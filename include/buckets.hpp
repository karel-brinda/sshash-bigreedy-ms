#pragma once

#include "util.hpp"
#include "bit_vector_iterator.hpp"
#include "ef_sequence.hpp"

namespace sshash {

/*
    Strings layout.

    SPSS input (as in SSHash v3): the strings are concatenated; pieces = [0, e_0, ..., e_{p-1}]
    holds their endpoints, so string i is [pieces[i], pieces[i+1]), and the k-mer at offset x of
    string i has id x - i * (k - 1). cum is empty.

    Masked-superstring input: the superstring is stored once and the strings are its maximal
    intervals [b_i, e_i) of represented k-mers (see builder/masked_superstring.hpp), which are
    sorted and may overlap by less than k - 1 bases. pieces = [0, e_0, ..., e_{p-1}] holds their
    ends, and cum = [0, c_1, ..., c_p] their cumulative numbers of k-mers, so that
        b_i = e_i - (k - 1) - (c_{i+1} - c_i)
    and the k-mer at offset x of interval i has id c_i + (x - b_i), i.e., ids are the ranks of the
    uppercase positions of the superstring. Since a super-k-mer lies in the interval that generated
    it, and every k-mer occurrence lies in at most one interval, the interval of a super-k-mer
    starting at offset x is the one with the smallest end >= x + k.
*/
struct buckets {
    bool masked_superstring() const { return cum.size() != 0; }

    /* Masked superstring only: [begin, end) and number of k-mers of interval i. */
    std::tuple<uint64_t, uint64_t, uint64_t> interval(uint64_t i, uint64_t k) const {
        assert(masked_superstring());
        uint64_t end = pieces.access(i + 1);
        uint64_t num_kmers = cum.access(i + 1) - cum.access(i);
        return {end - (k - 1) - num_kmers, end, num_kmers};
    }

    std::pair<lookup_result, uint64_t> offset_to_id(uint64_t offset, uint64_t k) const {
        if (masked_superstring()) return offset_to_id_masked_superstring(offset, k);
        auto [pos, contig_begin, contig_end] = pieces.locate(offset);

        /* The following two facts hold. */
        assert(pieces.access(pos) == contig_end);
        assert(contig_end >= offset);

        bool shift = contig_end > offset;
        uint64_t contig_id = pos - shift;

        if (contig_end == offset) {
            assert(shift == false);
            assert(pos + 1 < pieces.size());
            contig_begin = contig_end;
            contig_end = pieces.access(pos + 1);
        }

        /* Now, the following facts hold. */
        assert(offset >= contig_id * (k - 1));
        assert(contig_begin <= offset);
        assert(offset < contig_end);

        uint64_t absolute_kmer_id = offset - contig_id * (k - 1);
        uint64_t relative_kmer_id = offset - contig_begin;
        uint64_t contig_length = contig_end - contig_begin;
        assert(contig_length >= k);
        uint64_t contig_size = contig_length - k + 1;

        lookup_result res;
        res.kmer_id = absolute_kmer_id;
        res.kmer_id_in_contig = relative_kmer_id;
        res.contig_id = contig_id;
        res.contig_size = contig_size;

        return {res, contig_end};
    }

    std::pair<lookup_result, uint64_t> offset_to_id_masked_superstring(uint64_t offset,
                                                                       uint64_t k) const {
        auto [pos, contig_end] = pieces.next_geq(offset + k);
        assert(pos > 0 and pos < pieces.size());
        assert(contig_end >= offset + k);
        uint64_t contig_id = pos - 1;
        uint64_t first_kmer_id = cum.access(contig_id);
        uint64_t contig_size = cum.access(contig_id + 1) - first_kmer_id;
        uint64_t contig_begin = contig_end - (k - 1) - contig_size;
        assert(contig_begin <= offset);

        lookup_result res;
        res.kmer_id_in_contig = offset - contig_begin;
        res.kmer_id = first_kmer_id + res.kmer_id_in_contig;
        res.contig_id = contig_id;
        res.contig_size = contig_size;

        return {res, contig_end};
    }

    uint64_t contig_length(uint64_t contig_id, uint64_t k) const {
        if (masked_superstring()) return cum.access(contig_id + 1) - cum.access(contig_id) + k - 1;
        uint64_t length = pieces.access(contig_id + 1) - pieces.access(contig_id);
        return length;
    }

    uint64_t contig_prefix(uint64_t contig_id, uint64_t k) const {
        uint64_t contig_begin = masked_superstring() ? std::get<0>(interval(contig_id, k))
                                                     : pieces.access(contig_id);
        bit_vector_iterator bv_it(strings, 2 * contig_begin);
        return bv_it.read(2 * (k - 1));
    }

    uint64_t contig_suffix(uint64_t contig_id, uint64_t k) const {
        uint64_t contig_end = pieces.access(contig_id + 1);
        bit_vector_iterator bv_it(strings, 2 * (contig_end - k + 1));
        return bv_it.read(2 * (k - 1));
    }

    std::pair<uint64_t, uint64_t> locate_bucket(uint64_t bucket_id) const {
        uint64_t begin = num_super_kmers_before_bucket.access(bucket_id) + bucket_id;
        uint64_t end = num_super_kmers_before_bucket.access(bucket_id + 1) + bucket_id + 1;
        assert(begin < end);
        return {begin, end};
    }

    lookup_result lookup(uint64_t bucket_id, uint64_t target_kmer, uint64_t k, uint64_t m) const {
        auto [begin, end] = locate_bucket(bucket_id);
        return lookup(begin, end, target_kmer, k, m);
    }

    lookup_result lookup(uint64_t begin, uint64_t end, uint64_t target_kmer, uint64_t k,
                         uint64_t m) const {
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            auto res = lookup_in_super_kmer(super_kmer_id, target_kmer, k, m);
            if (res.kmer_id != constants::invalid_uint64) {
                assert(is_valid(res));
                return res;
            }
        }
        return lookup_result();
    }

    lookup_result lookup_in_super_kmer(uint64_t super_kmer_id, uint64_t target_kmer, uint64_t k,
                                       uint64_t m) const {
        uint64_t offset = offsets.access(super_kmer_id);
        auto [res, contig_end] = offset_to_id(offset, k);
        bit_vector_iterator bv_it(strings, 2 * offset);
        uint64_t window_size = std::min<uint64_t>(k - m + 1, contig_end - offset - k + 1);
        for (uint64_t w = 0; w != window_size; ++w) {
            uint64_t read_kmer = bv_it.read_and_advance_by_two(2 * k);
            if (read_kmer == target_kmer) {
                res.kmer_id += w;
                res.kmer_id_in_contig += w;
                assert(is_valid(res));
                return res;
            }
        }
        return lookup_result();
    }

    lookup_result lookup_canonical(uint64_t bucket_id, uint64_t target_kmer,
                                   uint64_t target_kmer_rc, uint64_t k, uint64_t m) const {
        auto [begin, end] = locate_bucket(bucket_id);
        return lookup_canonical(begin, end, target_kmer, target_kmer_rc, k, m);
    }

    lookup_result lookup_canonical(uint64_t begin, uint64_t end, uint64_t target_kmer,
                                   uint64_t target_kmer_rc, uint64_t k, uint64_t m) const {
        for (uint64_t super_kmer_id = begin; super_kmer_id != end; ++super_kmer_id) {
            uint64_t offset = offsets.access(super_kmer_id);
            auto [res, contig_end] = offset_to_id(offset, k);
            bit_vector_iterator bv_it(strings, 2 * offset);
            uint64_t window_size = std::min<uint64_t>(k - m + 1, contig_end - offset - k + 1);
            for (uint64_t w = 0; w != window_size; ++w) {
                uint64_t read_kmer = bv_it.read_and_advance_by_two(2 * k);
                if (read_kmer == target_kmer) {
                    res.kmer_id += w;
                    res.kmer_id_in_contig += w;
                    res.kmer_orientation = constants::forward_orientation;
                    assert(is_valid(res));
                    return res;
                }
                if (read_kmer == target_kmer_rc) {
                    res.kmer_id += w;
                    res.kmer_id_in_contig += w;
                    res.kmer_orientation = constants::backward_orientation;
                    assert(is_valid(res));
                    return res;
                }
            }
        }
        return lookup_result();
    }

    uint64_t id_to_offset(uint64_t id, uint64_t k) const {
        if (masked_superstring()) {
            uint64_t contig_id = cum.prev_leq(id);  // c_i <= id < c_{i+1}
            assert(contig_id + 1 < cum.size());
            auto [begin, end, num_kmers] = interval(contig_id, k);
            (void)end;
            (void)num_kmers;
            return begin + (id - cum.access(contig_id));
        }
        constexpr uint64_t linear_scan_threshold = 8;
        uint64_t lo = 0;
        uint64_t hi = pieces.size() - 1;
        assert(pieces.access(0) == 0);
        while (lo < hi) {
            if (hi - lo <= linear_scan_threshold) {
                for (; lo < hi; ++lo) {
                    uint64_t val = pieces.access(lo) - lo * (k - 1);
                    if (val > id) break;
                }
                break;
            }
            uint64_t mid = lo + (hi - lo) / 2;
            uint64_t val = pieces.access(mid);
            assert(val >= mid * (k - 1));
            if (id <= val - mid * (k - 1)) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        if (lo < pieces.size() and pieces.access(lo) - lo * (k - 1) > id) --lo;
        return id + lo * (k - 1);
    }

    void access(uint64_t kmer_id, char* string_kmer, uint64_t k) const {
        uint64_t offset = id_to_offset(kmer_id, k);
        bit_vector_iterator bv_it(strings, 2 * offset);
        uint64_t read_kmer = bv_it.read(2 * k);
        util::uint64_to_string_no_reverse(read_kmer, string_kmer, k);
    }

    struct iterator {
        iterator() {}

        iterator(buckets const* ptr, uint64_t kmer_id, uint64_t k, uint64_t num_kmers)
            : m_buckets(ptr), m_kmer_id(kmer_id), m_k(k), m_num_kmers(num_kmers) {
            bv_it = bit_vector_iterator(m_buckets->strings, -1);
            offset = m_buckets->id_to_offset(m_kmer_id, k);
            if (m_buckets->masked_superstring()) {
                /* The interval of the k-mer at offset (see offset_to_id_masked_superstring). */
                auto [pos, piece_end] = m_buckets->pieces.next_geq(offset + k);
                (void)piece_end;
                contig_id = pos - 1;
                next_offset = piece_end;
                start_piece();
            } else {
                auto [pos, piece_end] = m_buckets->pieces.next_geq(offset);
                if (piece_end == offset) pos += 1;
                pieces_it = m_buckets->pieces.at(pos);
                next_piece();
            }
            ret.second.resize(k, 0);
        }

        bool has_next() const { return m_kmer_id != m_num_kmers; }

        std::pair<uint64_t, std::string> next() {
            if (offset == next_offset - m_k + 1) {
                if (m_buckets->masked_superstring()) {
                    /* Consecutive intervals may overlap: restart at the begin of the next one. */
                    auto [begin, end, num_kmers] = m_buckets->interval(++contig_id, m_k);
                    (void)num_kmers;
                    offset = begin;
                    next_offset = end;
                    start_piece();
                } else {
                    offset = next_offset;
                    next_piece();
                }
            }

            while (offset != next_offset - m_k + 1) {
                ret.first = m_kmer_id;
                if (clear) {
                    util::uint64_to_string_no_reverse(read_kmer, ret.second.data(), m_k);
                } else {
                    memmove(ret.second.data(), ret.second.data() + 1, m_k - 1);
                    ret.second[m_k - 1] = util::uint64_to_char(last_two_bits);
                }
                clear = false;
                read_kmer >>= 2;
                last_two_bits = bv_it.get_next_two_bits();
                read_kmer += last_two_bits << (2 * (m_k - 1));
                ++m_kmer_id;
                ++offset;
                return ret;
            }

            return next();
        }

    private:
        std::pair<uint64_t, std::string> ret;
        buckets const* m_buckets;
        uint64_t m_kmer_id, m_k, m_num_kmers;
        uint64_t offset;
        uint64_t next_offset;
        bit_vector_iterator bv_it;
        ef_sequence<true>::iterator pieces_it;

        uint64_t contig_id;  // masked superstring only

        uint64_t read_kmer;
        uint64_t last_two_bits;
        bool clear;

        void next_piece() {
            next_offset = pieces_it.next();
            start_piece();
        }

        void start_piece() {
            bv_it.at(2 * offset);
            assert(next_offset > offset);
            read_kmer = bv_it.take(2 * m_k);
            clear = true;
        }
    };

    iterator at(uint64_t kmer_id, uint64_t k, uint64_t size) const {
        return iterator(this, kmer_id, k, size);
    }

    uint64_t num_bits() const {
        return pieces.num_bits() + cum.num_bits() + num_super_kmers_before_bucket.num_bits() +
               8 * (offsets.bytes() + strings.bytes());
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(pieces);
        visitor.visit(cum);
        visitor.visit(num_super_kmers_before_bucket);
        visitor.visit(offsets);
        visitor.visit(strings);
    }

    ef_sequence<true> pieces;
    ef_sequence<true> cum;  // empty unless the input is a masked superstring
    ef_sequence<false> num_super_kmers_before_bucket;
    pthash::compact_vector offsets;
    pthash::bit_vector strings;

private:
    bool is_valid(lookup_result res) const {
        return (res.contig_size == constants::invalid_uint32 or
                res.kmer_id_in_contig < res.contig_size) and
               (res.contig_id == constants::invalid_uint32 or res.contig_id < pieces.size());
    }
};

}  // namespace sshash