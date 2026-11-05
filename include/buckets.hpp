#pragma once

#include "external/pthash/external/bits/include/elias_fano.hpp"

#include "util.hpp"
#include "kmer_iterator.hpp"

namespace sshash {

/*
    Strings layout.

    SPSS input (as in SSHash v4): the strings are concatenated; pieces = [0, e_0, ..., e_{p-1}]
    holds their endpoints, so string i is [pieces[i], pieces[i+1]), and the k-mer at offset x of
    string i has id x - i * (k - 1). cum is empty.

    Masked-superstring input: the superstring is stored once and the strings are its maximal
    intervals [b_i, e_i) of represented k-mers (see builder/masked_superstring.hpp), which are
    sorted and may overlap by less than k - 1 bases. pieces = [0, e_0, ..., e_{p-1}] holds their
    ends, and cum = [0, c_1, ..., c_p] their cumulative numbers of k-mers, so that
        b_i = e_i - (k - 1) - (c_{i+1} - c_i)
    and the k-mer at offset x of interval i has id c_i + (x - b_i), i.e., ids are the ranks of the
    uppercase positions of the superstring. Every k-mer occurrence lies in at most one interval:
    the k-mer at offset x is represented if and only if x >= b_i for the interval i with the
    smallest end >= x + k. Minimizer positions stay absolute: a position shared by the
    super-k-mers of two overlapping intervals is stored once, which is enough since a lookup
    derives the k-mer offset from the query and then finds its interval as above.
*/
template <class kmer_t>
struct buckets  //
{
    bool masked_superstring() const { return cum.size() != 0; }

    lookup_result offset_to_id(const uint64_t offset, const uint64_t k) const {
        if (masked_superstring()) return offset_to_id_masked_superstring(offset, k);
        auto p = pieces.locate(offset);
        uint64_t contig_id = p.first.pos;
        uint64_t contig_begin = p.first.val;
        uint64_t contig_end = p.second.val;

        /* The following facts hold. */
        assert(offset >= contig_id * (k - 1));
        assert(contig_begin <= offset);
        assert(offset < contig_end);
        /****************************/

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
        assert(contig_begin == res.contig_begin(k));
        assert(contig_end == res.contig_end(k));

        return res;
    }

    /* Return an invalid result if the k-mer at offset is not represented. */
    lookup_result offset_to_id_masked_superstring(const uint64_t offset, const uint64_t k) const {
        auto [pos, contig_end] = pieces.next_geq(offset + k);
        if (contig_end < offset + k) return lookup_result();  // past the last interval
        assert(pos > 0);
        uint64_t contig_id = pos - 1;
        uint64_t contig_size = cum.diff(contig_id);
        uint64_t contig_begin = contig_end - (k - 1) - contig_size;
        if (offset < contig_begin) return lookup_result();  // not in any interval

        lookup_result res;
        res.kmer_id_in_contig = offset - contig_begin;
        res.kmer_id = cum.access(contig_id) + res.kmer_id_in_contig;
        res.contig_id = contig_id;
        res.contig_size = contig_size;
        return res;
    }

    /* Whether res, returned by offset_to_id(offset, k), is the k-mer at offset. */
    bool is_kmer_at(lookup_result const& res, const uint64_t offset, const uint64_t k) const {
        if (masked_superstring()) return res.kmer_id != constants::invalid_uint64;
        return offset + k - 1 < res.contig_end(k);
    }

    /* Where the k-mer of a valid lookup result begins in strings. */
    uint64_t kmer_offset(lookup_result const& res, const uint64_t k) const {
        if (masked_superstring()) {
            return contig_offsets(res.contig_id, k).first + res.kmer_id_in_contig;
        }
        return res.kmer_id + res.contig_id * (k - 1);
    }

    uint64_t contig_first_kmer_id(const uint64_t contig_id, const uint64_t k) const {
        if (masked_superstring()) return cum.access(contig_id);
        return pieces.access(contig_id) - contig_id * (k - 1);
    }

    /* Return where the contig begins and ends in strings. */
    std::pair<uint64_t, uint64_t>  // [begin, end)
    contig_offsets(const uint64_t contig_id, const uint64_t k) const {
        uint64_t end = pieces.access(contig_id + 1);
        uint64_t begin = masked_superstring() ? end - (k - 1) - cum.diff(contig_id)
                                              : pieces.access(contig_id);
        assert(end > begin);
        return {begin, end};
    }

    kmer_t contig_prefix(const uint64_t contig_id, const uint64_t k) const {
        uint64_t contig_begin = contig_offsets(contig_id, k).first;
        return util::read_kmer_at<kmer_t>(strings, k - 1, kmer_t::bits_per_char * contig_begin);
    }

    kmer_t contig_suffix(const uint64_t contig_id, const uint64_t k) const {
        uint64_t contig_end = pieces.access(contig_id + 1);
        return util::read_kmer_at<kmer_t>(strings, k - 1,
                                          kmer_t::bits_per_char * (contig_end - k + 1));
    }

    std::pair<uint64_t, uint64_t> locate_bucket(const uint64_t bucket_id) const {
        uint64_t begin = bucket_sizes.access(bucket_id) + bucket_id;
        uint64_t end = bucket_sizes.access(bucket_id + 1) + bucket_id + 1;
        assert(begin < end);
        return {begin, end};
    }

    lookup_result lookup(uint64_t begin, uint64_t end, kmer_t kmer, minimizer_info mini_info,
                         const uint64_t k, const uint64_t m) const  //
    {
        { /* check minimizer first */
            uint64_t pos_in_seq = offsets.access(begin);
            uint64_t read_mmer = uint64_t(
                util::read_kmer_at<kmer_t>(strings, m, kmer_t::bits_per_char * pos_in_seq));
            if (read_mmer != mini_info.minimizer) {
                auto res = lookup_result();
                res.minimizer_found = false;
                return res;
            }
        }

        for (uint64_t i = begin; i != end; ++i) {
            auto res = lookup(i, kmer, mini_info, k, m);
            if (res.kmer_id != constants::invalid_uint64) {
                assert(is_valid(res));
                return res;
            }
        }

        return lookup_result();
    }

    lookup_result lookup(uint64_t i, kmer_t kmer, minimizer_info mini_info, const uint64_t k,
                         const uint64_t m) const  //
    {
        (void)m;
        uint64_t pos_in_seq = offsets.access(i);
        if (pos_in_seq >= mini_info.pos_in_kmer) {
            uint64_t offset = pos_in_seq - mini_info.pos_in_kmer;
            auto res = offset_to_id(offset, k);
            if (is_kmer_at(res, offset, k)) {
                auto read_kmer =
                    util::read_kmer_at<kmer_t>(strings, k, kmer_t::bits_per_char * offset);
                if (read_kmer == kmer) {
                    assert(is_valid(res));
                    return res;
                }
            }
        }
        return lookup_result();
    }

    lookup_result lookup_canonical(uint64_t begin, uint64_t end, kmer_t kmer, kmer_t kmer_rc,
                                   minimizer_info mini_info, const uint64_t k,
                                   const uint64_t m) const  //
    {
        { /* check minimizer first */
            uint64_t pos_in_seq = offsets.access(begin);
            uint64_t read_mmer = uint64_t(
                util::read_kmer_at<kmer_t>(strings, m, kmer_t::bits_per_char * pos_in_seq));
            auto tmp = kmer_t(mini_info.minimizer);
            tmp.reverse_complement_inplace(m);
            uint64_t minimizer_rc = uint64_t(tmp);
            if (read_mmer != mini_info.minimizer and read_mmer != minimizer_rc) {
                auto res = lookup_result();
                res.minimizer_found = false;
                return res;
            }
        }

        for (uint64_t i = begin; i != end; ++i) {
            auto res = lookup_canonical(i, kmer, kmer_rc, mini_info, k, m);
            if (res.kmer_id != constants::invalid_uint64) {
                assert(is_valid(res));
                return res;
            }
        }

        return lookup_result();
    }

    lookup_result lookup_canonical(uint64_t i, kmer_t kmer, kmer_t kmer_rc,
                                   minimizer_info mini_info, const uint64_t k,
                                   const uint64_t m) const  //
    {
        const uint64_t pos_in_seq = offsets.access(i);
        uint64_t pos_in_kmer = mini_info.pos_in_kmer;
        auto res = check_position(pos_in_seq, pos_in_kmer, kmer, kmer_rc, k);
        if (res.kmer_id != constants::invalid_uint64) {
            assert(is_valid(res));
            return res;
        }
        pos_in_kmer = k - m - mini_info.pos_in_kmer;
        return check_position(pos_in_seq, pos_in_kmer, kmer, kmer_rc, k);
    }

    lookup_result check_position(const uint64_t pos_in_seq, const uint64_t pos_in_kmer,  //
                                 const kmer_t kmer, const kmer_t kmer_rc,                //
                                 const uint64_t k) const                                 //
    {
        if (pos_in_seq >= pos_in_kmer) {
            uint64_t offset = pos_in_seq - pos_in_kmer;
            auto res = offset_to_id(offset, k);
            if (is_kmer_at(res, offset, k)) {
                auto read_kmer =
                    util::read_kmer_at<kmer_t>(strings, k, kmer_t::bits_per_char * offset);
                if (read_kmer == kmer) {
                    assert(is_valid(res));
                    return res;
                }
                if (read_kmer == kmer_rc) {
                    assert(is_valid(res));
                    res.kmer_orientation = constants::backward_orientation;
                    return res;
                }
            }
        }
        return lookup_result();
    }

    uint64_t id_to_offset(const uint64_t id, const uint64_t k) const {
        if (masked_superstring()) {
            auto [contig_id, first_kmer_id] = cum.prev_leq(id);  // c_i <= id < c_{i+1}
            return contig_offsets(contig_id, k).first + (id - first_kmer_id);
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

    void access(const uint64_t kmer_id, char* string_kmer, const uint64_t k) const {
        uint64_t offset = id_to_offset(kmer_id, k);
        auto read_kmer = util::read_kmer_at<kmer_t>(strings, k, kmer_t::bits_per_char * offset);
        util::uint_kmer_to_string(read_kmer, string_kmer, k);
    }

    struct iterator {
        iterator() {}

        iterator(buckets const* ptr,                                        //
                 const uint64_t begin_kmer_id, const uint64_t end_kmer_id,  // [begin,end)
                 const uint64_t k)
            : m_buckets(ptr)
            , m_begin_kmer_id(begin_kmer_id)
            , m_end_kmer_id(end_kmer_id)
            , m_k(k)
            , m_it(ptr->strings, m_k)  //
        {
            m_offset = m_buckets->id_to_offset(m_begin_kmer_id, k);
            if (m_buckets->masked_superstring()) {
                /* the interval of the k-mer at m_offset (see offset_to_id_masked_superstring) */
                auto [pos, piece_end] = m_buckets->pieces.next_geq(m_offset + k);
                m_contig_id = pos - 1;
                m_next_offset = piece_end;
                start_piece();
            } else {
                auto [pos, piece_end] = m_buckets->pieces.next_geq(m_offset);
                if (piece_end == m_offset) pos += 1;
                m_pieces_it = m_buckets->pieces.get_iterator_at(pos);
                next_piece();
            }
            m_ret.second.resize(m_k, 0);
        }

        bool has_next() const { return m_begin_kmer_id != m_end_kmer_id; }

        std::pair<uint64_t, std::string> next() {
            if (m_offset == m_next_offset - m_k + 1) {
                if (m_buckets->masked_superstring()) {
                    /* Consecutive intervals may overlap: restart at the begin of the next one. */
                    auto [begin, end] = m_buckets->contig_offsets(++m_contig_id, m_k);
                    m_offset = begin;
                    m_next_offset = end;
                    start_piece();
                } else {
                    m_offset = m_next_offset;
                    next_piece();
                }
            }
            m_ret.first = m_begin_kmer_id;
            if (m_clear) {
                util::uint_kmer_to_string(m_it.get(), m_ret.second.data(), m_k);
                assert(kmer_t::bits_per_char * m_offset == m_it.position());
                m_it.at(kmer_t::bits_per_char * (m_offset + m_k));
            } else {
                memmove(m_ret.second.data(), m_ret.second.data() + 1, m_k - 1);
                m_ret.second[m_k - 1] = kmer_t::uint64_to_char(m_it.get_next_char());
            }
            m_clear = false;
            ++m_begin_kmer_id;
            ++m_offset;
            return m_ret;
        }

    private:
        std::pair<uint64_t, std::string> m_ret;
        buckets const* m_buckets;
        uint64_t m_begin_kmer_id, m_end_kmer_id;
        uint64_t m_k;
        uint64_t m_offset;
        uint64_t m_next_offset;
        kmer_iterator<kmer_t> m_it;
        bits::elias_fano<true, false>::iterator m_pieces_it;
        uint64_t m_contig_id;  // masked superstring only
        bool m_clear;

        void next_piece() {
            m_next_offset = m_pieces_it.value();
            m_pieces_it.next();
            start_piece();
        }

        void start_piece() {
            m_it.at(kmer_t::bits_per_char * m_offset);
            assert(m_next_offset > m_offset);
            m_clear = true;
        }
    };

    iterator at(const uint64_t begin_kmer_id, const uint64_t end_kmer_id, const uint64_t k) const {
        return iterator(this, begin_kmer_id, end_kmer_id, k);
    }

    uint64_t num_bits() const {
        return 8 * (pieces.num_bytes() + cum.num_bytes() + bucket_sizes.num_bytes() +
                    offsets.num_bytes() + strings.num_bytes());
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

    bits::elias_fano<true, false> pieces;
    bits::elias_fano<true, true> cum;  // empty unless the input is a masked superstring
    bits::elias_fano<false, false> bucket_sizes;
    bits::compact_vector offsets;
    bits::bit_vector strings;

private:
    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.pieces);
        visitor.visit(t.cum);
        visitor.visit(t.bucket_sizes);
        visitor.visit(t.offsets);
        visitor.visit(t.strings);
    }

    bool is_valid(lookup_result res) const {
        return res.contig_size != constants::invalid_uint64 and  //
               res.kmer_id_in_contig < res.contig_size and       //
               res.contig_id < pieces.size();                    //
    }
};

}  // namespace sshash