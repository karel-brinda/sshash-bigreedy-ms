#pragma once

#include "external/pthash/external/bits/include/compact_vector.hpp"
#include "external/pthash/external/bits/include/endpoints_sequence.hpp"

namespace sshash {

struct num_bits {
    num_bits() : per_absolute_offset(0), per_relative_offset(0), per_string_id(0) {}
    uint64_t per_absolute_offset;
    uint64_t per_relative_offset;
    uint64_t per_string_id;
};

template <typename Seq>
struct offsets  //
{
    struct builder {
        builder() {}

        void reserve(uint64_t n) { m_v.reserve(n); }
        void push_back(uint64_t val) { m_v.push_back(val); }

        uint64_t operator[](uint64_t i) {
            assert(i < m_v.size());
            return m_v[i];
        }

        uint64_t front() const { return m_v.front(); }
        uint64_t back() const { return m_v.back(); }
        uint64_t size() const { return m_v.size(); }
        uint64_t num_bytes() const { return sizeof(m_nb) + essentials::vec_bytes(m_v); }

        void set_num_bits(num_bits nb) { m_nb = nb; }
        void set_k(uint64_t /* k */) {}

        uint64_t num_strings() const {
            assert(m_v.size() >= 2);
            return m_v.size() - 1;
        }
        uint64_t string_begin(uint64_t i) const { return m_v[i]; }
        uint64_t string_end(uint64_t i) const { return m_v[i + 1]; }

    protected:
        num_bits m_nb;
        std::vector<uint64_t> m_v;
    };

    uint64_t string_begin(uint64_t string_id) const { return m_seq.access(string_id); }
    uint64_t string_end(uint64_t string_id) const { return m_seq.access(string_id + 1); }

    /* First kmer-id of the string: strings partition `strings`, so string j
       begins at offset access(j) and j * (k - 1) offsets do not start a kmer. */
    uint64_t first_kmer_id(const uint64_t string_id, const uint64_t k) const {
        return m_seq.access(string_id) - string_id * (k - 1);
    }

    /* Fill the kmer-ids of res from its already-located string fields. */
    void set_ids(lookup_result& res, const uint64_t k) const {
        res.kmer_id = res.kmer_offset - res.string_id * (k - 1);
        res.kmer_id_in_string = res.kmer_offset - res.string_begin;
    }

    std::pair<uint64_t, uint64_t> id_to_offset(const uint64_t kmer_id, const uint64_t k) const  //
    {
        constexpr uint64_t linear_scan_threshold = 32;
        uint64_t lo = 0;
        uint64_t hi = m_seq.size() - 1;
        assert(m_seq.access(0) == 0);
        while (hi - lo > linear_scan_threshold) {
            uint64_t mid = lo + (hi - lo) / 2;
            uint64_t val = m_seq.access(mid);
            assert(val >= mid * (k - 1));
            if (kmer_id <= val - mid * (k - 1)) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        assert(lo < hi);
        assert(hi < m_seq.size());
        for (auto it = m_seq.get_iterator_at(lo); lo <= hi; ++lo, ++it) {
            uint64_t val = *it - lo * (k - 1);
            if (val > kmer_id) break;
        }
        assert(lo > 0);
        /* return [string_id, offset] */
        return {lo - 1, kmer_id + (lo - 1) * (k - 1)};
    }

    uint64_t access(uint64_t i) const {
        assert(i < size());
        return m_seq.access(i);
    }

    uint64_t size() const { return m_seq.size(); }

    uint64_t num_bytes() const {
        return sizeof(m_num_bits_per_relative_offset) + m_seq.num_bytes();
    }

    struct iterator {
        iterator() {}
        iterator(offsets const* e, uint64_t pos) { m_it = e->m_seq.get_iterator_at(pos); }

        uint64_t value() const { return *m_it; }
        void next() { ++m_it; }

    private:
        typename Seq::iterator m_it;
    };

    iterator get_iterator_at(uint64_t pos) const {
        assert(pos < size());
        return {this, pos};
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

protected:
    Seq m_seq;
    uint64_t m_num_bits_per_relative_offset;

    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_num_bits_per_relative_offset);
        visitor.visit(t.m_seq);
    }
};

struct decoded_offsets
    : public offsets<bits::endpoints_sequence<>>  //
{
    struct builder
        : public offsets<bits::endpoints_sequence<>>::builder  //
    {
        uint64_t num_bits_per_offset() const { return m_nb.per_absolute_offset; }

        uint64_t encode(uint64_t offset, uint64_t, uint64_t) { return offset; }

        void build(decoded_offsets& e) {
            assert(std::is_sorted(m_v.begin(), m_v.end()));
            e.m_seq.encode(m_v.begin(), m_v.size(), m_v.back());
            std::vector<uint64_t>().swap(m_v);
        }
    };

    struct decoded_offset {
        uint64_t absolute_offset;
    };

    decoded_offset decode(const uint64_t encoded_offset) const { return {encoded_offset}; }

    void offset_to_id(lookup_result& res, decoded_offset /* p */, const uint64_t k) const  //
    {
        assert(res.kmer_offset != constants::invalid_uint64);

        auto q = m_seq.locate(res.kmer_offset);
        res.string_id = q.first.pos;
        res.string_begin = q.first.val;
        res.string_end = q.second.val;
        res.kmer_id = res.kmer_offset - res.string_id * (k - 1);     // absolute kmer id
        res.kmer_id_in_string = res.kmer_offset - res.string_begin;  // relative kmer id

        assert(res.string_id < m_seq.size());
        assert(res.string_begin < res.string_end);
        assert(res.kmer_offset >= res.string_id * (k - 1));
        assert(res.string_begin <= res.kmer_offset);
        assert(res.kmer_offset < res.string_end);
    }
};

struct encoded_offsets
    : public offsets<bits::compact_vector>  //
{
    struct builder
        : public offsets<bits::compact_vector>::builder  //
    {
        uint64_t num_bits_per_offset() const {
            return m_nb.per_string_id + m_nb.per_relative_offset;
        }

        uint64_t encode(uint64_t offset, uint64_t begin, uint64_t string_id) {
            /* encode offset as | string-id | relative offset | */
            assert(string_id < m_v.size());
            assert(offset >= begin);
            assert((offset - begin) < (1ULL << m_nb.per_relative_offset));
            uint64_t relative_offset = offset - begin;
            return (string_id << m_nb.per_relative_offset) + relative_offset;
        }

        void build(encoded_offsets& e) {
            assert(std::is_sorted(m_v.begin(), m_v.end()));
            e.m_seq.build(m_v.begin(), m_v.size(), m_nb.per_absolute_offset);
            e.m_num_bits_per_relative_offset = m_nb.per_relative_offset;
            std::vector<uint64_t>().swap(m_v);
        }
    };

    struct decoded_offset {
        uint64_t absolute_offset;
        uint64_t relative_offset;
        uint64_t string_id;
        uint64_t string_begin;
        uint64_t string_end;
    };

    decoded_offset decode(const uint64_t encoded_offset) const {
        uint64_t relative_offset = encoded_offset & ((1ULL << m_num_bits_per_relative_offset) - 1);
        uint64_t string_id = encoded_offset >> m_num_bits_per_relative_offset;
        assert(string_id + 1 < m_seq.size());
        uint64_t begin = m_seq.access(string_id);
        uint64_t end = m_seq.access(string_id + 1);
        return {begin + relative_offset, relative_offset, string_id, begin, end};
    }

    void offset_to_id(lookup_result& res, decoded_offset p, const uint64_t k) const  //
    {
        assert(res.kmer_offset != constants::invalid_uint64);

        res.string_id = p.string_id;
        res.string_begin = p.string_begin;
        res.string_end = p.string_end;
        res.kmer_id = res.kmer_offset - res.string_id * (k - 1);     // absolute kmer id
        res.kmer_id_in_string = res.kmer_offset - res.string_begin;  // relative kmer id

        assert(res.string_id < m_seq.size());
        assert(res.string_begin < res.string_end);
    }
};

/*
    Offsets for strings that are embedded, possibly overlapping, intervals of
    one single text, as for a BiGreedy masked superstring, where string j is
    the j-th maximal run of represented kmer starts: string j spells
    text[begin_j..end_j), intervals are sorted by begin, have length at least
    k, and consecutive intervals overlap by fewer than k characters, so every
    kmer occurrence of the text lies in at most one interval. Only the kmers
    of the intervals are indexed, and global kmer-ids are assigned in interval
    order. A plain SPSS is the degenerate case begin_{j+1} == end_j.
    Minimizer offsets are encoded as | string-id | relative offset |, like for
    encoded_offsets, so that the same text position remains distinguishable
    when it anchors kmers of two adjacent overlapping intervals.
*/
struct interval_offsets  //
{
    struct builder {
        builder() : m_k(0) {}

        void reserve(uint64_t n) { m_endpoints.reserve(n); }

        void set_k(uint64_t k) { m_k = k; }

        /* SPSS input: the begin offset of each string, plus the total length
           at the end, exactly as for the other builders. */
        void push_back(uint64_t val) {
            assert(m_begins.empty() and m_ends.empty());
            m_endpoints.push_back(val);
        }

        /* Masked-superstring input: one explicit interval per string. */
        void push_interval(uint64_t begin, uint64_t end) {
            assert(m_endpoints.empty());
            assert(end > begin);
            assert(m_begins.empty() or begin > m_begins.back());
            assert(m_ends.empty() or end > m_ends.back());
            m_begins.push_back(begin);
            m_ends.push_back(end);
        }

        uint64_t num_strings() const {
            return m_endpoints.empty() ? m_begins.size() : m_endpoints.size() - 1;
        }
        uint64_t string_begin(uint64_t i) const {
            return m_endpoints.empty() ? m_begins[i] : m_endpoints[i];
        }
        uint64_t string_end(uint64_t i) const {
            return m_endpoints.empty() ? m_ends[i] : m_endpoints[i + 1];
        }

        uint64_t operator[](uint64_t i) { return string_begin(i); }

        uint64_t front() const {
            return m_endpoints.empty() ? m_begins.front() : m_endpoints.front();
        }
        uint64_t back() const { return m_endpoints.empty() ? m_ends.back() : m_endpoints.back(); }
        uint64_t size() const { return m_endpoints.empty() ? m_begins.size() : m_endpoints.size(); }

        uint64_t num_bytes() const {
            return sizeof(m_nb) + sizeof(m_k) + essentials::vec_bytes(m_endpoints) +
                   essentials::vec_bytes(m_begins) + essentials::vec_bytes(m_ends);
        }

        void set_num_bits(num_bits nb) { m_nb = nb; }

        uint64_t num_bits_per_offset() const {
            return m_nb.per_string_id + m_nb.per_relative_offset;
        }

        /* encode offset as | string-id | relative offset | */
        uint64_t encode(uint64_t offset, uint64_t begin, uint64_t string_id) {
            assert(offset >= begin);
            uint64_t relative_offset = offset - begin;
            assert(m_nb.per_relative_offset == 64 or
                   relative_offset < (1ULL << m_nb.per_relative_offset));
            return (string_id << m_nb.per_relative_offset) + relative_offset;
        }

        void build(interval_offsets& e) {
            assert(m_k > 0);
            const uint64_t p = num_strings();
            assert(p > 0);

            if (!m_endpoints.empty()) { /* SPSS input: degenerate intervals */
                assert(m_endpoints.front() == 0);
                assert(std::is_sorted(m_endpoints.begin(), m_endpoints.end()));
                m_begins.assign(m_endpoints.begin(), m_endpoints.end() - 1);
                m_ends.assign(m_endpoints.begin() + 1, m_endpoints.end());
                std::vector<uint64_t>().swap(m_endpoints);
            }
            assert(m_begins.size() == p and m_ends.size() == p);

            {
                const uint64_t width =
                    std::max<uint64_t>(bits::util::ceil_log2_uint64(m_ends.back() + 1), 1);
                bits::compact_vector::builder cvb(m_begins.begin(), p, width);
                cvb.build(e.m_begins);
            }
            {
                /* prepend a zero, as endpoints_sequence requires a leading
                   zero: string j then ends at position j + 1 */
                std::vector<uint64_t> ends;
                ends.reserve(p + 1);
                ends.push_back(0);
                ends.insert(ends.end(), m_ends.begin(), m_ends.end());
                e.m_ends.encode(ends.begin(), ends.size(), ends.back());
            }
            {
                /* cum[j] = number of kmers of the strings before string j */
                std::vector<uint64_t> cum;
                cum.reserve(p + 1);
                cum.push_back(0);
                for (uint64_t j = 0; j != p; ++j) {
                    assert(m_ends[j] - m_begins[j] >= m_k);
                    cum.push_back(cum.back() + (m_ends[j] - m_begins[j] - m_k + 1));
                }
                e.m_cum.encode(cum.begin(), cum.size(), cum.back());
            }

            e.m_num_bits_per_relative_offset = m_nb.per_relative_offset;

            std::vector<uint64_t>().swap(m_begins);
            std::vector<uint64_t>().swap(m_ends);
        }

    private:
        num_bits m_nb;
        uint64_t m_k;
        std::vector<uint64_t> m_endpoints, m_begins, m_ends;
    };

    struct decoded_offset {
        uint64_t absolute_offset;
        uint64_t relative_offset;
        uint64_t string_id;
        uint64_t string_begin;
        uint64_t string_end;
    };

    uint64_t num_strings() const { return m_begins.size(); }
    uint64_t size() const { return m_begins.size() + 1; }

    uint64_t string_begin(uint64_t string_id) const {
        assert(string_id < num_strings());
        return m_begins.access(string_id);
    }
    uint64_t string_end(uint64_t string_id) const {
        assert(string_id < num_strings());
        return m_ends.access(string_id + 1);
    }

    uint64_t first_kmer_id(const uint64_t string_id, const uint64_t /* k */) const {
        return m_cum.access(string_id);
    }

    decoded_offset decode(const uint64_t encoded_offset) const {
        uint64_t relative_offset = encoded_offset & ((1ULL << m_num_bits_per_relative_offset) - 1);
        uint64_t string_id = encoded_offset >> m_num_bits_per_relative_offset;
        assert(string_id < num_strings());
        uint64_t begin = string_begin(string_id);
        uint64_t end = string_end(string_id);
        return {begin + relative_offset, relative_offset, string_id, begin, end};
    }

    /* Fill the kmer-ids of res from its already-located string fields. */
    void set_ids(lookup_result& res, const uint64_t /* k */) const {
        assert(res.kmer_offset >= res.string_begin);
        res.kmer_id_in_string = res.kmer_offset - res.string_begin;  // relative kmer id
        res.kmer_id = m_cum.access(res.string_id) + res.kmer_id_in_string;  // absolute kmer id
    }

    /*
        Locate the string that fully contains [res.kmer_offset,
        res.kmer_offset + k): since intervals are sorted by begin and overlap
        by fewer than k characters, it can only be the string with the
        smallest end >= res.kmer_offset + k. If no string contains the kmer,
        the string fields are invalidated so that the containment checks of
        the caller fail.
    */
    void offset_to_id(lookup_result& res, decoded_offset /* p */, const uint64_t k) const  //
    {
        assert(res.kmer_offset != constants::invalid_uint64);

        if (res.kmer_offset + k > m_ends.back()) {
            res.string_id = constants::invalid_uint64;
            res.string_begin = constants::invalid_uint64;
            res.string_end = constants::invalid_uint64;
            return;
        }

        auto q = m_ends.locate(res.kmer_offset + k - 1);
        assert(q.second.pos > 0);
        res.string_id = q.second.pos - 1;  // account for the leading zero
        res.string_begin = string_begin(res.string_id);
        res.string_end = q.second.val;
        assert(res.string_end == string_end(res.string_id));
        assert(res.kmer_offset + k <= res.string_end);

        if (res.kmer_offset >= res.string_begin) set_ids(res, k);
    }

    /* Return [string_id, offset] for the kmer with the given id. */
    std::pair<uint64_t, uint64_t> id_to_offset(const uint64_t kmer_id, const uint64_t k) const {
        if (kmer_id >= m_cum.back()) {
            /* out-of-range id: an empty iterator can still be constructed at
               the end of the last string (never dereferenced) */
            const uint64_t last = num_strings() - 1;
            return {last, string_end(last) - k + 1};
        }
        auto q = m_cum.locate(kmer_id);
        const uint64_t string_id = q.first.pos;
        assert(q.first.val <= kmer_id);
        assert(string_id < num_strings());
        return {string_id, string_begin(string_id) + (kmer_id - q.first.val)};
    }

    uint64_t access(uint64_t i) const { /* i-th string endpoint, as for the other flavours */
        return i == 0 ? (num_strings() ? string_begin(0) : 0) : string_end(i - 1);
    }

    uint64_t num_bytes() const {
        return sizeof(m_num_bits_per_relative_offset) + m_begins.num_bytes() +
               m_ends.num_bytes() + m_cum.num_bytes();
    }

    template <typename Visitor>
    void visit(Visitor& visitor) const {
        visit_impl(visitor, *this);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visit_impl(visitor, *this);
    }

private:
    uint64_t m_num_bits_per_relative_offset;
    bits::compact_vector m_begins;
    bits::endpoints_sequence<> m_ends;
    bits::endpoints_sequence<> m_cum;

    template <typename Visitor, typename T>
    static void visit_impl(Visitor& visitor, T&& t) {
        visitor.visit(t.m_num_bits_per_relative_offset);
        visitor.visit(t.m_begins);
        visitor.visit(t.m_ends);
        visitor.visit(t.m_cum);
    }
};

}  // namespace sshash
