#pragma once

/*
    Masked-superstring input, ported from 11_sshash-lite-bigreedy-ms (include/builder/bigreedy_ms.hpp).

    A masked superstring is one unwrapped FASTA record over ACGTacgt: an uppercase base marks the
    start of a represented k-mer. The strings indexed are the maximal intervals [begin, end) whose
    k-mers are all represented, i.e., each maximal uppercase run [first, last] gives the interval
    [first, last + k). Intervals are sorted, have length >= k, and consecutive ones overlap by
    less than k - 1 bases, so each k-mer occurrence of the superstring lies in at most one
    interval. K-mer ids are the ranks of the uppercase positions in text order.
*/

#include <cstdint>
#include <istream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sshash {

struct embedding_interval {
    uint64_t begin;
    uint64_t end;
};

inline bool is_masked_superstring_base(char base) {
    return base == 'A' or base == 'C' or base == 'G' or base == 'T' or base == 'a' or
           base == 'c' or base == 'g' or base == 't';
}

inline std::string read_masked_superstring(std::istream& is) {
    std::string header, sequence, extra;
    if (!std::getline(is, header) or header.empty() or header.front() != '>') {
        throw std::runtime_error("masked-superstring input must start with one FASTA header");
    }
    if (!std::getline(is, sequence) or sequence.empty()) {
        throw std::runtime_error("masked-superstring input must contain one unwrapped DNA sequence");
    }
    if (std::getline(is, extra)) {
        throw std::runtime_error(
            "masked-superstring input must contain exactly one unwrapped FASTA record");
    }
    for (char base : sequence) {
        if (!is_masked_superstring_base(base)) {
            throw std::runtime_error("masked-superstring input contains invalid DNA characters");
        }
    }
    return sequence;
}

inline std::vector<embedding_interval> masked_superstring_intervals(
    std::string const& masked_superstring, uint64_t k) {
    uint64_t n = masked_superstring.size();
    if (n < k) throw std::runtime_error("masked-superstring input length must be at least k");

    for (uint64_t pos = n - k + 1; pos != n; ++pos) {
        char base = masked_superstring[pos];
        if (base == 'A' or base == 'C' or base == 'G' or base == 'T') {
            throw std::runtime_error(
                "masked-superstring input has an uppercase character in its final k-1 positions");
        }
    }

    std::vector<embedding_interval> intervals;
    for (uint64_t pos = 0; pos <= n - k;) {
        if (masked_superstring[pos] >= 'a' and masked_superstring[pos] <= 'z') {
            ++pos;
            continue;
        }
        uint64_t begin = pos;
        do {
            ++pos;
        } while (pos <= n - k && masked_superstring[pos] >= 'A' &&
                 masked_superstring[pos] <= 'Z');
        intervals.push_back({begin, pos + k - 1});
    }
    if (intervals.empty()) {
        throw std::runtime_error("masked-superstring input has no represented k-mer starts");
    }
    return intervals;
}

}  // namespace sshash
