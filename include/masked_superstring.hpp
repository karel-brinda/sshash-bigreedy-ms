#pragma once

#include <istream>
#include <stdexcept>
#include <string>
#include <vector>

namespace sshash {

/*
    Support for BiGreedy masked superstrings in the KmerCamel format
    (as output by 'kmercamel compute'): one FASTA record with a single,
    unwrapped sequence line, where an uppercase base marks the start of a
    represented kmer (mask bit 1) and a lowercase base an unrepresented
    start (mask bit 0). The represented kmers are exactly the kmers that
    start at uppercase positions; with the default min-one mask, each
    represented kmer has exactly one uppercase occurrence, as required
    for the dictionary to assign unique kmer-ids.
*/

struct masked_superstring_interval {
    uint64_t begin;
    uint64_t end;
};

template <typename Kmer>
std::string read_masked_superstring(std::istream& is)  //
{
    std::string header, sequence, extra;
    if (!std::getline(is, header) or header.empty() or header.front() != '>') {
        throw std::runtime_error("masked-superstring input must start with one FASTA header");
    }
    if (!std::getline(is, sequence) or sequence.empty()) {
        throw std::runtime_error(
            "masked-superstring input must contain one unwrapped DNA sequence");
    }
    if (std::getline(is, extra)) {
        throw std::runtime_error(
            "masked-superstring input must contain exactly one unwrapped FASTA record");
    }
    for (char c : sequence) {
        if (!Kmer::is_valid(c)) {
            throw std::runtime_error("masked-superstring input contains invalid characters");
        }
    }
    return sequence;
}

[[maybe_unused]] static bool is_on_mask_char(char c) { return c >= 'A' and c <= 'Z'; }

/*
    Return the maximal runs of represented kmer starts, as intervals
    [begin, end) over the superstring: a run of uppercase positions
    [first, last] yields the interval [first, last + k), whose kmers are
    exactly the kmers starting at the run's positions. Intervals are
    sorted by begin, have length at least k, and consecutive intervals
    overlap by at most k - 2 characters (their uppercase runs are
    separated by at least one lowercase position), so each kmer
    occurrence of the superstring lies in at most one interval.
*/
[[maybe_unused]] static std::vector<masked_superstring_interval> masked_superstring_intervals(
    std::string const& masked_superstring, const uint64_t k)  //
{
    const uint64_t n = masked_superstring.size();
    if (n < k) throw std::runtime_error("masked-superstring input is shorter than k");

    for (uint64_t pos = n - k + 1; pos != n; ++pos) {
        if (is_on_mask_char(masked_superstring[pos])) {
            throw std::runtime_error(
                "masked-superstring input has an uppercase character in its final k-1 positions");
        }
    }

    std::vector<masked_superstring_interval> intervals;
    for (uint64_t pos = 0; pos <= n - k;) {
        if (!is_on_mask_char(masked_superstring[pos])) {
            ++pos;
            continue;
        }
        uint64_t begin = pos;
        do {
            ++pos;
        } while (pos <= n - k and is_on_mask_char(masked_superstring[pos]));
        intervals.push_back({begin, pos + k - 1});
    }
    if (intervals.empty()) {
        throw std::runtime_error("masked-superstring input has no represented kmer starts");
    }
    return intervals;
}

}  // namespace sshash
