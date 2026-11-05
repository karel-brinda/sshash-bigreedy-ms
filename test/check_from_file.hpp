#pragma once

#include <algorithm>  // for std::transform

#include "external/gz/zip_stream.hpp"
#include "include/builder/masked_superstring.hpp"

namespace sshash {

template <class kmer_t, input_file_type fmt>
bool check_correctness_lookup_access(std::istream& is, dictionary<kmer_t> const& dict) {
    const uint64_t k = dict.k();
    std::string sequence;
    uint64_t num_kmers = 0;
    uint64_t num_sequences = 0;
    lookup_result prev;
    prev.contig_id = 0;

    std::string got_kmer_str(k, 0);
    std::string expected_kmer_str(k, 0);

    std::cout << "checking correctness of access and positive lookup..." << std::endl;

    while (!is.eof())  //
    {
        if constexpr (fmt == input_file_type::cf_seg) {
            std::getline(is, sequence, '\t');  // skip '\t'
            std::getline(is, sequence);        // DNA sequence
        } else {
            static_assert(fmt == input_file_type::fasta);
            std::getline(is, sequence);  // header sequence
            std::getline(is, sequence);  // DNA sequence
        }

        if (sequence.length() < k) continue;

        /* transform 50% of the read nucleotides into lower-case letters
           (assuming the input is upper-case):
           lower-case kmers must be found anyway in the index */
        if ((num_sequences & 1) == 0) {
            std::transform(sequence.begin(), sequence.end(), sequence.begin(),
                           [](char c) { return std::tolower(c); });
        }
        ++num_sequences;

        for (uint64_t i = 0; i + k <= sequence.length(); ++i) {
            assert(util::is_valid<kmer_t>(sequence.data() + i, k));

            kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(sequence.data() + i, k);
            auto orientation = constants::forward_orientation;

            if (num_kmers != 0 and num_kmers % 5000000 == 0) {
                std::cout << "checked " << num_kmers << " kmers" << std::endl;
            }

            /* transform 50% of the kmers into their reverse complements */
            if ((num_kmers & 1) == 0) {
                uint_kmer.reverse_complement_inplace(k);
                orientation = constants::backward_orientation;
            }

            util::uint_kmer_to_string(uint_kmer, expected_kmer_str.data(), k);
            uint64_t id = dict.lookup(expected_kmer_str.c_str());

            /*
                Since we assume that we stream through the file from which the index was built,
                ids are assigned sequentially to kmers, so it must be id == num_kmers.
            */
            if (id != num_kmers) std::cout << "wrong id assigned" << std::endl;

            if (id == constants::invalid_uint64) {
                std::cout << "kmer '" << expected_kmer_str << "' not found!" << std::endl;
            }
            assert(id != constants::invalid_uint64);

            auto curr = dict.lookup_advanced(expected_kmer_str.c_str());
            assert(curr.kmer_id == id);

            if (curr.kmer_orientation != orientation) {
                std::cout << "ERROR: got orientation " << int(curr.kmer_orientation)
                          << " but expected " << int(orientation) << std::endl;
            }
            assert(curr.kmer_orientation == orientation);

            if (num_kmers == 0) {
                if (curr.contig_id != 0) {
                    std::cout << "contig_id " << curr.contig_id << " but expected 0" << std::endl;
                }
                assert(curr.contig_id == 0);  // at the beginning, contig_id must be 0
            } else {
                if (curr.kmer_id != prev.kmer_id + 1) {
                    std::cout << "ERROR: got curr.kmer_id " << curr.kmer_id << " but expected "
                              << prev.kmer_id + 1 << std::endl;
                }
                assert(curr.kmer_id == prev.kmer_id + 1);  // kmer_id must be sequential

                if (curr.kmer_id_in_contig >= curr.contig_size) {
                    std::cout << "ERROR: got curr.kmer_id_in_contig " << curr.kmer_id_in_contig
                              << " but expected something < " << curr.contig_size << std::endl;
                }
                assert(curr.kmer_id_in_contig <
                       curr.contig_size);  // kmer_id_in_contig must always be < contig_size

                if (curr.contig_id == prev.contig_id) {
                    /* same contig */
                    if (curr.contig_size != prev.contig_size) {
                        std::cout << "ERROR: got curr.contig_size " << curr.contig_size
                                  << " but expected " << prev.contig_size << std::endl;
                    }
                    assert(curr.contig_size == prev.contig_size);  // contig_size must be same
                    if (curr.kmer_id_in_contig != prev.kmer_id_in_contig + 1) {
                        std::cout << "ERROR: got curr.kmer_id_in_contig " << curr.kmer_id_in_contig
                                  << " but expected " << prev.kmer_id_in_contig + 1 << std::endl;
                    }
                    assert(curr.kmer_id_in_contig ==
                           prev.kmer_id_in_contig + 1);  // kmer_id_in_contig must be sequential
                } else {
                    /* we have changed contig */
                    if (curr.contig_id != prev.contig_id + 1) {
                        std::cout << "ERROR: got curr.contig_id " << curr.contig_id
                                  << " but expected " << prev.contig_id + 1 << std::endl;
                    }
                    assert(curr.contig_id ==
                           prev.contig_id + 1);  // contig_id must be sequential since we stream
                    if (curr.kmer_id_in_contig != 0) {
                        std::cout << "ERROR: got curr.kmer_id_in_contig " << curr.kmer_id_in_contig
                                  << " but expected 0" << std::endl;
                    }
                    assert(curr.kmer_id_in_contig ==
                           0);  // kmer_id_in_contig must be 0 when we change contig
                }
            }

            /* check also contig_size() */
            uint64_t contig_size = dict.contig_size(curr.contig_id);
            if (contig_size != curr.contig_size) {
                std::cout << "ERROR: got contig_size " << contig_size << " but expected "
                          << curr.contig_size << std::endl;
            }
            assert(contig_size == curr.contig_size);

            prev = curr;

            // check access
            dict.access(id, got_kmer_str.data());
            kmer_t got_uint_kmer = util::string_to_uint_kmer<kmer_t>(got_kmer_str.data(), k);
            kmer_t got_uint_kmer_rc = got_uint_kmer;
            got_uint_kmer_rc.reverse_complement_inplace(k);
            if (got_uint_kmer != uint_kmer and got_uint_kmer_rc != uint_kmer) {
                std::cout << "ERROR: got '" << got_kmer_str << "' but expected '"
                          << expected_kmer_str << "'" << std::endl;
            }
            ++num_kmers;
        }
    }
    std::cout << "checked " << num_kmers << " kmers" << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;

    return check_correctness_negative_lookup(dict);
}

template <class kmer_t, input_file_type fmt>
bool check_correctness_navigational_kmer_query(std::istream& is, dictionary<kmer_t> const& dict) {
    const uint64_t k = dict.k();
    std::string sequence;
    uint64_t num_kmers = 0;

    std::cout << "checking correctness of navigational queries for kmers..." << std::endl;

    while (!is.eof())  //
    {
        if constexpr (fmt == input_file_type::cf_seg) {
            std::getline(is, sequence, '\t');  // skip '\t'
            std::getline(is, sequence);        // DNA sequence
        } else {
            static_assert(fmt == input_file_type::fasta);
            std::getline(is, sequence);  // header sequence
            std::getline(is, sequence);  // DNA sequence
        }
        for (uint64_t i = 0; i + k <= sequence.length(); ++i) {
            assert(util::is_valid<kmer_t>(sequence.data() + i, k));
            if (num_kmers != 0 and num_kmers % 5000000 == 0) {
                std::cout << "checked " << num_kmers << " kmers" << std::endl;
            }
            neighbourhood<kmer_t> curr = dict.kmer_neighbours(sequence.data() + i);
            if (i + k < sequence.length()) {
                char next_nuc = sequence[i + k];
                bool next_nuc_not_found = curr.forward[kmer_t::char_to_uint(next_nuc)].kmer_id ==
                                          constants::invalid_uint64;
                if (next_nuc_not_found) {
                    std::cout << "expected forward[" << next_nuc << "]" << std::endl;
                }
                assert(!next_nuc_not_found);
            }

            if (i != 0) {
                char prev_nuc = sequence[i - 1];
                bool prev_nuc_not_found = curr.backward[kmer_t::char_to_uint(prev_nuc)].kmer_id ==
                                          constants::invalid_uint64;
                if (prev_nuc_not_found) {
                    std::cout << "expected backward[" << prev_nuc << "]" << std::endl;
                }
                assert(!prev_nuc_not_found);
            }

            ++num_kmers;
        }
    }
    std::cout << "checked " << num_kmers << " kmers" << std::endl;

    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

template <class kmer_t>
bool check_correctness_weights(std::istream& is, dictionary<kmer_t> const& dict) {
    uint64_t k = dict.k();
    std::string line;
    uint64_t kmer_id = 0;

    if (!dict.weighted()) {
        std::cerr << "ERROR: the dictionary does not store weights" << std::endl;
        return false;
    }

    std::cout << "checking correctness of weights..." << std::endl;

    while (!is.eof()) {
        std::getline(is, line);  // header line
        if (line.empty()) break;

        uint64_t i = 0;
        i = line.find_first_of(' ', i);
        assert(i != std::string::npos);

        i += 1;
        i += 5;  // skip "LN:i:"
        uint64_t j = line.find_first_of(' ', i);
        assert(j != std::string::npos);

        char* end;
        uint64_t seq_len = std::strtoull(line.data() + i, &end, 10);
        i = j + 1;
        i += 5;  // skip "ab:Z:"

        for (uint64_t j = 0; j != seq_len - k + 1; ++j, ++kmer_id) {
            uint64_t expected = std::strtoull(line.data() + i, &end, 10);
            i = line.find_first_of(' ', i) + 1;
            uint64_t got = dict.weight(kmer_id);
            if (expected != got) {
                std::cout << "ERROR for kmer_id " << kmer_id << ": expected " << expected
                          << " but got " << got << std::endl;
            }
            if (kmer_id != 0 and kmer_id % 5000000 == 0) {
                std::cout << "checked " << kmer_id << " weights" << std::endl;
            }
        }

        std::getline(is, line);  // skip DNA sequence
    }

    std::cout << "checked " << kmer_id << " weights" << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

/*
   The input file must be the one the index was built from.
   Throughout the code, we assume the input does not contain any duplicate.
*/
template <class kmer_t>
bool check_correctness_lookup_access(dictionary<kmer_t> const& dict, std::string const& filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        if (util::ends_with(filename, ".cf_seg.gz")) {
            good = check_correctness_lookup_access<kmer_t, input_file_type::cf_seg>(zis, dict);
        } else {
            good = check_correctness_lookup_access<kmer_t, input_file_type::fasta>(zis, dict);
        }
    } else {
        if (util::ends_with(filename, ".cf_seg")) {
            good = check_correctness_lookup_access<kmer_t, input_file_type::cf_seg>(is, dict);
        } else {
            good = check_correctness_lookup_access<kmer_t, input_file_type::fasta>(is, dict);
        }
    }
    is.close();
    return good;
}

/*
   The input file must be the one the index was built from.
   Throughout the code, we assume the input does not contain any duplicate.
*/
template <class kmer_t>
bool check_correctness_navigational_kmer_query(dictionary<kmer_t> const& dict,
                                               std::string const& filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        if (util::ends_with(filename, ".cf_seg.gz")) {
            good = check_correctness_navigational_kmer_query<kmer_t, input_file_type::cf_seg>(zis,
                                                                                              dict);
        } else {
            good = check_correctness_navigational_kmer_query<kmer_t, input_file_type::fasta>(zis,
                                                                                             dict);
        }
    } else {
        if (util::ends_with(filename, ".cf_seg")) {
            good = check_correctness_navigational_kmer_query<kmer_t, input_file_type::cf_seg>(is,
                                                                                              dict);
        } else {
            good =
                check_correctness_navigational_kmer_query<kmer_t, input_file_type::fasta>(is, dict);
        }
    }
    is.close();
    return good;
}

/*
   The input file must be the one the index was built from.
   Only for FASTA files since CUTTLEFISH does not store kmer weights.
*/
template <class kmer_t>
bool check_correctness_weights(dictionary<kmer_t> const& dict, std::string const& filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        good = check_correctness_weights(zis, dict);
    } else {
        good = check_correctness_weights(is, dict);
    }
    is.close();
    return good;
}

/*
    Masked-superstring input (port of check_correctness_lookup_access_masked_superstring of
    12_sshash-bigreedy-ms, via 14_sshash-v3-bigreedy-ms). K-mer ids must be the ranks of the
    uppercase positions in text order: this is what detects masks that are not min-one, which
    construction does not detect (a k-mer represented twice gets one id only). Also checked:
    lookup_advanced (contig = interval), access, and that every k-mer starting at a lowercase
    position is found if and only if it is represented (the junction negatives of the overlapping
    intervals).
*/
template <class kmer_t>
bool check_correctness_lookup_access_masked_superstring(std::istream& is,
                                                        dictionary<kmer_t> const& dict) {
    const uint64_t k = dict.k();
    std::string sequence = read_masked_superstring(is);
    auto intervals = masked_superstring_intervals(sequence, k);

    std::cout << "checking correctness of lookup and access over the masked superstring..."
              << std::endl;

    auto canonical = [&](kmer_t x) {
        kmer_t x_rc = x;
        x_rc.reverse_complement_inplace(k);
        return std::min(x, x_rc);
    };
    std::vector<kmer_t> represented;
    represented.reserve(dict.size());

    uint64_t expected_kmer_id = 0;
    std::string kmer_str(k, 0);
    std::string got_kmer_str(k, 0);
    for (uint64_t contig_id = 0; contig_id != intervals.size(); ++contig_id) {
        auto interval = intervals[contig_id];
        uint64_t contig_size = interval.end - interval.begin - k + 1;
        for (uint64_t pos = interval.begin; pos + k <= interval.end; ++pos, ++expected_kmer_id) {
            if (expected_kmer_id != 0 and expected_kmer_id % 5000000 == 0) {
                std::cout << "checked " << expected_kmer_id << " kmers" << std::endl;
            }
            kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(sequence.data() + pos, k);
            represented.push_back(canonical(uint_kmer));
            kmer_t queried = uint_kmer;
            int64_t orientation = constants::forward_orientation;
            if (expected_kmer_id & 1) {
                /* transform 50% of the queries into their reverse complements */
                queried.reverse_complement_inplace(k);
                orientation = constants::backward_orientation;
            }
            util::uint_kmer_to_string(queried, kmer_str.data(), k);

            uint64_t id = dict.lookup(kmer_str.c_str());
            if (id != expected_kmer_id) {
                std::cout << "ERROR: expected kmer_id " << expected_kmer_id << " but got " << id
                          << " (at position " << pos << ")" << std::endl;
                return false;
            }
            auto res = dict.lookup_advanced(kmer_str.c_str());
            if (res.kmer_id != expected_kmer_id or res.contig_id != contig_id or
                res.kmer_id_in_contig != pos - interval.begin or res.contig_size != contig_size or
                (!dict.canonical() and res.kmer_orientation != orientation)) {
                std::cout << "ERROR: at position " << pos << ", expected (kmer_id "
                          << expected_kmer_id << ", contig " << contig_id << ", kmer_id_in_contig "
                          << pos - interval.begin << ", contig_size " << contig_size
                          << ", orientation " << orientation << ") but got (" << res.kmer_id
                          << ", " << res.contig_id << ", " << res.kmer_id_in_contig << ", "
                          << res.contig_size << ", " << res.kmer_orientation << ")" << std::endl;
                return false;
            }
            if (dict.contig_size(contig_id) != contig_size) {
                std::cout << "ERROR: expected contig_size " << contig_size << " but got "
                          << dict.contig_size(contig_id) << std::endl;
                return false;
            }
            dict.access(id, got_kmer_str.data());
            if (!(util::string_to_uint_kmer<kmer_t>(got_kmer_str.data(), k) == uint_kmer)) {
                std::cout << "ERROR: accessed kmer '" << got_kmer_str << "' at id " << id
                          << " does not match the kmer at position " << pos << std::endl;
                return false;
            }
        }
    }
    if (expected_kmer_id != dict.size()) {
        std::cout << "ERROR: checked " << expected_kmer_id << " kmers but the dictionary has "
                  << dict.size() << std::endl;
        return false;
    }
    std::cout << "checked " << expected_kmer_id << " kmers" << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;

    std::cout << "checking correctness of lookup at the non-represented positions..." << std::endl;
    std::sort(represented.begin(), represented.end());
    uint64_t num_checked = 0, num_found = 0;
    uint64_t interval_id = 0;
    for (uint64_t pos = 0; pos + k <= sequence.size(); ++pos) {
        while (interval_id != intervals.size() and intervals[interval_id].end < pos + k) {
            ++interval_id;
        }
        if (interval_id != intervals.size() and intervals[interval_id].begin <= pos) continue;
        kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(sequence.data() + pos, k);
        bool expected =
            std::binary_search(represented.begin(), represented.end(), canonical(uint_kmer));
        util::uint_kmer_to_string(uint_kmer, kmer_str.data(), k);
        bool got = dict.lookup(kmer_str.c_str()) != constants::invalid_uint64;
        if (got != expected) {
            std::cout << "ERROR: the kmer at the non-represented position " << pos << " is "
                      << (got ? "found" : "not found") << " but it is "
                      << (expected ? "" : "not ") << "represented" << std::endl;
            return false;
        }
        ++num_checked;
        num_found += got;
    }
    std::cout << "checked " << num_checked << " kmers (" << num_found
              << " of them represented elsewhere)" << std::endl;
    std::cout << "EVERYTHING OK!" << std::endl;

    std::cout << "checking correctness of negative lookup with random kmers..." << std::endl;
    uint64_t num_lookups = std::min<uint64_t>(1000000, dict.size());
    for (uint64_t i = 0; i != num_lookups; ++i) {
        random_kmer(kmer_str.data(), k);
        kmer_t uint_kmer = util::string_to_uint_kmer<kmer_t>(kmer_str.data(), k);
        bool expected =
            std::binary_search(represented.begin(), represented.end(), canonical(uint_kmer));
        bool got = dict.lookup(kmer_str.c_str()) != constants::invalid_uint64;
        if (got != expected) {
            std::cout << "ERROR: random kmer '" << kmer_str << "' is " << (got ? "" : "not ")
                      << "found but it is " << (expected ? "" : "not ") << "represented"
                      << std::endl;
            return false;
        }
    }
    std::cout << "EVERYTHING OK!" << std::endl;
    return true;
}

template <class kmer_t>
bool check_correctness_lookup_access_masked_superstring(dictionary<kmer_t> const& dict,
                                                        std::string const& filename) {
    std::ifstream is(filename.c_str());
    if (!is.good()) throw std::runtime_error("error in opening the file '" + filename + "'");
    bool good = true;
    if (util::ends_with(filename, ".gz")) {
        zip_istream zis(is);
        good = check_correctness_lookup_access_masked_superstring<kmer_t>(zis, dict);
    } else {
        good = check_correctness_lookup_access_masked_superstring<kmer_t>(is, dict);
    }
    is.close();
    return good;
}

}  // namespace sshash