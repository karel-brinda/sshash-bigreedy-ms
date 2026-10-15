#pragma once

#include "dictionary.hpp"
#include "offsets.hpp"
#include "kmer.hpp"

namespace sshash {

/* interval_offsets subsumes a plain SPSS (degenerate intervals) and supports
   BiGreedy masked superstrings, whose strings are the maximal runs of
   represented kmer starts, embedded in one single text. */
using dictionary_type = dictionary<default_kmer_t, interval_offsets>;
// using dictionary_type = dictionary<default_kmer_t, decoded_offsets>;
// using dictionary_type = dictionary<default_kmer_t, encoded_offsets>;

}  // namespace sshash
