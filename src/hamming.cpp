#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// Nucleotide-to-integer conversion map
static std::map<char, int> nt2int = {{'A',0}, {'C', 1}, {'G', 2}, {'T', 3}, {'N', 4}, {'-', 5}};
static std::map<char, int> nt2intext = {
  {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}, {'N', 4}, {'-', 5},
  {'Y', 6}, {'R', 7}, {'W', 8}, {'S', 9}, {'K', 10}, {'M', 11}, {'D', 12},
  {'V', 13}, {'H', 14}, {'B', 15}, {'X', 16}
};

// [[Rcpp::export]]
int hamming (std::string s1, std::string s2) {
  // Calculate hamming distance between two sequences of equal length
  int s1_i, s2_i;
  int mismatch = 0;
  size_t len = s1.length();
  // Length check
  if (len != s2.length()) {
    Rprintf("Error: aligned strings are not the same length.\n");
    return 0;
  }
  // Main loop with comparisons
  for (size_t i=0; i<len; i++) {
    s1_i = nt2int[s1.at(i)];
    s2_i = nt2int[s2.at(i)];
    if (s1_i != s2_i) mismatch++;
  }
  return mismatch;
}

// [[Rcpp::export]]
IntegerVector compare_alignment (std::string s1, std::string s2) {
  // alignment strings should not have - in both sequences (like from multiple aln)
  // in the same position.
  int s1_i, s2_i, s1_p, s2_p;
  int match = 0;
  int mismatch = 0;
  int gapopen = 0;
  int gapextend = 0;
  size_t len = s1.length();
  // Length check
  if (len != s2.length()) {
    Rprintf("Error: aligned strings are not the same length.\n");
    return R_NilValue;
  }
  // Keep track of previous symbols for gap extension.
  // Initialize them to any other than 5.
  s1_p = 0;
  s2_p = 0;
  // Main loop with comparisons
  for (size_t i=0; i<len; i++) {
    s1_i = nt2int[s1.at(i)];
    s2_i = nt2int[s2.at(i)];
    if ((s1_i == s2_i) && (s1_i != 4)) { // match except N-N
      match++;
    } else if (s1_i == 5) { // letter is -
      if (s1_p == 5) gapextend++;
      else gapopen++;
    } else if (s2_i == 5) {
      if (s2_p == 5) gapextend++;
      else gapopen++;
    } else { // the only thing left is a mismatch
      mismatch++;
    }
    s1_p = s1_i;
    s2_p = s2_i;
  }
  return IntegerVector::create(match, mismatch, gapopen, gapextend);
}

// [[Rcpp::export]]
IntegerVector compare_alignment_ext (std::string s1, std::string s2) {
  // alignment strings should not have - in both sequences (like from multiple aln)
  // in the same position.
  int s1_i, s2_i, s1_p, s2_p;
  int match = 0;
  int mismatch = 0;
  int gapopen = 0;
  int gapextend = 0;
  size_t len = s1.length();
  // Length check
  if (len != s2.length()) {
    Rprintf("Error: aligned strings are not the same length.\n");
    return R_NilValue;
  }
  // Keep track of previous symbols for gap extension.
  // Initialize them to any other than 5.
  s1_p = 0;
  s2_p = 0;
  // Main loop with comparisons
  for (size_t i=0; i<len; i++) {
    s1_i = nt2int[s1.at(i)];
    s2_i = nt2int[s2.at(i)];
    if (s1_i == 5) { // letter is -
      if (s1_p == 5) gapextend++;
      else gapopen++;
    } else if (s2_i == 5) { // letter is -
      if (s2_p == 5) gapextend++;
      else gapopen++;
    } else if (
        ( (s1_i == s2_i) ) ||
        (  s1_i == 4) || (s2_i == 4) || (s1_i == 16) || (s2_i == 16) ||
        ( ((s1_i == 6) && ((s2_i == 1) || (s2_i == 3))) ||
          ((s2_i == 6) && ((s1_i == 1) || (s1_i == 3))) ) ||
        ( ((s1_i == 7) && ((s2_i == 0) || (s2_i == 2))) ||
          ((s2_i == 7) && ((s1_i == 0) || (s1_i == 2))) ) ||
        ( ((s1_i == 8) && ((s2_i == 0) || (s2_i == 3))) ||
          ((s2_i == 8) && ((s1_i == 0) || (s1_i == 3))) ) ||
        ( ((s1_i == 9) && ((s2_i == 1) || (s2_i == 2))) ||
          ((s2_i == 9) && ((s1_i == 1) || (s1_i == 2))) ) ||
        ( ((s1_i == 10) && ((s2_i == 2) || (s2_i == 3))) ||
          ((s2_i == 10) && ((s1_i == 2) || (s1_i == 3))) ) ||
        ( ((s1_i == 11) && ((s2_i == 0) || (s2_i == 1))) ||
          ((s2_i == 11) && ((s1_i == 0) || (s1_i == 1))) ) ||
        ( ((s1_i == 12) && ((s2_i == 0) || (s2_i == 2) || (s2_i == 3) || (s2_i == 7) || (s2_i == 8) || (s2_i == 10))) ||
          ((s2_i == 12) && ((s1_i == 0) || (s1_i == 2) || (s1_i == 3) || (s1_i == 7) || (s1_i == 8) || (s1_i == 10))) ) ||
        ( ((s1_i == 13) && ((s2_i == 0) || (s2_i == 1) || (s2_i == 2) || (s2_i == 1) || (s2_i == 7) || (s2_i == 9))) ||
          ((s2_i == 13) && ((s1_i == 0) || (s1_i == 1) || (s1_i == 2) || (s1_i == 1) || (s1_i == 7) || (s1_i == 9))) ) ||
        ( ((s1_i == 14) && ((s2_i == 0) || (s2_i == 1) || (s2_i == 3) || (s2_i == 11) || (s2_i == 8) || (s2_i == 6))) ||
          ((s2_i == 14) && ((s1_i == 0) || (s1_i == 1) || (s1_i == 3) || (s1_i == 11) || (s1_i == 8) || (s1_i == 6))) ) ||
        ( ((s1_i == 15) && ((s2_i == 1) || (s2_i == 2) || (s2_i == 3) || (s2_i == 9) || (s2_i == 6) || (s2_i == 10))) ||
          ((s2_i == 15) && ((s1_i == 1) || (s1_i == 2) || (s1_i == 3) || (s1_i == 9) || (s1_i == 6) || (s1_i == 10)))
        )
    ) {
      // if symbols match or they match N, call it a match
      match++;
    } else { // the only thing left is a mismatch
      mismatch++;
    }
    s1_p = s1_i;
    s2_p = s2_i;
  }
  return IntegerVector::create(match, mismatch, gapopen, gapextend);
}
