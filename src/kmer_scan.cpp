/* k-mer scanner
 * uses only A,C,G,T
 * ref. Bioinformatics algorithms I p.39-42
 */
#include <Rcpp.h>
#include <map>
#include <cmath>
#include <sstream>
#include <fstream>
#include <boost/algorithm/string.hpp>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]

static std::map<char, int> letters = {{'A',0},{'C',1}, {'G',2}, {'T',3}, {'N', 4} };
static std::map<int, char> lettersInverse = {{ 0,'A'},{ 1,'C'}, { 2,'G'}, { 3,'T'}, {4, 'N'} };

int symbol2num (const char sym) {
  switch (sym) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  case 'N':
    return 4;
  default:
    Rprintf("symbol2num: invalid character in input:%c.\n", sym);
    return 5;
  }
}

char num2symbol (int n) {
  switch (n) {
  case 1:
    return 'A';
  case 2:
    return 'C';
  case 3:
    return 'G';
  case 4:
    return 'T';
  case 5:
    return 'N';
  default:
    Rprintf("num2symbol: invalid number in input.");
    return '\0';
  }
}

/* Conver char sequence to integer */
int* seq_to_int(const char* in_seq) {
  int i, len = sizeof(in_seq)/sizeof(char);
  int* out_seq = new int[len];
  for (i=0; i<len; i++, in_seq++) {
    switch (*in_seq) {
    case 'A':
      out_seq[i] = 0;
      break;
    case 'C':
      out_seq[i] = 1;
      break;
    case 'G':
      out_seq[i] = 2;
      break;
    case 'T':
      out_seq[i] = 3;
      break;
    case 'N':
      out_seq[i] = 4;
      break;
    }
  }
  return out_seq;
}


int kmer2num1 (std::string kmer) {
  if (kmer.empty()) return 0;
  int symnum = letters[kmer.back()];
  kmer.pop_back();
  return 4 * kmer2num1(kmer) + symnum;
}

int kmer2num (std::string kmer) {
  // Convert kmer to number
  if (kmer.empty()) return 0;
  int symnum = symbol2num(kmer.back());
  kmer.pop_back();
  return 4*kmer2num(kmer) + symnum;
}

// Convert a single kmer to number.
int kmer_to_num (int* kmer, size_t first, size_t last) {
  if (last == first) return 0;
  return 4*kmer_to_num(kmer, first, last-1) + kmer[last];
}


std::string num2kmer (int n, int k) {
  // Convert number n to kmer
  if (k == 1) return std::string(1, num2symbol(n));
  return num2kmer(n/4, k-1) +  num2symbol(n%4);
}

// [[Rcpp::export]]
std::set<int> kmer_extract (std::string seq, int k) {
  /*
   * If N is found, set n=k, then n-- each round until n=0.
   * n tells how many more iterations we need to skip.
   */
  int k_n;
  std::set<int> kmers;
  // Copy the entire string to C char*, so we don't make
  // copies of each substring in the for loop. 
  int len = seq.length();
  const char* seqc = seq.c_str();
  int* seqn = seq_to_int(seqc);
  
  for (int i=0; i<len; i++) {
    std::cout << i << ": " << seqn[i] << "\n";
  }
  // Locations of the first and last kmer element
  size_t i0=0;
  int n=0;
  // Check if there's N in the first kmer. If yes, start from
  // the first position after N.
  for (size_t i=0; i<k; ++i) {
     if (seqn[i] == 4) n++;
  }
  for (size_t i=0; i<len-k+1; ++i) {
    if (n==0) { // this window does not contains Ns, generate kmer...
      k_n = kmer_to_num(seqn, i, i+k-1);
      kmers.insert(k_n);
    } else { // still have Ns, keep skipping
      n--;
    }
  }
  delete[] seqn;
  return kmers;
}

// std::map<int, int> kmers_from_string (std::string seq, int k) {
//   // counts is a map that maps an arbitrary kmer index (which is int)
//   // to the number of times we saw it (int)
//   std::map<int, int> counts;
//   int len = seq.size(); // seq length
//   int k_n;  // numerical representation of kmer in the loop
//   // map.find(key) != map.end()
//   for (size_t i=0; i<len-k+1; ++i) {
//     k_n = kmer2num(seq.substr(i, k));
//     // counts[k_n] = counts.find(k_n) == counts.end() ? 1 : counts[k_n]++;
//     if (counts.find(k_n) == counts.end()) { // key not found
//       counts[k_n] = 1;
//     } else { // key found
//       counts[k_n]++;
//     }
//   }
//   return counts;
// }

// std::vector<std::map<int, int> > kmers_from_strings (std::string seq, int k) {
//   std::vector<std::map<int, int> > counts;
//   
// }



int *kmer_counts (std::string seq, int k) {
  // Count the kmers in seq
  int len = seq.size();
  // Generate empty kmer count array
  int n = pow(4, k);
  int k_n = 0;
  int *counts = new int[n];
  for (int i=0; i<n; i++) counts[i] = 0;
  
  for (size_t i=0; i<len-k+1; ++i) {
    k_n = kmer2num(seq.substr(i, k));
    counts[k_n]++;
  }
  return counts;
}

std::vector<unsigned int> all_max (unsigned int *array) {
  // Return the array containing the indices of all the maximum values in 'array'
  std::vector<unsigned int> ids = {0};
  unsigned int val = array[0];
  for (unsigned int i=1; i < sizeof(array)/sizeof(array[0]); i++) {
    if (array[i] > val) {
      ids = {i};
      val = array[i];
    } else if (array[i] == val) {
      ids.push_back(i);
    }
  }
  return ids;
}

// [[Rcpp::export]]
void test_map (std::string kmer)
{
  int i = kmer2num1(kmer);
  // std::cout << kmer << ": " << kmer2num(kmer) << "\n";
}

// [[Rcpp::export]]
void test_switch (std::string kmer)
{
  int i = kmer2num(kmer);
  //std::cout << kmer << ": " << kmer2num2(kmer) << "\n";
}

// [[Rcpp::export]]
void test2() {
  /*  Here is defined an array of 10  int size elements:                         */
  
  int array[10] = { 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
  int last = 0;
  
  /*  The array elements are numbered (indexed) 0 thru 9.
   ie: the 1st element is index 0 , the 2nd is index 1 etc.
   
   So the first element will always be at index 0.
   To find the index of the last element we use the "sizeof" operator to 
   find both the size of the whole array  AND  the size of one element:
   */
  
  printf("\nsizeof(array): %lu \nsizeof(array[0]): %lu \nsizeof(int): %lu \n", 
         sizeof(array), sizeof(array[0]), sizeof(int) );
  
  /* So if the size is 10 then the last index is 9. */
  
  last = ( sizeof(array) / sizeof(array[0]) ) - 1;
  
  printf("The last element is at index %d and the value stored there is %d. \n\n", last, array[last] );
  
}

// [[Rcpp::export]]
std::string read_file (std::string filename) {
  // Load filename into string stream. String is obtained with buffer.str()
  // https://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
  std::ifstream t(filename);
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

// [[Rcpp::export]]
std::vector<std::string> split_string (std::string text, std::string delim) {
  std::vector<std::string> strs;
  boost::split(strs, text, boost::is_any_of(delim));
  return strs;
}

// [[Rcpp::export]]
std::vector<std::string> read_file_to_vector (std::string filename, std::string delim = "\n") {
  // Load filename into a vector of strings
  // https://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
  std::ifstream t(filename);
  std::stringstream buffer;
  buffer << t.rdbuf();
  std::vector<std::string> strs;
  std::string text = buffer.str();
  boost::split(strs, text, boost::is_any_of(delim));
  return strs;
}


// Prototyping
/*** R
# x = "igor segota"
kmer  = 'AGTACTGCATCGTCGCAAATGCC'
kmer1 = 'AAAAC'
kmer2 = 'AGTTCTGCATCGACGCAAATACC'

# source('/data1/igor/himap/r/read_fastx.R')
# query_list = fasta_reader(query_fasta_filename)
# system.time(sapply(query_list[['seqs']], test_map))
# system.time(sapply(query_list[['seqs']], test_switch))
# fasta_data = fasta_reader('/data1/igor/himap/data/V3-V4_337F-805R_hang25_sequences.fasta')
*/
