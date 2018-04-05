// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <RcppParallel.h>
using namespace RcppParallel;

#define POSTERIOR_MATCH_FILE "/data1/igor/himap/himap_mergepairs_match_qs.txt"
#define POSTERIOR_MISMATCH_FILE "/data1/igor/himap/himap_mergepairs_mismatch_qs.txt"
#define PHRED_OFFSET 33

// Precomputer posterior q scores for matching and mismatching bases
std::vector< std::vector<int> > post_q_match;
std::vector< std::vector<int> > post_q_mismatch;

// Load pre-computed posterior quality scores for matching and mismatching
// bases.
std::vector< std::vector<int> > load_posterior (std::string filename) {
  // Load pre-computer posterior Q scores into a vector of vector of ints.
  std::vector< std::vector<int> > out;
  std::vector<int> row;
  std::ifstream fin(filename.c_str());
  int prev_q1 = 0, q1, q2, q;
  while (fin >> q1 >> q2 >> q) {
    if (q1 == prev_q1) {
      // same first index
      row.push_back(q);
    } else {
      // new first index
      out.push_back(row);
      row.clear();
      row.push_back(q);
    }
    prev_q1 = q1;
  }
  return(out);
}


// Modified dada2 nwalign_endsfree to allow Ns and 
void nt2int(char *oseq, const char *iseq) {
  int i, len = strlen(iseq);
  for (i = 0; i < len; i++, iseq++, oseq++) {
    switch (*iseq) {
    case 'A':
      *oseq = 1;
      break;
    case 'C':
      *oseq = 2;
      break;
    case 'G':
      *oseq = 3;
      break;
    case 'T':
      *oseq = 4;
      break;
    case 'N':
      *oseq = 5;
      break;
    case '-':
      *oseq = '-';
      break;
    default:
      //Rprintf("invalid character in input:%c.\n",*iseq);
      *oseq = '\0';
    }
  }
  *oseq = '\0';
  return;
}

void int2nt(char *oseq, const char *iseq) {
  int i, len = strlen(iseq);
  for (i = 0; i < len; i++, iseq++, oseq++) {
    switch (*iseq) {
    case 1:
      *oseq = 'A';
      break;
    case 2:
      *oseq = 'C';
      break;
    case 3:
      *oseq = 'G';
      break;
    case 4:
      *oseq = 'T';
      break;
    case 5:
      *oseq = 'N';
      break;
    case '-':
      *oseq = '-';
      break;
    default:
      break;
    }
  }
  *oseq = '\0';
  return;
}

char **himap_nwalign_endsfree(const char *s1, const char *s2, 
                              const char *q1, const char *q2,
                              int score[5][5], int gap_p) {
  static size_t nnw = 0;
  int i, j;
  int l, r;
  int band = -1;
  unsigned int len1 = strlen(s1);
  unsigned int len2 = strlen(s2);
  int diag, left, up;
  
  unsigned int nrow = len1+1;
  unsigned int ncol = len2+1;
  int *d = (int *) malloc(nrow * ncol * sizeof(int)); //E
  int *p = (int *) malloc(nrow * ncol * sizeof(int)); //E
  if(d == NULL || p == NULL) Rcpp::stop("Memory allocation failed.");
  
  // Fill out left columns of d, p.
  for (i = 0; i <= len1; i++) {
    d[i*ncol] = 0; // ends-free gap
    p[i*ncol] = 3;
  }
  
  // Fill out top rows of d, p.
  for (j = 0; j <= len2; j++) {
    d[j] = 0; // ends-free gap
    p[j] = 2;
  }
  
  // Calculate left/right-bands in case of different lengths
  int lband, rband;
  if(len2 > len1) {
    lband = band;
    rband = band+len2-len1;
  } else if(len1 > len2) {
    lband = band+len1-len2;
    rband = band;
  } else {
    lband = band;
    rband = band;
  }
  
  // Fill out band boundaries of d.
  if(band>=0 && (band<len1 || band<len2)) {
    for(i=0;i<=len1;i++) {
      if(i-lband-1 >= 0) { d[i*ncol + i-lband-1] = -9999; }
      if(i+rband+1 <= len2) { d[i*ncol + i+rband+1] = -9999; }
    }
  }
  
  // Fill out the body of the DP matrix.
  for (i = 1; i <= len1; i++) {
    if(band>=0) {
      l = i-lband; if(l < 1) { l = 1; }
      r = i+rband; if(r>len2) { r = len2; }
    } else { l=1; r=len2; }
    
    for (j = l; j <= r; j++) {
      // Score for the left move.
      if (i == len1) {
        left = d[i*ncol + j-1]; // Ends-free gap.
      } else {
        left = d[i*ncol + j-1] + gap_p;
      }
      
      // Score for the up move.
      if (j == len2) {
        up = d[(i-1)*ncol + j]; // Ends-free gap.
      } else {
        up = d[(i-1)*ncol + j] + gap_p;
      }
      
      // Score for the diagonal move.
      diag = d[(i-1)*ncol + j-1] + score[s1[i-1]-1][s2[j-1]-1];
      
      // Break ties and fill in d,p.
      if (up >= diag && up >= left) {
        d[i*ncol + j] = up;
        p[i*ncol + j] = 3;
      } else if (left >= diag) {
        d[i*ncol + j] = left;
        p[i*ncol + j] = 2;
      } else {
        d[i*ncol + j] = diag;
        p[i*ncol + j] = 1;
      }
    }
  }
  
  char *al0 = (char *) malloc((len1+len2+1) * sizeof(char));
  char *al1 = (char *) malloc((len1+len2+1) * sizeof(char));
  if(al0 == NULL || al1 == NULL) Rcpp::stop("Memory allocation failed.");
  
  char *qs0 = (char *) malloc((len1+len2+1) * sizeof(char));
  char *qs1 = (char *) malloc((len1+len2+1) * sizeof(char));
  if(qs0 == NULL || qs1 == NULL) Rcpp::stop("Memory allocation failed.");
  
  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;  
  
  while ( i > 0 || j > 0 ) {
    //    Rprintf("(%i, %i): p=%i, d=%i\n", i, j, p[i*ncol + j], d[i*ncol + j]);
    switch ( p[i*ncol + j] ) {
    case 1:
      i = i - 1;
      j = j - 1;
      al0[len_al] = s1[i];
      al1[len_al] = s2[j];
      qs0[len_al] = q1[i];
      qs1[len_al] = q2[j];
      break;
    case 2:
      j = j - 1;
      al0[len_al] = '-';
      al1[len_al] = s2[j];
      qs0[len_al] = ' ';
      qs1[len_al] = q2[j];
      break;
    case 3:
      i = i - 1;
      al0[len_al] = s1[i];
      al1[len_al] = '-';
      qs0[len_al] = q1[i];
      qs1[len_al] = ' ';
      break;
    default:
      Rcpp::stop("N-W Align out of range.");
    }
    len_al++;
  }
  al0[len_al] = '\0';
  al1[len_al] = '\0';
  qs0[len_al] = '\0';
  qs1[len_al] = '\0';
  
  // Now merge the two alignments
  for (i=0; i<len_al; i++) {
    // "letter" (int) of the first alignment string
    if (al0[len_al-i-1] == al1[len_al-i-1]) {
      // matching letters in the alignment
    }
  }
  
  // Allocate memory to alignment strings.
  char **al = (char **) malloc( 4 * sizeof(char *) ); //E
  if (al == NULL)  Rcpp::stop("Memory allocation failed.");
  al[0] = (char *) malloc(len_al+1); //E
  al[1] = (char *) malloc(len_al+1); //E
  al[2] = (char *) malloc(len_al+1); //E
  al[3] = (char *) malloc(len_al+1); //E
  char *merged_seq = (char *) malloc(len_al+1);
  char *merged_qual = (char *) malloc(len_al+1);
  if (al[0] == NULL || al[1] == NULL)  Rcpp::stop("Memory allocation failed.");
  
  // Reverse the alignment strings (since traced backwards).
  for (i=0;i<len_al;i++) {
    al[0][i] = al0[len_al-i-1];
    al[1][i] = al1[len_al-i-1];
    al[2][i] = qs0[len_al-i-1];
    al[3][i] = qs1[len_al-i-1];
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';
  al[2][len_al] = '\0';
  al[3][len_al] = '\0';
  
  // Free allocated memory
  free(d);
  free(p);
  free(al0);
  free(al1);
  free(qs0);
  free(qs1);
  
  nnw++;
  return al;
}

std::vector<std::string> *himap_merge_alignment(char** al) {
  // al[0] seq1, al[1] seq2, al[2] qual1, al[3] qual 2
  int i, len = strlen(al[0]);
  int q1, q2, q;
  // char** out = (char **) malloc(2*sizeof(char *));
  // out[0] = (char *) malloc(len+1); // merged sequence
  // out[1] = (char *) malloc(len+1); // merged qual scores
  std::vector<std::string> *out;
  
  for (i=0; i<len; i++) {
    q1 = int(al[2][i]) - PHRED_OFFSET;
    q2 = int(al[3][i]) - PHRED_OFFSET;
    if (al[0][i] == al[1][i]) { // matching letters
      if (al[0][i] == 'N') {
        out[0][i] = 'N';
        out[1][i] = '\"';
      } else {
        out[0][i] = al[0][i];
        q = post_q_match[q1][q2] + PHRED_OFFSET;
        out[1][i] = (char)q;
      }
    } else { // mismatch
      if (al[0][i] == '-') { // fwd is missing, use just reverse 
        out[0][i] = al[1][i];
        out[1][i] = al[3][i];
      } else if (al[1][i] == '-') { // rev is missing, use forward
        out[0][i] = al[0][i];
        out[1][i] = al[2][i];
      } else { // mismatch and not indel
        q = post_q_mismatch[q1][q2] + PHRED_OFFSET;
        out[1][i] = (char)q;
        if (q1 >= q2) {
          out[0][i] = al[0][i];
        } else {
          out[0][i] = al[1][i];
        }
      }
    }
  }
  return(out);
}

// [[Rcpp::export]]
Rcpp::CharacterVector C_nwalign(std::string s1, std::string s2, 
                                std::string q1, std::string q2,
                                int match, int mismatch, int gap_p) {
  // q1 and q2 are quality scores
  // Load pre-computed posterior quality scores
  post_q_match = load_posterior(POSTERIOR_MATCH_FILE);
  post_q_mismatch = load_posterior(POSTERIOR_MISMATCH_FILE);
  
  // Make  c-style 2d score array
  int i, j;
  int c_score[5][5];
  for(i=0;i<5;i++) {
    for(j=0;j<5;j++) {
      if(i==j) { 
        if (i==4) {
          // both i and j are N, declare it mismatch
          c_score[i][j] = mismatch;
        } else {
          c_score[i][j] = match; 
        }
      }
      else { 
        c_score[i][j] = mismatch; 
      }
    }
  }
  
  // Make integer-ized c-style sequence strings
  char **al;
  std::vector<std::string> merged;
  char *seq1 = (char *) malloc(s1.size()+1); //E
  char *seq2 = (char *) malloc(s2.size()+1); //E
  if (seq1 == NULL || seq2 == NULL)  Rcpp::stop("Memory allocation failed.");
  nt2int(seq1, s1.c_str());
  nt2int(seq2, s2.c_str());
  const char *qs1 = q1.c_str();
  const char *qs2 = q2.c_str();
  
  // Perform alignment and convert back to ACGT
  al = himap_nwalign_endsfree(seq1, seq2, qs1, qs2, c_score, gap_p);
  int2nt(al[0], al[0]);
  int2nt(al[1], al[1]);
  
  // Now merge the aligned sequences
  merged = himap_merge_alignment(al);
  
  // Generate R-style return vector
  Rcpp::CharacterVector rval;
  rval.push_back(std::string(merged[0]));
  rval.push_back(std::string(merged[1]));
  // rval.push_back(std::string(al[0]));
  // rval.push_back(std::string(al[1]));
  // rval.push_back(std::string(al[2]));
  // rval.push_back(std::string(al[3]));
  // Clean up
  free(seq1);
  free(seq2);
  free(al[0]);
  free(al[1]);
  free(al[2]);
  free(al[3]);
  free(merged[0]);
  free(merged[1]);
  free(al);
  return(rval);
}

struct AlignMerge : public RcppParallel::Worker {
  // input sequences and quality scores
  std::string s1, s2, q1, q2;
  // output merged sequence and quality scores are in an array
  // merged[0] = sequence, merged[1] = quality scores
  char **merged = (char **) malloc(2*sizeof(char *));
  
  // Initialize an object
  AlignMerge(std::string s1, std::string s2, std::string q1, std::string q2,
             char **merged) : 
    s1(s1), s2(s2), q1(q1), q2(q2), merged(merged) {}
};


struct SquareRoot : public RcppParallel::Worker
{
  // source matrix
  RMatrix<double> input;
  
  // destination matrix
  RMatrix<double> output;
  
  // initialize with source and destination
  SquareRoot(const Rcpp::NumericMatrix input, Rcpp::NumericMatrix output) 
    : input(input), output(output) {}
  
  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    std::transform(input.begin() + begin, 
                   input.begin() + end, 
                   output.begin() + begin, 
                   ::sqrt);
  }
};


// [[Rcpp::export]]
Rcpp::CharacterVector C_parallel_mergepairs(std::vector<std::string> s1_v,
  std::vector<std::string> s2_v,
  std::vector<std::string> q1_v,
  std::vector<std::string> q2_v,
  int match, int mismatch, int gap_p) {

  // Load pre-computed posterior quality scores
  post_q_match = load_posterior(POSTERIOR_MATCH_FILE);
  post_q_mismatch = load_posterior(POSTERIOR_MISMATCH_FILE);
  
  // Make  c-style 2d score array
  int i, j;
  int c_score[5][5];
  for (i=0; i<5; i++) {
    for (j=0; j<5; j++) {
      if (i==j) { 
        if (i==4) {
          // both i and j are N, declare it mismatch
          c_score[i][j] = mismatch;
        } else {
          c_score[i][j] = match; 
        }
      }
      else { 
        c_score[i][j] = mismatch; 
      }
    }
  }
  
}


// [[Rcpp::export]]
Rcpp::CharacterVector C_mergepairs(
    std::vector<std::string> s1_v, std::vector<std::string> s2_v, 
    std::vector<std::string> q1_v, std::vector<std::string> q2_v,
    int match, int mismatch, int gap_p) {
  // q1 and q2 are quality scores
  // Make  c-style 2d score array
  int i, j;
  int c_score[5][5];
  for(i=0;i<5;i++) {
    for(j=0;j<5;j++) {
      if(i==j) { 
        if (i==4) {
          // both i and j are N, declare it mismatch
          c_score[i][j] = mismatch;
        } else {
          c_score[i][j] = match; 
        }
      }
      else { 
        c_score[i][j] = mismatch; 
      }
    }
  }
  
  // Make integer-ized c-style sequence strings
  std::vector<char **> al;
  std::vector<char *> seq1, seq2, qs1, qs2;
  std::vector<std::string> :: iterator it;
  
  for (it=s1_v.begin(); it != s1_v.end(); ++it) {
    seq1.push_back((char *) malloc((*it).size()+1));
  }
  for (it=s2_v.begin(); it != s2_v.end(); ++it) {
    seq2.push_back((char *) malloc((*it).size()+1));
  }
  for (it=q1_v.begin(); it != q1_v.end(); ++it) {
    qs1.push_back((char *) malloc((*it).size()+1));
  }
  for (it=q2_v.begin(); it != q2_v.end(); ++it) {
    qs2.push_back((char *) malloc((*it).size()+1));
  }
  
  // if (seq1 == NULL || seq2 == NULL)  Rcpp::stop("Memory allocation failed.");
  int k = 0;
  for (it=s1_v.begin(); it != s1_v.end(); ++it) {
    nt2int(seq1.at(k++), (*it).c_str());
  }
  k = 0;
  for (it=s2_v.begin(); it != s2_v.end(); ++it) {
    nt2int(seq2.at(k++), (*it).c_str());
  }
  
  // Now do multithreaded alignment
  
  // Perform alignment and convert back to ACGT
  // al = himap_nwalign_endsfree(seq1, seq2, qs1, qs2, c_score, gap_p);
  // int2nt(al[0], al[0]);
  // int2nt(al[1], al[1]);
  // Generate R-style return vector
  Rcpp::CharacterVector rval;
  // rval.push_back(std::string(al[0]));
  // rval.push_back(std::string(al[1]));
  // // Clean up
  // free(seq1);
  // free(seq2);
  // free(al[0]);
  // free(al[1]);
  // free(al);
  return(rval);
}

