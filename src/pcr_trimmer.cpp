#include <Rcpp.h>
#include "himap.h"
using namespace std;
// [[Rcpp::plugins(cpp11)]]

// Nucleotide-to-integer conversion map
// static map<char, int> nt2int = {{'A',0}, {'C', 1}, {'G', 2}, {'T', 3}, {'N', 4}, {'-', 5}, {' ', 5}};
// static map<int, char> int2nt = {{0, 'A'}, {1, 'C'}, {2, 'G'}, {3, 'T'}, {4, 'N'}, {5, '-'}};


int max_index(int *x, int size) {
  int id = 0;
  for (int i=1; i<size; i++) {
    if (x[i] > x[id]) id = i;
  }
  return id;
}


/* nt2int takes an input string (iseq) and converts it to numeric indices,
which are stored in the char output array (oseq).
Need memory at oseq equal to strlen of iseq, which is not checked currently.
nt2int(seq,seq) gives previous functionality.
*/
//
// void nt2int(char *oseq, const char *iseq) {
//   int i, len = strlen(iseq);
//
//   for (i = 0; i < len; i++, iseq++, oseq++) {
//     switch (*iseq) {
//     case 'A':
//       *oseq = 1;
//       break;
//     case 'C':
//       *oseq = 2;
//       break;
//     case 'G':
//       *oseq = 3;
//       break;
//     case 'T':
//       *oseq = 4;
//       break;
//     case 'N':
//       *oseq = 5;
//       break;
//     case '-':
//       *oseq = '-';
//       break;
//     default:
//       Rprintf("invalid character in input:%c.\n",*iseq);
//     *oseq = '\0';
//     }
//   }
//   *oseq = '\0';
//   return;
// }
//
// /* Converts seq in index form back to nts. */
// void int2nt(char *oseq, const char *iseq) {
//   int i, len = strlen(iseq);
//   for (i = 0; i < len; i++, iseq++, oseq++) {
//     switch (*iseq) {
//     case 1:
//       *oseq = 'A';
//       break;
//     case 2:
//       *oseq = 'C';
//       break;
//     case 3:
//       *oseq = 'G';
//       break;
//     case 4:
//       *oseq = 'T';
//       break;
//     case 5:
//       *oseq = 'N';
//       break;
//     case '-':
//       *oseq = '-';
//       break;
//     default:
//       break;
//     }
//   }
//   *oseq = '\0';
//   return;
// }



// Faster fitting alignment with a single indel penalty
// [[Rcpp::export]]
Rcpp::CharacterVector fitting_alignment_fast(std::string s1s, std::string s2s,
                                   int match=1, int mismatch=-1, int indel=-1) {

  // Output variables
  Rcpp::CharacterVector out;
  std::string s1_align;
  std::string s2_align;
  // FittingAlignment fitaln;
  // Conver std::string to char*
  char *s1 = (char *) malloc(s1s.size()+1);
  char *s2 = (char *) malloc(s2s.size()+1);

  // Initialize scoring matrix and backtracking matrix
  //int min_score = -500;

  nt2int(s1, s1s.c_str());
  nt2int(s2, s2s.c_str());

  // First index 0,1,2 middle, lower, upper scoring layers
  int len1 = strlen(s1);
  int len2 = strlen(s2);
  int nrow = len1+1;
  int ncol = len2+1;
  int hw = nrow*ncol;
  int i, j, k;
  // ri, ci;
  // i_next;
  // unsigned int s_start = 0;
  // unsigned int s_end = ncol-2;

  // Alignment scoring matrix
  int c_score[5][5];
  for (i=0; i<5; i++) {
    for (j=0; j<5; j++) {
      if (i==j) { c_score[i][j] = match; }
      else { c_score[i][j] = mismatch; }
    }
  }

  // Keep track of maximum scores
  int right, diag, down;
  // int score_last[ncol+1];

  //int max_dir[3]; // store cummulative score for each directional step in middle layer
  // int max_dir_u[2], max_dir_l[2]; // scores for lower and upper layers
  // int max_dir_last[ncol+1]; // last node in the middle layer is special for this aln

  // Convert std::string into character array
  //int diag; // diagonal score. mismatch or match
  // int max_direction;

  // Initialize scoring matrix (flattened to 1D array)
  int *s = (int *) malloc(hw*sizeof(int));
  s[0] = 0; // first element

  // Initialize Backtracking matrices
  // int *bt = new int[hw];
  // int *bt = (int *) malloc(hw*sizeof(int));
  int *b = (int *) malloc(hw*sizeof(int));
  // bt[0] = -1;
  b[0] = -1;
  // backtrack: 1=left, 2=diag, 3=up
  // special for first and last rows:
  //  if x > 3, then its a shortcut: x-3 tells us how many steps left we're going.

  /* ---------- ROW 0 (w/o first element) -------- */
  for (j=1; j<ncol; j++) {
    s[j] = 0; // zero-cost for any gap from beginning
    // bt[j] = j-1;
    b[j] = 1;
  }

  /* --------- COLUMN 0 (w/o first element) --------- */
  for (i=1; i<nrow; i++) {
    s[i*ncol] = i*indel;
    // bt[i*ncol] = (i-1)*ncol;
    b[i*ncol] = 3;
  }

  for (i=1; i<nrow; i++) {
    for (j=1; j<ncol; j++) {
      k = i*ncol + j;
      right = s[k-1] + (i==len1 ? 0 : indel); // 0 weights in last row
      // diag = s[i-ncol-1] + (s1[r-1] == s2[c-1] ? match : mismatch);
      diag = s[k-ncol-1] + c_score[s1[i-1]][s2[j-1]];
      down = s[k-ncol] + indel;
      if (down >= diag && down >= right) { // down wins
        // bt[k] = k-ncol;
        b[k] = 3;
        s[k] = down;
      } else if (right >= diag) {
        // bt[k] = k-1;
        b[k] = 1;
        s[k] = right;
      } else {
        // bt[k] = k-ncol-1;
        b[k] = 2;
        s[k] = diag;
      }
    }
  }

  char *al0 = (char *) malloc((len1+len2+1) * sizeof(char));
  char *al1 = (char *) malloc((len1+len2+1) * sizeof(char));
  if(al0 == NULL || al1 == NULL) Rcpp::stop("Memory allocation failed.");

  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;

    while (i > 0 || j > 0) {
    switch ( b[i*ncol + j] ) {
    case 2:
      al0[len_al] = s1[--i];
      al1[len_al] = s2[--j];
      break;
    case 1:
      al0[len_al] = '-';
      al1[len_al] = s2[--j];
      break;
    case 3:
      al0[len_al] = s1[--i];
      al1[len_al] = '-';
      break;
    default:
      Rcpp::stop("N-W Align out of range.");
    }
    len_al++;
  }

  al0[len_al] = '\0';
  al1[len_al] = '\0';

  // Allocate memory to alignment strings.
  char **al = (char **) malloc( 2 * sizeof(char *) ); //E
  if (al == NULL)  Rcpp::stop("Memory allocation failed.");
  al[0] = (char *) malloc(len_al+1); //E
  al[1] = (char *) malloc(len_al+1); //E
  if (al[0] == NULL || al[1] == NULL)  Rcpp::stop("Memory allocation failed.");

  // Reverse the alignment strings (since traced backwards).
  for (i=0;i<len_al;i++) {
    al[0][i] = al0[len_al-i-1];
    al[1][i] = al1[len_al-i-1];
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';

  int2nt(al[0], al[0]);
  int2nt(al[1], al[1]);

  // Store the best alignment score
  // fitaln.score = s[hw-1];
  out.push_back(std::string(al[0]));
  out.push_back(std::string(al[1]));
  // Return fitaln struct

  free(s);
  // free(bt);
  free(s1);
  free(s2);
  free(al0);
  free(al1);
  free(al);
  return out;
  //return patch::to_string(fitaln.score);
}

// [[Rcpp::export]]
Rcpp::List fit_align_fast(Rcpp::CharacterVector q, Rcpp::CharacterVector s, int match, int mismatch,
               int indel) {
  Rcpp::CharacterVector fa = fitting_alignment_fast(Rcpp::as<std::string>(q), Rcpp::as<std::string>(s),
                                          match, mismatch, indel);
  return Rcpp::List::create(
    Rcpp::_["query"] = fa[0],
    Rcpp::_["subject"] = fa[1]
  );

}


char **nw_fitting_align(const char *s1, const char *s2, int score[5][5],
                        int indel) {
  static size_t nnw = 0;
  int i, j;
  // int l, r;
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
    //d[i*ncol] = 0; // ends-free gap
    d[i*ncol] = i*indel;
    p[i*ncol] = 3;
  }

  // Fill out top rows of d, p.
  for (j = 0; j <= len2; j++) {
    d[j] = 0; // ends-free gap
    p[j] = 2;
  }

  // Fill out the body of the DP matrix.
  for (i = 1; i <= len1; i++) {
    for (j = 1; j <= len2; j++) {
      // Score for the left move.
      if (i == len1) {
        left = d[i*ncol + j-1]; // Ends-free gap.
      } else {
        left = d[i*ncol + j-1] + indel;
      }

      // Score for the up move.
      if (j == len2) {
        up = d[(i-1)*ncol + j] + indel; // Ends-free gap.
      } else {
        up = d[(i-1)*ncol + j] + indel;
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

  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;

  while ( i > 0 || j > 0 ) {
    //    Rprintf("(%i, %i): p=%i, d=%i\n", i, j, p[i*ncol + j], d[i*ncol + j]);
    switch ( p[i*ncol + j] ) {
    case 1:
      al0[len_al] = s1[--i];
      al1[len_al] = s2[--j];
      break;
    case 2:
      al0[len_al] = '-';
      al1[len_al] = s2[--j];
      break;
    case 3:
      al0[len_al] = s1[--i];
      al1[len_al] = '-';
      break;
    default:
      Rcpp::stop("N-W Align out of range.");
    }
    len_al++;
  }
  al0[len_al] = '\0';
  al1[len_al] = '\0';


  // Allocate memory to alignment strings.
  char **al = (char **) malloc( 2 * sizeof(char *) ); //E
  if (al == NULL)  Rcpp::stop("Memory allocation failed.");
  al[0] = (char *) malloc(len_al+1); //E
  al[1] = (char *) malloc(len_al+1); //E
  if (al[0] == NULL || al[1] == NULL)  Rcpp::stop("Memory allocation failed.");

  // Reverse the alignment strings (since traced backwards).
  for (i=0;i<len_al;i++) {
    al[0][i] = al0[len_al-i-1];
    al[1][i] = al1[len_al-i-1];
  }
  al[0][len_al] = '\0';
  al[1][len_al] = '\0';

  // Free allocated memory
  free(d);
  free(p);
  free(al0);
  free(al1);

  nnw++;
  return al;
}



// [[Rcpp::export]]
Rcpp::CharacterVector C_nwalign(std::string s1, std::string s2,
                                int match, int mismatch, int indel) {
  int i, j;
  char **al;
  // Make integer-ized c-style sequence strings
  char *seq1 = (char *) malloc(s1.size()+1); //E
  char *seq2 = (char *) malloc(s2.size()+1); //E
  if (seq1 == NULL || seq2 == NULL)  Rcpp::stop("Memory allocation failed.");
  nt2int(seq1, s1.c_str());
  nt2int(seq2, s2.c_str());
  // Make  c-style 2d score array
  int c_score[5][5];
  for(i=0;i<5;i++) {
    for(j=0;j<5;j++) {
      if(i==j) { c_score[i][j] = match; }
      else { c_score[i][j] = mismatch; }
    }
  }
  // Perform alignment and convert back to ACGT
  al = nw_fitting_align(seq1, seq2, c_score, indel);
  int2nt(al[0], al[0]);
  int2nt(al[1], al[1]);
  // Generate R-style return vector
  Rcpp::CharacterVector rval;
  rval.push_back(std::string(al[0]));
  rval.push_back(std::string(al[1]));
  // Clean up
  free(seq1);
  free(seq2);
  free(al[0]);
  free(al[1]);
  free(al);
  return(rval);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
//
/*** R
# # Smaller data
s_seq = "GTAGGCTTAAGGCCTA"
q_seq = "TAGATACAGG"
# # Bigger data
s_seq = "CAATCACCCCAATCCCTCAATCCTGGCCCCACGCATAGGCTAATGCCAATCGCGGCCAGGGTATAACCGCCATAACTGTGGGTCAGAAGGGATAAGTTCCACAATCCTATTTTCCTCGAGGCGCTTCGATGCGTTAACGCGTACACTCTGTCGGCCAACCGTGTGGGAGCCGAATTGGCTGGGCTGTTGAACATTCTATCAGTAGATAAACGAAGGTACATCCGAGGTTGTCGATCGACCGCGGGGTCGTAGCGCGTGCATGTTCCTTTCAGGCCCACATACTCCGGAACGGTTCATATCACGACTATTCTTGCACAATCGGACAACGGTGTACCATGGTGGACACCGTAGGAGACCAATACTGCGTAAATCATAAGCATTGGAGAGTGGACTGCTAGCGAGGCTCACCATGGAGTCTCGGTCGGCATCTCCTGACTGCTGTTCCATCGCGTTTTTCTTTTACTCACGCAATAAATCAATACCCCCTAACACAGGCCTGCTCCAGCCTTATTAAGGCCATAGTAGCTCTACATGTAGACCGAACGGAAGCACAGTTTGGTAGAAATTCTTAATCGACTATGGTCCGTGCAGGCCAAAAAAGGAATAATCTTCGAATTCTCACGCCTTCATTAGGGCGCACATGGTGGGGTAAATCACTGCACTCTGTTCGCAGTTAAGCGTTGCAATCAATATCGGCAGAACTCGGAGTCCGTATAAAGCCGCCTCAGCGTGCACACGCCCGTGCGGCACGTCATTAGACGAGGATTCCGGGGGACTGGCCTGTTCGTAATCCACTAAAACAATGGTCCTACCATCTAAAACGCACCGTGTTCCCCTCTACGGGAACCCCCTAGAT"
q_seq = "AGAGCGCAGAGAAGTCATTAGAACATGTAGCACATCGCTTATTAAGGGTCAATACCTAAAGGGCCTAACTATACGCCACACGGAACAGCTC"
fit_align_fast(q_seq, s_seq, match=1, mismatch=-1, indel=-1)
*/
