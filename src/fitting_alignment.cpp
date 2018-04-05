#include <Rcpp.h>
using namespace Rcpp;


struct FittingAlignment {
  std::string qseq, sseq;
  int score;
  unsigned int s_start, s_end, no_mismatches;
  unsigned int lefthang, righthang, no_matches;
  unsigned int q_no_gapopen, q_no_gapext, s_no_gapopen, s_no_gapext;
};

struct Backtracker {
  // Backtracking matrix
  int *r, *c;
  // r = row index, c = column index, l = layer index
};

FittingAlignment fitting_alignment(std::string s1, std::string s2,
                                   int match, int mismatch,
                                   int gap_open, int gap_extend) {
  
  // Output variables
  std::string s1_align;
  std::string s2_align;
  FittingAlignment fitaln;
  unsigned int left_hang = 0;
  unsigned int right_hang = 0;
  unsigned int no_matches = 0;
  unsigned int no_mismatches = 0;
  unsigned int q_no_gapopen = 0;
  unsigned int q_no_gapext = 0;
  unsigned int s_no_gapopen = 0;
  unsigned int s_no_gapext = 0;
  bool left_h = TRUE;
  bool right_h = TRUE;
  
  // Initialize scoring matrix and backtracking matrix
  int **scores, **scores_l, **scores_u;
  int min_score = -500;
  // First index 0,1,2 middle, lower, upper scoring layers
  int nrow = s1.length()+1;
  int ncol = s2.length()+1;
  unsigned int s_start = 0;
  unsigned int s_end = ncol-2;
  
  int max_dir[3]; // store cummulative score for each directional step in middle layer
  int max_dir_u[2], max_dir_l[2]; // scores for lower and upper layers
  int max_dir_last[ncol-1+1]; // last node in the middle layer is special for this aln
  
  // Convert std::string into character array
  int diag; // diagonal score. mismatch or match
  int max_direction; 

  // Define scoring matrix
  scores = new int*[nrow];
  scores_u = new int*[nrow];
  scores_l = new int*[nrow];
  // Declare and initialize scoring and backtracking matrices
  int ****btr;
  btr = new int***[3];
  for (int l=0; l<3; ++l) { // layer l
    btr[l] = new int**[nrow];
    for (int r=0; r<nrow; ++r) { // row r
      btr[l][r] = new int*[ncol];
      for (int c=0; c<ncol; ++c) { // column c
        btr[l][r][c] = new int[3];
        btr[l][r][c][0] = -1; // target row
        btr[l][r][c][1] = -1; // target col
        btr[l][r][c][2] = -1; // target layer
      }
    }
  }
  
  for (int r=0; r<nrow; ++r) {
    scores_l[r] = new int[ncol];
    scores_u[r] = new int[ncol];
    scores[r] = new int[ncol];
    for (int c=0; c<ncol; ++c) {
      scores_l[r][c] = 0;
      scores_u[r][c] = 0;
      scores[r][c] = 0;
    }
  }
  scores_l[0][0] = min_score;
  scores_u[0][0] = min_score;
  
  // 0 = lower layer, 1 = middle layer, 2 = upper layer
  
  /* ---------- ROW 0 -------- */
  for (int w=1; w<ncol; w++) {
    // Update lower layer
    scores_l[0][w] = min_score; // nowhere to backtrack. we should never be here.

    // Update upper layer
    max_dir_u[0] = scores_u[0][w-1] + gap_extend; // gap extension in upper
    max_dir_u[1] = scores[0][w-1] + gap_open; // gap open from middle
    btr[2][0][w][0] = 0;
    btr[2][0][w][1] = w-1;
    if (max_dir_u[1] > max_dir_u[0]) { // middle -> upper layer jump
      btr[2][0][w][2] = 1;
      scores_u[0][w] = max_dir_u[1];
    } else { // upper -> upper layer transition
      // always prefer extension over open even if penalty is the same
      btr[2][0][w][2] = 2;
      scores_u[0][w] = max_dir_u[0];
    }
    // For fitting alignment set these to zero. If we are here then
    // the only thing left is the left gap which should not be counter.
    scores_u[0][w] = 0;
    btr[2][0][w][0] = 0;
    btr[2][0][w][1] = 0;
    btr[2][0][w][2] = 1;

    // Update middle layer
    // scores[0][w] = scores_u[0][w]; // zero-cost for going upper --> middle (close gap)
    // btr[1][0][w][0] = 0;
    // btr[1][0][w][1] = w;
    // btr[1][0][w][2] = 2;
    scores[0][w] = 0; // these are just shortcuts from the 0,0 node, take this since we
    // should never have positive gap opening penalty.
    btr[1][0][w][0] = 0;
    btr[1][0][w][1] = 0;
    btr[1][0][w][2] = 1;

  }

  /* --------- COLUMN 0 --------- */
  for (int h=1; h<nrow; h++) {
    // Lower layer
    max_dir_l[0] = scores_l[h-1][0] + gap_extend;
    max_dir_l[1] = scores[h-1][0] + gap_open;
    btr[0][h][0][0] = h-1;
    btr[0][h][0][1] = 0;
    if (max_dir_l[1] > max_dir_l[0]) { // Opening better
      btr[0][h][0][2] = 1;
      scores_l[h][0] = max_dir_l[1];
    } else { // Extending better or equal to opening
      btr[0][h][0][2] = 0;
      scores_l[h][0] = max_dir_l[0];
    }

    // Upper layer: should never be able to get here
    scores_u[h][0] = min_score; // leave backtrackers undetermined (-1)

    // Middle layer: the only way to get here is from lower layer
    scores[h][0] = scores_l[h][0]; // zero-cost for lower-->middle
    btr[1][h][0][0] = h;
    btr[1][h][0][1] = 0;
    btr[1][h][0][2] = 0;
  }

  

  // Process scores 
  
  // Fill the rest of the elements
  for (int h=1; h<nrow; h++) { // h = row index
    for (int w=1; w<ncol; w++) { // w = column index
      // Update scoring matrix
      if (h == nrow-1 && w == ncol-1) {
        // last element. add shortcuts from left steps
        //std::cout << "processing last step...\n";
        for (int k=0; k<ncol-1; k++) {
          max_dir_last[k] = scores[h][k];
        }
        diag = s1[h-1] == s2[w-1] ? match : mismatch;
        max_dir_last[ncol-1] = scores[h-1][w-1] + diag; // diag move to last node
        max_dir_last[ncol-1+1] = scores_u[h][w]; // stop down to last node
        max_direction = std::distance(max_dir_last,
                                      std::max_element(max_dir_last,
                                                       max_dir_last + sizeof(max_dir_last) / sizeof(int)));

        // Update backtracking matrix
        if (max_direction < ncol-1) { // shortcuts to last node (free)
          btr[1][h][w][0] = h;
          btr[1][h][w][1] = max_direction;
          btr[1][h][w][2] = 1;
        } else if (max_direction == ncol-1) { // diag move
          btr[1][h][w][0] = h-1;
          btr[1][h][w][1] = w-1;
          btr[1][h][w][2] = 1;
        } else if (max_direction == ncol) { // step down
          btr[1][h][w][0] = h-1;
          btr[1][h][w][1] = w;
          btr[1][h][w][2] = 1;
        }
        // Update last node in the scoring matrix
        scores[h][w] = max_dir_last[max_direction];
      } else { // element that is not the last one

        // Update lower layer first (it depends on row-1 )
        max_dir_l[0] = scores_l[h-1][w] + gap_extend;
        max_dir_l[1] = scores[h-1][w] + gap_open;
        btr[0][h][w][0] = h-1;
        btr[0][h][w][1] = w;
        if (max_dir_l[1] > max_dir_l[0]) { // Gap open wins
          btr[0][h][w][2] = 1;
          scores_l[h][w] = max_dir_l[1];
        } else { // Gap extend wins
          btr[0][h][w][2] = 0;
          scores_l[h][w] = max_dir_l[0];
        }
        // Update upper layer (it depends on col-1 from mid and up)
        max_dir_u[0] = scores_u[h][w-1] + gap_extend;
        max_dir_u[1] = scores[h][w-1] + gap_open;
        btr[2][h][w][0] = h;
        btr[2][h][w][1] = w-1;
        if (max_dir_u[1] > max_dir_u[0]) { // Gap open better
          btr[2][h][w][2] = 1;
          scores_u[h][w] = max_dir_u[1];
        } else { // Gap extend better
          btr[2][h][w][2] = 2;
          scores_u[h][w] = max_dir_u[0];
        }
        // Update middle layer
        max_dir[0] = scores_u[h][w]; // upper->middle
        max_dir[1] = scores_l[h][w]; // lower->middle
        diag = s1[h-1] == s2[w-1] ? match : mismatch; // diagonal move
        max_dir[2] = scores[h-1][w-1] + diag;
        // Get the index of the maximum value
        max_direction = std::distance(
          max_dir,
          std::max_element(max_dir, max_dir + sizeof(max_dir) / sizeof(int)));
        scores[h][w] = max_dir[max_direction];
        if (max_direction == 2) { // diagonal move
          btr[1][h][w][0] = h-1;
          btr[1][h][w][1] = w-1;
          btr[1][h][w][2] = 1;
        } else { // same row, col, just different layer
          btr[1][h][w][0] = h;
          btr[1][h][w][1] = w;
          if (max_direction == 0) { // upper->middle
            btr[1][h][w][2] = 2;
          } else { // max_direction == 1, lower ->middle
            btr[1][h][w][2] = 0;
          }
        }

      }
    }
  }


  // Assemble the aligned sequences by backtracking
  int r = nrow-1; // row index
  int c = ncol-1;  // col index
  int l = 1; // start from last node in middle layer
  int temp[] = {0, 0, 0};

  std::vector<char> s1a, s2a;
  s1a.reserve(ncol-1);
  s2a.reserve(ncol-1);

  while (!(r==0 && c==0)) {
    // Which way to backtrack?
    // std::cout << "btr: " << l << ", " << r << ", " << c << "\n";
    if (btr[l][r][c][0] < r && btr[l][r][c][1] < c) { // diagonal reverse ...[2] == l here
      // std::cout << "diag\n";
      s1a.push_back(s1[r-1]);
      s2a.push_back(s2[c-1]);
      if (s1a.back() != s2a.back()) no_mismatches++;
      else no_matches++;
    } else { // non-diagonal moves
      if (l == 2 && btr[l][r][c][2] == 1) { // middle -> upper gap open
        q_no_gapopen++;
        s1a.push_back('-');
        s2a.push_back(s2[c-1]);
      } else if (l == 0 && btr[l][r][c][2] == 1) { // middle -> lower gap open
        s_no_gapopen++;
        s1a.push_back(s1[r-1]);
        s2a.push_back('-');
      } else if (l == 0 && btr[l][r][c][2] == 0) { // lower->lower gap extend
        s_no_gapext++;
        s1a.push_back(s1[r-1]);
        s2a.push_back('-');
      } else if (l == 2 && btr[l][r][c][2] == 2) { // upper-> upper gap extend
        q_no_gapext++;
        s1a.push_back('-');
        s2a.push_back(s2[c-1]);
      } else if (l == 1 && btr[l][r][c][2] == 1 && btr[l][r][c][1] < c) { // shortcut jump
        for (int s=c; s > btr[l][r][c][1]; --s) {
          // std::cout << "processing shortcut...\n";
          s1a.push_back('-');
          // std::cout << "query adding " << '-' << "\n";
          s2a.push_back(s2[s-1]);
          // std::cout << "subject adding " << s2[s-1] << "\n";
          if (r == nrow-1) { // shortcut for gap at the end
            s_end--;
            right_hang++;
          } else { // shortcut for gap at the beginning
            s_start++;
            left_hang++;
          }
        }
      }
    }
    // Update row and col indexes
    temp[0] = btr[l][r][c][0];
    temp[1] = btr[l][r][c][1];
    temp[2] = btr[l][r][c][2];
    r = temp[0];
    c = temp[1];
    l = temp[2];
  }

  // Now assemble alignment strings by poping the elements
  // of the stack in reverse.
  while (!s1a.empty()) {
    s1_align.push_back(s1a.back());
    s1a.pop_back();
  }
  while (!s2a.empty()) {
    s2_align.push_back(s2a.back());
    s2a.pop_back();
  }

  // Store the best alignment score
  fitaln.score = scores[nrow-1][ncol-1];
  fitaln.qseq = s1_align;
  fitaln.sseq = s2_align;
  fitaln.no_mismatches = no_mismatches;
  fitaln.lefthang = left_hang;
  fitaln.righthang = right_hang;
  fitaln.no_matches = no_matches;
  fitaln.q_no_gapopen = q_no_gapopen;
  fitaln.q_no_gapext = q_no_gapext;
  fitaln.s_no_gapopen = s_no_gapopen;
  fitaln.s_no_gapext = s_no_gapext;
  // fitaln.s_start = s_start;
  // fitaln.s_end = s_end;
  // std::cout << "Best alignment score: " << fitaln.score << "\n";
  
  // Clean-up all raw arrays with delete command to release memory
  
  // Return fitaln struct
  return fitaln;
}

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
      Rprintf("invalid character in input:%c.\n",*iseq);
    *oseq = '\0';
    }
  }
  *oseq = '\0';
  return;
}

/* Converts seq in index form back to nts. */
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



// Faster fitting alignment with a single indel penalty
// [[Rcpp::export]]
CharacterVector fitting_alignment_fast(std::string s1s, std::string s2s,
                                   int match=1, int mismatch=-1, int indel=-1) {
  
  // Output variables
  CharacterVector out;
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
  int i, j, k, ri, ci, i_next;
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
  int score_last[ncol+1];
  
  //int max_dir[3]; // store cummulative score for each directional step in middle layer
  // int max_dir_u[2], max_dir_l[2]; // scores for lower and upper layers
  int max_dir_last[ncol+1]; // last node in the middle layer is special for this aln
  
  // Convert std::string into character array
  //int diag; // diagonal score. mismatch or match
  int max_direction; 
  
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
  

  // // Iterate over every other element except for the last
  // for (j=0; j < (nrow-1)*(ncol-1) - 1; j++) { // flattened inner index
  //   // convert to flattened outer index
  //   ri = j/(ncol-1) + 1; // outer row
  //   ci = j%(ncol-1) + 1; // outer col
  //   i = ri*ncol + ci;
  //   // std::cout << "j=" << j << ", i=" << i << ", row: " << ri << ", col: " << ci << "\n";
  //   max_dir[0] = s[i-ncol] + indel; // step down from prev row
  //   max_dir[1] = s[i-1] + indel;     // step right from prev col
  //   max_dir[2] = s[i-ncol-1] + (s1[ri-1] == s2[ci-1] ? match : mismatch); // diag move
  //   // std::cout << "max_dir: " << max_dir[0] << ", " << max_dir[1] << ", " << max_dir[2] << "\n";
  //   max_direction = max_index(max_dir, 3);
  //   // std::cout << "max_direction: " << max_direction << "\n";
  //   s[i] = max_dir[max_direction];
  //   if (max_direction == 0) { // step down
  //     bt[i] = i-ncol;
  //   } else if (max_direction == 1) {
  //     bt[i] = i-1;
  //   } else {
  //     bt[i] = i-ncol-1;
  //   }
  // }
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
  
  
  // // Finalize the last element. Allow 0 cost jump from any last row element.
  // for (j=0; j<ncol-1; j++) {
  //   max_dir_last[j] = s[j+(nrow-1)*ncol];
  // }
  // // Last two elements of this array are diagonal and step down
  // max_dir_last[ncol-1] = s[hw-1-ncol-1] + (s1[nrow-1] == s2[ncol-1] ? match : mismatch);
  // max_dir_last[ncol] = s[hw-1-ncol] + indel;
  // max_direction = max_index(max_dir_last, ncol+1);
  // s[hw-1] = max_dir_last[max_direction];
  // if (max_direction < ncol-1) { // one of the shortcuts is the best
  //   bt[hw-1] = max_direction + (nrow-1)*ncol;
  // } else if (max_direction == ncol-1) {
  //   bt[hw-1] = j-ncol-1;
  // } else {
  //   bt[hw-1] = j-ncol;
  // }
  
  // Print backtracking and scoring matrices
  // for (int r=0; r<nrow; r++) {
  //   for (int c=0; c<ncol; c++) {
  //     std::cout << "score [row " << r << ", col " << c << "] = " << s[i] << "\n";
  //     std::cout << "bt [row " << r << ", col " << c << "] = " << bt[i] << "\n";
  //   }
  // }
  

  // // Assemble the aligned sequences by backtracking
  // std::vector<char> s1a, s2a;
  // s1a.reserve(ncol-1);
  // s2a.reserve(ncol-1);
  // 
  // i = hw-1;
  // while (i > 0) {
  //   i_next = bt[i];
  //   ri = i / ncol - 1;
  //   ci = i % ncol - 1;
  //   // std::cout << "row: " << ri << ", col: " << ci << "\n";
  //   if (i_next == i-ncol) { // step up
  //     s1a.push_back(s1[ri]);
  //     s2a.push_back('-');
  //   } else if (i_next == i-ncol-1) { // step diag
  //     s1a.push_back(s1[ri]);
  //     s2a.push_back(s2[ci]);
  //   } else { // left move (do all of them in case of shortcuts)
  //     for (j=0; j<i-i_next; j++) {
  //       s1a.push_back('-');
  //       s2a.push_back(s2[ci-j]);
  //     }
  //   }
  //   i = i_next;
  // }

  char *al0 = (char *) malloc((len1+len2+1) * sizeof(char));
  char *al1 = (char *) malloc((len1+len2+1) * sizeof(char));
  if(al0 == NULL || al1 == NULL) Rcpp::stop("Memory allocation failed.");

  // Trace back over p to form the alignment.
  size_t len_al = 0;
  i = len1;
  j = len2;

  // i = hw-1;
  // while (i > 0) {
  //   i_next = bt[i];
  //   ri = i / ncol - 1;
  //   ci = i % ncol - 1;
  //   // std::cout << "row: " << ri << ", col: " << ci << "\n";
  //   if (i_next == i-ncol) { // step up
  //     al0[len_al] = s1[ri];
  //     al1[len_al] = '-';
  //     len_al++;
  //   } else if (i_next == i-ncol-1) { // step diag
  //     al0[len_al] = s1[ri];
  //     al1[len_al] = s2[ci];
  //     len_al++;
  //   } else { // left move (do all of them in case of shortcuts)
  //     for (j=0; j<i-i_next; j++) {
  //       al0[len_al] = '-';
  //       al1[len_al] = s2[ci-j];
  //       len_al++;
  //     }
  //   }
  //   i = i_next;
  // }
  // 1 = left, 2 = diag, 3 = up
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
  
  // Now assemble alignment strings by poping the elements
  // of the stack in reverse.
  // while (!s1a.empty()) {
  //   s1_align.push_back(s1a.back());
  //   s1a.pop_back();
  // }
  // while (!s2a.empty()) {
  //   s2_align.push_back(s2a.back());
  //   s2a.pop_back();
  // }
  
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


// // C_nwalign
// // [[Rcpp::export]]
// Rcpp::CharacterVector C_nwalign(std::string s1, std::string s2, int match, int mismatch, int indel);
// RcppExport SEXP _himap_C_nwalign(SEXP s1SEXP, SEXP s2SEXP, SEXP matchSEXP, SEXP mismatchSEXP, SEXP indelSEXP) {
//   BEGIN_RCPP
//   Rcpp::RObject rcpp_result_gen;
//   Rcpp::RNGScope rcpp_rngScope_gen;
//   Rcpp::traits::input_parameter< std::string >::type s1(s1SEXP);
//   Rcpp::traits::input_parameter< std::string >::type s2(s2SEXP);
//   Rcpp::traits::input_parameter< int >::type match(matchSEXP);
//   Rcpp::traits::input_parameter< int >::type mismatch(mismatchSEXP);
//   Rcpp::traits::input_parameter< int >::type indel(indelSEXP);
//   rcpp_result_gen = Rcpp::wrap(fitting_alignment_fast(s1, s2, match, mismatch, indel));
//   return rcpp_result_gen;
//   END_RCPP
// }

// [[Rcpp::export]]
List fit_align(CharacterVector q, CharacterVector s, int match, int mismatch,
               int gap_open, int gap_extend) {
  FittingAlignment fa = fitting_alignment(as<std::string>(q), as<std::string>(s),
                                          match, mismatch, gap_open, gap_extend);
  return List::create(
    _["query"] = fa.qseq,
    _["subject"] = fa.sseq,
    _["score"] = fa.score,
    _["no_mismatches"] = fa.no_mismatches,
    _["no_matches"] = fa.no_matches,
    _["left"] = fa.lefthang,
    _["right"] = fa.righthang,
    _["q_no_gapopen"] = fa.q_no_gapopen,
    _["q_no_gapext"] = fa.q_no_gapext,
    _["s_no_gapopen"] = fa.s_no_gapopen,
    _["s_no_gapext"] = fa.s_no_gapext
    //_["s_start"] = fa.s_start+1,
    //_["s_end"] = fa.s_end+1
  );
    
}

// [[Rcpp::export]]
List fit_align_fast(CharacterVector q, CharacterVector s, int match, int mismatch,
               int indel) {
  CharacterVector fa = fitting_alignment_fast(as<std::string>(q), as<std::string>(s),
                                          match, mismatch, indel);
  return List::create(
    _["query"] = fa[0],
    _["subject"] = fa[1]
    // _["score"] = fa.score
    // _["no_mismatches"] = fa.no_mismatches,
    // _["no_matches"] = fa.no_matches,
    // _["left"] = fa.lefthang,
    // _["right"] = fa.righthang,
    // _["q_no_gapopen"] = fa.q_no_gapopen,
    // _["q_no_gapext"] = fa.q_no_gapext,
    // _["s_no_gapopen"] = fa.s_no_gapopen,
    // _["s_no_gapext"] = fa.s_no_gapext
    //_["s_start"] = fa.s_start+1,
    //_["s_end"] = fa.s_end+1
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
