#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// Patch bugged std::to_string() function.
// taken from the depths of StackOverflow...
namespace patch {
template < typename T > std::string to_string( const T& n ) {
  std::ostringstream stm ;
  stm << n ;
  return stm.str() ;
}
}

struct FittingAlignment {
  std::string qseq, sseq;
  int score;
  unsigned int s_start, s_end, no_mismatches, no_indels;
};

struct Backtracker {
  // Backtracking matrix
  int **r, **c;
  // r = row index, c = column index, l = layer index
};


struct Backtrack {
  int ****x;
};

FittingAlignment fitting_alignment(std::string s1, std::string s2,
                                   int match=2, int mismatch=-2, int indel=-1,
                                   int gap_open=-2, int gap_extend=-1) {
  
  // Output variables
  std::string s1_align;
  std::string s2_align;
  FittingAlignment fitaln;
  unsigned int no_mismatches = 0;
  unsigned int no_indels = 0;
  
  // Initialize scoring matrix and backtracking matrix
  int **scores, **scores_l, **scores_u;
  int min_score = -500;
  // First index 0,1,2 middle, lower, upper scoring layers
  int height = s1.length()+1;
  int width = s2.length()+1;
  unsigned int s_start = 0;
  unsigned int s_end = width-2;
  
  // Backtracking matrices
  Backtracker bt, bt_l, bt_u;

  int max_dir[3]; // store cummulative score for each directional step in middle layer
  int max_dir_u[2], max_dir_l[2]; // scores for lower and upper layers
  int max_dir_last[width-1+1]; // last node in the middle layer is special for this aln
  
  // Convert std::string into character array
  int diag; // diagonal score. mismatch or match
  int max_direction; 

  // Define scoring matrix
  scores = new int*[height];
  scores[0] = new int[width];
  scores[0][0] = 0;
  
  scores_u = new int*[height];
  scores_u[0] = new int[width];
  scores_u[0][0] = min_score;
  
  scores_l = new int*[height];
  scores_l[0] = new int[width];
  scores_l[0][0] = min_score;

  // Declare and initialize scoring and backtracking matrices
  int ****btr;
  btr = new int***[3];
  for (int l=0; l<3; ++l) { // layer l
    btr[l] = new int**[height];
    for (int r=0; r<height; ++r) { // row r
      btr[l][r] = new int*[width];
      for (int c=0; c<width; ++c) { // column c
        btr[l][r][c] = new int[3];
        btr[l][r][c][0] = -1; // target row
        btr[l][r][c][1] = -1; // target col
        btr[l][r][c][2] = -1; // target layer
      }
    }
  }
  // 0 = lower layer, 1 = middle layer, 2 = upper layer
  
  
  // std::vector<std::vector<std::vector<int> > > btr;
  // dim2: row index (height)
  // dim3: col index (width)
  // btr.r = new int**[3]; // dim4: 0 = row index (y), 1 = column index (y)
  // btr.c = new int**[3];
  
  // for (int i=0; i<3; ++i) {
  //   btr.r[i] = new int*[height];
  //   btr.r[i][0] = new int[width];
  //   btr.c[i] = new int*[height];
  //   btr.c[i][0] = new int[width];
  //   // btr[i][0][0] = new int[2];
  //   btr.r[i][0][0] = 0;
  //   btr.c[i][0][0] = 0;
  // }
  
  // Initialize Backtracking matrices
  bt.r = new int*[height];
  bt.c = new int*[height];
  bt.r[0] = new int[width];
  bt.c[0] = new int[width];
  bt.r[0][0] = -1; 
  bt.c[0][0] = -1;
  
  // bt_l.r = new int*[height];
  // bt_l.c = new int*[height];
  // bt_l.r[0] = new int[width];
  // bt_l.c[0] = new int[width];
  // bt_l.r[0][0] = -1;
  // bt_l.c[0][0] = -1;
  // 
  // bt_u.r = new int*[height];
  // bt_u.c = new int*[height];
  // bt_u.r[0] = new int[width];
  // bt_u.c[0] = new int[width];
  // bt_u.r[0][0] = -1;
  // bt_u.c[0][0] = -1;
  
  /* ---------- ROW 0 -------- */
  for (int w=1; w<width; w++) {
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
    
    // Update middle layer
    scores[0][w] = scores_u[0][w]; // zero-cost for going upper --> middle (close gap)
    btr[1][0][w][0] = 0;
    btr[1][0][w][1] = w;
    btr[1][0][w][2] = 2;
    
    // Deprecated. Remove later.
    bt.r[0][w] = 0;
    bt.c[0][w] = w-1;
  }
  
  /* --------- COLUMN 0 --------- */
  for (int h=1; h<height; h++) {
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
    scores[h] = new int[width];
    scores[h][0] = h*indel;
    
    bt.r[h] = new int[width];
    bt.c[h] = new int[width];
    bt.r[h][0] = h-1;
    bt.c[h][0] = 0;
  }
  
  // Process scores 
  
  // Fill the rest of the elements
  for (int h=1; h<height; h++) {
    for (int w=1; w<width; w++) {
      // Update scoring matrix
      if (h == height-1 && w == width-1) {
        // last element. add shortcuts from left steps
        //std::cout << "processing last step...\n";
        for (int k=0; k<width-1; k++) {
          max_dir_last[k] = scores[h][k];
        }
        diag = s1[h-1] == s2[w-1] ? match : mismatch;
        max_dir_last[width-1] = scores[h-1][w-1] + diag; // diag move to last node
        max_dir_last[width-1+1] = scores[h][w-1] + indel; // stop down to last node
        max_direction = std::distance(max_dir_last, 
                                      std::max_element(max_dir_last, 
                                                       max_dir_last + sizeof(max_dir_last) / sizeof(int)));
        
        // Update backtracking matrix
        if (max_direction < width-1) { // shortcuts to last node
          // bt[h][w] = bt[h][max_direction];
          bt.r[h][w] = h;
          bt.c[h][w] = max_direction;
          btr[0][h][w][0] = h;
          btr[0][h][w][1] = max_direction;
        } else if (max_direction == width-1) { // diag move
          // bt[h][w] = bt[h-1][w-1];
          bt.r[h][w] = h-1;
          bt.c[h][w] = w-1;
          btr[0][h][w][0] = h-1;
          btr[0][h][w][1] = w-1;
        } else if (max_direction == width) { // step down
          // bt[h][w] = bt[h-1][w];
          bt.r[h][w] = h-1;
          bt.c[h][w] = w;
          btr[0][h][w][0] = h-1;
          btr[0][h][w][1] = w;
        }
        // Update last node in the scoring matrix
        scores[h][w] = max_dir_last[max_direction];
      } else {
        max_dir[0] = scores[h-1][w] + indel; // step down
        max_dir[1] = scores[h][w-1] + indel; // step right
        diag = s1[h-1] == s2[w-1] ? match : mismatch;
        max_dir[2] = scores[h-1][w-1] + diag;
        // Get the index of the maximum value
        max_direction = std::distance(max_dir, 
                                      std::max_element(max_dir, 
                                                       max_dir + sizeof(max_dir) / sizeof(int)));
        scores[h][w] = max_dir[max_direction];
        if (max_direction == 0) { // step down
          bt.r[h][w] = h-1;
          bt.c[h][w] = w;
          // btr[0][h][w][0] = h-1;
          // btr[0][h][w][1] = w;
        } else if (max_direction == 1) { // step right
          bt.r[h][w] = h;
          bt.c[h][w] = w-1;
          // btr[0][h][w][0] = h;
          // btr[0][h][w][1] = w-1;
        } else if (max_direction == 2) { // diag move
          bt.r[h][w] = h-1;
          bt.c[h][w] = w-1;
          // btr[0][h][w][0] = h-1;
          // btr[0][h][w][1] = w-1;
        }
      }
    }
  }
  
  // DEBUG: print backtracking matrices
  for (int h=0; h<height; ++h) {
    for (int w=0; w<width; ++w) {
      std::cout << "bt[" << h << "][" << w << "] row, col: " << bt.r[h][w] << ", " << bt.c[h][w] << "\n";
    }
  }
  
  // Assemble the aligned sequences by backtracking
  int r = height-1; // row index
  int c = width-1;  // col index
  
  std::vector<char> s1a, s2a;
  s1a.reserve(width-1);
  s2a.reserve(width-1);
  
  while (r > 0 || c > 0) {
    // Which way to backtrack?
    std::cout << "row, col = " << r << ", " << c << "\n";
    if (bt.r[r][c] == r && bt.c[r][c] < c) { // left step
    //if (btr[0][r][c][0] == r && btr[0][r][c][1] < c) {
      std::cout << "left\n";
      for (int s=c; s > bt.c[r][c]; --s) {
        s1a.push_back('-');
        //std::cout << "query adding " << '-' << "\n";
        s2a.push_back(s2[s-1]);
        //std::cout << "subject adding " << s2[s-1] << "\n";
      }
    } else if (bt.r[r][c] < r && bt.c[r][c] == c) { // up step
    //} else if (btr[0][r][c][0] < r && btr[0][r][c][1] == c) {
      std::cout << "up\n";
      s1a.push_back(s1[r-1]);
      s2a.push_back('-');
    } else if (bt.r[r][c] < r && bt.c[r][c] < c) { // diag step
    //} else if (btr[0][r][c][0] < r && btr[0][r][c][1] < c) {
      std::cout << "diag\n";
      s1a.push_back(s1[r-1]);
      s2a.push_back(s2[c-1]);
    }
    // Update row and col indexes    
    r = bt.r[r][c];
    c = bt.c[r][c];
    // r = btr[0][r][c][0];
    // c = btr[0][r][c][1];
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
  
  // Print the scores matrix for debugging purposes
  // for (int i=0; i<height; ++i) {
  //   for (int j=0; j<width; ++j) {
  //     std::cout << "scores[" << i << "][" << j << "] = " << scores[i][j] << "\n";
  //   }
  // }
  
  // Store the best alignment score
  fitaln.score = scores[height-1][width-1];
  fitaln.qseq = s1_align;
  fitaln.sseq = s2_align;
  // std::cout << "Best alignment score: " << fitaln.score << "\n";
  
  // Clean-up all raw arrays with delete command to release memory
  
  // Return fitaln struct
  return fitaln;
  //return patch::to_string(fitaln.score);
}




// [[Rcpp::export]]
List fit_aln(CharacterVector q, CharacterVector s) {
  FittingAlignment fa = fitting_alignment(as<std::string>(q), as<std::string>(s));
  return List::create(
    _["query"] = fa.qseq,
    _["subject"] = fa.sseq,
    _["score"] = fa.score
  );
    
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
// 
/*** R
fit_aln("TAGATA", "GTAGGCTTAAGG")
*/
