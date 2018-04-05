int main() {
    
}


char **himap_nwalign_endsfree(const char *s1, const char *s2, int score[5][5], int gap_p, int band) {
    static size_t nnw = 0;
    int i, j;
    int l, r;
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