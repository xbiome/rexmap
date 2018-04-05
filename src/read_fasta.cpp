// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
// for mmap:
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h> // read
#include <boost/algorithm/string.hpp>

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



const char* map_file(const char* fname, size_t& length);

int main_function()
{
  size_t length;
  const char* f = map_file("/data1/igor/himap/16s_uniq_refseq_fullgen_bact_arch_2017-08-24.fasta", length);
  const char* l = f + length;
  
  uintmax_t m_numLines = 0;
  while (f && f!=l)
    if ((f = static_cast<const char*>(memchr(f, '\n', l-f))))
      m_numLines++, f++;
    
    std::cout << "m_numLines = " << m_numLines << "\n";
  return 0;
}

void handle_error(const char* msg) {
  perror(msg); 
  exit(255);
}

const char* map_file(const char* fname, size_t& length)
{
  int fd = open(fname, O_RDONLY);
  if (fd == -1)
    handle_error("open");
  
  // obtain file size
  struct stat sb;
  if (fstat(fd, &sb) == -1)
    handle_error("fstat");
  
  length = sb.st_size;
  
  const char* addr = static_cast<const char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
  if (addr == MAP_FAILED)
    handle_error("mmap");
  
  // TODO close fd at some point in time, call munmap(...)
  return addr;
}

// Copy-pasted from https://stackoverflow.com/questions/17925051/fast-textfile-reading-in-c
static uintmax_t wc(char const *fname)
{
  static const unsigned int BUFFER_SIZE = 16*1024;
  int fd = open(fname, O_RDONLY);
  if(fd == -1)
    handle_error("open");
  
  /* Advise the kernel of our access pattern.  */
  posix_fadvise(fd, 0, 0, 1);  // FDADVICE_SEQUENTIAL
  
  char buf[BUFFER_SIZE + 1];
  uintmax_t lines = 0;
  
  while(size_t bytes_read = read(fd, buf, BUFFER_SIZE))
  {
    if(bytes_read == (size_t)-1)
      handle_error("read failed");
    if (!bytes_read)
      break;
    
    for(char *p = buf; (p = (char*) memchr(p, '\n', (buf + bytes_read) - p)); ++p)
      ++lines;
  }
  
  return lines;
}

std::vector<std::string> wc2(const char* fname)
{
  static const unsigned int BUFFER_SIZE = 16*1024;
  std::vector<std::string> out;
  int fd = open(fname, O_RDONLY);
  if(fd == -1)
    handle_error("open");
  
  /* Advise the kernel of our access pattern.  */
  posix_fadvise(fd, 0, 0, 1);  // FDADVICE_SEQUENTIAL
  
  char buf[BUFFER_SIZE + 1];
  uintmax_t lines = 0;
  
  while(size_t bytes_read = read(fd, buf, BUFFER_SIZE))
  {
    if(bytes_read == (size_t)-1)
      handle_error("read failed");
    if (!bytes_read)
      break;
    
    for(char *p = buf; (p = (char*) memchr(p, '\n', (buf + bytes_read) - p)); ++p) {
      ++lines;
      // std::cout << p << " xx \n";
      out.push_back(p);
    }
      
  }
  
  return out;
}


// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp::export]]
CharacterVector test_fasta() {
  return main_function();
}
// [[Rcpp::export]]
std::vector<std::string> file_reader(std::string filename) {
  // is there a way to avoid this conversion?
  return wc2(filename.c_str());
}

std::string slurp(std::ifstream& in) {
  return static_cast<std::stringstream const&>(std::stringstream() << in.rdbuf()).str();
}

// [[Rcpp::export]]
std::vector<std::string> read_file(std::string filename) {
  std::ifstream in_file;
  std::string s;
  std::vector<std::string> out;
  in_file.open(filename.c_str());
  s = slurp(in_file);
  boost::split(out, s, boost::is_any_of("\n"));
  return out;
}
  

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

my.read.lines2=function(fname) {
   s = file.info( fname )$size
   buf = readChar( fname, s, useBytes=T)
   x = strsplit( buf,"\n",fixed=T,useBytes=T)[[1]]
}

filename = "/data1/igor/himap/16s_uniq_refseq_fullgen_bact_arch_2017-08-24.fasta"
testfq = '/data1/igor/app/pcr_pacbio/test.fastq'

# timesTwo(42)
# system.time(my.read.lines2("/data1/igor/himap/16s_uniq_refseq_fullgen_bact_arch_2017-08-24.fasta"))
system.time(test_fasta())
# tmp = file_reader(filename)

*/
