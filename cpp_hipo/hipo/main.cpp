/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "main.hpp"

using namespace std;

int main(int argc, char **argv) {
  if (argc == 2) {
    char *infilename = argv[1];
    char out[] = "out.root";
    datahandeler(infilename, out);
  } else if (argc == 3) {
    char *infilename = argv[1];
    char *outfilename = argv[2];
    datahandeler(infilename, outfilename);
  } else {
    std::cerr << RED << "Error: \n";
    std::cerr << BOLDRED << "\tNeed input file and output file\n";
    std::cerr << RESET << "Usage:\n\t";
    std::cerr << BOLDWHITE << argv[0] << " infile.root outfile.root\n\n";
    std::cerr << RESET << std::endl;
  }

  return 0;
}
