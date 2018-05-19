/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "glob_files.hpp"
#include "main.hpp"

int main(int argc, char **argv) {
  std::vector<std::string> files;
  if (argc >= 2) {
    files = glob(argv[1]);
    // for (auto f : files) std::cout << f << std::endl;
  } else if (argc < 2) {
    std::cerr << RED << "Error: \n";
    std::cerr << BOLDRED << "\tNeed input file and output file\n";
    std::cerr << RESET << "Usage:\n\t";
    std::cerr << BOLDWHITE << argv[0] << " infile.root outfile.root\n\n";
    std::cerr << RESET << std::endl;
    return 1;
  }

  std::string outfilename;
  if (argc == 2) {
    outfilename = "out.root";
  } else if (argc == 3) {
    outfilename = argv[2];
  }

  DataHandeler *dh = new DataHandeler(files, outfilename);
  dh->run();
  delete dh;
  return 0;
}
