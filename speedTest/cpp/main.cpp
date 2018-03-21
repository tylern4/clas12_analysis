/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "main.hpp"

using namespace std;

int main(int argc, char **argv) {
  gStyle->SetOptFit(1111);
  auto start = std::chrono::high_resolution_clock::now();
  if (argc == 2) {
    char *infilename = argv[1];
    test(infilename, "out.root");
  } else if (argc == 3) {
    char *infilename = argv[1];
    char *outfilename = argv[2];
    test(infilename, outfilename);
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << duration.count() << std::endl;

  return 0;
}
