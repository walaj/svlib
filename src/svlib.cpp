/* svlib - Tools for producing/benchmarking structural variants
 * Copyright 2016 Jeremiah Wala
 * Written by Jeremiah Wala (jwala@broadinstitute.org)
 * Released under the MIT license
 */

#include <iostream>
#include "seq2vcf.h"
#include "simulate.h"
#include "realigntest.h"

static const char *USAGE_MESSAGE =
"Program: svlib \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ]\n"
"Usage: snowman <command> [options]\n\n"
"Commands:\n"
"           seqtovcf          Convert an aligned long sequence BAM/SAM/CRAM to a VCF\n"
"           sim               Simulate indels and rearrangements on a genome\n"
"           realigntest       Simulate rearrangements and contigs and test realignment performance\n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

int main(int argc, char** argv) {

  if (argc <= 1) {
    std::cerr << USAGE_MESSAGE;
    return 0;
  } else {
    std::string command(argv[1]);
    if (command == "help" || command == "--help") {
      std::cerr << USAGE_MESSAGE;
      return 0;
    } else if (command == "seqtovcf") {
      runSeqToVCF(argc -1, argv + 1);
    } else if (command == "sim") {
      runSimulation(argc - 1, argv + 1);
    } else if (command == "realigntest") {
      runRealignTest(argc - 1, argv + 1);
    }
    else {
      std::cerr << USAGE_MESSAGE;
      return 0;
    }
  } 

  return 0;

}
