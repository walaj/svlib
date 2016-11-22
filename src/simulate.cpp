 #include "simulate.h"

 #include <getopt.h>
 #include <string>
 #include <sstream>
 #include <fstream>
 #include <iostream>
 #include <cstdlib>
 #include <algorithm>

 #include "SeqLib/GenomicRegion.h"
 #include "SeqLib/BamReader.h"

 #include "SeqFrag.h"
// #include "SimGenome.h"
 #include "PowerLawSim.h"
// #include "SimTrainerWalker.h"

 static SeqLib::BamReader bwalker;
 static faidx_t * findex;

 namespace opt {

   static std::string refgenome = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
   static uint32_t seed = 0;
   static std::string bam;
   static int nbreaks = 10;
   static double plaw = 1.0001;
   static int nindels = 10;
   static bool scramble = false; // add scrambled inserts

   static std::string string_id = "noid";
   static int viral_count = 0;
   static std::string blacklist;
 }

 static const char *SIMULATE_USAGE_MESSAGE =
 "Usage: svlib simulate\n\n"
 "  Description: Simulate indels and rearrangementss on a reference genome\n"
 "\n"
 "  Required input\n"
 "  -G, --reference-genome               Indexed ref genome for BWA-MEM. Default (Broad):\n"
 "  -b, --bam                            BAM file to train the simulation with\n"
 "  General options\n"
 "  -s, --seed                           Seed for the random number generator [0]\n"
 "  -a, --string-id                      String to prepend output files with [noid]\n"
 "  -l, --power-law                      Power law to simulate sizes from (x^-l). [1.0001]\n"
 "  -R, --num-rearrangements             Number of rearrangements to simulate [10]\n"
 "  -X, --num-indels                     Number of indels to simulate [10]\n"
 "      --blacklist                      BED file specifying blacklist regions not to put breaks in\n"
 "\n";

enum {
  OPT_BLACKLIST,
  OPT_SCRAMBLE
};

 static const char* shortopts = "hG:c:n:s:k:b:E:I:D:R:X:a:f:M:l:";
 static const struct option longopts[] = {
   { "help",                 no_argument, NULL, 'h' },
   { "reference-genome",     required_argument, NULL, 'G' },
   { "string-id",            required_argument, NULL, 'A' },
   { "seed",                 required_argument, NULL, 's' },
   { "bam",                  required_argument, NULL, 'b' },
   { "num-rearrangements",   required_argument, NULL, 'R' },
   { "num-indels",           required_argument, NULL, 'X' },
   { "viral-integration",    required_argument, NULL, 'M' },
   { "power-law",            required_argument, NULL, 'l' },
   { "blacklist",            required_argument, NULL, OPT_BLACKLIST},
   { "add-scrambled-inserts",no_argument, NULL, OPT_SCRAMBLE},
   { NULL, 0, NULL, 0 }
 };

// helper for opening files
template <typename T>
void fopen(const std::string& s, T& o) {
  o.open(s.c_str(), std::ios::out);
}

 void runSimulation(int argc, char** argv) {

   parseSimulationOptions(argc, argv);

   // seed the RNG
   srand(opt::seed);
   
   // load the reference
   findex = fai_load(opt::refgenome.c_str());  // load the reference

  // load the blacklist
  SeqLib::GRC blacklist;
  if (!opt::blacklist.empty()) {
    std::cerr << "...reading blacklist file " << opt::blacklist << std::endl;
    blacklist = SeqLib::GRC(opt::blacklist, bwalker.Header());
    blacklist.CreateTreeMap();
    std::cerr << "...read in " << blacklist.size() << " blacklist regions " << std::endl;
  }

  // open the bam file for the header
  SeqLib::BamReader rr;
  if (!rr.Open(opt::bam)) {
    std::cerr << "could not open BAM " << opt::bam << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cerr << "...opening output" << std::endl;
  std::ofstream outfasta;
  fopen(opt::string_id + ".contigs.fa", outfasta);
    
  std::ofstream events;
  fopen(opt::string_id + ".events.txt", events);
  
  PowerLawSim(rr.Header(), findex, opt::nbreaks, -opt::plaw, blacklist, outfasta, events);
  outfasta.close();
  events.close();

}

void parseSimulationOptions(int argc, char** argv) {
  
  bool die = false;
  
  if (argc < 2) 
    die = true;

  std::string t;

  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': die = true; break;
    case 'G': arg >> opt::refgenome; break;
    case 's': arg >> opt::seed; break;
    case 'b': arg >> opt::bam; break;
    case 'R': arg >> opt::nbreaks; break;
    case 'X': arg >> opt::nindels; break;
    case 'a': arg >> opt::string_id; break;
    case 'l': arg >> opt::plaw; break;
    case OPT_BLACKLIST: arg >> opt::blacklist; break;
    case 'M': arg >> opt::viral_count; break;
    default: die= true; 
    }
  }

  if (die) {
    std::cerr << "\n" << SIMULATE_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

}
