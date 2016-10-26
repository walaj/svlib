#include "vcf2bedpe.h"

#include "gzstream.h"
#include <string>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <cassert>
#include <unordered_map>

#include "SeqLib/GenomicRegion.h"

namespace opt {
  std::string input;
}

static const char* shortopts = "h";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { NULL, 0, NULL, 0 }
};

static const char *RUN_VCF2BEDPE_MESSAGE =
"Usage: svlib vcftobedpe <VCF> [OPTION]\n\n"
"  Description: Produce a BEDPE from a breakend or DEL/TRA/DUP/INV SV VCF\n"
"\n"
"  General options\n"
"  -h, --help                          Display this help and exit\n"
"\n";

// define a class to hold break and info
class BreakRegion : public SeqLib::GenomicRegion {

public:

  // the 1 is a dummy, store name in chr_name
  BreakRegion(const std::string c, const std::string p) : 
    SeqLib::GenomicRegion("1", p, p, SeqLib::BamHeader()), chr_name(c) {}

  std::string chr_name;
};

void runVCFToBEDPE(int argc, char** argv) {
  
  parseVCFToBEDPEOptions(argc, argv);

  // open the file
  igzstream vcf(opt::input.c_str(), std::ios::in);
  if (!vcf) {
    std::cerr << "Cannot open VCF " << opt::input << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;

  // hold the variants if in BND format
  std::unordered_map<std::string, std::vector<BreakRegion>> bks;

  while (std::getline(vcf, line, '\n')) {
    if (line.find("#") != std::string::npos) // skip comment lines
      continue;

    // breakend format
    if (line.find("BND")) {
      std::istringstream iss_r(line);
      std::string val;
      size_t count = 0;
      std::string pos, chr, id, ref, alt, qual, filter, info;
      std::vector<std::string> geno; 
      while (std::getline(iss_r, val, '\t')) {
	switch(count) {
	case 0: chr = val; break;
	case 1: pos = val; break;
	case 2: id = val; break;
        case 3: ref = val; break;
	case 4: alt = val; break;
        case 5: qual = val; break;
        case 6: filter = val; break;
	case 7: info = val; break;
	default: geno.push_back(val); break;
	}
	++count;
      }

      // make the break region
      BreakRegion br(chr, pos);
      br.chr_name = chr;
      
      // get the id
      std::string token_id = id.substr(0, id.find(":"));
      
      // insert the region for this id
      bks[token_id].push_back(br);
      assert(bks[token_id].size() <= 2);
    }

    // DUP etc format
    else {

    }

  }

  // assure they are all there
  for (auto& b : bks) {
    if (b.second.size() != 2) {
      std::cerr << "ERROR on id: " << b.first << " expecting 2 breaks, found " << b.second.size() << std::endl;
      continue;
    }

    BreakRegion *br1 = &b.second[0];
    BreakRegion *br2 = &b.second[1];
    if (*br2 < *br1)
      std::swap(*br1, *br2);
    
    // write it out
    std::cout << br1->chr_name << "\t" << br1->pos1 << "\t" << br1->pos1 << "\t" 
	      << br2->chr_name << "\t" << br2->pos1 << "\t" << br2->pos1 << std::endl;
  }
   

}

void parseVCFToBEDPEOptions(int argc, char** argv) {

  bool die = false;

  if (argc < 2) 
    die = true;
  else
    opt::input = std::string(argv[1]);

  bool help = false;
  std::stringstream ss;

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': help = true; break;
    }
  }

  if (die || help) {
    std::cerr << "\n" << RUN_VCF2BEDPE_MESSAGE;
    die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
  }

}
