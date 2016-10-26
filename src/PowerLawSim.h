#ifndef SVLIB_POWER_LAW_SIM_H
#define SVLIB_POWER_LAW_SIM_H

#include <fstream>
#include <vector>

#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BamHeader.h"

void PowerLawSim(const SeqLib::BamHeader& hdr, faidx_t* findex, int num_breaks, double power_law, SeqLib::GRC& grc, 
		 std::ofstream& outfasta, std::ofstream& events);
int weightedRandom(const std::vector<double>& cs);
std::vector<int> drawFromPower(double x0, double x1, double power, int n_draws);
void genRandomSequence(const SeqLib::BamHeader& hdr, std::string& s, SeqLib::GenomicRegion& gr, int width, faidx_t * findex, SeqLib::GRC& grc);

struct SVEvent {
  
  SeqLib::GenomicRegion reg1;
  SeqLib::GenomicRegion reg2;
  
  int break1;
  int break2;

  std::string ins;
  std::string r_ins;

  std::string etype;

  int number;
};

#endif
