#include "realigntest.h"

 #include <getopt.h>
 #include <string>
 #include <sstream>
 #include <iostream>
 #include <cstdlib>

#include "PowerLawSim.h"
#include "SeqLib/RefGenome.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamHeader.h"

 namespace opt {

   static std::string refgenome = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
   static std::string bam;
   static uint32_t seed = 0;
   static int num = 1000;
   static int min = 30;
   static int max = 500;
   static bool evaluate = false;
   static int error = 5;
 }

 static const char *REALIGN_USAGE_MESSAGE =
 "Usage: svlib realigntest -G <ref> -b <bam> -n <num> > sim.fasta \n\n"
 "  Description: Simulate contigs that span breakpoints, realign, and evaluate performance\n"
 "\n" 
 "  Mode:\n"
 "  -E, --evaluate                       Take a BAM of contig alignments and evaluate performance\n"
 "  General options\n"
 "  -G, --reference-genome               Indexed ref genome for BWA-MEM. Default (Broad): /seq/reference/...)\n"
 "  -s, --seed                           Seed for the random number generator [0]\n"
 "  -b, --bam                            BAM file with header to define the genome size\n"
 "  -n, --num-rearrangements             Number of contigs to simulate [1000]\n"
 "  -m, --min-size                       Minimum size of a contig fragment [30]\n"
 "  -M, --max-size                       Maximum size of a contig fragment [500]\n"
 "  -e, --error-num                      Number of SNP errors to place on each contig [1]\n"
 "\n";

static const char* shortopts = "hG:s:n:b:M:m:Ee:";
static const struct option longopts[] = {
  { "help",                 no_argument, NULL, 'h' },
  { "evaluate",             no_argument, NULL, 'E' },
  { "reference-genome",     required_argument, NULL, 'G' },
  { "seed",                 required_argument, NULL, 's' },
  { "bam",                  required_argument, NULL, 'b' },
  { "num-rearrangements",   required_argument, NULL, 'n' },
  { "min-size",             required_argument, NULL, 'm' },
  { "max-size",             required_argument, NULL, 'M' },
  { "error-num",            required_argument, NULL, 'e' },
  { NULL, 0, NULL, 0 }
};

void runRealignTest(int argc, char** argv) {
  
  parseRealignOptions(argc, argv);

  if (opt::evaluate) {
    evaluateRealignment();
    return;
  }
  
  // seed the RNG
  srand(opt::seed);

  SeqLib::RefGenome ref;
  if (!ref.LoadIndex(opt::refgenome)) {
    std::cerr << "Could not load reference genome: " << opt::refgenome << std::endl;
    exit(EXIT_FAILURE);
  }

  // read a bam to get the headers
  SeqLib::BamReader rr;
  if (!rr.Open(opt::bam)) {
    std::cerr << "Could not load BAM file: " << opt::bam << std::endl;
    exit(EXIT_FAILURE);
  }
  const SeqLib::BamHeader hdr = rr.Header();

  // set the weights for choosing random chromosome
  const std::vector<double> CHR_CUMSUM_WEIGHT_X = {0.08209014, 0.16218732, 0.22740558, 0.29036182, 0.34994586, 0.40630223, 0.45871420,
					     0.50691887, 0.55342720, 0.59806527, 0.64252937, 0.68661320, 0.72454415, 0.75989948,
					     0.79366797, 0.82342611, 0.85016757, 0.87588214, 0.89535614, 0.91611346, 0.93196494,
					     0.94886198, 1.00000000};

  const char TCGA[5] = "TCGA";
    
  int diff = opt::max - opt::min;
  assert(diff >= 0);
  size_t count = 0;
  // choose a set of random breaks
  for (int i = 0; i < opt::num; ++i) {
    SeqLib::GenomicRegion gr;
    gr.chr = weightedRandom(CHR_CUMSUM_WEIGHT_X);
    gr.pos1 = 1e6 + rand() % (int)(hdr.GetSequenceLength(gr.chr) - 2e6);
    gr.pos2 = gr.pos1 + (rand() % diff) + opt::min;
    gr.strand = '+';
    std::string seq1 = ref.QueryRegion(hdr.IDtoName(gr.chr), gr.pos1, gr.pos2);
    if (seq1.empty() || seq1.find("N") != std::string::npos)
      continue;
    if (rand() % 2) {
      gr.strand = '-';
      SeqLib::rcomplement(seq1);
    }
      
    SeqLib::GenomicRegion gr2;
    gr2.chr = weightedRandom(CHR_CUMSUM_WEIGHT_X);
    gr2.pos1 = 1e6 + rand() % (int)(hdr.GetSequenceLength(gr2.chr) - 2e6);
    gr2.pos2 = gr2.pos1 + (rand() % diff) + opt::min;
    gr2.strand = '+';
    std::string seq2 = ref.QueryRegion(hdr.IDtoName(gr2.chr), gr2.pos1, gr2.pos2);
    if (seq2.empty() || seq2.find("N") != std::string::npos)
      continue;
    if (rand() % 2) {
      gr2.strand = '-';
      SeqLib::rcomplement(seq2);
    }

    std::string seq = seq1 + seq2;
    for (int i = 0; i < opt::error; ++i) {
      int pos = rand() % seq.length();
      seq[pos] = TCGA[rand() % 4];
    }
    
    ++count;
    std::stringstream name;
    name << ">" << gr.ChrName(hdr) << "_" << gr.pos1 << "_" << gr.pos2 << "_" << gr.strand << "_" 
	 << gr2.ChrName(hdr) << "_" << gr2.pos1 << "_" << gr2.pos2 << "_" << gr2.strand; 
    std::cout << name.str() << std::endl << seq << std::endl;
  }

}

static void split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss;
  ss.str(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
}

void evaluateRealignment() {

  // read in the bam
  // read a bam to get the headers
  SeqLib::BamReader rr;
  if (!rr.Open(opt::bam)) {
    std::cerr << "Could not load BAM file: " << opt::bam << std::endl;
    exit(EXIT_FAILURE);
  }
  const SeqLib::BamHeader hdr = rr.Header();

  // group by contigs
  std::unordered_map<std::string, SeqLib::BamRecordVector> map;
  SeqLib::BamRecord r;
  while (rr.GetNextRecord(r)) 
    map[r.Qname()].push_back(r);
  
  // evaluate
  for (const auto& c : map) {

    // parse the name
    std::vector<std::string> tokens;
    split(c.first, '_', tokens);
    assert(tokens.size() == 8);
    SeqLib::GenomicRegion gr1(tokens[0], tokens[1], tokens[2], hdr);
    gr1.strand = tokens[3].at(0);
    SeqLib::GenomicRegion gr2(tokens[4], tokens[5], tokens[6], hdr);
    gr2.strand = tokens[7].at(0);

    bool hit1 = false;
    bool hit2 = false;
    int32_t max_mapq = 0;
    int32_t min_mapq = 10000;
    for (const auto& f : c.second) {
      if (gr1.GetOverlap(f.AsGenomicRegion()))
	hit1 = true;
      if (gr2.GetOverlap(f.AsGenomicRegion()))
	hit2 = true;
      max_mapq = std::max(f.MapQuality(), max_mapq);
      min_mapq = std::min(f.MapQuality(), min_mapq);
    }

    std::cout << gr1.ChrName(hdr) << "\t" << gr1.pos1 << "\t" << gr1.pos2 << "\t" << gr1.strand << "\t" 
	      << gr2.ChrName(hdr) << "\t" << gr2.pos1 << "\t" << gr2.pos2 << "\t" << gr2.strand << "\t" 
	      << hit1 << "\t" << hit2 << "\t" << c.second.size() << "\t" << gr1.Width() << "\t" << gr2.Width()
	      << "\t" << max_mapq << "\t" << min_mapq << std::endl; 

  }

}

void parseRealignOptions(int argc, char** argv) {
  
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
    case 'n': arg >> opt::num; break;
    case 'm': arg >> opt::min; break;
    case 'M': arg >> opt::max; break;
    case 'e': arg >> opt::error; break;
    case 'E': opt::evaluate = true;; break;
    default: die= true; 
    }
  }

  if (die) {
    std::cerr << "\n" << REALIGN_USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

}
