#include "seq2vcf.h"

#include <string>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <cassert>
#include <map>
#include <iomanip>

#include "AlignedContig.h"
#include "PlottedRead.h"

// trim this many bases from front and back of read when determining coverage
// this should be synced with the split-read buffer in BreakPoint2 for more accurate 
// representation of covearge of INFORMATIVE reads (eg ones that could be split)
#define INFORMATIVE_COVERAGE_BUFFER 4

namespace opt {
  static std::string bam;
  static std::map<std::string, std::string> short_bams;
  static int pad = 0;
  static int cores = 1;
}

static int sample_number = 0;

static const char* shortopts = "hp:P:n:t:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "case-bam",                required_argument, NULL, 't' },
  { "control-bam",             required_argument, NULL, 'n' },
  { "pad",                     required_argument, NULL, 'P' },
  { "cores",                   required_argument, NULL, 'p' },
  { NULL, 0, NULL, 0 }
};

static const char *RUN_SEQTOVCF_MESSAGE =
"Usage: svlib seqtovcf <BAM> [OPTION]\n\n"
"  Description: Produce SV calls directly from alignment of long-reads to a genome\n"
"\n"
"  General options\n"
"  -h, --help                          Display this help and exit\n"
"  -p, --cores                         Number of cores to run [1]\n"
"  -P, --pad                           Distance to look outside of contig region for supporting reads [100]\n"
"\n";

std::string bamOptParse(std::map<std::string, std::string>& obam, std::istringstream& arg, int sample_number, const std::string& prefix) {
  std::stringstream ss;
  std::string bam;
  arg >> bam;
  ss.str(std::string());
  ss << prefix << std::setw(3) << std::setfill('0') << sample_number;
  obam[ss.str()] = bam;
  return bam;
}


void runSeqToVCF(int argc, char** argv) {

  parseSeqToVCFOptions(argc, argv);

  // open the input
  SeqLib::BamReader reader;
  if (!reader.Open(opt::bam)) {
    std::cerr << "could not open BAM: " << opt::bam << std::endl;
    exit(EXIT_FAILURE);
  }
  SeqLib::BamHeader hdr = reader.Header();

  // Create the queue and consumer (worker) threads
  WorkQueue<LongReaderWorkItem*>  queue; // queue of work items to be processed by threads
  std::vector<ConsumerThread<LongReaderWorkItem, LongReaderThreadItem>*> threadqueue;

  // loop through and add jobs to the queue
  size_t count = 0;
  std::string curr_name;
  int max_mapq = 0, max_len = 0;
  SeqLib::BamRecordVector brv;
  SeqLib::GRC regions;
  SeqLib::BamRecord r;
  while(reader.GetNextRecord(r)) {
    
    if (++count % 200000 == 0) 
      std::cerr << "...reading in contig " << SeqLib::AddCommas(count) << std::endl;

    // first one
    if (curr_name.empty())
      curr_name = r.Qname();
    
    // add the MC tag 
    if (r.ChrID() < hdr.NumSequences() && r.ChrID() >= 0) 
      r.AddZTag("MC", hdr.IDtoName(r.ChrID()));
    else 
      continue;

    // add to existing
    if (curr_name == r.Qname())
      add_contig(r, regions, brv, max_mapq, max_len);

    // we found all of the contig alignments for this qname
    // start processing it
    else {
      curr_name = r.Qname();
      assert(brv.size());
      queue.add(new LongReaderWorkItem(ContigElement(brv, regions)));
      brv.clear();
      regions.clear();

      // start the next one
      add_contig(r, regions, brv, max_mapq, max_len);
    }

  }

  std::cerr << "...starting sv parsing threads" << std::endl;
  // establish and start the threads
  for (int i = 0; i < opt::cores; i++) {
    LongReaderThreadItem * ti =  new LongReaderThreadItem(i, opt::short_bams); // create the thread-specific data
    ConsumerThread<LongReaderWorkItem, LongReaderThreadItem>* threadr = // establish new thread to draw from queue
      new ConsumerThread<LongReaderWorkItem, LongReaderThreadItem>(queue, ti);
    threadr->start(); 
    threadqueue.push_back(threadr); // add thread to the threadqueue
  }
  
  // wait for the threads to finish
  for (int i = 0; i < opt::cores; ++i)  {
    threadqueue[i]->join();
    std::cout << threadqueue[i]->GetThreadData()->bps.str() << std::endl;
  }
}

// process an individual contig
bool LongReaderWorkItem::runLR(LongReaderThreadItem* thread_item) {

  std::vector<AlignedContig> this_alc;

  // contruct the BWA index of the long seq, for read realignment
  SeqLib::UnalignedSequenceVector usv;
  assert(m_contig.brv[0].Qname().length());
  //assert(m_contig.brv[0].Sequence().find("N") == std::string::npos);

  // here we have to flip if it has (-) alignment to reference.
  // this is because the pipeline assumes it came from a 
  // de novo assembly, which is PRE alignment. Thus, we don't want 
  // to take the i.Sequeunce directly, as this has been reverse-complemented
  // BY THE ALIGNER.
  std::string sss = m_contig.brv[0].Sequence();
  if (m_contig.brv[0].ReverseFlag())
    SeqLib::rcomplement(sss);
  usv.push_back({m_contig.brv[0].Qname(), sss, std::string()});

  // make the index
  SeqLib::BWAWrapper bw;
  bw.ConstructIndex(usv);

  // add the prefixes
  std::set<std::string> prefixes;
  for (const auto& b : opt::short_bams)
    prefixes.insert(b.first);

  // setup coverages
  std::map<std::string, STCoverage> covs;
  for (const auto& c : thread_item->readers)
    covs.insert(std::pair<std::string, STCoverage>(c.first, STCoverage()));

  // grab the reads
  SeqLib::BamRecordVector reads;
  for (auto& reader : thread_item->readers) {
    grab_reads(reader.first, reader.second, m_contig.regions, reads,
	       covs[reader.first]);
  }
    

  // make the actual AligneContig
  this_alc.push_back(AlignedContig(m_contig.brv, prefixes));
  assert(this_alc.size() == 1);
  AlignedContig * ac = &this_alc.back();
  //for (auto& kk : ac->m_frag_v) // set max indel for all of the frags
  //  kk.m_max_indel = 20;

  // false says dont worry about "local"
  std::vector<BreakPoint> allbreaks = ac->getAllBreakPoints(false); 
  for (auto& i : allbreaks)
    i.repeatFilter();
  //for (auto& i : allbreaks)
  //  i.addCovs(covs, clip_covs);
  for (auto& i : allbreaks)
    i.scoreBreakpoint(8, 2.5, 7, 3, 1, 100);
  //for (auto& i : allbreaks)
  //  i.setRefAlt(ref_genomeAW, ref_genome_viral_dummy);


  // write the breakpoints to this thread
  bool no_reads = true;
  for (auto& i : allbreaks)
    thread_item->bps << i.toFileString(no_reads) << std::endl;

  return true;
}

void parseSeqToVCFOptions(int argc, char** argv) {

  bool die = false;

  if (argc < 2) 
    die = true;
  else
    opt::bam = std::string(argv[1]);

  bool help = false;
  std::stringstream ss;

  std::string tmp;
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'h': help = true; break;
    case 'p': arg >> opt::cores; break;
    case 'P': arg >> opt::pad; break;
    case 'n': 
      tmp = bamOptParse(opt::short_bams, arg, sample_number++, "n"); break;
    case 't':  
      tmp = bamOptParse(opt::short_bams, arg, sample_number++, "t"); break;

    }
  }

  if (die || help) {
    std::cerr << "\n" << RUN_SEQTOVCF_MESSAGE;
    die ? exit(EXIT_FAILURE) : exit(EXIT_SUCCESS);
  }

}

void add_contig(const SeqLib::BamRecord& r, SeqLib::GRC& regions, SeqLib::BamRecordVector& brv,
		int& max_mapq, int& max_len) {

  SeqLib::GenomicRegion gr = r.AsGenomicRegion(); // find the region covered by contig
  gr.Pad(opt::pad); // add to it
  regions.add(gr);
  brv.push_back(r);
  max_mapq = std::max(r.MapQuality(), max_mapq);
  max_len = std::max(r.Length(), max_len);
  
}

void grab_reads(const std::string& prefix, SeqLib::BamReader& reader, const SeqLib::GRC& regions, SeqLib::BamRecordVector& brv,
		STCoverage& cov) {

  if (!reader.SetMultipleRegions(regions)) {
    std::cerr << "...failed to set regions " << std::endl;
    for (const auto& rr : regions)
      std::cerr << "\t" << rr << std::endl;
  }

  const size_t hard_limit = 1e6;
  SeqLib::BamRecord r;
  size_t count = 0;
  while (reader.GetNextRecord(r) && ++count < hard_limit) {

    // dont even mess with them
    if (r.CountNBases())
      continue;

    cov.addRead(r, INFORMATIVE_COVERAGE_BUFFER, false); 
    
    std::string srn =  prefix + "_" + std::to_string(r.AlignmentFlag()) + "_" + r.Qname();
    assert(srn.length());
    r.AddZTag("SR", srn);

    // for memory conservation
    r.RemoveTag("BQ");
    r.RemoveTag("OQ");

    brv.push_back(r);
  }
  
}
