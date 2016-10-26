#include "seq2vcf.h"

#include <string>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <cassert>
#include <map>
#include <iomanip>

#include "vcf.h"
#include "gzstream.h"

#include "AlignedContig.h"
#include "PlottedRead.h"
#include "alignReadsToContigs.h"

// trim this many bases from front and back of read when determining coverage
// this should be synced with the split-read buffer in BreakPoint2 for more accurate 
// representation of covearge of INFORMATIVE reads (eg ones that could be split)
#define INFORMATIVE_COVERAGE_BUFFER 4

static pthread_mutex_t lock;

namespace opt {
  static bool no_unfiltered = false;
  static int max_coverage = 100;
  static std::string refgenome;
  static std::string bam;
  static std::string analysis_id = "s2v";
  static std::map<std::string, std::string> short_bams;
  static int pad = 0;
  static int cores = 1;
}

static std::string args = "svlib ";
static int sample_number = 0;

static SeqLib::BamHeader hdr;

static const char* shortopts = "hp:P:n:t:G:a:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "case-bam",                required_argument, NULL, 't' },
  { "control-bam",             required_argument, NULL, 'n' },
  { "pad",                     required_argument, NULL, 'P' },
  { "cores",                   required_argument, NULL, 'p' },
  { "reference-genome",        required_argument, NULL, 'G' },
  { "analysis-id",             required_argument, NULL, 'a' },
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

  // open the mutex
  if (pthread_mutex_init(&lock, NULL) != 0) {
    std::cerr << "\n mutex init failed\n";
    exit(EXIT_FAILURE);
  }

  // open the input
  SeqLib::BamReader reader;
  if (!reader.Open(opt::bam)) {
    std::cerr << "could not open BAM: " << opt::bam << std::endl;
    exit(EXIT_FAILURE);
  }

  hdr = reader.Header();

  // open the bps
  ogzstream bps_file, aln_file;
  std::string bps_name = opt::analysis_id + ".bps.txt.gz";
  std::string aln_name = opt::analysis_id + ".alignments.txt.gz";
  bps_file.open(bps_name.c_str(), std::ios::out);
  aln_file.open(aln_name.c_str(), std::ios::out);

  // Create the queue and consumer (worker) threads
  WorkQueue<LongReaderWorkItem*>  queue; // queue of work items to be processed by threads
  std::vector<ConsumerThread<LongReaderWorkItem, LongReaderThreadItem>*> threadqueue;

  std::cerr << "...starting sv parsing threads" << std::endl;
  // establish and start the threads
  for (int i = 0; i < opt::cores; i++) {
    LongReaderThreadItem * ti =  new LongReaderThreadItem(i, opt::short_bams); // create the thread-specific data
    ti->bps_file = &bps_file;
    ti->aln_file = &aln_file;
    if (!ti->LoadReference(opt::refgenome)) {
      std::cerr << "ERROR could not load reference genome " << opt::refgenome << std::endl;
      exit(EXIT_FAILURE);
    }
      
    ConsumerThread<LongReaderWorkItem, LongReaderThreadItem>* threadr = // establish new thread to draw from queue
      new ConsumerThread<LongReaderWorkItem, LongReaderThreadItem>(queue, ti);
    threadr->start(); 
    threadqueue.push_back(threadr); // add thread to the threadqueue
  }

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
      // only add if multi-part alignment or has gapped alignment
      if (brv.size() > 1 || brv[0].MaxDeletionBases() || brv[0].MaxInsertionBases()) { 
	LongReaderWorkItem * item = new LongReaderWorkItem(ContigElement(brv, regions));
	queue.add(item);
      }
      brv.clear();
      regions = SeqLib::GRC(); // create a new one

      // start the next one
      add_contig(r, regions, brv, max_mapq, max_len);
    }

  }
    

  // wait for the threads to finish
  for (int i = 0; i < opt::cores; ++i)  {
    threadqueue[i]->join();
  }

  // close the text files
  aln_file.close();
  bps_file.close();

  // put args into string for VCF later
  for (int i = 0; i < argc; ++i)
    args += std::string(argv[i]) + " ";

  // make the vcf header
  VCFHeader header;
  header.filedate = fileDateString();
  header.source = args;
  header.reference = opt::refgenome;

  if (!hdr.isEmpty())
    for (int i = 0; i < hdr.NumSequences(); ++i)
      header.addContigField(hdr.IDtoName(i),hdr.GetSequenceLength(i));

  for (auto& b : opt::short_bams) {
    std::string fname = b.second; //bpf.filename();
    header.addSampleField(fname);
    header.colnames += "\t" + fname; 
  }

  // primary VCFs
  if (SeqLib::read_access_test(bps_name)) {
    VCFFile vcf(bps_name, opt::analysis_id, hdr, header, !opt::no_unfiltered);
    
    if (!opt::no_unfiltered) {
      std::string basename = opt::analysis_id + ".seqtovcf.unfiltered.";
      vcf.include_nonpass = true;
      std::cerr << "...writing unfiltered VCFs" << std::endl;
      vcf.writeIndels(basename, false, opt::short_bams.size() == 1);
      vcf.writeSVs(basename, false, opt::short_bams.size() == 1);
    }

    std::cerr << "...writing filtered VCFs" << std::endl;
    std::string basename = opt::analysis_id + ".seqtovcf.";
    vcf.include_nonpass = false;
    vcf.writeIndels(basename, false, opt::short_bams.size() == 1);
    vcf.writeSVs(basename, false, opt::short_bams.size() == 1);

  } else {
    std::cerr << "ERROR: Failed to make VCF. Could not file bps file " + bps_name << std::endl;
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
  bw.SetGapOpen(16); // default 6
  bw.SetMismatchPenalty(9); // default 4
  bw.ConstructIndex(usv);

  // add the prefixes
  std::set<std::string> prefixes;
  for (const auto& b : opt::short_bams)
    prefixes.insert(b.first);

  // setup coverages
  std::unordered_map<std::string, STCoverage*> covs;
  for (const auto& c : thread_item->readers)
    covs.insert(std::pair<std::string, STCoverage*>(c.first, new STCoverage()));

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

  // align reads to contigs
  alignReadsToContigs(hdr, bw, usv, reads, *ac, &thread_item->ref);

  // score read alignments
  // Get contig coverage, discordant matching to contigs, etc
  // repeat sequence filter
  ac->assessRepeats();
  
  ac->splitCoverage();
  // now that we have all the break support, check that the complex breaks are OK
  ac->refilterComplex(); 
  // add discordant reads support to each of the breakpoints
  //ac->addDiscordantCluster(dmap);
  // add in the cigar matches
  //ac->checkAgainstCigarMatches(cigmap);

  // false says dont worry about "local"
  std::vector<BreakPoint> allbreaks = ac->getAllBreakPoints(false); 
  for (auto& i : allbreaks)
    i.repeatFilter();
  for (auto& i : allbreaks)
    i.addCovs(covs);
  for (auto& i : allbreaks)
    i.scoreBreakpoint(8, 2.5, 7, 3, 1, 100);
  for (auto& i : allbreaks)
    i.setRefAlt(&thread_item->ref, NULL);

  // write the breakpoints to this thread
  bool no_reads = true;
  for (auto& i : allbreaks)
    thread_item->bps << i.toFileString(no_reads) << std::endl;

  //write the alignments to file
  thread_item->aln << *ac << std::endl;

  // clear coverage data
  for (auto& c : covs)
    delete c.second;

  // print info
  ++thread_item->num_processed;
  if (thread_item->num_processed % 1000 == 0) {
    std::cerr << "...processed " << SeqLib::AddCommas(thread_item->num_processed) << " seqs on thread " 
	      << thread_item->thread_id << std::endl;
    thread_item->WriteFiles(&lock);
  }
  
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
    case 'G': arg >> opt::refgenome; break;
    case 'a': arg >> opt::analysis_id; break;
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
		STCoverage* cov) {

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

    cov->addRead(r, INFORMATIVE_COVERAGE_BUFFER, false); 
    
    std::string srn =  prefix + "_" + std::to_string(r.AlignmentFlag()) + "_" + r.Qname();
    assert(srn.length());
    r.AddZTag("SR", srn);

    // for memory conservation
    r.RemoveTag("BQ");
    r.RemoveTag("OQ");

    brv.push_back(r);
  }

  SeqLib::BamRecordVector new_reads;

  const int m_seed = 6;

  for (auto& r : brv) {
      double this_cov1 = cov->getCoverageAtPosition(r.ChrID(), r.Position());
      double this_cov2 = cov->getCoverageAtPosition(r.ChrID(), r.PositionEnd());
      double this_cov = std::max(this_cov1, this_cov2);
      double sample_rate = 1; // dummy, always set if max_coverage > 0
      if (this_cov > 0) 
	sample_rate = 1 - (this_cov - opt::max_coverage) / this_cov; // if cov->inf, sample_rate -> 0. if cov -> max_cov, sample_rate -> 1
      
      // this read should be randomly sampled, cov is too high
      if (this_cov > opt::max_coverage) {
#ifdef QNAME
	  if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1)) {
	    std::cerr << "subsampling because this_cov is " << this_cov << " and max cov is " << opt::max_coverage << " at position " << r.Position() << " and end position " << r.PositionEnd() << std::endl;
	    std::cerr << " this cov 1 " << this_cov1 << " this_cov2 " << this_cov2 << std::endl;
	  }
#endif
	  uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(r.Qname().c_str()) ^ m_seed);
	  if ((double)(k&0xffffff) / 0x1000000 <= sample_rate) // passed the random filter
	    new_reads.push_back(r);
      }
      else // didn't have a coverage problems
	{
	  new_reads.push_back(r);
	}
      
  }
  
  brv = new_reads;
  
}

std::string fileDateString() {
  // set the time string
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  std::stringstream month;
  std::stringstream mdate;
  if ( (now->tm_mon+1) < 10)
    month << "0" << now->tm_mon+1;
  else 
    month << now->tm_mon+1;
  mdate << (now->tm_year + 1900) << month.str() <<  now->tm_mday;
  return mdate.str();
}
