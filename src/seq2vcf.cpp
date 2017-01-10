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

SeqLib::BamWriter contig_writer;

static SeqLib::GRC black;
namespace opt {
  static bool read_tracking = false;
  static bool no_unfiltered = false;
  static int max_coverage = 1000;
  static std::string refgenome;
  static std::string bam;
  static std::string analysis_id = "s2v";
  static std::map<std::string, std::string> short_bams;
  static int pad = 0;
  static int cores = 1;
  static int readlen = 0;
  static std::string blacklist;
  std::unordered_set<std::string> contig_list;
  std::string list;
  bool write_extracted_contigs = false;
}

static std::string args = "svlib ";
static int sample_number = 0;

static SeqLib::BamHeader hdr;

static const char* shortopts = "hp:P:n:t:G:a:m:B:L:W";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "case-bam",                required_argument, NULL, 't' },
  { "contig-list",             required_argument, NULL, 'L' },
  { "control-bam",             required_argument, NULL, 'n' },
  { "pad",                     required_argument, NULL, 'P' },
  { "cores",                   required_argument, NULL, 'p' },
  { "reference-genome",        required_argument, NULL, 'G' },
  { "max-coverage",            required_argument, NULL, 'm' },
  { "analysis-id",             required_argument, NULL, 'a' },
  { "write-extracted-contigs", no_argument, NULL, 'W' },
  { "blacklist",               required_argument, NULL, 'B' },
  { NULL, 0, NULL, 0 }
};

static const char *RUN_SEQTOVCF_MESSAGE =
"Usage: svlib seqtovcf <BAM> -G <REF> [OPTION]\n\n"
"  Description: Produce SV calls directly from alignment of long-reads to a genome\n"
"  NOTE: The input BAM must be qname sorted\n"
"        The optinal input short-read BAMs must be coordinate sorted and indexed\n"
"\n"
"  Required input\n"
"  -G, --reference-genome              An indexed reference genome\n"
"  General options\n"
"  -h, --help                          Display this help and exit\n"
"  -p, --cores                         Number of cores to run [1]\n"
"  -t, --case-bam                      A case (eg tumor) BAM of short reads for scoring / genotyping. Accepts multiple\n"
"  -n, --control-bam                   A control (eg normal) BAM of short reads for scoring / genotyping. Accepts multiple\n"
"  -P, --pad                           Distance to look outside of contig region for supporting reads [100]\n"
"  -m, --max-coveage                   Downsample reads to this coverage when extracting for scoring/genotyping [100]\n"
"  -a, --analysis-id                   A unique analysis id to prepend output files with [s2v]\n"
"  -B, --blacklist                     BED file of regions not to consider further, e.g. centromeric\n"
"  -L, --contig-list                   List of sequence names to only include. Names not here will not be processes.\n"
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

  // setup a way to write which contigs were extracted
  if (opt::write_extracted_contigs) {
    std::string obam = opt::analysis_id + ".extracted.contigs.bam";
    std::cerr << "...will write extracted contigs to " << obam << std::endl;
    contig_writer = SeqLib::BamWriter(); 
    contig_writer.SetHeader(hdr);
    contig_writer.Open(obam);
    contig_writer.WriteHeader();
  }
  
  // open blacklist
  if (!opt::blacklist.empty()) {
    black.ReadBED(opt::blacklist, hdr);
    black.CreateTreeMap();
  }

  // open the contig list
  if (!opt::list.empty()) {
    std::ifstream fs(opt::list);
    std::string cc;
    while(fs >> cc)
      opt::contig_list.insert(cc);
    std::cerr << "...read in " << SeqLib::AddCommas(opt::contig_list.size()) << " contig names to keep " << std::endl;
  }

  // open the bps
  ogzstream bps_file, aln_file;
  std::string bps_name = opt::analysis_id + ".bps.txt.gz";
  std::string aln_name = opt::analysis_id + ".alignments.txt.gz";
  bps_file.open(bps_name.c_str(), std::ios::out);

  // give header
  bps_file << BreakPoint::header() << std::endl;

  //aln_file.open(aln_name.c_str(), std::ios::out);

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
    
    if (opt::contig_list.size() && !opt::contig_list.count(r.Qname()))
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
      assign_contig(brv, queue, regions);

      brv.clear();
      regions = SeqLib::GRC(); // create a new one

      // start the next one
      add_contig(r, regions, brv, max_mapq, max_len);
    }

  }

  // assign the last one
  assign_contig(brv, queue, regions);

  // close the bam
  if (contig_writer.IsOpen())
    contig_writer.Close();
    
  std::cerr << "...processing " << SeqLib::AddCommas(queue.size()) << " contigs" << std::endl;
  opt::cores = opt::cores > queue.size() ? queue.size() : opt::cores; 

  std::cerr << "...starting sv parsing threads" << std::endl;
  // establish and start the threads
  for (int i = 0; i < opt::cores; i++) {
    LongReaderThreadItem * ti =  new LongReaderThreadItem(i, opt::short_bams); // create the thread-specific data
    ti->bps_file = &bps_file;
    //ti->aln_file = &aln_file;
    if (!ti->LoadReference(opt::refgenome)) {
      std::cerr << "ERROR could not load reference genome " << opt::refgenome << std::endl;
      exit(EXIT_FAILURE);
    }
      
    ConsumerThread<LongReaderWorkItem, LongReaderThreadItem>* threadr = // establish new thread to draw from queue
      new ConsumerThread<LongReaderWorkItem, LongReaderThreadItem>(queue, ti);
    threadr->start(); 
    threadqueue.push_back(threadr); // add thread to the threadqueue
  }

  // wait for the threads to finish
  for (int i = 0; i < opt::cores; ++i)  {
    threadqueue[i]->join();
  }

  // write whatever remains
  for (const auto& c : threadqueue) 
    c->GetThreadData()->WriteFiles(&lock);

  //std::vector<BreakPoint> bpall;
  //for (const auto& c : threadqueue)
  //  bpall.insert(bpall.end(), c->GetThreadData()->bp_glob.begin(), c->GetThreadData()->bp_glob.end());

  //std::cerr << "...parsed " << bpall.size() << " breakpoints " << std::endl;

  // de duplicate the breakpoints
  //std::sort(bpall.begin(), bpall.end());
  //bpall.erase( std::unique( bpall.begin(), bpall.end() ), bpall.end() );

  // send breakpoints to file
  //for (auto& i : bpall) 
  //  bps_file << i.toFileString(!opt::read_tracking) << std::endl;

  // close the text files
  //aln_file.close();
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
  AlignedContig ac(m_contig.brv, prefixes);
  
  // get the breaks, check that one is valid
  std::vector<BreakPoint> tmpbreaks = ac.getAllBreakPoints(false); 
  size_t blist = 0;
  for (const auto& b : tmpbreaks)
    blist += (black.CountOverlaps(b.b1.gr) | black.CountOverlaps(b.b2.gr)) ? 1 : 0;
  if (blist == tmpbreaks.size()) {
    //std::cerr << " skipping breaks " << std::endl;
    //for (const auto& b : tmpbreaks)
    //  std::cerr << "    " << b << std::endl;
    return true;
  }

  // align reads to contigs
  alignReadsToContigs(hdr, bw, usv, reads, ac, &thread_item->ref);

  // score read alignments
  // Get contig coverage, discordant matching to contigs, etc
  // repeat sequence filter
  ac.assessRepeats();
  
  ac.splitCoverage();
  // now that we have all the break support, check that the complex breaks are OK
  ac.refilterComplex(); 
  // add discordant reads support to each of the breakpoints
  //ac.addDiscordantCluster(dmap);
  // add in the cigar matches
  //ac.checkAgainstCigarMatches(cigmap);

  // false says dont worry about "local"
  std::vector<BreakPoint> allbreaks = ac.getAllBreakPoints(false); 

  for (auto& i : allbreaks)
    i.readlen = opt::readlen;
  for (auto& i : allbreaks)
    i.repeatFilter();
  for (auto& i : allbreaks)
    i.addCovs(covs);
  for (auto& i : allbreaks)
    i.scoreBreakpoint(8, 2.5, 7, 3, 1, 100);
  for (auto& i : allbreaks)
    i.setRefAlt(&thread_item->ref, NULL);

  // write the breakpoints to this thread
  //for (auto& b : allbreaks)
  //  b.lite(); // clear out string data to make easier storage

  // de duplicate the breakpoints
  std::sort(allbreaks.begin(), allbreaks.end());
  allbreaks.erase( std::unique( allbreaks.begin(), allbreaks.end() ), allbreaks.end() );

  // add them to the global
  for (const auto& b : allbreaks)
    if (b.hasMinimal() || !opt::short_bams.size())
      thread_item->bp_glob.push_back(b);
  
  //write the alignments to file
  //thread_item->aln << ac << std::endl;

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
    case 'm': arg >> opt::max_coverage; break;
    case 'B': arg >> opt::blacklist; break;
    case 'W': opt::write_extracted_contigs = true; break;
    case 'L': arg >> opt::list; break;
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

    // update the readlength
    opt::readlen = std::max(opt::readlen, (int)r.Length());

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

void assign_contig(const SeqLib::BamRecordVector& brv, WorkQueue<LongReaderWorkItem*>& queue, 
		   const SeqLib::GRC& regions) {
  
  // only add if multi-part alignment or has gapped alignment
  int max_mapq = 0;
  for (const auto& a : brv)
    max_mapq = a.MapQuality() > max_mapq ? a.MapQuality() : max_mapq;
  
  SeqLib::BamRecordVector brv2;
  
  for (const auto& b : brv)
    if (b.NumMatchBases() > 100) //(brv.size() > 1) || brv[0].MaxDeletionBases() >=20  || brv[0].MaxInsertionBases() >= 20)
      brv2.push_back(b);

  //if (brv[0].Qname()=="30893") {
  //  std::cerr << " HERE " << brv[0] << std::endl;
  //  std::cerr << " HERE " << brv2.size() << std::endl;
  //}
  
  if (!brv2.size()) // no pieces had a good match
    return;	

  if ((brv2.size() > 1) || brv2[0].MaxDeletionBases() > 0  || brv2[0].MaxInsertionBases() > 0) {
    //if (max_mapq >= -1) 
    queue.add(new LongReaderWorkItem(ContigElement(brv2, regions)));
    for (const auto& b : brv2)
      contig_writer.WriteRecord(b);
  }
}
