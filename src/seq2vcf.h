#ifndef SVLIB_SEQ2VCF_H__
#define SVLIB_SEQ2VCF_H__

#include "gzstream.h"
#include <ctime>
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/RefGenome.h"
#include <map>
#include "SeqLib/GenomicRegionCollection.h"
#include "BreakPoint.h"
#include "STCoverage.h"
#include "pthread-lite.h"

void runSeqToVCF(int argc, char** argv);
void parseSeqToVCFOptions(int argc, char** argv);
void add_contig(const SeqLib::BamRecord& r, SeqLib::GRC& regions, SeqLib::BamRecordVector& brv,
		int& max_mapq, int& max_len);
void grab_reads(const std::string& prefix, SeqLib::BamReader& reader, const SeqLib::GRC& regions, SeqLib::BamRecordVector& brv,
		STCoverage * cov);

// hold a single contig and all of its alignments
struct ContigElement {

  SeqLib::BamRecordVector brv;
  SeqLib::GRC regions;

  // construct a new contig with all of its alignments
  ContigElement(const SeqLib::BamRecordVector& b, const SeqLib::GRC& r) : brv(b), regions(r) {}
};

// skeleton for any thread-specific data
struct LongReaderThreadItem {

LongReaderThreadItem(size_t i, const std::map<std::string, std::string>& map) : thread_id(i) {
  // open the short read readers
  for (const auto& b : map) {
    readers[b.first] = SeqLib::BamReader();
    if (!readers[b.first].Open(b.second)) {
      std::cerr << " failed to open BAM: " << b.second << std::endl;
      exit(EXIT_FAILURE);
    }
  }
    
}

  /** Load a reference genome for getting ref bases */
  bool LoadReference(const std::string& r) {
    return (ref.LoadIndex(r));
  }

  /** Write the files store in this thread and clear to save mem */
  void WriteFiles(pthread_mutex_t* lock) {
    
    // de duplicate the breakpoints
    std::sort(bp_glob.begin(), bp_glob.end());
    bp_glob.erase( std::unique( bp_glob.begin(), bp_glob.end() ), bp_glob.end() );
    
    pthread_mutex_lock(lock);
    for (auto& b : bp_glob)
      (*bps_file) << b.toFileString(true) << std::endl;
    
    //bps.str(std::string());
    //(*aln_file) << aln.str();
    //aln.str(std::string());
    pthread_mutex_unlock(lock);
    bp_glob.clear();
  }

  size_t thread_id = 0; // unique id for this thread
  
  std::stringstream bps; // stream to store BPS output
  std::stringstream aln; // stream to store ALN output
  
  std::vector<BreakPoint> bp_glob; // store all breakpoints

  size_t num_processed = 0; // count number processed 

  // each thread has unique reader
  std::map<std::string, SeqLib::BamReader> readers; 

  SeqLib::RefGenome ref; // pointer to index reference genome

  ogzstream * bps_file; // pointer to a file to dump the bps info
  ogzstream * aln_file; // pointer to a file to dump the aln info

};

// work item for a single long reader
class LongReaderWorkItem {
 
 public:
  
 LongReaderWorkItem(const ContigElement& c) : m_contig(c) {}
  
  // must take a thread item, even if not needed (for consistentcy with thread template)
  bool run(LongReaderThreadItem* thread_data) { 
    return runLR(thread_data); // first element should always be thread data
  } 
  
  // private:
  
  // data 
  ContigElement m_contig; // a set of aligned contigs with same qname
  
  bool runLR(LongReaderThreadItem* thread_item);
};

std::string fileDateString();

void assign_contig(const SeqLib::BamRecordVector& brv, WorkQueue<LongReaderWorkItem*>& queue,
		   const SeqLib::GRC& regions);

#endif
