#ifndef SVLIB_SEQFRAG_H__
#define SVLIB_SEQFRAG_H__

#include <string>

#include "SeqLib/BWAWrapper.h"

using SeqLib::GenomicRegion;

struct Indel {
  
Indel() : len(0), type('N') {}

Indel(size_t l, char t, const std::string& rseq, const std::string& aseq, const std::string& lseq) : len(l), type(t) {
  ref_seq = rseq;
  alt_seq = aseq;
  lead_base = lseq;
  assert(lseq.length() == 1);
  assert(t == 'D' || alt_seq.length() == l); 
  assert(t == 'I' || ref_seq.length() == l); 
}

  size_t len;
  char type;
  std::string ref_seq, alt_seq, lead_base;
  SeqLib::GenomicRegion gr;
  int frag_id;

  friend std::ostream& operator<<(std::ostream& out, const Indel& i);


};

class SeqFrag {

 public:

 SeqFrag(const GenomicRegion& gr, faidx_t * findex) : m_gr(gr), m_index(findex) {}

  void getSeqFromRef(faidx_t * findex);

  int getLeftSide() const;

  int getRightSide() const;

  friend std::ostream& operator<<(std::ostream& out, const SeqFrag& s);

  std::string m_seq;

  void addScrambledEnds(size_t left_len, size_t right_len);

  char getStrand() const { return m_gr.strand; }

  void addIndels(size_t n);

  void addIns();

  void spikeMicrobe();

  std::vector<Indel> m_indels;  

  int frag_id;

  std::string left_scramble;
  std::string right_scramble;
 
  GenomicRegion m_gr;   
  
  int32_t phage_site = -1;
  std::string phage_string;

 private:

  faidx_t * m_index;

};


#endif
