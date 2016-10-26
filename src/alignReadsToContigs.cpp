#include "alignReadsToContigs.h"

void alignReadsToContigs(const SeqLib::BamHeader& hdr, const SeqLib::BWAWrapper& bw, 
			 const SeqLib::UnalignedSequenceVector& usv, 
			 SeqLib::BamRecordVector& bav_this, 
			 AlignedContig& this_alc, 
			 const SeqLib::RefGenome * rg) {
   
  if (!usv.size())
    return;
 
  // get the reference info
  SeqLib::GRC g;
  for (auto& i : this_alc.getAsGenomicRegionVector()) {
    i.Pad(100);
    g.add(i);
  }
  g.MergeOverlappingIntervals();
 
  // get the reference sequence
  std::vector<std::string> ref_alleles;
  for (auto& i : g)
    if (i.chr < 24) //1-Y
      try {
	std::string tmpref = rg->QueryRegion(i.ChrName(hdr), i.pos1, i.pos2);
	ref_alleles.push_back(tmpref); 
      } catch (...) {
	std::cerr << "Caught exception for ref_allele on " << i << std::endl;
      }
  // make the reference allele BWAWrapper
  SeqLib::BWAWrapper bw_ref;
  SeqLib::UnalignedSequenceVector usv_ref;
  int aa = 0;
  for (auto& i : ref_alleles) {
    if (!i.empty())
      usv_ref.push_back({std::to_string(aa++), i, std::string()}); // name, seq, qual
  }
  if (!usv_ref.size())
    bw_ref.ConstructIndex(usv_ref);
   
  // set up custom alignment parameters, mean
  bw_ref.SetGapOpen(16); // default 6
  bw_ref.SetMismatchPenalty(9); // default 2
 
  for (auto i : bav_this) {
     
    SeqLib::BamRecordVector brv, brv_ref;
 
    // try the corrected seq first
    std::string seqr = i.GetZTag("KC");
    if (seqr.empty())
      seqr = i.QualitySequence();
     
    bool hardclip = false;
    assert(seqr.length());
    bw.AlignSequence(seqr, i.Qname(), brv, hardclip, 0.60, 10000);
 
    if (brv.size() == 0) 
      continue;
 
    // get the maximum non-reference alignment score
    int max_as = 0;
    for (auto& r : brv)
      max_as = std::max(max_as, r.GetIntTag("AS"));
 
    // align to the reference alleles
    if (!bw_ref.IsEmpty())
      bw_ref.AlignSequence(seqr, i.Qname(), brv_ref, hardclip, 0.60, 10);
 
    // get the maximum reference alignment score
    int max_as_r = 0;
    for (auto& r : brv_ref)
      max_as_r = std::max(max_as_r, r.GetIntTag("AS"));
     
    // reject if better alignment to reference
    if (max_as_r > max_as) {
      //std::cerr << " Alignment Rejected for " << max_as_r << ">" << max_as << "  " << i << std::endl;
      //std::cerr << "                        " << max_as_r << ">" << max_as << "  " << brv_ref[0] << std::endl;
      continue;
    }
     
    // make sure we have only one alignment per contig
    std::set<std::string> cc;
 
    // check which ones pass
    SeqLib::BamRecordVector bpass;
    for (auto& r : brv) {
 
      // make sure alignment score is OK
      if ((double)r.NumMatchBases() * 0.5 > r.GetIntTag("AS")/* && i.GetZTag("SR").at(0) == 't'*/)
        continue;
       
      bool length_pass = (r.PositionEnd() - r.Position()) >= ((double)seqr.length() * 0.75);
 
      if (length_pass && !cc.count(usv[r.ChrID()].Name)) {
	bpass.push_back(r);
	cc.insert(usv[r.ChrID()].Name);
      }
    }
 
    // annotate the original read
    for (auto& r : bpass) {
      if (r.ReverseFlag())
	i.SmartAddTag("RC","1");
      else
	i.SmartAddTag("RC","0");
 
      i.SmartAddTag("SL", std::to_string(r.Position()));
      i.SmartAddTag("SE", std::to_string(r.PositionEnd()));
      i.SmartAddTag("TS", std::to_string(r.AlignmentPosition()));
      i.SmartAddTag("TE", std::to_string(r.AlignmentEndPosition()));
      i.SmartAddTag("SC", r.CigarString());
      i.SmartAddTag("CN", usv[r.ChrID()].Name);
 
      if (this_alc.getContigName() != usv[r.ChrID()].Name)
	continue;
      this_alc.AddAlignedRead(i);
      
    } // end passing bwa-aligned read loop 
  } // end main read loop
}
 
