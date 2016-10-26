#ifndef SVLIB_READ_TO_CONTIG_H
#define SVLIB_READ_TO_CONTIG_H

#include "SeqLib/BWAWrapper.h"
#include "SeqLib/RefGenome.h"
#include "AlignedContig.h"

void alignReadsToContigs(const SeqLib::BamHeader& hdr, const SeqLib::BWAWrapper& bw, 
			 const SeqLib::UnalignedSequenceVector& usv, 
			 SeqLib::BamRecordVector& bav_this, 
			 AlignedContig& this_alc, 
			 const SeqLib::RefGenome *  rg);

#endif
