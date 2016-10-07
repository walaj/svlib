#ifndef SVLIB_BREAKS_H__
#define SVLIB_BREAKS_H__

#include "SeqLib/GenomicRegion.h"

struct ReducedBreakEnd {
  
  ReducedBreakEnd() {}
  
  ReducedBreakEnd(const SeqLib::GenomicRegion& g, int mq, const std::string & chr_n);
  
  std::string chr_name;
  SeqLib::GenomicRegion gr;
  int32_t mapq:8, dummy:8, nm:16;
  
};

struct ReducedBreakPoint {

   char * ref;
   char * alt;
   char * cname;

};

#endif
