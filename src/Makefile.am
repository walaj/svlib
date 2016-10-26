bin_PROGRAMS = svlib

svlib_CPPFLAGS = \
     -I$(top_srcdir)/SeqLib \
     -I$(top_srcdir)/SeqLib/htslib -Wno-sign-compare

svlib_LDADD = \
	$(top_builddir)/SeqLib/src/libseqlib.a \
	$(top_builddir)/SeqLib/bwa/libbwa.a \
	$(top_builddir)/SeqLib/fermi-lite/libfml.a \
	$(top_builddir)/SeqLib/htslib/libhts.a 

svlib_SOURCES = svlib.cpp gzstream.C alignReadsToContigs.cpp vcf.cpp seq2vcf.cpp AlignmentFragment.cpp AlignedContig.cpp BreakPoint.cpp \
	DiscordantCluster.cpp STCoverage.cpp simulate.cpp PowerLawSim.cpp SeqFrag.cpp realigntest.cpp

