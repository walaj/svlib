[![Build Status](https://travis-ci.org/walaj/svlib.svg?branch=master)](https://travis-ci.org/walaj/svlib)

## *svlib* - tools for benchmarking structural variations
========================================================================

**License:** [GPLv3][license]

Table of contents
=================

  * [Installation](#gh-md-toc)
  * [Description](#description)
    * [vcftobedpe](#vcftobedpe)
    * [sim](#sim)
    * [realigntest](#realigntest)
    * [seqtovcf](#seqtovcf)
  * [Attributions](#attributions)

Installation
============
I have succesfully built on Unix with GCC-4.8, 4.9 and 5.1

```
git clone --recursive https://github.com/walaj/svlib.git
cd svlib
./configure
make 
```

Description
===========

### vcftobedpe
Convert a BND formatted SV VCF into a BEDPE. Output format will be:
```bash
chr1 pos1 pos1 chr2 pos2 pos2 id 0 str1 str2 info genotypes
```
Examples:
```bash
## convert a snowman vcf and sort with bedtools
svlib vcftobedpe a.snowman.sv.vcf | sortBed -i stdin
```

### sim
Simulate rearrangements and indels
```bash
svlib sim -G $REFHG19 -s 42 -R 1000 -X 10000 -A mysim
```

### realigntest
Test the ability of an aligner to realign an SV contig
```bash
svlib realigntest -G $REFHG19 -b some.bam > sim.fasta
## align sim.fasta with some aligner (eg BWA) to something like sim.bam
svlib realigntest -G $REFHG19 -E sim.bam > results
```

### seqtovcf
Call SVs and indels from a qname sorted BAM of long sequence alignments
```bash
CORES=10
TUM=tumor.shortreads.bam ## eg standard Illumina reads for scoring/genotyping
NORM=normal.shortreads.bam ## should be coordinate sorted
svlib seqtovcf qsorted.contigs.bam -p $CORES -t $TUM -n $NORM -a output_id -G $REFHG19
```

##### Using seqtovcf with noisy BAMs
``seqtovcf`` can optionally use exclusion lists to avoid processing variants that 
fall in noise regions (e.g. centromeres), or to avoid re-running the same variants
supported by different long-sequences.

```bash
## first pass to find contigs which support a unique variant (without reads for speed)
## -W flag will write a *.extracted.contigs.bam, so you can see which ones were chosen
## in the case that multiple contigs support the same variant, svlib only uses the highest quality contig
svlib seqtovcf qsorted.contigs.bam -a first_pass -G $REFHG19 -W
gunzip -c first_pass.bps.txt.gz | cut -f20 | sort | uniq > contigs.to.process
## second pass (with reads for genotyping)
svlib seqtovcf qsorted.contigs.bam -L contigs.to.process -B excluded_regions.bed -a second_pass -t $TUM -n $NORM
```

Attributions
============

*svlib* is developed and maintained by Jeremiah Wala (jwala@broadinstitute.org) --  Rameen Berkoukhim's lab -- Dana Farber Cancer Institute, Boston, MA

This project was developed in collaboration with the Cancer Genome Analysis team at the Broad Institute. Particular thanks to:
* Cheng-Zhong Zhang (Matthew Meyerson Lab)
* Marcin Imielinski

[license]: https://github.com/walaj/svlib/blob/master/LICENSE
