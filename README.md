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
Convert a BND formatted SV VCF into a BEDPE
```bash
## convert a snowman vcf and sort with bedtools
svlib vcftobedpe a.snowman.sv.vcf | sortBed -i stdin
```

### sim
Simulate rearrangements and indels
```bash
svlib sim 
```

### realigntest
Test the ability of an aligner to realign an SV contig
```bash
svlib realigntest
```

### seqtovcf
Call SVs and indels from an qname sorted BAM of long sequence alignments
```bash
svlib seqtovcf
```

Attributions
============

svlib is developed and maintained by Jeremiah Wala (jwala@broadinstitute.org) --  Rameen Berkoukhim's lab -- Dana Farber Cancer Institute, Boston, MA

This project was developed in collaboration with the Cancer Genome Analysis team at the Broad Institute. Particular thanks to:
* Cheng-Zhong Zhang (Matthew Meyerson Lab)
* Marcin Imielinski

[license]: https://github.com/walaj/svlib/blob/master/LICENSE
