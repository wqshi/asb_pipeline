## The pipeline to extract allele specific binding events from CHiP-seq datasets.

Author: Wenqiang Shi.

Email:wqshi.nudt@gmail.com

Date: 2017 June 25th

## Introduction

This pipeline tries to extract allele specific binding events from ChIP_seq datasets. The pipeline takes as input both ChIP-Seq and called genotype data from the same cell. The called genotypes are used to create the associated personalized genome to map ChIP-Seq reads. The next step is to extract the number of mapped reads for each allele at known heterozygous sites within ChIP-seq peaks, and then an ASB event is called if the number of mapped ChIP-seq reads on one allele is significantly different from the other allele (Binomial test).

More details please see the citation.


## Example

Please see the test_script.sh in the ./python/



## Software used in the testing:

* Python 2.7.3

* pandas-0.12.0

* VCFtools (v0.1.11), making sure it's in the $PATH.

* Novoalign V3.01.00, including novoutil and novoindex. Should be in the $PATH. 

* R 3.1.3

* OS: CentOS 5

## To do:
* Add the scripts for read simulation.
* Adjust asb calling based on copy number variants

## Citation:
Wenqiang Shi, Oriol Fornes, Anthony Mathelier, Wyeth W. Wasserman; Evaluating the impact of single nucleotide variants on transcription factor binding. Nucleic Acids Res 2016; 44 (21): 10106-10116. doi: 10.1093/nar/gkw691






















