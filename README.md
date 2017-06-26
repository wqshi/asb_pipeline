##The pipeline to extract allele specific binding events from CHiP-seq datasets.
Author: Wenqiang Shi.
Email:wqshi.nudt@gmail.com
Date: 2017 June 25th

##Introduction
This pipeline tries to extract allele specific binding events from ChIP_seq datasets. The pipeline takes as input both ChIP-Seq and called genotype data from the same cell. The called genotypes are used to create the associated personalized genome to map ChIP-Seq reads. The next step is to extract the number of mapped reads for each allele at known heterozygous sites within ChIP-seq peaks, and then an ASB event is called if the number of mapped ChIP-seq reads on one allele is significantly different from the other allele (Binomial test).

More details please see the citation.


##Example
Please the test_script.sh in the ./python/



##Software needed:
Python 2.7.3
pandas-0.12.0
VCFtools (v0.1.11)
Novoalign V3.01.00
R 3.1.3

##To do:
Add the scripts for read simulation.
Consider the copy number variants

#Citation:
Shi, Wenqiang, et al. "Evaluating the impact of single nucleotide variants on transcription factor binding." Nucleic acids research 44.21 (2016): 10106-10116.
