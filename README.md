# _C. inopinata_ population genomics

Here are brief descriptions of the code and data associated with Woodruff et al. 2023 manuscript "Patterns of genomic diversity in a fig-associated close relative of _C. elegans_."

Downstream data used for analysis and figures can be found in the folder "data." Various additional files needed to run some commands can be found in the folder "additional_files."

FASTQ and BAM files have been deposited in the SRA (Bioproject PRJNA769443). VCF files have been deposited at Figshare (https://figshare.com/projects/C_inopinata_population_genomics_2021/123973).


Specifics about software versions can be found here.
```
software_and_versions.txt
```


This was used to process raw reads.
```
process_and_demultiplex_reads.sh 
```


This was used to align/filter reads, call genotypes, and estimate population genetic statistics (Pi, FST, FIS).
```
align_genotype_pop_gen.sh 
```


This was used to process previously-available _C. elegans_ alignments.
```
elegans.sh 
```


This was used to extract summary statistics/effect sizes and to perform PCA and DAPC.
```
statistics.R
```

This was used to make figures.
```
figures.R
```


These were used to extract site coverage information.
```
coverage.sh
```
```
coverage.R
```

This was used to make the last supplemental figure concerning nucleotide diversity estimates in RAD data.
```
 elegans_subsample_RAD.sh
```

Code for a previous version of this manuscript can be found here.
```
2021 code
```

Code for work performed towards revisions following peer review in 2023 can be here. This includes: comparisons of inbreeding coefficients among the X and autosomes; work on LD decay; and the examination of potentially hyperdivergent regions.
```
GBE_revisions_2023
```

If there are any questions about this please contact me at gcwoodruff@ou.edu.
