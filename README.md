# _C. inopinata_ population genomics

Here are brief descriptions of the code and data associated with Woodruff et al. 2021 manuscript "Alignment of genetic differentiation across trophic levels in a fig community."

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
process_and_demultiplex_reads.sh 
```


This was used to process previously-available _C. elegans_ alignments.
```
elegans.sh 
```


This was used to extract summary statistics/effect sizes and to perform PCA and DAPC.
```
statistics.R
```

These were used to extract site coverage information.
```
coverage.sh
```
```
coverage.R
```

If there are any questions about this please contact me at gcwoodruff@ou.edu.
