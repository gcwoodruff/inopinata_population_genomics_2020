#these were the slurm scripts used to:
	#Get differences between the two sites
	#Extract only those sites were the difference is less than or equal to 50kb.
	#Did it for each chromosome


#here's a slurm header, originally split in parallel to crunch chromosomes independently

#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --output=LD_I_%J_stdout.txt
#SBATCH --error=LD_I_%J_stderr.txt
#SBATCH --time=12:00:00
#SBATCH --job-name=LD_I
#
#################################################

#make some folders
mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD 

mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023

mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/01_awk_diff/

mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/02_awk_abs_50k/



#these were previously done analyses-- starting on line 1788 in align_genotype_pop_gen.sh 
cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/28_VCFTools_LD/01_VCFTools_LD/

#get difference between the two sites in bp
awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $0, $3-$2}' Sp34_Chr1_out.geno.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/01_awk_diff/Sp34_Chr1_out.geno_diff.ld

awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $0, $3-$2}' Sp34_Chr2_out.geno.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/01_awk_diff/Sp34_Chr2_out.geno_diff.ld

awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $0, $3-$2}' Sp34_Chr3_out.geno.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/01_awk_diff/Sp34_Chr3_out.geno_diff.ld

awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $0, $3-$2}' Sp34_Chr4_out.geno.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/01_awk_diff/Sp34_Chr4_out.geno_diff.ld

awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $0, $3-$2}' Sp34_Chr5_out.geno.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/01_awk_diff/Sp34_Chr5_out.geno_diff.ld

awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $0, $3-$2}' inopinata_females_X_out.geno.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/01_awk_diff/Sp34_ChrX_out.geno_diff.ld


#get the absolute value of that difference and remove differences >50kb
cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/01_awk_diff/

awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($6) < 50001) print $0}' Sp34_Chr1_out.geno_diff.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/02_awk_abs_50k/Sp34_Chr1_out.geno_diff_less_50kb.ld

awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($6) < 50001) print $0}' Sp34_Chr2_out.geno_diff.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/02_awk_abs_50k/Sp34_Chr2_out.geno_diff_less_50kb.ld

awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($6) < 50001) print $0}' Sp34_Chr3_out.geno_diff.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/02_awk_abs_50k/Sp34_Chr3_out.geno_diff_less_50kb.ld

awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($6) < 50001) print $0}' Sp34_Chr4_out.geno_diff.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/02_awk_abs_50k/Sp34_Chr4_out.geno_diff_less_50kb.ld

awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($6) < 50001) print $0}' Sp34_Chr5_out.geno_diff.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/02_awk_abs_50k/Sp34_Chr5_out.geno_diff_less_50kb.ld

awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($6) < 50001) print $0}' Sp34_ChrX_out.geno_diff.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/02_awk_abs_50k/Sp34_ChrX_out.geno_diff_less_50kb.ld

#next, workflow.R