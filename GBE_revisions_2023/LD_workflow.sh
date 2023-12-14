#For looking at decay of LD
#made folder /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD from other machine

mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023

mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/00_workflow_R

cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/00_workflow_R




#going to try with bash

#test it


cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/28_VCFTools_LD/01_VCFTools_LD/

head -1000 Sp34_Chr1_out.geno.ld > Sp34_Chr1_out.geno_top_1000.ld

awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $0, $3-$2}' Sp34_Chr1_out.geno_top_1000.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/Sp34_Chr1_out.geno_top_1000_diff.ld

cd /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/


awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($6) < 50001) print $0}' Sp34_Chr1_out.geno_top_1000_diff.ld > Sp34_Chr1_out.geno_top_1000_diff_less_50kb.ld


awk -F '\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($6) < 100001) print $0}' Sp34_Chr1_out.geno_top_1000_diff.ld > Sp34_Chr1_out.geno_top_1000_diff_less_100kb.ld

#ok, I'll think it will work

mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/01_awk_diff

mkdir /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/02_awk_abs_50k


#okay, trying again.....

awk 'BEGIN {FS="\t"}  {OFS="\t"} {print $0, $3-$2}' Sp34_Chr1_out.geno.ld > /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/01_awk_diff/Sp34_Chr1_out.geno_diff_copy.ld


#number of lines of chr 1 ld file, 820307261 Sp34_Chr1_out.geno.ld

# 820,307,261 Sp34_Chr1_out.geno.ld


#scp gcwoodruff@dtn2.oscer.ou.edu:/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/02_awk_abs_50k/Sp34_Chr1_out.geno_diff_less_50kb.ld /Users/gavin/genome/pop_gen_revisions_9-2023/LD/Sp34_Chr1_out.geno_diff_less_50kb.ld

#scp -r gcwoodruff@dtn2.oscer.ou.edu:/ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/2023/02_awk_abs_50k/ /Users/gavin/genome/pop_gen_revisions_9-2023/LD/
