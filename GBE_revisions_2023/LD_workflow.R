# /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/28_VCFTools_LD/01_VCFTools_LD

# /ourdisk/hpc/figwormlab/gcwoodruff/dont_archive/pop_gen_revisions_2023/LD/28_VCFTools_LD/01_VCFTools_LD/Sp34_Chr1_out.geno.ld
#load some libraries
library(ggplot2)
library(cowplot)
library(patchwork)
library(lemon)

#set working directory
setwd("/Users/gavin/genome/pop_gen_revisions_9-2023/LD/02_awk_abs_50k")
#get data in there
Idat <- read.table("Sp34_Chr1_out.geno_diff_less_50kb.ld",sep="\t",header=TRUE)
IIdat <- read.table("Sp34_Chr2_out.geno_diff_less_50kb.ld",sep="\t",header=TRUE)
IIIdat <- read.table("Sp34_Chr3_out.geno_diff_less_50kb.ld",sep="\t",header=TRUE)
IVdat <- read.table("Sp34_Chr4_out.geno_diff_less_50kb.ld",sep="\t",header=TRUE)
Vdat <- read.table("Sp34_Chr5_out.geno_diff_less_50kb.ld",sep="\t",header=TRUE)
Xdat <- read.table("Sp34_ChrX_out.geno_diff_less_50kb.ld",sep="\t",header=TRUE)

#okay, what are the sample sizes, really
mean(Idat$N_INDV)
#[1] 19.85234
mean(IIdat$N_INDV)
#[1] 19.74694
mean(IIIdat$N_INDV)
#[1] 19.80741
mean(IVdat$N_INDV)
#[1] 19.73034
mean(Vdat$N_INDV)
#[1] 19.71218
mean(Xdat$N_INDV)
#[1] 16.62216

#okay, about 20 for autosomes and 17 for X
one_over_N_aut <- 1/20

one_over_N_X <- 1/17

#ggplot(dat, aes(x=BP_diff,y=R.2)) + geom_hline(yintercept=one_over_N,linetype="dotted") + geom_point() + geom_smooth() + theme_cowplot()
	#okay, way way way too many points
	#let's try aggregate

#aggregating mean R2 by the distance between sites
#X0 is the absolute value of the difference between sites

Iaggdatmean <- aggregate(R.2 ~ X0, FUN=mean, data=Idat)
IIaggdatmean <- aggregate(R.2 ~ X0, FUN=mean, data=IIdat)
IIIaggdatmean <- aggregate(R.2 ~ X0, FUN=mean, data=IIIdat)
IVaggdatmean <- aggregate(R.2 ~ X0, FUN=mean, data=IVdat)
Vaggdatmean <- aggregate(R.2 ~ X0, FUN=mean, data=Vdat)
Xaggdatmean <- aggregate(R.2 ~ X0, FUN=mean, data=Xdat)

#add chr
Iaggdatmean$chr <- "I"
IIaggdatmean$chr <- "II"
IIIaggdatmean$chr <- "III"
IVaggdatmean$chr <- "IV"
Vaggdatmean$chr <- "V"
Xaggdatmean$chr <- "X"

#combine all
all_agg_dat <- rbind(Iaggdatmean,IIaggdatmean,IIIaggdatmean,IVaggdatmean,Vaggdatmean,Xaggdatmean)

#get bp in kb
all_agg_dat$kb_diff <- all_agg_dat$X0/1000
#just get sites 20kb or less
all_agg_dat_20kb <- all_agg_dat[all_agg_dat$kb_diff<20.001,]

#this is the supplementlal figure
ggplot(all_agg_dat_20kb,aes(x=kb_diff,y=R.2)) + facet_rep_wrap(~chr) + geom_point(alpha=0.05,size=0.25) + geom_smooth(size=0.5) + geom_hline(yintercept=one_over_N_aut,linetype="dotted",colour="#1b9e77") + geom_hline(yintercept=one_over_N_X,linetype="dashed",colour="#d95f02")  + theme_cowplot() + theme(strip.background = element_blank()) + scale_y_continuous(limits=c(0,0.5),breaks=seq(0,0.5,0.05)) + scale_x_continuous(limits=c(0,20),breaks=seq(0,20,2)) + labs(y = expression(italic("r"^2)), x="Distance between sites (kb)")
#save it
ggsave("LD_decay.png",bg="white",units="in",height=7,width=10)


#where is the asymptote of these fits??

p <- ggplot(all_agg_dat_20kb,aes(x=kb_diff,y=R.2)) + facet_rep_wrap(~chr) + geom_point(alpha=0.05,size=0.25) + geom_smooth(size=0.5) + geom_hline(yintercept=one_over_N_aut,linetype="dotted",colour="#1b9e77") + geom_hline(yintercept=one_over_N_X,linetype="dashed",colour="#d95f02")  + theme_cowplot() + theme(strip.background = element_blank()) + scale_y_continuous(limits=c(0,0.5),breaks=seq(0,0.5,0.05)) + scale_x_continuous(limits=c(0,20),breaks=seq(0,20,2)) + labs(y = expression(italic("r"^2)), x="Distance between sites (kb)")

fitdat <- ggplot_build(p)$data[[2]]

aggregate(y ~ PANEL, FUN=min, data=fitdat)

#  PANEL          y
#1     1 0.07110303
#2     2 0.07515192
#3     3 0.08220562
#4     4 0.07516240
#5     5 0.07027378
#6     6 0.07863917











 