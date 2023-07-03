library(ggplot2)
library(ggforce)
library(reshape2)
library(cowplot)
library(plyr)
library(rstatix)

#coverage before genotype filter

before_dat <- read.table("all_depth_before_genotype_filter.tsv",sep="\t",header=TRUE)

as.data.frame(tapply(before_dat$Cov,before_dat$sample,mean))

#A03.bam                                                18.944008
#A05.bam                                                11.483445
#B03.bam                                                 3.631566
#dec_2016_D01.bam                                       39.609752
#dec_2016_D02.bam                                       23.867116
#dec_2016_D03.bam                                       40.501326
#dec_2016_D04.bam                                       43.581105
#dec_2016_D05.bam                                       44.880955
#dec_2016_D10.bam                                       40.891739
#may_2017_A11.bam                                        8.523264
#may_2017_A12.bam                                        5.723795
#may_2017_B11.bam                                       42.292507
#may_2017_B12.bam                                        3.336316
#may_2017_C11.bam                                       42.495331
#may_2017_C12.bam                                       34.896262
#may_2017_D11.bam                                       14.077731
#may_2017_D12.bam                                       21.346691
#may_2017_E11.bam                                       41.879821
#may_2017_E12.bam                                       25.178702
#may_2017_F11.bam                                       49.765125
#may_2017_G12.bam                                       21.203651
#may_2017_H10.bam                                        5.714281
#may_2017_H11.bam                                       44.251051
#may_2017_H12.bam                                       18.649488

#after genotype filter
dat <- read.table("site_coverage.tsv",sep="\t",header=TRUE)

#rowmeans(dat[,-2])

dat_melt <- melt(dat, id=c("Chr","BP"))

names(dat_melt)[names(dat_melt) == 'variable'] <- 'sample'
names(dat_melt)[names(dat_melt) == 'value'] <- 'coverage'

as.data.frame(tapply(dat_melt$coverage,dat_melt$sample,mean))

#including 0's (that is, this includes sites where there were no reads for a given sample as 0x coverage.)

#A03                                                 86.072013
#A05                                                 42.039726
#B03                                                  4.558154
#dec_2016_D01                                       148.943142
#dec_2016_D02                                       100.611371
#dec_2016_D03                                       142.468046
#dec_2016_D04                                       148.582910
#dec_2016_D05                                       151.459927
#dec_2016_D10                                       142.462295
#may_2017_A11                                        35.830220
#may_2017_A12                                        11.741385
#may_2017_B11                                       148.069186
#may_2017_B12                                         5.708617
#may_2017_C11                                       147.480875
#may_2017_C12                                       130.539169
#may_2017_D11                                        65.963778
#may_2017_D12                                        90.936673
#may_2017_E11                                       146.058346
#may_2017_E12                                       106.973225
#may_2017_F11                                       152.052098
#may_2017_G12                                        93.502914
#may_2017_H10                                        19.535643
#may_2017_H11                                       150.428023
#may_2017_H12                                        71.630273

dat_melt_b <- dat_melt

dat_melt_b[, 4][dat_melt_b[, 4] == 0] <- NA

#checking stats of retained coverages-- this should have nothing less than 15x!!!
summary(dat_melt_b$coverage)

#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's
#      15       55       99      112      160      358 15099368

#yay



as.data.frame(tapply(dat_melt_b$coverage,dat_melt_b$sample,mean,na.rm=TRUE))

#coverage for each sample at the sites!!!!

#A03                                                                    86.58202
#A05                                                                    49.38279
#B03                                                                    21.55718
#dec_2016_D01                                                          148.95871
#dec_2016_D02                                                          101.00682
#dec_2016_D03                                                          142.79411
#dec_2016_D04                                                          148.74074
#dec_2016_D05                                                          151.70852
#dec_2016_D10                                                          143.25006
#may_2017_A11                                                           44.95167
#may_2017_A12                                                           31.34714
#may_2017_B11                                                          148.08684
#may_2017_B12                                                           28.70805
#may_2017_C11                                                          147.51059
#may_2017_C12                                                          131.37059
#may_2017_D11                                                           67.00712
#may_2017_D12                                                           91.51727
#may_2017_E11                                                          146.09961
#may_2017_E12                                                          107.46245
#may_2017_F11                                                          152.35357
#may_2017_G12                                                           93.96749
#may_2017_H10                                                           31.60188
#may_2017_H11                                                          150.44184
#may_2017_H12                                                           73.52293
#library(plyr)

dat_melt_c <- na.omit(dat_melt_b)


count(dat_melt_c$sample)

#number of sites included in each sample

#              x    freq
#1           A03 4945568
#2           A05 4235124
#3           B03 1051911
#4  dec_2016_D01 4974352
#5  dec_2016_D02 4955395
#6  dec_2016_D03 4963512
#7  dec_2016_D04 4969593
#8  dec_2016_D05 4966720
#9  dec_2016_D10 4947514
#10 may_2017_A11 3965387
#11 may_2017_A12 1863388
#12 may_2017_B11 4974279
#13 may_2017_B12  989257
#14 may_2017_C11 4973870
#15 may_2017_C12 4943387
#16 may_2017_D11 4897410
#17 may_2017_D12 4943311
#18 may_2017_E11 4973467
#19 may_2017_E12 4952224
#20 may_2017_F11 4965028
#21 may_2017_G12 4950276
#22 may_2017_H10 3075365
#23 may_2017_H11 4974415
#24 may_2017_H12 4846807

#there are 4974872 total sites in the set.

coun_dat <- count(dat_melt_c$sample)




coun_dat$perc_sites_retained <-  (coun_dat$freq/4974872)*100

cov_dat <- as.data.frame(tapply(dat_melt_b$coverage,dat_melt_b$sample,mean,na.rm=TRUE))

coun_dat$mean_coverage_sites_retained <- tapply(dat_melt_b$coverage,dat_melt_b$sample,mean,na.rm=TRUE)

coun_dat

#              x    freq perc_sites_retained mean_coverage_sites_retained
#1           A03 4945568            99.41096                     86.58202
#2           A05 4235124            85.13031                     49.38279
#3           B03 1051911            21.14448                     21.55718
#4  dec_2016_D01 4974352            99.98955                    148.95871
#5  dec_2016_D02 4955395            99.60849                    101.00682
#6  dec_2016_D03 4963512            99.77165                    142.79411
#7  dec_2016_D04 4969593            99.89389                    148.74074
#8  dec_2016_D05 4966720            99.83614                    151.70852
#9  dec_2016_D10 4947514            99.45008                    143.25006
#10 may_2017_A11 3965387            79.70832                     44.95167
#11 may_2017_A12 1863388            37.45600                     31.34714
#12 may_2017_B11 4974279            99.98808                    148.08684
#13 may_2017_B12  989257            19.88507                     28.70805
#14 may_2017_C11 4973870            99.97986                    147.51059
#15 may_2017_C12 4943387            99.36712                    131.37059
#16 may_2017_D11 4897410            98.44293                     67.00712
#17 may_2017_D12 4943311            99.36559                     91.51727
#18 may_2017_E11 4973467            99.97176                    146.09961
#19 may_2017_E12 4952224            99.54475                    107.46245
#20 may_2017_F11 4965028            99.80213                    152.35357
#21 may_2017_G12 4950276            99.50560                     93.96749
#22 may_2017_H10 3075365            61.81797                     31.60188
#23 may_2017_H11 4974415            99.99081                    150.44184
#24 may_2017_H12 4846807            97.42576                     73.52293


library(cowplot)

library(dplyr)

#samp_after_dat <- sample_n(dat_melt_b, 300000)

EcoRI_sites <- read.table("inopinata_EcoRI_sites.tsv",sep="\t",header=TRUE)


EcoRI_sites$BP_start.plusone <- EcoRI_sites$BP_start+1

EcoRI_sites$chr.bp_start <- as.factor(paste(EcoRI_sites$Chr,EcoRI_sites$BP_start.plusone))

EcoRI_sites_ii <- EcoRI_sites[EcoRI_sites$Chr != "CSP34.Sp34_mitochondrion",]

EcoRI_sites_ii$EcoRI_cut_sites <- as.numeric(gsub('\\.', 0, EcoRI_sites_ii$EcoRI_cut_sites))

summary(EcoRI_sites_ii$EcoRI_cut_sites)

#figure

before_dat$BP_start.plusone <- before_dat$BP_start+1

before_dat$chr.bp_start <- as.factor(paste(before_dat$Chr,before_dat$BP_start.plusone))

summary_before_dat <- before_dat %>% group_by(chr.bp_start) %>% get_summary_stats(Cov, type="common")

mean_before_dat <- data.frame(chr.bp_start = summary_before_dat$chr.bp_start, cov= summary_before_dat$mean)

out <- do.call(rbind,strsplit(as.character(mean_before_dat$chr.bp_start),' '))

mean_before_dat$Chr <- as.factor(out[,1])
mean_before_dat$BP <- as.numeric(out[,2])


mean_before_dat$MB <- mean_before_dat$BP/1000000

mean_before_dat_ii <- mean_before_dat[mean_before_dat$Chr != "CSP34.Sp34_mitochondrion",]


fig_df <- merge(x = mean_before_dat_ii, y = EcoRI_sites_ii, by = "chr.bp_start", all = TRUE)

fig_df_ii <- data.frame(Chr = fig_df$Chr.x, MB = fig_df$MB, Cov= fig_df$cov, EcoRI_cut_sites= fig_df$EcoRI_cut_sites)

fig_df_ii$coverage_per_cut_site = fig_df$cov/fig_df$EcoRI_cut_sites


fig_df_melt <- melt(fig_df_ii, id=c("Chr","MB"))



levels(fig_df_melt$variable)[levels(fig_df_melt$variable)=="Cov"] <- "Coverage"
levels(fig_df_melt$variable)[levels(fig_df_melt$variable)=="EcoRI_cut_sites"] <- "EcoRI cut sites"
levels(fig_df_melt$variable)[levels(fig_df_melt$variable)=="coverage_per_cut_site"] <- "Coverage per cut site"


fig_df_melt$variable <-factor(fig_df_melt$variable, levels = c("Coverage", "EcoRI cut sites", "Coverage per cut site"))


library(cowplot)
library(lemon)


fig_df_melt[is.na(fig_df_melt$value),]

fig_df_ii[is.infinite(fig_df_ii$coverage_per_cut_site),]




ggplot(fig_df_melt, aes(x = MB, y = value)) + geom_point(alpha=0.12, size=0.25,colour="black")  + stat_smooth(size=0.5, se=FALSE) + facet_rep_grid(variable ~ Chr,scales="free_y") + theme_cowplot() + theme(strip.background = element_blank()) + scale_x_continuous(breaks=c(0,10,20),minor_breaks=c(5,15))


#ggsave("coverage.png", width=7, height=6,units = "in")  








