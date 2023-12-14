#load libraries 
library(ggplot2)
library(cowplot)
library(lemon)
library(ggforce)
library(patchwork)
library(car)

#set working directory
setwd("/Users/gavin/genome/pop_gen_revisions_9-2023/X_FIS")
#get fis data in there
FISdat <- read.table("fis_elegans_and_inopinata.tsv", header=T,sep="\t")
#species, chr as factor
FISdat$species <- as.factor(FISdat$species)
FISdat$Chr <- as.factor(FISdat$Chr)

#chromosome names

Chr <-c("I","II","III","IV","V","X")
###chromosome sizes
endsCI<-c(20594552,20116996,19436473,21017644,23638155,18190508)
endsCE<-c(15072434,15279421,13783801,17493829,20924180,17718942)

#normalized chr position
FISdat$NORMPOS<-888

for (i in 1:6){

  FISdat[FISdat$species=="elegans" & FISdat$Chr==Chr[i],]$NORMPOS <- FISdat[FISdat$species=="elegans" & FISdat$Chr==Chr[i],]$BP/endsCE[i]
  FISdat[FISdat$species=="inopinata" & FISdat$Chr==Chr[i],]$NORMPOS <- FISdat[FISdat$species=="inopinata" & FISdat$Chr==Chr[i],]$BP/endsCI[i]

}

#get pi data in there


PIdat <- read.csv("ino_elg_pi.csv", header=T)
#species, chr as factor
PIdat$species <- as.factor(PIdat$species)
PIdat$scaffold <- as.factor(PIdat$scaffold)

levels(PIdat$species)[match("C. inopinata",levels(PIdat$species))] <- "inopinata"
levels(PIdat$species)[match("C. elegans",levels(PIdat$species))] <- "elegans"


#normalized chr position
PIdat$NORMPOS<-888

for (i in 1:6){

  PIdat[PIdat$species=="elegans" & PIdat$scaffold==Chr[i],]$NORMPOS <- PIdat[PIdat$species=="elegans" & PIdat$scaffold==Chr[i],]$start/endsCE[i]
  PIdat[PIdat$species=="inopinata" & PIdat$scaffold==Chr[i],]$NORMPOS <- PIdat[PIdat$species=="inopinata" & PIdat$scaffold==Chr[i],]$start/endsCI[i]

}
#get a column of chromosome and base pair position
FISdat$chr_pos <- paste(FISdat$Chr,FISdat$BP)
#subtract one to align the two data types
PIdat$start_m_one <- PIdat$start-1
PIdat$chr_pos <- paste(PIdat$scaffold, PIdat$start_m_one)

#subset by species
inopPI <- PIdat[PIdat$species == "inopinata",]
inopFIS <- FISdat[FISdat$species == "inopinata",]
elegPI  <- PIdat[PIdat$species == "elegans",]
elegFIS <- FISdat[FISdat$species == "elegans",]


#join FIS and PI data to same sites
inop_merge <- merge(inopFIS,inopPI, by="chr_pos")

eleg_merge <- merge(elegFIS,elegPI, by="chr_pos")


#column for "is this site X or autosome" ?

inop_merge$Xaut <- ifelse(inop_merge$Chr=="X", "X", "Autosome")
eleg_merge$Xaut <- ifelse(eleg_merge$Chr=="X", "X", "Autosome")


#linear model, is there a relationship between Fis and Pi in inopinata
summary(lm(Fis ~ pi_all, data=inop_merge))

#Call:
#lm(formula = Fis ~ pi_all, data = inop_merge)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-0.81176 -0.11146 -0.01980  0.08295  0.83137
#
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.184517   0.003653   50.51   <2e-16 ***
#pi_all      -3.695177   0.266000  -13.89   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.1728 on 7040 degrees of freedom
#Multiple R-squared:  0.02668, Adjusted R-squared:  0.02654
#F-statistic:   193 on 1 and 7040 DF,  p-value: < 2.2e-16

#remove X chromosome
inop_merge_no_X <- inop_merge[inop_merge$Chr != "X",]

#linear model, is there a relationship between Fis and Pi in inopinata when just considering autosomes
summary(lm(Fis ~ pi_all, data=inop_merge_no_X))

#Call:
#lm(formula = Fis ~ pi_all, data = inop_merge_no_X)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-0.77017 -0.08713 -0.00188  0.07917  0.87347
#
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.127475   0.003121  40.851  < 2e-16 ***
#pi_all      -1.578513   0.220039  -7.174 8.13e-13 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.1378 on 6304 degrees of freedom
#Multiple R-squared:  0.008097,  Adjusted R-squared:  0.00794
#F-statistic: 51.46 on 1 and 6304 DF,  p-value: 8.131e-13

#look only at the X chromosome
inop_merge_no_X <- inop_merge[inop_merge$Chr != "X",]

#linear model, is there a relationship between Fis and Pi in inopinata when just considering X
inop_merge_X <- inop_merge[inop_merge$Chr == "X",]

summary(lm(Fis ~ pi_all, data=inop_merge_X))

#Call:
#lm(formula = Fis ~ pi_all, data = inop_merge_X)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-0.58685 -0.11480  0.00236  0.11853  0.56150
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  0.44762    0.01279  35.004   <2e-16 ***
#pi_all      -2.12189    1.41693  -1.498    0.135
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.1903 on 734 degrees of freedom
#Multiple R-squared:  0.003046,  Adjusted R-squared:  0.001688
#F-statistic: 2.243 on 1 and 734 DF,  p-value: 0.1347


#put both species data frames together
all_merge <- rbind(inop_merge,eleg_merge)


#Ne estimates
mean(sapply(inopPI$pi_all, function(th) { return(th/(2*4*2.3e-09))}))
#[1] 613505.9

mean(sapply(elegPI$pi_all, function(th) { return(th/(2*4*2.3e-09))}))
#[1] 117894.3

(mean(inopPI$pi_all))/(2*4*2.3e-09)
#613505.9

(mean(elegPI$pi_all))/(2*4*2.3e-09)
#117894.3

#okay, pi on the X for inopinata

inop_merge_X <- inop_merge[inop_merge$Xaut == "X",]
inop_merge_aut <- inop_merge[inop_merge$Xaut != "X",]


summary(inop_merge_X$pi_all)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.000200 0.004100 0.006750 0.007545 0.009800 0.041100

summary(inop_merge_aut$pi_all)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.00020 0.00630 0.01010 0.01179 0.01530 0.08500



summary(inop_merge_X$Fis)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-0.1753  0.3175  0.4314  0.4316  0.5456  1.0000

summary(inop_merge_aut$Fis)

#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#-0.65422  0.02251  0.10646  0.10887  0.18765  1.00000

#observed aut mean fis = 0.10887
#expected mean pi of x,
mean(inop_merge_aut$Fis)*(4/3)
#0.1451582
  #way higher

#what about elegans?

#X and aut for elegans
eleg_merge_X <- eleg_merge[eleg_merge$Xaut == "X",]
eleg_merge_aut <- eleg_merge[eleg_merge$Xaut != "X",]


summary(lm(Fis ~ pi_all, data=eleg_merge))
#Call:
#lm(formula = Fis ~ pi_all, data = eleg_merge)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-1.22760  0.02515  0.04267  0.04916  0.56469
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.961659   0.001711  562.17   <2e-16 ***
#pi_all      -10.820870   0.343224  -31.53   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.1472 on 9235 degrees of freedom
#Multiple R-squared:  0.09717, Adjusted R-squared:  0.09707
#F-statistic:   994 on 1 and 9235 DF,  p-value: < 2.2e-16
 #y=-10.82x+0.96, p<0.0001, r2=0.097
 
summary(lm(Fis ~ pi_all, data=eleg_merge_aut))
#
#Call:
#lm(formula = Fis ~ pi_all, data = eleg_merge_aut)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-1.22646  0.03015  0.04381  0.04990  0.53238
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.960243   0.001939   495.3   <2e-16 ***
#pi_all      -10.138130   0.356954   -28.4   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.1508 on 7567 degrees of freedom
#Multiple R-squared:  0.09633, Adjusted R-squared:  0.09621
#F-statistic: 806.7 on 1 and 7567 DF,  p-value: < 2.2e-16
 #y=-10.14x+0.96, p<0.0001, r2=0.15

#okay, pi on the X for elegans


summary(eleg_merge_X$pi_all)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.000100 0.000600 0.001000 0.001251 0.001500 0.036700

summary(eleg_merge_aut$pi_all)

#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.000000 0.000500 0.001000 0.002432 0.002200 0.059600
#observed mean pi = 0.001251
#expected mean pi of x,
mean(eleg_merge_aut$pi_all)*0.75
#0.001824178
  #lower


summary(eleg_merge_aut$Fis)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-0.3692  0.9637  1.0000  0.9356  1.0000  1.0000

summary(eleg_merge_X$Fis)

#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#-0.04679  0.96040  1.00000  0.94702  1.00000  1.00000


#C. remanei data from paper Teterina at al. "Genomic diversity landscapes in outcrossing and selfing Caenorhabditis nematodes" https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010879 ; https://figshare.com/articles/dataset/Supplementary_files_for_Genomic_diversity_landscapes_in_outcrossing_and_selfing_Caenorhabditis_nematodes_/23826486
#FIS
tetFISdat <- read.table("diversity_stats_Fis_C.elegans_C.remanei.txt", header=T,sep="\t",stringsAsFactors=T)

#get just remanei
remFISdat <- tetFISdat[tetFISdat$Species == "C.remanei",]
#autosomes and X
remFISdat_aut <- remFISdat[remFISdat$chrom != "X",]

remFISdat_X <- remFISdat[remFISdat$chrom == "X",]

summary(remFISdat_aut$Fis)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.000008 0.280368 0.382170 0.386261 0.480395 1.000000

summary(remFISdat_X$Fis)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.03365 0.22594 0.32835 0.34214 0.43584 0.85457

wilcox.test(remFISdat_aut$Fis,remFISdat_X$Fis)
#  Wilcoxon rank sum test with continuity correction
#
#data:  remdat_aut$Fis and remdat_X$Fis
#W = 127494, p-value = 5.277e-05
#alternative hypothesis: true location shift is not equal to 0


#ok, pi

tetdivdat <- read.table("diversity_stats_diploSHIC_BETA_C.elegans_C.remanei.txt", header=T,sep="\t",stringsAsFactors=T)


#get just remanei
remdat <- tetdivdat[tetdivdat$Species == "C.remanei",]
#autosomes and X
remdat_aut <- remdat[remdat$chrom != "X",]

remdat_X <- remdat[remdat$chrom == "X",]

summary(remdat_aut$pi)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.004095 0.011284 0.015018 0.015409 0.018810 0.031992

summary(remdat_X$pi)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.002572 0.010123 0.012330 0.013329 0.015306 0.032918

#observed mean pi = 0.013329
#expected mean pi of x,
mean(remdat_aut$pi)*0.75
#0.01155641
  #higher

#make a chr & pos column
remdat$chr_pos <- paste(remdat$chrom,remdat$classifiedWinStart.x)

remFISdat$chr_pos <- paste(remFISdat$chrom,remFISdat$start)
#join remanei FIS and Pi data by site
rem_merge <- merge(remdat,remFISdat, by="chr_pos")

#compare X/autosome ratio of population genetic statistics with some null hypothetical value

#eleg pi X:Aut ratio
eleg_X_aut_ratio_pi <- mean(eleg_merge_X$pi_all)/mean(eleg_merge_aut$pi_all)
#0.5142014

#bootstrap to get 95% CI

# function to obtain the mean 
Bmean <- function(data, indices) {
  d <- data[indices] # allows boot to select sample 
    return(mean(d))
} 

elegpibootX <- boot(data=eleg_merge_X$pi_all, statistic=Bmean, R=10000)
# get 95% confidence interval 
elegpibootXci <- confint( elegpibootX, level=.95, type='perc')  
elegpibootaut <- boot(data=eleg_merge_aut$pi_all, statistic=Bmean, R=10000)
# get 95% confidence interval 
elegpibootautci <- confint(elegpibootaut, level=.95, type='perc')
#make a row of X:aut ratio and 95% CI  
eleg_pi_row <- c(eleg_X_aut_ratio_pi,elegpibootXci[1]/elegpibootautci[1],elegpibootXci[2]/elegpibootautci[2], "Pi", "C. elegans")

#inop pi X:Aut ratio
inop_X_aut_ratio_pi <- mean(inop_merge_X$pi_all)/mean(inop_merge_aut$pi_all)
#0.6401462
#bootstrap to get 95% CI
inoppibootX <- boot(data=inop_merge_X$pi_all, statistic=Bmean, R=10000)

# get 95% confidence interval 
inoppibootXci <- confint( inoppibootX, level=.95, type='perc')  
inoppibootaut <- boot(data=inop_merge_aut$pi_all, statistic=Bmean, R=10000)

# get 95% confidence interval 
inoppibootautci <- confint( inoppibootaut, level=.95, type='perc')  
#make a row of X:aut ratio and 95% CI  
inop_pi_row <- c(inop_X_aut_ratio_pi,inoppibootXci[1]/inoppibootautci[1],inoppibootXci[2]/inoppibootautci[2], "Pi", "C. inopinata")


#rem pi X:aut ratio
rem_X_aut_ratio_pi <- mean(remdat_X$pi)/mean(remdat_aut$pi)
#0.865063
#bootstrap to get 95% CI
rempibootX <- boot(data=remdat_X$pi, statistic=Bmean, R=10000)

# get 95% confidence interval 
rempibootXci <- confint( rempibootX, level=.95, type='perc')  
rempibootaut <- boot(data=remdat_aut$pi, statistic=Bmean, R=10000)

# get 95% confidence interval 
rempibootautci <- confint( rempibootaut, level=.95, type='perc')  

#make a row of X:aut ratio and 95% CI 
rem_pi_row <- c(rem_X_aut_ratio_pi,rempibootXci[1]/rempibootautci[1],rempibootXci[2]/rempibootautci[2], "Pi", "C. remanei")

#eleg fis X:aut ratio
eleg_X_aut_ratio_fis <- mean(eleg_merge_X$Fis)/mean(eleg_merge_aut$Fis)
#1.012219
#bootstrap to get 95% CI
elegfisbootX <- boot(data=eleg_merge_X$Fis, statistic=Bmean, R=10000)
# get 95% confidence interval 
elegfisbootXci <- confint( elegfisbootX, level=.95, type='perc')  
elegfisbootaut <- boot(data=eleg_merge_aut$Fis, statistic=Bmean, R=10000)
# get 95% confidence interval 
elegfisbootautci <- confint(elegfisbootaut, level=.95, type='perc')  


#make a row of X:aut ratio and 95% CI 
eleg_fis_row <- c(eleg_X_aut_ratio_fis,elegfisbootXci[1]/elegfisbootautci[1],elegfisbootXci[2]/elegfisbootautci[2], "Fis", "C. elegans")



#inop fis X:aut ratio
inop_X_aut_ratio_Fis <- mean(inop_merge_X$Fis)/mean(inop_merge_aut$Fis)
#bootstrap to get 95% CI
inopFisbootX <- boot(data=inop_merge_X$Fis, statistic=Bmean, R=10000)

# get 95% confidence interval 
inopFisbootXci <- confint( inopFisbootX, level=.95, type='perc')  
inopFisbootaut <- boot(data=inop_merge_aut$Fis, statistic=Bmean, R=10000)

# get 95% confidence interval 
inopFisbootautci <- confint( inopFisbootaut, level=.95, type='perc')  


#make a row of X:aut ratio and 95% CI 

inop_fis_row <- c(inop_X_aut_ratio_Fis,inopFisbootXci[1]/inopFisbootautci[1],inopFisbootXci[2]/inopFisbootautci[2], "Fis", "C. inopinata")


#rem fis X:aut ratio
rem_X_aut_ratio_Fis <- mean(remFISdat_X$Fis)/mean(remFISdat_aut$Fis)
#0.8857847
#bootstrap to get 95% CI
remFisbootX <- boot(data=remFISdat_X$Fis, statistic=Bmean, R=10000)

# get 95% confidence interval 
remFisbootXci <- confint( remFisbootX, level=.95, type='perc')  
remFisbootaut <- boot(data=remFISdat_aut$Fis, statistic=Bmean, R=10000)

# get 95% confidence interval 
remFisbootautci <- confint( remFisbootaut, level=.95, type='perc')  

#make a row of X:aut ratio and 95% CI 
rem_fis_row <- c(rem_X_aut_ratio_Fis,remFisbootXci[1]/remFisbootautci[1],remFisbootXci[2]/remFisbootautci[2], "Fis", "C. remanei")
#combine all the X:aut ratio rows to make a df for plotting
rat_dat <- as.data.frame(rbind(eleg_pi_row,inop_pi_row,rem_pi_row,eleg_fis_row,inop_fis_row,rem_fis_row))
#as numeric
rat_dat$V1 <- as.numeric(rat_dat$V1)
rat_dat$V2 <- as.numeric(rat_dat$V2)
rat_dat$V3 <- as.numeric(rat_dat$V3)
ggplot(rat_dat,aes(x=V5, y=V1)) + geom_bar(stat="identity",fill="lightblue") + geom_errorbar(aes(ymin=V2, ymax=V3), colour="black", width=.1, position=position_dodge(0.1)) + facet_rep_wrap(~ V4)
#set aside by statistic
rat_dat_pi <- rat_dat[rat_dat$V4 == "Pi",]

rat_dat_fis <- rat_dat[rat_dat$V4 == "Fis",]

#this is a supplemental figure for X:Autosome ratio
a <- ggplot(rat_dat_pi,aes(x=V5, y=V1)) + geom_bar(stat="identity",fill="lightblue") + geom_errorbar(aes(ymin=V2, ymax=V3), colour="black", width=.1, position=position_dodge(0.1)) + geom_hline(yintercept=0.75,linetype="dotted") + theme_cowplot() + theme(axis.text.x = element_text(angle = 45,hjust=1)) + scale_y_continuous(limits=c(0,1),breaks=c(0,0.25,0.50,0.75,1.0)) + ylab("PiX/Piautosome") + xlab("Species")

b <- ggplot(rat_dat_fis,aes(x=V5, y=V1)) + geom_bar(stat="identity",fill="lightblue") + geom_errorbar(aes(ymin=V2, ymax=V3), colour="black", width=.1, position=position_dodge(0.1)) + geom_hline(yintercept=1.333,linetype="dotted") + theme_cowplot() + theme(axis.text.x = element_text(angle = 45,hjust=1)) + scale_y_continuous(limits=c(0,4.0),breaks=seq(0,4,0.5)) + ylab("Fx/F") + xlab("Species")

a+b
#ggsave("x_aut_ratios.pdf", units="in",height=5,width=6)


#this is a supplemental figure ; scatterplot of FIS and Pi

#column for "is this site X or autosome" ?

rem_merge$Xaut <- ifelse(rem_merge$chrom.x=="X", "X", "Autosome")

rem_merge_aut <- rem_merge[rem_merge$Xaut == "Autosome",]


a <- ggplot(inop_merge, aes(x=pi_all,y=Fis)) + geom_point(alpha=0.25, aes(colour=Xaut)) + geom_smooth(method="lm", linetype="dotted", se=FALSE) + scale_colour_manual(values=c("black","red")) + theme_cowplot() + xlab("Pi") + theme(legend.position="none")



b <- ggplot(eleg_merge, aes(x=pi_all,y=Fis)) + geom_point(alpha=0.25, aes(colour=Xaut)) + geom_smooth(method="lm", linetype="dotted", se=FALSE)  + scale_colour_manual(values=c("black","red")) + theme_cowplot()  + xlab("Pi") + theme(legend.position="none")


c <- ggplot(rem_merge, aes(x=pi,y=Fis)) + geom_point(alpha=0.25, aes(colour=Xaut)) + scale_colour_manual(values=c("black","red")) + geom_smooth(method="lm", linetype="dotted", se=FALSE) + theme_cowplot()  + xlab("Pi")


d <- ggplot(inop_merge_aut, aes(x=pi_all,y=Fis)) + geom_point(alpha=0.25) + geom_smooth(method="lm", linetype="dotted", se=FALSE) + theme_cowplot()  + xlab("Pi")

e <- ggplot(eleg_merge_aut, aes(x=pi_all,y=Fis)) + geom_point(alpha=0.25) + geom_smooth(method="lm", linetype="dotted", se=FALSE) + theme_cowplot()  + xlab("Pi")

f <- ggplot(rem_merge_aut, aes(x=pi,y=Fis)) + geom_point(alpha=0.25) + geom_smooth(method="lm", linetype="dotted", se=FALSE)  + theme_cowplot()  + xlab("Pi")

(a+b+c)/(d+e+f)


ggsave("fis_by_pi.pdf", units="in",height=7,width=10)


#FIS genomic landscape , supplemental figure

FISdat$MB <- FISdat$BP/1000000


ggplot(FISdat, aes(x = MB, y = Fis)) + geom_point(alpha=0.12, size=0.25,colour="black")  + stat_smooth(size=0.5, se=FALSE) + facet_rep_grid(species ~ Chr) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face="italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab(expression(italic(F[IS]))) + theme(plot.title = element_text(size=12)) + scale_x_continuous(breaks=c(0,10,20),minor_breaks=c(5,15)) + scale_y_continuous(limits=c(-0.76,1.1),breaks=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))



#ggsave("eleg_inop_fis_genomic_landscape.pdf", units="in",height=7,width=10)
#ggsave("eleg_inop_fis_genomic_landscape.png", units="in",bg="white",height=7,width=10)



#supplemental figure, sina plots for FIS and Pi by chromosome

a <- ggplot(inop_merge, aes(x=Chr,y=pi_all))+ geom_sina(size=0.2, alpha=0.2,scale="width") + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + theme_cowplot() +ylim(0,0.06)

b <- ggplot(eleg_merge, aes(x=Chr,y=pi_all))+ geom_sina(size=0.2, alpha=0.2,scale="width") + stat_summary(aes(group=Chr),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + theme_cowplot()  +ylim(0,0.06)

c <- ggplot(remdat, aes(x=chrom,y=pi)) + geom_sina(size=0.2, alpha=0.2,scale="width") + stat_summary(aes(group=chrom),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + theme_cowplot()  +ylim(0,0.06)

d <- ggplot(inop_merge, aes(x=Chr,y=Fis))+ geom_sina(size=0.2, alpha=0.2,scale="width") + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + theme_cowplot()  +ylim(-0.5,1)


e <- ggplot(eleg_merge, aes(x=Chr,y=Fis))+ geom_sina(size=0.2, alpha=0.2,scale="width") + stat_summary(aes(group=Chr),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + theme_cowplot()  +ylim(-0.5,1)


f<- ggplot(remFISdat, aes(x=chrom,y=Fis)) + geom_sina(size=0.2, alpha=0.2,scale="width") + stat_summary(aes(group=chrom),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + theme_cowplot()  +ylim(-0.5,1)

(a+b+c)/(d+e+f)


#ggsave("pi_fis_chrom_sina.pdf", units="in",height=7,width=10)
#ggsave("pi_fis_chrom_sina.png", units="in",bg="white",height=7,width=10)

