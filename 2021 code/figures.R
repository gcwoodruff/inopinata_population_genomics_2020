
#figures for Woodruff et al. 2021 C. inopinata population genomics paper.


library(ggplot2)
library(lemon)
library(ggforce)
library(patchwork)
library(RColorBrewer)
library(rstatix)

#pi, nucleotide diversity figures

dat <- read.csv("ino_elg_pi.csv", header=TRUE)


levels(dat$species)[levels(dat$species)=="elegans"] <- "C. elegans"
levels(dat$species)[levels(dat$species)=="inopinata"] <- "C. inopinata"


dat$species <-factor(dat$species, levels = c("C. elegans", "C. inopinata"))

dat$MB <- dat$start/1000000


ggplot(dat, aes(x = MB, y = pi_all)) + geom_point(alpha=0.12, size=0.25,colour="black")  + stat_smooth(size=0.5, se=FALSE) + facet_rep_grid(species ~ scaffold) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face="italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab("π (Nucleotide diversity)") + theme(plot.title = element_text(size=12)) + scale_x_continuous(breaks=c(0,10,20),minor_breaks=c(5,15))



#ggsave("inopinata_elegans_pi_genome.pdf", width=6.5, height=4,units = "in", useDingbats=FALSE)  
#ggsave("inopinata_elegans_pi_genome.png", width=6.5, height=4,units = "in")  


#normalize by chromosome position

ele <- dat[dat$species == "C. elegans",]
ino <- dat[dat$species == "C. inopinata",]

ele_I <- ele[ele$scaffold == "I",]
ele_II <- ele[ele$scaffold == "II",]
ele_III <- ele[ele$scaffold == "III",]
ele_IV <- ele[ele$scaffold == "IV",]
ele_V <- ele[ele$scaffold == "V",]
ele_X <- ele[ele$scaffold == "X",]

ele_I$norm_dist_center <- abs((7536217-ele_I$start)/7536217)/2
ele_II$norm_dist_center <- abs((7639711-ele_II$start)/7639711)/2
ele_III$norm_dist_center <- abs((6891901-ele_III$start)/6891901)/2
ele_IV$norm_dist_center <- abs((8746915-ele_IV$start)/8746915)/2
ele_V$norm_dist_center <- abs((10462090-ele_V$start)/10462090)/2
ele_X$norm_dist_center <- abs((8859471-ele_X$start)/8859471)/2

ino_I <- ino[ino$scaffold == "I",]
ino_II <- ino[ino$scaffold == "II",]
ino_III <- ino[ino$scaffold == "III",]
ino_IV <- ino[ino$scaffold == "IV",]
ino_V <- ino[ino$scaffold == "V",]
ino_X <- ino[ino$scaffold == "X",]

ino_I$norm_dist_center <- abs((10297276-ino_I$start)/10297276)/2
ino_II$norm_dist_center <- abs((10058498-ino_II$start)/10058498)/2
ino_III$norm_dist_center <- abs((9718237-ino_III$start)/9718237)/2
ino_IV$norm_dist_center <- abs((10508822-ino_IV$start)/10508822)/2
ino_V$norm_dist_center <- abs((11819078-ino_V$start)/11819078)/2
ino_X$norm_dist_center <- abs((9095254-ino_X$start)/9095254)/2


pi_dat <- rbind(ele_I,ele_II,ele_III,ele_IV,ele_V,ele_X,ino_I,ino_II,ino_III,ino_IV,ino_V,ino_X)

#define chromosome arms and centers
pi_dat$chr_str_type <- ifelse(pi_dat$norm_dist_center >= 0.25,"Arms", "Centers")



ggplot(pi_dat, aes(x = species, y = pi_all)) + geom_sina(aes(colour=chr_str_type),size=0.2, alpha=0.2,scale="width") + stat_summary(aes(group=chr_str_type),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12, face="italic"), axis.text.y = element_text(colour="black", size=12),legend.text=element_text(colour="black", size=12),legend.title=element_text(colour="black", size=13),axis.title=element_text(size=14), strip.background = element_rect(colour="white", fill="white"), strip.text = element_text(size=12, colour="black"),legend.key=element_blank(),legend.position = "none") + xlab("Species") + ylab("π (Nucleotide diversity)") + labs(colour = "Chromosome\nregion",size=12) + guides(colour = guide_legend(override.aes = list(size=1,alpha=1))) + scale_y_continuous(limits=c(0,0.08),breaks=c(0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08))



#ggsave("inopinata_elegans_pi_sina_no_legend_8-7-20.pdf", width=4, height=3,units = "in", useDingbats=FALSE)  


#ggsave("inopinata_elegans_pi_sina_no_legend_8-7-20.png", width=4, height=3,units = "in")  



#fstat figures



dat <- read.csv("all_fst_combined.csv", header=TRUE)


dat$MB <- dat$BP_start/1000000


dat$pop_pair <- factor(dat$pop_pair, levels = c("ishi-irio","irio-yona","ishi-yona"))





#fis
new_dat <- read.table("fis_elegans_and_inopinata.tsv", header = TRUE, sep = "\t")


new_dat$MB <- new_dat$BP/1000000



levels(new_dat$species)[levels(new_dat$species)=="elegans"] <- "C. elegans"
levels(new_dat$species)[levels(new_dat$species)=="inopinata"] <- "C. inopinata"


ino_fis_dat <- new_dat[new_dat$species == "C. inopinata",]

ele_fis_dat <- new_dat[new_dat$species == "C. elegans",]

ggplot(new_dat, aes(x = MB, y = Fis)) + geom_point(alpha=0.12, size=0.25,colour="black")  + stat_smooth(size=0.5, se=FALSE) + facet_rep_grid(species ~ Chr) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12, face="italic"), strip.background = element_blank(), axis.ticks = element_line(colour = "black")) + xlab("Position (MB)") + ylab(expression(italic(F[IS]))) + theme(plot.title = element_text(size=12)) + scale_x_continuous(breaks=c(0,10,20),minor_breaks=c(5,15)) + scale_y_continuous(limits=c(-0.76,1.1),breaks=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))



#ggsave("inopinata_elegans_fis_genome.pdf", width=6.5, height=4,units = "in", useDingbats=FALSE)  

#ggsave("inopinata_elegans_fis_genome.png", width=6.5, height=4,units = "in")  

#prep data for composite plot


ino_fis_prep <- data.frame(Chr=ino_fis_dat$Chr,BP=ino_fis_dat$BP,F=ino_fis_dat$Fis,MB=ino_fis_dat$MB)
ino_fis_prep$pop_pair <- "fis"

ino_fst_prep <- data.frame(Chr=dat$Chr,BP=dat$BP_start,F=dat$FST,MB=dat$MB,pop_pair=dat$pop_pair)

f_stats <- rbind(ino_fst_prep,ino_fis_prep)



x_labels <- c(expression(paste("Ir-Is ",italic(F[ST]))),expression(paste("Ir-Yo ",italic(F[ST]))),expression(paste("Is-Yo ",italic(F[ST]))),expression(italic(F[IS])))


ggplot(f_stats, aes(x = MB, y = F)) + geom_point(aes(colour=pop_pair),alpha=0.25, size=0.15) + geom_smooth(aes(colour=pop_pair), size=0.75,se=FALSE) + facet_rep_wrap( ~ Chr, nrow=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.title=element_text(size = 12), legend.text = element_text(size = 11), legend.key = element_blank(),legend.position = "none") + xlab("Position (MB)") + ylab(expression(italic(F)))  + scale_colour_manual(values=c("#a6bddb","#2b8cbe","#045a8d","red")) + scale_y_continuous(limits=c(-0.76,1.1),breaks=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)) + scale_x_continuous(breaks=c(0,10,20))



ggsave("f_genome_plot.png", width=7.25, height=4.75, units = "in")


#composite plot




a <- ggplot(f_stats, aes(x = MB, y = F)) + geom_point(aes(colour=pop_pair),alpha=0.25, size=0.15) + geom_smooth(aes(colour=pop_pair), size=0.75,se=FALSE) + facet_rep_wrap( ~ Chr, nrow=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=14), axis.text.x=element_text(colour="black", size=12),axis.text.y=element_text(colour="black", size=12),strip.text.x = element_text(size=12), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.title=element_text(size = 12), legend.text = element_text(size = 11), legend.key = element_blank(),legend.position = "none") + xlab("Position (MB)") + ylab(expression(italic(F)))  + scale_colour_manual(values=c("#a6bddb","#2b8cbe","#045a8d","red")) + scale_y_continuous(limits=c(-0.76,1.1),breaks=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)) + scale_x_continuous(breaks=c(0,10,20))



b <- ggplot(f_stats, aes(x = pop_pair, y = F)) + geom_sina(size=0.15, alpha=0.25,scale="width",aes(colour=pop_pair)) + scale_colour_manual(values=c("#a6bddb","#2b8cbe","#045a8d","red")) + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=13), axis.text.y = element_text(colour="black", size=13),legend.position = "none",axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"), strip.text = element_text(size=13, colour="black")) + xlab(expression(paste(italic(F)," statistic"))) + ylab(expression(italic(F))) + scale_x_discrete(labels= x_labels) + scale_y_continuous(limits=c(-0.76,1.1),breaks=c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1))


a/b

#ggsave("f_composite_plot.png", width=7, height=5.5, units = "in")

#ggsave("f_composite_plot.pdf", width=7, height=5.5, units = "in", useDingbats=FALSE)


#natural history worm histograms


#get data in there; data from 2018 woodruff and phillips 2018 bmc ecology paper
ecol_dat <- read.table("okinawa_septica_stage_worm_presence_data_R_edit_5-3-18.csv", sep="\t", header=T)

#getting data types right
ecol_dat$Latitude <- as.numeric(as.matrix(ecol_dat$Latitude))
ecol_dat$Longitude <- as.numeric(as.matrix(ecol_dat$Longitude))
ecol_dat$Fig_stage <- as.factor(ecol_dat$Fig_stage)
ecol_dat$Plant <- as.factor(ecol_dat$Plant)
ecol_dat$worm_presence <- as.factor(ecol_dat$worm_presence)

#set aside pollinated and unpollinated

ecol_dat$pollination[ecol_dat$Fig_stage != "1"] <- "yes"
ecol_dat$pollination <- as.factor(ecol_dat$pollination)
levels(ecol_dat$pollination) <- c("yes", "no")
ecol_dat$pollination[ecol_dat$Fig_stage != "1"] <- "yes"
ecol_dat$pollination[ecol_dat$Fig_stage == "1"] <- "no"

ecol_dat[,c(8,18)]

pollinated <- ecol_dat[ecol_dat$pollination == "yes",]
unpollinated <- ecol_dat[ecol_dat$pollination == "no",]


#hm, some missing data here, let's assume that if the thing is unpollinated, then there are no foundresses...

ecol_dat_b <- ecol_dat

ecol_dat_b$foundresses[ecol_dat_b$pollination == "no"] <- 0


levels(ecol_dat_b$pollination) <- c("Pollinated figs","Unpollinated figs")

pollinated_figs_foundress_data <- ecol_dat_b[ecol_dat_b$pollination == "Pollinated figs",]$foundresses

worms_per_wasp_data <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,2,2,3,3,5,6)
	#this is actually the data, see 2018 woodruff and phillips 2018 bmc ecology paper


a <- ggplot() + aes(pollinated_figs_foundress_data) + geom_histogram(binwidth=1, fill="#111A1A",colour="black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.position="none", strip.text = element_text(colour="black", size=12), strip.background = element_rect(colour="white", fill="white"),axis.title=element_text(size=13), axis.ticks = element_line(colour = "black")) + xlab("Number of wasps per fig") + ylab("Frequency") + scale_x_continuous(breaks=c(0:11)) + scale_y_continuous(limits=c(0,60),breaks=c(0,10,20,30,40,50,60))


b <- ggplot() + aes(worms_per_wasp_data) + geom_histogram(binwidth=1, fill="#FF1A00",colour="black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=12), axis.text.y = element_text(colour="black", size=12),legend.position="none", strip.text = element_text(colour="black", size=12), strip.background = element_rect(colour="white", fill="white"),axis.title=element_text(size=13), axis.ticks = element_line(colour = "black")) + xlab("Number of worms per wasp") + ylab("Frequency") + scale_x_continuous(breaks=c(0:6)) + scale_y_continuous(limits=c(0,20),breaks=c(0,5,10,15,20))



a/b

#ggsave("natural_history_histograms.png",height=5,width=4,units="in")
#ggsave("natural_history_histograms.pdf",height=5,width=4,units="in")


#comparison of worm F stats with fig wasp F stats


dat <- read.table("bisculatus_data.tsv", header = TRUE, sep = "\t")


ggplot(dat, aes(x = f_stat, y = value)) + geom_hline(yintercept=0.14531,size=0.8,linetype="dashed") + geom_hline(yintercept=0.016015,size=0.8,linetype="dotted") + geom_sina(aes(colour=paper),size=1, alpha=1,scale="width") + stat_summary(aes(group=paper),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="black",position = position_dodge(width = 0.9)) + scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=13, face="italic"), axis.text.y = element_text(colour="black", size=13),legend.position = "none",axis.title=element_text(size=15), strip.background = element_rect(colour="white", fill="white"), strip.text = element_text(size=13, colour="black"),legend.key=element_blank()) + scale_y_continuous(limits=c(-0.1,0.5),breaks=c(-0.1,0,0.1,0.2,0.3,0.4,0.5)) + ylab(expression(italic(F))) + xlab("Pollinating fig wasp F statistic")

#ggsave("bisculatus_sina_no_legend.pdf",height=5,width=4,useDingbats=FALSE)

#supplemental figure , coverage


before_dat <- read.table("all_depth_before_genotype_filter.tsv",sep="\t",header=TRUE)

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

fig_df_melt[is.na(fig_df_melt$value),]

fig_df_ii[is.infinite(fig_df_ii$coverage_per_cut_site),]




ggplot(fig_df_melt, aes(x = MB, y = value)) + geom_point(alpha=0.12, size=0.25,colour="black")  + stat_smooth(size=0.5, se=FALSE) + facet_rep_grid(variable ~ Chr,scales="free_y") + theme_cowplot() + theme(strip.background = element_blank()) + scale_x_continuous(breaks=c(0,10,20),minor_breaks=c(5,15))


#ggsave("coverage.png", width=7, height=6,units = "in")  



#supplemental figure, nucleotide diversity across Caenorhabditis genus


dat<- read.table("pi_estimates_caenorhabditis_with_inopinata_for_presentation_10-16-20.tsv", sep="\t", header=T)


species_order <- c("C. brenneri","C. sinica","C. remanei","C. latens","C. japonica","C. inopinata","C. briggsae","C. elegans","C. tropicalis")

dat$species <- factor(dat$species, levels = c("C. brenneri","C. sinica","C. remanei","C. latens","C. japonica","C. inopinata","C. briggsae","C. elegans","C. tropicalis"))


ggplot(dat,aes(x=reorder(species,-PI), y=PI)) + geom_bar(stat="identity",aes(fill=reproductive_mode)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(colour="black", size=14,face="italic",angle = 45, hjust =1), axis.text.y = element_text(colour="black", size=15),legend.position="none", strip.text.x = element_text(size=15),axis.title=element_text(size=17), strip.background = element_rect(colour="white", fill="white"),axis.ticks = element_line(colour = "black")) + xlab("Species") + ylab("π (nucleotide diversity)") + scale_fill_brewer(breaks=c("gon", "self"), labels=c("Male/Female", "Selfer")) + annotate("segment", x=0.5, xend=6.5, y=.17, yend=.17, arrow=arrow(ends="both", angle=90, length=unit(.2,"cm"))) + annotate("segment", x=6.6, xend=9.6, y=.17, yend=.17, arrow=arrow(ends="both", angle=90, length=unit(.2,"cm"))) + annotate("text", x=3.5, y=0.177, label="Obligate outcrossers", size=6) + annotate("text", x="C. briggsae", y=0.176, label="Selfers", size=6)


#ggsave("caeno_pi_WITH_inopinata.png", width=7, height=6, units = "in")





#Supplemental figure, PCA


pca_df <- read.table("24_inopinata_worms_pca.tsv", header = TRUE, sep = "\t")


#plots


ggplot(pca_df, aes(x = PC1, y = PC2)) + geom_point(aes(colour=Island)) + scale_colour_manual(values=c(c("#111A1A","#FF0000","#005AB3"))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=14, colour="black"),axis.ticks = element_line(colour = "black"), axis.title=element_text(size=16), legend.position = "none") + xlab("PC1 (14%)") + ylab("PC2 (9%)") + scale_x_continuous(limits=c(-110,50),breaks=c(-100,-75,-50,-25,0,25,50))


#ggsave("pca_plot_6-30-21.pdf", width=3.5, height=3,units = "in",useDingbats=FALSE)  




#Supplemental figure, DAPC clustering


plotdf <- read.table("dapc_clusters_for_figure.tsv", header = TRUE, sep = "\t")

plotdf$K <- as.factor(plotdf$K)
plotdf$Group <- as.factor(plotdf$Group)

#best K by BIC = 12
k12 <- plotdf[plotdf$K == "12",]


colourCount = length(unique(k12$K))
getPalette = colorRampPalette(brewer.pal(12, "Set1"))


plotdf$fig_id <- factor(plotdf$fig_id ,levels = c("89-4","89-3","89-2","87-3","72-17","72-16","44-2","44-1","22-2"))

plotdf$sample_id <-factor(plotdf$sample_id, levels = c("may_2017_H11.bam","may_2017_F11.bam","A05.bam","dec_2016_D10.bam","B03.bam","A03.bam","may_2017_D11.bam","may_2017_C11.bam","may_2017_A12.bam","may_2017_E11.bam","dec_2016_D02.bam","dec_2016_D01.bam","may_2017_H10.bam","may_2017_E12.bam","may_2017_D12.bam","may_2017_C12.bam","dec_2016_D05.bam","dec_2016_D04.bam","dec_2016_D03.bam","may_2017_H12.bam","may_2017_G12.bam","may_2017_B12.bam","may_2017_B11.bam","may_2017_A11.bam"))


ggplot(plotdf, aes(x = sample_id, y = Posterior, fill = Group)) + geom_bar(stat = "identity") + scale_fill_manual(values=getPalette(12)) + facet_rep_wrap( ~ K, ncol=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(), axis.title=element_text(size=16), axis.text.x=element_blank(),axis.text.y=element_text(colour="black", size=15),strip.text.x = element_text(size=15), strip.background = element_blank(), axis.ticks = element_line(colour = "black"),legend.position = "none") + xlab("Sample") + ylab("Posterior") +scale_y_continuous(limits=c(0,1),breaks=c(0,1)) 

#ggsave("DAPC_cluster_plot_no_labels.pdf",height=6,width=8,units="in",useDingbats=FALSE)


#supplemental figure for BIC plot following DAPC in file "statistics.R"

