
args = commandArgs(trailingOnly=TRUE)

gt_dat <- read.table(args[1], sep="\t",header=TRUE)

dp_dat <- read.table(args[2], sep="\t",fill=TRUE,header=TRUE)


new_dp_mat <-NULL

for(i in 1:nrow(gt_dat)) {
	new_dp_row <- NULL
	new_dp_row_count <- 1
	gt_row_count <- 0
    for (j in 1:24){
		gt_row_count <- gt_row_count + 1
		new_dp_row <- cbind(new_dp_row,ifelse(gt_dat[i,gt_row_count] != "./.:.",dp_dat[i,new_dp_row_count],0))
		new_dp_row_count <- ifelse(gt_dat[i,gt_row_count] != "./.:.",new_dp_row_count+1,new_dp_row_count)}
	new_dp_mat <- rbind(new_dp_mat,new_dp_row)
}

write.table(new_dp_mat,file = args[3],row.names=FALSE,quote=FALSE,sep="\t")
