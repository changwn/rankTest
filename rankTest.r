

rankTest <- function(data_ttt, tMarker=NULL, cytoMarker=NULL, window_size=100){

#load("C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_try_TIMER/data/RNAseq/coadRNAseq.RData")

T_markers<-c("CD3E" ,"CD2","CD3G", "ITK" , "CD3D","SIRPG","CD48","CD6" ,"TIGIT") #R1
Cyto_markers<-c("PDCD1","CD8A","SLA2","NKG7","PRF1","GZMA","GZMH","GZMB") #R2

if(tMarker == NULL){tMarker <- T_markers}
if(cytoMarker == NULL){cytoMarker <- Cyto_markers}

data1 <- data_ttt[tMarker, ]
data2 <- data_ttt[cytoMarker, ]
data_bind <- rbind(data1, data2)

ddd1 <- svd(data1)
v <- ddd1$v

#use v to sort data_bind
#cor(t(data1), v[, 1])
v_base <- -v[, 1]

position <- order(v_base)

#source("C:/Users/wnchang/Documents/F/PhD_Research/2018_05_07_CIBERSORT_BCV/BCV_v2_20180518.R")
library(bcv)
source("./BCV_ttest2.R")


#-----------------------------
#title : bcv test based immune cell activity prediction
#-----------------------------

N <- ncol(data_bind)
rank_storage <- matrix(NA, N-window_size+1, 4)	#4 column to storage data
for(i in 1:(N - window_size + 1) ){
	left <- i
	right <- i + window_size - 1
	colNum <- position[left:right]
	tmp <- data_bind[,colNum]

	#method 1 : myself
	#hhhh <- BCV_v2(tmp, 2, 2, maxrank=5, norm_type="F") #version2
		#hhhh$msep
		#print.bcv_my(hhhh)
	#rank_i <- hhhh$min_rank
	#rank_i <- bcv_Ttest(hhhh, maxrank = 5)
	#rank_storage[i, 1] <- rank_i

	#method 2: from bcv paper
	#bcv_result1 <- cv.svd.gabriel(tmp, 2, 2, maxrank = 5)
	#rank_ii <- summary(bcv_result1)$rank.1se	#rank.best  do not work
	#rank_storage[i, 2] <- rank_ii	

	#method 3: from Zhang
	#vv <- c()
	#for(j in 1:10){
	#	pp <- BCV_ttest2(tmp, rounds=50)
	#	vv <- rbind(vv, pp)		#10 times ttest for smoothing curve 
	#}
	#bb <- apply(vv, 2, median)

	bb <- BCV_ttest2(tmp, rounds=50)
	rank_storage[i, 3] <- bb[1]
	rank_storage[i, 4] <- bb[2]
}

pdf("Rank_test(cytotoxicity_infiltrating_level).pdf")
par(mfrow = c(2,1))

#matplot(rank_storage, type='l', col=c("grey","blue") ,main="rank curve as load increase")
#matplot(rank_storage[, 1], type='l', col=c("red") ,main="(my-ttest)rank curve as load increase")
#matplot(rank_storage[, 2], type='l', col=c("blue") ,main="(original-bcv)rank curve as load increase")
#rank_storage
matplot(rank_storage[, 3], type='l', col=c("red") ,main="p value of rank1")
matplot(rank_storage[, 4], type='l', col=c("blue") ,main="p value of rank2")

dev.off()

}

