library(SiZer)
rr2 <- read.table("original_data_100grains/rarf100grPPE_18sites.txt", check.names = F, stringsAsFactors = F)
ma2 <- rr2[,6:ncol(rr2)]
ma2[ma2>0] <- 1
di2 <- data.frame(rr2[,c(2,4)],ric=rowSums(ma2))
di2 <- di2[order(di2$e., di2$age),]

si <- read.table("original_data_100grains/sites.txt", sep=",", header=T, stringsAsFactors = F)
si <- si[!(si$e. %in% c(1093, 209, 1177, 301, 1103, 1195, 1082)),]
Encoding(si$sigle) <-"UTF-8"
ec <- unique(si$e.)

tb <-di2
i=16
smm <- 100 
mx <- max(tb[,3])
i=1
tb <- tb[tb$age<=11900,]

par(mfrow=c(3,8), mar=c(2,4,2.5,0.5), oma=c(3,3,1,1))
for(i in 1:NROW(ec))
{
vib <- tb[tb$e.==ec[i],]
y=vib$ric
x=vib$age
assign(si[i,3], SiZer(x, y, h=c(500,5000), x.grid = seq(-100, 11900, 300),  derv=1))
}


i=1
tabl <- data.frame( si[i,3],get(si[i,3])$h.grid,   get(si[i,3])$slope)
for(i in 2:nrow(si)){
  tabl <- rbind(tabl,data.frame( si[i,3],get(si[i,3])$h.grid,   get(si[i,3])$slope))
}

colnames(tabl) <- c("file", "bandwidth",get(si[i,3])$x.grid)
write.table(tabl, "ready_data/SiZer_output.csv", sep=";", row.names = F, quote=F)
