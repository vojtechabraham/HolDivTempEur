library(SiZer)
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

rr2 <- read.table("original_data_100grains/rarf100grPPE_18sites.txt", check.names = F, stringsAsFactors = F)
ma2 <- rr2[,6:ncol(rr2)]
ma2[ma2>0] <- 1
di2 <- data.frame(rr2[,c(2,4)],ric=rowSums(ma2))
di2 <- di2[order(di2$e., di2$age),]
tb <-di2

tabl <- read.table( "ready_data/SiZer_output.csv", sep=";", header=T,check.names = F)
bandy<- read.table("ready_data/sites_order_Fig3.txt", sep="\t", stringsAsFactors = F)
Encoding(bandy$V4) <- "UTF-8"


#tiff("Fig3.tiff",  compression = "lzw", height = 1300, width = 1000, res=300, units = "mm", pointsize = 55)
cairo_pdf("Fig3.pdf",  height = 26, width = 20, pointsize = 30)

nf<- layout(matrix(nrow = 11, ncol = 8, data=c(1:10,45,
                                               1:6,11:14,45,
                                               15:20,11:14,45,
                                               15:20,21:24,45,
                                               25:30,21:24,45,
                                               25:30,31:34,45,
                                               35:40,31:34,45,
                                               35:40,41:44,45)
                   , byrow = F,),
            height = c(rep(1,10), 0.3))

layout.show(nf)
i=1 
for(i in 1:nrow(bandy)){   
  
  #par(mfrow=c(2,1))
  par(mar=c(0,4,3,0.5), mgp=c(2.5, 0.5, 0) )
  b <- bandy[i,]
  mt <- t(tabl[tabl$file==b$V1,3:ncol(tabl)])
  h <- tabl[tabl$file==b$V1,2]
  mt <- mt[nrow(mt):1,]
  grix <- as.integer(colnames(tabl)[3:ncol(tabl)])
  
  image((mt), axes=F, ylab = "resolution" , xlab="", col = c("#d7191c",  "gray", "#0571b0","white"))
  axis(side=2, at=seq(0,1, 1/(ncol(mt)-1))[seq(1,NROW(h),3)], labels =round(h[seq(1,NROW(h),3)]), cex.axis=0.7, las=1 ,  tcl=par("tcl")*0.6)
  abline(h=(1/(NROW(h)-1))*12, lwd=3)
  title(main=b[1,"V4"], line=1, cex.main=1.3)
  vib <- tb[tb$e.==b[1,"V3"],]
  cs <- cbind(rev(grix), mt[,as.integer(b[1,2])])
  y=vib[,"ric"]
  x=vib[,"age"]
  xg <- grix
  rozestup <- (xg[-1]-xg[-NROW(xg)])[1]
  model2 <- locally.weighted.polynomial(x,y,h= (h[as.integer(b[1,2])]))
  par(mar=c(2,4,0,0.5), mgp=c(2.5, 0.5, 0) )
  plot(model2, use.ess = F, xaxs="i",main="", pch=16,ylab="pollen richness", xlab="age BP", cex.axis=0.7, las=1, cex=0.8, ylim=c(7,37), xlim=c(12050, -250), xaxt="n")
  axis(side=2, at=seq(10,35,5), labels =seq(10,35,5), cex.axis=0.7, las=1 ,  tcl=par("tcl")*0.6)
  axis(1, at=seq(0,12000,2000), labels=seq(0,12,2), cex=0.7, tcl=par("tcl")*0.6)
  if(NROW(cs[cs[,2]== -1,1])>0){
    rect(xleft = cs[cs[,2]==-1,1]-(rozestup/2), xright = cs[cs[,2]==-1,1]+(rozestup/2), ybottom = 0, ytop = 100, border = "transparent", col = add.alpha("#d7191c", 0.5))
  }
  if(NROW(cs[cs[,2]== 1,1])>0){
    rect(xleft = cs[cs[,2]== 1,1]-(rozestup/2), xright = cs[cs[,2]== 1,1]+(rozestup/2), ybottom = 0, ytop = 100, border = "transparent", col = add.alpha("#0571b0", 0.5))
  }
  
  if(i %in% c(3,18)){
    plot(1:2,1:2, axes=F, bty="n", col="transparent", xlab="", ylab="")
    plot(1:2,1:2, axes=F, bty="n", col="transparent", xlab="", ylab="")
    plot(1:2,1:2, axes=F, bty="n", col="transparent", xlab="", ylab="")
    plot(1:2,1:2, axes=F, bty="n", col="transparent", xlab="", ylab="")}
  
}


plot(1:2,1:2, axes=F, bty="n", col="transparent", xlab="", ylab="")
text(x=1.5, y=0, labels = "age (kyr BP)", cex = 1.5, xpd=T)
dev.off()






