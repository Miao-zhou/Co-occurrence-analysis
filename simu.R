dir <- "/Users/xizhou/Nutstore\ Files/paper_threshold/data/"
setwd(dir)
set.seed(100)
nt <- 100
nf <- 200
pf <- runif(nf,min=0.05,max=0.2)
pt1 <- runif(nt,min=0.05,max=0.2)
pt2 <- runif(nt,min=0.05,max=0.2)

rho01 <- rep(0.1,nt)
r200_01 <- simu(rho01,pf,pt1,pt2,nt,nf,200)
r500_01 <- simu(rho01,pf,pt1,pt2,nt,nf,500)
r1000_01 <- simu(rho01,pf,pt1,pt2,nt,nf,1000)
r2000_01 <- simu(rho01,pf,pt1,pt2,nt,nf,2000)

rho02 <- rep(0.2,nt)
r200_02 <- simu(rho02,pf,pt1,pt2,nt,nf,200)
r500_02 <- simu(rho02,pf,pt1,pt2,nt,nf,500)
r1000_02 <- simu(rho02,pf,pt1,pt2,nt,nf,1000)
r2000_02 <- simu(rho02,pf,pt1,pt2,nt,nf,2000)









res1 <- rbind(r200_01,r500_01,r1000_01,r2000_01)
res1 <- cbind(c(200,500,1000,2000),res1)
res2 <- rbind(r200_02,r500_02,r1000_02,r2000_02)
res <- cbind(res1,res2)

rownames(res) <- NULL
colnames(res) <- c("nc","Type I error","Type II error","Type I error","Type II error")
#library(xtable)
library(kableExtra)
kable(res,"latex") %>%
  add_header_above(c(" ", "rho=0.1" = 2, "rho=0.2" = 2))

#library(ROCR)

#df <- data.table(pred=1-as.vector(vh),labels=as.vector(v2))
#df <- df[!is.na(pred),]
#df[,truth:=0]
#df[labels=="+",truth:=1]
#df[,labels:=NULL]
#pred <- prediction(df$pred, df$truth)
#perf <- performance(pred,"tpr","fpr")
#png("roc_sim.png",width=1000,height=1000)
#par(mar=c(5.1,5.1,4.1,2.1))
#par("cex.axis"=2)
#plot(perf,cex.lab=3)
#dev.off()

my_palette <- colorRampPalette(c("gray", "white","orange"))(n = 100)
png("fig3.png",width=2000,height=1500)
par("cex.axis"=3)
heatmap.2(v1,cellnote=v3,notecex=5,notecol="blue",Rowv=FALSE,Colv=FALSE,labRow=FALSE,labCol=FALSE,
   trace="none",key.par=list(mar=c(3.5,0,3,0)),lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(0.5, 5), lwid=c(1, 10, 1),margins=c(3,0),col=my_palette,key.title=NA,density.info="none",key.xlab="")
dev.off()