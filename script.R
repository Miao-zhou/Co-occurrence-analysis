##"pelvic organ prolapse"[Mesh] AND ("2007/01/01"[PDAT] : "2016/12/31"[PDAT])##

library(pubMR)
library(data.table)
library(tidyr)
library(corrplot)
library(igraph)

dir <- "/Users/xizhou/Nutstore\ Files/Nutstore/mydoc/paper_threshold/data/"
setwd(dir)
obj <- AB(input="zuo.xml")
obj1=data.table(PMID=obj@PMID,MS=obj@MS)
MS <- obj1[,MS]
idx <- sapply(MS,is.null)
obj1 <- obj1[!idx,]
obj1 = obj1 %>% unnest(MS) 
v <- table(obj1[,c("MS","PMID")])
v <- crossprod(t(v))
v1 <- v
diag(v1) <- NA
p=hyp(v,length(obj@PMID))
diag(p) <- NA
idx <- which(rowSums(p==0,na.rm=TRUE)>0)
vr <- v1[idx,idx]
rownames(vr)

pr <- p[idx,idx]
pr[is.na(pr)] <- 0
s <- 1-pr
s1 <- s
diag(s1) <- 0
pdf("co.pdf",w=10,h=10)
corrplot(s,diag=F,type="upper",tl.srt=45,tl.col=1,tl.cex=0.5,cl.lim=c(0,1))
dev.off()


g <- graph.adjacency(s, mode = "undirected", weighted =T, diag = F)
E(g)$width <- as.numeric(cut(E(g)$weight,4))



#png("fig1.png",w=2000,h=1000)
#par(mfrow=c(1,2))
#set.seed(100)
#corrplot(s,diag=F,type="upper",tl.srt=45,tl.col=1,tl.cex=0.5,cl.lim=c(0,1))
#text(-2,65,"a", cex=3)
#g <- graph.adjacency(s, mode = "undirected", weighted =T, diag = F)
#E(g)$width <- as.numeric(cut(E(g)$weight,4))
#E(g)$width[E(g)$width==1] <- 0.25
#E(g)$width[E(g)$width==2] <- 1
#plot(g)
#text(-1.3,1,"b", cex=3)
#dev.off()


png("fig1.png",w=2000,h=2000)
s1 <- s
s1[s<=0.95] <- 0.95
corrplot(s1,diag=F,p.mat=1-s,type="upper",method="circle",tl.srt=45,tl.col=1,tl.cex=1.5,cl.length=5,cl.lim=c(0.95,1),is.corr=FALSE,cl.cex=2,pch.cex=4,pch.col="green",insig="label_sig",sig.level=1e-16)
dev.off()
#text(-2,65,"a", cex=3)
#g <- graph.adjacency(s, mode = "undirected", weighted =T, diag = F)
#E(g)$width <- as.numeric(cut(E(g)$weight,4))
#E(g)$width[E(g)$width==1] <- 0.25
#E(g)$width[E(g)$width==2] <- 1
#plot(g)
#text(-1.3,1,"b", cex=3)
#dev.off()



png("fig2.png",w=2000,h=2000)
set.seed(100)
g <- graph.adjacency(s, mode = "undirected", weighted =T, diag = F)
E(g)$width <- as.numeric(cut(E(g)$weight,c(0,0.95,0.99,1,2),right=FALSE))
E(g)$width[E(g)$width==1] <- 0.75
E(g)$width[E(g)$width==2] <- 1.5
E(g)$color <- ifelse(E(g)$width==3,"blue","gray")
plot(g,vertex.size=4,vertex.label.cex=2)
dev.off()

#s1 <- s
#diag(s1) <- 0
#x1 <- c(rownames(s1)[which(s1==1, arr.ind=T)[,1]],colnames(s1)[which(s1==1, arr.ind=T)[,2]])
#x1 <- unique(x1)
#cat(x1,sep="\n")
#x2 <- c(rownames(s1)[which(s1>0.99&s1<=1, arr.ind=T)[,1]],colnames(s1)[which(s1>0.99&s1<=1, arr.ind=T)[,2]])
#x2 <- unique(x2)
#cat(x2,sep="\n")

idx <- which(rowSums(p==0,na.rm=TRUE)>0)
vr <- v1[idx,idx]
cat(rownames(vr),sep="\n")


f <- fread("/Users/xizhou/Nutstore\ Files/Nutstore/mydoc/paper_threshold/data/venn.csv")
library(VennDiagram)
library(scales)

library(VennDiagram)
venn.plot <- draw.pairwise.venn(
area1 = length(f[!htest=="",htest]), 
area2 = length(f[!zuo=="",zuo]),
cross.area = sum(f[!htest=="",htest]%in%f[!zuo=="",zuo]),
category = c("h-test", "d-method"),
fill = c("blue", "yellow"), 
lty = "blank",
cex = 2, 
cat.cex = 3,
cat.pos = c(190, 0), 
cat.dist = c(0.07,0.06), 
cat.just = list(c(0, 0), c(0, 0)),
ext.pos = 0, 
ext.dist = -0.05,
ext.length = 0.85, 
ext.line.lwd = 2,
ext.line.lty = "dashed",
alpha=0.3,
euler.d=T)

png(filename="venn.png",width=500,height=300)
grid.draw(venn.plot);
dev.off()



