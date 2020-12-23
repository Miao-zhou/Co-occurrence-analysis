

Co-occurrence analysis
==========
It is a probabilistic model based on the Hypergeometric distribution  for detecting statistically significant pairwise co-occurred MeSH terms.

## Authors

[周晓北] (Zhou Xiaobei)  
[崔雷] (Cui Lei)  
[黄德生] (Huang Desheng)    
[周淼]  (Zhou Miao)

## Download data and code
```
 $ git clone https://github.com/Miao-zhou/Co-occurrence-analysis.git
```

The detailed methods are introduced in our paper.

## Analogue simulation 
We set up a simulation to reflect the reality of MeSH term co-occurrence data and we simulate different variables. Detailed information are elaborated in our article.
The code is following to realize this process:

```r
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
```
The function "simu()" is in "code.R".

- We calculated type I and II errors to evaluate the performance of our model:



"![Image text](https://raw.githubusercontent.com/Miao-zhou/Co-occurrence-analysis/main/simulation%20result.png)"


## Real data application
We chose a bibliometric study about pelvic organ prolapse (POP) as our real data application to illustrate how the probabilistic model can be used for MeSH word co-occurrence analysis and compared our new result with  the original one.

### Visualizing the result 

**Visualizing p-value matrix**

We get a high dimensional p-value matrix (3192×3192) produced by the model, then we visualized metrics for clearly interpreting this result.
The visualizing code are following:

```r
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

png("fig1.png",w=2000,h=2000)
s1 <- s
s1[s<=0.95] <- 0.95
corrplot(s1,diag=F,p.mat=1-s,type="upper",method="circle",tl.srt=45,tl.col=1,tl.cex=1.5,cl.length=5,cl.lim=c(0.95,1),is.corr=FALSE,cl.cex=2,pch.cex=4,pch.col="green",insig="label_sig",sig.level=1e-16)
dev.off()
```
ALL function are defined in the file "code.R".


"![Image text](https://raw.githubusercontent.com/Miao-zhou/Co-occurrence-analysis/main/fig1.png)"


**Network of MeSH terms**
- Then, we built a network structure of MeSH terms (0 < pval < 0.05) of POP dataset  to explore the relationship among these MeSH terms.

The code are following:

```r
png("fig2.png",w=2000,h=2000)
set.seed(100)
g <- graph.adjacency(s, mode = "undirected", weighted =T, diag = F)
E(g)$width <- as.numeric(cut(E(g)$weight,c(0,0.95,0.99,1,2),right=FALSE))
E(g)$width[E(g)$width==1] <- 0.75
E(g)$width[E(g)$width==2] <- 1.5
E(g)$color <- ifelse(E(g)$width==3,"blue","gray")
plot(g,vertex.size=4,vertex.label.cex=2)
dev.off()
```


"![Image text](https://raw.githubusercontent.com/Miao-zhou/Co-occurrence-analysis/main/fig2.png)"


**Comparison of results**

We compare our results of our method with the results of original method. Then, we visualize the results of comparison.

```r
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
```

We can get the result like this:


"![Image text](https://raw.githubusercontent.com/Miao-zhou/Co-occurrence-analysis/main/venn.png)"




