.ghyp <- function(k,n1,n2,n)
{
   max_sp <- max(n1,n2)
   min_sp <- min(n1,n2)
   psite <- as.numeric(n+1)   
   p_gt <- 0
   prob_share_site <- rep(x=0,times=psite)
   all.probs <- phyper(0:min_sp,min_sp,n-min_sp,max_sp)
   prob_share_site[1] <- all.probs[1]
   for (jj in 2:length(all.probs)) 
      prob_share_site[jj] <- all.probs[jj]-all.probs[jj-1]
   for (jj in 0:n) 
   {
      if(jj >= k) 
         p_gt <- prob_share_site[(jj+1)]+ p_gt
   }
   p_gt
}

hyp <- function(x,n)
{
   nr <- nrow(x)
   nc <- ncol(x)
   co <- diag(x)
   psite <- as.numeric(n+1)   
   p <- matrix(0,nrow=nr,ncol=nc)
   colnames(p) <- colnames(x)
   rownames(p) <- rownames(x)
   for(i in seq(nr))
   {   
      for(j in seq(nc))
      {       
         sp1 <- co[i]
         sp2 <- co[j]
         if(x[i,j]==0)
            p[i,j] <- 1
         else
            p[i,j] <- .ghyp(x[i,j],sp1,sp2,n)
      }
      cat("i",i, "\n")  
   }
   p
}


#nt <- 100
#nf <- 200
#nc <- 200
#rho <- runif(nt,min=0.1,max=0.2)
#pf <- runif(nf,min=0.05,max=0.2)
#pt1 <- runif(nt,min=0.05,max=0.2)
#pt2 <- runif(nt,min=0.05,max=0.2)

simu <- function(rho,pf,pt1,pt2,nt,nf,nc)
{
   require(mipfp)
   require(gplots)
   require(data.table)
   require(dendextend)
   nr <- 2*nt+nf
   yt <- list()
   for(i in seq(nt))
   {
      corr <- matrix(c(1,rho[i],rho[i],1),2,2)
      d <- ObtainMultBinaryDist(cor=corr,marg.probs=c(pt1[i],pt2[i]))
      y <- RMultBinary(n=nc,mult.bin.dist=d)
      yt[[i]] <- t(y$binary.sequences)
   }
   yt <- do.call("rbind",yt)
   yf <- list()
   for(i in seq(nf))
   {
      tmp <- rbinom(nc,1,pf[i])
      if(all(tmp==0))
         tmp[sample(nc,1)] <- 1
      yf[[i]] <- tmp  
   }   
   yf <- do.call("rbind",yf)
   y <- rbind(yt,yf)
   idx <- sample(nr,nr)
   y <- y[idx,]
   idx_true <- match(seq(2*nt),idx)
   rownames(y) <- paste0("word_",seq(nr))
   colnames(y) <- paste0("pmid_",seq(nc))
   v <- crossprod(t(y))
   v1 <- v
   diag(v1) <- NA
   x <- v1
   x <- sweep(x,1L,rowMeans(x,na.rm=TRUE),check.margin=FALSE)
   sx <- apply(x,1L,sd,na.rm=TRUE)
   x <- sweep(x,1L,sx,"/",check.margin=FALSE)
   s <- idx_true
   idxr <- s[seq(s)%%2==1]
   idxc <- s[seq(s)%%2==0]
   v2 <- matrix(NA, nrow=nrow(v),ncol=ncol(v))
   for(i in seq(idxr))
      v2[idxr[i],idxc[i]] <- "+"
   for(i in seq(idxr))
      v2[idxc[i],idxr[i]] <- "+"

   vh=hyp(v,nc)
   diag(vh) <- NA

   df <- data.table(pred=as.vector(vh),labels=as.vector(v2))
   df <- df[!is.na(pred),]
   df[,truth:=0]
   df[labels=="+",truth:=1]
   df[,labels:=NULL]
   type1_error <- df[truth==0,mean(pred<0.05)]
   type2_error <- df[truth==1,mean(pred>=0.05)]
   c(type1_error,type2_error)
}