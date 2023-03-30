###################################################
### chunk number 1: defnormmat
###################################################
#line 11 "mds.Rnw"
normmat <- function(mat){
  matmax <- apply(mat,2,max)
  nmat <- -matrix(rep(matmax,dim(mat)[1]),nrow=dim(mat)[1],byrow=T)+mat
  nmat
}


###################################################
### chunk number 2: defreadneut
###################################################
#line 21 "mds.Rnw"
readneut <- function(dat,logtransform=TRUE){
  ncol <- dim(dat)[2]
  nrow <- dim(dat)[1]
  dunmat <- matrix(0,nrow,ncol)
  cenmat <- matrix(0,nrow,ncol)
  ltmat <-  matrix(0,nrow,ncol)
  gtmat <- matrix(0,nrow,ncol)
  for(i in 1:ncol){
    thiscol <- as.character(dat[,i])
    thiscen <- rep(1,length(thiscol))
    lt <- grep("<",thiscol)
    thiscen[lt] <- 0
    gt <- grep(">",thiscol)
    thiscen[gt] <- 2
    if(logtransform){
      thisneut <- log2(as.double(gsub("<|>","",thiscol)))
    }
    else{
      thisneut <- as.double(gsub("<|>","",thiscol))
    }
    if(length(lt)>0){
      minlt <- min(thisneut[lt])
      ltvec <- rep(minlt,nrow)
      ltvec[lt] <- thisneut[lt]
    }
    else{
      ltvec <- rep(NA,nrow)
    }
    if(length(gt)>0){
      maxgt <- max(thisneut[gt])
      gtvec <- rep(maxgt,nrow)
      gtvec[lt] <- thisneut[gt]
    }
    else{
      gtvec <- rep(NA,nrow)
    }
    ltmat[,i] <- ltvec
    gtmat[,i] <- gtvec
    dunmat[,i] <- thisneut
    cenmat[,i] <- thiscen
  }
  duncmat <- dunmat
  duncmat[cenmat!=1] <- NA
  list(d=dunmat,dc=duncmat,cen=cenmat,lt=ltmat,gt=gtmat)
}


###################################################
### chunk number 3: defprepdata
###################################################
#line 71 "mds.Rnw"
prepdata <- function(datlist){
  N1 <- dim(datlist$dc)[1] # rows, no of viruses
  N2 <- dim(datlist$dc)[2] # cols, no of sera
  rw <- rep(seq(1,N1),N2)
  cl <- rep(seq(1,N2),each=N1)
  dc <- as.vector(datlist$dc)
  cen <- as.vector(datlist$cen)
  lt <- as.vector(datlist$lt)
  gt <- as.vector(datlist$gt)
  list(d=dc,lod=cbind(lt,gt),cen=cen,rw=rw,cl=cl,M=N1*N2,N1=N1,N2=N2)
}


###################################################
### chunk number 4: defprepdatamcmcglmm
###################################################
#line 87 "mds.Rnw"
prepdata.mcmcglmm <- function(datlist){
  N1 <- dim(datlist$dc)[1] # rows, no of viruses
  N2 <- dim(datlist$dc)[2] # cols, no of sera
  rw <- rep(seq(1,N1),N2)
  cl <- rep(seq(1,N2),each=N1)
  cen <- as.vector(datlist$cen)
  lt <- as.vector(datlist$lt)
  gt <- as.vector(datlist$gt)
  d1 <- as.vector(datlist$d)
  d1[cen==0] <- -Inf
  d2 <- as.vector(datlist$d)
  data.frame(d1=d1,d2=d2,rw=rw,cl=cl)
}


###################################################
### chunk number 5: defmapdatamds
###################################################
#line 105 "mds.Rnw"
mapdata.mds <- function(datlist,output,idx){
  N <- datlist$N1
  M <- datlist$N2
  r <- 2
  nidx <- length(idx)
  xmat <- array(0,dim=c(N+M,r,nidx))
  for(i in 1:nidx){
      n <- idx[i]
      x1 <- as.matrix(output$x1)[n,]
      y1 <- as.matrix(output$y1)[n,]
      x2 <- as.matrix(output$x2)[n,]
      y2 <- as.matrix(output$y2)[n,]
      xmat[1:N,1,i] <- x1
      xmat[(N+1):(N+M),1,i] <- y1
      xmat[1:N,2,i] <- x2
      xmat[(N+1):(N+M),2,i] <- y2
  }
  xmat
}


###################################################
### chunk number 6: defmapdatapca
###################################################
#line 129 "mds.Rnw"
mapdata.pca <- function(datlist,output,idx){
  N <- datlist$N1
  M <- datlist$N2
  r <- 2
  nidx <- length(idx)
  xmat <- array(0,dim=c(N,r,nidx))
  for(i in 1:nidx){
      n <- idx[i]
      x <- as.matrix(output$x)[n,]
      xmat[1:N,1:r,i] <- x
  }
  xmat
}


###################################################
### chunk number 7: defbmdsnonorm
###################################################
#line 149 "mds.Rnw"
bmds.nonorm.bug <- "
model{
  for(i in 1:M){
       cen[i] ~ dinterval(d[i],lod[i,1:2])
       d[i] ~ dnorm(delta[i],tau.d)T(0,)
       delta[i] <- sqrt(pow(x1[rw[i]],2)+pow(x2[rw[i]],2)+pow(y1[cl[i]],2)+pow(y2[cl[i]],2)-2*(x1[rw[i]]*y1[cl[i]]+x2[rw[i]]*y2[cl[i]]))
  }
  for(i in 1:N1){
       x1[i] ~ dnorm(0,tau.x0[1])
       x2[i] ~ dnorm(0,tau.x0[2])
  }
  for(i in 1:N2){
       y1[i] ~ dnorm(0,tau.x0[1])
       y2[i] ~ dnorm(0,tau.x0[2])
  }
  for(i in 1:2){
       tau.x0[i] ~ dgamma(0.01,0.01)
  }
  tau.x[1:2] <- sort(tau.x0)
  tau.d ~ dgamma(0.01,0.01)
}
"
write(bmds.nonorm.bug,'bmds.nonorm.bug')


###################################################
### chunk number 8: defbmdsnonormunif
###################################################
#line 177 "mds.Rnw"
bmds.nonorm.unif.bug <- "
model{
  for(i in 1:M){
       cen[i] ~ dinterval(d[i],lod[i,1:2])
       d[i] ~ dnorm(delta[i],tau.d)T(0,)
       delta[i] <- sqrt(pow(x1[rw[i]],2)+pow(x2[rw[i]],2)+pow(y1[cl[i]],2)+pow(y2[cl[i]],2)-2*(x1[rw[i]]*y1[cl[i]]+x2[rw[i]]*y2[cl[i]]))
  }
  for(i in 1:N1){
       x1[i] ~ dunif(-z,z)
       x2[i] ~ dunif(-z,z)
  }
  for(i in 1:N2){
       y1[i] ~ dunif(-z,z)
       y2[i] ~ dunif(-z,z)
  }
  tau.d ~ dgamma(0.01,0.01)
}
"
write(bmds.nonorm.unif.bug,'bmds.nonorm.unif.bug')


###################################################
### chunk number 9: defbpca
###################################################
#line 201 "mds.Rnw"
bpca.bug <- '
model{
  for(i in 1:M){
    cen[i] ~ dinterval(d[i],lod[i,1:2])
    d[i] ~ dnorm(wx[rw[i],cl[i]]+mu[cl[i]],tau.y)
  }
  wx <- x%*%t(W)
  for(j in 1:N2){
    mu[j] ~ dnorm(0,0.001)
  }
  for(i in 1:N1){
    for(j in 1:2){
       x[i,j] ~ dnorm(0,1)
    }
  }
  for(i in 1:N2){
    for(j in 1:2){
      W[i,j] ~ dnorm(0,tau.W)
    }
  }
  tau.y ~ dgamma(0.01,0.01)
  tau.W ~ dgamma(0.01,0.01)
}
'
write(bpca.bug,'bpca.bug')


