##############################################
#Multi-graphon estimator by S. Chandna and P.A. Maugis,
#'Nonparametric regression for multiple heterogeneous networks'
#arXiv preprint arXiv:2001.04938
#Real data implementation: Templeton255
##############################################
require(reshape2);
require(ggplot2);
require(grid);
require(loe);
require(vegan);
require(mgcv);
require(Matrix);
require(rgl);
require(GoFKernel)
require(fields)

load("Templeton255.Rdata")###adjacencies and network covariates
m <- length(G); n <- ncol(G[[1]]);
N=m;
zcci=cci/max(cci);##network-level covariate

graph_smpl_dist_base <- function(G){#Vertex distances based on neighborhood similarity
  ## Input:
  # G : a graph sample; i.e., a list of adjacency matrices of simple graphs where vertices can be identified across realization.
  # 
  ## Output:
  # A n xn matrix dst, where n is the order of the entris of G. dst_{ij} is the distance between 
  # vertices i and j computed using G. The metric is that of Airoldi et. al. (2013) simplified
  # to address undirected graphs, and randomized for stability. Assuming G is an iid sample drawn
  # from a fixed graphon f and that each vertex i is associated with a latent x_i, this metric estimates
  # consistentli \int_{[0,1]} (f(x_i,x)-f(x_j,x))dx and therefore measures how the neighborhood of i
  # differs from that of j.
  #
  ## Bibliography:
  # E. M. Airoldi, T. B. Costa and S. H. Chan, (2013) "Stochastic blockmodel approximation of a graphon: Theory and consistent estimation", NIPS
  #
  #
  ## Allocating variables
  m <- length(G);
  N=m;
  n <- ncol(G[[1]]);
  Gg   <- array(unlist(G),c(n,n,m));
  # rijk <- array(0,c(n,n,n));
  #
  ## Core loop
  R <- sample.int(m,floor(m/2)); #Randomization, conrasting with Airoldi et. al. (2013)
  
  #view(R)
  
  rij <- (apply(Gg[,,R],c(1,2),mean)%*%apply(Gg[,,-R],c(1,2),mean))/(n-2)
  # for (i in 1:n){ for (j in 1:n){ for (k in 1:n){ if (!(k %in% c(i,j))){
  # R <- sample.int(N,floor(N/2)); #Randomization, conrasting with Airoldi et. al. (2013)
  # rijk[i,j,k] <- sum(Gg[i,k,R])*sum(Gg[j,k,-R]);
  # }}}};# rijk[i,j,k] corresponds to \hat r_{ij}^k in eq. (4) of Airoldi et. al. (2013)
  # rij <- apply(rijk,c(1,2),sum);
  #
  ## Metric computation
  rhohat <- sum(Gg)/(n*(n-1)*m)
  dij <- (outer(diag(rij),diag(rij),'+')-rij-t(rij))/(rhohat^2);#(N^2*(n-2)); #Only the first 4 terms of eq. (5) of Airoldi et. al. (2013)
  dij[which(dij<0)] <- 0;
  #
  ## Output
  return(list(dij=dij,rij=rij))
}
#
graph_smpl_dist <- function(G,r=20,gram=F,low_rank=0){#Further randomizes vertex distance
  ## Input:
  # G : a graph sample; i.e., a list of adjacency matrices of simple graphs where vertices can be identified across realization.
  # r : number of independent replications
  #
  ## Output:
  # A n xn matrix dst, where n is the order of the entris of G. dst_{ij} is the average distance
  # between vertices i and j computed using graph_smpl_dist_base in r independent replications. 
  #
  #
  ## Main loop (vectorized)
  if (!gram){
    dst <- replicate(r,graph_smpl_dist_base(G)$dij,simplify=F)
  } else {
    dst <- replicate(r,graph_smpl_dist_base(G)$rij,simplify=F)
  }
  #
  ## Allocating and averaging
  dst <- array(unlist(dst),c(ncol(G[[1]]),ncol(G[[1]]),length(G)));
  dst <- apply(dst,c(1,2),mean);
  ## Eigenvalue thresholding
  if (low_rank>0){
    edc <- eigen(dst);
    sbs <- if (gram) {head(1:nrow(dst),low_rank)} else {c(1,tail(1:nrow(dst),low_rank))};
    dst <- edc$vectors[,sbs]%*%diag(edc$values[sbs])%*%t(edc$vectors[,sbs]);
  }
  #
  ## Output
  return(dst)
}

graph_smpl_order <- function(G,type='i',ini='rand',...){# Ordering of vertices using metric or ordinal embedding method
  ## Input:
  # G : a graph sample; i.e., a list of adjacency matrices of simple graphs where vertices can be identified across realization.
  # ... : aditional arguments to pass on to graph_smpl_dist
  #
  # Output:
  # A permutaion perm of the vertices in G such that adjacent vertices in the ordering are close
  # according the the distance matrix computed with graph_smpl_dist. Uses the SOE function from the
  # 'loe' package.
  #
  #
  ## Computing distance matrix dst if not provided
  if (!('dst' %in% names(list(...)))){
    dst  <- graph_smpl_dist(G,...);
  };
  #
  if (type == 'i'){
    ## Calling Isomap function
    perm <- order(scores(isomap(dst,k=ncol(G[[1]]),ndim=1)));
  } else {
    ## Calling SOE function
    comp_dst <- get.order(dst)
    res_SOE <- (SOE(comp_dst,ncol(G[[1]]),p=1,iniX=ini,rnd=nrow(comp_dst)));
    
  }
  #
  ## Output:
  return(res_SOE)
}

graph_smpl_dist_plot <- function(dst){#Plot distance matrix
  ## Input:
  # dst: distance matrix
  #
  ## Output:
  # a ggplot2 object with entries of dst plotted as stat_countour geaom.
  ## Reshaping data
  dst.copy <- melt(dst);
  names(dst.copy) <- c("x","y","z");
  #
  ## Output
  ggplot(dst.copy,aes(x,y))+
    stat_contour(aes(z=z,colour=..level..))+
    scale_x_continuous(expand = c(0, 0), limits = c(0,1))+
    scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
    theme(axis.text = element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          plot.margin = unit(c(0,0,0,-.5),"cm"),
          aspect.ratio=1)
}

gram <- graph_smpl_dist(G,r=150,gram=T);
lat_v <- eigen(gram)$vectors[,1:3]%*%diag(sqrt(eigen(gram)$values[1:3]))


PP=25##user-specified number of repetitions to obtain the minimal stress embedding [Terada and Luxburg, 2014]

perm_pp<-matrix(NA,n,PP);
xihatt_pp<-matrix(NA,n,PP);
stressfn<-matrix(NA,PP,1);

for (p in 1:PP){
  res_SOE<- graph_smpl_order(G,type='j')
  xihatt_pp[,p]=res_SOE$X  
  perm_pp[,p]<-order(res_SOE$X)
  stressfn[p]<-res_SOE$str
  print(stressfn[p])
}

minstress=which.min(stressfn);
xihatt=xihatt_pp[,minstress]
perm<-perm_pp[,minstress];

#save(xihatt,perm, file = "xihatt_Templeton_PP25.RData")
#load("xihatt_Templeton_PP25.RData")

load("permSOEBrain2_PP25_jcgs.RData")

sxihatt=sort(xihatt)

Gshat<-lapply(1:m, matrix, data= NA, nrow=n, ncol=n)
for(i in 1:m)
{
  Gi=G[[i]];
  
  Gshat[[i]]=Gi[sort.list(xihatt),sort.list(xihatt)];
}
image(Gshat[[2]])##visualize adjacency with nodes sorted by xihatt
brain_nodes_sorted=brain_nodes[sort.list(xihatt)]

N=m;

v1<-rep(sxihatt,times=n)
v111<-rep(v1,times=N);

v2<-rep(sxihatt,each=n);
v211<-rep(v2,times=N);


v3<-rep(zcci,each=n*n);
v311<-v3;


ys1<-unlist(Gshat)
y1=cbind(c(ys1));

xxx11_multigr <- bam(y1~s(v111,v211,v311,fx=TRUE),family=gaussian(link = "identity"))
rhoHat=matrix(0,N,1);
for (j in 1:N)
{
  rhoHat[j]=sum((Gshat[[j]]))/(n*n);
}
rhoHat1 =sum(rhoHat)/N


##########################get predictions for unique net-covariate values
lun=length(unique(cci));
indcc=matrix(0,lun,1);
ucci=unique(cci);

for (jj in 1:lun){
  indcc_vec=which(cci==ucci[jj])
  indcc[jj]=indcc_vec[1];
}
# 

preds11_uniP=array(0,dim=c(n,n,lun))

for (jj in 1:lun) {
  loc1_unP<- data.frame(v111=rep(sxihatt,each=n),v211=rep(sxihatt,times=n),v311=zcci[indcc[jj]])#aa instead of 1:n
  preds11_uniP[,,jj] <- matrix(predict.gam(xxx11_multigr,newdata=loc1_unP,type='response'),n)
}

#writeMat("preds_Templeton255.mat",preds11_uniP=preds11_uniP,cci=cci,zcci=zcci,indcc=indcc,rhoHat=rhoHat)
#writeMat("xihatt_Templeton255.mat",xihatt=xihatt,perm=perm)

##############Example multi-graphon plot
zind=2#specify any cci index between 1 to lun
po=2;
image.plot(abs(preds11_uniP[,,zind]/rhoHat1)^(1/po))##multigraphon corresponding to cci[zzind]



#Part II: Graphon estimation assuming iid nets:
Gshatavg=matrix(NA,n,n);
Gs=array(data=NA,dim=c(n,n,N));
preds11_gr=array(0,dim=c(n,n));

for(i in 1:N)
{
  Gs[,,i]=Gshat[[i]];
}
Gshatavg=apply(Gs,c(1,2),mean);

y1<-cbind(c(Gshatavg))


xxx11_gr <- bam(y1~s(v1,v2,fx=TRUE),family=gaussian(link = "identity"))
loc1_gr<- data.frame(v1=rep(sxihatt,each=n),v2=rep(sxihatt,times=n))
preds11_gr <- matrix(predict.gam(xxx11_gr,newdata=loc1_gr,type='response'),n)

image.plot(abs(preds11_gr/rhoHat1))
