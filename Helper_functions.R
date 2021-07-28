####################################
# Helper function ROC
####################################

vec_2_norm = function(x){
  sqrt(sum(x^2))
}


pagerank = function(A, eps=1e-08, d=0.85, E=rep(1/nrow(A), nrow(A))){
  v = runif(ncol(A))
  v = v/sum(v)
  last_v = rep(1,ncol(A))*100
  A = as.matrix(A)
  
  for(i in 1:ncol(A)){
    s = sum(A[,i])
    if(s!=0){
      A[,i] = A[,i] / s
    }
  }
  
  while( vec_2_norm(v-last_v) > eps ){
    last_v = v
    v = d* A %*% v + (1-d)* E
  }
  as.numeric(v)
}


add.periphery = function(P, size=nrow(P), p=mean(P), degree.seed=apply(P, 1, sum), type=c('configuration', 'constant')){
  if( type=='constant' ){
    retval = matrix(p, nrow(P)+size, nrow(P)+size)
  }else{
    d1 = apply(P, 1, sum)
    d2 = sample(degree.seed, size, replace=TRUE)
    d2 = d2/mean(d2)*sqrt(p*sum(d1))
    retval = c(d1,d2) %*% t(c(d1,d2)) / sum(d1) # so that the core part has the same degree distributio
  }
  retval[1:nrow(P), 1:nrow(P)] = P
  retval
}


P.to.A = function(P){
  n = nrow(P)
  A = matrix(runif(n^2), n, n) < P
  A[lower.tri(A, diag=TRUE)] = 0
  A = A + t(A)
  return(A)
}


ROC = function(score, pos=c(rep(TRUE,round(length(score)/2)), rep(FALSE,length(score)-round(length(score)/2))), plot=FALSE, correct=FALSE){
  # score: numeric vector, the coreness score of each node
  # pos: boolean vector, the label of each node
  # plot: if true, a ROC curve will be plotted; otherwise, only AUC is returned.
  # correct: if true, the AUC is always greater than 0.5
  score = as.numeric(score)
  n1 = sum(pos)
  n2 = length(pos) - n1
  
  o = order(score, decreasing=TRUE)
  score = score[o]
  pos = pos[o]
  
  TPR = rep(0, length(score)+1)
  FPR = rep(0, length(score)+1)
  
  for(i in 1:length(score)){
    TPR[i+1] = sum(pos[1:i]) / n1
    FPR[i+1] = sum(!pos[1:i]) / n2
  }
  
  auc = sum(diff(FPR) * TPR[-length(TPR)])
  
  if(correct & auc<0.5){
    auc = 1 - auc
    TPR = 1 - TPR[length(TPR):1]
    FPR = 1 - FPR[length(FPR):1]
  }
  
  if(plot){
    plot(FPR, TPR, type='l')
    text(0.75,0.2,paste('AUC =',auc))
  }
  
  list('AUC'=auc, 'FPR_TPR'=cbind(FPR, TPR))
}


clu_FPR_TPR = function(clu, truth){
  # Returns the (False Positive Rate, True Positive Rate) of a clustering
  TPR = sum(truth[clu==clu[1]])/n1
  FPR = sum(!truth[clu==clu[1]])/n2
  if( TPR < FPR ){
    TPR = sum(truth[clu!=clu[1]])/n1
    FPR = sum(!truth[clu!=clu[1]])/n2
  }
  c(FPR, TPR)
}


k_core = function(A, k=3){
  keep.id = rep(TRUE, ncol(A))
  n = -1
  
  while( n != sum(keep.id) ){
    n = sum(keep.id)
    if(n==0){
      break
    }
    keep.id = apply(A[, keep.id,drop=FALSE], 1, sum) >= k
  }
  
  return(keep.id)
}


getmode = function(v){
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

####################################
# Graphon functions
####################################
g1 = function(x, K = floor(log(length(x)))){
  n = length(x)
  l = floor(x*K)
  P = sapply(l, function(y){retval=(l+1)/(K+1); retval[y!=l]=0.3/(K+1); retval})
  P
}

g2 = function(x){
  P = sapply(x, function(y){sin(5*pi*(x+y-1)+1)/2+0.5})
}

g3 = function(x){
  P = sapply(x, function(y){1/(1+exp(15*(0.8*abs(y-x))^(4/5)-0.1))})
}







