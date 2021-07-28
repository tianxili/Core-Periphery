library(irlba)

Ours.ER.function = function(A, r=3){
  r = max(r, 2)
  svd.fit.A = irlba(A, nv=r+3, nu=r+3, maxit=100000)
  
  Sigma = diag(svd.fit.A$d[1:r], nrow=r)
  u = svd.fit.A$u[,1:r,drop=FALSE]
  v = svd.fit.A$v[,1:r,drop=FALSE]
  
  v_ = t( t(v) - apply(v, 2, mean) )
  
  temp.mat = Sigma %*% t(v_) %*% v_ %*% Sigma
  
  retval = sapply(1:nrow(A), function(i){u[i,,drop=FALSE] %*% temp.mat %*% t(u[i,,drop=FALSE])})
  
  return(sqrt(retval))
}


Ours.Config.function = function(A, d=3){
  d = max(d, 2)
  svd.fit = irlba(A, nv=min(d+3, nrow(A)), maxit = 10000)
  A.recon = svd.fit$u[,1:d,drop=FALSE] %*% diag(svd.fit$d[1:d], nrow=d) %*% t( svd.fit$v[,1:d,drop=FALSE] )
  A.recon.normal = A.recon %*% diag( 1/rowSums(A) )
  apply(A.recon.normal, 1, sd) * sqrt(nrow(A)-1)
}



