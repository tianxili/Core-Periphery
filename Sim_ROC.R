library(irlba)
library(igraph)
library(randnet)

source('Helper_functions.R')
source('Methods.R')


####################################
# Simulation Parameters
####################################
n1 = 1000 # core size
n2 = 1000 # periphery size
K = 10 # Number of repetitions to average

sim.type = 'configuration'
# sim.type = 'constant'

lambda = 0.02
tt = c(-0.4, -0.2, 0)

P1 = lambda*(1-tt) # Core density
P2 = lambda*(1+tt) # Periphery density

epsilon = .Machine$double.eps


####################################
# Simulation
####################################

fun_idxs = c(1, 2, 3)
fun_idx = 0
for(g in c(g1, g2, g3)){
  fun_idx = fun_idx+1
  print(paste('g',fun_idxs[fun_idx],sep=''))
  
  if(sim.type=='constant'){
    pdf(paste('g',fun_idxs[fun_idx],'_ROC','_N1=',n1,'_N2=',n2,'_Const.pdf', sep=''), width=12, height=2.8)
  }else{
    pdf(paste('g',fun_idxs[fun_idx],'_ROC','_N1=',n1,'_N2=',n2,'_Config.pdf', sep=''), width=12, height=2.8)
  }
  layout(matrix(1:(length(P1)+1), nrow=1, byrow=TRUE))
  
  par(mar = c(0.5, 0.5, 0.5, 1.5))
  x = seq(0.00001, 1-0.00001, length.out = n1)
  P0 = g(x)
  P0 = P0*0.025/mean(P0)
  P0 = P0[seq(1, length(x), length.out=100), seq(1, length(x), length.out=100)]
  image(-P0[,ncol(P0):1], col=hcl.colors(500), xaxt= "n", yaxt= "n")
  
  
  for(n in 1:length(P1)){
    mean.p.c = P1[n]
    mean.p.p = P2[n]
    
    retval.degree = matrix(0, n1+n2+1, 2)
    retval.pagerank = matrix(0, n1+n2+1, 2)
    retval.ours = matrix(0, n1+n2+1, 2)
    retval.locc = matrix(0, n1+n2+1, 2)
    retval.submat_local = matrix(0, n1+n2+1, 2)
    
    ours.clu = c(0, 0)
    oursk.clu = c(0, 0)
    ASE.clu = c(0, 0)
    K.core.clu = matrix(0, 2*mean.p.c*(n1+n2), 2)
    
    for(k in 1:K){
      
      truth = c(rep(TRUE, n1), rep(FALSE, n2))
      
      x = runif(n1)
      P0 = g(x)
      print(mean(P0))
      P0 = P0*mean.p.c/mean(P0)
      P = add.periphery(P0, n2, p=mean.p.p, type=sim.type)
      rm(P0)
      
      diag(P) = 0
      P[P<0] = 0
      P[P>1] = 1
      #A = P
      A = P.to.A(P)
      # A = Matrix(A)
      
      degree.score.A = rowSums(A)
      
      
      d_select_vec = rep(NA, 5)
      for(i_rank in 1:length(d_select_vec)){
        d_select_vec[i_rank] = ECV.Rank(A, ceiling((n1+n2)^(1/3)) , weighted=FALSE, mode="undirected")$auc.rank
      }
      print(d_select_vec)
      d_select = getmode(d_select_vec)
      # d_select = d_select_vec[fun_idx]
      print(paste('Rank =', d_select))
      
      
      # Scores
      pagerank.score.A = pagerank(A)
      
      # g = graph.adjacency(A, mode='undirected', diag=FALSE)
      locc.score.A = transitivity(graph.adjacency(A, mode='undirected', diag=FALSE), type='localundirected', isolates='zero')
      
      svd.fit.A = irlba(A, nv=d_select+1, maxit = 10000)
      
      submat_local.score = abs(svd.fit.A$u[,1])
      
      A.embed = svd.fit.A$u[,1:d_select,drop=FALSE] %*% diag(sqrt(svd.fit.A$d[1:d_select]), nrow=d_select)
      
      
      if(sim.type=='constant'){
        A_recon.score = spectral.filter.score(A, d_select)
      }else{
        A_recon.score = Ours.Config.function(A, d_select)
      }
      
      
      retval.degree = retval.degree + ROC(degree.score.A, truth)$FPR_TPR
      retval.pagerank = retval.pagerank + ROC(pagerank.score.A, truth)$FPR_TPR
      retval.ours = retval.ours + ROC(A_recon.score, truth)$FPR_TPR
      retval.locc = retval.locc + ROC(locc.score.A, truth)$FPR_TPR
      retval.submat_local = retval.submat_local + ROC(submat_local.score, truth)$FPR_TPR

      
      if(sim.type=='constant'){
        rho_hat = sum(A)/nrow(A)/(nrow(A)-1)
        threshold.value = sqrt( rho_hat^(0.99) * log(nrow(A)) )
      }else{
        rho_hat = sum(A)/nrow(A)/(nrow(A)-1)
        threshold.value = sqrt( log(nrow(A)) / rho_hat^(1.01) ) / nrow(A)
      }
      
      ours.clu = ours.clu + clu_FPR_TPR(A_recon.score > threshold.value, truth)
      oursk.clu = oursk.clu + clu_FPR_TPR(kmeans(log(A_recon.score+epsilon), centers=2, iter.max=20)$cluster, truth)
      ASE.clu = ASE.clu + clu_FPR_TPR(kmeans(A.embed, centers=2, iter.max=20)$cluster, truth)
      
      for(i in 1:nrow(K.core.clu)){
        k.core.score = k_core(A, i)
        K.core.clu[i,] = K.core.clu[i,] + clu_FPR_TPR(k.core.score, truth)
      }
      
    }
    
    retval.degree = retval.degree/K
    retval.pagerank = retval.pagerank/K
    retval.ours = retval.ours/K
    retval.locc = retval.locc/K
    retval.submat_local = retval.submat_local/K

    ASE.clu = ASE.clu/K
    ours.clu = ours.clu/K
    oursk.clu = oursk.clu/K
    K.core.clu = K.core.clu/K
    K.core.clu = K.core.clu[apply(K.core.clu, 1, sum)!=2,]
    
    
    par(mar = c(3.1, 6, 2.1, 0.5))
    if(TRUE){
      plot(c(0,1), c(0,1), type='n', xlab='', ylab='', main=bquote(paste('d'['core'], ' = ',.(round(mean(P[1:n1,])*(n1+n2))),' / ','d'['periphery'], ' = ',.(round(mean(P[(n1+1):(n1+n2),])*(n1+n2))))))
    }else{
      plot(c(0,1), c(0,1), type='n', xlab='', ylab='')
    }
    title(xlab='False Positive', ylab='True Positive', line=2, cex.lab=1)
    
    lines(retval.degree, col=1, lwd=2, lty=2)
    lines(retval.pagerank, col=3, lwd=2, lty=3)
    lines(retval.submat_local, col=4, lwd=2, lty=4) #
    lines(retval.locc, col=5, lwd=2, lty=5) #
    lines(retval.ours, col='red', lwd=2, lty=1) #
    
    points(K.core.clu, col='Brown', pch=2, lwd=2, cex=1.5)
    points(ASE.clu[1], ASE.clu[2], col=7, pch=6, lwd=2, cex=1.5)
    points(ours.clu[1], ours.clu[2], col='red', pch=8, lwd=2, cex=1.5) #
    points(oursk.clu[1], oursk.clu[2], col='red', pch=3, lwd=2, cex=1.5) #
    
    legend('bottomright', legend=c('Degree', 'PageRank', 'EigenVec', 'Local CC', 'ASE', 'k-core', 'Ours'), col=c(1, 3, 4, 5, 7, 'Brown', 'red'), lty=c(2, 3, 4, 5, 0, 0, 1), pch=c(-1, -1, -1, -1, 6, 2, -1), lwd=2, cex=0.82, bg='white')
    
  }
  
  dev.off()
}











