library(quantreg)

# NEEDS THIS FUNCTION
ranks_ext = function(v, score = "wilcoxon", tau = 0.5, trim = NULL)
{
  if(score %in% c('wilcoxon', 'normal', 'sign', 'normalscale', 'halfnormalscale', 'lehmann', 'interquartile'))
  {
    return(ranks(v=v, score=score, tau=tau, trim=trim))
  }
  
  else if(score=="invLehmann")
  {
    J <- ncol(v$sol)
    taus <- v$sol[1, ]
    dt <- taus[2:J] - taus[1:(J - 1)]
    taus <- taus[2:(J - 1)]
    phi <- c(0, taus*log(taus), 0)
    dphi <- diff(phi)
    ranks <- (((v$dsol[, 2:J] - v$dsol[, 1:(J - 1)])))
    ranks <- as.vector(ranks %*% (dphi/dt))
    return(list(ranks = ranks, A2 = 1))
  }
}

# ccomb = function(x) pcauchy(mean(tan((0.5-x)*pi)), lower.tail = FALSE)


myfun = function(Ystar, Gstar, B)
{
  n = length(Gstar)
  Gstar = matrix(Gstar, ncol=1)
  R = matrix(c(1,0,0,2), nrow=2)
  S = t(Gstar) %*% B
  M = R * (sum(Gstar^2))
  M_inv = chol2inv(chol(M))
  T_rs = S %*% M_inv %*% t(S)
  pGREAT = pchisq(T_rs, length(S), lower.tail=FALSE)
  return(pGREAT)
}

# myfun = function(Ystar, Gstar, B)
# {
#   n = length(Gstar)
#   Gstar = matrix(Gstar, ncol=1)
#   M = sum(Gstar^2)/n
#   S = t(Gstar) %*% B/sqrt(n)
#   T_rs = (S^2/M) / rank_w$A2
#   pGREAT = pchisq(T_rs, 1, lower.tail=FALSE)
#   pGREAT = c(pGREAT, ccomb(pGREAT))
#   return(pGREAT)
# }








#rank_weight = c('wilcoxon', 'normal', 'lehmann', 'invLehmann')


# X=covar[,c("Sex","PC1","PC2","PC3","genocohort2","genocohort3","genocohort4","drug_resistant")]
# Y=data.frame(Y=covar$BIS_irnt)
# dati=data.frame(Y, X)
# dati=na.omit(dati)
# 
# 
# # NULL MODEL:
# # HERE, CHANGE THE COVARIATES
# mod0 = rq(paste0("Y~",paste0(colnames(X),collapse="+")), tau=-1, data=dati)
# A = NULL
# B = matrix(NA, nrow=nrow(dati), ncol=length(rank_weight))
# 
# for(w in 1:length(rank_weight))
# {
#   rank_w = ranks_ext(mod0, score=rank_weight[w])
#   B[,w] = rank_w$ranks
#   A = c(A, rank_w$A2)
# }
# 
# Ystar = lm(as.matrix(Y)~as.matrix(X))$residuals
# 
# for(i in 1:nsnp)
# {
#   G=data.frame(G=as.numeric(t(g)))
#   Gstar = lm(as.matrix(G)~as.matrix(X))$residuals
#   pval = myfun(Ystar, Gstar, B)
#   names(pval) <- c(rank_weight,"ccomb")
# }
# 
