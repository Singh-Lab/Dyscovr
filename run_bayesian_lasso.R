############################################################
### Run Bayesian Models
### Written By: Amey Pasarkar, Fall 2022; Modified by Sara Geraghty; Princeton University
############################################################

#library(gglasso)
#library(statmod)
#library(mvtnorm)
#library(MCMCpack)
#library(SuppDists)
#library(mvnfast)
#library(MBSGS)
#library(MASS)
#library(mgcv)
#library(mnormt)
#library(truncnorm)

# Link to Git Repo with Amey's original code: https://github.com/Singh-Lab/BayesianGroupLasso
# This file contains two parts, each with an implementation of a Bayesian group 
# LASSO method. See README.md in above repo for in-depth documentation on these 
# methods. 
# PART 1: BGL Method (Kyung et al.): https://projecteuclid.org/journals/bayesian-analysis/volume-5/issue-2/Penalized-regression-standard-errors-and-Bayesian-lassos/10.1214/10-BA607.full
# PART 2: BGLSS Method (Xu and Ghosh): https://projecteuclid.org/journals/bayesian-analysis/volume-10/issue-4/Bayesian-Variable-Selection-and-Estimation-for-Group-Lasso/10.1214/14-BA929.full

################################################################################
### PART 1: BGL IMPLEMENTATION
################################################################################
#' Implementation of Bayesian Group Lasso from Kyung et. al 2010
#' Original source code is very minimal/did not run/was unreadable, heavy 
#' modifications made
#' Returns: model list object. Attributes are accessed by model$attr 
#' Attributes are:
#' Sparse_Beta: Vector of size p. The final posterior mean estimate for each 
#' regressor.
#' Groups that contain 0 in their 98% confidence interval are set to 0
#' To get the coefficients for a particular group, use 
#' model$Sparse_Beta[,which(groups==k)] to get group k
#' Confidence: Vector of size # of groups. Contains confidence score for each 
#' group P(Beta_k<0). Equivalent to finding tail-size from 0
#' Betas: (niterations-burnin) x p regressosrs matrix. Contains the gibbs 
#' sampling estimate for each regressor at each iteration after convergence. 
#' Each row describes output at an iteration. To get the history of coefficients 
#' for group k, use model$Betas[,which(groups==k)] to get group k.
#' Lambda: (niterations-burnin) length vector. Contains the gibbs sampling 
#' estimate for lambda at each iteration. Use mean(model$Lambda) to get the 
#' exact posterior estimate.
#' Sigma2: (niterations-burnin) length vector. Contains the gibbs sampling 
#' estimate for sigma^2 at each iteration. Use mean(model$Sigma2) to get the 
#' exact posterior estimate. 
#' Tau2: (niteations-burnin) x number of groups matrix. Contains the gibbs 
#' sampling estimate for Tau2 for each group at each iteration
#' Tau2 is a less interpretable parameter, it mainly helps with creating an 
#' efficient mathematical formulation. Each row describes the Tau2 for each 
#' group at a given iteration
#' @param X n (number of datapoints) x p (number of features) matrix. 
#' Note: X MUST BE SCALED AND CENTERED. use function X <- scale(X) to do so
#' @param Y vector of response variable. Y must be centered. use Y<-Y-mean(Y) 
#' to do so
#' @param groups vector of length p. Describes which group each regressor 
#' belongs to.
#' @param scale a TRUE/FALSE value indicating whether or not X and Y should be 
#' scaled and centered
#' @param niterations how long to run the Gibbs sampler
#' @param burnin how many of the first iterations to throw out (effectively by 
#' when do you expect model to converge)
runBGL <- function(X, Y, groups, scale, niterations=10000, burnin=5000){
  
  if(scale) {
    X <- scale(X)
    Y <- Y - mean(Y)
  }
  
  n <- nrow(X)
  p <- ncol(X)
  G <- max(groups)
  
  # Define orthonormal basis transformations
  # See Slide 11 in the othonormal_basis_justification.pdf for brief explanation 
  # of why this is needed
  # This is important for allowing us to update each beta parameter one at a time 
  # At each iteration, we will use these transformation to sample the betas, and 
  # then transform back to the original space
  ei  <- eigen(t(X) %*% X) 
  Pt <- diag(ei$values^(1/2))%*%t(ei$vectors) # used to bring beta vector into the orthonormal basis
  Pt_i <- ei$vectors %*% diag(ei$values^(-1/2)) # brings beta vector out of the basis
  
  nX <- X%*%ei$vectors%*%diag(ei$values^(-1/2)) # nX is our orthonormal design matrix
  
  nXX <- t(nX)%*%nX #this will be an identity matrix. (linear algebra property)
  
  # Initialize the parameters by sampling from the prior
  a=1
  b=0.1
  lambda.g.sq  <- rgamma(1,shape=a, rate=b) #sample lambda^2, these a & b yield relatively flat prior
  sig.sqg   <- runif(1,0.1,10)
  tau.sqg   <- beta.g <- c()
  for (g in 1:G){
    mgg         <- length(which(groups==g))
    tau.g       <- rgamma(1,(mgg+1)/2, rate=lambda.g.sq/2)
    mean.g      <- rep(0,mgg)
    cov.g       <- matrix(0,nrow=mgg, ncol=mgg)
    diag(cov.g) <- sig.sqg*tau.g
    betag       <- mvtnorm::rmvnorm(1,mean=mean.g,sigma=cov.g)
    tau.sqg[g]  <- tau.g
    beta.g      <- c(beta.g,t(betag))
  }
  
  # Define matrices/vectors that will contain the information of each iteration 
  # of our gibbs sampler. Save original prior terms
  p <- ncol(X)
  sigsqg.post <- lambdag.post <- NULL 
  tausqg.post <- rbind( tau.sqg,matrix(rep(NA,niterations*G),ncol=G) )
  beta.g      <- rbind( t(beta.g),matrix(rep(NA,niterations*p),ncol=p) )
  
  # Some pre-loading of matrices:
  Xg_ALL = vector(mode='list', length=G)
  Xg_not_ALL = vector(mode='list', length=G)
  for(g in 1:G){
    Xg_ALL[[g]] <- nX[,which(groups==g)] 
    Xg_not_ALL[[g]] <- nX[,which(groups!=g)] 
  }
  # Gibbs
  print("made it to Gibbs sampling")
  for (M in 1:niterations)  {
    # if(M %% 1000 == 0){
    # print(M)
    # }
    # Full Conditional Posteriors  
    gam       <- c()
    beta.p    <- c()
    beta.norm <- c()
    
    # scale all the beta's to orthonormal basis, get the new beta in this basis, 
    # then send back up
    scale_beta <- Pt %*% beta.g[M,]
    for(g in 1:G){
      Xg           <- Xg_ALL[[g]]#nX[,which(groups==g)] 
      Xg_not       <- Xg_not_ALL[[g]]#nX[,which(groups!=g)] 
      beta_not     <- scale_beta[which(groups!=g)]  # get latest estimates of all 
                                                    # the non-beta terms   
      
      # Define covariance matrix for posterior of group g
      dtau.inv     <- diag(1/tau.sqg[g],nrow=length(which(groups==g)),
                           ncol=length(which(groups==g)))
      # because of standardization, t(Xg)%*% Xg is the identity 
      XXg <- diag(1,nrow=length(which(groups==g)),ncol=length(which(groups==g))) 
      
      cov.be<- sig.sqg * solve(XXg+dtau.inv) 
      
      #Define mean for posterior of group g
      subt <- rep(0,nrow(X))
      subt <- Xg_not %*% beta_not
      mean.be    <- 1/sig.sqg * cov.be %*% t(Xg) %*% (Y - 1/2*subt)
      
      # sample from multivariate normal
      # Original:
      # betag      <- mvtnorm::rmvnorm(1,mean=mean.be,sigma=cov.be)
      # New:
      betag      <- mvnfast::rmvn(1,mu=mean.be,sigma=cov.be)
      beta.p     <- c(beta.p, t(betag))
      
      # Update Tau2 for group g, as described in paper
      repeat{
        gam[g]  <- rinvGauss(1, nu=sqrt(lambda.g.sq*sig.sqg/sum(betag^2)), 
                             lambda=lambda.g.sq)
        if (gam[g] > 0) break    	
      }
      tau.sqg[g] <- 1/gam[g]
      #beta.norm[g] <- sum(betag^2)
    }
    
    # As mentioned earlier, transform betas out of the orthonormal basis
    beta.g[M+1,] <- t(Pt_i %*% beta.p)
    
    # Calculating the norm of each group, used for updating other 
    # parameter posteriors
    beta.norm <- rep(0, G)
    for(i in 1:G){
      beta.norm[g] <- sum(beta.g[M+1,which(groups==i)]^2)
    }
    tausqg.post[M+1,] <- tau.sqg 	
    
    # Sigma2 update, in variable sig.sq 
    sh.sig      <- (n-1+p)/2
    sc.sig      <- 1/2*t(Y-X%*%beta.g[M+1,])%*%(Y-X%*%beta.g[M+1,])+ 
      1/2*sum(beta.norm/tau.sqg)
    sig.sqg     <- rinvgamma(1, shape=sh.sig, scale=sc.sig)
    sigsqg.post <- c(sigsqg.post, sig.sqg)
    
    # Lambda update
    sh.lam       <- (p + G)/2 + a 
    sc.lam       <- 1/2*sum(tau.sqg) + b 
    lambda.g.sq  <- rgamma(1, shape=sh.lam, rate=sc.lam)
    lambdag.post <- c(lambdag.post, lambda.g.sq)
  }
  
  params <- list('Betas'=beta.g[(burnin+1):(1+niterations), ], 
                 'Lambda'=lambdag.post[(burnin+1):(1+niterations)], 
                 'Sigma2'=sigsqg.post[(burnin+1):(1+niterations)], 
                 'Tau2'=tausqg.post[(burnin+1):(1+niterations), ])
  
  # set groups to 0 based on 98% interval
  final_beta <- selectGroups(params$Betas, groups)
  params$Sparse_Beta <- final_beta
  params$Confidence <- getConfidence(params$Betas, groups)
  
  return(params)
  
}

selectGroups <- function(posteriorBeta, groups){
  #print(groups)
  betas <- colMeans(posteriorBeta)
  G <- max(groups)
  group_counts <- rep(0, G)
  for(i in 1:length(betas)){
    g <- groups[i]
    quantiles <- quantile(posteriorBeta[,i], probs=c(0.01, 0.99))
    #print(quantiles)
    if(quantiles[['1%']]<0 && quantiles[['99%']]>0){
      group_counts[g] <- group_counts[g]+1
    }
    if(group_counts[g]==length(which(groups==g))){ 
      betas[which(groups==g)] <- rep(0, length(which(groups==g)))
    }
  }
  return(betas)
}

getConfidence <- function(posteriorBeta, groups){
  # For each group, estimate the probability that a regressor crosses over 0
  # Effectively measuring the size of the tail after 0 for each regressor in a group. 
  # Similar intuition to how p-value works
  
  betas <- colMeans(posteriorBeta)
  
  zero_probs <- rep(0, max(groups))
  for(i in 1:max(groups)){
    means <- betas[which(groups==i)]
    
    # Use two different cases for multivariate case (size of group>1) or 
    # univariate (size of group==1)
    if(length(means)>1){
      #estimate cov from the gibbs sampler
      covariance <- cov(posteriorBeta[,which(groups==i)])
      
      #Define lower and upper bounds of integration for each regressor
      lower <- rep(0, length(means))
      upper <- rep(0, length(means))
      for(j in 1:length(means)){
        if(means[j]<0){
          lower[j]<- 0
          upper[j]<- Inf
        }
        else{
          lower[j]<- -Inf
          upper[j]<- 0
        }
      }
      #Numerically integrate
      zero_probs[i] <- pmvnorm(lower=lower, upper=upper, mean=means, 
                               sigma=covariance,keepAttr=FALSE)
    }
    else{
      variance <- var(posteriorBeta[,which(groups==i)])
      if(means<0){
        # want the probability beta>0 if mean<0, so do 1-pnorm
        zero_probs[i] <- 1-pnorm(0, mean=means, sd=variance^0.5) 
      }
      if(means>0){
        zero_probs[i] <- pnorm(0, mean=means, sd=variance^0.5)
      }
    }
    
    
  }
  return(zero_probs)
}


################################################################################
### PART 2: BGLSS IMPLEMENTATION
################################################################################
#' Function for hyperparameter tuning & training of BGLSS method
#' BGLSS paper provided code that runs well-hyperparameter tuning was added by 
#' Amey, some slight modifications to main BGLSS function were made
#' Returns: model list object. Attributes are accessed by model$attr 
#' Attributes are:
#' pos_median: Vector of size p. The final posterior MEDIAN estimate for each 
#' regressor.
#' Confidence: Vector of size # of groups. Contains confidence score for each 
#' group P(Beta_k=0). Lower is better
#' Betas: (niterations-burnin) x p regressosrs matrix. Contains the gibbs 
#' sampling estimate for each regressor at each iteration after convergence. 
#' Each row describes output at an iteration
#' To get the history of coefficients for group k, use 
#' model$Betas[,which(groups==k)] to get group k
#' Lambda2: Choice of lambda^2 in the model - it is fit via own method so it's 
#' not in the gibbs sampler
#' model_type: What were the hyperparameters chosen
#' @param X n (number of datapoints) x p (number of features) matrix.
#' Note: X MUST BE SCALED AND CENTERED. use function X <- scale(X) to do so
#' @param Y vector of response variable. Y must be centered. use Y<-Y-mean(Y) 
#' to do so
#' @param group_sizes vector of length #number of groups. Describes the size 
#' of each group. This assumes that all columns of a group are together. 
#' I.e. if you say group_sizes=c(2,3,3), group membership would be (1,1,2,2,2,3,3,3)
#' @param scale a TRUE/FALSE value indicating whether or not X and Y should be 
#' scaled and centered
#' @param niterations how long to run the Gibbs sampler
#' @param burnin how many of the first iterations to throw out (effectively by
#' when do you expect model to converge)
#' @param nfolds how many folds to use for doing the sparsity tuning. 3 is 
#' sufficient. Set nfolds=0 if you do not want to do hyperparameter tuning (still 
#' tends to perform well)
#' Note: Other options for running BGLSS include using a pre-specified lambda 
#' from CVglasso. This would always do worse, so I only test two types of 
#' hyperparameter settings
runBGLSS <- function(X, Y, group_sizes, scale, niterations=10000, burnin=5000, 
                     nfolds=3){
  
  set.seed(1)
  
  if(scale) {
    X <- scale(X)
    Y <- Y - mean(Y)
  }
  
  #Create folds
  if(nfolds!=0){
    folds <- cut(seq(1,nrow(X)),breaks=nfolds,labels=FALSE)
    
    best_model <- NA #will hold something about the settings
    best_error <- -1
    groups <- convert_grouptype(group_sizes)
    #Perform nfold cross validation
    #Model1: Defaults. Model2: Adjust prior on pi0 to reflect observed sparsity. 
    print('Tuning Hyperparameters')
    for(modeltype in 1:2){
      error <- 0
      for(i in 1:nfolds){
        #print(i)
        testIndexes <- which(folds==i,arr.ind=TRUE)
        testX <- X[testIndexes, ]
        trainX <- X[-testIndexes, ]
        trainY <- Y[-testIndexes]
        testY <- Y[testIndexes]
        
        gr_cv<-cv.gglasso(trainX, trainY, group=groups, loss=c('ls'), 
                          pred.loss='L2', intercept=FALSE, nfolds=5)
        coefs <- coef(gr_cv$gglasso.fit, s = gr_cv$lambda.1se)
        betas <- coefs[2:length(coefs)]
        sparsity <- calculate_sparsity(betas, groups)
        
        lambda <- gr_cv$lambda.1se
        if(modeltype==1){
          res<-ModBGLSS(trainY, trainX, group_size=group_sizes, 
                        niter=floor(niterations/2), burnin=floor(burnin/2))
        }
        else if(modeltype==2){
          res<-ModBGLSS(trainY, trainX, group_size=group_sizes, a=2*sparsity, 
                        b=2*sparsity,niter=floor(niterations/2), 
                        burnin=floor(burnin/2))
        }
        
        betas <- res$pos_median
        
        error <- error + mean((testY-testX %*% betas)**2)
        
      }
      if(error<best_error || best_error<0){
        best_model <- modeltype
        best_error <- error
      }
      gr_cv<-cv.gglasso(X, Y, group=groups, loss=c('ls'), pred.loss='L2', 
                        intercept=FALSE, nfolds=5)
      coefs <- coef(gr_cv$gglasso.fit, s = gr_cv$lambda.1se)
      betas <- coefs[2:length(coefs)]
      sparsity <- calculate_sparsity(betas, groups)
      lambda <- gr_cv$lambda.1se
      #now train model on the full dataset using given settings
    }
    print(sprintf('Training Main Model: %i', best_model))
    if(best_model==1){
      model <- ModBGLSS(Y, X, group_size=group_sizes, niter=niterations, 
                        burnin=burnin)
    }
    else if(best_model==2){
      model<-ModBGLSS(Y, X, group_size=group_sizes, niter=niterations, 
                      burnin=burnin, a=2*sparsity, b=2*sparsity)
    }
    model$model_type <- best_model
    return(model)
  }
  else{
    model <- ModBGLSS(Y, X, group_size=group_sizes, niter=niterations, 
                      burnin=burnin)
    return(model)
  }
}

ModBGLSS = function(Y, X, niter = 10000, burnin = 5000, group_size, a=1, b=1,
                    num_update = 100, niter.update =100, verbose = FALSE, 
                    alpha=1e-1, gamma=1e-1, pi_prior=TRUE, pi=0.5, 
                    update_tau=TRUE, option.weight.group=FALSE, 
                    option.update="global", lambda2_update=NULL)
{
  ####################################
  # Create and Initialize parameters #
  ####################################
  n = length(Y)
  p = dim(X)[2]
  ngroup = length(group_size)
  # Initialize parameters
  tau2 = rep(1, ngroup)
  sigma2 = 1
  l = rep(0, ngroup)
  beta = vector(mode='list', length=ngroup)
  for(i in 1:ngroup) beta[[i]]=rep(0, group_size[i])
  Z = rep(0, ngroup)
  
  ##################################
  # Compute lambda2 via EM         
  # Authors suggest doing 100 iterations of Gibbs/EM combo to fit lambda2, 
  # then leaving lambda2 fixed for rest of model
  ##################################
  if (update_tau==TRUE){
    #fit_for_lambda2 = BGLSS_EM_lambda(Y, X, num_update = num_update, 
    #                                  niter = niter.update, 
    #                                  group_size = group_size,
    #                                  option.update=option.update,
    #                                  option.weight.group=option.weight.group)
    fit_for_lambda2 = fitlambdaEM(Y, X, num_update = num_update, 
                                  niter = niter.update, group_size = group_size,
                                  option.update=option.update, 
                                  option.weight.group=option.weight.group)
    lambda2 = apply(fit_for_lambda2$lambda2_path,2,tail,1)} else
    {
      lambda2 <- lambda2_update
    }
  
  #print(c("lambda2",option.update))
  
  ###############################
  # avoid duplicate computation #
  ###############################
  
  YtY = t(Y) %*% Y
  XtY = t(X) %*% Y
  XtX = t(X) %*% X
  
  XktY = vector(mode='list', length=ngroup)
  XktXk = vector(mode='list', length=ngroup)
  XktXmk = vector(mode='list', length=ngroup)
  
  begin_idx = 1
  for(i in 1:ngroup)
  {
    end_idx = begin_idx + group_size[i] - 1
    Xk = X[,begin_idx:end_idx]
    XktY[[i]] = t(Xk) %*% Y
    XktXk[[i]] = t(Xk) %*% Xk
    XktXmk[[i]] = t(Xk) %*% X[,-(begin_idx:end_idx)]
    begin_idx = end_idx + 1
  }
  
  #####################
  # The Gibbs Sampler #
  #####################
  
  # burnin
  coef = array(0, dim=c(p, niter-burnin))
  all_lg = array(0, dim=c(ngroup, niter-burnin))
  coef_tau = array(0, dim=c(ngroup, niter))
  for (iter in 1:niter)
  {
    # Update beta's
    for(i in 1:ngroup)
    {
      bmk = c()
      for(j in 1:ngroup)
      {
        if(j!=i) bmk = c(bmk, beta[[j]])
      }
      f1 = XktY[[i]] - XktXmk[[i]] %*% bmk
      f2 = XktXk[[i]]+1/tau2[i]*diag(nrow=group_size[i])
      f2_inverse = solve(f2)
      mu = f2_inverse %*% f1
      ### Main part
      l[i] = pi/(pi+(1-pi)*(tau2[i])^(-group_size[i]/2)*det(f2)^(-1/2)*
                   exp(t(f1)%*%mu/(2*sigma2)))
      maxf <- max(f2)
      trythis <- (-group_size[i]/2)*log(tau2[i]) + (-1/2)*log(det(f2/maxf)) + 
        (-dim(f2)[1]/2)*log(maxf) + t(f1)%*%mu/(2*sigma2)
      l[i] = pi/(pi+(1-pi)*exp(trythis))
      
      if(runif(1)<l[i])
      {
        beta[[i]] = rep(0, group_size[i])
        Z[i] = 0
      }
      else
      {
        beta[[i]] = mvnfast::rmvn(1, mu=mu, sigma=sigma2 * f2_inverse)
        Z[i] = 1
      }
    }
    
    # Update tau2's
    if(update_tau) {
      
      if (option.weight.group== FALSE){
        for(i in 1:ngroup)
        {
          if(Z[i]==0){tau2[i] = rgamma(1, shape=(group_size[i]+1)/2, 
                                       rate=lambda2[i]/2)}
          else{tau2[i] = 1/rig(1, mean=sqrt(lambda2[i]*sigma2/sum(beta[[i]]^2)), 
                               scale = 1/(lambda2[i]))}
        }
      }else{
        for(i in 1:ngroup)
        {
          if(Z[i]==0){tau2[i] = rgamma(1, shape=(group_size[i]+1)/2, 
                                       rate=lambda2[i]*(group_size[i])/2)}
          else{tau2[i] = 1/rig(1, mean=sqrt(lambda2[i]*group_size[i]*
                                              sigma2/sum(beta[[i]]^2)), 
                               scale = 1/(group_size[i]*lambda2[i]))}
          coef_tau[i,iter]  = tau2[i]
          
        }
        
      }
    }
    
    
    # Update sigma2
    s=0
    for(i in 1:ngroup)
    {
      s = s + sum(beta[[i]]^2)/tau2[i]
    }
    beta_vec = c()
    for(j in 1:ngroup) beta_vec = c(beta_vec, beta[[j]])
    if(iter > burnin)
      coef[,iter-burnin] = beta_vec
    all_lg[,iter-burnin] = l # this is the letter l, not the number 1. It's a list of each group's zero probability
    sigma2 = rinvgamma(1, shape=(n-1)/2 + sum(Z*group_size)/2 + alpha,
                       scale=(YtY-2*t(beta_vec)%*%XtY+t(beta_vec)%*%
                                XtX%*%beta_vec+s)/2 + gamma)
    
    # Update pi
    if(pi_prior==TRUE)
      pi = rbeta(1, shape1=a+ngroup-sum(Z), shape2=b+sum(Z))
  }
  
  # output the posterior mean and median as our estimator
  pos_mean = apply(coef, 1, mean)
  pos_median = apply(coef, 1, median)
  
  zeroprobs <- rep(0, ngroup)
  for(g in 1:ngroup){
    zeroprobs[g]=mean(all_lg[g,])
  }
  list(pos_median = pos_median, Betas = coef, Confidence=zeroprobs, 
       Lambda2=lambda2, zeroProb=all_lg)
}

fitlambdaEM <- function (Y, X, num_update = 100, niter = 100, group_size, a = 1, 
                         b = 1, verbose = FALSE, delta = 0.001, alpha = 0.1, 
                         gamma = 0.1, pi_prior = TRUE, pi = 0.5, 
                         option.update = "global", option.weight.group = FALSE) 
{
  n = length(Y)
  p = dim(X)[2]
  ngroup = length(group_size)
  tau2 = rep(1, ngroup)
  sigma2 = 1
  lambda2 = 1
  matlambda2 = rep(1, ngroup)
  lambda2_path = rep(-1, num_update)
  matlambda2_path = matrix(-1, ncol = ngroup, nrow = num_update)
  l = rep(0, ngroup)
  beta = vector(mode = "list", length = ngroup)
  for (i in 1:ngroup) beta[[i]] = rep(0, group_size[i])
  Z = rep(0, ngroup)
  YtY = t(Y) %*% Y
  XtY = t(X) %*% Y
  XtX = t(X) %*% X
  XktY = vector(mode = "list", length = ngroup)
  XktXk = vector(mode = "list", length = ngroup)
  XktXmk = vector(mode = "list", length = ngroup)
  begin_idx = 1
  for (i in 1:ngroup) {
    end_idx = begin_idx + group_size[i] - 1
    Xk = X[, begin_idx:end_idx]
    XktY[[i]] = t(Xk) %*% Y
    XktXk[[i]] = t(Xk) %*% Xk
    XktXmk[[i]] = t(Xk) %*% X[, -(begin_idx:end_idx)]
    begin_idx = end_idx + 1
  }
  for (update in 1:num_update) {
    coef = array(0, dim = c(p, niter))
    tau2_each_update = array(0, dim = c(ngroup, niter))
    for (iter in 1:niter) {
      if (verbose == TRUE) {
        #print(iter)
      }
      for (i in 1:ngroup) {
        bmk = c()
        for (j in 1:ngroup) {
          if (j != i) 
            bmk = c(bmk, beta[[j]])
        }
        f1 = XktY[[i]] - XktXmk[[i]] %*% bmk
        f2 = XktXk[[i]] + 1/tau2[i] * diag(nrow = group_size[i])
        f2_inverse = solve(f2)
        mu = f2_inverse %*% f1
        l[i] = pi/(pi + (1 - pi) * (tau2[i])^(-group_size[i]/2) * 
                     det(f2)^(-1/2) * exp(t(f1) %*% mu/(2 * sigma2)))
        maxf <- max(f2)
        trythis <- (-group_size[i]/2) * log(tau2[i]) + 
          (-1/2) * log(det(f2/maxf)) + (-dim(f2)[1]/2) * 
          log(maxf) + t(f1) %*% mu/(2 * sigma2)
        l[i] = pi/(pi + (1 - pi) * exp(trythis))
        if (runif(1) < l[i]) {
          beta[[i]] = rep(0, group_size[i])
          Z[i] = 0
        }
        else {
          beta[[i]] = mvnfast::rmvn(1, mu=mu, sigma=sigma2 * f2_inverse)
          Z[i] = 1
        }
      }
      if (option.weight.group == FALSE) {
        for (i in 1:ngroup) {
          if (Z[i] == 0) {
            tau2[i] = rgamma(1, shape = (group_size[i] + 
                                           1)/2, rate = matlambda2[i]/2)
          }
          else {
            tau2[i] = 1/rig(1, mean = sqrt(matlambda2[i] * 
                                             sigma2/sum(beta[[i]]^2)), 
                            scale = 1/(matlambda2[i]))
          }
        }
      }
      else {
        for (i in 1:ngroup) {
          if (Z[i] == 0) {
            tau2[i] = rgamma(1, shape = (group_size[i] + 
                                           1)/2, rate = matlambda2[i] * 
                               (group_size[i])/2)
          }
          else {
            tau2[i] = 1/rig(1, mean = sqrt(matlambda2[i] * 
                                             group_size[i] * sigma2/sum(beta[[i]]^2)), 
                            scale = 1/(group_size[i] * matlambda2[i]))
          }
        }
      }
      tau2_each_update[, iter] = tau2
      s = 0
      for (i in 1:ngroup) {
        s = s + sum(beta[[i]]^2)/tau2[i]
      }
      beta_vec = c()
      for (j in 1:ngroup) beta_vec = c(beta_vec, beta[[j]])
      coef[, iter] = beta_vec
      sigma2 = rinvgamma(1, shape = (n - 1)/2 + 
                           sum(Z * group_size)/2 + alpha, 
                         scale = (YtY - 2 * t(beta_vec) %*% XtY + t(beta_vec) %*% 
                                    XtX %*% beta_vec + s)/2 + gamma)
      if (pi_prior == TRUE) 
        pi = rbeta(1, shape1 = a + ngroup - sum(Z), 
                   shape2 = b + sum(Z))
    }
    tau2_mean = apply(tau2_each_update, 1, mean)
    matlambda2 = (group_size + 1)/(tau2_mean * group_size)
    lambda2 = (p + ngroup)/sum(tau2_mean * group_size)
    if (option.update == "global") 
      matlambda2 <- rep(lambda2, ngroup)
    if (option.weight.group == FALSE) 
      matlambda2 <- rep((p + ngroup)/sum(tau2_mean), ngroup)
    matlambda2_path[update, ] = matlambda2
  }
  pos_mean = apply(coef, 1, mean)
  pos_median = apply(coef, 1, median)
  list(pos_mean = pos_mean, pos_median = pos_median, coef = coef, 
       lambda2_path = matlambda2_path)
}

# HELPER FUNCTIONS
calculate_sparsity <- function(betas, group_assignments){
  # given a list of betas and their respective group assignments, determine how 
  # many groups are 0.
  groups <- max(group_assignments)
  sparsity <- 0
  # the entire group should be pushed to 0, so just look to see if the group 
  # contains a 0
  for(i in 1:groups){ 
    if(0 %in% betas[which(group_assignments==i)]){
      sparsity <- sparsity + 1
    }
  }
  return(sparsity/groups)
}
convert_grouptype <- function(group_sizes){
  # convert a list of group sizes to a list of group assignments by indices 
  # (what gglasso uses)
  groups <- vector(mode='numeric', length=sum(group_sizes))
  group_num <- 1
  index<-1
  for(j in group_sizes){
    for(i in index:(index+j-1)){
      groups[i]<-group_num
    }
    index <- index+j
    group_num<-group_num+1
  }
  return(groups)
}