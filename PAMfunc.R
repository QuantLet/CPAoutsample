# -------------------------------------------------------------------------------
# One-step SCAD with mutliplier bootstrap
# -------------------------------------------------------------------------------

# Define function for computing derivative of the SCAD penalty
pen.prime = function(beta, lambda, a){
  indicator  = ifelse(abs(beta) <= lambda, 1, 0) 
  tmp        = numeric(0)
  for (jpp in 1:length(beta)){
    tmp[jpp]   = max(a * lambda - abs(beta[jpp]), 0)
  }
  pen.value  = (lambda * indicator) + ((tmp/(a - 1)) * (1 - indicator))
  pen.value
}

# Define function for computing the SCAD penalty
pen.scad = function(beta, lambda, a){
  pen.value = numeric(0)
  for (jps in 1:length(beta)){
    if (abs(beta[jps]) <= lambda){
      pen.value[jps] = lambda * abs(beta[jps])
    } else if (lambda < abs(beta[jps]) && abs(beta[jps]) <= a * lambda){
      pen.value[jps] = - ((abs(beta[jps])^2 - 2 * a * lambda * abs(beta[jps]) + lambda^2)/
                            (2 * (a - 1)))
    } else {
      pen.value[jps] =  ((a + 1) * lambda^2)/2
    }
  }
  pen.value
}

# Define function for submatrix X_U
XU.fct = function(X, beta, lambda, a){
  pen   = pen.prime(beta, lambda, a)       # Evaluate penalty function
  if (min(pen) == 0){
    U.index  = which(pen == 0)
    X.U      = as.matrix(X[, U.index])
  } else {
    X.U      = matrix(c(-Inf, -Inf, -Inf, -Inf), nrow = 2, ncol = 2)
  }
  X.U
}

# Define function for submatrix X_V
XV.fct = function(X, beta, lambda, a){
  pen     = pen.prime(beta, lambda, a)    # Evaluate penalty function
  if (max(pen) == 0){
    X.V      = matrix(c(-Inf, -Inf, -Inf, -Inf), nrow = 2, ncol = 2)
    beta.new = numeric(0)
  } else {
    V.index  = which(pen != 0)
    pen.new  = pen[V.index]
    beta.new = beta[V.index]
    X.V      = as.matrix(X[, V.index])
    for (jxv in 1:length(V.index)){
      X.V[, jxv] = X.V[, jxv] * lambda/ (pen.new[jxv])
    }
  }
  values = list(X.V, beta.new)
  names(values) = c("Xmatrix", "betapar")
  return(values)
}

# Define function for columns' order in  X_U and X_V
order.fct = function(X, beta, lambda, a){
  pen       = pen.prime(beta, lambda, a)       # Evaluate penalty function
  order     = seq(1, dim(X)[2], 1)
  U.index   = which(pen == 0)
  V.index   = which(pen != 0)
  U.order   = order[U.index]
  V.order   = order[V.index]
  order.new = c(U.order, V.order)
  
  order.new
}

# Define function for finding real X.V coefficients
V.coeff.fct = function(coeff, beta, lambda, a){
  pen       = pen.prime(beta, lambda, a) # Evaluate penalty function 
  V.index   = which(pen != 0)
  pen.new   = pen[V.index]
  V.coeff   = numeric(0)
  for (jvc in 1:length(V.index)){
    V.coeff[jvc] = coeff[jvc] * lambda / (pen.new[jvc]) # Transform coeff.'s back
  }
  V.coeff
}

# Define function to select a grid of lambda values
grid.fct = function(n.obs, Y, X){
  
  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X)/n.obs
  X.cent  = scale(X, X.mean, FALSE)
  X.norm  = sqrt(drop(one %*% (X.cent^2)/n.obs))
  X.maxl  = scale(X.cent, FALSE, X.norm)
  
  Y.mean  = drop(one %*% Y)/n.obs
  Y.maxl  = scale(Y, Y.mean, FALSE)
  
  max.lambda.tmp = numeric(0)
  for (j in 1:dim(X.maxl)[2]){
    max.lambda.tmp = c(max.lambda.tmp, t(as.vector(Y.maxl)) %*% X.maxl[, j])
  }
  max.lambda = max.lambda.tmp[which.max(max.lambda.tmp)]
  
  lambda.grid = c(max.lambda/(n.obs^seq(0.49, 0.01, -0.03)*1000)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)*100)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)*10)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03))
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)/10)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)/100))
  lambda.grid
}


# Define function to select a grid of lambda values for the prolonged interval 
# under homogeneity
grid.fct.T = function(n.obs, Y, X){
  
  n.obs.L = length(Y)
  one     = rep(1, n.obs.L)
  X.mean  = drop(one %*% X)/n.obs.L
  X.cent  = scale(X, X.mean, FALSE)
  X.norm  = sqrt(drop(one %*% (X.cent^2)/n.obs.L))
  X.maxl  = scale(X.cent, FALSE, X.norm)
  
  Y.mean  = drop(one %*% Y)/n.obs.L
  Y.maxl  = scale(Y, Y.mean, FALSE)
  
  max.lambda.tmp = numeric(0)
  for (j in 1:dim(X.maxl)[2]){
    max.lambda.tmp = c(max.lambda.tmp, t(as.vector(Y.maxl)) %*% X.maxl[, j])
  }
  max.lambda = max.lambda.tmp[which.max(max.lambda.tmp)]
  
  lambda.grid = c(max.lambda/(n.obs^seq(0.49, 0.01, -0.03)*1000)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)*100)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)*10)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03))
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)/10)
                  , max.lambda/(n.obs^seq(0.49, 0.01, -0.03)/100)) 
  lambda.grid
}

# Define function to compute One-step SCAD algorithm, Li & Zou (2008)
Onestep.SCAD = function(X, Y, a, n.obs){

  beta.tmp    = list()
  act.set     = list()
  a0.tmp      = list()
  
  lambda.grid = grid.fct(n.obs, Y, X)
  
  # Normalize design X (mean = 0, var = 1) and response Y (mean = 0)
  X.orig  = X
  Y.orig  = Y
  
  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X)/n.obs
  X.cent  = scale(X, X.mean, FALSE)
  X.norm0 = sqrt(drop(one %*% (X.cent^2)/n.obs))
  X.norm  = scale(X.cent, FALSE, X.norm0)
  
  X.norm2 = scale(X.orig, FALSE, X.norm0)
  
  Y.mean  = drop(one %*% Y)/n.obs
  Y.norm  = scale(Y, Y.mean, FALSE)
  
  X       = X.norm 
  Y       = Y.norm
  
  # 0. STEP: Find initial values of coefficients - unpenalized fit
  object.ols    = glmnet(X.norm2, Y.orig, family = "gaussian", alpha = 0, lambda = 0,
                         standardize = FALSE, intercept = TRUE)

  # Beta with normalisation - to be used in the next step
  beta.fit.ols  = as.vector(object.ols$beta)
  
  # Original unpenalized coefficients
  beta.0        = as.vector(object.ols$beta)/X.norm0
  a.0           = object.ols$a0
  
  mu0.hat       = X.norm %*% beta.fit.ols
  
  # 1. STEP: Create working data
  # 1.a)
  Y.star = mu0.hat   # Create response according to the inital fit (not centralised)
  X.star = as.matrix(X.norm) # Scaled but nor centralised matrix X 
  
  # Compute the best couple (lambda, beta(lambda)) for the algorithm
  for (llam in 1:length(lambda.grid)){
    
    lambda1  = lambda.grid[llam]
    
    # 1.b)
    X.U      = XU.fct(X.star, beta.fit.ols, lambda1, a)
    XV.tmp   = XV.fct(X.star, beta.fit.ols, lambda1, a)
    X.V      = XV.tmp$Xmatrix
    V.beta   = XV.tmp$betapar
    ord      = order.fct(X.star, beta.fit.ols, lambda1, a)
    
    # 1.c)
    if (X.V[1, 1] != -Inf){
      if (X.U[1, 1] != -Inf){
        Y.2star      = Y.star - X.U %*% ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star
        X.V2star     = X.V - X.U %*% ginv(t(X.U) %*% X.U) %*% t(X.U) %*% X.V 
      } else {
        Y.2star      = Y.star
        X.V2star     = X.V
      }
      
      # 2. STEP: Lasso estimation using CV
      if (dim(X.V)[2] == 1){
        V.coeff.tmp   = ifelse(abs(V.beta) <= 2 * lambda1, 
                               sign(V.beta) * max((abs(V.beta) - lambda1), 0), 
                               ((a - 1) * V.beta - sign(V.beta) * a * lambda1)/(a - 2))
        
        V.coeff       = V.coeff.tmp
        
      } else {
        object        = glmnet(X.V2star, Y.2star, family = "gaussian", alpha = 1, 
                               lambda = lambda1, standardize = FALSE, intercept = FALSE) 
        
        V.coeff.tmp   = as.vector(object$beta)  # Extract coefficients from the fit
        V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit.ols, lambda1, a)
      }
      
      # 3. STEP: Coefficients associated with X_U and X_V
      if (X.U[1, 1] != -Inf){
        U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) 
                            %*% (Y.star - X.V %*% V.coeff.tmp)) 
      } else {
        U.coeff = numeric(0)
      }
      
      coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
      coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
      
      act.set.tmp   = sum(coeff != 0) # No. of nonzero coefficients
      beta.fit      = coeff/X.norm0          
    } else if (X.V[1, 1] == -Inf){
      # 2. STEP: Skip - no V coefficients
      # 3. STEP: Coefficients associated with X_U and X_V
      U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
      V.coeff = numeric(0)
      
      coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
      coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
      
      act.set.tmp   = sum(coeff != 0) # No. of nonzero coefficients
      beta.fit      = coeff/X.norm0                   
    }
    
    act.set[[llam]]     = act.set.tmp
    beta.tmp[[llam]]    = beta.fit
    a0.tmp[[llam]]      = Y.mean - (X.mean %*% beta.tmp[[llam]])
  }
  
  bic    = numeric(0)
  for (lbic in 1:length(lambda.grid)){
    bic[lbic] = (log(n.obs) * act.set[[lbic]]/n.obs * max(log(log(n.par)), 
                                                          sqrt(n.obs)/n.par)
                  + log((t(Y.orig - rep(a0.tmp[[lbic]], n.obs) - 
                             X.orig %*% beta.tmp[[lbic]]) %*% 
                           (Y.orig - rep(a0.tmp[[lbic]], n.obs) 
                            - X.orig %*% beta.tmp[[lbic]]))/n.obs))
  }
  
  index         = which.min(bic)
  
  values        = list(a.0, a0.tmp[[index]], beta.0, beta.tmp[[index]], 
                       lambda.grid[index], act.set[[index]], index, bic[index])
  names(values) = c("a.0", "a", "beta.0", "beta", "lambda", "act.set", "index", "bic")
  return(values)
}

# Define function to compute One-step SCAD algorithm combined with multiplier bootstrap
Onestep.SCAD.MB = function(X, Y, a, n.obs, lambda.value, multipl){
  
  beta.tmp = list()
  act.set  = list()
  
  # Normalize design X (mean = 0, var = 1) and response Y (mean = 0)
  X.orig  = X
  Y.orig  = Y
  
  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X)/n.obs
  X.cent  = scale(X, X.mean, FALSE)
  X.norm0 = sqrt(drop(one %*% (X.cent^2)/n.obs))
  X.norm  = scale(X.cent, FALSE, X.norm0)
  
  Y.mean  = drop(one %*% Y)/n.obs
  Y.norm  = scale(Y, Y.mean, FALSE)
  
  X.mb.mean.vector = t((t(X.orig) %*% multipl))/sum(multipl)
  X.mb.mean.matrix = matrix(rep(X.mb.mean.vector, n.obs), 
                            ncol= length(X.mb.mean.vector), byrow = TRUE)
  
  X.mb       = (X.orig - X.mb.mean.matrix) * sqrt(multipl)
  Y.mb       = (Y.orig - ((t(Y.orig) %*% multipl)/sum(multipl))) * sqrt(multipl)
  
  X.mbnorm   = scale(X.mb, FALSE, X.norm0)
  
  object.ols = glmnet(X.mbnorm, Y.mb, family = "gaussian", alpha = 0, lambda = 0,
                       standardize = FALSE, intercept = FALSE)
  
  # Beta with normalisation - to be used in the next step
  beta.fit.ols  = as.vector(object.ols$beta)
  
  # Original OLS coefficients
  beta.0        = as.vector(object.ols$beta)/X.norm0
  a.0           = (((t(Y.orig) %*% multipl)/sum(multipl)) 
                   - ((t((t(X.orig) %*% multipl)) %*% beta.0)/sum(multipl)))
  
  mu0.hat       = X.mbnorm %*% beta.fit.ols # Initial fit
  
  # 1. STEP: Create working data
  # 1.a)
  Y.star   = mu0.hat   # Create response according to the inital fit
  X.star   = X.mbnorm
  
  lambda   = lambda.value
  
  # 1.b)
  X.U      = XU.fct(X.star, beta.fit.ols, lambda, a)
  XV.tmp   = XV.fct(X.star, beta.fit.ols, lambda, a)
  X.V      = XV.tmp$Xmatrix
  V.beta   = XV.tmp$betapar
  ord      = order.fct(X.star, beta.fit.ols, lambda, a)
  
  # 1.c)
  if (X.V[1, 1] != -Inf){
    if (X.U[1, 1] != -Inf){
      Y.2star      = Y.star - X.U %*% ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star
      X.V2star     = X.V - X.U %*% ginv(t(X.U) %*% X.U) %*% t(X.U) %*% X.V 
    } else {
      Y.2star      = Y.star
      X.V2star     = X.V
    }
    
    # 2. STEP: Lasso estimation using CV
    if (dim(X.V)[2] == 1){
      V.coeff.tmp   = ifelse(abs(V.beta) <= 2 * lambda, 
                             sign(V.beta) * max((abs(V.beta) - lambda), 0), 
                             ((a - 1) * V.beta - sign(V.beta) * a * lambda)/(a - 2))
      
      V.coeff       = V.coeff.tmp
      
    } else {
      object        = glmnet(X.V2star, Y.2star, family = "gaussian", alpha = 1, 
                             lambda = lambda.value, standardize = FALSE, 
                             intercept = FALSE)
      
      V.coeff.tmp   = as.vector(object$beta)  # Extract coefficients from the fit
      V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit.ols, lambda, a)
    }
    
    # 3. STEP: Coefficients associated with X_U and X_V
    if (X.U[1, 1] != -Inf){
      U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U)
                          %*% (Y.star - X.V %*% V.coeff.tmp)) 
    } else {
      U.coeff = numeric(0)
    }
    
    coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
    coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
    
    act.set       = sum(coeff != 0) # No. of nonzero coefficients
    beta.fit      = coeff/X.norm0                 
  } else if (X.V[1, 1] == -Inf){
    # 2. STEP: Skip - no V coefficients
    # 3. STEP: Coefficients associated with X_U and X_V
    U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
    V.coeff = numeric(0)
    
    coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
    coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
    
    act.set       = sum(coeff != 0) # No. of nonzero coefficients
    beta.fit      = coeff/X.norm0                   
  }
  
  lambda.fit      = lambda.value
  
  a0.tmp          = (((t(Y.orig) %*% multipl)/sum(multipl))
                     - ((t((t(X.orig) %*% multipl)) %*% beta.fit)/sum(multipl)))
  
  bic    = (log(n.obs) * act.set/n.obs * max(log(log(n.par)), sqrt(n.obs)/n.par)
            + log((t(Y.mb - X.mb %*% beta.fit) %*% (Y.mb - X.mb %*% beta.fit))/n.obs))
  
  values        = list(a.0, a0.tmp, beta.0, beta.fit, lambda.fit, act.set, Y.mean, bic)
  names(values) = c("a.0", "a", "beta.0", "beta", "lambda", "act.set", "intercept", "bic")
  return(values)
}

# Define function to compute One-step SCAD algorithm for the prolonged interval 
# under homogeneity
Onestep.SCAD.MB.shift = function(X.tmp.L, Y.tmp.L, X.tmp.R, Y.tmp.R, a, lambda.L, 
                                 lambda.R, multipl.L, multipl.R, beta.L, beta.R, 
                                 a.L, a.R){
  
  beta.tmp  = list()
  n.obs.L   = length(Y.tmp.L)
  n.obs.R   = length(Y.tmp.R)
  
  # Normalize design X.T (mean = 0, var = 1) and response Y.T (mean = 0)
  Y.orig.T  = c(Y.tmp.L, Y.tmp.R)
  X.orig.T  = rbind(X.tmp.L, X.tmp.R)
  multipl.T = c(multipl.L, multipl.R)
  
  n.obs.T   = n.obs.L + n.obs.R
  one.T     = rep(1, n.obs.T)
  X.mean.T  = drop(one.T %*% X.orig.T)/n.obs.T
  X.cent.T  = scale(X.orig.T, X.mean.T, FALSE)
  X.norm0.T = sqrt(drop(one.T %*% (X.cent.T^2)/n.obs.T))
  X.norm.T  = scale(X.cent.T, FALSE, X.norm0.T)
  
  X.mb.mean.vector = t((t(X.orig.T) %*% multipl.T))/sum(multipl.T)
  X.mb.mean.matrix = matrix(rep(X.mb.mean.vector, n.obs.T), 
                            ncol= length(X.mb.mean.vector), byrow = TRUE)
  
  X.mb.T           = (X.orig.T - X.mb.mean.matrix) * sqrt(multipl.T)
  X.mbnorm.T       = scale(X.mb.T, FALSE, X.norm0.T)
  
  Y.T       = c(Y.tmp.L, (Y.tmp.R - (rep((a.R - a.L), n.obs.R) 
                                     + X.tmp.R %*% (beta.R - beta.L))))
  Y.mb.T    = (Y.T - ((t(Y.T) %*% multipl.T)/sum(multipl.T))) * sqrt(multipl.T)
  
  # 0. STEP: Find initial values of coefficients - penalized fit
  object.ols    = glmnet(X.mbnorm.T, Y.mb.T, family = "gaussian", alpha = 0, lambda = 0,
                         standardize = FALSE, intercept = FALSE)
  
  # Beta with normalisation - to be used in the next step
  beta.fit.ols  = as.vector(object.ols$beta)
  
  # Original OLS coefficients
  beta.0        = as.vector(object.ols$beta)/X.norm0.T
  a.0           = (((t(Y.T) %*% multipl.T)/sum(multipl.T)) 
                   - ((t((t(X.orig.T) %*% multipl.T)) %*% beta.0)/sum(multipl.T)))
  
  mu0.hat       = X.mbnorm.T %*% beta.fit.ols # Initial fit
  
  # 1. STEP: Create working data
  # 1.a)
  Y.star     = mu0.hat   # Create response according to the inital fit
  X.star     = X.mbnorm.T
  
  # Compute the best couple (lambda, beta(lambda)) for the algorithm
  index      = which(grid.fct(length(Y.tmp.L), Y.tmp.L, X.tmp.L) == lambda.L)
  lambda1    = grid.fct.T(length(Y.star), Y.tmp.L, X.tmp.L)[index]
  
  # 1.b)
  X.U        = XU.fct(X.star, beta.fit.ols, lambda1, a)
  XV.tmp     = XV.fct(X.star, beta.fit.ols, lambda1, a)
  X.V        = XV.tmp$Xmatrix
  V.beta     = XV.tmp$betapar
  ord        = order.fct(X.star, beta.fit.ols, lambda1, a)
  
  # 1.c)
  if (X.V[1, 1] != -Inf){
    if (X.U[1, 1] != -Inf){
      Y.2star      = Y.star - X.U %*% ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star
      X.V2star     = X.V - X.U %*% ginv(t(X.U) %*% X.U) %*% t(X.U) %*% X.V 
    } else {
      Y.2star      = Y.star
      X.V2star     = X.V
    }
    
    # 2. STEP: Lasso estimation using CV
    if (dim(X.V)[2] == 1){
      V.coeff.tmp   = ifelse(abs(V.beta) <= 2 * lambda1, 
                             sign(V.beta) * max((abs(V.beta) - lambda1), 0), 
                             ((a - 1) * V.beta - sign(V.beta) * a * lambda1)/(a - 2))
      
      V.coeff       = V.coeff.tmp
      
    } else {
      object        = glmnet(X.V2star, Y.2star, family = "gaussian", alpha = 1, 
                             lambda = lambda1, standardize = FALSE, intercept = FALSE)
      
      V.coeff.tmp   = as.vector(object$beta)  # Extract coefficients from the fit
      V.coeff       = V.coeff.fct(V.coeff.tmp, beta.fit.ols, lambda1, a)
    }
    
    # 3. STEP: Coefficients associated with X_U and X_V
    if (X.U[1, 1] != -Inf){
      U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) 
                          %*% (Y.star - X.V %*% V.coeff.tmp)) 
    } else {
      U.coeff = numeric(0)
    }
    
    coeff.tmp2    = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
    coeff         = coeff.tmp2[2, order(coeff.tmp2[1, ])]
    
    beta.fit      = coeff/X.norm0.T        
  } else if (X.V[1, 1] == -Inf){
    # 2. STEP: Skip - no V coefficients
    # 3. STEP: Coefficients associated with X_U and X_V
    U.coeff = as.vector(ginv(t(X.U) %*% X.U) %*% t(X.U) %*% Y.star) 
    V.coeff = numeric(0)
    
    coeff.tmp2  = rbind(ord, c(U.coeff, V.coeff))           # Fitted coefficients
    coeff       = coeff.tmp2[2, order(coeff.tmp2[1, ])]
    
    beta.fit    = coeff/X.norm0.T                 
  }
  
  lambda.fit    = lambda1
  
  a0.tmp        = (((t(Y.T) %*% multipl.T)/sum(multipl.T)) 
                   - ((t((t(X.orig.T) %*% multipl.T)) %*% beta.fit)/sum(multipl.T)))
  act.set       = sum(beta.fit != 0)
  
  values        = list(a.0, a0.tmp, beta.0, beta.fit, lambda.fit, act.set)
  names(values) = c("a.0", "a", "beta.0", "beta", "lambda", "act.set")
  return(values)
}

# Define function to compute bootstrapped penalized log-likelihood function
loglik.pen.MB = function(a0, beta.0, beta, lambda, X, Y, multipl){
  
  n.obs      = length(Y)
  
  one        = rep(1, n.obs)
  X.mean     = drop(one %*% X)/n.obs
  X.cent     = scale(X, X.mean, FALSE)
  X.norm0    = sqrt(drop(one %*% (X.cent^2)/n.obs))
  
  X          = X * sqrt(multipl)
  Y          = (Y - a0) * sqrt(multipl)
  
  var.eps    = 1 
  
  loglik1    = (- (1/(n.obs * 2 * var.eps)  *  t(Y - X %*% beta) %*% (Y - X %*% beta))
                - (t(pen.prime(beta.0 * X.norm0, lambda, a)) %*% abs(beta * X.norm0)))
  
  loglik1
}

# Define function to compute bootstrapped unpenalized log-likelihood function
loglik.MB = function(a0, beta, X, Y, multipl){
  
  n.obs   = length(Y)
  
  X       = X * sqrt(multipl)
  Y       = (Y - a0) * sqrt(multipl)
  
  loglik1 = (- (1/(n.obs * 2)  *  t(Y -  X %*% beta) %*% (Y - X %*% beta)))
  
  loglik1
}

# Define function to compute penalized log-likelihood function
loglik.pen = function(a0, beta.0, beta, lambda, X, Y){
  
  n.obs   = length(Y)
  
  one     = rep(1, n.obs)
  X.mean  = drop(one %*% X)/n.obs
  X.cent  = scale(X, X.mean, FALSE)
  X.norm0 = sqrt(drop(one %*% (X.cent^2)/n.obs))
  
  loglik1    = (- (1/(n.obs * 2)  *  
                    (t(Y - rep(a0, n.obs) -  X %*% beta) 
                   %*% (Y - rep(a0, n.obs) - X %*% beta)))
                - (t(pen.prime(beta.0 * X.norm0, lambda, a)) %*% abs(beta * X.norm0)))
  
  loglik1
}

# Define function to compute unpenalized log-likelihood function
loglik = function(a0, beta, X, Y){
  n.obs   = length(Y)
  
  loglik1 = (- (1/(n.obs * 2 )  *  t(Y - rep(a0, n.obs) 
                                    - (X %*% beta)) %*% (Y - rep(a0, n.obs) 
                                                         - (X %*% beta))))
  
  loglik1
}
