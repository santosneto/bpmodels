#' Probability Density Function
#'
#'
#'
#'@export
dBP <- function(x,mu=1,sigma=1,log=FALSE)
{
  if (any(mu < 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(x <= 0))
    stop(paste("x must be positive", "\n", ""))

  a <- mu*(1+sigma)
  b <- 2 + sigma

  fy <- dbetapr(x, shape1 = a, shape2 = b, scale = 1, log = log)
  fy

}
#' Cumulative Distribution Function
#'
#'
#'
#'@export

pBP <-  function(q,mu=1,sigma=1, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu < 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(q < 0))
    stop(paste("q must be positive", "\n", ""))

  a <- mu*(1+sigma)
  b <- 2 + sigma

  cdf <- pbetapr(q, shape1 = a, shape2 = b, scale=1, lower.tail = lower.tail,
                 log.p = log.p)
  cdf
}

#' Random Numbers
#'
#'
#'
#'@export

rBP <- function(n,mu=1,sigma=1)
{
  if (any(mu < 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(n <= 0))
    stop(paste("n must be a positive integer", "\n", ""))

  n <- ceiling(n)

  a <- mu*(1+sigma)
  b <- 2 + sigma

  r <- rbetapr(n,shape1=a,shape2=b,scale=1)

  r
}


#' Quantile Function
#'
#'
#'
#'@export


qBP <- function(p, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu < 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(p <= 0) | any(p >= 1))
    stop(paste("p must be between 0 and 1", "\n", ""))

  a <- mu*(1+sigma)
  b <- 2 + sigma

  q <- qbetapr(p, shape1 = a, shape2 = b,scale=1, lower.tail = lower.tail, log.p = log.p)

  q
}



#' gamlss BP
#'
#'
#'
#'@export
BP <- function (mu.link = "log", sigma.link = "log")
{
  mstats <- checklink("mu.link", "Beta Prime", substitute(mu.link), c("log", "identity", "sqrt"))
  dstats <- checklink("sigma.link", "Beta Prime",substitute(sigma.link), c("log", "identity", "sqrt"))
  structure(list(family = c("BP", "Beta Prime"),
                 parameters = list(mu = TRUE,sigma = TRUE), nopar = 2, type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,


                 dldm = function(y, mu, sigma){
                   a <- mu*(1+sigma)
                   b <- mu*(1+sigma)+sigma+2
                   Phi <-  (1+sigma)
                   yast <- log(y) - log(1+y)
                   muast <- digamma(a) - digamma(b)
                   dldm <- Phi*(yast - muast)

                   dldm
                 },
                 d2ldm2 = function(mu, sigma){
                   Phi2 <- (1+sigma)^2
                   a <- mu*(1+sigma)
                   b <- mu*(1+sigma)+sigma+2
                   d2dldm2 <- -Phi2*(trigamma(a) - trigamma(b))

                   d2dldm2
                 },
                 dldd = function(y, mu, sigma){
                   Phi <-  (1+sigma)
                   a <- mu*(1+sigma)
                   b <- mu*(1+sigma)+sigma+2
                   ystar <- mu*log(y) - (1+mu)*log(1+y)
                   mustar <- mu*digamma(a) - (1+mu)*digamma(b) + digamma(Phi+1)

                   dldd <- ystar - mustar

                   dldd
                 },
                 d2ldd2 = function(mu,sigma){
                   Phi <-  (1+sigma)
                   a <- mu*(1+sigma)
                   b <- mu*(1+sigma)+sigma+2

                   d2ldd2 <- -(mu^2)*trigamma(a) + ((1+mu)^2)*trigamma(b) - trigamma(Phi+1)

                   d2ldd2

                 },
                 d2ldmdd = function(mu,sigma){

                   a <- mu*(1+sigma)
                   b <- mu*(1+sigma)+sigma+2
                   Phi <-  (1+sigma)
                   gammaast <- Phi*(trigamma(b) + mu*(trigamma(b)-trigamma(a)))

                   d2ldmdd <- gammaast

                   d2ldmdd

                 },
                 G.dev.incr = function(y, mu, sigma,...){-2*dBP(y, mu, sigma, log = TRUE)},
                 rqres = expression(rqres(pfun = "pBP", type = "Continuous", y = y, mu = mu, sigma = sigma)),
                 mu.initial = expression({mu <- mean(y)}),
                 sigma.initial = expression({sigma <-  mean(y)*(1+mean(y))/var(y) }),
                 mu.valid = function(mu) all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),
                 y.valid = function(y) all(y > 0)),
            class = c("gamlss.family","family"))
}


#' GEE BP
#'
#'
#'
#'@export
geeBP = function(formula, data, id, tol = 0.001, maxiter = 25, corstr = "independence", linkmu = "log", print=FALSE){
  
  namescor = c("independence", "unstructured", "exchangeable", "AR-1", "one-dependent",
               "one-dependent-stat","two-dependent","two-dependent-stat")
  if(all(namescor != corstr)){
    stop("the correlation structure is not defined")
  }
  nameslink = c("log", "identity")
  if(all(nameslink != linkmu)){
    stop("the link function is not defined")
  }
  
  formula = as.formula(formula)
  call <- match.call()
  X = as.matrix(model.matrix(formula, data = data)) # Matriz de especificação
  p = ncol(X) # Número de parâmetros
  y = model.frame(formula, data = data)[,1] # Variável resposta
  nr = table(id) # Número de repetições
  n = max(id) # Número de unidades experimentais
  N = nrow(X)
  
  if(linkmu == "log"){
    mod0 = gamlss(formula,family = BP(mu.link = "log"), trace = FALSE, data = data)
  }
  else{
    mod0 = gamlss(formula,family = BP(mu.link = "identity"), trace = FALSE, data = data)
  }
  
  if(corstr == "independence"){
    cat("a gamlss object was returned")
    return(mod0)
  }
  
  beta = mod0$mu.coefficients # Chute inicial para beta
  phi = max(0.01,mod0$sigma.coefficients) # Chute inicial para phi
  
  # Modelo sob suposição de dependência
  cont = 1
  repeat{
    
    eta = X%*%beta
    if(linkmu == "log"){
      mu = as.vector(exp(eta)) # mi para a ligação logarítmica
    } else{
      mu = as.vector(eta) # mi para a ligação logarítmica
    }
    # Cálculo da função de variância
    vmu = mu^2 + mu
    
    # y estrela
    ys = log(y/(y+1))
    
    # mi estrela
    mus = digamma(mu*(phi+1))-digamma(mu*(phi+1)+phi+2)
    
    # Variância de b_ij
    vmus = trigamma(mu*(phi+1))-trigamma(mu*(phi+1)+phi+2)
    
    # Vetor b_i
    u = ys-mus
    
    #Matrizes utilizadas para o cálculo da equação de estimação
    if(linkmu == "log"){
      G = diag(as.vector(mu)) # G para a ligação logarítimica
    } else{
      G = diag(1,N,N) # G para a ligação logarítimica
    }
    A = diag(as.vector(vmus))
    Lambda = (phi+1)*G%*%A
    
    uc = split(u,id)
    scomb = matrix(u,n,nr[1],byrow = TRUE)
    if(corstr == "unstructured"){
      R = matrix(0,nr[1],nr[1])
      for(j in 1:nr[1]){
        for(l in j:nr[1]){
          num = sum(scomb[,j]*scomb[,l])
          den1 = sqrt(sum(scomb[,j]^2))
          den2 = sqrt(sum(scomb[,l]^2))
          R[j,l] = num/(den1*den2)
          R[l,j] = R[j,l]
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    } else if(corstr == "AR-1"){
      cnum = cden1 = cden2 = 0
      for(i in 1:n){
        for(j in 1:(nr[i]-1)){
          cnum = cnum + uc[[i]][j]*uc[[i]][j+1]
          cden1 = cden1 + (uc[[i]][j]^2)
          cden2 = cden2 + (uc[[i]][j+1]^2)
        }
      }
      alpha = cnum/sqrt(cden1*cden2)
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = list(NULL)
      for(i in 1:n){
        R[[i]] = matrix(0,nr[i],nr[i])
        for(j in 1:nr[i]){
          for(l in 1:nr[i]){
            R[[i]][j,l] = alpha^(abs(j-l))
          }
        }
      }
      # Matriz de correlação AR-1
      Rm = as.matrix(bdiag(R))
      R=R[[1]]
    } else if(corstr == "exchangeable"){
      cnum = cden = 0
      for(i in 1:n){
        aux = uc[[i]]%*%t(uc[[i]])
        cnum = cnum + sum(aux[upper.tri(aux)])*(2/(nr[i]-1))
        for(j in 1:(nr[i])){
          cden = cden + (uc[[i]][j]^2)
        }
      }
      alpha = (cnum/cden)
      Rm = matrix(0,N,N)
      R = list(NULL)
      for(i in 1:n){
        R[[i]] = matrix(alpha,nr[i],nr[i])
        diag(R[[i]]) = 1
      }
      # Matriz de correlação Uniforme
      Rm = as.matrix(bdiag(R))
      R = R[[1]]
    } else if(corstr == "one-dependent"){
      alpha = 0
      den = 0
      for(i in 1:n){
        for(j in 1:nr[1]){
          den = den + (scomb[i,j]^2)
        }
      }
      
      for(j in 1:(nr[1]-1)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+1])
        }
        alpha[j] = num/den
      }
      alpha = (N/n)*alpha
      
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      
      R = matrix(0,nr[1],nr[1])
      for(i in 1:nr[1]){
        for(j in 1:nr[1]){
          if(j == (i+1)){
            R[i,j] = alpha[i]
            R[j,i] = R[i,j]
          }
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    } else if(corstr == "two-dependent"){
      alpha1 = 0
      den = 0
      for(i in 1:n){
        for(j in 1:nr[1]){
          den = den + (scomb[i,j]^2)
        }
      }
      for(j in 1:(nr[1]-1)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+1])
        }
        alpha1[j] = num/den
      }
      alpha1 = (N/n)*alpha1
      alpha2 = 0
      for(j in 1:(nr[1]-2)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+2])
        }
        alpha2[j] = num/den
      }
      alpha2 = (N/n)*alpha2
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = matrix(0,nr[1],nr[1])
      for(i in 1:nr[1]){
        for(j in 1:nr[1]){
          if(j==(i+1)){
            R[i,j] = alpha1[i]
            R[j,i] = R[i,j]
          }
          if(j == (i+2)){
            R[i,j] = alpha2[i]
            R[j,i] = R[i,j]
          }
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    } else if(corstr == "one-dependent-stat"){
      alpha = 0
      den = 0
      for(i in 1:n){
        for(j in 1:nr[1]){
          den = den + (scomb[i,j]^2)
        }
      }
      for(j in 1:(nr[1]-1)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+1])
        }
        alpha[j] = num/den
      }
      alpha = (N/n)*alpha
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = matrix(0,nr[1],nr[1])
      for(i in 1:nr[1]){
        for(j in 1:nr[1]){
          if(j == (i+1)){
            R[i,j] = sum(alpha)/(nr[1]-1)
            R[j,i] = R[i,j]
          }
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    } else{
      alpha1 = 0
      den = 0
      for(i in 1:n){
        for(j in 1:nr[1]){
          den = den + (scomb[i,j]^2)
        }
      }
      for(j in 1:(nr[1]-1)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+1])
        }
        alpha1[j] = num/den
      }
      alpha1 = (N/n)*alpha1
      
      alpha2 = 0
      for(j in 1:(nr[1]-2)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+2])
        }
        alpha2[j] = num/den
      }
      alpha2 = (N/n)*alpha2
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      
      R = matrix(0,nr[1],nr[1])
      for(i in 1:nr[1]){
        for(j in 1:nr[1]){
          if(j == (i+1)){
            R[i,j] = sum(alpha1)/(nr[1]-1)
            R[j,i] = R[i,j]
          }
          if(j == (i+2)){
            R[i,j] = sum(alpha2)/(nr[1]-1)
            R[j,i] = R[i,j]
          }
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    }
    
    Omega = sqrt(A)%*%Rm%*%sqrt(A)
    W = Lambda%*%solve(Omega)%*%Lambda
    z = eta + solve(Lambda)%*%u
    
    #Novo valor de beta
    beta1 = solve(t(X)%*%W%*%X)%*%(t(X)%*%W%*%z)
    
    # Verificar se convergiu: beta1 é aproximadamente beta
    dif = abs(beta1-beta)
    if(sum(dif)<=(tol*p)){
      beta = beta1
      if(print == TRUE) cat("The algorithm converged")
      converg = 1
      break
    }
    
    # Se não convergir em 50 iterações o algoritmo para
    if(cont == maxiter){
      if(print == TRUE) cat("Maximum number of iterations reached")
      converg = 0
      break
    }
    beta = beta1
    
    # Resíduo de Pearson
    r = (y-mu)*(1/sqrt(vmu))
    
    # Cálculo do novo phi
    phi = 1/(sum(r^2)/(N-p))
    cont = cont + 1
  }
  
  # Matriz de sensibilidade
  S = -t(X)%*%W%*%X
  invOmega = solve(Omega)
  
  # Covariância de beta
  VarBeta = solve(S)%*%t(X)%*%Lambda%*%invOmega%*%u%*%t(u)%*%invOmega%*%Lambda%*%X%*%solve(S)
  
  # Estimativa do erro padrão de beta
  SEbeta = sqrt(diag(VarBeta))
  
  # A função retorna na primeira coluna as estimativas de beta e na segunda coluna o erro padrão
  # respectivo
  fit <- list()
  attr(fit, "class") <- c("geeBP")
  fit$title <- "geeBP:  BETA PRIME GENERALIZED ESTIMATING EQUATIONS"
  fit$model <- list()
  fit$model$link <- linkmu
  fit$model$varfun <- "mu(1+mu)"
  fit$model$corstr <- corstr
  fit$call <- call
  fit$formula <- formula
  fit$nclusters = n
  fit$clusters = nr
  fit$nobs <- N
  fit$iterations <- cont
  fit$coefficients <- beta
  eta <- as.vector(X %*% fit$coefficients)
  fit$linear.predictors <- eta
  mu <- as.vector(mu)
  fit$fitted.values <- mu
  fit$residuals <- r
  fit$family <- "Beta prime"
  fit$y <- as.vector(y)
  fit$id <- as.vector(id)
  fit$max.id <- max(nr)
  fit$working.correlation <- R
  fit$scale <- phi
  fit$robust.variance <- VarBeta
  fit$robust.se = SEbeta
  if(corstr == "unstructured"){
    fit$alpha = fit$working.correlation[upper.tri(fit$working.correlation)]
  }
  if(corstr == "AR-1"||corstr == "exchangeable"||corstr == "one-dependent"){
    fit$alpha = alpha
  }
  if(corstr == "two-dependent"){
    fit$alpha = list(diag1 = alpha1, diag2 = alpha2)
  }
  if(corstr == "one-dependent-stat"){
    fit$alpha = R[1,2]
  }
  if(corstr == "two-dependent-stat"){
    fit$alpha = c(R[1,2],R[1,3])
  }
  #fit$xnames <- colnames(model.matrix(formula, data = data))
  #dimnames(fit$robust.variance) <- list(fit$xnames, fit$xnames)
  #dimnames(fit$naive.variance) <- list(fit$xnames, fit$xnames)
  fit$comp$X = X
  fit$comp$W = W
  fit$comp$u = u
  fit$comp$Lambda = Lambda
  fit$comp$vmu = vmu
  fit$comp$A = A
  fit$comp$G = G
  fit$comp$Rm = Rm
  fit$comp$Omega = Omega
  # QIC
  Q = (1/phi)*(y*(log(mu)-log(y))+(y-1)*(log(1+mu)-log(1+y)))
  psiq = sum(Q)
  VI = A
  oi = t(X)%*%Lambda%*%solve(VI)%*%Lambda%*%X
  QIC = -2*psiq + 2*sum(diag(oi%*%VarBeta))
  CIC = sum(diag(oi%*%VarBeta))
  fit$QIC = QIC
  fit$CIC = CIC
  devi = -2*phi*Q
  D = sum(devi)
  phia = D/n
  fit$EQIC = (1/phia)*D + sum(log(2*pi*phia*(diag(A)+1/6))) + 2*CIC
  return(fit)
}




