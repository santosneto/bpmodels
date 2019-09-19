
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


#' Residuals
#'
#'
#'
#'@export

residuals.pearson <- function(model)
{

  mu <- fit$mu.fv
  phi <- fit$sigma.fv
  y <- model$y

  resP <- (sqrt(phi)*(y-mu))/sqrt(mu*(1+mu))


  resP
}



#' Envelope BP
#'
#'
#'
#'@export

envelope.BP <- function(model,k=100,link=c("log","identity"),color = "grey50", xlabel = "Theorical Quantile",ylabel = "Empirical Quantile",font="serif")
{

  n=model$N
  td  = model$residuals
  sigma = model$sigma.fv
  mu = model$mu.fv
  re <- matrix(0,n,k)
  mu.link <- model$mu.link
  sigma.link <- model$sigma.link

  if(length(model$mu.coefficients)==1) x <- model$mu.x else x <- model$mu.x[,-1];
  if(length(model$sigma.coefficients)==1) z <- model$sigma.x else z <- model$sigma.x[,-1];

  for(i in 1:k)
  {
    y1 <- mapply(rBP,1,mu,sigma)

    if(length(model$mu.coefficients)==1) form.x <- y1 ~ x-1 else  form.x <- y1 ~ x;
    if(length(model$sigma.coefficients)==1) form.z <- y1 ~ z-1 else  form.z <- y1 ~ z;

    conh0 = gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)

    if(mu.link=="log" & sigma.link == "identity"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="log",sigma.link="identity"),method=RS(),control = conh0)
    } else if(mu.link=="log" & sigma.link == "log"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="log",sigma.link="log"),method=RS(),control = conh0)
    } else if(mu.link=="log" & sigma.link == "sqrt"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="log",sigma.link="sqrt"),method=RS(),control = conh0)
    } else if(mu.link=="identity" & sigma.link == "sqrt"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="identity",sigma.link="sqrt"),method=RS(),control = conh0)
    } else if(mu.link=="identity" & sigma.link == "log"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="identity",sigma.link="log"),method=RS(),control = conh0)
    } else if(mu.link=="identity" & sigma.link == "identity"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="identity",sigma.link="identity"),method=RS(),control = conh0)
    } else if(mu.link=="sqrt" & sigma.link == "sqrt"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="sqrt",sigma.link="identity"),method=RS(),control = conh0)
    } else if(mu.link=="sqrt" & sigma.link == "log"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="sqrt",sigma.link="log"),method=RS(),control = conh0)
    } else{
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="sqrt",sigma.link="identity"),method=RS(),control = conh0)
    }

    rdf <- model1$residuals
    re[,i] <- sort(rdf)
  }
  e10 <- numeric(n)
  e20 <- numeric(n)
  e11 <- numeric(n)
  e21 <- numeric(n)
  e12 <- numeric(n)
  e22 <- numeric(n)

  for(l in 1:n)
  {
    eo = sort(re[l,])
    e10[l] = eo[ceiling(k*0.01)]
    e20[l] = eo[ceiling(k*(1 - 0.01))]
    e11[l] = eo[ceiling(k*0.05)]
    e21[l] = eo[ceiling(k*(1 - 0.05))]
    e12[l] = eo[ceiling(k*0.1)]
    e22[l] = eo[ceiling(k*(1 - 0.1))]
  }

  a <- qqnorm(e10, plot.it = FALSE)$x
  r <- qqnorm(td, plot.it = FALSE)$x
  xb = apply(re, 1, mean)
  rxb <- qqnorm(xb, plot.it = FALSE)$x

  df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
  ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel)+theme(text=element_text(size=30,family=font))
}


diag.BP <- function(model,mu.link = "log",sigma.link = "log",scheme="case.weight")
{

  x <- model$mu.x
  z  <- model$sigma.x
  y <- model$y
  p<-ncol(x)
  q<-ncol(z)

  linkstr <- mu.link
  linkobj <- make.link(linkstr)
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta  <- linkobj$mu.eta

  sigma_linkstr <- sigma.link
  sigma_linkobj <- make.link(sigma_linkstr)
  sigma_linkfun <- sigma_linkobj$linkfun
  sigma_linkinv <- sigma_linkobj$linkinv
  sigma_mu.eta  <- sigma_linkobj$mu.eta


  B=function(Delta,I,M)
  {
    B=(t(Delta)%*%(I-M)%*%Delta)
    return(B)
  }

  loglik <- function(vP)
  {
    betab = vP[1:p]
    alpha = vP[-(1:p)]
    eta   = as.vector(x%*%betab)
    tau   = as.vector(z%*%alpha)
    mu    = linkinv(eta)
    sigma = sigma_linkinv(tau)

    a <- mu*(1+sigma)
    b <- 2 + sigma

    fy <- (a-1)*log(y) - (a+b)*log(1+y) - lbeta(a,b)

    return(sum(fy))
  }

  muest <- model$mu.coefficients
  sigmaest <- model$sigma.coefficients
  x0<- c(muest,sigmaest)
  h0 <- hessian(loglik,x0)

  Ldelta= h0[(p+1):(p+q),(p+1):(p+q)]
  Lbeta=h0[1:p,1:p]
  b11=cbind(matrix(0, p, p), matrix(0, p, q))
  b12=cbind(matrix(0, q, p), solve(Ldelta))
  B1= rbind(b11, b12)  #parameter beta
  b211 =cbind(solve(Lbeta), matrix(0, p, q))
  b212= cbind(matrix(0, q, p), matrix(0, q, q))
  B2=rbind(b211,b212)  # parameter delta

  b311 =cbind(matrix(0, p, p), matrix(0, p, q))
  b312= cbind(matrix(0, q, p), matrix(0, q, q))
  B3=rbind(b311,b312)  # parameter theta

  if(scheme=="case.weight")
  {
    ############################Case Weight####################################

    mu <- model$mu.fv
    sigma <- model$sigma.fv
    eta <- linkfun(mu)
    ai <- mu.eta(eta)
    a <- mu*(1+sigma)
    b <- mu*(1+sigma)+sigma+2
    Phi <-  (1+sigma)
    yast <- log(y) - log(1+y)
    muast <- digamma(a) - digamma(b)
    dldm <- Phi*(yast - muast)
    Deltamu <- crossprod(x,diag(ai*dldm))


    tau <- sigma_linkfun(sigma)
    bi <- sigma_mu.eta(tau)
    ystar <- mu*log(y) - (1+mu)*log(1+y)
    mustar <- mu*digamma(a) - (1+mu)*digamma(b) + digamma(Phi+1)
    dldd <- ystar - mustar

    Deltasigma <- crossprod(z,diag(bi*dldd))

    Delta <- rbind(Deltamu,Deltasigma)

    ##################theta#########################
    BT<-B(Delta,solve(h0),B3)
    autovmaxthetaPC<- eigen(BT,symmetric=TRUE)$val[1]
    vetorpcthetaPC<- eigen(BT,symmetric=TRUE)$vec[,1]
    dmaxG.theta<-abs(vetorpcthetaPC)
    vCithetaPC<-2*abs(diag(BT))
    Cb0<-vCithetaPC
    Cb.theta<-Cb0/sum(Cb0)
    ######################betas########################
    BM<-B(Delta,solve(h0),B1)
    autovmaxbetaPC<-eigen(BM,symmetric=TRUE)$val[1]
    vetorpcbetaPC<-eigen(BM,symmetric=TRUE)$vec[,1]
    dmaxG.beta<-abs(vetorpcbetaPC)
    vCibetaPC<-2*abs(diag(BM))
    Cb1<-vCibetaPC
    Cb.beta<-Cb1/sum(Cb1)
    ####################alphas#########################
    BD<-B(Delta,solve(h0),B2)
    autovmaxdeltaPC<-eigen(BD,symmetric=TRUE)$val[1]
    vetordeltaPC<-eigen(BD,symmetric=TRUE)$vec[,1]
    dmaxG.alpha=abs(vetordeltaPC)
    vCideltaPC=2*abs(diag(BD))
    Cb2=vCideltaPC
    Cb.alpha=Cb2/sum(Cb2)

    result <- list(dmax.beta = dmaxG.beta,
                   dmax.alpha = dmaxG.alpha,
                   dmax.theta = dmaxG.theta,
                   Ci.beta = Cb.beta,
                   Ci.alpha = Cb.alpha,
                   Ci.theta = Cb.theta)
    return(result)
  }

  if(scheme=="response")
  {
    ############################Response####################################
    mu <- model$mu.fv
    sigma <- model$sigma.fv
    eta <- linkfun(mu)
    ai <- mu.eta(eta)
    tau <- sigma_linkfun(sigma)
    bi <- sigma_mu.eta(tau)
    sy<- sqrt((mu*(1+mu))/sigma)
    Phi <-  (1+sigma)

    dymu <- Phi*(1/(y*(1+y)))
    Deltamu <- crossprod(x,diag(ai*dymu*sy))
    p<-ncol(x)
    q<-ncol(z)
    dysigma <- mu*(1/(y*(1+y))) - 1/(1+y)
    Deltasigma <- crossprod(z,diag(bi*dysigma*sy))
    Delta <- rbind(Deltamu,Deltasigma)

    ###############thetas###########################
    BT<-B(Delta,solve(h0),B3)
    autovmaxthetaPC<- eigen(BT,symmetric=TRUE)$val[1]
    vetorthetaRP<- eigen(BT,symmetric=TRUE)$vec[,1]
    dmaxG.theta<-abs(vetorthetaRP)
    vCithetaRP<-2*abs(diag(BT))
    Cb0<-vCithetaRP
    Cb.theta<-Cb0/sum(Cb0)

    #################betas##########################
    BM=B(Delta,solve(h0),B1)
    autovmaxbetaRP <- eigen(BM,symmetric=TRUE)$val[1]
    vetorbetaRP <- eigen(BM,symmetric=TRUE)$vec[,1]
    dmaxG.beta <- abs(vetorbetaRP)
    vCibetaRP <- 2*abs(diag(BM))
    Cb1 <- vCibetaRP
    Cb.beta <- Cb1/sum(Cb1)
    ####################alpha#######################
    BD=B(Delta,solve(h0),B2)
    autovmaxdeltaRP <- eigen(BD,symmetric=TRUE)$val[1]
    vetordeltaRP <- eigen(BD,symmetric=TRUE)$vec[,1]
    dmaxG.alpha <- abs(vetordeltaRP)
    vCideltaRP <- 2*abs(diag(BD))
    Cb2 <- vCideltaRP
    Cb.alpha <- Cb2/sum(Cb2)


    result <- list(dmax.beta = dmaxG.beta,
                   dmax.alpha = dmaxG.alpha,
                   dmax.theta = dmaxG.theta,
                   Ci.beta = Cb.beta,
                   Ci.alpha = Cb.alpha,
                   Ci.theta = Cb.theta)
    return(result)
  }

  # if(scheme=="location")
  # {
  #   ############################Location Predictor####################################
  #   l <- lx
  #   mu <- model$mu.fv
  #   sigma <- model$sigma.fv
  #   eta <- linkfun(mu)
  #   ai <- mu.eta(eta)
  #   tau <- sigma_linkfun(sigma)
  #   bi <- sigma_mu.eta(tau)
  #   bl <- coef(model,what="mu")[l]
  #   sxl <- sd(x[,l])
  #   dmu <- (-1/(2*mu)) + sigma/((y*sigma) + y + (sigma*mu)) +  ((sigma+1)*y)/(4*(mu^2)) - (sigma^2)/(4*y*(sigma+1))
  #   ai1 <- rep(0,length(mu))
  #   dmu2 <- 1/(2*mu*mu) - (sigma^2)/(((y*sigma) + y + (sigma*mu))^2) - (y*(sigma+1))/(2*mu*mu*mu)
  #   ci <- dmu2*(ai^2) + dmu*ai1*ai
  #   xaux <- matrix(0,nrow=nrow(x),ncol=ncol(x))
  #   xaux[,l] <- ones(nrow(x),1)
  #   x0 <- xaux
  #
  #   Deltamu <- bl*sxl*crossprod(x,diag(ci)) + sxl*crossprod(x0,diag(dmu*ai))
  #
  #   dmusigma <- y/(((y*sigma) + y + (sigma*mu))^2) + y/(4*mu*mu) - (sigma*(sigma+2))/(4*(sigma+1)*(sigma+1)*y)
  #   mi <- dmusigma*bi*ai
  #   Deltasigma <- bl*sxl*crossprod(z,diag(mi))
  #
  #   Delta <- rbind(Deltamu,Deltasigma)
  #
  #
  #   ###############thetas###########################
  #   BT<-B(Delta,solve(h0),B3)
  #   autovmaxthetaPX<- eigen(BT,symmetric=TRUE)$val[1]
  #   vetorthetaPX<- eigen(BT,symmetric=TRUE)$vec[,1]
  #   dmaxG.theta<-abs(vetorthetaPX)
  #   vCithetaPX<-2*abs(diag(BT))
  #   Cb0<-vCithetaPX
  #   Cb.theta<-Cb0/sum(Cb0)
  #
  #   #################betas##########################
  #   BM=B(Delta,solve(h0),B1)
  #   autovmaxbetaPX <- eigen(BM,symmetric=TRUE)$val[1]
  #   vetorbetaPX <- eigen(BM,symmetric=TRUE)$vec[,1]
  #   dmaxG.beta <- abs(vetorbetaPX)
  #   vCibetaPX <- 2*abs(diag(BM))
  #   Cb1 <- vCibetaPX
  #   Cb.beta <- Cb1/sum(Cb1)
  #   ####################alpha#######################
  #   BD=B(Delta,solve(h0),B2)
  #   autovmaxdeltaPX <- eigen(BD,symmetric=TRUE)$val[1]
  #   vetordeltaPX <- eigen(BD,symmetric=TRUE)$vec[,1]
  #   dmaxG.alpha <- abs(vetordeltaPX)
  #   vCideltaPX <- 2*abs(diag(BD))
  #   Cb2 <- vCideltaPX
  #   Cb.alpha <- Cb2/sum(Cb2)
  #
  #   result <- list(dmax.beta = dmaxG.beta,
  #                  dmax.alpha = dmaxG.alpha,
  #                  dmax.theta = dmaxG.theta,
  #                  Ci.beta = Cb.beta,
  #                  Ci.alpha = Cb.alpha,
  #                  Ci.theta = Cb.theta)
  #   return(result)
  # }
  #
  # if(scheme=="precision")
  # {
  #   ############################Precision Predictor####################################
  #   mu <- model$mu.fv
  #   sigma <- model$sigma.fv
  #   eta <- linkfun(mu)
  #   ai <- mu.eta(eta)
  #   tau <- sigma_linkfun(sigma)
  #   bi <- sigma_mu.eta(tau)
  #   k <- lz
  #   ak <- coef(model,what="sigma")[k]
  #   szk <- sd(z[,k])
  #   dsigma <- (y+ mu)/((sigma*y) + y + (sigma*mu)) - y/(4*mu) - (sigma*(sigma+2)*mu)/(4*(sigma+1)*(sigma+1)*y) + sigma/(2*(sigma+1))
  #   Deltasigma <- crossprod(z,diag(bi*dsigma))
  #
  #   bi1 <- rep(2,length(sigma))
  #   dsigma2 <- 1/(2*(sigma+1)*(sigma+1)) - ((y+mu)^2)/(((y*sigma) + y + (sigma*mu))^2) - mu/(2*(sigma+1)*(sigma+1)*(sigma+1)*y)
  #   wi <- dsigma2*(bi^2) + dsigma*bi1*bi
  #   zaux <- matrix(0,nrow=nrow(z),ncol=ncol(z))
  #   zaux[,k] <- ones(nrow(z),1)
  #   z0 <- zaux
  #
  #   Deltasigma <- ak*szk*crossprod(z,diag(wi)) + szk*crossprod(z0,diag(dsigma*bi))
  #
  #   dmusigma <- y/(((y*sigma) + y + (sigma*mu))^2) + y/(4*mu*mu) - (sigma*(sigma+2))/(4*(sigma+1)*(sigma+1)*y)
  #   mi <- dmusigma*bi*ai
  #   Deltamu <- ak*szk*crossprod(x,diag(mi))
  #
  #   Delta <- rbind(Deltamu,Deltasigma)
  #
  #
  #   ###############thetas###########################
  #   BT<-B(Delta,solve(h0),B3)
  #   autovmaxthetaPZ<- eigen(BT,symmetric=TRUE)$val[1]
  #   vetorthetaPZ<- eigen(BT,symmetric=TRUE)$vec[,1]
  #   dmaxG.theta<-abs(vetorthetaPZ)
  #   vCithetaPZ<-2*abs(diag(BT))
  #   Cb0<-vCithetaPZ
  #   Cb.theta<-Cb0/sum(Cb0)
  #
  #   #################betas##########################
  #   BM=B(Delta,solve(h0),B1)
  #   autovmaxbetaPZ <- eigen(BM,symmetric=TRUE)$val[1]
  #   vetorbetaPZ <- eigen(BM,symmetric=TRUE)$vec[,1]
  #   dmaxG.beta <- abs(vetorbetaPZ)
  #   vCibetaPZ <- 2*abs(diag(BM))
  #   Cb1 <- vCibetaPZ
  #   Cb.beta <- Cb1/sum(Cb1)
  #   ####################alpha#######################
  #   BD=B(Delta,solve(h0),B2)
  #   autovmaxdeltaPZ <- eigen(BD,symmetric=TRUE)$val[1]
  #   vetordeltaPZ <- eigen(BD,symmetric=TRUE)$vec[,1]
  #   dmaxG.alpha <- abs(vetordeltaPZ)
  #   vCideltaPZ <- 2*abs(diag(BD))
  #   Cb2 <- vCideltaPZ
  #   Cb.alpha <- Cb2/sum(Cb2)
  #
  #   result <- list(dmax.beta = dmaxG.beta,
  #                  dmax.alpha = dmaxG.alpha,
  #                  dmax.theta = dmaxG.theta,
  #                  Ci.beta = Cb.beta,
  #                  Ci.alpha = Cb.alpha,
  #                  Ci.theta = Cb.theta)
  #   return(result)
  # }
  ############################Joint Predictor####################################

}

envelope.GA <- function(model,k=100,link=c("log","identity"),color = "grey50", xlabel = "Theorical Quantile",ylabel = "Empirical Quantile",font="serif")
{

  n=model$N
  td  = model$residuals
  sigma = model$sigma.fv
  mu = model$mu.fv
  re <- matrix(0,n,k)
  mu.link <- model$mu.link
  sigma.link <- model$sigma.link

  if(length(model$mu.coefficients)==1) x <- model$mu.x else x <- model$mu.x[,-1];
  if(length(model$sigma.coefficients)==1) z <- model$sigma.x else z <- model$sigma.x[,-1];

  for(i in 1:k)
  {
    y1 <- mapply(rGA,1,mu,sigma)

    if(length(model$mu.coefficients)==1) form.x <- y1 ~ x-1 else  form.x <- y1 ~ x;
    if(length(model$sigma.coefficients)==1) form.z <- y1 ~ z-1 else  form.z <- y1 ~ z;

    conh0 = gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)

    if(mu.link=="log" & sigma.link == "identity"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=GA(mu.link="log",sigma.link="identity"),method=RS(),control = conh0)
    } else if(mu.link=="log" & sigma.link == "log"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=GA(mu.link="log",sigma.link="log"),method=RS(),control = conh0)
    } else if(mu.link=="log" & sigma.link == "sqrt"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=GA(mu.link="log",sigma.link="sqrt"),method=RS(),control = conh0)
    } else if(mu.link=="identity" & sigma.link == "sqrt"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=GA(mu.link="identity",sigma.link="sqrt"),method=RS(),control = conh0)
    } else if(mu.link=="identity" & sigma.link == "log"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=GA(mu.link="identity",sigma.link="log"),method=RS(),control = conh0)
    } else if(mu.link=="identity" & sigma.link == "identity"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=GA(mu.link="identity",sigma.link="identity"),method=RS(),control = conh0)
    } else if(mu.link=="sqrt" & sigma.link == "sqrt"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=GA(mu.link="sqrt",sigma.link="identity"),method=RS(),control = conh0)
    } else if(mu.link=="sqrt" & sigma.link == "log"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=GA(mu.link="sqrt",sigma.link="log"),method=RS(),control = conh0)
    } else{
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=GA(mu.link="sqrt",sigma.link="identity"),method=RS(),control = conh0)
    }

    rdf <- model1$residuals
    re[,i] <- sort(rdf)
  }
  e10 <- numeric(n)
  e20 <- numeric(n)
  e11 <- numeric(n)
  e21 <- numeric(n)
  e12 <- numeric(n)
  e22 <- numeric(n)

  for(l in 1:n)
  {
    eo = sort(re[l,])
    e10[l] = eo[ceiling(k*0.01)]
    e20[l] = eo[ceiling(k*(1 - 0.01))]
    e11[l] = eo[ceiling(k*0.05)]
    e21[l] = eo[ceiling(k*(1 - 0.05))]
    e12[l] = eo[ceiling(k*0.1)]
    e22[l] = eo[ceiling(k*(1 - 0.1))]
  }

  a <- qqnorm(e10, plot.it = FALSE)$x
  r <- qqnorm(td, plot.it = FALSE)$x
  xb = apply(re, 1, mean)
  rxb <- qqnorm(xb, plot.it = FALSE)$x

  df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
  ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel)+theme(text=element_text(size=30,family=font))
}


envelope.IG <- function(model,k=100,link=c("log","identity"),color = "grey50", xlabel = "Theorical Quantile",ylabel = "Empirical Quantile",font="serif")
{

  n=model$N
  td  = model$residuals
  sigma = model$sigma.fv
  mu = model$mu.fv
  re <- matrix(0,n,k)
  mu.link <- model$mu.link
  sigma.link <- model$sigma.link

  if(length(model$mu.coefficients)==1) x <- model$mu.x else x <- model$mu.x[,-1];
  if(length(model$sigma.coefficients)==1) z <- model$sigma.x else z <- model$sigma.x[,-1];

  for(i in 1:k)
  {
    y1 <- mapply(rIG,1,mu,sigma)

    if(length(model$mu.coefficients)==1) form.x <- y1 ~ x-1 else  form.x <- y1 ~ x;
    if(length(model$sigma.coefficients)==1) form.z <- y1 ~ z-1 else  form.z <- y1 ~ z;

    conh0 = gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)

    if(mu.link=="log" & sigma.link == "identity"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=IG(mu.link="log",sigma.link="identity"),method=RS(),control = conh0)
    } else if(mu.link=="log" & sigma.link == "log"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=IG(mu.link="log",sigma.link="log"),method=RS(),control = conh0)
    } else if(mu.link=="log" & sigma.link == "sqrt"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=IG(mu.link="log",sigma.link="sqrt"),method=RS(),control = conh0)
    } else if(mu.link=="identity" & sigma.link == "sqrt"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=IG(mu.link="identity",sigma.link="sqrt"),method=RS(),control = conh0)
    } else if(mu.link=="identity" & sigma.link == "log"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=IG(mu.link="identity",sigma.link="log"),method=RS(),control = conh0)
    } else if(mu.link=="identity" & sigma.link == "identity"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=IG(mu.link="identity",sigma.link="identity"),method=RS(),control = conh0)
    } else if(mu.link=="sqrt" & sigma.link == "sqrt"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=IG(mu.link="sqrt",sigma.link="identity"),method=RS(),control = conh0)
    } else if(mu.link=="sqrt" & sigma.link == "log"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=IG(mu.link="sqrt",sigma.link="log"),method=RS(),control = conh0)
    } else{
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=IG(mu.link="sqrt",sigma.link="identity"),method=RS(),control = conh0)
    }

    rdf <- model1$residuals
    re[,i] <- sort(rdf)
  }
  e10 <- numeric(n)
  e20 <- numeric(n)
  e11 <- numeric(n)
  e21 <- numeric(n)
  e12 <- numeric(n)
  e22 <- numeric(n)

  for(l in 1:n)
  {
    eo = sort(re[l,])
    e10[l] = eo[ceiling(k*0.01)]
    e20[l] = eo[ceiling(k*(1 - 0.01))]
    e11[l] = eo[ceiling(k*0.05)]
    e21[l] = eo[ceiling(k*(1 - 0.05))]
    e12[l] = eo[ceiling(k*0.1)]
    e22[l] = eo[ceiling(k*(1 - 0.1))]
  }

  a <- qqnorm(e10, plot.it = FALSE)$x
  r <- qqnorm(td, plot.it = FALSE)$x
  xb = apply(re, 1, mean)
  rxb <- qqnorm(xb, plot.it = FALSE)$x

  df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
  ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel)+theme(text=element_text(size=30,family=font))
}


envelope <- function(model, k = 100, precision = c("fixed","varying"), dist = RBS(mu.link = "identity", sigma.link = "identity"),color = "grey50", xlabel = "Theorical quantile",ylabel = "Empirical quantile",font="serif")
{
  if (precision != "fixed") {
    alfa1 <- ceiling(k * alpha)
    alfa2 <- ceiling(k * (1 - alpha))
    n <- model$N
    td <- model$residuals
    sigma <- model$sigma.fv
    mu <- model$mu.fv
    re <- matrix(0, n, k)
    X <- model$mu.x
    Z <- model$sigma.x
    p <- ncol(X)
    q <- ncol(Z)
    for (i in 1:k) {
      y1 <- mapply(rRBS, n = 1, mu = mu, sigma = sigma)
      nresp <- y1
      if (p == 1) form <- nresp ~ 0 + X else form <- nresp ~ X[, -1]
      if (q == 1) form1 <- nresp ~ 0 + Z else form1 <- nresp ~ Z[, -1]
      conh0 = gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
      model1 <- gamlss(formula = form, sigma.formula = form1, family = dist, method = RS(), control = conh0)
      rdf <- model1$residuals
      re[, i] = sort(rdf)
    }
    e10 <- numeric(n)
    e20 <- numeric(n)
    e11 <- numeric(n)
    e21 <- numeric(n)
    e12 <- numeric(n)
    e22 <- numeric(n)
    for (l in 1:n) {
      eo = sort(re[l,])
      e10[l] = eo[ceiling(k*0.01)]
      e20[l] = eo[ceiling(k*(1 - 0.01))]
      e11[l] = eo[ceiling(k*0.05)]
      e21[l] = eo[ceiling(k*(1 - 0.05))]
      e12[l] = eo[ceiling(k*0.1)]
      e22[l] = eo[ceiling(k*(1 - 0.1))]
    }


    a <- qqnorm(e10, plot.it = FALSE)$x
    r <- qqnorm(td, plot.it = FALSE)$x
    xb = apply(re, 1, mean)
    rxb <- qqnorm(xb, plot.it = FALSE)$x

    df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
    ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel)+theme(text=element_text(size=30,family=font))
  }
  else {

    n <- model$N
    td <- model$residuals
    sigma <- model$sigma.fv[1]
    mu <- model$mu.fv
    alpha <- sqrt(2/sigma)
    bet.a <- (sigma * mu)/(sigma + 1)
    re <- matrix(0, n, k)
    X <- model$mu.x
    p <- ncol(X)
    mulink <- model$mu.link
    for (i in 1:k) {
      y1 <- mapply(rRBS, n = 1, mu = mu, sigma = sigma)
      nresp <- y1
      if (p == 1) form <- nresp ~ 0 + X else form <- nresp ~ X[, -1]

      conh0 <- gamlss.control(trace = FALSE, autostep = FALSE, save = TRUE)
      model1 <- gamlss(formula = form, family = dist, method = RS(), control = conh0)
      rdf <- model1$residuals
      re[, i] <- sort(rdf)
    }
    e10 <- numeric(n)
    e20 <- numeric(n)
    e11 <- numeric(n)
    e21 <- numeric(n)
    e12 <- numeric(n)
    e22 <- numeric(n)
    for (l in 1:n) {
      eo = sort(re[l,])
      e10[l] = eo[ceiling(k*0.01)]
      e20[l] = eo[ceiling(k*(1 - 0.01))]
      e11[l] = eo[ceiling(k*0.05)]
      e21[l] = eo[ceiling(k*(1 - 0.05))]
      e12[l] = eo[ceiling(k*0.1)]
      e22[l] = eo[ceiling(k*(1 - 0.1))]
    }
    a <- qqnorm(e10, plot.it = FALSE)$x
    r <- qqnorm(td, plot.it = FALSE)$x
    xb = apply(re, 1, mean)
    rxb <- qqnorm(xb, plot.it = FALSE)$x

    df <- data.frame(r=r,xab=a,emin=cbind(e10,e11,e12),emax=cbind(e20,e21,e22),xb=xb,td=td,rxb=rxb)
    ggplot(df,aes(r,td))+geom_ribbon(aes(x=xab, ymin=emin.e10, ymax=emax.e20),fill=color,alpha=0.5)  + geom_ribbon(aes(x=xab, ymin=emin.e11, ymax=emax.e21),fill=color,alpha=0.5) + geom_ribbon(aes(x=xab, ymin=emin.e12, ymax=emax.e22),fill=color,alpha=0.5) +scale_fill_gradient(low = "grey25", high = "grey75")+ geom_point() + geom_line(aes(rxb,xb),lty=2)+xlab(xlabel)+ylab(ylabel)+theme(text=element_text(size=30,family=font))
  }
}


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



geeBP = function(formula, data, id, tol = 0.001, maxiter = 80, corstr = "independence", linkmu = "log"){

  formula = as.formula(formula)
  X = as.matrix(model.matrix(formula, data = data)) # Matriz de especificação
  p = ncol(X) # Número de parâmetros
  y = model.frame(formula, data = data)[,1] # Variável resposta
  nt = table(id) # Número de repetições
  n = max(id) # Número de unidades experimentais
  N = nrow(X)
  mod0 = gamlss(formula,family = BP(mu.link = "log"), trace = FALSE, data = data,method = CG())
  beta = mod0$mu.coefficients # Chute inicial para beta
  phi = mod0$sigma.coefficients # Chute inicial para phi
  if(phi<0){
    phi = 1
  }

  # Modelo sob suposição de dependência
  cont = 1
  repeat{
    eta = X%*%beta
    if(linkmu == "log"){
      mu = as.vector(exp(eta)) # mi para a ligação logarítmica
    }
    if(linkmu == "identity"){
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
    }
    if(linkmu == "identity"){
      G = diag(N) # G para a ligação logarítimica
    }
    A = diag(as.vector(vmus))
    Lambda = (phi+1)*G%*%A

    uc = split(u,id)
    scomb = matrix(u,n,nt[1],byrow = TRUE)
    if(corstr == "unstructured"){
      R = matrix(0,nt[1],nt[1])
      for(j in 1:nt[1]){
        for(l in j:nt[1]){
          num = sum(scomb[,j]*scomb[,l])
          den1 = sqrt(sum(scomb[,j]^2))
          den2 = sqrt(sum(scomb[,l]^2))
          R[j,l] = num/(den1*den2)
          R[l,j] = R[j,l]
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    }
    if(corstr == "AR-1"){
      cnum = cden1 = cden2 = 0
      for(i in 1:n){
        for(j in 1:(nt[i]-1)){
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
        R[[i]] = matrix(0,nt[i],nt[i])
        for(j in 1:nt[i]){
          for(l in 1:nt[i]){
            R[[i]][j,l] = alpha^(abs(j-l))
          }
        }
      }
      # Matriz de correlação AR-1
      Rm = as.matrix(bdiag(R))
    }
    if(corstr == "exchangeable"){
      cnum = cden = 0
      for(i in 1:n){
        aux = uc[[i]]%*%t(uc[[i]])
        cnum = cnum + sum(aux[upper.tri(aux)])*(2/(nt[i]-1))
        for(j in 1:(nt[i])){
          cden = cden + (uc[[i]][j]^2)
        }
      }
      alpha = (cnum/cden)
      Rm = matrix(0,N,N)
      R = list(NULL)
      for(i in 1:n){
        R[[i]] = matrix(alpha,nt[i],nt[i])
        diag(R[[i]]) = 1
      }
      # Matriz de correlação Uniforme
      Rm = as.matrix(bdiag(R))
    }
    if(corstr == "one-dependent"){
      alpha = 0
      den = 0
      for(i in 1:n){
        for(j in 1:t[1]){
          den = den + (scomb[i,j]^2)
        }
      }
      for(j in 1:(nt[1]-1)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+1])
        }
        alpha[j] = num/den
      }
      alpha = (N/n)*alpha
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = matrix(0,nt[1],nt[1])
      for(i in 1:nt[1]){
        for(j in 1:nt[1]){
          if(j==(i+1)){
            R[i,j] = alpha[i]
            R[j,i] = R[i,j]
          }
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    }
    if(corstr == "two-dependent"){
      alpha1 = 0
      den = 0
      for(i in 1:n){
        for(j in 1:nt[1]){
          den = den + (scomb[i,j]^2)
        }
      }
      for(j in 1:(nt[1]-1)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+1])
        }
        alpha1[j] = num/den
      }
      alpha1 = (N/n)*alpha1
      alpha2 = 0
      for(j in 1:(nt[1]-2)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+2])
        }
        alpha2[j] = num/den
      }
      alpha2 = (N/n)*alpha2
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = matrix(0,nt[1],nt[1])
      for(i in 1:nt[1]){
        for(j in 1:nt[1]){
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
    }
    if(corstr == "one-dependent-stat"){
      alpha = 0
      den = 0
      for(i in 1:n){
        for(j in 1:nt[1]){
          den = den + (scomb[i,j]^2)
        }
      }
      for(j in 1:(nt[1]-1)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+1])
        }
        alpha[j] = num/den
      }
      alpha = (N/n)*alpha
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = matrix(0,nt[1],nt[1])
      for(i in 1:nt[1]){
        for(j in 1:nt[1]){
          if(j==(i+1)){
            R[i,j] = sum(alpha)/(nt[1]-1)
            R[j,i] = R[i,j]
          }
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    }
    if(corstr == "two-dependent-stat"){
      alpha1 = 0
      den = 0
      for(i in 1:n){
        for(j in 1:nt[1]){
          den = den + (scomb[i,j]^2)
        }
      }
      for(j in 1:(nt[1]-1)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+1])
        }
        alpha1[j] = num/den
      }
      alpha1 = (N/n)*alpha1
      alpha2 = 0
      for(j in 1:(nt[1]-2)){
        num = 0
        for(i in 1:n){
          num = num + (scomb[i,j]*scomb[i,j+2])
        }
        alpha2[j] = num/den
      }
      alpha2 = (N/n)*alpha2
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = matrix(0,nt[1],nt[1])
      for(i in 1:nt[1]){
        for(j in 1:nt[1]){
          if(j==(i+1)){
            R[i,j] = sum(alpha1)/(nt[1]-1)
            R[j,i] = R[i,j]
          }
          if(j == (i+2)){
            R[i,j] = sum(alpha2)/(nt[1]-1)
            R[j,i] = R[i,j]
          }
        }
      }
      diag(R) = 1
      Rm = kronecker(diag(n),R)
    }


    Omega = sqrt(A)%*%Rm%*%sqrt(A)
    W = Lambda%*%solve(Omega,tol=1e-150)%*%Lambda
    z = eta + solve(Lambda,tol=1e-150)%*%u

    #Novo valor de beta
    beta1 = solve(t(X)%*%W%*%X,tol=1e-150)%*%(t(X)%*%W%*%z)

    # Verificar se convergiu: beta1 é aproximadamente beta
    dif = abs(beta1-beta)
    if(sum(dif)<=(tol*p)){
      beta = beta1
      paste0("The algorithm converged")
      converg = 1
      break
    }

    # Se não convergir em 50 iterações o algoritmo para
    if(cont == maxiter){
      paste0("Maximum number of iterations reached")
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

  # Matriz de sensibilidade considerando a estrutura de independência
  if(corstr == "independence"){
    # Matriz de sensibilidade
    S = -t(X)%*%W%*%X
    invOmega = solve(Omega)

    # Covariância de beta
    VarBeta = solve(S)%*%t(X)%*%Lambda%*%invOmega%*%A%*%invOmega%*%Lambda%*%X%*%solve(S)
  }
  else
  {
    # Matriz de sensibilidade
    S = -t(X)%*%W%*%X
    invOmega = solve(Omega)

    # Covariância de beta
    VarBeta = solve(S)%*%t(X)%*%Lambda%*%invOmega%*%u%*%t(u)%*%invOmega%*%Lambda%*%X%*%solve(S)
  }

  # Estimativa do erro padrão de beta
  SEbeta = sqrt(diag(VarBeta))

  # A função retorna na primeira coluna as estimativas de beta e na segunda coluna o erro padrão
  # respectivo
  fit <- list()
  fit$X = X
  fit$W = W
  fit$Lam = Lambda
  fit$u = u
  fit$coef = beta
  fit$EP = SEbeta
  fit$rp = r
  fit$rcov = VarBeta
  fit$vmu = vmu
  fit$phi = phi
  fit$mu = mu
  fit$A = A
  fit$G = G
  fit$R = R
  fit$Rm = Rm
  fit$alpha = alpha
  fit$Omega = Omega
  return(fit)
}


#'Simulations
#'
#'
#'
#'@export
sim.gee.bp <- function(n=10,r=2,phi=2,REP=10,alpha=0.5)
{

  theta <- c(1.5,-0.5,0.1)
  nt <- r
  # Cluster
  idvec = rep(1:n, each = nt)
  # Tamanho da amostra total
  N = n*nt

  set.seed(13)
  cx = runif(N,0,1)
  set.seed(19)
  cw = rgamma(N,1,2)
  dias = rep(1:nt,n)
  # Valores de mi
  mi = exp(theta[1]+theta[2]*cx+theta[3]*cw)
  mi.fit <- matrix(mi,n,byrow=T)
  alfa <- alpha
  phin <- phi

  for(Matcor in c("AR-1","exchangeable"))
  {

    betaar1 = matrix(0,REP,3)
    sear1 = matrix(0,REP,3)
    betaex = matrix(0,REP,3)
    seex = matrix(0,REP,3)
    betaun = matrix(0,REP,3)
    seun = matrix(0,REP,3)
    betaind = matrix(0,REP,3)
    seind = matrix(0,REP,3)
    betagam = matrix(0,REP,3)
    segam = matrix(0,REP,3)
    betaod = matrix(0,REP,3)
    seod = matrix(0,REP,3)

    s = 0
    set.seed(1317291619)
    while(s < REP){
      s = s + 1
      # index

      #
      ym = matrix(0,n,nt)
      for(k in 1:n){
        mu <- rep(0,nt)
        Sigma = matrix(0,nt,nt)
        if(Matcor == "AR-1"){
          for(i in 1:nt){
            for(j in 1:nt){
              Sigma[i,j] = alfa^(abs(i-j))
            }
          }
        }
        if(Matcor == "exchangeable"){
          Sigma = matrix(alfa,nt,nt)
          diag(Sigma) = 1
        }
        rawvars <- mvrnorm(n=1, mu=mu, Sigma=Sigma)
        pvars <- pnorm(rawvars)
        ym[k,] <- qBP(pvars, mi.fit[k,], phin)
      }
      y = vec(t(ym))

      original = data.frame(y,cx,cw)


      resul1 = geeBP(y~cx+cw, data = original, id = idvec, corstr = "AR-1")
      resul2 = geeBP(y~cx+cw, data = original, id = idvec, corstr = "exchangeable")
      resul3 = geeBP(y~cx+cw, data = original, id = idvec, corstr = "unstructured")
      resul4 = gamlss(y~cx+cw,family = BP(mu.link = "log"), trace = FALSE, method = CG())
      if(Matcor == "AR-1"){

        resul5 =geeglm(y~cx+cw, id = idvec, data = original, corstr = "ar1", family = Gamma("log"))

      }
      if(Matcor == "exchangeable"){
        resul5 = geeglm(y~cx+cw, id = idvec, data = original, corstr = "ex", family = Gamma("log"))
      }
      resul6 = try(geeBP(y~cx+cw, data = original, id = idvec, corstr = "one-dependent-stat"))
      if(inherits(resul6, "try-error")){

        s = s - 1
        next
      }
      betaar1[s,] = resul1$coef
      sear1[s,] = resul1$EP
      betaex[s,] = resul2$coef
      seex[s,] = resul2$EP
      betaun[s,] = resul3$coef
      seun[s,] = resul3$EP
      betaind[s,] = resul4$mu.coefficients
      seind[s,] = sqrt(diag(vcov(resul4))[1:3])
      betagam[s,] = resul5$coefficients
      segam[s,] = sqrt(diag(vcov(resul5)))
      betaod[s,] = resul6$coef
      seod[s,] = resul6$EP
    }


    mtheta <- matrix(rep(theta,REP), REP, byrow = TRUE)
    bar1 <- colMeans((betaar1-mtheta))
    bex <- colMeans((betaex-mtheta))
    bun <- colMeans((betaun-mtheta))
    bind <- colMeans((betaind-mtheta))
    bgam <- colMeans((betagam-mtheta))
    bod <- colMeans((betaod-mtheta))

    # RMSE

    rsmear1 <- sqrt(colMeans((betaar1-mtheta)^2))
    rsmeex <- sqrt(colMeans((betaex-mtheta)^2))
    rsmeun <- sqrt(colMeans((betaun-mtheta)^2))
    rsmeind <- sqrt(colMeans((betaind-mtheta)^2))
    rsmegam <- sqrt(colMeans((betagam-mtheta)^2))
    rsmeod <- sqrt(colMeans((betaod-mtheta)^2))


    cat("Sampling unit size:",n,"\n")
    cat("Number of replicates:",nt,"\n")
    cat("Correlation:",alfa,"\n")
    cat("Precision:",phin,"\n")
    cat("Type of correlation matrices:",Matcor,"\n")
    cat("Monte Carlo replications:", REP,"\n")
    cat("------Correlation structures - AR1------\n")
    cat(" Mean: \n",sep=" | ", round(colMeans(betaar1),5),"\n EP: \n",round(colMeans(sear1),5) ,"\n SD: \n", round(colSds(betaar1),5),"\n Bias: \n",round(bar1,5),"\n RMSE: \n",round(rsmear1,5))
    cat("\n------Correlation structures - Exchangeable------\n")
    cat(" Mean: \n",sep=" | ", round(colMeans(betaex),5),"\n EP: \n",round(colMeans(seex),5) ,"\n SD: \n", round(colSds(betaex),5),"\n Bias: \n",round(bex,5),"\n RMSE: \n",round(rsmeex,5))
    cat("\n------Correlation structures - Unstructured------\n")
    cat(" Mean: \n",sep=" | ", round(colMeans(betaun),5),"\n EP: \n",round(colMeans(seun),5) ,"\n SD: \n", round(colSds(betaun),5),"\n Bias: \n",round(bun,5),"\n RMSE: \n",round(rsmeun,5))
    cat("\n------Correlation structures - Independence------\n")
    cat(" Mean: \n",sep=" | ", round(colMeans(betaind),5),"\n EP: \n",round(colMeans(seind),5) ,"\n SD: \n", round(colSds(betaind),5),"\n Bias: \n",round(bind,5),"\n RMSE: \n",round(rsmeind,5))
    cat("\n------Correlation structures - One-dependent-stat------\n")
    cat(" Mean: \n",sep=" | ", round(colMeans(betaod),5),"\n EP: \n",round(colMeans(seod),5) ,"\n SD: \n", round(colSds(betaod),5),"\n Bias: \n",round(bod,5),"\n RMSE: \n",round(rsmeod,5))
    cat("\n------GEE Gamma------\n")
    cat(" Mean: \n",sep=" | ", round(colMeans(betagam),5),"\n EP: \n",round(colMeans(segam),5) ,"\n SD: \n", round(colSds(betagam),5),"\n Bias: \n",round(bgam,5),"\n RMSE: \n",round(rsmegam,5),"\n")

  } #End
}
