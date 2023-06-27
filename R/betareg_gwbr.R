#' @title Global Beta Regression Model
#'
#' @description Fits a global regression model using the beta distribution, recommended for rates and proportions, via maximum likelihood using a parametrization with mean (transformed by the link function) and precision parameter (called phi). For more details see Ferrari and Cribari-Neto (2004).
#'
#' @param yvar A vector with the response variable name.
#' @param xvar A vector with descriptive variable(s) name(s).
#' @param data A data set object with \code{yvar} and \code{xvar}.
#' @param link The link function used in modeling. The options are: \code{"logit"}, \code{"probit"}, \code{"loglog"} or \code{"cloglog"}. The default is \code{"logit"}.
#' @param maxint A Maximum number of iterations to numerically maximize the log-likelihood function in search of the estimators. The default is \code{maxint=100}.
#'
#' @return A list that contains:
#'
#' \itemize{
#' \item \code{parameter_estimates} - Parameter estimates.
#' \item \code{phi} - Precision parameter estimate.
#' \item \code{residuals} - Table with observed values (\code{y}), estimated values in classical regression (\code{yhatcl}), pure residual in classical regression (\code{ecl}), estimated values (\code{yhat}), the link function applied in the estimated values (\code{eta}), pure residual (\code{res}), standardized residual (\code{resstd}), standardized weighted residual 2 (\code{resstd2}), residual deviance (\code{resdeviance}), Cooks distance (\code{cookD}) and generalized leverage (\code{glbp}).
#' \item \code{log_likelihood} - Log-likelihood of the fitted model.
#' \item \code{aicc} - Corrected Akaike information criterion.
#' \item \code{r2} - Pseudo R2 and adjusted pseudo R2 statistics.
#' \item \code{bp_test} - Breusch-Pagan test for heteroscedasticity.
#' \item \code{link_function} - The link function used in modeling.
#' \item \code{n_iter} - Number of iterations used in convergence.
#' }
#'
#' @examples
#' data(saopaulo)
#' output_list=betareg_gwbr("prop_landline",c("prop_urb","prop_poor"),saopaulo)
#'
#' ## Parameters
#' output_list$parameter_estimates
#'
#' ## R2 and AICc
#' output_list$r2
#' output_list$aicc
#' @export

betareg_gwbr=function(yvar,xvar,data,link=c("logit", "probit", "loglog", "cloglog"), maxint=100){
  x=data[,xvar]
  y=data[,yvar]
  n <- length(y)
  for(i in 1:n){
    if(y[i]>=0.99999 | y[i]<=0.00001){
      y[i] <- (y[i]*(n-1)+.5)/n
    }
  }
  x <- as.matrix(cbind(x=matrix(1, n, 1), x))
  k <- as.numeric(ncol(x))

  if(length(link)==4){
    link=c("logit")
  }

  if(length(link)>1){
    message('ERROR: Link Function should be one of logit, loglog, cloglog or probit.')
    stop()
  }

  if(toupper(link)=="LOGIT"){yc <- log(y/(1-y))}else{
    if(toupper(link)=="CLOGLOG"){yc <- log(-log(1-y))}else{
      if(toupper(link)=="LOGLOG"){yc <- -log(-log(y))}else{
        if(toupper(link)=="PROBIT"){yc <- qnorm(y)}else{
          message('ERROR: Link Function should be one of logit, loglog, cloglog or probit.')
          stop()
        }
      }
    }
  }

  betai <- solve(t(x)%*%x)%*%t(x)%*%yc
  e <- yc-x%*%betai

  if(toupper(link)=="LOGIT"){
    linkf <- function(eta){
      ilink <- exp(eta)/(1+exp(eta))
      ilink <- ifelse(ilink==1,.99,ilink)
      dlink <- (1/ilink)+(1/(1-ilink))
      dlink2 <- 1/((1-ilink)*(1-ilink))-(1/(ilink*ilink))
      links <- as.matrix(cbind(ilink, dlink, dlink2))
      return(links)
    }
  }

  if(toupper(link)=="PROBIT"){
    linkf <- function(eta){
      ilink <- pnorm(eta)
      dlink <- dnorm(ilink)
      dlink2 <- -(1/sqrt(2*acos(-1)))*ilink*exp(-ilink*ilink/2)
      links <- as.matrix(cbind(ilink, dlink, dlink2))
      return(links)
    }
  }

  if(toupper(link)=="LOGLOG"){
    linkf <- function(eta){
      ilink <- exp(-exp(-eta))
      ilink <- ifelse(ilink<=0,.01,ilink)
      dlink <- (1/log(ilink))*(-1/ilink)
      dlink2 <- (1/(log(ilink)*log(ilink)))*(1/(ilink*ilink)) + (1/log(ilink))*(1/(ilink*ilink))
      links <- as.matrix(cbind(ilink, dlink, dlink2))
      return(links)
    }
  }

  if(toupper(link)=="CLOGLOG"){
    linkf <- function(eta){
      ilink <- 1-exp(-exp(eta))
      ilink <- ifelse(ilink>=.99999,.99,ilink)
      dlink <- (-1/log(1-ilink))*(1/(1-ilink))
      dlink2 <- (-1/(log(1-ilink)*log(1-ilink))*(1/(1-ilink)))*(1/(1-ilink)) + (-1/log(1-ilink))*(1/((1-ilink)*(1-ilink)))
      links <- as.matrix(cbind(ilink, dlink, dlink2))
      return(links)
    }
  }

  etai <- x%*%betai
  mu <- linkf(etai)[,1]
  gmu <- linkf(etai)[,2]
  sigma2 <- as.numeric(t(e)%*%e)/(t((n-k)%*%gmu)*gmu)

  phi <- 0
  for(i in 1:n){
    phi <- phi+mu[i]%*%(1-mu[i])/(sigma2[i]%*%n)
  }
  phi <- ifelse(phi<1, phi, phi-1)

  param <- t(rbind(betai,phi))

  it <- 0

  max_like <- function(param){
    it <- it+1
    beta <- param[1:ncol(param)-1]
    phi <- param[ncol(param)]
    eta <- x%*%beta
    mu <- linkf(eta)[,1]
    lgamma1 <- lgamma(phi%*%mu)
    arg <- (1-mu)*phi
    arg <- ifelse(arg<=0,1E-23,arg)
    lgamma2 <- lgamma(arg)
    lgamma3 <- (phi%*%mu-1)*log(y)
    lgamma4 <- ((1-mu)*phi-1)*log(1-y)
    lnl <- 0
    n <- length(y)
    for(i in 1:n){
      lnl <- lnl+lgamma(phi)-lgamma1[i]-lgamma2[i]+lgamma3[i]+lgamma4[i]
    }
    return(lnl)
  }
  param <- rbind(betai,phi)
  optn <- as.matrix(1)
  con <- rbind(cbind(matrix(NA, 1, k),.01), matrix(NA, 1, k+1))
  it <- 0
  parami <- t(rbind(betai,phi))
  dif <- 1
  etai <- x%*%betai

  phi <- as.numeric(phi)
  repeat{
    mu <- linkf(etai)[,1]
    gmu <- linkf(etai)[,2]
    ye <- log(y/(1-y))
    mue <- digamma(mu*phi)-digamma((1-mu)*phi)
    t <- 1/gmu
    c <- (trigamma(mu*phi)*mu-trigamma((1-mu)*phi)*(1-mu))*phi
    z <- (trigamma(mu*phi)*phi+trigamma((1-mu)*phi)*phi)/(gmu*gmu)
    d <- trigamma(mu*phi)*mu*mu + trigamma((1-mu)*phi)*(1-mu)*(1-mu)-as.numeric(trigamma(phi))

    ub <- (t(x)%*%(t*(ye-mue)))*phi
    up <- 0
    for(i in 1:n){
      up <- up+mu[i]%*%(ye[i]-mue[i])+log(1-y[i])-digamma((1-mu[i])*phi)+digamma(phi)
    }

    kbb <- as.numeric(phi)*t(x)%*%(as.numeric(z)*x)
    kbp <- t(x)%*%(t*c)
    kpp <- sum(d)
    km <- rbind(cbind(kbb,kbp),cbind(t(kbp),kpp))

    param <- t(parami)+solve(km)%*%rbind(ub,up)
    b <- t(param)
    b[ncol(b)] <- ifelse(b[ncol(b)]<=0,.01,b[ncol(b)])
    mlike1 <- max_like(b)
    mlike2 <- max_like(parami)
    dif <- mlike1-mlike2

    parami <- b
    betai <- parami[1:ncol(parami)-1]
    etai <- x%*%betai
    phi <- parami[ncol(parami)]
    xr <- parami
    it <- it+1
    if(abs(dif)<=0.00000001 | it>=maxint){
      break
    }
  }

  kinv <- solve(km)
  dp <- sqrt(diag(kinv))

  param <- xr
  beta <- param[1:ncol(x)]
  phi <- param[ncol(x)+1]
  eta <- x%*%beta
  mu <- linkf(eta)[,1]
  gmu <- linkf(eta)[,2]

  tstat <- c(beta,phi)/dp
  probt <- 2*(1-pt(abs(tstat),n-k))
  vec <- cbind(c(beta,phi),dp)

  res <- y-mu

  lgamma1t <- lgamma(phi*y)
  arg <- (1-y)*phi
  arg <- ifelse(arg<=0,1e-23,arg)
  lgamma2t <- lgamma((1-y)*phi)
  lgamma3t <- (phi*y-1)*log(y)
  lgamma4t <- ((1-y)*phi-1)*log(1-y)

  lgamma1 <- lgamma(phi*mu)
  arg <- (1-mu)*phi;
  arg <- ifelse(arg<=0,1E-23,arg)
  lgamma2 <- lgamma((1-mu)*phi)
  lgamma3 <- (phi*mu-1)*log(y)
  lgamma4 <- ((1-mu)*phi-1)*log(1-y)

  lnlmutil <- lgamma(phi)-lgamma1t-lgamma2t+lgamma3t+lgamma4t
  lnl <- lgamma(phi)-lgamma1-lgamma2+lgamma3+lgamma4
  resdeviance <- sign(y-mu)*sqrt(2*abs(lnlmutil-lnl))

  vary <- mu*(1-mu)/(1+phi)
  resstd <- (y-mu)/sqrt(vary)

  ye <- log(y/(1-y))
  mue <- digamma(mu*phi)-digamma((1-mu)*phi)
  m <- 1/(y*(1-y))
  gmu2 <- linkf(eta)[,3]
  q <- (phi*(trigamma(mu*phi)+trigamma((1-mu)*phi))+(ye-mue)*gmu2/gmu)/(gmu*gmu)
  f <- c-(ye-mue)
  g <- sum(d)-(1/phi)*(t(f)*t(t))%*%x%*%solve(t(x)%*%(q*x))%*%t(x)%*%(t*f)
  b <- -(y-mu)/(y*(1-y))
  glb <- (t*x)%*%solve(t(x)%*%(q*x))%*%t(x*t*m)
  glbp <- diag(glb+as.numeric(1/(g*phi))*(t*x%*%solve(t(x)%*%(q*x))%*%t(x)%*%(t*f))%*%((t(f)*t(t))%*%x%*%solve(t(x)%*%(q*x))%*%t(x*t*m)-t(b)))

  h <- diag(sqrt(z)*x%*%solve(t(x)%*%(z*x))%*%t(x*sqrt(z)))
  cookD <- h*(resstd*resstd)/(k%*%(1-h)*(1-h))

  resstd2 <- (ye-mue)/sqrt((trigamma(mu*phi)+trigamma((1-mu)*phi))*(1-h))

  eta <- x%*%beta
  mat <- cbind(eta,yc)
  pseudor2 <- (cor(mat)%*%cor(mat)-1)[1,1]
  adjpr2 <- 1-((n-1)/(n-ncol(x)-1))*(1-pseudor2)

  aic <- 2*(k+1)-2*sum(lnl)
  aicc <- aic+2*(k+1)*(k+2)/(n-(k+1)-1)

  sseb <- t(res)%*%res
  gbp <- matrix(0, n, 1)
  fbp <- sseb/n
  for(i in 1:n){
    tmp <- res[i]*res[i]
    gbp[i] <- tmp/fbp-1
  }
  lm <- .5*t(gbp)%*%x%*%solve(t(x)%*%x)%*%t(x)%*%gbp
  probbp <- (1-pchisq(abs(lm),k-1))
  vecbp <- c(lm, k-1, probbp)

  odds <- exp(beta[2:length(beta)])
  oddsl <- c(0,odds)

  betacl <- solve(t(x)%*%x)%*%t(x)%*%y
  yhatcl <- x%*%betacl
  ecl <- y-yhatcl
  sse <- t(ecl)%*%ecl
  ssr <- (t(betacl)%*%t(x)%*%y-(1/n)%*%t(y)%*%matrix(1,n,n)%*%y)
  mse <- sse/(n-k)
  ssto <- t(y)%*%y-(1/n)%*%t(y)%*%matrix(1,n,n)%*%y
  r2 <- ssr/ssto
  adjr2 <- 1-((n-1)/(n-k))*(1-r2)
  vbeta <- as.numeric(mse)*solve(t(x)%*%x)
  dpcl <- sqrt(diag(vbeta))
  veccl <- cbind(betacl,dpcl)
  tstatcl <- betacl/dpcl
  probtcl <- 2*(1-pt(abs(tstatcl),n-k))

  gbp <- matrix(0,n,1)
  fbp <- sse/n
  for(i in 1:n){
    tmp <- ecl[i]^2
    gbp[i] <- tmp/fbp-1
  }
  lm <- .5*t(gbp)%*%x%*%solve(t(x)%*%x)%*%t(x)%*%gbp

  sigmahatcl <- t(ecl)%*%ecl/n
  lnlcl <- -n*log(sigmahatcl)/2-n*log(2*acos(-1))/2-sum((y-yhatcl)*(y-yhatcl))/(2*sigmahatcl)

  aicl <- 2*k-2*lnlcl
  aiccl <- aicl+2*k*(k+1)/(n-k-1)

  res=data.frame(y,
                 yhatcl=as.vector(yhatcl),
                 ecl,
                 yhat=mu,
                 eta=as.vector(eta),
                 res,
                 resstd,
                 resstd2,
                 resdeviance,
                 cookD=as.vector(cookD),
                 glbp=as.vector(ecl))

  parameter_estimates=data.frame(cbind(vec,tstat,probt,c(oddsl,NA)), row.names = c("Intercept",colnames(x)[-1],"Phi"))
  names(parameter_estimates)=c("Estimate", "Std. Error", "t Value","Pr>|t|", "Odds Ratio")

  result <- list(parameter_estimates=parameter_estimates,
                 phi=phi,
                 residuals=res,
                 log_likelihood=sum(lnl),
                 aicc=aicc,
                 r2=data.frame('Pseudo R2'=pseudor2, 'Adj pseudo R2'=adjpr2),
                 bp_test=data.frame('Statistic Value'=vecbp[1],'df'=round(vecbp[2]),'p-value'=vecbp[3]),
                 link_function=link,
                 n_iter=it-2)

  return(result)
}
