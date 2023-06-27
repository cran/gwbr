#' @title Geographically Weighted  Beta Regression
#'
#' @description Fits a local regression model for each location using the beta distribution, recommended for rates and proportions, using a parametrization with mean (transformed by the link function) and precision parameter (called phi). For more details see Da Silva and Lima (2017).
#'
#' @param yvar A vector with the response variable name.
#' @param xvar A vector with descriptive variable(s) name(s).
#' @param lat A vector with the latitude variable name.
#' @param long A vector with the longitude variable name.
#' @param h The bandwidth parameter.
#' @param data A data set object with \code{yvar} and \code{xvar}.
#' @param xglobal A vector with descriptive variable(s) name(s) with global effect.
#' @param grid A data set with the location variables. Only used when the location variable are in another data set, different from data set used in parameter \code{data}. Variable name \code{"lat"} is expected for latitude and \code{"long"} for longitude.
#' @param method The kernel function used. The options are: \code{"fixed_g"}, \code{"fixed_bsq"} or \code{"adaptive_bsq"}. The default is \code{"fixed_g"}.
#' @param link The link function used in modeling. The options are: \code{"logit"}, \code{"probit"}, \code{"loglog"} or \code{"cloglog"}. The default is \code{"logit"}.
#' @param distancekm Logical. If \code{TRUE} use the distance in kilometers otherwise, use the Euclidean distance. The default is \code{TRUE}.
#' @param global Logical. If \code{TRUE} return to global model, giving the results from \code{betareg_gwbr} function. The default is \code{FALSE}.
#' @param maxint A maximum number of iterations to numerically maximize the log-likelihood function in search of the parameter estimates. The default is \code{maxint=100}.
#'
#' @return A list that contains:
#'
#' \itemize{
#' \item \code{parameter_estimates_qtls} - Parameter estimates quartiles and interquartile range.
#' \item \code{parameter_estimates_desc} - Parameter estimates mean, minimum and maximum.
#' \item \code{std_qtls} - Standard deviation quartiles and interquartile range.
#' \item \code{std_desc} - Standard deviation mean, minimum and maximum.
#' \item \code{est_n_parameters} - Number of parameters.
#' \item \code{est_gwr_parameters} - Effective number of parameters in the local model.
#' \item \code{phi} - Vector of precision parameter estimates.
#' \item \code{global_parameter} - Global parameter estimates, when existing.
#' \item \code{global_phi} - Global scale parameter estimate, when existing.
#' \item \code{global_parameter_tab} - Global parameter estimates table, when existing.
#' \item \code{residuals} - Table with observed values (\code{y}), estimated values (\code{yhat}), the link function applied in the estimated values (\code{eta}), pure residual (\code{res}), standardized residual (\code{resstd}), standardized weighted residual 2 (\code{resstd2}), residual deviance (\code{resdeviance}), Cooks distance (\code{cookD}), generalized leverage (\code{glbp}) and number of iterations (\code{iteration}).
#' \item \code{log_likelihood} - Log-likelihood of the fitted model.
#' \item \code{aicc} - Corrected Akaike information criterion.
#' \item \code{r2} - Pseudo R2 and adjusted pseudo R2 statistics.
#' \item \code{bp_test} - Breusch-Pagan test for heteroscedasticity.
#' \item \code{w} - Matrix of weights.
#' \item \code{parameters} - Table with parameter estimates of each model.
#' \item \code{significance} - Significance level of each model.
#' \item \code{bandwidth} - Bandwidth used.
#' \item \code{link_function} - The link function used in modeling.
#' }
#'
#' @examples
#' \donttest{
#' data(saopaulo)
#' output_list=gwbr("prop_landline",c("prop_urb", "prop_poor"),"y","x",116.3647,saopaulo)
#'
#' ## Descriptive statistics of the parameter estimates
#' output_list$parameter_estimates_desc
#'
#' ## Table with all parameter estimates and your respective statistics
#' output_list$parameters
#' }
#' @export

gwbr <- function(yvar, xvar, lat, long, h, data, xglobal=NA_character_, grid=data.frame(), method=c("fixed_g", "fixed_bsq", "adaptative_bsq"), link=c("logit", "probit", "loglog", "cloglog"), distancekm=TRUE, global=FALSE, maxint=100){
if(global==TRUE){
  result=betareg_gwbr(yvar=yvar, xvar=xvar, data=data, link=link, maxint=maxint)
}else{
    if(!is.numeric(h) | h<=0){
      message("ERROR: 'h' must be numeric and greater than 0")
      stop()
    }

    y=as.matrix(data[,yvar])
    x=as.matrix(data[,xvar])

    if(!is.na(xglobal)){
      xf=data[,xglobal]
      nf=ncol(xf)
    }

    n=length(y)
    x=cbind(matrix(1,n,1),x)
    for(i in 1:n){
      if(y[i]>=.99999 | y[i]<=.00001){
        y[i]=(y[i]*(n-1)+.5)/n
      }
    }
    yhat=matrix(0,n,1)
    k=ncol(x)

    if(length(link)==4){
      link=c("logit")
    }

    if(length(link)>1){
      message('ERROR: Link Function should be one of logit, loglog, cloglog or probit.')
      stop()
    }

    if(toupper(link)=="LOGIT"){
      yc=log(y/(1-y))
      linkf=function(eta){
        eta=ifelse(eta>700,700,eta)
        ilink=exp(eta)/(1+exp(eta))
        ilink=ifelse(ilink>.9999,.9999,ilink)
        ilink=ifelse(ilink<=1E-10,.00001,ilink)
        dlink=(1/ilink)+(1/(1-ilink))
        dlink2=1/((1-ilink)*(1-ilink))-(1/(ilink*ilink))
        links=as.matrix(cbind(ilink, dlink, dlink2))
        return(links)
      }
    }

    if(toupper(link)=="PROBIT"){
      yc=qnorm(y)
      linkf=function(eta){
        ilink=pnorm(eta)
        dlink=dnorm(ilink)
        dlink2=-(1/sqrt(2*acos(-1)))*ilink*exp(-ilink*ilink/2)
        links=as.matrix(cbind(ilink, dlink, dlink2))
        return(links)
      }
    }

    if(toupper(link)=="LOGLOG"){
      yc=-log(-log(y))
      linkf=function(eta){
        ilink=exp(-exp(-eta))
        ilink=ifelse(ilink<=0,.01,ilink)
        dlink=(1/log(ilink))*(-1/ilink)
        dlink2=(1/(log(ilink)*log(ilink)))*(1/(ilink*ilink)) + (1/log(ilink))*(1/(ilink*ilink))
        links=as.matrix(cbind(ilink, dlink, dlink2))
        return(links)
      }
    }

    if(toupper(link)=="CLOGLOG"){
      yc=log(-log(1-y))
      linkf=function(eta){
        ilink=1-exp(-exp(eta))
        ilink=ifelse(ilink>=.99999,.99,ilink)
        dlink=(-1/log(1-ilink))*(1/(1-ilink))
        dlink2=(-1/(log(1-ilink)*log(1-ilink))*(1/(1-ilink)))*(1/(1-ilink)) + (-1/log(1-ilink))*(1/((1-ilink)*(1-ilink)))
        links=as.matrix(cbind(ilink, dlink, dlink2))
        return(links)
      }
    }

    if(sum(toupper(link)==c("LOGIT", "PROBIT", "PROBIT", "LOGLOG", "CLOGLOG"))==0){
      message('ERROR: Link Function should be one of logit, loglog, cloglog or probit.')
    }

    if(length(method)==3){
      method=c("fixed_g")
    }

    if(sum(toupper(method)==c("FIXED_G", "FIXED_BSQ", "ADAPTIVE_BSQ"))==0){
      message('ERROR: Method should be one of fixed_g, fixed_bsq or adaptive_bsq.')
      stop()
    }

    if(length(method)>1){
      message('ERROR: Method should be one of fixed_g, fixed_bsq or adaptive_bsq.')
      stop()
    }

    coord=cbind(data[,long],data[,lat])

    if(nrow(grid)==0){
      points=cbind(data[,long],data[,lat])
    }else{
      points=cbind(grid[,long],grid[,lat])
    }

    m=nrow(points)
    dist_ <- as.matrix(dist(coord))
    seq=as.matrix(1:n)

    bi=matrix(0,ncol(x)*m,4)
    rsqri=matrix(0,m,1)
    sumwi=matrix(0,m,1)
    varbi=matrix(0,ncol(x)*m,1)
    varbigg=matrix(0,ncol(x)*m,1)
    varbis=matrix(0,ncol(x)*m,ncol(x))
    s=matrix(0,m,1)
    s_=matrix(0,m,m)
    s2=matrix(0,m,1)
    bit=matrix(0,m,ncol(x)+1)
    ss=matrix(0,m,1)

    rnd=1

    for(i in 1:m){
      seqi=matrix(i,n,1)
      dist=cbind(seqi,seq,as.matrix(dist_[,i]))
      if(distancekm==TRUE){
        dist[,3]=dist[,3]*111
      }
      u=nrow(dist)
      w=matrix(0,u,1)

      for(jj in 1:u){
        if(toupper(method)=="FIXED_G"){
          w[jj]=exp(-0.5*(dist[jj,3]/h)^2)
        }
        if(toupper(method)=="FIXED_BSQ"){
          w[jj]=(1-(dist[jj,3]/h)^2)^2
        }
      }
      if(toupper(method)=="ADAPTIVE_BSQ"){
        dist=dist[order(dist[,3]),]
        dist=cbind(dist,1:nrow(dist))
        w=matrix(0,n,2)
        hn=dist[round(h),3]
        for(jj in 1:(n-1)){
          if(dist[jj,4]<=h){
            w[jj,1]=(1-(dist[jj,3]/hn)^2)^2
          }else{
            w[jj,1]=0
          }
          w[jj,2]=dist[jj,2]
        }
        w=w[,1]
      }

      w=as.vector(w)
      if(det(t(x)%*%as.matrix(w*x))==0){
        b=matrix(0,ncol(x),1)
      }else{
        b=solve(t(x)%*%as.matrix(w*x))%*%t(x)%*%(w*yc)
      }
      m1=(i-1)*ncol(x)+1
      m2=m1+(ncol(x)-1)
      bi[m1:m2,1]=i
      bi[m1:m2,2]=b
      bi[m1:m2,3]=points[i,1]
      bi[m1:m2,4]=points[i,2]
      yhat[i]=x[i,]%*%b
      ss[i]=(x[i,]%*%solve(t(x)%*%as.matrix(w*x))%*%t(x))[1]
    }
    e=yc-yhat

    betai_=matrix(t(bi[,1:2]),m, byrow = T)
    i=seq(2,ncol(betai_),2)
    betai_=betai_[,i]

    eta=as.matrix(x)%*%t(betai_)

    mu=linkf(eta)[,1:m]
    gmu=linkf(eta)[,(m+1):(2*m)]

    sigma2=as.numeric(t(e)%*%e)/((n-sum(ss))*(gmu*gmu))
    phi=matrix(0,1*m,1)
    for(i in 1:m){
      for(j in 1:m){
        phi[i]=phi[i]+mu[j,i]*(1-mu[j,i])/(sigma2[j,i]*n)
      }
    }
    phi=ifelse(phi<1,phi,phi-1)
    parameters=cbind(betai_,phi)

    max_like=function(param){
      betai2=t(param[1:length(param)-1])
      phii2=param[length(param)]
      etai2=as.matrix(x)%*%t(betai2)
      mu2=as.matrix(linkf(etai2)[,1])
      lgamma1=lgamma(phii2*mu2)
      arg=(1-mu2)*phii2
      arg=ifelse(arg<=0,1E-23,arg)
      lgamma2=lgamma(arg)
      lgamma3=as.vector(phii2*mu2-1)*log(y)
      lgamma4=((1-mu2)*phii2-1)*log(1-y)
      lnl=0
      for(j in 1:length(y)){
        lnl=lnl+(lgamma(phii2)-lgamma1[j]-lgamma2[j]+lgamma3[j]+lgamma4[j])*w[j]
      }
      return(lnl)
    }

    beta=matrix(0,m,ncol(x)+3)
    yhat=matrix(0,m,1)
    iteration=matrix(0,m,1)
    phi=matrix(0,m,1)
    stdb=matrix(0,m,ncol(x))
    stdphi=matrix(0,m,1)
    w1=matrix(0,m,m)
    bb=matrix(0,ncol(x)*n,n)

    for(i in 1:m){
      seqi=matrix(i,n,1)
      dist=cbind(seqi,seq,as.matrix(dist_[,i]))
      if(distancekm==TRUE){
        dist[,3]=dist[,3]*111
      }
      u=nrow(dist)
      w=matrix(0,u,1)
      for(jj in 1:u){
        if(toupper(method)=="FIXED_G"){
          w[jj]=exp(-0.5*(dist[jj,3]/h)^2)
        }
        if(toupper(method)=="FIXED_BSQ"){
          w[jj]=ifelse(dist[jj,3]<h,
                       (1-(dist[jj,3]/h)^2)^2,
                       0)
        }
      }
      if(toupper(method)=="ADAPTIVE_BSQ"){
        dist=dist[order(dist[,3]),]
        dist=cbind(dist,1:nrow(dist))
        w=matrix(0,n,2)
        hn=dist[round(h),3]
        for(jj in 1:n){
          if(dist[jj,4]<=h){
            w[jj,1]=(1-(dist[jj,3]/hn)^2)^2
          }else{
            w[jj,1]=0
          }
          w[jj,2]=dist[jj,2]
        }
        w=w[,1]
      }

      parami=t(parameters[i,])
      optn=matrix(1)
      con=rbind(cbind(matrix(NA,1,k),.01),matrix(NA,1,k+1))
      it=0

      dif=1
      parami[length(parami)]=ifelse(parami[length(parami)]<=0,.01,parami[length(parami)])
      betai=t(parami[1:length(parami)-1])
      phii=parami[length(parami)]
      etaini=as.matrix(x)%*%t(betai)
      while(abs(dif)>0.00000001 & it<maxint){
        mu=as.matrix(linkf(etaini)[,1])
        mu=ifelse(mu<1e-7,1e-7,mu)
        mu=ifelse(mu>=0.9999999,0.99999,mu)
        if(sum(mu>0.993)==length(mu)){
          it=maxint
        }else{
          gmu=as.matrix(linkf(etaini)[,2])
          ye=log(y/(1-y))
          mue=digamma(mu*phii)-digamma((1-mu)*phii)
          t=1/gmu
          c=phii*(trigamma(mu*phii)*mu - trigamma((1-mu)*phii)*(1-mu))
          z=(phii*trigamma(mu*phii)+phii*trigamma((1-mu)*phii))/(gmu*gmu)
          d=w*(trigamma(mu*phii)*mu*mu + trigamma((1-mu)*phii)*(1-mu)*(1-mu)-trigamma(phii))

          w=as.vector(w)

          ub= phii*(t((x*w))%*%(t*(ye-mue)))
          up=0
          for(ii in 1:length(mu)){
            up=up+(mu[ii]*(ye[ii]-mue[ii])+log(1-y[ii])-digamma((1-mu[ii])%*%phii)+digamma(phii))*w[ii]
          }
          z=as.vector(z)

          kbb=phii*t(x)%*%as.matrix(x*w*z)
          kbp=t(x*w)%*%(t*c)
          kpp=sum(d)
          km=rbind(cbind(kbb,kbp),cbind(t(kbp),kpp))

          if(det(km)>0){
            paramf=t(parami)+solve(km, tol=0)%*%rbind(ub,up)
          }else{
            paramf=parami
          }
          b=t(paramf)
          b[ncol(b)]=ifelse(b[ncol(b)]<=0,.01,b[ncol(b)])
          mlike1=max_like(b)
          mlike2=max_like(parami)
          dif=mlike1-mlike2

          parami=b
          betai=t(parami[1:length(parami)-1])
          phii=parami[length(parami)]
          etaini=as.matrix(x)%*%t(betai)
          it=it+1
        }
      }
      xr=parami
      beta[i,ncol(x)+1]=i
      beta[i,ncol(x)+2]=coord[i,1]
      beta[i,ncol(x)+3]=coord[i,2]
      beta[i,1:ncol(x)]=t(xr[1:ncol(x)])
      yhat[i]=x[i,]%*%as.matrix(xr[1:ncol(x)])
      iteration[i]=it
      phi[i]=xr[ncol(xr)]
      s_[i,]=(x[i,]%*%solve(t(x)%*%as.matrix(x*w*z))%*%t(x*w*z))

      w1[,i]=w
    }
    vv2=sum(diag(t(s_)%*%s_))
    vv1=sum(diag(s_))

    v1=sum(s)
    v1=2*vv1-vv2

    mu=linkf(yhat)[,1]
    gmu=linkf(yhat)[,2]

    t=1/gmu
    c=phi*(trigamma(mu*phi)*mu - trigamma((1-mu)*phi)*(1-mu))
    z=(phi*trigamma(mu*phi)+phi*trigamma((1-mu)*phi))/(gmu*gmu)

    for(i in 1:m){
      d=w1[,i]*(trigamma(mu*phi)*mu*mu + trigamma((1-mu)*phi)*(1-mu)*(1-mu)-trigamma(phi))
      kbb=(t(as.vector(phi)*x)%*%as.matrix(x*as.vector(z)*w1[,i]))
      kbp=t(x*w1[,i])%*%(t*c)
      kpp=sum(d)
      km=rbind(cbind(kbb,kbp),cbind(t(kbp),kpp))
      if(det(km)>0){
        kinv=solve(km, tol=0)
      }else{
        kinv=matrix(1E10,nrow(km),ncol(km))
      }
      dp=sqrt(diag(kinv))
      stdb[i,]=t(dp[1:ncol(x)])
      stdphi[i]=dp[length(dp)]
    }

    if(!is.na(xglobal)){
      if(toupper(link)=="LOGIT"){yc <- log(y/(1-y))}else{
        if(toupper(link)=="CLOGLOG"){yc <- log(-log(1-y))}else{
          if(toupper(link)=="LOGLOG"){yc <- -log(-log(y))}else{
            if(toupper(link)=="PROBIT"){yc <- qnorm(y)}
          }
        }
      }
      yc <- as.matrix(yc)

      is=diag(n)-s_
      bf=solve(t(xf)%*%t(is)%*%is%*%xf)%*%t(xf)%*%t(is)%*%is%*%yc
      ef=is%*%yc-(is%*%xf)%*%bf

      is2=t(is)%*%is
      etai=is%*%xf%*%bf
      mu=linkf(etai)[,1]
      gmu=linkf(etai)[,2]
      sigma2=as.numeric(t(ef)%*%ef)/((n-k)*(gmu*gmu))

      global_phi=0
      for(i in 1:n){
        global_phi=global_phi+mu[i]%*%(1-mu[i])/(sigma2[i]%*%n)
      }
      global_phi=ifelse(global_phi<1,global_phi,global_phi-1)

      param=rbind(bf,global_phi)

      max_likegf <- function(param){
        it=it+1
        beta=param[1:ncol(param)-1]
        phi=param[ncol(param)]
        eta=is%*%xf%*%beta
        mu=linkf(eta)[,1]
        lgamma1=t(lgamma(phi%*%mu))
        arg=(1-mu)*phi
        arg=ifelse(arg<=0,1E-23,arg)
        lgamma2=lgamma(arg)
        lgamma3=(phi*mu-1)*log(y)
        lgamma4=((1-mu)*phi-1)*log(1-y)
        lnl=0
        n=nrow(y)
        for(i in 1:n){
          lnl= lnl+lgamma(phi)-lgamma1[i]-lgamma2[i]+lgamma3[i]+lgamma4[i]
        }
        return(lnl)
      }
      param=rbind(bf,global_phi)
      it=0

      parami=t(rbind(bf,global_phi))
      dif=1
      etai=is%*%xf%*%bf

      while(abs(dif)>0.00000001 & it<maxint){
        mu=linkf(etai)[,1]
        gmu= linkf(etai)[,2]
        ye=log(y/(1-y))
        mue=digamma(mu%*%global_phi)-digamma((1-mu)%*%global_phi)
        t=1/gmu
        c=t(global_phi%*%(as.numeric(trigamma(mu%*%global_phi))*mu - as.numeric(trigamma((1-mu)%*%global_phi))*(1-mu)))
        z=(as.numeric(global_phi)*trigamma(mu%*%global_phi)+as.numeric(global_phi)*trigamma((1-mu)%*%global_phi))/(gmu*gmu)
        d=as.numeric(trigamma(mu%*%global_phi))*mu*mu + as.numeric(trigamma((1-mu)%*%global_phi))*(1-mu)*(1-mu)-as.numeric(trigamma(global_phi))


        ub=global_phi%*%(t(xf)%*%t(is)%*%(t*(ye-as.numeric(mue))))
        up=0
        for(i in 1:n){
          up=up+as.numeric(as.numeric(mu[i]%*%(ye[i]-mue[i]))+log(1-y[i])-digamma((1-mu[i])*global_phi)+digamma(global_phi))
        }
        kbb=global_phi%*%t(xf)%*%t(is)%*%(z*(is%*%xf))
        kbp=t(xf)%*%t(is)%*%(t*c)
        kpp=sum(d)
        km=rbind(cbind(kbb,kbp),cbind(t(kbp),kpp))

        param=t(parami)+solve(km, tol=0)%*%rbind(ub,up)
        b=t(param)
        b[ncol(b)]=ifelse(b[ncol(b)]<=0,.01,b[ncol(b)])
        mlike1=max_likegf(b)
        mlike2=max_likegf(parami)
        dif=mlike1-mlike2

        parami=b
        betai=parami[1:ncol(parami)-1]
        etai=is%*%xf%*%betai
        phi=parami[ncol(parami)]
        xr=parami
        it=it+1
      }
      bf=xr[1]
      phif=xr[2]

      global_tab=data.frame(xglobal,bf,phif)

      rnd_=2
    }

    res=y-mu

    beta_=cbind(beta[,1:ncol(x)],phi)
    qntl=sapply(as.data.frame(beta_), function(x) quantile(x, probs = c(.25,.5,.75), type = 2))
    qntl=rbind(qntl,(qntl[3,]-qntl[1,]))
    descriptb=rbind(sapply(as.data.frame(beta_), mean),
                    sapply(as.data.frame(beta_), min),
                    sapply(as.data.frame(beta_), max))
    rownames(qntl)=c("P25", "P50", "P75", "IQR")
    rownames(descriptb)=c("Mean", "Min", "Max")
    colnames(qntl)=colnames(descriptb)=c("Intercept", xvar, "Phi")

    stdbeta_=cbind(stdb,stdphi)
    qntls=sapply(as.data.frame(stdbeta_), function(x) quantile(x, probs = c(.25,.5,.75), type = 2))
    qntls=rbind(qntls,(qntls[3,]-qntls[1,]))
    descripts=rbind(sapply(as.data.frame(stdbeta_), mean),
                    sapply(as.data.frame(stdbeta_), min),
                    sapply(as.data.frame(stdbeta_), max))
    rownames(qntls)=c("P25", "P50", "P75", "IQR")
    rownames(descripts)=c("Mean", "Min", "Max")
    colnames(qntls)=colnames(descripts)=c("Intercept", xvar, "Phi")

    tstat=beta[,1:ncol(x)]/stdb
    probt=2*(1-pt(abs(tstat),m-k))
    tstatp=phi/stdphi
    probtp=2*(1-pt(abs(tstatp),m-k))
    bistdt_=cbind(beta[,(ncol(x)+1):ncol(beta)],
                  beta[,1:ncol(x)],
                  stdb,tstat,probt,phi,stdphi,tstatp,probtp)
    colname1_=c("Intercept", xvar)
    label_=rbind(matrix(c("std_"),ncol(x)),
                 matrix(c("tstat_"),ncol(x)),
                 matrix(c("probt_"),ncol(x)))
    colname_=c(c("id", "x", "y"),
               colname1_,
               t(paste0(label_, matrix(t(colname1_),ncol(x)))),
               c("phi", "std_phi", "tstat_phi", "probt_phi"))


    lgamma1t=lgamma(phi*y)
    arg=(1-y)*phi
    arg=ifelse(arg<=0,1E-23,arg)
    lgamma2t=lgamma((1-y)*phi)
    lgamma3t=(phi*y-1)*log(y)
    lgamma4t=((1-y)*phi-1)*log(1-y)

    lgamma1=lgamma(phi*mu)
    arg=(1-mu)*phi
    arg=ifelse(arg<=0,1E-23,arg)
    lgamma2=lgamma((1-mu)*phi)
    lgamma3=(phi*mu-1)*log(y)
    lgamma4=((1-mu)*phi-1)*log(1-y)

    lnlmutil=lgamma(phi)-lgamma1t-lgamma2t+lgamma3t+lgamma4t
    lnl= lgamma(phi)-lgamma1-lgamma2+lgamma3+lgamma4
    resdeviance=sign(y-mu)*sqrt(2*abs(lnlmutil-lnl))

    vary=mu*(1-mu)/(1+phi)
    resstd=(y-mu)/sqrt(vary)

    ye=log(y/(1-y))
    mue=digamma(mu*phi)-digamma((1-mu)*phi)
    m=1/(y*(1-y))
    gmu2=linkf(beta[,1:ncol(x)])[,3]
    q=(phi*(trigamma(mu*phi)+trigamma((1-mu)*phi))+(ye-mue)*gmu2/gmu)/(gmu*gmu)
    f=c-(ye-mue)
    g=sum(d)-t((1/phi)*f*t)%*%as.matrix(x)%*%solve(t(x)%*%as.matrix(as.vector(q)*x))%*%t(x)%*%(t*f)
    b=-(y-mu)/(y*(1-y))
    glb=as.matrix(t*x)%*%solve(t(x)%*%as.matrix(as.vector(q)*x))%*%t(x*t*as.vector(m))
    glbp=diag(glb+((1/as.vector(as.numeric(g)*phi))*(t*as.matrix(x)%*%solve(t(x)%*%as.matrix(as.vector(q)*x))%*%t(x)%*%(t*f)%*%((t(f)*t(t))%*%as.matrix(x)%*%solve(t(x)%*%as.matrix(as.vector(q)*x))%*%t(x*t*as.vector(m))-t(b)))))

    z=as.vector(z)
    H=diag(as.matrix(sqrt(z)*x)%*%solve(t(x)%*%as.matrix(z*x))%*%t(x*sqrt(z)))
    cookD=H*(resstd*resstd)/(k*(1-H)*(1-H))

    resstd2=(ye-mue)/sqrt((trigamma(mu*phi)+trigamma((1-mu)*phi))*(1-H))

    eta=yhat
    mat=cbind(eta,yc)
    pseudor2=(cor(mat)%*%t(cor(mat)) -1)[1,1]
    adjpr2=1-((n-1)/(n-v1))*(1-pseudor2)

    sseb=t(res)%*%res
    gbp=matrix(0,n,1)
    fbp=sseb/n
    for(i in 1:n){
      tmp=res[i]**2
      gbp[i]=tmp/fbp - 1
    }
    lm=.5*t(gbp)%*%as.matrix(x)%*%solve(t(x)%*%as.matrix(x))%*%t(x)%*%gbp
    probbp=(1-pchisq(abs(lm),k-1))
    vecbp=cbind(lm,k-1,probbp)

    AIC= 2*(v1) - 2*sum(lnl)
    AICC=AIC+2*(v1)*(v1+1)/(n-v1-1)

    yhat=mu
    res_=data.frame(y,yhat,eta,res,resstd,resstd2,resdeviance,cookD,glbp,iteration)
    parameters2_=data.frame(bistdt_)
    names(parameters2_)=colname_

    sig_=matrix("not significant at 90%",n,ncol(x))
    for(i in 1:n){
      for(j in 1:ncol(x)){
        if(probt[i,j]<0.01*((ncol(x))/v1)){
          sig_[i,j]="significant at 99%"
        }else{
          if(probt[i,j]<0.05*((ncol(x))/v1)){
            sig_[i,j]="significant at 95%"
          }else{
            if(probt[i,j]<0.1*((ncol(x))/v1)){
              sig_[i,j]="significant at 90%"
            }else{
              sig_[i,j]="not significant at 90%"
            }
          }
        }
      }
    }

    sigp_=matrix("not significant at 90%",n,1)
    for(i in 1:n){
      if(probtp[i]<0.01*((ncol(x))/v1)){
        sigp_[i]="significant at 99%"
      }else{
        if(probtp[i]<0.05*((ncol(x))/v1)){
          sigp_[i]="significant at 95%"
        }else{
          if(probtp[i]<0.1*((ncol(x))/v1)){
            sigp_[i]="significant at 90%"
          }else{
            sigp_[i]="not significant at 90%"
          }
        }
      }
    }

    sig_=as.data.frame(cbind(sig_,sigp_))
    names(sig_)=c(paste0("sig_",c(colname1_,"Phi")))

    if(!is.na(xglobal)){
      global_par=rbind(xglobal,bf)
    }else{
      global_par=global_phi=global_tab=c("There's no global parameters in the model")

    }

    result <- list(parameter_estimates_qtls=as.data.frame(qntl),
                   parameter_estimates_desc=as.data.frame(descriptb),
                   std_qtls=as.data.frame(qntls),
                   std_desc=as.data.frame(descripts),
                   est_n_parameters=c(vv1=vv1,vv2=vv2),
                   est_gwr_parameters=as.numeric(v1),
                   phi=as.vector(phi),
                   global_parameter=global_par,
                   global_phi=global_phi,
                   global_parameter_tab=global_tab,
                   residuals=res_,
                   log_likelihood=sum(lnl),
                   aic=AIC,
                   aicc=AICC,
                   r2=c(`Pseudo R2`=pseudor2, `Adj. Pseudo R2`=adjpr2),
                   bp_test=c(`Statistic Value`=vecbp[1],`df`=vecbp[2],`p-value`=vecbp[3]),
                   w=as.matrix(w1),
                   parameters=parameters2_,
                   significance=sig_,
                   bandwidth=as.numeric(h),
                   link_function=as.character(link))

}
return(result)
}

