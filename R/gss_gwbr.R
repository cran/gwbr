#' @title Golden Section Search Algorithm
#'
#' @description The Golden Section Search (GSS) algorithm is used in searching for the best bandwidth for geographically weighted regression. For more details see Da Silva and Mendes (2018).
#'
#' @param yvar A vector with the response variable name.
#' @param xvar A vector with descriptive variable(s) name(s).
#' @param lat A vector with the latitude variable name.
#' @param long A vector with the longitude variable name.
#' @param data A data set object with \code{yvar} and \code{xvar}.
#' @param method Kernel function used to set bandwidth parameter. The options are: \code{"fixed_g"}, \code{"fixed_bsq"} or \code{"adaptive_bsq"}. The default is \code{"fixed_g"}.
#' @param link The link function used in modeling. The options are: \code{"logit"}, \code{"probit"}, \code{"loglog"} or \code{"cloglog"}. The default is \code{"logit"}.
#' @param type Can be \code{"cv"}, when the Cross-Validation function is used to estimate the bandwidth or \code{"aic"}, when the AIC function is used. The default is \code{"cv"}.
#' @param globalmin Logical. If \code{TRUE} search for the global minimum. The default is \code{TRUE}.
#' @param distancekm Logical. If \code{TRUE} use the distance in kilometers otherwise, use the Euclidean distance. The default is \code{TRUE}.
#' @param maxint A maximum number of iterations to numerically maximize the log-likelihood function in search of parameter estimates. The default is \code{maxint=100}.
#'
#' @return A list that contains:
#'
#' \itemize{
#' \item \code{global_min} - Global minimum of the function, giving the best bandwidth (\code{h}).
#' \item \code{local_mins} - Local minimums of the function.
#' \item \code{type} - Function used to estimate the bandwidth.
#' }
#'
#' @examples
#' \donttest{
#' data(saopaulo)
#' output_list=gss_gwbr("prop_landline",c("prop_urb","prop_poor"),"y","x",saopaulo,"fixed_g")
#'
#' ## Best bandwidth
#' output_list$global_min
#' }
#' @export

gss_gwbr <- function(yvar, xvar, lat, long, data, method=c("fixed_g", "fixed_bsq", "adaptive_bsq"), link=c("logit", "probit", "loglog", "cloglog"), type=c("cv", "aic"), globalmin=TRUE, distancekm=TRUE, maxint=100){

  y=as.matrix(data[,yvar])
  x=as.matrix(data[,xvar])

  n=length(y)
  x=cbind(matrix(1,n,1),x)
  k=ncol(x)
  yhat=matrix(0,n,1)
  ss=matrix(0,n,1)

  PHI="local"

  for(i in 1:n){
    if(y[i]>=0.99999 | y[i]<=0.00001){
      y[i]=(y[i]*(n-1)+.5)/n
    }
  }

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

  if(sum(toupper(link)==c("LOGIT", "PROBIT", "LOGLOG", "CLOGLOG"))==0){
    message('ERROR: Link Function should be one of LOGIT, LOGLOG, CLOGLOG or PROBIT.')
    stop()
  }

  if(length(method)==3){
    method=c("fixed_g")
  }

  if(sum(toupper(method)==c("FIXED_G", "FIXED_BSQ", "ADAPTIVE_BSQ"))==0){
    message('ERROR: Method should be one of FIXED_G, FIXED_BSQ or ADAPTIVE_BSQ.')
    stop()
  }

  if(length(method)>1){
    message('ERROR: Method should be one of fixed_g, fixed_bsq or adaptive_bsq.')
    stop()
  }

  if(length(type)==2){
    type=c("cv")
  }

  if(sum(toupper(type)==c("CV", "AIC"))==0){
    message('ERROR: Type should be one of CV or AIC.')
    stop()
  }

  max_like=function(param){
    betai2=param[1:length(param)-1]
    phii2=param[length(param)]
    etai2=as.matrix(x)%*%betai2
    mu2=linkf(etai2)[,1]
    lgamma1=lgamma(phii2%*%mu2)
    arg=(1-mu2)*phii2
    arg=ifelse(arg<=0,1E-23,arg)
    lgamma2=lgamma(arg)
    lgamma3=(phii2*mu2-1)*log(y)
    lgamma4=((1-mu2)*phii2-1)*log(1-y)
    lnl=0
    n=length(y)
    for(ll in 1:n){
      lnl=lnl+(lgamma(phii2)-lgamma1[ll]-lgamma2[ll]+lgamma3[ll]+lgamma4[ll])%*%w[ll]
    }
    return(lnl)
  }

  coord=cbind(data[,long],data[,lat])

  dist_=as.matrix(dist(coord))
  seq=c(1:n)

  if(toupper(method)=="FIXED_G"){
    if(toupper(type)=="CV"){
      ax=0
      bx=trunc(max(dist_)+1)
      if(distancekm==TRUE){
        bx=bx*111
      }

      r=0.61803399
      tol=0.001
      h0=ax
      h3=bx
      h1=bx-r*(bx-ax)
      h2=ax+r*(bx-ax)

      bi=matrix(0,ncol(x)*n,2)
      for(i in 1:n){
        seqi=matrix(i,n,1)
        dist=cbind(seqi,seq,dist_[,i])
        if(distancekm==TRUE){
          dist[,3]=dist[,3]*111
        }
        u=nrow(dist)
        w=matrix(0,u,1)

        for(jj in 1:u){
          w[jj]=exp(-0.5*(dist[jj,3]/h1)**2)
        }
        w[i]=0
        w <- as.vector(w)
        if(det(t(x)%*%as.matrix(w*x))==0){
          b=matrix(0,ncol(x),1)
        }else{
          b=solve(t(x)%*%as.matrix(w*x))%*%t(x)%*%as.matrix(w*yc)
        }
        m1=(i-1)*ncol(x)+1
        m2=m1+(ncol(x)-1)
        bi[m1:m2,1]=i
        bi[m1:m2,2]=b
        yhat[i]=x[i,]%*%b
        ss[i]=(x[i,]%*%solve(t(x)%*%as.matrix(w*x))%*%t(x))[1]
      }

      e=yc-yhat

      betai_=matrix(t(bi[,1:2]),n, byrow = T)
      ii=seq(2,ncol(betai_),2)
      betai_=betai_[,ii]

      eta=as.matrix(x)%*%t(betai_)

      mu=linkf(eta)[,1:n]
      gmu=linkf(eta)[,(n+1):(2*n)]
      sigma2=as.numeric(t(e)%*%e)/((n-sum(ss))*(gmu*gmu))
      phi=matrix(0,1*n,1)
      for(ii in 1:n){
        for(j in 1:n){
          phi[ii]=phi[ii]+mu[j,ii]*(1-mu[j,ii])/(sigma2[j,ii]*n)
        }
      }
      phi=ifelse(phi<1,phi,phi-1)
      parameters=cbind(betai_,(phi-1))
      yhat=matrix(0,n,1)

      cv=function(h){
        phi=matrix(0,1*n,1)
        s=matrix(0,n,1)
        s_=matrix(0,n,n)
        w1=matrix(0,n,n)
        w11=matrix(0,n,n)

        for(i in 1:n){
          seqi=matrix(i,n,1)
          dist=cbind(seqi,seq,dist_[,i])
          if(distancekm==TRUE){
            dist[,3]=dist[,3]*111
          }
          u=nrow(dist)
          w=matrix(0,u,1)

          for(jj in 1:u){
            w[jj]=exp(-0.5*(dist[jj,3]/h)**2)
            w[i]=0
          }

          parami=t(parameters[i,])
          it=0
          dif=1
          parami[ncol(parami)]=ifelse(parami[ncol(parami)]<=0,.01,parami[ncol(parami)])
          betai=t(parami[1:ncol(parami)-1])
          phii=parami[ncol(parami)]
          etaini=as.matrix(x)%*%t(betai)

          while(abs(dif)>0.00000001 & it<maxint){
            mu=linkf(etaini)[,1]
            mu=ifelse(mu<1e-7,1e-7,mu)
            mu=ifelse(mu>=0.9999999,0.99999,mu)
            if(sum(mu>0.993)==length(mu)){
              it=maxint
            }else{
              gmu=linkf(etaini)[,2]
              ye=log(y/(1-y))
              mue=digamma(mu*phii)-digamma((1-mu)*phii)
              t=1/gmu
              c=t(phii%*%(trigamma(mu*phii)*mu - trigamma((1-mu)*phii)*(1-mu)))
              z=t((phii%*%trigamma(mu*phii)+phii%*%trigamma((1-mu)*phii))/(gmu*gmu))
              d=w*(trigamma(mu*phii)*mu*mu + trigamma((1-mu)*phii)*(1-mu)*(1-mu)-trigamma(phii))
              w <- as.vector(w)
              ub= phii*(t((x*w))%*%(t*(ye-mue)))
              up=0
              for(ii in 1:length(mu)){
                up=up+(mu[ii]*(ye[ii]-mue[ii])+log(1-y[ii])-digamma((1-mu[ii])%*%phii)+digamma(phii))*w[ii]
              }
              z <- as.vector(z)
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
          b=t(xr[1:ncol(x)])
          s_[i,]=(x[i,]%*%solve(t(x)%*%as.matrix(x*w*z))%*%t(x*w*z))[1,]
          yhat[i]=x[i,]%*%t(b)
        }
        cv=t(y-yhat)%*%(y-yhat)

        v1=NULL


        res=c(cv,v1)
        return(res)
      }

      ax=0
      bx=trunc(max(dist_)+1)
      if(distancekm==TRUE){
        bx=bx*111
      }
      r=0.61803399
      tol=0.1

      if(globalmin==FALSE){
        lower=ax
        upper=bx
        xmin=matrix(0,1,2)
        end_iter=1
      }else{
        lower=cbind(ax,(1-r)*bx,r*bx)
        upper=cbind((1-r)*bx,r*bx,bx)
        xmin=matrix(0,3,2)
        end_iter=3
      }

      min_bandwidth_=NULL
      for(gmy in 1:end_iter){
        ax1=lower[gmy]
        bx1=upper[gmy]
        h0=ax1
        h3=bx1
        h1=bx1-r*(bx1-ax1)
        h2=ax1+r*(bx1-ax1)

        rnd_=1
        res1=cv(h1)
        cv1=res1[1]
        res2=cv(h2)
        cv2=res2[1]

        dd_output=c(gmy,h1,cv1,h2,cv2)

        intcv=1
        while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & intcv<200){
          if(cv2<cv1){
            h0=h1
            h1=h3-r*(h3-h0)
            h2=h0+r*(h3-h0)
            cv1=cv2
            res2=cv(h2)
            cv2=res2[1]
          }else{
            h3=h2
            h1=h3-r*(h3-h0)
            h2=h0+r*(h3-h0)
            cv2=cv1
            res1=cv(h1)
            cv1=res1[1]
          }
          dd_output=c(gmy,h1,cv1,h2,cv2)
          intcv=intcv+1
        }
        if(cv1<cv2){
          golden=cv1
          xmin[gmy,1]=golden
          xmin[gmy,2]=h1
          npar=res1[1]
        }else{
          golden=cv2
          xmin[gmy,1]=golden
          xmin[gmy,2]=h2
          npar=res2[1]
        }
        min_bandwidth_=xmin
      }

      colnames(min_bandwidth_)=c('golden', 'bandwidth')
      if(globalmin==TRUE){
        h=xmin[which.min(xmin[,1]),2]
      }else{
        h=xmin[,2]
      }

      result <- list(global_min=h,
                     local_mins=min_bandwidth_,
                     type)

    }

    if(toupper(type)=="AIC"){
      ax=0
      cx=trunc(max(as.matrix(dist(coord)))+1)
      bx=cx/2

      r=0.61803399
      tol=0.001
      cc=1-r
      h0=ax
      h3=cx


      if(abs(cx-bx)>abs(bx-ax)){
        h1=bx
        h2=bx+cc*(cx-bx)
      }else{
        h2=bx
        h1=bx-cc*(bx-ax)
      }

      bi=matrix(0,ncol(x)*n,2)

      for(i in 1:n){
        d=matrix(0,1,3)
        dist=d
        for(j in 1:n){
          if(distancekm==TRUE){
            dif=abs(coord[i,1]-coord[j,1])
            raio=acos(-1)/180
            argument=sin(coord[i,2]*raio)*sin(coord[j,2]*raio)+cos(coord[i,2]*raio)*cos(coord[j,2]*raio)*cos(dif*raio)
            if(argument>=1){
              arco=0
            }else{
              arco=acos(sin(coord[i,2]*raio)*sin(coord[j,2]*raio)+cos(coord[i,2]*raio)*cos(coord[j,2]*raio)*cos(dif*raio))
            }
            d1=arco*6371
            if(d1<=1e-3){d1=0}
          }else{
            d1=sqrt((coord[i,1]-coord[j,1])**2+(coord[i,2]-coord[j,2])**2)
          }
          if(d1!=0){
            d[1]=i
            d[2]=j
            if(distancekm==TRUE){
              d[3]=arco*6371
            }else{
              d[3]=sqrt((coord[i,1]-coord[j,1])**2+(coord[i,2]-coord[j,2])**2)
            }
            dist=rbind(dist,d)
          }
        }
        u=nrow(dist)
        w=matrix(0,u,1)
        x1=x[i,]
        y1=yc[i,]
        for(jj in 2:u){
          w[jj]=exp(-(dist[jj,3]/h1)**2)
          x1=rbind(x1,x[dist[jj,2],])
          y1=rbind(y1,yc[dist[jj,2],])
        }
        w <- as.vector(w)

        if(det(t(x1)%*%as.matrix(w*x1))==0){
          b=matrix(0,ncol(x),1)
        }else{
          b=solve(t(x1)%*%as.matrix(w*x1))%*%t(x1)%*%as.matrix(w*y1)
        }
        m1=(i-1)*ncol(x)+1
        m2=m1+(ncol(x)-1)
        bi[m1:m2,1]=i
        bi[m1:m2,2]=b
        yhat[i]=x[i,]%*%b
        ss[i]=(x1%*%solve(t(x1)%*%as.matrix(w*x1))%*%t(x1))[1]
      }
      e=yc-yhat

      betai_=matrix(t(bi[,1:2]),n, byrow = T)
      ii=seq(2,ncol(betai_),2)
      betai_=betai_[,ii]

      eta=as.matrix(x)%*%t(betai_)

      x1=x
      y1=y
      mu=linkf(eta)[,1:n]
      gmu=linkf(eta)[,(n+1):(2*n)]
      sigma2=as.numeric(t(e)%*%e)/((n-sum(ss))*(gmu*gmu))
      phi=matrix(0,1*n,1)

      for(ii in 1:n){
        for(j in 1:n){
          phi[ii]=phi[ii]+mu[j,ii]*(1-mu[j,ii])/(sigma2[j,ii]*n)
        }
      }
      phi=ifelse(phi<1,phi,phi-1)
      parameters=cbind(betai_,(phi-1))
      yhat=matrix(0,n,1)

      cv=function(h){
        phi=matrix(0,1*n,1)
        s=matrix(0,n,1)
        s_=matrix(0,n,n)
        w1=matrix(0,n,n)
        w11=matrix(0,n,n)

        for(i in 1:n){
          d=matrix(0,1,3)
          dist=d
          for(j in 1:n){
            if(distancekm==TRUE){
              dif=abs(coord[i,1]-coord[j,1])
              raio=acos(-1)/180
              argument=sin(coord[i,2]*raio)*sin(coord[j,2]*raio)+cos(coord[i,2]*raio)*cos(coord[j,2]*raio)*cos(dif*raio)
              if(argument>=1){
                arco=0
              }else{
                arco=acos(sin(coord[i,2]*raio)*sin(coord[j,2]*raio)+cos(coord[i,2]*raio)*cos(coord[j,2]*raio)*cos(dif*raio))
              }
              d1=arco*6371
              if(d1<=1e-3){d1=0}
            }else{
              d1=sqrt((coord[i,1]-coord[j,1])**2+(coord[i,2]-coord[j,2])**2)
            }
            if(d1!=0){
              d[1]=i
              d[2]=j
              if(distancekm==TRUE){
                d[3]=arco*6371
              }else{
                d[3]=sqrt((coord[i,1]-coord[j,1])**2+(coord[i,2]-coord[j,2])**2)
              }
              dist=rbind(dist,d)
            }
          }
          u=nrow(dist)
          w=matrix(0,u,1)
          x1=x[i,]
          y1=y[i,]
          for(jj in 2:u){
            w[jj]=exp(-(dist[jj,3]/h)**2)
            x1=rbind(x1,x[dist[jj,2],])
            y1=rbind(y1,y[dist[jj,2],])
          }
          parami=t(parameters[i,])
          it=0
          dif=1
          parami[ncol(parami)]=ifelse(parami[ncol(parami)]<=0,.01,parami[ncol(parami)])
          betai=t(parami[1:ncol(parami)-1])
          phii=parami[ncol(parami)]
          etaini=as.matrix(x1)%*%t(betai)

          while(abs(dif)>0.00000001 & it<maxint){
            mu=linkf(etaini)[,1]
            mu=ifelse(mu<1e-7,1e-7,mu)
            mu=ifelse(mu>=0.9999999,0.99999,mu)
            if(sum(mu>0.993)==length(mu)){
              it=maxint
            }else{
              gmu=linkf(etaini)[,2]
              ye=log(y1/(1-y1))
              mue=digamma(mu*phii)-digamma((1-mu)*phii)
              t=1/gmu
              c=t(phii%*%(trigamma(mu*phii)*mu - trigamma((1-mu)*phii)*(1-mu)))
              z=t((phii%*%trigamma(mu*phii)+phii%*%trigamma((1-mu)*phii))/(gmu*gmu))
              d=w*(trigamma(mu*phii)*mu*mu + trigamma((1-mu)*phii)*(1-mu)*(1-mu)-trigamma(phii))
              w <- as.vector(w)
              ub= phii*(t((x1*w))%*%(t*(ye-mue)))
              up=0
              for(ii in 1:length(mu)){
                up=up+(mu[ii]*(ye[ii]-mue[ii])+log(1-y1[ii])-digamma((1-mu[ii])%*%phii)+digamma(phii))*w[ii]
              }

              z <- as.vector(z)
              kbb=phii*t(x1)%*%as.matrix(x1*w*z)
              kbp=t(x1*w)%*%(t*c)
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
              etaini=as.matrix(x1)%*%t(betai)
              it=it+1
            }
          }
          xr=parami
          b=t(xr[1:ncol(x)])

          phi[i]=xr[ncol(xr)]
          s[i]=((sqrt(z)*as.matrix(x1))%*%solve(t(x1)%*%as.matrix(x1*w*z))%*%t(x1*w*sqrt(z)))[1,1]

          s_[i,]=((sqrt(z)*as.matrix(x1))%*%solve(t(x1)%*%as.matrix(x1*w*z))%*%t(x1*w*sqrt(z)))[1,]
          if(i==1){
            s1_=s_
          }else{
            if(i==n){
              s1_[i,]=c(s_[i,2:i],s_[i,1])
            }else{
              s1_[i,]=c(s_[i,2:i],s_[i,1],s_[i,(i+1):nrow(s_)])
            }
          }
          if(length(w)<n){
            w1[,i]=rbind(w,rep(0,n-nrow(w)))
          }else{
            w1[,i]=w
          }
          if(i==1){
            w11[,i]=w1[,i]
          }else{
            if(i==n){
              w11[,i]=c(w1[2:i,i],w1[1,i])
            }else{
              w11[,i]=c(w1[2:i,i],w1[1,i],w1[(i+1):nrow(w1),i])
            }
          }
          yhat[i]=x[i,]%*%t(b)
        }
        cv=t(y-yhat)%*%(y-yhat)

        s_=s1_
        vv2=sum(diag(t(s_)%*%s_))
        vv1=sum(diag(s_))
        v1=sum(s)+1
        v1=2*vv1-vv2+1
        arg=phi*mu
        arg=ifelse(arg<=0,1E-23,arg)
        lgamma1=lgamma(arg)
        arg=(1-mu)*phi
        arg=ifelse(arg<=0,1E-23,arg)
        lgamma2=lgamma(arg)
        lgamma3=(phi*mu-1)*log(y)
        lgamma4=((1-mu)*phi-1)*log(1-y)
        lnl= lgamma(phi)-lgamma1-lgamma2+lgamma3+lgamma4
        aic= 2*(v1) - 2*sum(lnl)
        aicc=aic+2*(v1)*(v1+1)/(n-v1-1)
        cv=aicc

        return(cv)
      }

      cv1=cv(h1)
      cv2=cv(h2)

      while(abs(h3-h0)>tol*(abs(h1)+abs(h2))){
        if(cv2<cv1){
          h0=h1
          h1=h2
          h2=r*h1+cc*h3
          cv1=cv2
          cv2=cv(h2)
        }else{
          h3=h2
          h2=h1
          h1=r*h2+cc*h0
          cv2=cv1
          cv1=cv(h1)
        }
      }
      dd_output=c(h1,cv1,h2,cv2)

      if(cv1<cv2){
        golden=cv1
        xmin=h1
      }else{
        golden=cv2
        xmin=h2
      }
      result=list(global_min=xmin,
                  local_mins="Only global minimum available for type=cv.",
                  type)
    }

  }

  if(toupper(method)=="FIXED_BSQ"){
    type="CV"
    ax=0
    bx=trunc(max(dist_)+1)
    if(distancekm==TRUE){
      bx=bx*111
    }

    r=0.61803399
    tol=0.1
    h0=ax
    h3=bx
    h1=bx-r*(bx-ax)
    h2=ax+r*(bx-ax)

    bi=matrix(0,ncol(x)*n,2)
    for(i in 1:n){
      seqi=matrix(i,n,1)
      dist=cbind(seqi,seq,dist_[,i])
      if(distancekm==TRUE){
        dist[,3]=dist[,3]*111
      }
      u=nrow(dist)
      w=matrix(0,u,1)
      for(jj in 1:u){
        w[jj]=(1-(dist[jj,3]/h1)^2)^2
      }

      w[i]=0
      w <- as.vector(w)

      w[dist[,3]<=h1]=0

      if(det(t(x)%*%as.matrix(w*x))==0){
        b=matrix(0,ncol(x),1)
      }else{
        b=solve(t(x)%*%as.matrix(w*x))%*%t(x)%*%as.matrix(w*yc)
      }
      m1=(i-1)*ncol(x)+1
      m2=m1+(ncol(x)-1)
      bi[m1:m2,1]=i
      bi[m1:m2,2]=b
      yhat[i]=x[i,]%*%b
      ss[i]=(x[i,]%*%solve(t(x)%*%as.matrix(w*x))%*%t(x))[1]
    }

    e=yc-yhat

    betai_=matrix(t(bi[,1:2]),n, byrow = T)
    ii=seq(2,ncol(betai_),2)
    betai_=betai_[,ii]

    eta=as.matrix(x)%*%t(betai_)

    mu=linkf(eta)[,1:n]
    gmu=linkf(eta)[,(n+1):(2*n)]
    sigma2=as.numeric(t(e)%*%e)/((n-sum(ss))*(gmu*gmu))
    phi=matrix(0,1*n,1)
    for(ii in 1:n){
      for(j in 1:n){
        phi[ii]=phi[ii]+mu[j,ii]*(1-mu[j,ii])/(sigma2[j,ii]*n)
      }
    }
    phi=ifelse(phi<1,phi,phi-1)
    parameters=cbind(betai_,(phi-1))
    yhat=matrix(0,n,1)

    cv <- function(h){
      phi=matrix(0,1*n,1)
      s=matrix(0,n,1)
      s_=matrix(0,n,n)
      w1=matrix(0,n,n)
      w11=matrix(0,n,n)

      for(i in 1:n){
        seqi=matrix(i,n,1)
        dist=cbind(seqi,seq,dist_[,i])
        if(distancekm==TRUE){
          dist[,3]=dist[,3]*111
        }
        u=nrow(dist)
        w=matrix(0,u,1)
        for(jj in 1:u){
          w[jj]=(1-(dist[jj,3]/h)^2)^2
          w[i]=0
        }
        position=which(dist[,3]>h)
        w[position]=0
        parami=t(parameters[i,])
        it=0
        dif=1
        parami[ncol(parami)]=ifelse(parami[ncol(parami)]<=0,.01,parami[ncol(parami)])
        betai=t(parami[1:ncol(parami)-1])
        phii=parami[ncol(parami)]
        etaini=as.matrix(x)%*%t(betai)

        while(abs(dif)>0.00000001 & it<maxint){
          mu=linkf(etaini)[,1]
          mu=ifelse(mu<1e-7,1e-7,mu)
          mu=ifelse(mu>=0.9999999,0.99999,mu)
          if(sum(mu>0.993)==length(mu)){
            it=maxint
          }else{
            gmu=linkf(etaini)[,2]
            ye=log(y/(1-y))
            mue=digamma(mu*phii)-digamma((1-mu)*phii)
            t=1/gmu
            c=t(phii%*%(trigamma(mu*phii)*mu - trigamma((1-mu)*phii)*(1-mu)))
            z=t((phii%*%trigamma(mu*phii)+phii%*%trigamma((1-mu)*phii))/(gmu*gmu))
            d=w*(trigamma(mu*phii)*mu*mu + trigamma((1-mu)*phii)*(1-mu)*(1-mu)-trigamma(phii))
            w <- as.vector(w)
            ub= phii*(t((x*w))%*%(t*(ye-mue)))
            up=0
            for(ii in 1:length(mu)){
              up=up+(mu[ii]*(ye[ii]-mue[ii])+log(1-y[ii])-digamma((1-mu[ii])%*%phii)+digamma(phii))*w[ii]
            }
            z <- as.vector(z)
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
        b=t(xr[1:ncol(x)])

        s_[i,]=(x[i,]%*%solve(t(x)%*%as.matrix(x*w*z))%*%t(x*w*z))[1,]

        yhat[i]=x[i,]%*%t(b)
      }
      cv=t(y-yhat)%*%(y-yhat)
      v1=NULL

      res=c(cv,v1)
      return(res)
    }

    ax=0
    bx=trunc(max(dist_)+1)
    if(distancekm==TRUE){
      bx=bx*111
    }
    r=0.61803399
    tol=0.1

    if(globalmin==FALSE){
      lower=ax
      upper=bx
      xmin=matrix(0,1,2)
      end_iter=1
    }else{
      lower=cbind(ax,(1-r)*bx,r*bx)
      upper=cbind((1-r)*bx,r*bx,bx)
      xmin=matrix(0,3,2)
      end_iter=3
    }

    min_bandwidth_=NULL
    for(gmy in 1:end_iter){
      ax1=lower[gmy]
      bx1=upper[gmy]
      h0=ax1
      h3=bx1
      h1=bx1-r*(bx1-ax1)
      h2=ax1+r*(bx1-ax1)

      rnd_=1
      res1=cv(h1)
      cv1=res1[1]
      res2=cv(h2)
      cv2=res2[1]

      dd_output=c(gmy,h1,cv1,h2,cv2)

      intcv=1
      while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & intcv<200){
        if(cv2<cv1){
          h0=h1
          h1=h3-r*(h3-h0)
          h2=h0+r*(h3-h0)
          cv1=cv2
          res2=cv(h2)
          cv2=res2[1]
        }else{
          h3=h2
          h1=h3-r*(h3-h0)
          h2=h0+r*(h3-h0)
          cv2=cv1
          res1=cv(h1)
          cv1=res1[1]
        }
        dd_output=c(gmy,h1,cv1,h2,cv2)
        intcv=intcv+1
      }
      if(cv1<cv2){
        golden=cv1
        xmin[gmy,1]=golden
        xmin[gmy,2]=h1
        npar=res1[1]
      }else{
        golden=cv2
        xmin[gmy,1]=golden
        xmin[gmy,2]=h2
        npar=res2[1]
      }
      min_bandwidth_=xmin
    }

    colnames(min_bandwidth_)=c('golden', 'bandwidth')
    if(globalmin==TRUE){
      h=xmin[which.min(xmin[,1]),2]
    }else{
      h=xmin[,2]
    }

    result <- list(global_min=h,
                   local_mins=min_bandwidth_,
                   type)
  }

  if(toupper(method)=="ADAPTIVE_BSQ"){
    if(toupper(type)=="CV"){
      ax=5
      bx=n
      r=0.61803399
      tol=0.1
      h0=ax
      h3=bx
      h1=bx-r*(bx-ax)
      h2=ax+r*(bx-ax)

      bi=matrix(0,ncol(x)*n,2)
      for(i in 1:n){
        seqi=matrix(i,n,1)
        dist=cbind(seqi,seq,dist_[,i])
        if(distancekm==TRUE){
          dist[,3]=dist[,3]*111
        }
        u=nrow(dist)
        w=matrix(0,u,1)
        for(jj in 1:u){
          w[jj]=exp(-0.5*(dist[jj,3]/h1)^2)
        }

        w[i]=0

        dist <- dist[order(dist[,3]),]
        dist <- cbind(dist, c(1:nrow(dist)))
        w=matrix(0,n,2)
        hn=dist[h1,3]
        for(jj in 1:n){
          if(dist[jj,4]<=h1){
            w[jj,1]=(1-(dist[jj,3]/hn)^2)^2
          }else{
            w[jj,1]=0
          }
          w[jj,2]=dist[jj,2]
        }

        w[which(w[,2]==i)]=0
        w=w[order(w[,2]),]
        w=w[,1]
        if(det(t(x)%*%as.matrix(w*x))==0){
          b=matrix(0,ncol(x),1)
        }else{
          b=solve(t(x)%*%as.matrix(w*x))%*%t(x)%*%as.matrix(w*yc)
        }
        m1=(i-1)*ncol(x)+1
        m2=m1+(ncol(x)-1)
        bi[m1:m2,1]=i
        bi[m1:m2,2]=b
        yhat[i]=x[i,]%*%b
        ss[i]=(x[i,]%*%solve(t(x)%*%as.matrix(w*x))%*%t(x))[1]
      }
      e=yc-yhat

      betai_=matrix(t(bi[,1:2]),n, byrow = T)
      ii=seq(2,ncol(betai_),2)
      betai_=betai_[,ii]

      eta=as.matrix(x)%*%t(betai_)

      mu=linkf(eta)[,1:n]
      gmu=linkf(eta)[,(n+1):(2*n)]
      sigma2=as.numeric(t(e)%*%e)/((n-sum(ss))*(gmu*gmu))
      phi=matrix(0,1*n,1)
      for(ii in 1:n){
        for(j in 1:n){
          phi[ii]=phi[ii]+mu[j,ii]*(1-mu[j,ii])/(sigma2[j,ii]*n)
        }
      }
      phi=ifelse(phi<1,phi,phi-1)
      parameters=cbind(betai_,(phi-1))
      yhat=matrix(0,n,1)

      cv=function(h){
        phi=matrix(0,1*n,1)
        s=matrix(0,n,1)
        s_=matrix(0,n,n)
        w1=matrix(0,n,n)
        w11=matrix(0,n,n)

        for(i in 1:n){
          seqi=matrix(i,n,1)
          dist=cbind(seqi,seq,dist_[,i])

          if(distancekm==TRUE){
            dist[,3]=dist[,3]*111
          }
          u=nrow(dist)
          w=matrix(0,u,1)
          for(jj in 1:u){
            w[jj]=exp(-0.5*(dist[jj,3]/h)^2)
            w[i]=0
          }
          dist <- dist[order(dist[,3]),]
          dist <- cbind(dist, c(1:nrow(dist)))
          w=matrix(0,n,2)
          hn=dist[h,3]
          for(jj in 1:n){
            if(dist[jj,4]<=h){
              w[jj,1]=(1-(dist[jj,3]/hn)^2)^2
            }else{
              w[jj,1]=0
            }
            w[jj,2]=dist[jj,2]
          }
          position2=w[which(w[,1]>0),2]

          w[which(w[,2]==i)]=0

          w=w[order(w[,2]),]
          w=w[,1]
          parami=t(parameters[i,])
          it=0
          dif=1
          parami[ncol(parami)]=ifelse(parami[ncol(parami)]<=0,.01,parami[ncol(parami)])
          betai=t(parami[1:ncol(parami)-1])
          phii=parami[ncol(parami)]
          etaini=as.matrix(x)%*%t(betai)

          while(abs(dif)>0.00000001 & it<maxint){
            mu=linkf(etaini)[,1]
            mu=ifelse(mu<1e-7,1e-7,mu)
            mu=ifelse(mu>=0.9999999,0.99999,mu)
            if(sum(mu>0.993)==length(mu)){
              it=maxint
            }else{
              gmu=linkf(etaini)[,2]
              ye=log(y/(1-y))
              mue=digamma(mu*phii)-digamma((1-mu)*phii)
              t=1/gmu
              c=t(phii%*%(trigamma(mu*phii)*mu - trigamma((1-mu)*phii)*(1-mu)))
              z=t((phii%*%trigamma(mu*phii)+phii%*%trigamma((1-mu)*phii))/(gmu*gmu))
              d=w*(trigamma(mu*phii)*mu*mu + trigamma((1-mu)*phii)*(1-mu)*(1-mu)-trigamma(phii))
              w <- as.vector(w)
              ub= phii*(t((x*w))%*%(t*(ye-mue)))
              up=0
              for(ii in 1:length(mu)){
                up=up+(mu[ii]*(ye[ii]-mue[ii])+log(1-y[ii])-digamma((1-mu[ii])%*%phii)+digamma(phii))*w[ii]
              }
              z <- as.vector(z)
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
          b=t(xr[1:ncol(x)])

          s_[i,]=(x[i,]%*%solve(t(x)%*%as.matrix(x*w*z))%*%t(x*w*z))[1,]

          yhat[i]=x[i,]%*%t(b)
        }
        cv=t(y-yhat)%*%(y-yhat)

        v1=NULL


        res=c(cv,v1)
        return(res)
      }

      ax=5
      bx=n
      r=0.61803399
      tol=0.1

      if(globalmin==FALSE){
        lower=ax
        upper=bx
        xmin=matrix(0,1,2)
        end_iter=1
      }else{
        lower=cbind(ax,(1-r)*bx,r*bx)
        upper=cbind((1-r)*bx,r*bx,bx)
        xmin=matrix(0,3,2)
        end_iter=3
      }

      min_bandwidth_=NULL
      for(gmy in 1:end_iter){
        ax1=lower[gmy]
        bx1=upper[gmy]
        h0=ax1
        h3=bx1
        h1=bx1-r*(bx1-ax1)
        h2=ax1+r*(bx1-ax1)

        rnd_=1
        res1=cv(h1)
        cv1=res1[1]
        res2=cv(h2)
        cv2=res2[1]

        dd_output=c(gmy,h1,cv1,h2,cv2)

        intcv=1
        while(abs(h3-h0) > tol*(abs(h1)+abs(h2)) & intcv<200){
          if(cv2<cv1){
            h0=h1
            h1=h3-r*(h3-h0)
            h2=h0+r*(h3-h0)
            cv1=cv2
            res2=cv(h2)
            cv2=res2[1]
          }else{
            h3=h2
            h1=h3-r*(h3-h0)
            h2=h0+r*(h3-h0)
            cv2=cv1
            res1=cv(h1)
            cv1=res1[1]
          }
          dd_output=c(gmy,h1,cv1,h2,cv2)
          intcv=intcv+1
        }
        if(cv1<cv2){
          golden=cv1
          xmin[gmy,1]=golden
          xmin[gmy,2]=h1
          npar=res1[1]
          xmin[gmy,2]=floor(h1)
        }else{
          golden=cv2
          xmin[gmy,1]=golden
          xmin[gmy,2]=h2
          npar=res2[1]
          xmin[gmy,2]=floor(h2)
        }
        min_bandwidth_=xmin
      }
      colnames(min_bandwidth_)=c('golden', 'bandwidth')

      if(globalmin==TRUE){
        h=xmin[which.min(xmin[,1]),2]
      }else{
        h=xmin[,2]
      }
      result <- list(global_min=h,
                     local_mins=min_bandwidth_,
                     type)
    }

    if(toupper(type)=="AIC"){
      ax=5
      cx=nrow(coord)
      bx=trunc(cx/2)

      r=0.61803399
      tol=0.001
      cc=1-r
      h0=ax
      h3=cx

      if(abs(cx-bx)>abs(bx-ax)){
        h1=bx
        h2=bx+cc*(cx-bx)
      }else{
        h2=bx
        h1=bx-cc*(bx-ax)
      }

      bi=matrix(0,ncol(x)*n,2)

      for(i in 1:n){
        d=matrix(0,1,3)
        dist=d
        for(j in 1:n){
          if(distancekm==TRUE){
            dif=abs(coord[i,1]-coord[j,1])
            raio=acos(-1)/180
            argument=sin(coord[i,2]*raio)*sin(coord[j,2]*raio)+cos(coord[i,2]*raio)*cos(coord[j,2]*raio)*cos(dif*raio)
            if(argument>=1){
              arco=0
            }else{
              arco=acos(sin(coord[i,2]*raio)*sin(coord[j,2]*raio)+cos(coord[i,2]*raio)*cos(coord[j,2]*raio)*cos(dif*raio))
            }
            d1=arco*6371
            if(d1<=1e-3){d1=0}
          }else{
            d1=sqrt((coord[i,1]-coord[j,1])**2+(coord[i,2]-coord[j,2])**2)
          }
          if(d1!=0){
            d[1]=i
            d[2]=j
            if(distancekm==TRUE){
              d[3]=arco*6371
            }else{
              d[3]=sqrt((coord[i,1]-coord[j,1])**2+(coord[i,2]-coord[j,2])**2)
            }
            dist=rbind(dist,d)
          }
        }
        u=nrow(dist)
        w=matrix(0,u,1)
        x1=x[i,]
        y1=yc[i,]

        for(jj in 2:u){
          w[jj]=exp(-(dist[jj,3]/h1)**2)
          x1=rbind(x1,x[dist[jj,2],])
          y1=rbind(y1,yc[dist[jj,2],])
        }
        w <- as.numeric(w)

        x1=x[i,]
        y1=yc[i,]


        dist <- dist[order(dist[,3]),]
        dist <- cbind(dist, c(1:nrow(dist)))
        w=matrix(0,n,2)
        hn=dist[h1,3]
        for(jj in 2:(n-1)){
          if(dist[jj,4]<=h1){
            w[jj,1]=(1-(dist[jj,3]/hn)**2)**2
          }else{
            w[jj,1]=0
          }
          w[jj,2]=dist[jj,2]
        }
        position=w[which(w[,1]>0),2]
        w=c(0,w[position,1])
        x1=rbind(x1,x[position,])
        y1=c(y1,yc[position,])

        if(det(t(x1)%*%as.matrix(w*x1))==0){
          b=matrix(0,ncol(x),1)
        }else{
          b=solve(t(x1)%*%as.matrix(w*x1))%*%t(x1)%*%as.matrix(w*y1)
        }
        m1=(i-1)*ncol(x)+1
        m2=m1+(ncol(x)-1)
        bi[m1:m2,1]=i
        bi[m1:m2,2]=b
        yhat[i]=x[i,]%*%b
        ss[i]=(x1%*%solve(t(x1)%*%as.matrix(w*x1))%*%t(x1))[1]
      }

      e=yc-yhat

      betai_=matrix(t(bi[,1:2]),n, byrow = T)
      ii=seq(2,ncol(betai_),2)
      betai_=betai_[,ii]

      eta=as.matrix(x)%*%t(betai_)

      x1=x
      y1=y
      mu=linkf(eta)[,1:n]
      gmu=linkf(eta)[,(n+1):(2*n)]
      sigma2=as.numeric(t(e)%*%e)/((n-sum(ss))*(gmu*gmu))
      phi=matrix(0,1*n,1)
      for(ii in 1:n){
        for(j in 1:n){
          phi[ii]=phi[ii]+mu[j,ii]*(1-mu[j,ii])/(sigma2[j,ii]*n)
        }
      }
      phi=ifelse(phi<1,phi,phi-1)
      parameters=cbind(betai_,(phi-1))
      yhat=matrix(0,n,1)

      cv=function(h){

        max_like=function(param){
          betai2=param[1:length(param)-1]
          phii2=param[length(param)]
          etai2=as.matrix(x1)%*%betai2
          mu2=linkf(etai2)[,1]
          lgamma1=t(lgamma(phii2%*%mu2))
          arg=(1-mu2)*phii2
          arg=ifelse(arg<=0,1E-23,arg)
          lgamma2=lgamma(arg)
          lgamma3=(phii2*mu2-1)*log(y1)
          lgamma4=((1-mu2)*phii2-1)*log(1-y1)
          lnl=0
          n=length(y1)
          for(ll in 1:n){
            lnl=lnl+(lgamma(phii2)-lgamma1[ll]-lgamma2[ll]+lgamma3[ll]+lgamma4[ll])%*%w[ll]
          }
          return(lnl)
        }

        phi=matrix(0,1*n,1)
        s=matrix(0,n,1)
        s_=matrix(0,n,n)
        w1=matrix(0,n,n)
        w11=matrix(0,n,n)

        for(i in 1:n){
          d=matrix(0,1,3)
          dist=d
          for(j in 1:n){
            if(distancekm==TRUE){
              dif=abs(coord[i,1]-coord[j,1])
              raio=acos(-1)/180
              argument=sin(coord[i,2]*raio)*sin(coord[j,2]*raio)+cos(coord[i,2]*raio)*cos(coord[j,2]*raio)*cos(dif*raio)
              if(argument>=1){
                arco=0
              }else{
                arco=acos(sin(coord[i,2]*raio)*sin(coord[j,2]*raio)+cos(coord[i,2]*raio)*cos(coord[j,2]*raio)*cos(dif*raio))
              }
              d1=arco*6371
              if(d1<=1e-3){d1=0}
            }else{
              d1=sqrt((coord[i,1]-coord[j,1])**2+(coord[i,2]-coord[j,2])**2)
            }
            if(d1!=0){
              d[1]=i
              d[2]=j
              if(distancekm==TRUE){
                d[3]=arco*6371
              }else{
                d[3]=sqrt((coord[i,1]-coord[j,1])**2+(coord[i,2]-coord[j,2])**2)
              }
              dist=rbind(dist,d)
            }
          }
          u=nrow(dist)
          w=matrix(0,u,1)
          x1=x[i,]
          y1=y[i,]
          for(jj in 2:u){
            w[jj]=exp(-(dist[jj,3]/h)**2)
            x1=rbind(x1,x[dist[jj,2],])
            y1=rbind(y1,y[dist[jj,2],])
          }
          x1=x[i,]
          y1=y[i,]
          dist <- dist[order(dist[,3]),]
          dist <- cbind(dist, c(1:nrow(dist)))

          w=matrix(0,n,2)
          hn=dist[h,3]
          for(jj in 2:(n-1)){
            if(dist[jj,4]<=h){
              w[jj,1]=(1-(dist[jj,3]/hn)**2)**2
            }else{
              w[jj,1]=0
            }
            w[jj,2]=dist[jj,2]
          }

          position=w[which(w[,1]>0),2]
          position2=c(i,position)
          w=c(0,w[position,1])

          x1=rbind(x1,x[position,])
          y1=c(y1,y[position,])
          parami=t(parameters[i,])
          it=0
          dif=1
          parami[ncol(parami)]=ifelse(parami[ncol(parami)]<=0,.01,parami[ncol(parami)])
          betai=t(parami[1:ncol(parami)-1])
          phii=parami[ncol(parami)]
          etaini=as.matrix(x1)%*%t(betai)

          while(abs(dif)>0.00000001 & it<maxint){
            mu=linkf(etaini)[,1]
            mu=ifelse(mu<1e-7,1e-7,mu)
            mu=ifelse(mu>=0.9999999,0.99999,mu)
            if(sum(mu>0.993)==length(mu)){
              it=maxint
            }else{
              gmu=linkf(etaini)[,2]
              ye=log(y1/(1-y1))
              mue=digamma(mu*phii)-digamma((1-mu)*phii)
              t=1/gmu
              c=t(phii%*%(trigamma(mu*phii)*mu - trigamma((1-mu)*phii)*(1-mu)))
              z=t((phii%*%trigamma(mu*phii)+phii%*%trigamma((1-mu)*phii))/(gmu*gmu))
              d=w*(trigamma(mu*phii)*mu*mu + trigamma((1-mu)*phii)*(1-mu)*(1-mu)-trigamma(phii))

              w <- as.vector(w)
              ub= phii*(t((x1*w))%*%(t*(ye-mue)))
              up=0
              for(ii in 1:length(mu)){
                up=up+(mu[ii]*(ye[ii]-mue[ii])+log(1-y1[ii])-digamma((1-mu[ii])%*%phii)+digamma(phii))*w[ii]
              }

              z <- as.vector(z)
              kbb=phii*t(x1)%*%as.matrix(x1*w*z)
              kbp=t(x1*w)%*%(t*c)
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
              etaini=as.matrix(x1)%*%t(betai)
              it=it+1
            }
          }
          xr=parami
          b=t(xr[1:ncol(x)])

          phi[i]=xr[ncol(xr)]
          s[i]=((sqrt(z)*as.matrix(x1))%*%solve(t(x1)%*%as.matrix(x1*w*z), tol=0)%*%t(x1*w*sqrt(z)))[1,1]
          s_p_=t(((sqrt(z)*as.matrix(x1))%*%solve(t(x1)%*%as.matrix(x1*w*z), tol=0)%*%t(x1*w*sqrt(z)))[1,])
          for(ksp in 1:ncol(s_p_)){
            s_[i,position2[ksp]]=s_p_[ksp]
          }
          if(length(w)<n){
            w1[,i]=c(w,rep(0,n-length(w)))
          }else{
            w1[,i]=w
          }
          if(i==1){
            w11[,i]=w1[,i]
          }else{
            if(i==n){
              w11[,i]=c(w1[2:i,i],w1[1,i])
            }else{
              w11[,i]=c(w1[2:i,i],w1[1,i],w1[(i+1):nrow(w1),i])
            }
          }
          yhat[i]=x[i,]%*%t(b)
        }

        mu2=mu
        mu=matrix(0,n,1)
        for(ksp in 1:ncol(s_p_)){
          mu[position2[ksp]]=mu2[ksp]
        }

        vv2=sum(diag(t(s_)%*%s_))
        vv1=sum(diag(s_))
        v1=sum(s)+1
        v1=2*vv1-vv2+1
        arg=phi*mu
        arg=ifelse(arg<=0,1E-23,arg)
        lgamma1=lgamma(arg)
        arg=(1-mu)*phi
        arg=ifelse(arg<=0,1E-23,arg)
        lgamma2=lgamma(arg)
        lgamma3=(phi*mu-1)*log(y)
        lgamma4=((1-mu)*phi-1)*log(1-y)
        lnl= lgamma(phi)-lgamma1-lgamma2+lgamma3+lgamma4
        aic= 2*(v1) - 2*sum(lnl)
        aicc=aic+2*(v1)*(v1+1)/(n-v1-1)

        return(aicc)
      }

      cv1=cv(h1)
      cv2=cv(h2)


      while(abs(h3-h0)>tol*(abs(h1)+abs(h2))){
        if(cv2<cv1){
          h0=h1
          h1=h2
          h2=r*h1+cc*h3
          cv1=cv2
          cv2=cv(h2)
        }else{
          h3=h2
          h2=h1
          h1=r*h2+cc*h0
          cv2=cv1
          cv1=cv(h1)
        }
      }
      dd_output=c(h1,cv1,h2,cv2)

      if(cv1<cv2){
        golden=cv1
        xmin=floor(h1)
      }else{
        golden=cv2
        xmin=floor(h2)
      }

      result=list(global_min=xmin,
                  local_mins="Only global minimum available for type=cv.",
                  type)
    }
  }

  return(result)
}

