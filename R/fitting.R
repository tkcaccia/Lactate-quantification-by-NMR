
voigt =   function(A, x, x0, fwhm, eta){
    G= exp( -(x-x0)^2 * 4 * log(2)  / (fwhm^2))
    L= (0.5 * fwhm)^2 / ((x-x0)^2 + (0.5 * fwhm)^2)
    A*((1-eta) * G + eta * L)
  }

doublet =
  function(x1,height,ppm,fwhm,eta,J,p1=1,p2=1,...){
    fit=voight(p1*height,x1,ppm-J/2,fwhm,eta)+
      voight(p2*height,x1,ppm+J/2,fwhm,eta)
    fit
  }

TSP.500 =
  function(x1,height,ppm,fwhm,eta) {
    fitting=voigt(height,x1,ppm,fwhm,eta)
    fitting=fitting+voigt(height/35,x1,ppm-0.00653,fwhm,eta)
    fitting=fitting+voigt(height/35,x1,ppm+0.00653,fwhm,eta)
    fitting
  }


pick.peaks =   function (x, span) 
  {
    span.width <- span * 2 + 1
    loc.max <- span.width + 1 - apply(embed(x, span.width), 1, which.max)
    loc.max[loc.max == 1 | loc.max == span.width] <- NA
    pks <- loc.max + 0:(length(loc.max) - 1)
    unique(pks[!is.na(pks)])
  }

search.siglet.2 = function(x1,y1,...){
    pp=pick.peaks(y1,5)
    if(length(pp)>0){
      ppm=x1[pp][which.max(y1[pp])]
    }else{
      ppm=NA
    }
    ppm
  }


search.doublet =   function(x1,y1,J,error,simila=0.5,...){
    ll=length(x1)
    if(ll>1000){
      
      l1=predict(loess(y1~x1,span = 0.02),x1)
      
      pp=pick.peaks(l1,50)
    }
    if(ll>300 & ll<=1000){
      l1=predict(loess(y1~x1,span = 0.1),x1)
      
      pp=pick.peaks(l1,50)
    }
    if(ll>200 & ll<=300){
      l1=predict(loess(y1~x1,span = 0.1),x1)
      
      pp=pick.peaks(l1,40)
    }
    if(ll>100 & ll<=200){
      l1=predict(loess(y1~x1,span = 0.1),x1)
      
      pp=pick.peaks(l1,20)
    }
    if(ll>50 & ll<=100){
      l1=predict(loess(y1~x1,span = 0.2),x1)
      
      pp=pick.peaks(l1,10)
    }
    if(ll>0 & ll<=50){
      l1=predict(loess(y1~x1,span = 0.2),x1)
      
      pp=pick.peaks(l1,5)
    }
    
    if(J<0.005){
      pp=pick.peaks(l1,5)
    }
    if(length(pp)!=0){
      M=as.matrix(dist(x1[pp]))
      W=which(M>(J-error) & M<(J+error) & M!=0,arr.ind=T)
      similar=abs(log(abs(y1[pp[W[,1]]]/y1[pp[W[,2]]])))<simila
      if(sum(similar)>1){
        W=W[similar,]
        W=W[which.max(rowSums(cbind(y1[pp[W[,1]]],y1[pp[W[,2]]]))),]
        ppm=mean(x1[pp[W]])
      }else{
        ppm=NA
      }
    }else{
      ppm=NA
    }
    ppm
}





TSP_ref.500 = 
  function(x,y){
    
    
    fr1 <- function(z) {   
      fwhmopt <- z[1]
      etaopt <- z[2]
      a=y1-TSP.500(x1,height,ppm,fwhmopt,etaopt)
      sum(a[sel]*a[sel])
    }
    fr2 <- function(z) {   
      heightopt <- z[1]
      ppmopt <- z[2]
      a=y1-TSP.500(x1,heightopt,ppmopt,fwhm,eta)
      sum(a[sel]*a[sel])
    }
    
    range=c(-0.1,0.1)
    sel=x>range[1] & x<range[2]
    x1=x[sel]
    y1=y[sel]
    pp=which.max(y1)
    
    ppm=x1[pp]
    height=y1[pp]
    
    range=c(-0.02,0.02)+ppm
    sel=x>range[1] & x<range[2]
    x1=x[sel]
    y1=y[sel]
    pp=which.max(y1)
    
    ppm=x1[pp]
    height=y1[pp]
    
    sel=x1>(ppm-0.02) & x1<(ppm+0.02)
    
    
    o=optim(c(0.002,1), fr1,method = "L-BFGS-B",
            lower = c(0.0001,0.8), 
            upper = c(0.01,2))    
    
    
    fwhm=abs(o$par[1])
    eta=o$par[2]
    #  o=optim(c(height,ppm), fr2)
    
    o=optim(c(height,ppm), fr2,method = "L-BFGS-B",
            lower = c(height,ppm)-c(height*0.005,0.0001), 
            upper = c(height,ppm)+c(height*0.05,0.0001))
    height=abs(o$par[1])
    ppm=o$par[2]
    
    area=sum(TSP.500(x1,height,ppm,fwhm,eta))
    RMSD=sqrt(o$value/length(x1))
    return(list(ppm=ppm,height=height,fwhm=fwhm,eta=eta,area=area,RMSD=RMSD))
    
  }




