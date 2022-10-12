
sumFunction =
  function(FUN, xv, height, ppm, fwhm, eta, J) {
    s=rep(0,length(xv))
    for(i in length(FUN)){
      s=s+FUN[[i]](xv, height[i], ppm[i], fwhm[i], eta[i], J[i])
    }
    s
  }

conta =
  function(v){
    l=NULL
    for(i in 1:length(v))
      l[i]=length(v[[i]])
    l
  }

param_deconvolution =
  function(xv,zv,size,ppm,FUN,fwhm_min,fwhm_max,
           eta_min,eta_max,...){
    link=as.numeric(as.factor(link))
    functions=unique(link)
    n_fun=length(functions)
    
    sel=rep(FALSE,length(xv))
    height=NULL
    fwhm=NULL
    eta=NULL
    nnn=length(FUN)
    
    a1o=ppm+size
    a2o=ppm-size
    
    sel_i=list()
    
    for(i in 1:nnn){
      
      sel_temp=rep(FALSE,length(xv))
      if(!is.na(a1o[i]) & !is.na(a2o[i])){
        sel_i[[i]]=(xv<a1o[i] & xv>a2o[i])
        sel_temp=sel_temp | sel_i[[i]]
        
      }
      
      height[i]=abs(max(zv[sel_i[[i]]]))
      fwhm[i]=fwhm_ref
      eta[i]=eta_ref
      sel=sel | sel_temp
    }
    height[is.na(height)]=0
    sel=which(sel)
    fwhm_min=fwhm*fwhm_min
    fwhm_max=fwhm*fwhm_max
    fwhm=apply(rbind(fwhm,fwhm_min),2,max)
    fwhm=apply(rbind(fwhm,fwhm_max),2,min)
    
    eta_min=eta*eta_min
    eta_max=eta*eta_max
    eta=apply(rbind(eta,eta_min),2,max)
    eta=apply(rbind(eta,eta_max),2,min)
    
    height_min=height*0.2
    height_max=height*1.5
    
    ppm_min=ppm-shift
    ppm_max=ppm+shift
    
    
    npeaks=conta(ppm)
    
    
    
    for(iii in 1:5){
      
      who=unlist((lapply(ppm,function(x) !any(is.na(x))))) |   height_max!=0
      wwho=which(who)
      lwwho=length(wwho)
      fr <- function(zzz) {   
        vv=rep(0,length(xv))
        for(i in 1:lwwho){
          
          vv=vv+sumFunction(FUN[wwho][i],
                            xv, 
                            protons[wwho][i]*height[wwho][i], 
                            ppm[wwho][i], zzz[i], zzz[i+lwwho],
                            J[i])
          
        }
        a=zv-vv
        a=a*a
        sum(abs(a[sel]))
      }
      if(!all(unlist((lapply(ppm,function(x) any(is.na(x))))))){ 
        oo=optim(c(fwhm[wwho],eta[wwho]), fr,method = "L-BFGS-B",
                 lower = c(fwhm_min[wwho],eta_min[wwho]), 
                 upper = c(fwhm_max[wwho],eta_max[wwho])+0.000001)
        
        fwhm[wwho]=abs(oo$par[1:length(wwho)])
        eta[wwho]=oo$par[length(wwho)+(1:length(wwho))]
        ########################
      }
      
      
      
      
      
      fr <- function(z) {   
        vv=rep(0,length(xv))
        for(i in 1:lwwho){
          vv=vv+sumFunction(FUN[wwho][i],
                            xv,  
                            height[wwho][i]*protons[wwho][i], 
                            z[i], fwhm[i], eta[i],J[i])
        }
        a=zv-vv
        a=a*a*zv
        sum(abs(a[sel]))  
      }
      if(!all(unlist((lapply(ppm,function(x) any(is.na(x))))))){ 
        oo=optim(ppm[wwho], fr,method = "L-BFGS-B",
                 lower = ppm_min[wwho], 
                 upper = ppm_max[wwho]+0.000001)
        for(it in 1:length(ppm[wwho])){
          start=1-sum(npeaks[wwho][length(npeaks[wwho]):it])+sum(npeaks[wwho])
          end=sum(npeaks[wwho][1:it])
          ppm[wwho][[it]]=oo$par[start:end]
          
        }
        
      }
      
      
      
      fr <- function(zzz) {   
        vv=rep(0,length(xv))
        for(j in 1:sum(peak.selection[wwho])){
          pos=(1:lwwho)[peak.selection[wwho]]
          i=pos[j]
          wi=which(link[i]==link)
          if(peak.selection[i]){
            vv=vv+sumFunction(FUN[wwho][wi],xv, zzz[j]*protons[wwho][wi], ppm[wwho][wi], fwhm[wi], eta[wi],J[wi])
          }
        }
        a=zv-vv
        a=a*a*zv
        
        sum(abs(a[sel]))    
      }
      if(!all(unlist((lapply(ppm,function(x) all(is.na(x))))))){ 
        
        height[height<=0]=0.01
        height_max[height_max<=0]=0.01
        
        oo=optim(c(height[wwho][peak.selection[wwho]]), fr,method = "L-BFGS-B",
                 lower = c(height_min[wwho][peak.selection[wwho]]), 
                 upper = c(height_max[wwho][peak.selection[wwho]])+0.000001) 
        height[wwho][peak.selection[wwho]]=oo$par
      }
      
      for(i in 1:nnn){
        if(!is.na(ppm[i])){
          if(!peak.selection[i]){
            height[i]=NA
            height[i]=mean(height[link[i]==link],na.rm=TRUE)
          }
        }
      }
      
      ##################VINCOLARE ###########################
      
      if(!all(unlist((lapply(ppm,function(x) any(is.na(x))))))){ 
        height[height<=0]=0.01
        temp_height=height
        for(hh in 1:nnn){
          if(!all(is.na(ppm[[hh]]))){
            fitto=rep(0,length(xv))
            for(i in (1:nnn)[-hh])
              
              if(!all(is.na(ppm[i]))){
                
                
                fitto=fitto+sumFunction(FUN[wwho][i],xv, temp_height[wwho][i], 
                                        ppm[wwho][i], fwhm[wwho][i], eta[wwho][i],J[wwho][i])
              }
            
            
            FF=sumFunction(FUN[wwho][hh],xv, 
                           temp_height[wwho][hh], 
                           ppm[wwho][hh], fwhm[wwho][hh], 
                           eta[wwho][hh],J[wwho][hh])
            
            
            pp_FF=pick.peaks(FF,30)
            if(length(pp_FF)>2) {
              pp_FF=pp_FF[order(FF[pp_FF],decreasing = T)[1:2]]
              
            }
            if(length(pp_FF)==0) {
              pp_FF=which.max(FF)
              
            }
            ln=length(xv)
            njj=nrow(as.matrix(pp_FF))
            temp_pp=matrix(nrow=njj,ncol=5)
            temp_pp[,1]=pp_FF
            temp_pp[,2]=pp_FF-1
            temp_pp[,3]=pp_FF-2
            temp_pp[,4]=pp_FF+1
            temp_pp[,5]=pp_FF+2
            temp_pp[temp_pp<1]=1
            temp_pp[temp_pp>ln]=ln
            
            mi = NULL
            for( jjj in 1:njj)
              mi[jjj]=min((max(zv[temp_pp[jjj,]])-max(fitto[temp_pp[jjj,]]))/max(FF[temp_pp[jjj,]]))
            mi=min(mi)
            
            ##############################
            
            #mi=min((zv[pp_FF]-fitto[pp_FF])/FF[pp_FF])
            mi[mi<0]=0.01
            height_max[hh]=temp_height[hh]*mi*1.3    ########URINE 1.1
            height_min=height_max*0.4
            
            height[hh]=temp_height[hh]*mi
            
          }
        }
      }
    }
    fit=matrix(0,nrow=nnn,ncol=length(xv))
    area=rep(NA,nnn)
    
    for(i in 1:nnn)
      if(!is.na(ppm[i])){
        if(!peak.selection[i]){
          height[i]=NA
          height[i]=mean(height[link[i]==link],na.rm=TRUE)
        }
        fit[i,]=FUN[[i]](xv, height[i], ppm[[i]], fwhm[i], eta[i],J[i])
        area[i]=sum(fit[i,])
      }
    
    
    signal=rowSums(fit[,sel,drop=FALSE])/length(sel)    #################ERRORE
    
    colSums_fit=colSums(fit,na.rm = TRUE)
    
    acc=signal^2/SN
    
    
    error1=rep(NA,nnn)
    error2=rep(NA,nnn)
    error3=rep(NA,nnn)
    t_a=zv
    tt=colSums_fit<zv
    t_a[tt]=colSums_fit[tt]
    t_a[t_a<0]=0
    t_b=abs(zv-colSums_fit)
    a=zv-colSums_fit
    e1=a; e1[e1<0]=0
    e2=-a; e2[e2<0]=0
    for(i in 1:nnn){
      if(!all(is.na(ppm[[i]]))){
        
        
        
        
        error1[i]=sum(e1[sel_i[[i]]])/sum(fit[i,sel_i[[i]],drop=FALSE])
        error2[i]=sum(e2[sel_i[[i]]])/sum(fit[i,sel_i[[i]],drop=FALSE])
        error3[i]=sum( t_a[sel_i[[i]]])/sum( t_b[sel_i[[i]]])
        
      }
    }
    
    return(list(height=height,ppm=ppm,fwhm=fwhm,eta=eta,area=area,fit=fit,signal=signal/SN,acc=acc,error1=error1,error2=error2,error3=error3))
  }
