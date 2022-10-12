
vect=as.numeric(as.vector(rr[,"ensemble"]))
WHO=ave(1:length(vect), vect, FUN = seq_along)
WHO[WHO!=1]="NO"
WHO[WHO==1]="YES"
rr=cbind(rr,WHO)

colorandom=as.numeric(as.factor(rr[,"link"]))
colo=sample(rainbow(max(rr[,"link"]),alpha = 0.5))

super_ppm=matrix(ncol=nrow(rr),nrow=addr$n)
super_area=matrix(ncol=nrow(rr),nrow=addr$n)
super_fwhm=matrix(ncol=nrow(rr),nrow=addr$n)

colnames(super_ppm)=as.vector(rr[,"metabolite"])
rownames(super_ppm)=addr$Names
colnames(super_area)=as.vector(rr[,"metabolite"])
rownames(super_area)=addr$Names
colnames(super_fwhm)=as.vector(rr[,"metabolite"])
rownames(super_fwhm)=addr$Names
r.baseline=0
l.baseline=0


nmarks=8
total_marker=matrix(1,nrow=addr$n,ncol=9)

k=1
for(k in k:addr$n){
  print(k)
  ppm_ref=tsp[[k]]$ppm
  eta_ref=tsp[[k]]$eta   ###   
  eta_ref=median(eta_ref)
  fwhm_ref=tsp[[k]]$fwhm  #####  
  fwhm_ref=median(fwhm_ref) *0.7
  height_ref=tsp[[k]]$height
  
  
  x=spectra[[k]]$Spectra[[1]][[2]]
  y=spectra[[k]]$Spectra[[1]][[1]]
  z=spectra[[k]]$resolved[[1]]
  b=spectra[[k]]$baseline[[1]]
  
  
  
  sel=x>7.5 & x<11.00
  yy=y[sel]
  xx=x[sel]
  ll=predict(loess(yy~xx),xx)
  SN = sd(yy-ll)
  
  ppm_all=NULL
  LL=NULL
  RR=NULL
  right_vis=NULL
  left_vis=NULL
  for(im in 1:nrow(rr)){
    
    center=rr[im,"coeff"]
    right=center-rr[im,"width"]
    left=center+rr[im,"width"]
    left=min(left,rr[im,"left"])
    right=max(right,rr[im,"right"])
    
    right_vis[im]=right-0.03
    left_vis[im]=left+0.03
    J=rr[im,"J"]
    JJ=J
    JJ[is.na(JJ)]=0
    JJ=max(JJ)
    JJ=JJ/2
    
    RR[im]=right
    LL[im]=left
    error=rr[im,"error"]
    simila=rr[im,"simila"]
    
    sel=x>right & x<left
    if(sum(sel)>5){
      x1=x[sel]
      y1=y[sel]
      z1=z[sel]
      b1=b[sel]
      
      ppm_all[im]=do.call(as.vector(rr[im,"search"]), list(x=x,
                                                           y=z,
                                                           left=left,
                                                           right=right,
                                                           x1=x1,
                                                           y1=z1,
                                                           J=J,
                                                           error=error,
                                                           simila=simila,
                                                           fwhm_ref=tsp[[k]]$fwhm,
                                                           eta_ref=tsp[[k]]$eta))
    }else{
      ppm_all[im]=NA
    }
  }
  
  for(i in 1:length(ppm_all)){
    link=rr[,"link"]
    ws=which(link[i]==link)
    if(any(is.na(ppm_all[ws]))){
      ppm_all[ws]=NA
    }
  }
  super_ppm[k,]=ppm_all
  
  
  
  ss=rr[,"WHO"]=="YES" &   ave(1:nrow(rr),   rr[,"ensemble"], FUN = seq_along)==1
  ss=(1:nrow(rr))[ss]
  
  for(im in ss){
    
    
    metavis=as.matrix(which(as.vector(rr[im,"ensemble"])==as.vector(rr[,"ensemble"])))
    
    xxlim=matrix(nrow=1,ncol=2)
    yylim=matrix(nrow=1,ncol=2)
    sel=rep(FALSE,length(x))
    for(jk in 1:1){
      
      xxlim[jk,]=c(max(left_vis[metavis[,jk]]),min(right_vis[metavis[,jk]]))
      
      temp_selection=(x>xxlim[jk,2] & x<xxlim[jk,1])
      
      sel=sel | temp_selection
      
      ytemp=z[temp_selection]                                      
      yylim[jk,]=c(min(c(0,min(ytemp))),max(max(ytemp)*2,10000))
      
      
    }
    sel=x>xxlim[,2] & x<xxlim[,1]
    xv=x[sel]
    yv=y[sel]
    zv=z[sel]
    bv=0 
    
    ww=which(as.vector(rr[im,"ensemble"])==as.vector(rr[,"ensemble"]))
    
    J=rr[ww,"J"]
    JJ=J
    JJ[is.na(JJ)]=0
    JJ=max(JJ)
    JJ=JJ/2
    size=rr[ww,"size"]
    fwhm_max=rr[ww,"fwhm_max"]
    fwhm_min=rr[ww,"fwhm_min"]
    eta_max=rr[ww,"eta_max"]
    eta_min=rr[ww,"eta_min"]
    shift=rr[ww,"shift"]
    
    ppm=ppm_all[ww]
    size=rr[ww,"size"]
    link=rr[ww,"link"]
    protons=rr[ww,"protons"]
    peak.selection=rr[ww,"peak.selection"]
    
    if(!all(unlist(lapply(ppm,function(x) any(is.na(x)))))){
      FUN=list()
      for(ii in 1:length(ww)){
        v=as.vector(rr[ww[ii],"FUN"])
        eval(parse(text = v))
        FUN[[ii]]=ff
      }
      
      
      sel_ww=unlist((lapply(ppm,function(x) !any(is.na(x)))))
      
      ww=ww[sel_ww]
      
      J=J[sel_ww]
      fwhm_min=fwhm_min[sel_ww]
      fwhm_max=fwhm_max[sel_ww]
      eta_min=eta_min[sel_ww]
      eta_max=eta_max[sel_ww]
      FUN=FUN[sel_ww]
      ppm=ppm[sel_ww]
      size=size[sel_ww]
      link=link[sel_ww]
      protons=protons[sel_ww]
      peak.selection=peak.selection[sel_ww]
      shift=shift[sel_ww]
      
      pp=do.call(as.vector(rr[im,"fit"]), list(xv=xv,zv=zv,size=size,ppm=ppm,J=J,
                                               FUN=FUN,fwhm_min=fwhm_min,fwhm_max=fwhm_max,
                                               eta_min=eta_min,eta_max=eta_max,SN=SN,
                                               link=link,protons=protons,peak.selection=peak.selection,
                                               shift=shift))
      
      second.b=rep(0,length(xv))
      
      for(im_temp in   which(rr[,"metabolite"]=="Lactate")){
        nom=paste(Metabolites_address,k,"-",as.vector(rr[im_temp,"metabolite"]),".png")
        png(nom,width = 2400,height = 1000)
        
        par(cex=6,mfrow=c(1,1))
        for(jk in 1:1){
          temp_selection=(xv>rr[im_temp,"right"] & xv<rr[im_temp,"left"])   #nuovo
          
          ytemp=zv[temp_selection]
          yylim[jk,2]=max(max(ytemp)*2,10000)
          
          plot(xv,zv,type="l",xlim=xxlim[jk,],ylim=yylim[jk,],cex.lab=2,cex.axis=2,
               main=as.vector(rr[im_temp,"metabolite"]),lwd=2,
               xlab="Chemical shift (ppm)",ylab="Intensity (a.u.)")
          
          abline(v=rr[im_temp,"left_vis"],col=2,lty=1,lwd=2)
          abline(v=rr[im_temp,"right_vis"],col=2,lty=1,lwd=2)
          
          
          if(any(ww==im_temp)){
            abline(v=ppm[ww==im_temp]+size[ww==im_temp],col=4,lty=1,lwd=2)
            abline(v=ppm[ww==im_temp]-size[ww==im_temp],col=4,lty=1,lwd=2)
            
            abline(v=ppm[ww==im_temp],col=3,lty=2,lwd=2)
            sel.range=(pp$fit[ww==im_temp,]/max(pp$fit[ww==im_temp,]))>0.03
            polygon(c(xv[sel.range],rev(xv[sel.range])),
                    c(bv[sel.range]+pp$fit[ww==im_temp,sel.range]+second.b[sel.range],
                      rev(bv[sel.range]+second.b[sel.range])),border=4,lwd=3)
          }
          
          mtext(addr$Names[k],cex=2)
          drangex=as.numeric(dist(xxlim[jk,]))
          drangey=as.numeric(dist(yylim[jk,]))
          
          
          if(!all(unlist(lapply(ppm,function(x) all(is.na(x)))))){
            super_area[k,ww]=pp$area
            super_ppm[k,ww]=pp$ppm
            super_fwhm[k,ww]=pp$fwhm
            points(xv,bv+second.b,col=4,type="l",lty=2,lwd=2)
            points(xv,colSums(pp$fit)+bv+second.b,type="l",col="#3333ff88",lwd=3)
            abline(v=r.baseline,col=5,lwd=2)
            abline(v=l.baseline,col=5,lwd=2)
          }
          for(h in 1:nrow(pp$fit)){
            poly.y=c(bv+pp$fit[h,]+second.b,rev(bv+second.b))
            poly.y[poly.y>(yylim[jk,2]*1.05)]=(yylim[jk,2]*1.05)
            poly.x=c(xv,rev(xv))
            polygon(poly.x,poly.y,col=colo[colorandom[ww[h]]],border=2)
            if(length(which(im_temp==ww))>0 & ww[h]==im_temp)
              points(super_ppm[k,im_temp],yylim[jk,1],cex=3,pch=21,bg=2)
          }
        }
        
        
        
        
        
        dev.off()   
      }
    }
  }
}
setwd(origin)
colnames(super_area)=as.vector(rr[,"metabolite"])
rownames(super_area)=addr$Names
