



eta_ref=rep(NA,addr$n)
fwhm_ref=rep(NA,addr$n)
area_ref=rep(NA,addr$n)
ppm_ref=rep(NA,addr$n)
RMSD_ref=rep(NA,addr$n)
Ans_ref=rep(NA,addr$n)
ns_ref=rep(NA,addr$n)
x_ref=1:addr$n
tsp=list()
for(k in 1:addr$n){
  
  u=spectra[[k]]$Spectra[[1]]
  x=u[[2]]
  y=u[[1]]
  z=spectra[[k]]$resolved[[1]]
  
  tsp[[k]]=TSP_ref.500(x,z)
  eta_ref[k]=tsp[[k]]$eta
  fwhm_ref[k]=tsp[[k]]$fwhm
  area_ref[k]=tsp[[k]]$area
  ppm_ref[k]=tsp[[k]]$ppm
  RMSD_ref[k]=tsp[[k]]$RMSD/spectra[[k]]$Spectra[[1]]$ns
  Ans_ref[k]=tsp[[k]]$area/spectra[[k]]$Spectra[[1]]$ns
  ns_ref[k]=spectra[[k]]$Spectra[[1]]$ns
}



for(k in 1:addr$n){
  u=spectra[[k]]$Spectra[[1]]
  x=u[[2]]
  y=u[[1]]
  z=spectra[[k]]$resolved[[1]]
  
  
  nom=paste(TSP_address,k,"-TSP.png",sep="")
  
  sel=x>(-0.03) & x<0.03
  x1=x[sel]
  y1=z[sel]
  
  x2=seq(-0.03,0.03,length.out = 3001)
  fit=TSP.500(x2,
              tsp[[k]]$height,
              tsp[[k]]$ppm,
              tsp[[k]]$fwhm,
              tsp[[k]]$eta)
  
  png(nom,width = 2400,height = 1000)
  
  layout(matrix(c(1,1,2,3,1,1,4,5,1,1,6,7), 3, 4, byrow = TRUE),widths=c(2,1,1), heights=c(1,1,1))
  par(cex=1.5,mai=c(1.5,1.5,0.5,0.5))
  plot(x1,y1,xlim=c(0.03,-0.03),xlab="Chemical shift (ppm)",ylab="Intensity (a.u.)",main=addr$Names[k])
  points(x2,fit,type="l",col=2,cex=2)
  abline(v=tsp[[k]]$ppm+c(-0.002,0.002),col=2,cex=2)
  ############################################
  plot(fwhm_ref,type="n",main="FWHM",xlab="")
  points(x_ref[-k],fwhm_ref[-k],cex=2)
  points(x_ref[k],fwhm_ref[k],pch=21,bg=2,cex=2)
  ############################################
  plot(eta_ref,type="n",main="ETA",xlab="")
  points(x_ref[-k],eta_ref[-k],cex=2)
  points(x_ref[k],eta_ref[k],pch=21,bg=2,cex=2)
  ############################################
  plot(ppm_ref,type="n",main="ppm",xlab="")
  points(x_ref[-k],ppm_ref[-k],cex=2)
  points(x_ref[k],ppm_ref[k],pch=21,bg=2,cex=2)
  ############################################
  plot(RMSD_ref,type="n",main="RMSD",xlab="")
  points(x_ref[-k],RMSD_ref[-k],cex=2)
  points(x_ref[k],RMSD_ref[k],pch=21,bg=2,cex=2)
  ############################################
  plot(Ans_ref,type="n",main="Area/ns",xlab="")
  points(x_ref[-k],Ans_ref[-k],cex=2)
  points(x_ref[k],Ans_ref[k],pch=21,bg=2,cex=2)
  ############################################
  plot(ns_ref,type="n",main="number of scans",xlab="")
  points(x_ref[-k],ns_ref[-k],cex=2)
  points(x_ref[k],ns_ref[k],pch=21,bg=2,cex=2)
  
  dev.off()
  
  
  print(k)
  
}

