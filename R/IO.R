

read_address =
  function(Address,nick=0){
    name=nick
    if(is.list(Address)){
      if(Address$Partition=="")   stop("You have to declare the partition (e.g., C:, D:)")
      DIR <- paste(Address$Partition,"/",sep="")
      SNames=NULL
      
      if(Address$Data_directory!=""){
        if(Address$Data_directory=="*"){
          if(nick==0){
            name=1
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_1=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_1,"/",sep="")
            
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_1
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_1,sep="")
            }
            
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          DIR_1=Address$Data_directory
          DIR=paste(DIR,Address$Data_directory,"/",sep="")
          if(any(name==1)){
            if(is.null(SNames)){
              SNames=DIR_1
            }else{
              SNames=paste(SNames,"_",DIR_1,sep="")
            }
          }   
        }
      }
      
      if(Address$User_name!=""){
        if(Address$User_name=="*"){
          if(nick==0){
            name=2
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_2=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_2,"/",sep="")
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_2
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_2,sep="")
            }
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          
          DIR_2=Address$User_name
          DIR=paste(DIR,Address$User_name,"/",sep="")
          if(any(name==2)){
            if(is.null(SNames)){
              SNames=DIR_2
            }else{
              SNames=paste(SNames,"_",DIR_2,sep="")
            }
          }   
        }
      }
      
      
      if(Address$Spectroscopy!=""){
        if(Address$Spectroscopy=="*"){
          if(nick==0){
            name=3
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_3=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_3,"/",sep="")
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_3
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_3,sep="")
            }
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          DIR_3=Address$Spectroscopy
          DIR=paste(DIR,Address$Spectroscopy,"/",sep="")
          if(any(name==3)){
            if(is.null(SNames)){
              SNames=DIR_3
            }else{
              SNames=paste(SNames,"_",DIR_3,sep="")
            }
          }   
        }
      }
      
      if(Address$Dataset_name!=""){
        if(Address$Dataset_name=="*"){
          if(nick==0){
            name=4
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_4=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_4,"/",sep="")
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_4
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_4,sep="")
            }
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          DIR_4=Address$Dataset_name
          DIR=paste(DIR,Address$Dataset_name,"/",sep="")
          if(any(name==4)){
            if(is.null(SNames)){
              SNames=DIR_4
            }else{
              SNames=paste(SNames,"_",DIR_4,sep="")
            }
          }   
        }
      }
      
      if(Address$Esperiment_number!=""){
        if(Address$Esperiment_number=="*"){
          if(nick==0){
            name=5
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_5=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_5,"/",sep="")
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_5
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_5,sep="")
            }
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          DIR_5=Address$Esperiment_number
          DIR=paste(DIR,Address$Esperiment_number,"/",sep="")
          if(any(name==5)){
            if(is.null(SNames)){
              SNames=DIR_5
            }else{
              SNames=paste(SNames,"_",DIR_5,sep="")
            }
          }   
        }
      }
      FID <- DIR
      if(Address$Pdata!=""){
        if(Address$Pdata=="*"){
          if(nick==0){
            name=6
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_6=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_6,"/",sep="")
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_6
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_6,sep="")
            }
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          DIR_6=Address$Pdata
          DIR=paste(DIR,Address$Pdata,"/",sep="")
          if(any(name==6)){
            if(is.null(SNames)){
              SNames=DIR_6
            }else{
              SNames=paste(SNames,"_",DIR_6,sep="")
            }
          }   
        }
      }
      
      if(Address$Processing_number!=""){
        if(Address$Processing_number=="*"){
          if(nick==0){
            name=7
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_7=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_7,"/",sep="")
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_7
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_7,sep="")
            }
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          DIR_7=Address$Processing_number
          DIR=paste(DIR,Address$Processing_number,"/",sep="")
          if(any(name==7)){
            if(is.null(SNames)){
              SNames=DIR_7
            }else{
              SNames=paste(SNames,"_",DIR_7,sep="")
            }
          }   
        }
        
      }
    }else{
      DIR=Address
      naddr=length(Address)
      str=strsplit(Address,"/")
      
      if(nick==0){
        SNames=Address
      }else{
        
        SName=NULL
        for(jj in 1:naddr){
          
          SName[jj]=paste(str[[jj]][nick+1],collapse = "_")
          
        }
      }
      ii=1
      while(ii<length(str[[1]]) & !file.exists(paste(paste(str[[1]][1:ii],collapse = "/"),"/fid",sep=""))){
        ii=ii+1
      }
      FID=NULL
      for(jj in 1:naddr){
        
        FID[jj]=paste(paste(str[[jj]][1:ii],collapse = "/"),"/",sep="")
      }
      FID
    }
    
    
    
    Experiments <- DIR
    
    
    sel=file.exists(paste(DIR,"1r",sep=""))
    deleted=Experiments[!sel]
    Experiments <- Experiments[sel]
    FID <- FID[sel]
    SNames = SNames[sel]
    aa=strsplit(Experiments,"/")
    loc=t(matrix(unlist(aa),ncol=length(aa)))
    Spec=NULL
    
    return(list(Names=SNames, n=length(Experiments),Source=Experiments,loc=loc,deleted=deleted,FID=FID))
    
  }

read_spectra = 
  function(Address,name=0){
    
    read=read_address(Address,name)
    SNames=read$Names
    Experiments=read$Source
    FIDs=read$FID
    loc=read$loc
    deleted=read$deleted
    rbnmr <- function (i){
      DataR = NULL
      setwd (FIDs[i])
      con_acqu=file("acqus")
      Acqus <- readLines(con_acqu, n = -1)
      DATE=as.numeric(gsub(".*[=]", "",Acqus[pmatch("##$DATE= ",Acqus)]))
      NS <- as.numeric(gsub(".*[=]", "",Acqus[pmatch("##$NS=",Acqus)]))
      
      class(DATE) = c('POSIXt','POSIXct')
      
      setwd (Experiments[i])
      con_proc=file("procs")
      Procs <- readLines(con_proc, n = -1)
      BYTORD <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$BYTORDP=",Procs)]))
      
      if (BYTORD=="0") {ENDIAN="little"} else {ENDIAN="big"}
      NC_proc <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$NC_proc=",Procs)]))
      OFFSET <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$OFFSET=",Procs)]))
      SW_p <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$SW_p=",Procs)]))
      SF <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$SF=",Procs)]))
      SI <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$SI=",Procs)]))
      con = file("1r", open="rb")
      
      DataR$y <- readBin (con,what="int",n=SI,endian=ENDIAN)
      close (con)
      close(con_proc)
      close(con_acqu)
      DataR$x <- seq (OFFSET,OFFSET-SW_p/SF,length=SI)
      DataR$y <- DataR$y/(2^(-NC_proc))
      DataR$date=DATE
      DataR$ns=NS
      DataR
    }
    Spec <- lapply (1:length(Experiments),rbnmr)
    return(list(Spectra=Spec, 
                Names=SNames, 
                n=length(Experiments),
                Source=Experiments,
                loc=loc,
                deleted=deleted))
    
  }


baselinefun =
  function(spectra){
    
    spectra$baseline=list()
    spectra$resolved=list()
    
    k=1
    
    for(k in k:spectra$n){
      
      u=spectra$Spectra[[k]]
      x=u[[2]]
      y=u[[1]]
      
      b=rep(0,length(y))
      
      sel=matrix(nrow=200,ncol=2)  
      
      i=1;    sel[i,]=  c(8.6,8.58)
      i=i+1;  sel[i,]=  c(8.534,8.51)
      i=i+1;  sel[i,]=  c(8.49,8.46)
      i=i+1;  sel[i,]=  c(8.41,8.39)
      i=i+1;  sel[i,]=  c(8.33,8.25)
      i=i+1;  sel[i,]=  c(8.184,8.13)
      i=i+1;  sel[i,]=  c(7.95,7.80)
      i=i+1;  sel[i,]=  c(7.735,7.7)
      i=i+1;  sel[i,]=  c(7.23,7.18)
      i=i+1;  sel[i,]=  c(6.97,6.92)
      i=i+1;  sel[i,]=  c(6.55,6.50)
      i=i+1;  sel[i,]=  c(6.20,5.85)
      i=i+1;  sel[i,]=  c(5.4,5.21)
      i=i+1;  sel[i,]=  c(4.7013,4.673)
      i=i+1;  sel[i,]=  c(4.6,4.45)
      i=i+1;  sel[i,]=  c(4.44,4.30)
      i=i+1;  sel[i,]=  c(4.289,4.242)
      i=i+1;  sel[i,]=  c(4.235,4.10)
      i=i+1;  sel[i,]=  c(4.093,3.851)
      i=i+1;  sel[i,]=  c(3.85,3.838)
      i=i+1;  sel[i,]=  c(3.837,3.7355)
      i=i+1;  sel[i,]=  c(3.734,3.664)
      i=i+1;  sel[i,]=  c(3.663,3.648)
      i=i+1;  sel[i,]=  c(3.646,3.516)
      i=i+1;  sel[i,]=  c(3.508,3.17)
      i=i+1;  sel[i,]=  c(3.10,2.92)
      i=i+1;  sel[i,]=  c(2.85,2.63)
      i=i+1;  sel[i,]=  c(2.55,2.28)
      i=i+1;  sel[i,]=  c(2.22,1.85)
      i=i+1;  sel[i,]=  c(1.82,1.6)
      i=i+1;  sel[i,]=  c(1.517,1.447)
      i=i+1;  sel[i,]=  c(1.41,1.32)  #      i=i+1;  sel[i,]=  c(1.41,1.245)
      i=i+1;  sel[i,]=  c(1.30,1.245)
      i=i+1;  sel[i,]=  c(1.24,1.152)
      i=i+1;  sel[i,]=  c(1.085,0.7)
      i=i+1;  sel[i,]=  c(0.30,-0.250)
      
      sel=sel[1:i,,drop=FALSE]
      
      sel.v=rep(FALSE,length(b))
      
      
      ff=y
      ff2=ff
      for(i in 1:nrow(sel)){
        ss1=(x<sel[i,1] & x>sel[i,2])
        sel.v=sel.v | ss1
        ss2=(x<sel[i,1]+0.001 & x>sel[i,2]-0.001)
        xl=lm(y[which(ss2 & !ss1)]~x[which(ss2 & !ss1)])
        ff[ss1]=x[ss1]* xl$coefficients[2]+ xl$coefficients[1]
        
        
      }  
      
      
      
      sel=(x<4.705)
      
      b[sel]=predict(loess(ff[sel]~x[sel],span=0.05),x[sel])
      
      sel=(x>5.20)
      b[sel]=predict(loess(ff[sel]~x[sel],span=0.03),x[sel])
      
      
      sel=(x>4.705 & x <5.20) 
      b[sel]=y[sel]
      sel=(x>9) 
      b[sel]=y[sel]
      sel=(x <(-0.2)) 
      b[sel]=y[sel]
      
      ###################################
      
      z=y-b
      spectra$baseline[[k]]=b
      spectra$resolved[[k]]=z
    }
    spectra
  }


