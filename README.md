# Lactate quantification by NMR in prostate tissue extracts

![This is an image](https://github.com/tkcaccia/Figures/Lactate.png)

This is a companion repository to the paper "Cooperation of high-fat diet and MYC oncogenic drive promotes lactate accumulation and immunosuppression to accelerate prostate cancer progression". All code is in R. The raw NMR spectra acquired using a Bruker AVANCE III 500 spectrometr equipped with a 5 mm inverse triple resonance 1H/13C/15N TXI probe and x, y, z gradient coils are provvided in the folder Data.



```
Working_folder="./Lactate-quantification-by-NMR/"
Output_address="./Lactate-quantification-by-NMR/Output"

setwd(Working_folder)

source("R/fitting.R")
source("R/deconvolution.R")
source("R/IO.R")

load("Data/Fitting-parameters.RData")

dir.create(Output_address)
TSP_address=paste(Output_address, "TSP/",sep="/")
Metabolites_address=paste(Output_address, "Metabolites/",sep="/")
dir.create(TSP_address)
dir.create(Metabolites_address)


##################################################################

Address=list(Partition=Working_folder,
             Data_directory="Data",
             User_name="",
             Spectroscopy="*",
             Dataset_name="*",
             Esperiment_number="",
             Pdata="pdata",
             Processing_number="1")

addr=read_address(Address)


oo=1:addr$n
names(oo)=addr$Names
oo=oo[c("prostate extracts_54","prostate extracts_104","prostate extracts_24",
        "prostate extracts_134","prostate extracts_184",
        "prostate extracts_34","prostate extracts_124","prostate extracts_194")]

addr$Names=addr$Names[oo]
addr$Source=addr$Source[oo]
addr$loc=addr$loc[oo,]
addr$FID=addr$FID[oo]
addr$n=length(oo)




i=1
spectra=list()
for(i in i:addr$n){
  spectrum=addr$Source[i]
  spectra[[i]]=read_spectra(spectrum)
  spectra[[i]]=baselinefun(spectra[[i]])
  print(i)
}





##########################################################################


setwd(Working_folder)

source("R/TSP.R")

source("R/quantification.R")

lactate=super_area[,"Lactate"]
TSP=area_ref

# Divide by the number of protons in each signal
lactate=lactate/3
TSP=TSP/9

# TSP's concentration is 5.804841238 mM

lactate=5.804841238*lactate/TSP

Volume_selected=c(735,565,400,380,430,606,545,275)
Volume_polar_phase=c(1074.1,743.2,575.4,520.3,587.4,791.2,656.9,469.9)
DLP_weight=c(44.8,31,24,21.7,24.5,33,27.4,19.6)

labels=c("MYC_CTD","MYC_CTD","MYC_CTD","MYC_CTD","MYC_CTD","MYC_HFD","MYC_HFD","MYC_HFD")


lactate=(lactate*0.52*Volume_polar_phase/Volume_selected)/DLP_weight

lactate=lactate*1000


shapiro.test((lactate))

var.test((lactate)[labels=="MYC_CTD"],(lactate)[labels=="MYC_HFD"])
t.test((lactate)[labels=="MYC_CTD"],(lactate)[labels=="MYC_HFD"], var.equal=TRUE)
```
