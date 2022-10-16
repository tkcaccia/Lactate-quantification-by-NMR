# Lactate quantification by NMR in prostate tissue extracts

![This is an image](https://github.com/tkcaccia/Lactate-quantification-by-NMR/blob/main/Figures/Lactate.png)

This is a companion repository to the manuscript “*High-fat diet and MYC cooperation promotes lactate accumulation and tumour microenvironment remodelling in prostate cancer*” by Boufaied N, Chetta P, Hallal T *et al.* The code has been generated to calculate lactate concentration in murine prostate tissues from NMR spectra. The entire code is in R and does not require non-standard hardware or installation. The raw lactate NMR spectra were acquired using a Bruker AVANCE III 500 spectrometer equipped with a 5 mm triple resonance 1H/13C/15N TXI probe and x, y, z gradient coils. Detailed description of NMR experiment and analysis is provided in the main manuscript and supplementary information. Lactate NMR spectra are provided in the folder Data. The expected run time to reproduce data shown in the manuscript is about 5 min.

## Initialization
After downloading the repository, please follow the script below. Both functions and parameters for the fitting are uploaded in the R environment. Do not forget to modify the address of the Working_folder and Output_address as the following:

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
```


## Data loading
The lactate NMR spectra (from the folder Data) are uploaded in the R environment and ordered according to the **Table S9** provided in the Supplementary Material:

```
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
```


## Baseline correction
Spectra baseline is then corrected:

```
i=1
spectra=list()
for(i in i:addr$n){
  spectrum=addr$Source[i]
  spectra[[i]]=read_spectra(spectrum)
  spectra[[i]]=baselinefun(spectra[[i]])
  print(i)
}
```

## Fitting of TSP and lactate signals
The signals of sodium trimethylsilyl propionate (TSP) and lactate are quantified and the output of the fitting is saved in the folder defined by the variable `Output_address` for manual inspection. The signal of lactate is deconvoluted from an unidentified metabolite.

```
setwd(Working_folder)
source("R/TSP.R")
source("R/quantification.R")
lactate=super_area[,"Lactate"]
TSP=area_ref
```

## Lactate absolute quantification
The goal is to obtain a final absolute lactate concentration in nmol/mg prostate tissue. First you need to divide each signal for the number of protons in each signal and multiply the ratio of lactate/TSP signal for TSP concentration (5.805 mM). Since only a portion of the entire polar phase has been analyzed (here called polar phase analyzed, see methods), the final lactate concentration is then calculated for the entire polar phase (obtained from the entire prostate tissue extraction, see method) and this is then divided for the weight of the prostate tissue. The amount of polar phase analyzed for each samples is saved in the variable `polar_phase_analyzed`. The amount of polar phase for each samples is saved in the variable `Volume_polar_phase`.

```
lactate=lactate/3
TSP=TSP/9
lactate=5.804841238*lactate/TSP

polar_phase_analyzed=c(735,565,400,380,430,606,545,275)
Volume_polar_phase=c(1074.1,743.2,575.4,520.3,587.4,791.2,656.9,469.9)
DLP_weight=c(44.8,31,24,21.7,24.5,33,27.4,19.6)

lactate=(lactate*0.52*Volume_polar_phase/polar_phase_analyzed)/DLP_weight
lactate=lactate*1000
```
## Statistical analysis
To determine whether a significant difference in lactate concentration is observed between prostate tissues from high-fat diet fed/obese mice and control diet fed mice, a t-test was performed after default confirmation of normality data distribution and equal variance.

### Shapiro-Wilk Normality Test
```
shapiro.test((lactate))
```
### F Test to Compare Two Variances
```
labels=c("MYC_CTD","MYC_CTD","MYC_CTD","MYC_CTD","MYC_CTD","MYC_HFD","MYC_HFD","MYC_HFD")
var.test((lactate)[labels=="MYC_CTD"],(lactate)[labels=="MYC_HFD"])
```
### Student's t-Test where the pooled variance is used to estimate the variance 
```
t.test((lactate)[labels=="MYC_CTD"],(lactate)[labels=="MYC_HFD"], var.equal=TRUE)
```
