Compute per-barcode ACE2 binding affinity
================
Tyler Starr
6/9/2022

This notebook reads in per-barcode counts from `count_variants.ipynb`
for ACE2-binding Tite-seq experiments, computes functional scores for
RBD ACE2-binding affiniity, and does some basic QC on variant binding
functional scores.

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))

#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$Titeseq_Kds_dir)){
  dir.create(file.path(config$Titeseq_Kds_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /app/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas_haswellp-r0.3.7.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3     forcats_0.4.0     stringr_1.4.0     dplyr_0.8.3      
    ##  [5] purrr_0.3.3       readr_1.3.1       tidyr_1.0.0       tibble_3.0.2     
    ##  [9] ggplot2_3.3.0     tidyverse_1.3.0   data.table_1.12.8 yaml_2.2.0       
    ## [13] knitr_1.26       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.3       cellranger_1.1.0 pillar_1.4.5     compiler_3.6.2  
    ##  [5] dbplyr_1.4.2     tools_3.6.2      digest_0.6.23    lubridate_1.7.4 
    ##  [9] jsonlite_1.6     evaluate_0.14    lifecycle_0.2.0  gtable_0.3.0    
    ## [13] pkgconfig_2.0.3  rlang_0.4.7      reprex_0.3.0     cli_2.0.0       
    ## [17] rstudioapi_0.10  DBI_1.1.0        haven_2.2.0      xfun_0.11       
    ## [21] withr_2.1.2      xml2_1.2.2       httr_1.4.1       fs_1.3.1        
    ## [25] hms_0.5.2        generics_0.0.2   vctrs_0.3.1      grid_3.6.2      
    ## [29] tidyselect_1.1.0 glue_1.3.1       R6_2.4.1         fansi_0.4.0     
    ## [33] readxl_1.3.1     rmarkdown_2.0    modelr_0.1.5     magrittr_1.5    
    ## [37] backports_1.1.5  scales_1.1.0     ellipsis_0.3.0   htmltools_0.4.0 
    ## [41] rvest_0.3.5      assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.3   
    ## [45] munsell_0.5.0    broom_0.7.0      crayon_1.3.4

## Setup

First, we will read in metadata on our sort samples, the table giving
number of reads of each barcode in each of the sort bins, and the
barcode-variant lookup tables, and merge these tables together.

``` r
#read dataframe with list of barcode runs
barcode_runs <- read.csv(file=config$barcode_runs,stringsAsFactors=F); barcode_runs <- subset(barcode_runs, select=-c(R1))

#eliminate rows from barcode_runs that are not from an tite-seq experiment
barcode_runs <- barcode_runs[barcode_runs$sample_type == "TiteSeq",]

#read file giving count of each barcode in each sort partition
counts <- data.table(read.csv(file=config$variant_counts_file,stringsAsFactors=F))

#eliminate rows from counts that are not part of an titration bin sample
counts <- subset(counts, sample %in% barcode_runs[barcode_runs$sample_type=="TiteSeq","sample"])

#read in barcode-variant lookup tables
dt_Wh1 <- data.table(read.csv(file=config$codon_variant_table_file_Wuhan_Hu_1,stringsAsFactors=F))
dt_BA1 <- data.table(read.csv(file=config$codon_variant_table_file_BA1,stringsAsFactors=F))
dt_BA2 <- data.table(read.csv(file=config$codon_variant_table_file_BA2,stringsAsFactors=F))

#merge, eliminate barcodes duplicated within a library
dt <- rbind(dt_Wh1,dt_BA1,dt_BA2); setkey(dt,barcode,library)
duplicates <- dt[duplicated(dt,by=c("barcode","library")),.(library,barcode)] #the data.table duplciates function annoyingly only flags the first of each duplicate so doesn't intrinsically allow removal of both of the entries of the duplicate. So, flat what are duplciates, and then remove
dt[,duplicate:=FALSE]
for(i in 1:nrow(duplicates)){
  dt[library==duplicates[i,library] & barcode==duplicates[i,barcode],duplicate:=TRUE]
}
dt <- dt[duplicate==FALSE,]; dt[,duplicate:=NULL]

dt <- merge(counts, dt, by=c("library","barcode")); rm(dt_Wh1);rm(dt_BA1);rm(dt_BA2);rm(counts); rm(duplicates)

#append replicate identifier to the library name
dt[,library:=paste(library,replicate,sep=""),by=c("library","replicate")];dt[,replicate:=NULL]

#make tables giving names of Titeseq samples and the corresponding ACE2 incubation concentrations
samples_TiteSeq <- data.frame(sample=unique(paste(barcode_runs[barcode_runs$sample_type=="TiteSeq","sample_type"],formatC(barcode_runs[barcode_runs$sample_type=="TiteSeq","concentration"], width=2,flag="0"),sep="_")),conc=c(10^-6, 10^-7, 10^-8, 10^-9, 10^-10, 10^-11, 10^-12, 10^-13,0))
```

Convert from Illumina read counts to estimates of the number of cells
that were sorted into a bin, and add some other useful information to
our data frame.

``` r
#for each bin, normalize the read counts to the observed ratio of cell recovery among bins
for(i in 1:nrow(barcode_runs)){
  lib <- paste(as.character(barcode_runs$library[i]), as.character(barcode_runs$replicate[i]), sep="")
  bin <- as.character(barcode_runs$sample[i])
  ratio <- sum(dt[library==lib & sample==bin,"count"])/barcode_runs$number_cells[i]
  if(ratio<1){ #if there are fewer reads from a FACS bin than cells sorted
    dt[library==lib & sample==bin, count.norm := as.numeric(count)] #don't normalize cell counts, make count.norm the same as count
    print(paste("reads < cells for",lib,bin,", un-normalized (ratio",ratio,")")) #print to console to inform of undersampled bins
  }else{
    dt[library==lib & sample==bin, count.norm := as.numeric(count/ratio)] #normalize read counts by the average read:cell ratio, report in new "count.norm" column
    print(paste("read:cell ratio for",lib,bin,"is",ratio))
  }
}
```

    ## [1] "read:cell ratio for pool1A TiteSeq_01_bin1 is 1.11050274280336"
    ## [1] "read:cell ratio for pool1A TiteSeq_01_bin2 is 1.00813623099594"
    ## [1] "read:cell ratio for pool1A TiteSeq_01_bin3 is 1.37137149083914"
    ## [1] "read:cell ratio for pool1A TiteSeq_01_bin4 is 1.39169899847368"
    ## [1] "reads < cells for pool1A TiteSeq_02_bin1 , un-normalized (ratio 0.745538616468764 )"
    ## [1] "read:cell ratio for pool1A TiteSeq_02_bin2 is 1.23872617810293"
    ## [1] "read:cell ratio for pool1A TiteSeq_02_bin3 is 1.36868687136275"
    ## [1] "read:cell ratio for pool1A TiteSeq_02_bin4 is 1.51512712834341"
    ## [1] "read:cell ratio for pool1A TiteSeq_03_bin1 is 1.59207041452787"
    ## [1] "read:cell ratio for pool1A TiteSeq_03_bin2 is 1.44499939364318"
    ## [1] "read:cell ratio for pool1A TiteSeq_03_bin3 is 1.42740048680076"
    ## [1] "read:cell ratio for pool1A TiteSeq_03_bin4 is 1.35930550791485"
    ## [1] "read:cell ratio for pool1A TiteSeq_04_bin1 is 1.38180223356423"
    ## [1] "read:cell ratio for pool1A TiteSeq_04_bin2 is 1.39084816487766"
    ## [1] "read:cell ratio for pool1A TiteSeq_04_bin3 is 1.35712388128491"
    ## [1] "reads < cells for pool1A TiteSeq_04_bin4 , un-normalized (ratio 0.70288784970055 )"
    ## [1] "read:cell ratio for pool1A TiteSeq_05_bin1 is 1.32287021605954"
    ## [1] "read:cell ratio for pool1A TiteSeq_05_bin2 is 1.52277079597675"
    ## [1] "reads < cells for pool1A TiteSeq_05_bin3 , un-normalized (ratio 0.753941749698894 )"
    ## [1] "read:cell ratio for pool1A TiteSeq_05_bin4 is 1.31467854376452"
    ## [1] "read:cell ratio for pool1A TiteSeq_06_bin1 is 1.4332920250161"
    ## [1] "read:cell ratio for pool1A TiteSeq_06_bin2 is 1.1514268551807"
    ## [1] "reads < cells for pool1A TiteSeq_06_bin3 , un-normalized (ratio 0.330102315562736 )"
    ## [1] "read:cell ratio for pool1A TiteSeq_06_bin4 is 6.55014749262537"
    ## [1] "read:cell ratio for pool1A TiteSeq_07_bin1 is 1.25710813241139"
    ## [1] "read:cell ratio for pool1A TiteSeq_07_bin2 is 1.45975074698315"
    ## [1] "read:cell ratio for pool1A TiteSeq_07_bin3 is 1.20073891625616"
    ## [1] "reads < cells for pool1A TiteSeq_07_bin4 , un-normalized (ratio 0.352678571428571 )"
    ## [1] "read:cell ratio for pool1A TiteSeq_08_bin1 is 1.40013384587533"
    ## [1] "read:cell ratio for pool1A TiteSeq_08_bin2 is 1.3681367607412"
    ## [1] "reads < cells for pool1A TiteSeq_08_bin3 , un-normalized (ratio 0.310344827586207 )"
    ## [1] "reads < cells for pool1A TiteSeq_08_bin4 , un-normalized (ratio 0.295454545454545 )"
    ## [1] "read:cell ratio for pool1A TiteSeq_09_bin1 is 1.48183802015763"
    ## [1] "read:cell ratio for pool1A TiteSeq_09_bin2 is 1.22922186578234"
    ## [1] "reads < cells for pool1A TiteSeq_09_bin3 , un-normalized (ratio 0.306188925081433 )"
    ## [1] "reads < cells for pool1A TiteSeq_09_bin4 , un-normalized (ratio 0.213576158940397 )"
    ## [1] "read:cell ratio for pool2A TiteSeq_01_bin1 is 1.44191696301595"
    ## [1] "read:cell ratio for pool2A TiteSeq_01_bin2 is 1.20662507587724"
    ## [1] "read:cell ratio for pool2A TiteSeq_01_bin3 is 1.13604669765249"
    ## [1] "read:cell ratio for pool2A TiteSeq_01_bin4 is 1.33645409368899"
    ## [1] "read:cell ratio for pool2A TiteSeq_02_bin1 is 1.78922969650961"
    ## [1] "read:cell ratio for pool2A TiteSeq_02_bin2 is 1.21751633174173"
    ## [1] "read:cell ratio for pool2A TiteSeq_02_bin3 is 1.1713037521961"
    ## [1] "read:cell ratio for pool2A TiteSeq_02_bin4 is 1.44315531118955"
    ## [1] "read:cell ratio for pool2A TiteSeq_03_bin1 is 1.38421333159703"
    ## [1] "read:cell ratio for pool2A TiteSeq_03_bin2 is 1.28520231923036"
    ## [1] "read:cell ratio for pool2A TiteSeq_03_bin3 is 1.32895582903287"
    ## [1] "read:cell ratio for pool2A TiteSeq_03_bin4 is 1.46519993532319"
    ## [1] "read:cell ratio for pool2A TiteSeq_04_bin1 is 1.49516568843486"
    ## [1] "read:cell ratio for pool2A TiteSeq_04_bin2 is 1.38664772359428"
    ## [1] "read:cell ratio for pool2A TiteSeq_04_bin3 is 1.33152354674085"
    ## [1] "reads < cells for pool2A TiteSeq_04_bin4 , un-normalized (ratio 0.957967629131018 )"
    ## [1] "read:cell ratio for pool2A TiteSeq_05_bin1 is 1.27548147781028"
    ## [1] "read:cell ratio for pool2A TiteSeq_05_bin2 is 1.41979815124894"
    ## [1] "read:cell ratio for pool2A TiteSeq_05_bin3 is 14.3515067544163"
    ## [1] "reads < cells for pool2A TiteSeq_05_bin4 , un-normalized (ratio 0.0840840840840841 )"
    ## [1] "read:cell ratio for pool2A TiteSeq_06_bin1 is 1.36289253025947"
    ## [1] "read:cell ratio for pool2A TiteSeq_06_bin2 is 1.28604154527947"
    ## [1] "reads < cells for pool2A TiteSeq_06_bin3 , un-normalized (ratio 0.395541401273885 )"
    ## [1] "reads < cells for pool2A TiteSeq_06_bin4 , un-normalized (ratio 0.0534351145038168 )"
    ## [1] "read:cell ratio for pool2A TiteSeq_07_bin1 is 1.53485062988225"
    ## [1] "read:cell ratio for pool2A TiteSeq_07_bin2 is 1.27425615674748"
    ## [1] "reads < cells for pool2A TiteSeq_07_bin3 , un-normalized (ratio 0.13588850174216 )"
    ## [1] "reads < cells for pool2A TiteSeq_07_bin4 , un-normalized (ratio 0.0765027322404372 )"
    ## [1] "read:cell ratio for pool2A TiteSeq_08_bin1 is 1.51899972585826"
    ## [1] "read:cell ratio for pool2A TiteSeq_08_bin2 is 1.57944008077596"
    ## [1] "reads < cells for pool2A TiteSeq_08_bin3 , un-normalized (ratio 0.21256038647343 )"
    ## [1] "reads < cells for pool2A TiteSeq_08_bin4 , un-normalized (ratio 0.0957178841309824 )"
    ## [1] "read:cell ratio for pool2A TiteSeq_09_bin1 is 1.43692679360111"
    ## [1] "read:cell ratio for pool2A TiteSeq_09_bin2 is 1.12152192016726"
    ## [1] "reads < cells for pool2A TiteSeq_09_bin3 , un-normalized (ratio 0.335949764521193 )"
    ## [1] "reads < cells for pool2A TiteSeq_09_bin4 , un-normalized (ratio 0.0263157894736842 )"

``` r
#annotate each barcode as to whether it's a homolog variant, SARS-CoV-2 wildtype, synonymous muts only, stop, nonsynonymous, >1 nonsynonymous mutations
dt[,variant_class:=as.character(NA)]
dt[n_codon_substitutions==0, variant_class := "wildtype"]
dt[n_codon_substitutions > 0 & n_aa_substitutions==0, variant_class := "synonymous"]
dt[n_aa_substitutions>0 & grepl("*",aa_substitutions,fixed=T), variant_class := "stop"]
dt[n_aa_substitutions == 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := "1 nonsynonymous"]
dt[n_aa_substitutions > 1 & !grepl("*",aa_substitutions,fixed=T), variant_class := ">1 nonsynonymous"]

#cast the data frame into wide format
dt <- dcast(dt, library + barcode + target + variant_class + aa_substitutions + n_aa_substitutions ~ sample, value.var="count.norm")
```

## Calculating mean bin for each barcode at each sample concentration

Next, for each barcode at each of the ACE2 concentrations, calculate the
“mean bin” response variable. This is calculated as a simple mean, where
the value of each bin is the integer value of the bin (bin1=unbound,
bin4=highly bound) – because of how bins are defined, the mean
fluorescence of cells in each bin are equally spaced on a log-normal
scale, so mean bin correlates with simple mean fluorescence.

We do not use the fluorescence boundaries of the FACS bins in our
calculations here, but we provide them for posterity’s sake below. For
the library titration sorts, the fluorescence boundaries for bins 1-4
are as follows:

    (-288, 254), (255, 1633), (1634, 10348), (10349, 262143)

``` r
#function that returns mean bin and sum of counts for four bins cell counts. Includes cutoffs for bimodal sample splits to filter out
calc.meanbin <- function(vec, split13filter=0.4, split24filter=0.4, split14filter=0.4){
  total <- sum(vec)
  if(is.na(total) | (vec[1] > split13filter*total & vec[3] > split13filter*total) | (vec[2] > split24filter*total & vec[4] > split24filter*total) | (vec[1] > split14filter*total & vec[4] > split14filter*total)){
    return(list(NA,NA))
  }else{
    return( list((vec[1]*1+vec[2]*2+vec[3]*3+vec[4]*4)/(vec[1]+vec[2]+vec[3]+vec[4]), total) )
  }
}
  

#iterate through Titeseq samples, compute mean_bin and total_count for each barcode variant
for(i in 1:nrow(samples_TiteSeq)){ #iterate through titeseq sample (concentration)
  meanbin_out <- paste(samples_TiteSeq[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
  totalcount_out <- paste(samples_TiteSeq[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
  bin1_in <- paste(samples_TiteSeq[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
  bin2_in <- paste(samples_TiteSeq[i,"sample"],"_bin2",sep="")
  bin3_in <- paste(samples_TiteSeq[i,"sample"],"_bin3",sep="")
  bin4_in <- paste(samples_TiteSeq[i,"sample"],"_bin4",sep="")
  dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
}
```

## Fit titration curves

We will use nonlinear least squares regression to fit curves to each
barcode’s titration series. We will do weighted nls, using the empirical
variance estimates from above to weight each observation. We will also
include a minimum cell count that is required for a meanbin estimate to
be used in the titration fit, and a minimum number of concentrations
with determined meanbin that is required for a titration to be reported.

``` r
#For QC and filtering, output columns giving the average number of cells that were sampled for a barcode across the 9 sample concentrations, and a value for the number of meanbin estimates that were removed for being below the # of cells cutoff
cutoff <- 2
dt[,TiteSeq_avgcount := mean(c(TiteSeq_01_totalcount,TiteSeq_02_totalcount,TiteSeq_03_totalcount,TiteSeq_04_totalcount,
                                TiteSeq_05_totalcount,TiteSeq_06_totalcount,TiteSeq_07_totalcount,TiteSeq_08_totalcount,
                                TiteSeq_09_totalcount),na.rm=T),by=c("library","barcode")]

#number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
dt[,TiteSeq_min_cell_filtered := sum(c(c(TiteSeq_01_totalcount,TiteSeq_02_totalcount,TiteSeq_03_totalcount,TiteSeq_04_totalcount,
                                        TiteSeq_05_totalcount,TiteSeq_06_totalcount,TiteSeq_07_totalcount,TiteSeq_08_totalcount,
                                        TiteSeq_09_totalcount)<cutoff,is.na(c(TiteSeq_01_totalcount,TiteSeq_02_totalcount,TiteSeq_03_totalcount,TiteSeq_04_totalcount,
                                                                             TiteSeq_05_totalcount,TiteSeq_06_totalcount,TiteSeq_07_totalcount,TiteSeq_08_totalcount,
                                                                             TiteSeq_09_totalcount))),na.rm=T),by=c("library","barcode")]

#function that fits a nls regression to the titration series, including an option to filter below certain thresholds for average cells across all samples, and number of samples below a cutoff of cells
fit.titration <- function(y.vals,x.vals,count.vals,min.cfu=cutoff,
                          min.means=0.8,min.average=cutoff,Kd.start=1e-9,
                          a.start=3,a.lower=2,a.upper=3,
                          b.start=1,b.lower=1,b.upper=1.5){
  indices <- count.vals>min.cfu & !is.na(y.vals)
  y <- y.vals[indices]
  x <- x.vals[indices]
  if((length(y) < min.means*length(y.vals)) | (mean(count.vals,na.rm=T) < min.average)){ #return NAs if < min.means fraction of concentrations have above min.cfu counts or if the average count across all concentrations is below min.average
    return(list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA)))
  }else{
    fit <- nls(y ~ a*(x/(x+Kd))+b,
               start=list(a=a.start,b=b.start,Kd=Kd.start),
               lower=list(a=a.lower,b=b.lower,Kd=min(x.vals[x.vals>0])/100), #constrain Kd to be no lower than 1/100x the lowest concentration value
               upper=list(a=a.upper,b=b.upper,Kd=max(x.vals[x.vals>0])*10), #constrain Kd to be no higher than the 10x highest concentration value
               algorithm="port")
    y.pred <- predict(fit,newdata=list(x=x))
    resid <- y - y.pred
    resid.norm <- resid/as.numeric(summary(fit)$coefficients["a","Estimate"])
    nMSR <- mean((resid.norm)^2,na.rm=T)
    return(list(as.numeric(summary(fit)$coefficients["Kd","Estimate"]),
                as.numeric(summary(fit)$coefficients["Kd","Std. Error"]),
                as.numeric(summary(fit)$coefficients["a","Estimate"]),
                as.numeric(summary(fit)$coefficients["b","Estimate"]),
                as.numeric(nMSR)))
  }
}

#fit titration to huACE2 Titeseq data for each barcode
dt[,c("Kd_ACE2","Kd_SE_ACE2","response_ACE2","baseline_ACE2","nMSR_ACE2") :=
     tryCatch(fit.titration(y.vals=c(TiteSeq_01_meanbin,TiteSeq_02_meanbin,TiteSeq_03_meanbin,TiteSeq_04_meanbin,
                                     TiteSeq_05_meanbin,TiteSeq_06_meanbin,TiteSeq_07_meanbin,TiteSeq_08_meanbin,
                                     TiteSeq_09_meanbin),
                            x.vals=samples_TiteSeq$conc,
                            count.vals=c(TiteSeq_01_totalcount,TiteSeq_02_totalcount,TiteSeq_03_totalcount,TiteSeq_04_totalcount,
                                         TiteSeq_05_totalcount,TiteSeq_06_totalcount,TiteSeq_07_totalcount,TiteSeq_08_totalcount,TiteSeq_09_totalcount)),
              error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))}),by=c("library","barcode")]
```

## QC and sanity checks

We will do some QC to make sure we got good titration curves for most of
our library barcodes. We will also spot check titration curves from
across our measurement range, and spot check curves whose fit parameters
hit the different boundary conditions of the fit variables.

We successfully generated *K*<sub>D</sub> estimates for 216040 of our
pool1A barcodes (73.12%), and 220272 of our pool2A barcodes (72.71%).

Why were estimates not returned for some barcodes? The histograms below
show that many barcodes with unsuccessful titration fits have lower
average cell counts and more concentrations with fewer than the minimum
cutoff number of cells (cutoff=2) than those that were fit. Therefore,
we can see the the majority of unfit barcodes come from our minimum read
cutoffs, meaning there weren’t too many curves that failed to be fit for
issues such as nls convergence.

``` r
par(mfrow=c(2,2))
hist(log10(dt[library=="pool1A" & !is.na(Kd_ACE2),TiteSeq_avgcount]+0.5),breaks=20,xlim=c(0,5),main="pool1A",col="gray50",xlab="average cell count across concentration samples")
hist(log10(dt[library=="pool1A" & is.na(Kd_ACE2),TiteSeq_avgcount]+0.5),breaks=20,add=T,col="red")

hist(log10(dt[library=="pool2A" & !is.na(Kd_ACE2),TiteSeq_avgcount]+0.5),breaks=20,xlim=c(0,5),main="pool2A",col="gray50",xlab="average cell count across concentration samples")
hist(log10(dt[library=="pool2A" & is.na(Kd_ACE2),TiteSeq_avgcount]+0.5),breaks=20,add=T,col="red")

hist(dt[library=="pool1A" & !is.na(Kd_ACE2),TiteSeq_min_cell_filtered],breaks=5,main="pool1A",col="gray50",xlab="number of sample concentrations below cutoff cell number",xlim=c(0,10))
hist(dt[library=="pool1A" & is.na(Kd_ACE2),TiteSeq_min_cell_filtered],breaks=16,add=T,col="red")

hist(dt[library=="pool2A" & !is.na(Kd_ACE2),TiteSeq_min_cell_filtered],breaks=5,main="pool2A",col="gray50",xlab="number of sample concentrations below cutoff cell number",xlim=c(0,10))
hist(dt[library=="pool2A" & is.na(Kd_ACE2),TiteSeq_min_cell_filtered],breaks=16,add=T,col="red")
```

<img src="compute_binding_Kd_files/figure-gfm/avgcount-1.png" style="display: block; margin: auto;" />

Let’s checkout what the data looks like for some curves that didn’t
converge on a titration fit, different cutoffs, boudnary conditions,
etc. I define a function that take a row from the data table and plots
the meanbin estimates and the fit titration curve (if converged). This
allows for quick and easy troubleshooting and spot-checking of curves.

In the plots below for non-converging fits, we can see that the data
seem to have very low plateaus/signal over the concentration range and
perhaps some noise. I understand why they are difficult to fit, and I am
not worried by their exclusion, as I can’t by eye tell what their fit
should be hitting. My best guess is they would have a “response”
parameter lower than the minimum allowable, but that is also a hard Kd
then to estimate reliably so I’m ok not fitting these relatively small
number of curves.

To allow manual checks of what the data looks like for different curve
fits, I define functions that take a row from the dt table and the
corresponding table of fits, and plots the meanbin estimates and the fit
titration curve (if converged). This allows for quick and easy
troubleshooting and spot-checking of curves.

``` r
#make functions that allow me to plot a titration for any given row from the counts data frames, for spot checking curves
plot.titration <- function(row,output.text=F){
  y.vals <- c();for(sample in samples_TiteSeq$sample){y.vals <- c(y.vals,paste(sample,"_meanbin",sep=""))};y.vals <- unlist(dt[row,y.vals,with=F])
  x.vals <- samples_TiteSeq$conc
  count.vals <- c();for(sample in samples_TiteSeq$sample){count.vals <- c(count.vals,paste(sample,"_totalcount",sep=""))};count.vals <- unlist(dt[row,count.vals,with=F])
  if(dt[row,variant_class] %in% c("wildtype","synonymous")){
    title <- dt[row,target]
  }else{
    title <- paste(dt[row,target],dt[row,aa_substitutions])
  }
  indices <- count.vals>cutoff & !is.na(count.vals)
  y.vals <- y.vals[indices]
  x.vals <- x.vals[indices]
  plot(x.vals,y.vals,xlab="[ACE2] (M)",
       ylab="mean bin",log="x",ylim=c(1,4),xlim=c(1e-13,1e-6),pch=19,main=title)
  Kd_var <- "Kd_ACE2"
  fit <- nls(y.vals ~ a*(x.vals/(x.vals+Kd))+b,
             start=list(a=3,b=1,Kd=dt[row,get(Kd_var)]),
             lower=list(a=2,b=1,Kd=1e-15),
             upper=list(a=3,b=1.5,Kd=1e-5), #constrain Kd to be no higher than the 10x highest concentration value
             algorithm="port") 
  if(!is.na(dt[row,get(Kd_var)])){
    lines(10^c(seq(-13,-6,0.25)),predict(fit,newdata=list(x.vals=10^c(seq(-13,-6,0.25)))))
    legend("topleft",bty="n",cex=1,legend=paste("Kd",format(dt[row,get(Kd_var)],digits=3),"M"))
  }
  if(output.text==T){ #for troubleshooting and interactive work, output some info from the counts table for the given row
    vars <- c("library","barcode","target","variant_class","aa_substitutions","TiteSeq_avgcount","TiteSeq_min_cell_filtered","Kd_ACE2","Kd_SE_ACE2","baseline_ACE2","response_ACE2","nMSR_ACE2")
    return(dt[row,..vars])
  }
}
```

Distribution of Kd estimates, with wt/syn barcodes in purple:

``` r
par(mfrow=c(3,1))
hist(log10(dt[target=="Wuhan_Hu_1",Kd_ACE2]),col="gray40",breaks=60,xlab="log10(KD), ACE2 (M)",main="Wuhan_Hu_1",xlim=c(-13,-5))
hist(log10(dt[target=="Wuhan_Hu_1" & variant_class %in% (c("synonymous","wildtype")),Kd_ACE2]),col="#92278F",add=T,breaks=60)
hist(log10(dt[target=="BA1",Kd_ACE2]),col="gray40",breaks=60,xlab="log10(KD), ACE2 (M)",main="BA1",xlim=c(-13,-5))
hist(log10(dt[target=="BA1" & variant_class %in% (c("synonymous","wildtype")),Kd_ACE2]),col="#92278F",add=T,breaks=60)
hist(log10(dt[target=="BA2",Kd_ACE2]),col="gray40",breaks=60,xlab="log10(KD), ACE2 (M)",main="BA2",xlim=c(-13,-5))
hist(log10(dt[target=="BA2" & variant_class %in% (c("synonymous","wildtype")),Kd_ACE2]),col="#92278F",add=T,breaks=60)
```

<img src="compute_binding_Kd_files/figure-gfm/Kd_distribution-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$Titeseq_Kds_dir,"/hist_Kd-per-barcode.pdf",sep="")))
```

Some stop variants eked through our RBD+ selection, either perhaps
because of stop codon readthrough, improper PacBio sequence annotation,
or other weirdness. Either way, the vast majority of nonsense mutants
were purged before this step, and the remaining ones are biased toward
unreliable and so we remove them.

``` r
#remove stop variants, which even if they eke through, either a) still have low counts and give poor fits as a result, or b) seem to be either dubious PacBio calls (lower variant_call_support) or have late stop codons which perhaps don't totally ablate funciton. Either way, the vast majority were purged before this step and we don't want to deal with the remaining ones!
dt[variant_class == "stop",c("Kd_ACE2","Kd_SE_ACE2","response_ACE2","baseline_ACE2","nMSR_ACE2") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]
```

Let’s take a look at some of the curves with *K*<sub>D,app</sub> values
across this distribution to get a broad sense of how things look.

First, curves with *K*<sub>D,app</sub> fixed at the 10<sup>-5</sup>
maximum. We can see these are all flat-lined curves with no response.

``` r
par(mfrow=c(2,2))
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 9e-6)[1])
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 9e-6)[2])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 9e-6)[1])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 9e-6)[2])
```

<img src="compute_binding_Kd_files/figure-gfm/1e-5_Kd-1.png" style="display: block; margin: auto;" />

Next, with *K*<sub>D,app</sub> around 10<sup>-6</sup>

``` r
par(mfrow=c(2,2))
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 1e-6 & dt$Kd_ACE2 < 1.2e-6)[1])
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 1e-6 & dt$Kd_ACE2 < 1.2e-6)[2])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 1e-6 & dt$Kd_ACE2 < 1.2e-6)[1])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 1e-6 & dt$Kd_ACE2 < 1.2e-6)[2])
```

<img src="compute_binding_Kd_files/figure-gfm/1e-6_Kd-1.png" style="display: block; margin: auto;" />

With *K*<sub>D,app</sub> around 10<sup>-7</sup>, we seem to be picking
up more consistent binding signals, though there are some noisy curves.

``` r
par(mfrow=c(2,2))
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 1e-7 & dt$Kd_ACE2 < 1.2e-7)[1])
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 1e-7 & dt$Kd_ACE2 < 1.2e-7)[2])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 1e-7 & dt$Kd_ACE2 < 1.2e-7)[1])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 1e-7 & dt$Kd_ACE2 < 1.2e-7)[2])
```

<img src="compute_binding_Kd_files/figure-gfm/1e-7_Kd-1.png" style="display: block; margin: auto;" />

At *K*<sub>D,app</sub> of 10<sup>-8</sup>, we are likewise picking up
some signal, perhaps a bit less noise than the -8 curves

``` r
par(mfrow=c(2,2))
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 1e-8 & dt$Kd_ACE2 < 1.2e-8)[1])
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 1e-8 & dt$Kd_ACE2 < 1.2e-8)[2])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 1e-8 & dt$Kd_ACE2 < 1.2e-8)[1])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 1e-8 & dt$Kd_ACE2 < 1.2e-8)[2])
```

<img src="compute_binding_Kd_files/figure-gfm/1e-8_Kd-1.png" style="display: block; margin: auto;" />

Same at *K*<sub>D,app</sub> of 10<sup>-9</sup>.

``` r
par(mfrow=c(2,2))
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 1e-9 & dt$Kd_ACE2 < 1.2e-9)[1])
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 1e-9 & dt$Kd_ACE2 < 1.2e-9)[2])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 1e-9 & dt$Kd_ACE2 < 1.2e-9)[1])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 1e-9 & dt$Kd_ACE2 < 1.2e-9)[2])
```

<img src="compute_binding_Kd_files/figure-gfm/1e-9_Kd-1.png" style="display: block; margin: auto;" />

*K*<sub>D,app</sub> of 10<sup>-10</sup>

``` r
par(mfrow=c(2,2))
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 1e-10 & dt$Kd_ACE2 < 1.2e-10)[1])
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 1e-10 & dt$Kd_ACE2 < 1.2e-10)[2])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 1e-10 & dt$Kd_ACE2 < 1.2e-10)[1])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 1e-10 & dt$Kd_ACE2 < 1.2e-10)[2])
```

<img src="compute_binding_Kd_files/figure-gfm/1e-10_Kd-1.png" style="display: block; margin: auto;" />

*K*<sub>D,app</sub> \~ 10<sup>-11</sup>. This is higher affinity than
the main bulk of curves.

``` r
par(mfrow=c(2,2))
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 1e-11 & dt$Kd_ACE2 < 2e-11)[1])
plot.titration(which(dt$library=="pool1A" & dt$Kd_ACE2 > 1e-11 & dt$Kd_ACE2 < 2e-11)[2])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 1e-11 & dt$Kd_ACE2 < 2e-11)[1])
plot.titration(which(dt$library=="pool2A" & dt$Kd_ACE2 > 1e-11 & dt$Kd_ACE2 < 2e-11)[2])
```

<img src="compute_binding_Kd_files/figure-gfm/1e-11_Kd-1.png" style="display: block; margin: auto;" />

## Data filtering by fit quality

Next, let’s filter out poor fits using the value we previously computed,
the *normalized* mean square residual (nMSR). This metric computes the
residual between the observed response variable and that predicted from
the titration fit, normalizes this residual by the response range of the
titration fit (which is allowed to vary between 1.5 and 3 per the
titration fits above), and computes the mean-square of these normalized
residuals.

Look at nMSR metric versus avgcoutn value, and layer on value of nMSR
filtering based on 20x the global median (and percentage filtered from
each background). Filter to NA fits with nMSR above this cutoff

``` r
median.nMSR <- median(dt$nMSR_ACE2,na.rm=T)
threshold <- 40
par(mfrow=c(2,2))
for(bg in c("Wuhan_Hu_1","BA1","BA2")){
  plot(log10(dt[target==bg,TiteSeq_avgcount]),dt[target==bg,nMSR_ACE2],main=bg,pch=19,col="#00000010",xlab="average cell count (log10)",ylab="nMSR",xlim=c(0,6),ylim=c(0,0.5))
  abline(h=threshold*median.nMSR,col="red",lty=2)
  legend("topleft",bty="n",cex=1,legend=paste(format(100*nrow(dt[target==bg & nMSR_ACE2 > threshold*median.nMSR & !is.na(nMSR_ACE2),])/nrow(dt[target==bg & !is.na(nMSR_ACE2),]),digits=3),"%"))
}

dt[nMSR_ACE2 > threshold*median.nMSR,c("Kd_ACE2","Kd_SE_ACE2","response_ACE2","baseline_ACE2") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))]
```

<img src="compute_binding_Kd_files/figure-gfm/nMSR_v_cell_count-1.png" style="display: block; margin: auto;" />

Last, convert our *K*<sub>D,app</sub> to 1) a log<sub>10</sub>-scale,
and 2) *K*<sub>A,app</sub>, the inverse of *K*<sub>D,app</sub>, such
that higher values are associated with tighter binding, as is more
intuitive. (If we want to continue to discuss in terms of
*K*<sub>D,app</sub>, since people are often more familiar with
*K*<sub>D</sub>, we can refer to the
log<sub>10</sub>(*K*<sub>A,app</sub>) as
-log<sub>10</sub>(*K*<sub>D,app</sub>), which are identical.

``` r
dt[,log10Ka := -log10(Kd_ACE2),by=c("barcode","library")]
```

Let’s visualize the final binding measurements as violin plots for the
different wildtype targets. In next notebook, we’ll evaluate count depth
and possibly apply further filtering to remove low-count expression
estimates

``` r
p1 <- ggplot(dt[!is.na(log10Ka),],aes(x=variant_class,y=log10Ka))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("huACE2, log10(Ka)")+xlab("variant class")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
  facet_wrap(~target,nrow=4)

grid.arrange(p1,ncol=1)
```

<img src="compute_binding_Kd_files/figure-gfm/binding_distribution_vioplot-1.png" style="display: block; margin: auto;" />

``` r
#save pdf
invisible(dev.print(pdf, paste(config$Titeseq_Kds_dir,"/violin-plot_log10Ka-by-target.pdf",sep="")))
```

We have generated binding measurements for 72.68% of the barcodes in our
libraries.

## Data Output

Finally, let’s output our measurements for downstream analyses.

``` r
dt[,.(library,barcode,target,variant_class,aa_substitutions,n_aa_substitutions,
     TiteSeq_avgcount,log10Ka)] %>%
  mutate_if(is.numeric, round, digits=6) %>%
  write.csv(file=config$Titeseq_Kds_file, row.names=F)
```
