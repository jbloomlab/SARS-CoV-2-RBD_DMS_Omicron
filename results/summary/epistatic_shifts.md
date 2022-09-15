Shifts in mutation effects among variant backgrounds
================
Tyler Starr
6/9/2022

-   <a href="#setup" id="toc-setup">Setup</a>
-   <a
    href="#calculate-site-wise-jensen-shannon-divergence-a-metric-of-divergence-in-site-specific-mutational-profiles"
    id="toc-calculate-site-wise-jensen-shannon-divergence-a-metric-of-divergence-in-site-specific-mutational-profiles">Calculate
    site-wise Jensen-Shannon divergence, a metric of divergence in
    site-specific mutational profiles</a>
-   <a href="#line-plots-of-js-divergence-from-wh1-across-rbd-sites"
    id="toc-line-plots-of-js-divergence-from-wh1-across-rbd-sites">Line
    plots of JS divergence from WH1 across RBD sites</a>
-   <a href="#map-divergence-to-pdb-structure"
    id="toc-map-divergence-to-pdb-structure">Map divergence to pdb
    structure</a>

This notebook analyzes sites whose mutation effects deviate most
strongly among the variant RBD backgrounds.

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra","bio3d","ggridges","ggrepel","GGally")
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
if(!file.exists(config$epistatic_shifts_dir)){
  dir.create(file.path(config$epistatic_shifts_dir))
}

#make pdb output directory
if(!file.exists(paste(config$epistatic_shifts_dir,"/pdbs/",sep=""))){
  dir.create(file.path(paste(config$epistatic_shifts_dir,"/pdbs/",sep="")))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.6 LTS
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
    ##  [1] GGally_1.4.0      ggrepel_0.8.1     ggridges_0.5.1    bio3d_2.4-0      
    ##  [5] gridExtra_2.3     forcats_0.4.0     stringr_1.4.0     dplyr_0.8.3      
    ##  [9] purrr_0.3.3       readr_1.3.1       tidyr_1.0.0       tibble_3.0.2     
    ## [13] ggplot2_3.3.0     tidyverse_1.3.0   data.table_1.12.8 yaml_2.2.0       
    ## [17] knitr_1.26       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.0   xfun_0.11          haven_2.2.0        colorspace_1.4-1  
    ##  [5] vctrs_0.3.1        generics_0.0.2     htmltools_0.4.0    rlang_0.4.7       
    ##  [9] pillar_1.4.5       glue_1.3.1         withr_2.1.2        DBI_1.1.0         
    ## [13] RColorBrewer_1.1-2 dbplyr_1.4.2       modelr_0.1.5       readxl_1.3.1      
    ## [17] plyr_1.8.5         lifecycle_0.2.0    munsell_0.5.0      gtable_0.3.0      
    ## [21] cellranger_1.1.0   rvest_0.3.5        evaluate_0.14      parallel_3.6.2    
    ## [25] fansi_0.4.0        broom_0.7.0        Rcpp_1.0.3         scales_1.1.0      
    ## [29] backports_1.1.5    jsonlite_1.6       fs_1.3.1           hms_0.5.2         
    ## [33] digest_0.6.23      stringi_1.4.3      grid_3.6.2         cli_2.0.0         
    ## [37] tools_3.6.2        magrittr_1.5       crayon_1.3.4       pkgconfig_2.0.3   
    ## [41] ellipsis_0.3.0     xml2_1.2.2         reprex_0.3.0       lubridate_1.7.4   
    ## [45] reshape_0.8.8      assertthat_0.2.1   rmarkdown_2.0      httr_1.4.1        
    ## [49] rstudioapi_0.10    R6_2.4.1           compiler_3.6.2

Define colorblind-friendly palette

``` r
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")
# The palette with black
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")
```

## Setup

Read in tables of mutant measurements from current study and prior VOCs
DMS datasets.

``` r
dt <- data.table(read.csv(file=config$final_variant_scores_mut_file,stringsAsFactors=F))
```

## Calculate site-wise Jensen-Shannon divergence, a metric of divergence in site-specific mutational profiles

For each pair of backgrounds, at each site I want to compute the
Jensen-Shannon divergence between the profile of mutation effects of all
mutations at the site. (Remove any measurements determined for \<3 or
\<5 bc to avoid measurements with lower precision driving noise in the
global JSD metric.)

``` r
#define a minbc bind measurement that requires three barcodes be averaged for a final determination, otherwise change to NA
dt[,bind_min3bc := bind]
dt[n_bc_bind < 3, bind_min3bc := NA]

dt[,bind_min5bc := bind]
dt[n_bc_bind < 5, bind_min5bc := NA]

#define a minbc expr measurement that requires three barcodes be averaged for a final determination, otherwise change to NA
dt[,expr_min3bc := expr]
dt[n_bc_expr < 3, expr_min3bc := NA]

dt[,expr_min5bc := expr]
dt[n_bc_expr < 5, expr_min5bc := NA]

#define a function for computing J-S divergence/distance between two affinity vectors (where affinity is given currently as log10-Kd)
JS <- function(vec1,vec2){
  vec1_pair <- vec1[!is.na(vec1) & !is.na(vec2)]
  vec2_pair <- vec2[!is.na(vec1) & !is.na(vec2)]
  pi1 <- 10^(vec1_pair)/sum(10^(vec1_pair))
  pi2 <- 10^(vec2_pair)/sum(10^(vec2_pair))
  n <- 0.5 * (pi1+pi2)
  JS <- 0.5 * (sum(pi1*log(pi1/n)) + sum(pi2*log(pi2/n)))
  #return(sqrt(JS)) #if doing distance
  return(JS) #if doing divergence
}

#first, for bind measurements
#data table for storing difference in correlation in profiles between bg pairs at each site
#generate table with all combinations of bg_1 and bg_2 for each site
diffs_bind <- data.table(expand.grid(site=unique(dt$position),bg_2=c("Wuhan-Hu-1_v2","Omicron_BA1","Omicron_BA2","Wuhan-Hu-1_v1","Alpha","Beta","Delta","Eta"),bg_1=c("Wuhan-Hu-1_v2","Omicron_BA1","Omicron_BA2","Wuhan-Hu-1_v1","Alpha","Beta","Delta","Eta")))

#remove duplicates where bg_1 and _2 the same
diffs_bind <- diffs_bind[bg_1 != bg_2,]

#loop through and compute JSD for each site for each pair of bgs, for bind metric
diffs_bind$JSD <- as.numeric(NA) #jensen-shannon divergence, from raw bind values (lower limit 5)
diffs_bind$JSD_min3bc <- as.numeric(NA) #jensen-shannon divergence, require a minimum of 3 bcs averaged
diffs_bind$JSD_min5bc <- as.numeric(NA) #jensen-shannon divergence, require a minimum of 5 bcs averaged
for(i in 1:nrow(diffs_bind)){
  x_uncens <- dt[target==diffs_bind[i,bg_1] & position==diffs_bind[i,site],bind]
  y_uncens <- dt[target==diffs_bind[i,bg_2] & position==diffs_bind[i,site],bind]
  x_min3bc <- dt[target==diffs_bind[i,bg_1] & position==diffs_bind[i,site],bind_min3bc]
  y_min3bc <- dt[target==diffs_bind[i,bg_2] & position==diffs_bind[i,site],bind_min3bc]
  x_min5bc <- dt[target==diffs_bind[i,bg_1] & position==diffs_bind[i,site],bind_min5bc]
  y_min5bc <- dt[target==diffs_bind[i,bg_2] & position==diffs_bind[i,site],bind_min5bc]
  diffs_bind[i,JSD := JS(x_uncens,y_uncens)]
  diffs_bind[i,JSD_min3bc := JS(x_min3bc,y_min3bc)]
  diffs_bind[i,JSD_min5bc := JS(x_min3bc,y_min5bc)]
}

#repeat for expr measurements
#data table for storign difference in correlation in profiles between bg pairs at each site
#generate table with all combinations of bg_1 and bg_2 for each site
diffs_expr <- data.table(expand.grid(site=unique(dt$position),bg_2=c("Wuhan-Hu-1_v2","Omicron_BA1","Omicron_BA2","Wuhan-Hu-1_v1","Alpha","Beta","Delta","Eta"),bg_1=c("Wuhan-Hu-1_v2","Omicron_BA1","Omicron_BA2","Wuhan-Hu-1_v1","Alpha","Beta","Delta","Eta")))

#remove duplicates where either bg_1 and _2 the same
diffs_expr <- diffs_expr[bg_1 != bg_2,]

#loop through and compute JSD for each site for each pair of bgs, for expr metric
diffs_expr$JSD <- as.numeric(NA) #jensen-shannon divergence, from raw expr values
diffs_expr$JSD_min3bc <- as.numeric(NA) #jensen-shannon divergence, require a minimum of 3 bcs averaged
diffs_expr$JSD_min5bc <- as.numeric(NA) #jensen-shannon divergence, require a minimum of 5 bcs averaged
for(i in 1:nrow(diffs_expr)){
  x_uncens <- dt[target==diffs_expr[i,bg_1] & position==diffs_expr[i,site],expr]
  y_uncens <- dt[target==diffs_expr[i,bg_2] & position==diffs_expr[i,site],expr]
  x_min3bc <- dt[target==diffs_expr[i,bg_1] & position==diffs_expr[i,site],expr_min3bc]
  y_min3bc <- dt[target==diffs_expr[i,bg_2] & position==diffs_expr[i,site],expr_min3bc]
  x_min5bc <- dt[target==diffs_expr[i,bg_1] & position==diffs_expr[i,site],expr_min5bc]
  y_min5bc <- dt[target==diffs_expr[i,bg_2] & position==diffs_expr[i,site],expr_min5bc]
  diffs_expr[i,JSD := JS(x_uncens,y_uncens)]
  diffs_expr[i,JSD_min3bc := JS(x_min3bc,y_min3bc)]
  diffs_expr[i,JSD_min5bc := JS(x_min3bc,y_min5bc)]
}
```

Output file with the site-pair JS divergences.

``` r
diffs_bind[,.(bg_1,bg_2,site,JSD,JSD_min3bc,JSD_min5bc)] %>%
  mutate_if(is.numeric, round, digits=6) %>%
  write.csv(file=config$JSD_v_WH1_file, row.names=F,quote=F)
```

Output file with the site-pair JS divergences.

``` r
diffs_expr[,.(bg_1,bg_2,site,JSD,JSD_min3bc,JSD_min5bc)] %>%
  mutate_if(is.numeric, round, digits=6) %>%
  write.csv(file=config$JSD_v_WH1_expr_file, row.names=F,quote=F)
```

Plotting/visualizing:

Utility function: plot scatterplot showing affinity of each of the 20
amino acids in a pair of sites

``` r
plot_scatter <- function(site, bg1, bg2, JSD=F, JSD_min3bc=T, JSD_min5bc=F,n_bc_cutoff=3,phenotype="bind"){
  x <- dt[target==bg1 & position==site,get(phenotype)]
  x_n_bc <- dt[target==bg1 & position==site,get(paste("n_bc_",phenotype,sep=""))]
  x_ref <- dt[target==bg1 & position==site & as.character(mutant)==as.character(wildtype),get(phenotype)]
  y <- dt[target==bg2 & position==site,get(phenotype)]
  y_n_bc <- dt[target==bg2 & position==site,get(paste("n_bc_",phenotype,sep=""))]
  y_ref <- dt[target==bg2 & position==site & as.character(mutant)==as.character(wildtype),get(phenotype)]
  x_min3bc <- dt[target==bg1 & position==site,get(paste(phenotype,"_min3bc",sep=""))]
  y_min3bc <- dt[target==bg2 & position==site,get(paste(phenotype,"_min3bc",sep=""))]
  x_min5bc <- dt[target==bg1 & position==site,get(paste(phenotype,"_min5bc",sep=""))]
  y_min5bc <- dt[target==bg2 & position==site,get(paste(phenotype,"_min5bc",sep=""))]
  chars <- dt[target==bg1 & position==site,mutant]
  cols <- rep("black",20); cols[which(x_n_bc < n_bc_cutoff | y_n_bc < n_bc_cutoff)] <- "orange"
  plot(x,y, xlim=if(phenotype=="bind"){c(4.5,12)}else{c(5.5,11)},ylim=if(phenotype=="bind"){c(4.5,12)}else{c(5.5,11)},pch=chars,xlab=paste(bg1,phenotype),ylab=paste(bg2,phenotype),col=cols,main=paste("site",site))
  abline(v=x_ref,lty=2,col="red")
  abline(h=y_ref,lty=2,col="red")
  if(JSD==T){
    val <- JS(x,y)
    legend("topleft",bty="n",cex=1,legend=paste("JSD:",format(val,digits=3)))
  }else if(JSD_min3bc==T){
    val <- JS(x_min3bc,y_min3bc)
    legend("topleft",bty="n",cex=1,legend=paste("JSD:",format(val,digits=3)))
  }else if(JSD_min5bc==T){
    val <- JS(x_min5bc,y_min5bc)
    legend("topleft",bty="n",cex=1,legend=paste("JSD:",format(val,digits=3)))
  }
}
```

``` r
par(mfrow=c(4,3))
plot_scatter(site=403,"Wuhan-Hu-1_v2","Omicron_BA1")
plot_scatter(site=406,"Wuhan-Hu-1_v2","Omicron_BA1")
plot_scatter(site=419,"Wuhan-Hu-1_v2","Omicron_BA1")
plot_scatter(site=439,"Wuhan-Hu-1_v2","Omicron_BA1")
plot_scatter(site=449,"Wuhan-Hu-1_v2","Omicron_BA1")
plot_scatter(site=455,"Wuhan-Hu-1_v2","Omicron_BA1")
plot_scatter(site=493,"Wuhan-Hu-1_v2","Omicron_BA1")
plot_scatter(site=494,"Wuhan-Hu-1_v2","Omicron_BA1")
plot_scatter(site=496,"Wuhan-Hu-1_v2","Omicron_BA1")
plot_scatter(site=498,"Wuhan-Hu-1_v2","Omicron_BA1")
plot_scatter(site=501,"Wuhan-Hu-1_v2","Omicron_BA1")
plot_scatter(site=505,"Wuhan-Hu-1_v2","Omicron_BA1")
```

<img src="epistatic_shifts_files/figure-gfm/scatters_BA1_diffs-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/bg-scatters_Omicron_BA1_shifts.pdf",sep=""),useDingbats=F))
```

``` r
par(mfrow=c(4,3))
plot_scatter(site=403,"Wuhan-Hu-1_v2","Omicron_BA2")
plot_scatter(site=406,"Wuhan-Hu-1_v2","Omicron_BA2")
plot_scatter(site=419,"Wuhan-Hu-1_v2","Omicron_BA2")
plot_scatter(site=439,"Wuhan-Hu-1_v2","Omicron_BA2")
plot_scatter(site=449,"Wuhan-Hu-1_v2","Omicron_BA2")
plot_scatter(site=455,"Wuhan-Hu-1_v2","Omicron_BA2")
plot_scatter(site=493,"Wuhan-Hu-1_v2","Omicron_BA2")
plot_scatter(site=494,"Wuhan-Hu-1_v2","Omicron_BA2")
plot_scatter(site=496,"Wuhan-Hu-1_v2","Omicron_BA2")
plot_scatter(site=498,"Wuhan-Hu-1_v2","Omicron_BA2")
plot_scatter(site=501,"Wuhan-Hu-1_v2","Omicron_BA2")
plot_scatter(site=505,"Wuhan-Hu-1_v2","Omicron_BA2")
```

<img src="epistatic_shifts_files/figure-gfm/scatters_BA2_diffs-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/bg-scatters_Omicron_BA2_shifts.pdf",sep=""),useDingbats=F))
```

## Line plots of JS divergence from WH1 across RBD sites

Make lineplots showing JS-D across sites for each variant compared to
WH1.

Also add gray shading for sites of escape from antibodies from our large
panel of antibodies we’ve profiled w.r.t. WH1 escape, downloaded from:
<https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/main/processed_data/escape_data.csv>

First, define sites of substantial antibody escape

``` r
dt_mAb <- data.table(read.csv(file=config$mut_antibody_escape,stringsAsFactors = F))
dt_mAb <- unique(dt_mAb[condition_type=="antibody",.(condition, condition_type, condition_subtype, site, wildtype, site_total_escape)])

dt_mAb[,site_average_escape:=mean(site_total_escape,na.rm=T),by=c("site")]

site_escape <- unique(dt_mAb[,.(wildtype, site, site_average_escape)])

#define sites for labeling as those with an average of 0.05 normalized site-wise escape across all mAbs
sig_mAb_sites <- site_escape[site_average_escape>0.15, site]


#define some epitope classes for adding highlights
label_df <- data.frame(xmin=sig_mAb_sites-0.5,
                       xmax=sig_mAb_sites+0.5)
```

``` r
#define focal bg for others to compare to
bg <- "Wuhan-Hu-1_v2"
temp <- diffs_bind[bg_1==bg,]
temp$target <- as.character(temp$bg_2)

#define colors for each bg
group.colors <- c("Wuhan-Hu-1_v2" = cbPalette[1], "Omicron_BA1" = cbPalette[2], "Omicron_BA2" = cbPalette[8],"Wuhan-Hu-1_v1" = "black", "Alpha" = cbPalette[3], "Beta" = cbPalette[6], "Delta" = cbPalette[7], "Eta" = cbPalette[5])

#define order for plotting of bgs
temp$target <- factor(temp$target,levels=c("Alpha","Beta","Delta","Eta","Wuhan-Hu-1_v1","Omicron_BA1","Omicron_BA2"))


ggplot(data=temp, aes(x=site, y=JSD, color=target))+
  geom_rect(data=label_df, aes(x=NULL, y=NULL, color=NULL,xmin=xmin, xmax=xmax, ymin=0,ymax=1.1*max(temp$JSD,na.rm=T)), alpha=0.2)+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("JS divergence versus Wuhan-Hu-1")
```

<img src="epistatic_shifts_files/figure-gfm/line_plots_JSD_v_WH1-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/JSD_v_WH1.pdf",sep=""),useDingbats=F))
```

Same but require minimum 3 bc for a measurement

``` r
ggplot(data=temp, aes(x=site, y=JSD_min3bc, color=target))+
  geom_rect(data=label_df, aes(x=NULL, y=NULL, color=NULL,xmin=xmin, xmax=xmax, ymin=0,ymax=1.1*max(temp$JSD,na.rm=T)), alpha=0.2)+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("JS divergence versus Wuhan-Hu-1")+
  geom_text_repel(aes(label=ifelse(((JSD_min3bc > 0.15)),as.character(site),'')),size=3,color="gray40")
```

<img src="epistatic_shifts_files/figure-gfm/line_plots_JSD_v_WH1_min3bc-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/JSD_v_WH1_min3bc.pdf",sep=""),useDingbats=F))
```

Include just Beta from the prior VOC dataset, and Beta-associated WH1
for “control”

``` r
ggplot(data=temp[bg_2 %in% c("Omicron_BA1","Omicron_BA2","Beta","Wuhan-Hu-1_v1")], aes(x=site, y=JSD_min3bc, color=target))+
  geom_rect(data=label_df, aes(x=NULL, y=NULL, color=NULL,xmin=xmin, xmax=xmax, ymin=0,ymax=1.1*max(temp$JSD,na.rm=T)), alpha=0.2)+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("JS divergence versus Wuhan-Hu-1")+
  geom_text_repel(aes(label=ifelse(((JSD_min3bc > 0.15)),as.character(site),'')),size=3,color="gray40")
```

<img src="epistatic_shifts_files/figure-gfm/line_plots_JSD_v_WH1_min3bc_Beta-only-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/JSD_v_WH1_min3bc_Omi-and-Beta-only.pdf",sep=""),useDingbats=F))
```

Scatterplots to more clearly show the additional shifts in each Omicron
variant versus Beta?

``` r
temp2 <- dcast(temp[bg_2 %in% c("Omicron_BA1","Omicron_BA2","Beta")], site ~ bg_2, value.var="JSD_min3bc")

p1 <- ggplot(data=temp2, aes(x=Omicron_BA1, y=Omicron_BA2))+
  geom_point()+
  geom_text_repel(aes(label=ifelse((Omicron_BA2 > 0.1 | Omicron_BA1 > 0.1),as.character(site),'')),size=3)+
  theme_classic()+
  xlab("epistatic shift versus Wuhan-Hu-1, Omicron BA.1")+
  ylab("epistatic shift versus Wuhan-Hu-1, Omicron BA.2")

p2 <- ggplot(data=temp2, aes(x=Omicron_BA1, y=Beta))+
  geom_point()+
  geom_text_repel(aes(label=ifelse((Omicron_BA1 > 0.1 | Beta > 0.1),as.character(site),'')),size=3)+
  theme_classic()+
  xlab("epistatic shift versus Wuhan-Hu-1, Omicron BA.1")+
  ylab("epistatic shift versus Wuhan-Hu-1, Beta")

p3 <- ggplot(data=temp2, aes(x=Omicron_BA2, y=Beta))+
  geom_point()+
  geom_text_repel(aes(label=ifelse((Omicron_BA2 > 0.1 | Beta > 0.1),as.character(site),'')),size=3)+
  theme_classic()+
  xlab("epistatic shift versus Wuhan-Hu-1, Omicron BA.2")+
  ylab("epistatic shift versus Wuhan-Hu-1, Beta")

grid.arrange(p1,p2,p3,nrow=1)
```

<img src="epistatic_shifts_files/figure-gfm/scatterplots_epistatic-shift_Beta-BA1-BA2-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/JSD_v_WH1_min3bc_Omi-versus-Beta-scatters.pdf",sep=""),useDingbats=F))
```

Repeat for expression measurements

``` r
#define focal bg for others to compare to
bg <- "Wuhan-Hu-1_v2"
temp <- diffs_expr[bg_1==bg,]
temp$target <- as.character(temp$bg_2)

#define order for plotting of bgs
temp$target <- factor(temp$target,levels=c("Alpha","Beta","Delta","Eta","Wuhan-Hu-1_v1","Omicron_BA1","Omicron_BA2"))


ggplot(data=temp, aes(x=site, y=JSD, color=target))+
  geom_rect(data=label_df, aes(x=NULL, y=NULL, color=NULL,xmin=xmin, xmax=xmax, ymin=0,ymax=1.1*max(temp$JSD,na.rm=T)), alpha=0.2)+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("JS divergence versus Wuhan-Hu-1")
```

<img src="epistatic_shifts_files/figure-gfm/line_plots_JSD_v_WH1_expr-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/JSD_v_WH1_expr.pdf",sep=""),useDingbats=F))
```

Same but require minimum 3 bc for a measurement

``` r
ggplot(data=temp, aes(x=site, y=JSD_min3bc, color=target))+
  geom_rect(data=label_df, aes(x=NULL, y=NULL, color=NULL,xmin=xmin, xmax=xmax, ymin=0,ymax=1.1*max(temp$JSD,na.rm=T)), alpha=0.2)+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("JS divergence versus Wuhan-Hu-1, expression DMS")+
  geom_text_repel(aes(label=ifelse(((JSD_min3bc > 0.15)),as.character(site),'')),size=3,color="gray40")
```

<img src="epistatic_shifts_files/figure-gfm/line_plots_JSD_v_WH1_min3bc_expr-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/JSD_v_WH1_min3bc_expr.pdf",sep=""),useDingbats=F))
```

Include just Beta from the prior VOC dataset, and Beta-associated
Wuhan-Hu-1 measurements for “control”

``` r
ggplot(data=temp[bg_2 %in% c("Omicron_BA1","Omicron_BA2","Beta","Wuhan-Hu-1_v1")], aes(x=site, y=JSD_min3bc, color=target))+
  geom_rect(data=label_df, aes(x=NULL, y=NULL, color=NULL,xmin=xmin, xmax=xmax, ymin=0,ymax=1.1*max(temp$JSD,na.rm=T)), alpha=0.2)+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("JS divergence versus Wuhan-Hu-1")+
  geom_text_repel(aes(label=ifelse(((JSD_min3bc > 0.1)),as.character(site),'')),size=3,color="gray40")
```

<img src="epistatic_shifts_files/figure-gfm/line_plots_JSD_v_WH1_expr_min3bc_Beta-only-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/JSD_v_WH1_min3bc_expr_Omi-and-Beta-only.pdf",sep=""),useDingbats=F))
```

``` r
temp[bg_2=="Wuhan-Hu-1_v1",bg_2:="Wuhan.Hu.1_v1"]
temp2 <- dcast(temp[bg_2 %in% c("Omicron_BA1","Omicron_BA2","Beta","Wuhan.Hu.1_v1")], site ~ bg_2, value.var="JSD_min3bc")

p1 <- ggplot(data=temp2, aes(x=Omicron_BA1, y=Omicron_BA2))+
  geom_point()+
  geom_text_repel(aes(label=ifelse((Omicron_BA2 > 0.1 | Omicron_BA1 > 0.1),as.character(site),'')),size=3)+
  theme_classic()+
  scale_x_continuous(limits = c(0, 0.31))+
  scale_y_continuous(limits = c(0, 0.31))

p2 <- ggplot(data=temp2, aes(x=Omicron_BA1, y=Beta))+
  geom_point()+
  geom_text_repel(aes(label=ifelse((Omicron_BA1 > 0.1 | Beta > 0.1),as.character(site),'')),size=3)+
  theme_classic()+
  scale_x_continuous(limits = c(0, 0.31))+
  scale_y_continuous(limits = c(0, 0.31))

p3 <- ggplot(data=temp2, aes(x=Omicron_BA2, y=Beta))+
  geom_point()+
  geom_text_repel(aes(label=ifelse((Omicron_BA2 > 0.1 | Beta > 0.1),as.character(site),'')),size=3)+
  theme_classic()+
  scale_x_continuous(limits = c(0, 0.31))+
  scale_y_continuous(limits = c(0, 0.31))

p4 <- ggplot(data=temp2, aes(x=Omicron_BA1, y=Wuhan.Hu.1_v1))+
  geom_point()+
  geom_text_repel(aes(label=ifelse((Wuhan.Hu.1_v1 > 0.1 | Beta > 0.1),as.character(site),'')),size=3)+
  theme_classic()+
  scale_x_continuous(limits = c(0, 0.31))+
  scale_y_continuous(limits = c(0, 0.31))

grid.arrange(p1,p2,p3,p4,nrow=1)
```

<img src="epistatic_shifts_files/figure-gfm/scatterplots_epistatic-shift_Beta-BA1-BA2_expr-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/JSD_v_WH1_min3bc_expr_Omi-versus-Beta-scatters.pdf",sep=""),useDingbats=F))
```

## Map divergence to pdb structure

First, bind

``` r
pdb_wh1 <- read.pdb(file=config$pdb_6m0j)
```

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
#iterate through backgrounds, output a pdb comparing its divergence to WH1 (using min3bc)
for(s in c("Omicron_BA1","Omicron_BA2")){
  b <- rep(0, length(pdb_wh1$atom$b))
  for(i in 1:nrow(pdb_wh1$atom)){
    if(pdb_wh1$atom$chain[i]=="E"){
      res <- pdb_wh1$atom$resno[i]
      JSD <- diffs_bind[bg_1=="Wuhan-Hu-1_v2" & bg_2==s & site==res, JSD_min3bc]
      if(length(JSD)>0){
        b[i] <- JSD
      }
    }
  }
  write.pdb(pdb=pdb_wh1, file=paste(config$epistatic_shifts_dir,"/pdbs/",s,"_v_WH1_JSD-min3bc.pdb",sep=""), b=b)
}
```

repeat for expression measures

``` r
pdb_wh1 <- read.pdb(file=config$pdb_6m0j)
```

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
#iterate through backgrounds, output a pdb comparing its divergence to WH1 (using min3bc)
for(s in c("Omicron_BA1","Omicron_BA2")){
  b <- rep(0, length(pdb_wh1$atom$b))
  for(i in 1:nrow(pdb_wh1$atom)){
    if(pdb_wh1$atom$chain[i]=="E"){
      res <- pdb_wh1$atom$resno[i]
      JSD <- diffs_expr[bg_1=="Wuhan-Hu-1_v2" & bg_2==s & site==res, JSD_min3bc]
      if(length(JSD)>0){
        b[i] <- JSD
      }
    }
  }
  write.pdb(pdb=pdb_wh1, file=paste(config$epistatic_shifts_dir,"/pdbs/",s,"_v_WH1_JSD-min3bc_expr.pdb",sep=""), b=b)
}
```

``` r
par(mfrow=c(1,2))
temp <- dt[position==493 & target %in% c("Omicron_BA1","Wuhan-Hu-1_v2"),]

temp[,delta_bind_Q493ref := as.numeric(NA)]
for(i in 1:nrow(temp)){
  temp$delta_bind_Q493ref[i] <- temp$bind[i] - temp[target==temp$target[i] & mutant=="Q",bind] 
}

plot(temp[target=="Wuhan-Hu-1_v2",delta_bind_Q493ref],temp[target=="Omicron_BA1",delta_bind_Q493ref],pch=temp[target=="Wuhan-Hu-1_v2",mutant],xlab="Q493x mutation, WH1",ylab="Q493x mutation, BA1",xlim=c(-2.5,1.7),ylim=c(-2.5,1.7))
abline(a=0,b=1)
abline(v=0,lty=2)
abline(h=0,lty=2)

temp2 <- dt[position==493 & target %in% c("Omicron_BA2","Wuhan-Hu-1_v2"),]

temp2[,delta_bind_Q493ref := as.numeric(NA)]
for(i in 1:nrow(temp2)){
  temp2$delta_bind_Q493ref[i] <- temp2$bind[i] - temp2[target==temp2$target[i] & mutant=="Q",bind] 
}

plot(temp2[target=="Wuhan-Hu-1_v2",delta_bind_Q493ref],temp2[target=="Omicron_BA2",delta_bind_Q493ref],pch=temp2[target=="Wuhan-Hu-1_v2",mutant],xlab="Q493x mutation, WH1",ylab="Q493x mutation, BA1",xlim=c(-2.5,1.7),ylim=c(-2.5,1.7))
abline(a=0,b=1)
abline(v=0,lty=2)
abline(h=0,lty=2)
```

<img src="epistatic_shifts_files/figure-gfm/look more at 493-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/scatter_Q493x_mut.pdf",sep=""),useDingbats=F))
```

For each position that varies between WH1 and BA.1 and BA.2, make little
plots that illustrate epistasis of the substitution itself to show
entrenchment and anti-entrenchment

``` r
entrench_BA1 <- data.table(site=NA,WH1_wt=NA,BA1_wt=NA,WH1_wt_bind=NA,WH1_mut_bind=NA,BA1_revert_bind=NA,BA1_wt_bind=NA)

for(pos in unique(dt$position)){
  WH1_wt <- dt[target=="Wuhan-Hu-1_v2" & wildtype==mutant & position==pos,wildtype]
  BA1_wt <- dt[target=="Omicron_BA1" & wildtype==mutant & position==pos,wildtype]
  if(WH1_wt != BA1_wt){
    entrench_BA1 <- rbind(entrench_BA1, 
                          list(site=pos,
                               WH1_wt=WH1_wt,
                               BA1_wt=BA1_wt,
                               WH1_wt_bind=dt[target=="Wuhan-Hu-1_v2" & position==pos & mutant==WH1_wt,bind],
                               WH1_mut_bind=dt[target=="Wuhan-Hu-1_v2" & position==pos & mutant==BA1_wt,bind],
                               BA1_revert_bind=dt[target=="Omicron_BA1" & position==pos & mutant==WH1_wt,bind],
                               BA1_wt_bind=dt[target=="Omicron_BA1" & position==pos & mutant==BA1_wt,bind]
                               )
                          )
  }
}
entrench_BA1 <- entrench_BA1[-1,]

par(mfrow=c(5,3))
for(i in 1:nrow(entrench_BA1)){
  plot(1:4, 
       entrench_BA1[i,4:7],
       main=paste(entrench_BA1[i,c(2,1,3)],collapse=""),
       pch=16,cex=2,
       col=c(cbbPalette[1],cbbPalette[1],cbbPalette[2],cbbPalette[2]),
       xlim=c(0.75, 4.25), ylim=c(6,10.5),
       ylab="ACE2 affinity (-log10Kd)",
       xlab="amino acid",
       xaxt="n")
  axis(1, at=1:4, labels=c("WH1", "WH1+\nmut","BA.1+\nrevert","BA.1"))
  points(1:2, entrench_BA1[i,4:5],type="l",lwd=1.5)
  points(3:4, entrench_BA1[i,6:7],type="l",lwd=1.5,col=cbbPalette[2])
}
```

<img src="epistatic_shifts_files/figure-gfm/entrenchment_Omicron_BA1-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/entrenchment_diagrams_BA1.pdf",sep=""),useDingbats=F))
```

``` r
entrench_BA2 <- data.table(site=NA,WH1_wt=NA,BA2_wt=NA,WH1_wt_bind=NA,WH1_mut_bind=NA,BA2_revert_bind=NA,BA2_wt_bind=NA)

for(pos in unique(dt$position)){
  WH1_wt <- dt[target=="Wuhan-Hu-1_v2" & wildtype==mutant & position==pos,wildtype]
  BA2_wt <- dt[target=="Omicron_BA2" & wildtype==mutant & position==pos,wildtype]
  if(WH1_wt != BA2_wt){
    entrench_BA2 <- rbind(entrench_BA2, 
                          list(site=pos,
                               WH1_wt=WH1_wt,
                               BA2_wt=BA2_wt,
                               WH1_wt_bind=dt[target=="Wuhan-Hu-1_v2" & position==pos & mutant==WH1_wt,bind],
                               WH1_mut_bind=dt[target=="Wuhan-Hu-1_v2" & position==pos & mutant==BA2_wt,bind],
                               BA2_revert_bind=dt[target=="Omicron_BA2" & position==pos & mutant==WH1_wt,bind],
                               BA2_wt_bind=dt[target=="Omicron_BA2" & position==pos & mutant==BA2_wt,bind]
                               )
                          )
  }
}
entrench_BA2 <- entrench_BA2[-1,]

par(mfrow=c(6,3))
for(i in 1:nrow(entrench_BA2)){
  plot(1:4, 
       entrench_BA2[i,4:7],
       main=paste(entrench_BA2[i,c(2,1,3)],collapse=""),
       pch=16,cex=2,
       col=c(cbbPalette[1],cbbPalette[1],cbbPalette[8],cbbPalette[8]),
       xlim=c(0.75, 4.25), ylim=c(6,10.5),
       ylab="ACE2 affinity (-log10Kd)",
       xlab="",
       xaxt="n")
  axis(1, at=1:4, labels=c("WH1", "WH1+\nmut","BA.2+\nrevert","BA.2"))
  points(1:2, entrench_BA2[i,4:5],type="l",lwd=1.5)
  points(3:4, entrench_BA2[i,6:7],type="l",lwd=1.5,col=cbbPalette[8])
}

invisible(dev.print(pdf, paste(config$epistatic_shifts_dir,"/entrenchment_diagrams_BA2.pdf",sep=""),useDingbats=F))
```

<img src="epistatic_shifts_files/figure-gfm/entrenchment_Omicron_BA2-1.png" style="display: block; margin: auto;" />
