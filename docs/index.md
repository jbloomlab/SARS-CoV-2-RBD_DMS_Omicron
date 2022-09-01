---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: default
---

## Overview 

We performed deep mutational scans to measure the impact on ACE2 binding of all single amino-acid mutations in the Omicron BA.1 and BA.2 variant RBDs, and compare mutational effects to those in the ancestral Wuhan-Hu-1 RBD background. We combine this with other variant RBD backgrounds described in this [prior publication](https://www.science.org/doi/10.1126/science.abo7896).

Here, we link to two interactive visualizations that enable exploration of the data from these deep mutational scans across SARS_CoV-2 variants.


### Instructions 

We have made two tools to help visualize the data from our deep mutational scans:

1. A set of interactive heatmaps that you can use to explore the change in ACE2-binding affinity ($$\Delta$$log10 $$K_D$$) or the change in RBD expression (log10(MFI)) caused by mutations in each RBD variant. To use this tool, click [here]({{ site.baseurl }}{% link heatmaps.md %}).

2. An interactive widget that you can use to visualize epistatic shifts in mutational effects on ACE2-binding affinity (-log10 $$K_D$$) between variant RBDs. To use this tool, click [here]({{ site.baseurl }}{% link epistasis.md %}).  

### Data

If you are interested in the raw data from our study, you can find the ACE2-binding affinity (-log10 $$K_D$$) and RBD expression (log10(MFI)) measurements from each experiment [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_Omicron/tree/main/results/final_variant_scores). You can find the data used to plot the epistatic shifts between variant backgrounds [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_Omicron/blob/main/results/epistatic_shifts/JSD_versus_Wuhan1_by_target.csv). Data for pre-Omicron variants (Alpha, Beta, Delta, Eta, and associated Wuhan-Hu-1 measurements) are described in [this paper](https://www.science.org/doi/10.1126/science.abo7896) and found on [this repository](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_variants).