---
layout: epistasis
permalink: /epistatic-shifts/
---

---

### Overview 

You can use this tool to explore epistatic shifts in mutational effects on ACE2-binding affinity (-log10 $$K_D$$) between SARS-CoV-2 receptor-binding domain (RBD) variants. Once you've selected a comparison of interest, you can investigate the mutation-level epistatic shifts. 

#### Instructions

To use this tool, select two SARS-CoV-2 variants that you wish to compare between. To do this, simply select a *'comparator'* in the dropdown menu below the plot and select a *'variant'* by clicking on a variant name in the legend above the plot. Now you're comparing the epistatic shift (Jensen-Shannon divergence between ACE2 binding affinities) at each RBD site between the variant backgrounds you selected. 

Now that you've selected two variants to compare between, you can investigate site level differences by clicking on the points in the higlighted line plot. Simply click on a point – you should see the size of the point change indicating your selection – then you will see the differences in affinities of each of the 20 amino acids measured in each variant appear in the scatter plot on the right. Hover over individual mutations to see exact numerical details.

#### Technical Details

The epistatic shift is calculated as the Jensen-Shannon divergence in the set of Boltzmann-weighted affinities for all amino acids at each site. Mutation affinities were experimentally measured via high-throughput ACE2-binding titrations with yeast-displayed RBDs. The "Barcodes" quantity in per-mutation tooltips indicates the number of internally replicated barcodes with which a mutation was measured, where higher numbers indicate higher-confidence measurements. Mutations with fewer than 3 barcodes in either background are excluded from the Jensen-Shannon divergence calculation.

Data labeled "v1" (Alpha, Beta, Delta, Eta, and associated Wuhan-Hu-1 v1 measurements) are from a previously published study described [here](https://www.science.org/doi/10.1126/science.abo7896). The Wuhan-Hu-1 library was also spiked into the Omicron BA.1 and BA.2 libraries in a separate experiment. Although each Wuhan-Hu-1 dataset is closely correlated, it is most appropriate to compare each variant dataset to its matched Wuhan-Hu-1 dataset due to internal control.


Raw data (*NOTE: PRELIMINARY*) can be found [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_Omicron/blob/main/results/epistatic_shifts/JSD_versus_Wuhan1_by_target.csv) for a table of all pairwise RBD epistatic shifts, and [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_Omicron/blob/main/results/final_variant_scores/final_variant_scores.csv) for individual measurements of RBD mutant affinities. The code used to make these plots can be found [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_Omicron/blob/main/Epistatic-Shifts-Interactive-Visualization.ipynb). 
