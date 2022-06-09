# Count variants in each sample
This Python Jupyter notebook counts occurrences of each barcode in each sample from Illumina barcode sequencing, and adds these counts to the codon variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://plotnine.readthedocs.io/en/stable/) for ggplot2-like plotting.

The analysis relies heavily on the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package:


```python
import itertools
import multiprocessing
import multiprocessing.pool
import os
import warnings

import alignparse
import alignparse.targets

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.illuminabarcodeparser
import dms_variants.utils
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the gray-grid one defined in `dms_variants`:


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.2.4
    Using dms_variants version 0.8.9


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

### Parameters for notebook
Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

## Input variant tables
Initialize the table of barcode-variant pairs from the respective `process_ccs` notebooks for each background.


```python
variants = pd.read_csv(config['codon_variant_table_file_Wuhan_Hu_1'], na_filter=None)
variants = variants.append(pd.read_csv(config['codon_variant_table_file_BA1'], na_filter=None))
variants = variants.append(pd.read_csv(config['codon_variant_table_file_BA2'], na_filter=None))

variants = variants.reset_index(drop=True)

display(HTML(variants.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Wuhan_Hu_1</td>
      <td>pool1</td>
      <td>AAAAAAAAAAAGGAGA</td>
      <td>4</td>
      <td>GGT166ATG</td>
      <td>G166M</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan_Hu_1</td>
      <td>pool1</td>
      <td>AAAAAAAAAAATTTAA</td>
      <td>4</td>
      <td></td>
      <td></td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>Wuhan_Hu_1</td>
      <td>pool1</td>
      <td>AAAAAAAAAACGCGTA</td>
      <td>3</td>
      <td>GAA154ACT</td>
      <td>E154T</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan_Hu_1</td>
      <td>pool1</td>
      <td>AAAAAAAAAACTCCAA</td>
      <td>2</td>
      <td>TTT156ATG</td>
      <td>F156M</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan_Hu_1</td>
      <td>pool1</td>
      <td>AAAAAAAAACCGATTA</td>
      <td>2</td>
      <td>CAG84GAA</td>
      <td>Q84E</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


Are there any barcodes in the same library that are shared across targets?
If so, we need to get rid of those as they will be confounded in barcode parsing:


```python
dup_barcodes = (
    variants
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print('Here are duplicated barcodes:')
display(HTML(dup_barcodes.head().to_html(index=False)))

print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(variants)} barcodes:")
variants = (
    variants
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(variants)} barcodes.")
```

    Here are duplicated barcodes:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>duplicate_count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>AAAAAAACGTGTTAAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AACTAAAATCTAGTCT</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AACTATGAACCAATAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAGTAATAGAAAGACA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AATAAAAAAAATTTAC</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


    
    Removing the 51 duplicated barcodes.Started with 598496 barcodes:
    After removing duplicates, there are 598394 barcodes.


Pull out a target sequence for matching to the barcode and flanking sequence regions. Note, in this pipeline this is ok because our different backgrounds don't have differing flanks or other features within the actual N16 region covered in Illumina sequencing. If ever placing in-line barcodes here in the future, we would need to modify this.


```python
# get wildtype gene sequence for primary target
targets = alignparse.targets.Targets(seqsfile=config['amplicons_Wuhan_Hu_1'],
                                     feature_parse_specs=config['feature_parse_specs_Wuhan_Hu_1'])
```

## Setup to parse barcodes
Read data frame with list of all barcode runs.


```python
# barcode runs with R1 files by semicolon string split
barcode_runs = (pd.read_csv(config['barcode_runs'])
                .assign(R1=lambda x: x['R1'].str.split('; '))
                )

display(HTML(barcode_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>replicate</th>
      <th>sample</th>
      <th>sample_type</th>
      <th>sort_bin</th>
      <th>concentration</th>
      <th>date</th>
      <th>number_cells</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>220325</td>
      <td>1068797</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s01-b1_S1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>1.0</td>
      <td>220325</td>
      <td>1719961</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s01-b2_S2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>1.0</td>
      <td>220325</td>
      <td>2572303</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s01-b3_S3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>1.0</td>
      <td>220325</td>
      <td>7162967</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s01-b4_S4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_02_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>2.0</td>
      <td>220325</td>
      <td>423019</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s02-b1_S5_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_02_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>2.0</td>
      <td>220325</td>
      <td>2967228</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s02-b2_S6_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_02_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>2.0</td>
      <td>220325</td>
      <td>3774832</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s02-b3_S7_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_02_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>2.0</td>
      <td>220325</td>
      <td>5556904</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s02-b4_S8_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_03_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>3.0</td>
      <td>220325</td>
      <td>2678879</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s03-b1_S9_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_03_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>3.0</td>
      <td>220325</td>
      <td>2721170</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s03-b2_S10_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_03_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>3.0</td>
      <td>220325</td>
      <td>4196378</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s03-b3_S11_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_03_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>3.0</td>
      <td>220325</td>
      <td>3425525</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s03-b4_S12_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_04_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>4.0</td>
      <td>220325</td>
      <td>4923879</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s04-b1_S13_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_04_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>4.0</td>
      <td>220325</td>
      <td>5177366</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s04-b2_S14_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_04_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>4.0</td>
      <td>220325</td>
      <td>2427897</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s04-b3_S15_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_04_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>4.0</td>
      <td>220325</td>
      <td>73134</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s04-b4_S16_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_05_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>5.0</td>
      <td>220325</td>
      <td>9951655</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s05-b1_S17_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_05_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>5.0</td>
      <td>220325</td>
      <td>2592136</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s05-b2_S18_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_05_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>5.0</td>
      <td>220325</td>
      <td>36532</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s05-b3_S19_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_05_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>5.0</td>
      <td>220325</td>
      <td>5164</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s05-b4_S20_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_06_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>6.0</td>
      <td>220325</td>
      <td>11699823</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s06-b1_S21_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_06_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>6.0</td>
      <td>220325</td>
      <td>892312</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s06-b2_S22_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_06_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>6.0</td>
      <td>220325</td>
      <td>5571</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s06-b3_S23_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_06_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>6.0</td>
      <td>220325</td>
      <td>678</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s06-b4_S24_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_07_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>7.0</td>
      <td>220325</td>
      <td>11928098</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s07-b1_S25_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_07_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>7.0</td>
      <td>220325</td>
      <td>683081</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s07-b2_S26_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_07_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>7.0</td>
      <td>220325</td>
      <td>812</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s07-b3_S27_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_07_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>7.0</td>
      <td>220325</td>
      <td>672</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s07-b4_S28_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_08_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>8.0</td>
      <td>220325</td>
      <td>11925657</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s08-b1_S29_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_08_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>8.0</td>
      <td>220325</td>
      <td>731087</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s08-b2_S30_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_08_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>8.0</td>
      <td>220325</td>
      <td>783</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s08-b3_S31_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_08_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>8.0</td>
      <td>220325</td>
      <td>616</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s08-b4_S32_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_09_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>9.0</td>
      <td>220325</td>
      <td>12094579</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s09-b1_S33_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_09_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>9.0</td>
      <td>220325</td>
      <td>526503</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s09-b2_S34_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_09_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>9.0</td>
      <td>220325</td>
      <td>921</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s09-b3_S35_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_09_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>9.0</td>
      <td>220325</td>
      <td>604</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s09-b4_S36_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_01_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>220325</td>
      <td>1171478</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s10-b1_S37_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_01_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>1.0</td>
      <td>220325</td>
      <td>1573265</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s10-b2_S38_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_01_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>1.0</td>
      <td>220325</td>
      <td>2375280</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s10-b3_S39_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_01_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>1.0</td>
      <td>220325</td>
      <td>7676548</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s10-b4_S40_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_02_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>2.0</td>
      <td>220325</td>
      <td>420211</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s11-b1_S41_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220428_VH00699_111_AAATNWKM5/Unaligned/Project_tstarr/220325_s11-b1_S65_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_02_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>2.0</td>
      <td>220325</td>
      <td>2830684</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s11-b2_S42_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_02_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>2.0</td>
      <td>220325</td>
      <td>3442464</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s11-b3_S43_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_02_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>2.0</td>
      <td>220325</td>
      <td>5913534</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s11-b4_S44_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_03_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>3.0</td>
      <td>220325</td>
      <td>2641634</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s12-b1_S45_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_03_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>3.0</td>
      <td>220325</td>
      <td>2510488</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s12-b2_S46_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_03_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>3.0</td>
      <td>220325</td>
      <td>3939307</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s12-b3_S47_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_03_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>3.0</td>
      <td>220325</td>
      <td>3500482</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s12-b4_S48_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_04_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>4.0</td>
      <td>220325</td>
      <td>5412146</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s13-b1_S49_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_04_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>4.0</td>
      <td>220325</td>
      <td>5328949</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s13-b2_S50_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_04_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>4.0</td>
      <td>220325</td>
      <td>2570653</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s13-b3_S51_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_04_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>4.0</td>
      <td>220325</td>
      <td>72349</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s13-b4_S52_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_05_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>5.0</td>
      <td>220325</td>
      <td>10391185</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s14-b1_S53_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_05_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>5.0</td>
      <td>220325</td>
      <td>2179850</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s14-b2_S54_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_05_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>5.0</td>
      <td>220325</td>
      <td>28870</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s14-b3_S55_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220428_VH00699_111_AAATNWKM5/Unaligned/Project_tstarr/220325_s14-b3_S66_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_05_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>5.0</td>
      <td>220325</td>
      <td>4662</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s14-b4_S56_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_06_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>6.0</td>
      <td>220325</td>
      <td>11816796</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s15-b1_S57_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_06_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>6.0</td>
      <td>220325</td>
      <td>753443</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s15-b2_S58_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_06_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>6.0</td>
      <td>220325</td>
      <td>4710</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s15-b3_S59_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_06_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>6.0</td>
      <td>220325</td>
      <td>393</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s15-b4_S60_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_07_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>7.0</td>
      <td>220325</td>
      <td>11956441</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s16-b1_S61_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_07_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>7.0</td>
      <td>220325</td>
      <td>654607</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s16-b2_S62_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_07_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>7.0</td>
      <td>220325</td>
      <td>574</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s16-b3_S63_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_07_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>7.0</td>
      <td>220325</td>
      <td>366</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s16-b4_S64_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_08_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>8.0</td>
      <td>220325</td>
      <td>11858829</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s17-b1_S65_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_08_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>8.0</td>
      <td>220325</td>
      <td>721997</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s17-b2_S66_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_08_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>8.0</td>
      <td>220325</td>
      <td>621</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s17-b3_S67_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_08_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>8.0</td>
      <td>220325</td>
      <td>397</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s17-b4_S68_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_09_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>9.0</td>
      <td>220325</td>
      <td>12077783</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s18-b1_S69_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_09_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>9.0</td>
      <td>220325</td>
      <td>512747</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s18-b2_S70_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_09_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>9.0</td>
      <td>220325</td>
      <td>637</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s18-b3_S71_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>TiteSeq_09_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>9.0</td>
      <td>220325</td>
      <td>418</td>
      <td>[/shared/ngs/illumina/tstarr/220415_VH00699_102_AAAYMVFM5/Unaligned/Project_tstarr/220325_s18-b4_S72_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>220413</td>
      <td>609600</td>
      <td>[/shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib12-22-24-26_bin1_1_S69_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib12-22-24-26_bin1_2_S70_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib12-22-24-26_bin1_3_S71_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>220413</td>
      <td>4137873</td>
      <td>[/shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib12-22-24-26_bin2_1_S72_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib12-22-24-26_bin2_2_S73_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib12-22-24-26_bin2_3_S74_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>220413</td>
      <td>6652066</td>
      <td>[/shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib12-22-24-26_bin3_1_S75_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib12-22-24-26_bin3_2_S76_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib12-22-24-26_bin3_3_S77_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>A</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>220413</td>
      <td>6307756</td>
      <td>[/shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib12-22-24-26_bin4_1_S78_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib12-22-24-26_bin4_2_S79_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib12-22-24-26_bin4_3_S80_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>220413</td>
      <td>648890</td>
      <td>[/shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib13-23-25_bin1_1_S81_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib13-23-25_bin1_2_S82_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib13-23-25_bin1_3_S83_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>220413</td>
      <td>3823934</td>
      <td>[/shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib13-23-25_bin2_1_S84_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib13-23-25_bin2_2_S85_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib13-23-25_bin2_3_S86_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>220413</td>
      <td>6712514</td>
      <td>[/shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib13-23-25_bin3_1_S87_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib13-23-25_bin3_2_S88_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib13-23-25_bin3_3_S89_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>A</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>220413</td>
      <td>6284677</td>
      <td>[/shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib13-23-25_bin4_1_S90_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib13-23-25_bin4_2_S91_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/220427_VH00319_200_AAATTVLM5/Unaligned/Project_tstarr/220413_lib13-23-25_bin4_3_S92_R1_001.fastq.gz]</td>
    </tr>
  </tbody>
</table>


Make sure library / sample combinations are unique:


```python
assert len(barcode_runs) == len(barcode_runs.groupby(['library', 'replicate', 'sample']))
```

Make sure the the libraries for which we have barcode runs are all in our variant table:


```python
unknown_libs = set(barcode_runs['library']) - set(variants['library'])
if unknown_libs:
    raise ValueError(f"Libraries with barcode runs not in variant table: {unknown_libs}")
```

Now we initialize an [IlluminaBarcodeParser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) for each library.

First, get the length of the barcode from the alignment target after making sure the same length for all targets:


```python
bclen = len(targets.targets[0].get_feature('barcode').seq)

assert (bclen == len(target.get_feature('barcode').seq) for target in targets.targets)

print(f"Barcodes of length {bclen}")
```

    Barcodes of length 16


The other barcode parsing params come from the config file:


```python
parser_params = config['illumina_barcode_parser_params']

display(HTML(
    pd.Series(parser_params, name='value')
    .rename_axis(index='parameter')
    .reset_index()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>parameter</th>
      <th>value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>upstream</td>
      <td>GGCCGC</td>
    </tr>
    <tr>
      <td>downstream</td>
      <td></td>
    </tr>
    <tr>
      <td>minq</td>
      <td>20</td>
    </tr>
    <tr>
      <td>upstream_mismatch</td>
      <td>1</td>
    </tr>
    <tr>
      <td>downstream_mismatch</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


The parser needs to know the set of valid barcodes, which are stored in the variant table and are different for each library.
So we create a different parser for each library using these valid barcode sets:


```python
# create dict keyed by library, value is parser for library
parsers = {lib: dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    bclen=bclen,
                    valid_barcodes=variants.loc[variants['library']==lib]['barcode'],
                    **parser_params)
           for lib in set(variants['library'])}

print('Number of valid barcodes searched for by each parser:')
display(HTML(
    pd.DataFrame([(lib, len(p.valid_barcodes)) for lib, p in parsers.items()],
                 columns=['library', 'number of valid barcodes'])
    .to_html(index=False)
    ))
```

    Number of valid barcodes searched for by each parser:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>number of valid barcodes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool2</td>
      <td>302948</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>295446</td>
    </tr>
  </tbody>
</table>


## Parse barcodes
We now parse the barcodes.
Since this will take a while, we utilize multiple CPUs via the Python [multiprocessing](https://docs.python.org/3.6/library/multiprocessing.html) module.
First, determine how many CPUs to use.
We use the minimum of the user-specified number hardcoded below and the number actually available.
(If you are running *interactively* on the Hutch cluster, you may need to reduce the number below in order to avoid an error as there is an enforced CPU limit on the home `rhino` nodes):


```python
ncpus = min(config['max_cpus'], multiprocessing.cpu_count())
print(f"Using {ncpus} CPUs")
```

    Using 8 CPUs


Parse the barcodes in parallel via a [multiprocessing.Pool](https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.pool.Pool) using all the available CPUs to get a list of the data frames with barcode counts / fates for each sample:


```python
def process_func(parser, r1files, library, replicate, sample):
    """Convenience function to be starmapped to multiprocessing pool."""
    return parser.parse(r1files, add_cols={'library': library, 'replicate': replicate, 'sample': sample})

# parallel computation of list of data frames
with multiprocessing.pool.Pool(processes=ncpus) as pool:
    bclist = pool.starmap(
                process_func,
                [(parsers[run.library], run.R1, run.library, run.replicate, run.sample)
                  for run in barcode_runs.itertuples()],
                )
```

Now concatenate the list into data frames of barcode counts and barcode fates:


```python
counts = pd.concat([samplecounts for samplecounts, _ in bclist],
                   sort=False,
                   ignore_index=True)

print('First few lines of counts data frame:')
display(HTML(counts.head().to_html(index=False)))

fates = pd.concat([samplefates for _, samplefates in bclist],
                  sort=False,
                  ignore_index=True)

print('First few lines of fates data frame:')
display(HTML(fates.head().to_html(index=False)))
```

    First few lines of counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>barcode</th>
      <th>count</th>
      <th>library</th>
      <th>replicate</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>CAAGTTTATGGGGCAA</td>
      <td>654</td>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>GGAGATCAAATGATAG</td>
      <td>633</td>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>CAAGGATGTGTTTTAT</td>
      <td>631</td>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>ACGAGAGCCGCAGTAC</td>
      <td>621</td>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>TCTGAATGGTGGACTT</td>
      <td>562</td>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
  </tbody>
</table>


    First few lines of fates data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>fate</th>
      <th>count</th>
      <th>library</th>
      <th>replicate</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>valid barcode</td>
      <td>1186902</td>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>invalid barcode</td>
      <td>229075</td>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>low quality barcode</td>
      <td>222673</td>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>unparseable barcode</td>
      <td>19288</td>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>failed chastity filter</td>
      <td>0</td>
      <td>pool1</td>
      <td>A</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
  </tbody>
</table>


## Examine fates of parsed barcodes
First, we'll analyze the "fates" of the parsed barcodes.
These fates represent what happened to each Illumina read we parsed:
 - Did the barcode read fail the Illumina chastity filter?
 - Was the barcode *unparseable* (i.e., the read didn't appear to be a valid barcode based on flanking regions)?
 - Was the barcode sequence too *low quality* based on the Illumina quality scores?
 - Was the barcode parseable but *invalid* (i.e., not in our list of variant-associated barcodes in the codon variant table)?
 - Was the barcode *valid*, and so will be added to variant counts.
 
First, we just write a CSV file with all the barcode fates:


```python
fatesfile = os.path.join(config['counts_dir'], 'barcode_fates.csv')
print(f"Writing barcode fates to {fatesfile}")
fates.to_csv(fatesfile, index=False)
```

    Writing barcode fates to results/counts/barcode_fates.csv


Next, we tabulate the barcode fates in wide format:


```python
display(HTML(fates
             .pivot_table(columns='fate',
                          values='count',
                          index=['library', 'replicate', 'sample'])
             .to_html()
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>fate</th>
      <th>failed chastity filter</th>
      <th>invalid barcode</th>
      <th>low quality barcode</th>
      <th>unparseable barcode</th>
      <th>valid barcode</th>
    </tr>
    <tr>
      <th>library</th>
      <th>replicate</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="40" valign="top">pool1</th>
      <th rowspan="40" valign="top">A</th>
      <th>SortSeq_bin1</th>
      <td>0</td>
      <td>391891</td>
      <td>230548</td>
      <td>17653</td>
      <td>756376</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>0</td>
      <td>1291220</td>
      <td>1747805</td>
      <td>125122</td>
      <td>7991800</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>0</td>
      <td>2539923</td>
      <td>3880131</td>
      <td>262335</td>
      <td>17601390</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>0</td>
      <td>3253049</td>
      <td>4043956</td>
      <td>277282</td>
      <td>17681950</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin1</th>
      <td>0</td>
      <td>229075</td>
      <td>222673</td>
      <td>19288</td>
      <td>1186902</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin2</th>
      <td>0</td>
      <td>266511</td>
      <td>318459</td>
      <td>27488</td>
      <td>1733955</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin3</th>
      <td>0</td>
      <td>527451</td>
      <td>664751</td>
      <td>52569</td>
      <td>3527583</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin4</th>
      <td>0</td>
      <td>1857881</td>
      <td>1963186</td>
      <td>158366</td>
      <td>9968694</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin1</th>
      <td>0</td>
      <td>83361</td>
      <td>63227</td>
      <td>4761</td>
      <td>315377</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin2</th>
      <td>0</td>
      <td>773683</td>
      <td>709452</td>
      <td>58667</td>
      <td>3675583</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin3</th>
      <td>0</td>
      <td>908260</td>
      <td>965335</td>
      <td>78760</td>
      <td>5166563</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin4</th>
      <td>0</td>
      <td>1830091</td>
      <td>1624485</td>
      <td>138211</td>
      <td>8419416</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin1</th>
      <td>0</td>
      <td>865873</td>
      <td>813794</td>
      <td>65849</td>
      <td>4264964</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin2</th>
      <td>0</td>
      <td>741541</td>
      <td>760058</td>
      <td>63334</td>
      <td>3932089</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin3</th>
      <td>0</td>
      <td>1309536</td>
      <td>1180160</td>
      <td>94570</td>
      <td>5989912</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin4</th>
      <td>0</td>
      <td>1026755</td>
      <td>907480</td>
      <td>76328</td>
      <td>4656335</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin1</th>
      <td>0</td>
      <td>1568522</td>
      <td>1450343</td>
      <td>110885</td>
      <td>6803827</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin2</th>
      <td>0</td>
      <td>1658188</td>
      <td>1455265</td>
      <td>118071</td>
      <td>7200930</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin3</th>
      <td>0</td>
      <td>616488</td>
      <td>638250</td>
      <td>53776</td>
      <td>3294957</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin4</th>
      <td>0</td>
      <td>7558</td>
      <td>10347</td>
      <td>1095</td>
      <td>51405</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin1</th>
      <td>0</td>
      <td>2659504</td>
      <td>2597809</td>
      <td>211719</td>
      <td>13164748</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin2</th>
      <td>0</td>
      <td>701112</td>
      <td>759502</td>
      <td>63053</td>
      <td>3947229</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin3</th>
      <td>0</td>
      <td>4223</td>
      <td>5986</td>
      <td>546</td>
      <td>27543</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin4</th>
      <td>0</td>
      <td>1334</td>
      <td>5693</td>
      <td>105</td>
      <td>6789</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin1</th>
      <td>0</td>
      <td>3018184</td>
      <td>3275350</td>
      <td>266295</td>
      <td>16769263</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin2</th>
      <td>0</td>
      <td>180111</td>
      <td>197676</td>
      <td>16420</td>
      <td>1027432</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin3</th>
      <td>0</td>
      <td>322</td>
      <td>2864</td>
      <td>36</td>
      <td>1839</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin4</th>
      <td>0</td>
      <td>825</td>
      <td>933</td>
      <td>66</td>
      <td>4441</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin1</th>
      <td>0</td>
      <td>2718543</td>
      <td>2936944</td>
      <td>234953</td>
      <td>14994909</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin2</th>
      <td>0</td>
      <td>176919</td>
      <td>190654</td>
      <td>15368</td>
      <td>997128</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin3</th>
      <td>0</td>
      <td>232</td>
      <td>2414</td>
      <td>19</td>
      <td>975</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin4</th>
      <td>0</td>
      <td>88</td>
      <td>201</td>
      <td>3</td>
      <td>237</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin1</th>
      <td>0</td>
      <td>3028052</td>
      <td>3194912</td>
      <td>258529</td>
      <td>16697516</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin2</th>
      <td>0</td>
      <td>177248</td>
      <td>197755</td>
      <td>15956</td>
      <td>1000227</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin3</th>
      <td>0</td>
      <td>81</td>
      <td>683</td>
      <td>13</td>
      <td>243</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin4</th>
      <td>0</td>
      <td>239</td>
      <td>401</td>
      <td>6</td>
      <td>182</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin1</th>
      <td>0</td>
      <td>3225156</td>
      <td>3471837</td>
      <td>285726</td>
      <td>17922207</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin2</th>
      <td>0</td>
      <td>113845</td>
      <td>123887</td>
      <td>10382</td>
      <td>647189</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin3</th>
      <td>0</td>
      <td>73</td>
      <td>3272</td>
      <td>65</td>
      <td>282</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin4</th>
      <td>0</td>
      <td>76</td>
      <td>921</td>
      <td>6</td>
      <td>129</td>
    </tr>
    <tr>
      <th rowspan="40" valign="top">pool2</th>
      <th rowspan="40" valign="top">A</th>
      <th>SortSeq_bin1</th>
      <td>0</td>
      <td>429957</td>
      <td>257609</td>
      <td>17752</td>
      <td>910440</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>0</td>
      <td>1092897</td>
      <td>1629576</td>
      <td>115773</td>
      <td>7467403</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>0</td>
      <td>2283737</td>
      <td>3679347</td>
      <td>265883</td>
      <td>16998195</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>0</td>
      <td>2555512</td>
      <td>3600988</td>
      <td>245622</td>
      <td>16337991</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin1</th>
      <td>0</td>
      <td>284774</td>
      <td>316982</td>
      <td>25563</td>
      <td>1689174</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin2</th>
      <td>0</td>
      <td>266794</td>
      <td>348335</td>
      <td>29295</td>
      <td>1898341</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin3</th>
      <td>0</td>
      <td>383113</td>
      <td>498491</td>
      <td>42391</td>
      <td>2698429</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin4</th>
      <td>0</td>
      <td>1727860</td>
      <td>1932410</td>
      <td>155028</td>
      <td>10259354</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin1</th>
      <td>0</td>
      <td>162401</td>
      <td>149089</td>
      <td>10366</td>
      <td>751854</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin2</th>
      <td>0</td>
      <td>560700</td>
      <td>652061</td>
      <td>53302</td>
      <td>3446404</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin3</th>
      <td>0</td>
      <td>582773</td>
      <td>748394</td>
      <td>63233</td>
      <td>4032171</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin4</th>
      <td>0</td>
      <td>1462312</td>
      <td>1647911</td>
      <td>131334</td>
      <td>8534148</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin1</th>
      <td>0</td>
      <td>591526</td>
      <td>686006</td>
      <td>56241</td>
      <td>3656585</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin2</th>
      <td>0</td>
      <td>479785</td>
      <td>616020</td>
      <td>48876</td>
      <td>3226485</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin3</th>
      <td>0</td>
      <td>846222</td>
      <td>972457</td>
      <td>80602</td>
      <td>5235165</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin4</th>
      <td>0</td>
      <td>858453</td>
      <td>972663</td>
      <td>77588</td>
      <td>5128906</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin1</th>
      <td>0</td>
      <td>1309378</td>
      <td>1539502</td>
      <td>125979</td>
      <td>8092055</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin2</th>
      <td>0</td>
      <td>1235624</td>
      <td>1419216</td>
      <td>114601</td>
      <td>7389375</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin3</th>
      <td>0</td>
      <td>547515</td>
      <td>640520</td>
      <td>51205</td>
      <td>3422885</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin4</th>
      <td>0</td>
      <td>9077</td>
      <td>13206</td>
      <td>1349</td>
      <td>69308</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin1</th>
      <td>0</td>
      <td>2326857</td>
      <td>2557162</td>
      <td>208823</td>
      <td>13253764</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin2</th>
      <td>0</td>
      <td>500445</td>
      <td>583659</td>
      <td>46942</td>
      <td>3094947</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin3</th>
      <td>0</td>
      <td>51030</td>
      <td>77093</td>
      <td>8958</td>
      <td>414328</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin4</th>
      <td>0</td>
      <td>67</td>
      <td>185</td>
      <td>3</td>
      <td>392</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin1</th>
      <td>0</td>
      <td>2798842</td>
      <td>3087605</td>
      <td>248821</td>
      <td>16105023</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin2</th>
      <td>0</td>
      <td>149463</td>
      <td>182696</td>
      <td>15219</td>
      <td>968959</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin3</th>
      <td>0</td>
      <td>296</td>
      <td>913</td>
      <td>115</td>
      <td>1863</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin4</th>
      <td>0</td>
      <td>12</td>
      <td>131</td>
      <td>3</td>
      <td>21</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin1</th>
      <td>0</td>
      <td>3026157</td>
      <td>3543789</td>
      <td>284853</td>
      <td>18351351</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin2</th>
      <td>0</td>
      <td>133055</td>
      <td>158292</td>
      <td>12562</td>
      <td>834137</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin3</th>
      <td>0</td>
      <td>28</td>
      <td>421</td>
      <td>2</td>
      <td>78</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin4</th>
      <td>0</td>
      <td>5</td>
      <td>22</td>
      <td>0</td>
      <td>28</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin1</th>
      <td>0</td>
      <td>2992343</td>
      <td>3447493</td>
      <td>274694</td>
      <td>18013558</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin2</th>
      <td>0</td>
      <td>181611</td>
      <td>214662</td>
      <td>18031</td>
      <td>1140351</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin3</th>
      <td>0</td>
      <td>48</td>
      <td>302</td>
      <td>81</td>
      <td>132</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin4</th>
      <td>0</td>
      <td>12</td>
      <td>138</td>
      <td>0</td>
      <td>38</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin1</th>
      <td>0</td>
      <td>2796299</td>
      <td>3823618</td>
      <td>271588</td>
      <td>17354890</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin2</th>
      <td>0</td>
      <td>91194</td>
      <td>106228</td>
      <td>8958</td>
      <td>575057</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin3</th>
      <td>0</td>
      <td>1991</td>
      <td>663</td>
      <td>38</td>
      <td>214</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin4</th>
      <td>0</td>
      <td>3</td>
      <td>259</td>
      <td>13</td>
      <td>11</td>
    </tr>
  </tbody>
</table>


Now we plot the barcode-read fates for each library / sample, showing the bars for valid barcodes in orange and the others in gray.
We see that the largest fraction of barcode reads correspond to valid barcodes, and most of the others are invalid barcodes (probably because the map to variants that aren't present in our variant table since we didn't associate all variants with barcodes). The exception to this is lib2 Titeseq_03_bin3; the PCR for this sample in the original sequencing run failed, so we followed it up with a single MiSeq lane. We did not filter out the PhiX reads from this data before parsing, so these PhiX reads will deflate the fraction of valid barcode reads as expected, but does not indicate any problems.


```python
barcode_fate_plot = (
    ggplot(
        fates
        .assign(sample=lambda x: pd.Categorical(x['sample'],
                                                x['sample'].unique(),
                                                ordered=True),
                fate=lambda x: pd.Categorical(x['fate'],
                                              x['fate'].unique(),
                                              ordered=True),
                is_valid=lambda x: x['fate'] == 'valid barcode'
                ), 
        aes('fate', 'count', fill='is_valid')) +
    geom_bar(stat='identity') +
    facet_grid('sample ~ library + replicate') +
    facet_grid('sample ~ library + replicate') +
    scale_fill_manual(CBPALETTE, guide=False) +
    theme(figure_size=(1.4 * (1 + fates['library'].nunique() + fates['replicate'].nunique()),
                       1.7 * (1.2 + fates['sample'].nunique())),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank()
          ) +
    scale_y_continuous(labels=dms_variants.utils.latex_sci_not,
                       name='number of reads')
    )

_ = barcode_fate_plot.draw()
```


    
![png](count_variants_files/count_variants_42_0.png)
    


## Output csv of barcode counts in variant-barcode lookup table


```python
print(f"Writing variant counts to {config['variant_counts_file']}")
counts.to_csv(config['variant_counts_file'], index=False)
```

    Writing variant counts to results/counts/variant_counts.csv


The [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) has lots of nice functions that can be used to analyze the counts it contains.
However, we do that in the next notebook so we don't have to re-run this entire (rather computationally intensive) notebook every time we want to analyze a new aspect of the counts.


```python

```
