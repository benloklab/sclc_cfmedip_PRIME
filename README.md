# SCLC cfMeDIP PRIME
This repository describes the Peripheral blood leukocyte methylation (PRIME) subtraction algorithm used in the following publication: Cell-free DNA methylation-defined prognostic subgroups in small-cell lung cancer identified by leukocyte methylation subtraction
doi: 10.1016/j.isci.2022.105487. Code written by Sami Ul Haq, MSc.

In brief, to increase specificity of the cancer signal in the cell-free DNA, i.e. in our using small cell lung cancer (SCLC) liquid biopsy samples, compared to non-cancer noise ratio, we implemented a novel approach utilizing *paired* peripheral blood leukocyte (PBL) genomic DNA collected from the same plasma source material of patient with SCLC at identical timepoints. Comparison of total plasma cell-free DNA to PBL genomic DNA by principal component analysis revealed that PBLs exhibited a distinct methylation signal. We implemented our novel algorithmic filter, PRIME, to reduce the PBL methylation signals. We started with whole-genome windows (n = 9,603,445), removed ENCODE-blacklist regions, and then selected windows hypomethylated across PBLs (median beta per window < 0.2). Within these PBL-hypomethylated windows, we further selected for windows with a CG-density threshold >= 5 to account for the functionality of the 5-methylcytosine antibody. Thus, PRIME filtered out non-tumor noise in the cfDNA, reducing 9,603,445 whole-genome windows to 196,582 SCLC-specific windows.

While we utilized a specific criterion (median beta per window < 0.2 and CG-density >= 5) to select high-confidence SCLC-specific windows, our approach can be modified to fit looser filtering criterion.

## Software requirements
PRIME was developed and tested using R version 3.6.1 on a high-performance computing cluster with a minimum memory of 60GB. 

## Script order
The R scripts in this repository should be run in the following sequence:

- step0.medestrand.converter.R
- step1.data.splitter.R
- step2.CpG.per.window.R
- step3.PBL.depleted.windows.finder.R
- step4.selecting.PRIME.windows.R

## Running the PRIME pipeline
Each script should be run on an HPC node with **~24-48hrs** time allocated and **60GB of memory** in the following manner:

The first step involves converting cfMeDIP-seq counts data to EPIC/450K-style beta-values
```bash
Rscript step0.medestrand.converter.R DIRECTORY_WITH_BAM_FILES
```

The next step is breaking up the massive chr1-22 genome-wide 300bp matrix (obtained in the previous step) into smaller chunks to ease computation
```bash
Rscript step1.data.splitter.R
```

Subsequently, the CpG density for each 300bp window needs to be caculated. Here the smaller chunks from the previous steps need to be passed 
```bash
for EACHFILE in $(ls PBL.depleted.windows.part*.RData)
do
sbatch step2.CpG.per.window.R.sh $EACHFILE
done
```

Then, the windows in PBLs with beta-values cutoffs are identified and saved
```bash
Rscript step3.PBL.depleted.windows.finder.R
```

Lastly, the CpG density data and the beta-value cutoffs are both examined and the final PRIME windows selected
```bash
Rscript step4.selecting.PRIME.windows.R
```
