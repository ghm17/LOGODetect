# LOGODetect
LOGODetect (LOcal Genetic cOrrelation Dectector) is a powerful tool to identify small segments that harbor local genetic correlation between two traits. We have now updated the software, which can identify small regions with significant local genetic correlation across two populations.

# Before starting
* `LOGODetect` is developed using R and tested in Linux environments. The statistical computing software R (>=3.5.1) and the following R packages are required:
  * [snowfall](https://CRAN.R-project.org/package=snowfall)
  * [data.table](https://CRAN.R-project.org/package=data.table) (>=1.11.8)
  * [optparse](https://CRAN.R-project.org/package=optparse) (>=1.6.6)
  * [BEDMatrix](https://CRAN.R-project.org/package=BEDMatrix) (>=2.0.3)
  * [XPASS](https://github.com/YangLabHKUST/XPASS)

* `LOGODetect` requires `ldsc` to calculate stratified genetic covariance, installation steps see [here](https://github.com/bulik/ldsc). Make sure that `ldsc` goes well after command `conda activate ldsc` is run. 

# Single-population cross-trait analysis

`LOGODetect` requires the reference genotype data and the pre-computed LD score. Here are the command line to download these reference:
```bash
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/LOGODetect/LOGODetect_data.tar.gz
tar -zxvf LOGODetect_data.tar.gz
```

### Applying LOGODetect

```bash

conda activate ldsc

Rscript /LOGODetect.R \
--sumstats PATH_TO_SUMSTAT1,PATH_TO_SUMSTAT2 \
--n_gwas N1,N2 \
--ref_dir PATH_TO_REFERENCE \
--pop POPULATION \
--ldsc_dir PATH_TO_LDSC \
--block_partition PATH_TO_GENOME_PARTITION \
--out_dir PATH_TO_OUTFILE \
# The following flags are optional.
--max_nsnps CN \
--interval INTER \
--chr CHR \
--n_cores N_CORE

conda deactivate

```
where the inputs in order are

* `PATH_TO_SUMSTAT1` and `PATH_TO_SUMSTAT2` (required): Full paths to two GWAS summary statistics (in hg19 genome build), separated by comma.
```
    CHR    SNP           BP        A1    A2    BETA       P
    1      rs4074137     1026707   A     C     0.9942     0.7198
    1      rs11260590    1027070   A     G     1.0108     0.6765
    1      rs74048006    1027845   T     G     0.9909     0.6547
    ...
```
Or:
```
    CHR    SNP           BP        A1    A2    OR        P
    1      rs4074137     1026707   C     A     0.9985    0.8976
    1      rs11260590    1027070   G     A     1.0049    0.7808
    1      rs74048006    1027845   T     G     0.9993    0.9579
    ...
```
* `N1` and `N2` (required): Sample sizes of two GWAS summary statistics, in the same order of the GWAS summary statistics, separated by comma.
* `PATH_TO_REFERENCE` (required): Full path to the directory that contains the reference genotype data and the pre-computed LD score. 
* `POPULATION` (required): Population for GWAS sample, currently EUR is allowed. 
* `PATH_TO_LDSC` (required): Full path to the directory that contains the script of LDSC.
* `PATH_TO_GENOME_PARTITION` (required): Full path to the genome partition file. Sample data in `./LOGODetect/block_partition.txt` is provided.
* `PATH_TO_OUTFILE` (required): Full directory to the output regions (and the temporary files). 
* `CN` (optional): Maximal number of SNPs in a true signal region. Default is 2000 for GWAS with about 5000000 SNPs. 
* `INTER` (optional): Number of SNPs in an output region will be an integer multiplied by INTER. Default is 20 for GWAS with about 5000000 SNPs. When analyzing GWAS summary statistics with reduced SNP density, the user can re-specify the parameters CN and INTER to save computation time. For example, for analysis restricted to HapMap3 SNPs with about 1000000 SNPs, the user can set CN as 500 and INTER as 5. 
* `CHR` (optional): The chromosome number, perform chromosome-specific analysis if provided.
* `N_CORE` (optional): Number of cores used in parallel computing. Default is 20. 

### A concrete example
```bash
#!/bin/bash
## ---- Download the required reference panel and example data for LOGODetect ---- ##
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/LOGODetect/LOGODetect_data.tar.gz
tar -zxvf LOGODetect_data.tar.gz
rm -rf LOGODetect_data.tar.gz


## ---- Applying LOGODetect ---- ##
cd LOGODetect_data

mkdir ./results

conda activate ldsc

Rscript /LOGODetect/LOGODetect.R \
--sumstats ./sumstats/BIP.txt,./sumstats/SCZ.txt \
--n_gwas 51710,105318 \
--ref_dir ./LOGODetect_1kg_ref \
--pop EUR \
--ldsc_dir /LOGODetect/ldsc \
--block_partition /LOGODetect/block_partition.txt \
--out_dir ./results \
--n_cores 25

conda deactivate

```

### Output
LOGODetect outputs a whitespace-delimited text file `LOGODetect_regions.txt` in `PATH_TO_OUTFILE` specified by the user, with each row representing one small segment and the columns as such:
* `chr`: The chromosome. 
* `begin_pos`: The starting position of this detected small region.
* `stop_pos`: The stopping position of this detected small region.
* `stat`: The scan statistic value. Positive value means positive local genetic correlation between two traits. 
* `pval`: The p-value of this detected small region.
* `pval_adj`: The adjusted p-value of this detected small region.

We have prepared the example output file for you in `/LOGODetect_data/results/LOGODetect_regions.txt`. 


# Cross-population analysis

* A detailed description can be seen [here](https://github.com/qlu-lab/X-Wing). Below shows the short tutorial. `LOGODetect` requires the reference genotype data in [plink](https://zzz.bwh.harvard.edu/plink/) format and the pre-computed LD score using [LDSC](https://github.com/bulik/ldsc). We have provided these data for EUR, EAS, AFR, SAS, and AMR population using 1000 Genomes Project phase 3 samples. Here are the command line to download these reference:
  * AFR reference:
    ```bash
    wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LOGODetect/LOGODetect_1kg_AFR.tar.gz
    tar -zxvf LOGODetect_1kg_AFR.tar.gz
    ```
  * AMR reference:
    ```bash
    wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LOGODetect/LOGODetect_1kg_AMR.tar.gz
    tar -zxvf LOGODetect_1kg_AMR.tar.gz
    ```
  * EAS reference:
    ```bash
    wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LOGODetect/LOGODetect_1kg_EAS.tar.gz
    tar -zxvf LOGODetect_1kg_EAS.tar.gz
    ```
  * EUR reference:
    ```bash
    wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LOGODetect/LOGODetect_1kg_EUR.tar.gz
    tar -zxvf LOGODetect_1kg_EUR.tar.gz
    ```
  * SAS reference:
    ```bash
    wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LOGODetect/LOGODetect_1kg_SAS.tar.gz
    tar -zxvf LOGODetect_1kg_SAS.tar.gz

### Applying LOGODetect

```bash
Rscript LOGODetect.R \
--sumstats PATH_TO_SUMSTAT1,PATH_TO_SUMSTAT2 \
--n_gwas N1,N2 \
--ref_dir PATH_TO_REFERENCE \
--pop POPULATION1,POPULATION2 \
--block_partition PATH_TO_GENOME_PARTITION \
--gc_snp PATH_TO_SNPLIST \
--out_dir PATH_TO_OUTFILE \
# The following flags are optional.
--n_cores N_CORE
```
where the inputs in order are

* `PATH_TO_SUMSTAT1` and `PATH_TO_SUMSTAT2` (required): Full paths to two GWAS summary statistics, separated by comma.
* `N1` and `N2` (required): Sample sizes of two GWAS summary statistics, in the same order of the GWAS summary statistics, separated by comma.
* `PATH_TO_REFERENCE` (required): Full path to the directory that contains the reference genotype data in the in the binary `plink` `.bed/.bim/.fam` format, the covariate files (e.g. principal components) for the reference panel, and the pre-computed LD scores. 
* `POPULATION1` and `POPULATION2` (required): Populations for GWAS sample (could be the same in single-population cross-trait analysis), in the same order of the GWAS summary statistics, separated by comma. AFR, AMR, EAS, EUR, and SAS are allowed. 
* `PATH_TO_GENOME_PARTITION` (required): Full path to the genome partition file.
* `PATH_TO_SNPLIST` (required): Full path to HapMap3 SNPs file. Sample data in `./LOGODetect/1kg_hm3_snp.txt` is provided.
* `PATH_TO_OUTFILE` (required): Full directory to the output regions. 
* `N_CORE` (optional): Number of cores used in parallel computing. Default is 20. 

### A concrete example
```bash
#!/bin/bash
## ---- Download the example data ---- ##
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/example/X-Wing_example.tar.gz
tar -zxvf X-Wing_example.tar.gz
rm -rf X-Wing_example.tar.gz

## ---- Download the required reference panel ---- ##
cd example/data

### EUR reference panel
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LOGODetect/LOGODetect_1kg_EUR.tar.gz
tar -zxvf LOGODetect_1kg_EUR.tar.gz
rm -rf LOGODetect_1kg_EUR.tar.gz

### EAS reference panel
wget ftp://ftp.biostat.wisc.edu/pub/lu_group/Projects/XWING/ref/LOGODetect/LOGODetect_1kg_EAS.tar.gz
tar -zxvf LOGODetect_1kg_EAS.tar.gz
rm -rf LOGODetect_1kg_EAS.tar.gz

## ---- Applying LOGODetect ---- ##
cd example

mkdir ./results/LOGODetect

Rscript /LOGODetect/LOGODetect.R \
--sumstats ./data/sumstats/BMI_EUR.txt,./data/sumstats/BMI_EAS.txt \
--n_gwas 359983,158284 \
--ref_dir ./data/LOGODetect_1kg_ref \
--pop EUR,EAS \
--block_partition /LOGODetect/block_partition.txt \
--gc_snp /LOGODetect/1kg_hm3_snp.txt \
--out_dir ./results/LOGODetect \
--n_cores 20

```

### Output
* The output files include the regions with significant local genetic correlation and the annotation based on local genetic correlation as shown below
  ```bash
  head ./results/LOGODetect/LOGODetect_regions.txt

  chr    begin_pos    stop_pos     stat                pval                   pval_adj
  1      32159588     32175927     13.9007303828433    0.00479904019196161    0.0107599532725034
  1      74992278     75006328     23.659452286243     0.0001999600079984     0.00107756226532471
  1      107872936    107899979    13.5101005165975    0.0017996400719856     0.00492533072332902
  ```



# Citation
If you use the software of LOGODetect, please cite: 

Guo, H., Li, J. J., Lu, Q., Hou, L. [Detecting Local Genetic Correlations with Scan Statistics.](https://www.nature.com/articles/s41467-021-22334-6) Nature Communications, 2021.

Miao, J., Guo, H., Song, G., Zhao, Z., Hou, L., & Lu, Q. (2023). [Quantifying portable genetic effects and improving cross-ancestry genetic prediction with GWAS summary statistics.](https://www.nature.com/articles/s41467-023-36544-7) Nature Communications, 2023.

The genetic covariance estimation is adapted from `ldsc`, see Bulik-Sullivan, B., et al. [An Atlas of Genetic Correlations across Human Diseases and Traits.](https://www.nature.com/articles/ng.3406) Nature Genetics, 2015. 

The LD blocks partition is adapted from `LDetect`, see Berisa, Tomaz, and Joseph K. Pickrell. [Approximately independent linkage disequilibrium blocks in human populations.](https://academic.oup.com/bioinformatics/article/32/2/283/1743626/) Bioinformatics (2016).
