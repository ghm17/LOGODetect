# LOGODetect
LOGODetect (LOcal Genetic cOrrelation Dectector) is a powerful tool to identify small segments that harbor local genetic correlation between two traits/diseases.

## Before starting
LOGODetect is built upon `R`, make sure that your R-version is no less than 3.5.0, and the R-package `data.table` is required to speed up reading large files. You will also need to install `plink` (downloaded [here](https://www.cog-genomics.org/plink/1.9)) to calculate LD matrix. Besides, in order to estimate the global genetic correlation between two diseases, you will need to install `Python 3` and `ldsc` (see [here](https://github.com/bulik/ldsc) for instructions).

## Tutorial
First download LOGODetect and the corresponding data.
        
    git clone git@github.com:ghm17/LOGODetect.git

**All the following steps should be carried out under the `./LOGODetect` directory!** Suppose we would like to find which part of the genome contributes to the genetic correlation between two diseases, e.g. bipolar disorder and schizophrenia. We need to prepare the following files:
* Two GWAS summary statistics files: We have prepared the example data for you in the directory `./example/BIP.txt.gz` and `./example/SCZ.txt.gz`. The GWAS summary statistics files need to be transformed into the standard format with the exact column names by yourself. The first few lines should look like this:

      head BIP.txt
      
      chr SNP pos A1 A2 N Z pval MAF
      1 rs4074137 1026707 C A 51710 0.360047754483072 0.7198 0.408
      1 rs11260590 1027070 G A 51710 -0.418067309888113 0.6765 0.106393579578418
      1 rs74048006 1027845 T G 51710 -0.448339477300638 0.6547 0.158787159156836
      1 rs76994018 1027846 T C 51710 -0.453310875113276 0.6513 0.158787159156836
      1 rs12077244 1028259 C T 51710 -0.456491245373629 0.6476 0.105393579578418
      1 rs6689308 1029805 A G 51710 -0.346923433004747 0.7273 0.160787159156836

The column names here denote chromosome, SNP rs-ID, position, effect allele, non-effect allele, sample size, z-score, p-value and minor allele frequency respectively.
* Plink bfiles of reference panel: These are files .bed/.bim/.fam format. We have already prepared these data (of EUR population from 1000 Genomes Project, with rare variants whose MAF < 1% filtered out) in directory `./Data/LD_matrix`. In particular, the reference genotype data are partitioned into blocks, each block covers 1 centiMorgan(cM), and the Plink bfiles of every two adjacent blocks are provided. These data are used to approximate the LD matrix of the sample cohorts.
* LD blocks partition: We have already prepared these data in directory `./Data/LD_block`. It contains a file `ldblock_count.txt` recording the LD block counts for each chromosome, and 22 files `ldblock_merged_chr@.txt` recording LD blocks partitions, each line represent the start and stop position of each LD block. On average, each LD block covers about 15 Mb. We perform LOGODetect within each LD block, and then aggregate the result together to control for FDR.

### Step 0-Calculate LD matrix
If this is the first time you use LOGODetect, run the following commands:
    
    cd ./LOGODetect
    source ./Code/calculateLD.txt  
    for chr in $(seq 22)
    do
    	Rscript ./Code/random_vector_generation.R ${chr}
    	Rscript ./Code/aggregation.R ${chr}
    done

The `random_vector_generation.R` and `aggregation.R` should be run for all chromosomes 1-22 (can be run in parallel manually by yourself). This step will generate random vectors which are needed in calculating the scan statistic null distribution. In particular, this step cost much memory and time, but it only needs to be run for one time!

### Step 1-Data preprocessing
        
    cd ./LOGODetect
    Rscript ./Code/preprocess.R BIP.txt SCZ.txt output_dir

The two arguments `BIP.txt` and `SCZ.txt` specify two input GWAS summary statistics files, the last argument `output_dir` specify the output directory. The argument `output_dir` should be kept consistent in the following steps. This step performs quality control and outputs two summary stat files for whole genome `output_dir/Data_QC/dat1.txt` and `output_dir/Data_QC/dat2.txt` and 44 summary stat files partitioned by chromosomes `output_dir/Data_QC/dat1_chr@.txt` and `output_dir/Data_QC/dat2_chr@.txt`.

### Step 2-Use ldsc to calculate heritability and genetic correlation of two traits
This step follows the instruction of `ldsc` to estimate genetic correlation, details see [here](https://github.com/bulik/ldsc).   

    mkdir output_dir/ldsc
    source activate ldsc
    cd ./LOGODetect
    python ./Data/ldsc/munge_sumstats.py \
    --sumstats BIP.txt \
    --out output_dir/ldsc/dat1_reformated \
    --merge-alleles Data/ldsc/w_hm3.snplist

    python ./Data/ldsc/munge_sumstats.py \
    --sumstats SCZ.txt \
    --out output_dir/ldsc/dat2_reformated \
    --merge-alleles Data/ldsc/w_hm3.snplist

    python ./Data/ldsc/ldsc.py \
    --rg output_dir/ldsc/dat1_reformated.sumstats.gz,output_dir/ldsc/dat2_reformated.sumstats.gz \
    --ref-ld-chr ./Data/ldsc/eur_w_ld_chr/baseline. \
    --w-ld-chr ./Data/ldsc/weights_hm3_no_hla/weights. \
    --out output_dir/ldsc/ldsc_trait1_trait2
    source deactivate

### Step 3-Calculate the scan statistic distribution under the null
    
    cd ./LOGODetect
    Rscript ./Code/BiScan_null.R chr N1 N2 output_dir

The `N1` and `N2` arguments denote the sample sizes of two summary stat files. The `BiScan_null.R` should be run for all chromosomes 1-22 (It is highly recommended to run this in parallel manually by yourself). 

### Step 4-Identify small regions that harbor local genetic correlation under different \theta
        
    cd ./LOGODetect
    Rscript ./Code/BiScan.R output_dir

### step 5-Select the best $\theta$ using the aggregated genetic covariance of identified regions as the metric and obtain the final results
        
    cd ./LOGODetect/Data
    wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_plinkfiles.tgz
    wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/hapmap3_snps.tgz
    tar -zxvf 1000G_Phase3_plinkfiles.tgz
    tar -zxvf hapmap3_snps.tgz
    
    cd ./LOGODetect
    Rscript ./Code/S-LDSC.R output_dir

The two files `1000G_Phase3_plinkfiles.tgz` and `hapmap3_snps.tgz` are needed when performing stratified-LDSC to calculate the aggregated genetic covariance of identified regions. Finally, `S-LDSC.R` selects the best $\theta$ and returns the final result of LOGODetect in `output_dir/LOGODetect_result.txt`. 

## Output
After running all the above steps, LOGODetect outputs a whitespace-delimited text file `output_dir/LOGODetect_result.txt`, with each row representing one small segment and the columns as such:
* `chr`: The chromosome. 
* `begin_pos`: The starting position of this detected small region (unit: MB).
* `stop_pos`: The stopping position of this detected small region (unit: MB).
* `stat`: The scan statistic value. Positive value means positive local genetic correlation between two traits. 
* `pval`: The p-value of this detected small region.
* `qval`: The q-value of this detected small region.
We have prepared the example output file for you in the directory `./example/LOGODetect_result.txt`. 

# Citation
If you use the software of LOGODetect, please cite: 

[Guo, et al. Detecting Local Genetic Correlations with Scan Statistics. Nature Communications, 2021](https://www.nature.com/articles/s41467-021-22334-6).

The genetic correlation estimation is adapted from `ldsc`, see [Bulik-Sullivan, B., et al. An Atlas of Genetic Correlations across Human Diseases and Traits. Nature Genetics, 2015](https://www.nature.com/articles/ng.3406). 

The LD blocks partition is adapted from `LDetect`, see [Berisa, Tomaz, and Joseph K. Pickrell. Approximately independent linkage disequilibrium blocks in human populations. Bioinformatics (2016)](https://academic.oup.com/bioinformatics/article/32/2/283/1743626/).
