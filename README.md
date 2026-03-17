## IMPORTANT: This repository is a fork of the original ccSAG with several fixes and improvements.

# 1. Installation
```
conda create -n ccSAG
conda activaet ccSAG
mamba install bioconda::fastp
mamba install bioconda::spades
mamba install bioconda::bwa
mamba install bioconda::blanstn
```

# 2. Usage

At first, clone this repository and change the file permission.
```
git clone https://github.com/ynishimuraLv/ccSAG.git
cd ccSAG
chmod +x ccSAG_cross_reference_cleaning.sh
chmod +x ccSAG_clamping.sh
```

## 2-1 ccSAG_cross_reference_cleaning.sh  

#### In this script, following steps are conducted.  
* Quality control of paired-end reads by fastp
* Assembling using each paired-end read file by SPAdes  
* Cross-reference read cleaning using bwa  
* Assembling using all cleaned reads by SPAdes  

#### example SAG data is contained in "example" directory.  
#### try
```
./ccSAG_cross_reference_cleaning.sh -o test -s example/Rawdata -c cleaned_merge_contigs.fasta
```

- The parameters can be also provided via a configuration file specified with the `-c` option.
- An example is ccSAG_cross_reference_cleaning.config

#### As the output files,  
* Cleaned contig file from all cleaned reads is in the assigned output directory.  
* Cleaned read file of each SAG is in the "output directory/QC" named "*_cleaned.fastq".  

## 2-2 ccSGA_clamping.sh  
#### In this script, following steps are conducted.  
* Assembling contig clamps using non-cleaned reads by SPAdes  
* connecting cleaned contigs by contig clamps using blastn  

#### try after executing 2-1(ccSAG_cross_reference_cleaning.sh)
```
./ccSAG_clamping.sh  -o test -r example/Rawdata -c example/cleaned_merge_contigs.fasta -C clamped_contigs.fasta -t 20
```

- The parameters can be also provided via a configuration file specified with the `-c` option.
- An example is ccSAG_clamping.config. 

#### As the output files,  
* Clamped contig file is in the assigned output directory  
