# cDriver_tools

Here you can find wrapper script (*cDriver.run.R*) for cDriver R package ([link](https://github.com/hanasusak/cDriver)), with example input files. In this way you can run **driver genes analysis** as *black box*. 

The manuscript for cDriver is published at Scientific Reports and can be downloaded from: https://www.nature.com/articles/s41598-017-12888-1

### Included files:

  - README file (this one) 
    * File with explanation how to run cDriver from command line
  - cDriver.run.R
    * Script to be run from command line and perform cDriver analysis
  - CLL_maf_like_file.txt
    * MAF-alike file (specification [link](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification))
    with additional columns (VAF, Damage_score, Ploidy, CCF_CNV) for 385 patients diagnosed with chronic lymphocytic leukemia. From additional columns only VAF column is obligatory and others are recommended.  All other standard MAF columns names are the same as in MAF files specification, but not all columns from MAF specification are nessesary (importan columns have values instead of dots in this example file). Patients IDs don't have to follow MAF specification. This file is required for the cDriver analysis.  
    **VAF** is abbreviation for Variant Allele Frequency, calculated as ratio between alternative reads and total number of reads.   
    **Damage_score**  in this file is CADD score normalized between 0 and 1.  
    **Ploidy** is integer number, where 1 means there is a deletion and 3 or more means that in region of this variant there is a copy gain. 
    Ploidy equal to 2 is for normal state.  
    **CCF_SNV** is value between 0 and 1 and it represents in which fraction of cancer cells this ploidy is observed (for Ploidy 2, CCF_SNV is 1). 
  - CLL_Gold_Standard_Genes.txt  
    * File with gene names (HUGO symbol) in one column, without header, which are considered as gold standard for the cancer. For the CLL we compiled list of 22 genes which have previously been published in context to have increased number of somatic point mutations in CLL. 
    
### Requirements

  - installed R (>= 3.1.0)
  - sudo rights (or that you are allowed to install R packages and their dependencies)
  - MPFR C library version 3.0.0 or higher

### Instructions to run cDrirver 

Dowload files from this repository, and follow insturctions bellow:

  Open your shell and type  first:
  ```Shell
Rscript cDriver.run.R -h
```
Here you can get familiar with input variables for cDriver tool. Then you can run it with example CLL file:
  ```Shell
Rscript cDriver.run.R -m CLL_maf_like_file.txt -p 5e-05 -c CLL_Gold_Standard_Genes.txt 
```

### Output from cDriver tool

  - Log file (.txt extension)
  - Results tab separated file (.tsv extension)
  - Genereted plots in one PDF file (if option to plot was chosen)

# License
This project is licensed under the terms of the MIT license. 
