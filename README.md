# cDriver_tools

### Included files:

  - README file (this one) 
    * File with explanation how to run cDriver from command line
  - cDriver.run.R
    * Script to be run from command line and perform cDriver analysis
  - CLL_maf_like_file.txt
    * MAF-alike file (specification [link](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification))
    with additional columns (VAF, Damage_score, Ploidy, CCF_SNV) for 385 patients diagnosed with chronic lymphocytic leukemia. Othere column names are the same as in MAF files specification, but not all columns from MAF specification are nessesary (importan columns have values instead of dots in example file). Patients IDs don't have to follow MAF specification. This file is required for the cDriver analysis 
    (but only  additional column VAF is obligatory and others are recommended).  
    **VAF** is abbreviation for Variant Allele Frequency, calculated as ratio between alternative reads and total number of reads.  
    **Damage_score**  in this file is CADD score normalized between 0 and 1.  
    **Ploidy** is integer number, where 1 means there is a deletion and 3 or more means that in region of this variant there is a copy gain. 
    Ploidy equal to 2 is for normal state.  
    **CCF_SNV** is value between 0 and 1 and it represents in which fraction of cancer cells this ploidy is observed (for Ploidy 2, CCF_SNV is 1). 
  - CLL_Gold_Standard_Genes.txt  
    * File with gene names (HUGO symbol) in one column, without header
    
### Requirements

  - installed R (>= 3.1.0)
  - sudo rights (or that you are allowed to install R packages and their dependencies)

### Instructions to run cDrirver 

  Open your shell and type  first:
  ```Shell
Rscript cDriver.run.R -h
```
Here you can get familiar with input variables for cDriver tool. Then you can run it with example CLL file:
  ```Shell
Rscript cDriver.run.R -m CLL+SLL_385_noFilter.maf -p 5e-05 
```

### Output from cDriver tool
