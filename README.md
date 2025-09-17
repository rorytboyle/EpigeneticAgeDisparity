Probe enrichment analysis workflow steps:
1) Run DMP by genetic ancestry: to be added
2) KnowYourCG analysis: to be added

mQTL analysis workflow steps:
1) Get genomic location of identified CpGs: get_windows.R
2) Locate imputed genotype data: locate_genotype_data_LPC.R
     - EDIT SO THAT DIRECTORY IS MOUNTED FROM WITHIN R
3) Get ROIs from PLINK files for these data: extract_ROI_from_PLINK_file.R
     - Update code with MAF flag
     - Update code with LD flag
     - secure password prompt
4) Run mQTL analysis: to be added
     - EDIT SO THAT ASSOCIATIONS SIG AT FDR ARE CORRECTLY COUNTED AND PLOTTED
     - Color mQTL plots by ancestry group
     - run entirely within R (avoid command line steps)
     - secure password prompt
5) Sanity check for nearby common SNPs: to be added
