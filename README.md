mQTL analysis workflow

Steps
1) Get genomic location of identified CpGs: get_windows.R
2) Locate imputed genotype data: locate_genotype_data_LPC.R
     - EDIT SO THAT DIRECTORY IS MOUNTED FROM WITHIN R
3) Get ROIs from PLINK files for these data: extract_ROI_from_PLINK_file.R
     - Update code with MAF flag
4) Run mQTL analysis: TO BE ADDED
     - EDIT SO THAT ASSOCIATIONS SIG AT FDR ARE CORRECTLY COUNTED AND PLOTTED
     - Color mQTL plots by ancestry group
