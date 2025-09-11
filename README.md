mQTL analysis workflow

Steps
1) Get genomic location of identified CpGs: get_windows.R
2) Locate imputed genotype data: locate_genotype_data_LPC.R
  2a) EDIT SO THAT DIRECTORY IS MOUNTED FROM WITHIN R
4) Get ROIs from PLINK files for these data: extract_ROI_from_PLINK_file.R
  4a) Update code with MAF flag
5) Run mQTL analysis: TO BE ADDED
   5a) EDIT SO THAT ASSOCIATIONS SIG AT FDR ARE CORRECTLY COUNTED AND PLOTTED
   5b) Color mQTL plots by ancestry group
