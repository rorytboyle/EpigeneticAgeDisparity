# get summary table of allele frequencies

## CHANGE LINE 392 --recode A to --make-bed THIS WILL CREATE PLINK FILES WITH SNP INFO IN .BIM (TAKE COL 1 , 4 , 5, 6 from BIM FILE OR)
# DELIMIT .gnomAD .json file by - to have chr pos allele 1 allele  =2 into a df
# merge  columns by pmbb snp list (BIM file) and gnomAD df (two merges)
# merge 1 by chr pos and allele 1 pmbb = allele 1 gnomad and allele 2 pmbb = allele 2 gnomad
# merge 2 by chr pos and allele 1 pmbb = allele 2 gnomad and allele 2 pmbb = allele 1 gnomad
# Should then have a list of allele freqs (col for every snp with allele freq for each snp)
# make sure gnomAD has allele freqs for all the snps in genomic window (should be same in gnomAD as in pmbb query)

# then filter snps based on diff in allele freq btwn groups (maybe just apply easy threshold that allele freqs are not identifical between groups

# if there is a difference, than you run the latest plink files and remove SNPs with identical allele freqs between groups
# plink command to ouput a .bim file that contains the snps that have different frequencies THIS is then fed into the mQTL analysis
# --exclude .txt file with list of SNP names with identifical allele freqs (col 2) THEN run -- recode A --dash out