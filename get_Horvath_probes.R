# get Horvarth CpGs from https://github.com/isglobal-brge/methylclock/blob/main/data/coefHorvath.rda
load("coefHorvath.rda")
horvath_CpG <- coefHorvath$CpGmarker
horvath_CpG_weights <- coefHorvath$CoefficientTraining

df <- data.frame(CpG=horvath_CpG, Weights=horvath_CpG_weights)
# remove the frist row which is intercept
df <- df[-1, ]
saveRDS(df, file="20251124_cpg_Horvath_weights.rds")