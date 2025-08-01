# Note: This code is provided for illustrative purposes only and comes with
# ABSOLUTELY no guarantee or warranty.

library(mice)
library(nnet) # for fitting multinomial logistic regression model.

################################################################################
# Read in data
################################################################################

zlist <- list(death1yr=0,female=0,age=0,acpulmed=0,cshock=0,diabetes=0,highbp=0,
  smokhx="",cva=0,dyslip=0,famhxcad=0,angina=0,cancer=0,dementia=0,pud=0,
  prevmi=0,asthma=0,depres=0,perartdis=0,chf=0,hyprthyr=0,
  as=0,prevrevasc=0,sysbp=0,diasbp=0,hrtrate=0,resp=0,hgb=0,wbc=0,
  sod=0,pot=0,glucose=0,urea=0,cr=0)

effect1.df <- data.frame(scanami_effect1.txt",zlist))
# NOTE: The dataset is not available for public dissemination. Please do not
# contact the author requesting the dataset.

remove(zlist)

effect1.df$smokhx <- as.factor(effect1.df$smokhx)

effect1.df$smokhx <- relevel(effect1.df$smokhx,ref="NEVER")
# Change reference level for smoking history to "NEVER".

################################################################################
# Multiple imputation using MICE: Use PMM for smoking history.
################################################################################

meth <- make.method(effect1.df)
meth["smokhx"] <- "pmm"
nimp <- round(100 * nic(effect1.df)/nrow(effect1.df),0)

effect1.imp.pmm <- mice(effect1.df,method=meth,m=nimp,maxit=20,seed=12022025,
  printFlag=F)

remove(meth,nimp)

################################################################################
# Multiple imputation using MICE: use multinomial logistic regression for
# smoking history.
################################################################################

meth <- make.method(effect1.df)
nimp <- round(100 * nic(effect1.df)/nrow(effect1.df),0)
print(nimp)

effect1.imp.mlr <- mice(effect1.df,method=meth,m=nimp,maxit=20,seed=12022025,
  printFlag=F)

remove(meth,nimp)

################################################################################
# Fit logistic regression model in imputed datasets
# Adjust for smoking history and all other covariates.
################################################################################

outcome.model.imp.mlr <- with(effect1.imp.mlr,glm(death1yr ~ female + age + 
  acpulmed + 
  cshock + diabetes + highbp + cva + dyslip + famhxcad + angina + cancer + 
  dementia + pud + prevmi + asthma + depres + perartdis + chf + hyprthyr + as +
  prevrevasc + sysbp + diasbp + hrtrate + resp + hgb + wbc + sod + pot + 
  glucose + urea + cr + smokhx,family="binomial"))

outcome.model.mlr <- summary(pool(outcome.model.imp.mlr),conf.int=TRUE)

outcome.model.imp.pmm <- with(effect1.imp.pmm,glm(death1yr ~ female + age + 
  acpulmed + cshock + diabetes + highbp + cva + dyslip + famhxcad + angina + 
  cancer + dementia + pud + prevmi + asthma + depres + perartdis + chf + 
  hyprthyr + as + prevrevasc + sysbp + diasbp + hrtrate + resp + hgb + wbc + 
  sod + pot + glucose + urea + cr + smokhx,family="binomial"))

outcome.model.pmm <- summary(pool(outcome.model.imp.pmm),conf.int=TRUE)

################################################################################
# Extract odds ratios and 95% CIs.
################################################################################

Variables <- outcome.model.mlr[,"term"]

OR.mlr <- round(exp(outcome.model.mlr[,"estimate"]),3)
CI.lower.mlr <- round(exp(outcome.model.mlr[,"conf.low"]),3)
CI.upper.mlr <- round(exp(outcome.model.mlr[,"conf.high"]),3)

OR.pmm <- round(exp(outcome.model.pmm[,"estimate"]),3)
CI.lower.pmm <- round(exp(outcome.model.pmm[,"conf.low"]),3)
CI.upper.pmm <- round(exp(outcome.model.pmm[,"conf.high"]),3)

MLR.summary <- paste0(OR.mlr," (",CI.lower.mlr,",",CI.upper.mlr,")")
PMM.summary <- paste0(OR.pmm," (",CI.lower.pmm,",",CI.upper.pmm,")")

Table1.df <- data.frame(Variables,MLR.summary,PMM.summary)

write.table(Table1.df,file="table1.txt",append=F,quote=FALSE,sep=" , ",
  row.names=FALSE,col.names=TRUE)
