### R implementation of PENGUINS
options(stringsAsFactors = FALSE)

## Require dependencies, and install automatically if user doesn't already have them installed
repos <- "https://cloud.r-project.org"
if (!require(data.table)) {
  install.packages("data.table", repos = repos)
  library(data.table)
}
if (!require(optparse)) {
  install.packages("optparse", repos = repos)
  library(optparse)
}
if (!require(GenomicSEM)) {
  if (!require(devtools)) {
    install.packages("devtools", repos = repos)
    library(devtools)
  }
  install_github("GenomicSEM/GenomicSEM")
  library(GenomicSEM)
}

## Read the arguments into R
option_list <- list(
  # required
  make_option("--sumstat_files", action = "store", default = NULL, type = "character"),
  #make_option("--y_path", action = "store", default = NULL, type = "character"),
  #make_option("--x_path", action = "store", default = NULL, type = "character"),
  make_option("--output", action = "store", default = NULL, type = "character"),
  make_option("--version", action = "store", default = "individual", type = "character"),
  make_option("--N", action = "store", default = NULL, type = "character"),
  make_option("--Ns", action = "store", default = NULL, type = "integer"),
  # optional
  make_option("--hm3", action = "store", default = "./w_hm3.noMHC.snplist", type = "character"),
  make_option("--ld", action = "store", default = "./eur_w_ld_chr/", type = "character"),
  make_option("--wld", action = "store", default = "./eur_w_ld_chr/", type = "character")
)

### step 0. specify the inputs
opt <- parse_args(OptionParser(option_list = option_list))
# Required inputs
sumstat_files <- unlist(strsplit(opt$sumstat_files, split = ","))
output_path <- opt$output
version <- opt$version
N <- as.numeric(unlist(strsplit(opt$N, split = ",")))
Ns <- opt$Ns
# Optional inputs
hm3 <- opt$hm3
ld <- opt$ld
wld <- opt$wld
cat("\nOutput Path: ", output_path)
cat("\nN: ", N)
cat("\nhm3: ", hm3)
cat("\nld: ", ld)
cat("\nwld: ", wld)
N <- c(766345, 269867)
Ns <- 195653
output_path <- "./penguins/output"
ld <- "/z/Comp/lu_group/Software/ldsc/Inputs/EUR_1KGphase1/eur_w_ld_chr/"
wld <- ld
munged_files <- c(paste0(output_path, "/munge1"), paste0(output_path, "/munge2"))

### step 1: run LD score regression for inputed GWAS sumstat files using GenomicSEM
munge(files = sumstat_files, N = N,
      hm3 = hm3, trait.names = munged_files,
      log.name = paste0(output_path, "/munge"),
      parallel = TRUE, cores = 2)
## sample prev and population prev currently assuming continuous variables
LDSCoutput <- ldsc(traits = paste0(munged_files, ".sumstats.gz")
      , sample.prev = c(NA, NA)
      , population.prev = c(NA, NA)
      , ld = ld, wld = wld)
# extract data from LDSCoutput
#y.var <- LDSCoutput[["S"]][1]
x.var <- LDSCoutput[["S"]][4]
xy.cov <- LDSCoutput[["S"]][2]
k <- nrow(LDSCoutput$S)
se <- matrix(0, k, k)
se[lower.tri(se, diag = TRUE)] <- sqrt(diag(LDSCoutput$V))
se.var <- se[4]
se.cov <- se[2]
print(LDSCoutput$I)
alpha <- LDSCoutput$I[2]
alpha
## STEP 2a - calculate phenotypic covariance between the exposure and outcome on the individual-level phenotype data
# individual level data (phenotype, covariates, run lm to get y.res and x.res)
if (version == "individual") {
  phen <- fread(paste0("/z/Comp/lu_group/Members/zzhao/PRS_covariate_SG/UKB_EA_scan/UKB_phen/clean/LLC.txt"))
  head(phen)
  colnames(phen)[3] <- "Y"
  covariate <- fread(paste0("/z/Comp/lu_group/Members/zzhao/PRS_covariate_SG/UKB_EA_scan/UKB_phen/clean/covariates.txt"))
  colnames(covariate)[3] <- "X"
  head(covariate)
  phen <- base::merge(x=phen, y=covariate ,by.x = c("IID", "FID") ,by.y = c("IID", "FID"))
  phen <- na.omit(phen)

  y.res <- resid(lm(Y ~ SEX + YOB + SEXYOB + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = phen))
  x.res <- resid(lm(X ~ SEX + YOB + SEXYOB + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = phen))
  y.res <- scale(y.res)
  x.res <- scale(x.res)
  dat <- data.frame(Y = y.res, X = x.res)

  cov.covxy <- cov(dat$Y , dat$X) - xy.cov
  cov.varx <- var(dat$X) - x.var
  cov.beta <- cov.covxy / cov.varx

  mu_x <- cov.covxy
  mu_y <- cov.varx
  sigma2_x <-  (1 + cov(dat$Y, dat$X)^2) / (nrow(dat) - 1) + se.cov^2

  sigma2_y <-  2 / (nrow(dat) - 1) + se.var^2
  se.est <- sqrt(sigma2_x / (mu_y^2) + (mu_x^2) * sigma2_y / (mu_y^4))

  # calculate marginal regression results as a baseline to output if using individual level data
  lm.marg <- lm(Y ~ X, data = dat)
  se.marg <- summary(lm.marg)$coefficient[2, 2]
  beta.marg <- as.numeric(coefficients(lm.marg)[2])
  p.marg <- 2 * pnorm(abs(beta.marg / se.marg), lower.tail = FALSE)
} else {
  ## Step 2b - SUMSTATS ONLY
  ldsc.cov <- (sqrt(prod(N)) / Ns) * alpha

  cov.covxy <- ldsc.cov - xy.cov
  cov.varx <- 1 - x.var
  cov.beta <- cov.covxy / cov.varx

  mu_x <- cov.covxy
  mu_y <- cov.varx
  sigma2_x <-  (1 + ldsc.cov^2) / (Ns - 1) + se.cov^2

  sigma2_y <-  2 / (Ns - 1) + se.var^2
  se.est <- sqrt(sigma2_x / (mu_y^2) + (mu_x^2) * sigma2_y / (mu_y^4))
}

p.cov <- 2 * pnorm(abs(cov.beta / se.est), lower.tail = FALSE)

## Output data
out <- c("BETA" = cov.beta, "SE" = se.est, "P" = p.cov)
out <- c("BETA" = cov.beta, "SE" = se.est, "P" = p.cov, "BETA_Marginal_Regression" = beta.marg, "SE_Marginal_Regression" = se.marg, "P_Marginal_Regression" = p.marg)
out
# executable
#write.table(out,paste0("./output/result.txt"),col.names=T,row.names=F,sep="\t",quote=F)
