# ==============================================================================
# DTS for all kind of simuls. We want to obtain the DTS between two datasets:
# "real samples" and "better simuls". The "simuls" datasets will be:
#   - fixed percentage simuls t0
#   - fixed percentage simuls t1
# ==============================================================================
# .../boxplots + tablas_DTS - draft2022/readme   simul_vs_simul

# Each one of these will be compared to neutral simuls
inputs  <- list("fixed_perc"    = "analysis_results/const_neutral/",
                "tr1"           = "analysis_results/const_neutral_tr1/")

# We'll load the experimental data too (any folder does it)
input_real <- "analysis_results/const_neutral_tr1/"

dataset_name <- names(inputs)

output  <- "analysis_results/comparison_results/"

r_dig <- 3 # digits to round to
if (!file.exists(output)) {system(paste("mkdir -p", output))}
metrics  <- c("faith", "shannon", "pielou", "obs_rich")

# - KS -----------------
# https://www.itl.nist.gov/div898/handbook/eda/section3/eda35g.htm
# The K-S test is based on the maximum distance between these two curves. 
# Another advantage is that it is an exact test (the chi-square goodness-of-fit test depends on an adequate sample size for the approximations to be valid). Despite these advantages, the K-S test has several important limitations:
# 
# It only applies to continuous distributions.
# It tends to be more sensitive near the center of the distribution than at the tails.
# Perhaps the most serious limitation is that the distribution must be fully specified. That is, if location, scale, and shape parameters are estimated from the data, the critical region of the K-S test is no longer valid. It typically must be determined by simulation.

# - AD test ------------
# Several goodness-of-fit tests, such as the Anderson-Darling test and the Cramer 
# Von-Mises test, are refinements of the K-S test
## Anderson-Darling: gives more weight to the tails than does the K-S test

# - DTS ----------------
library("twosamples")
# Testing the hypothesis that two samples come from the same distribution using randomiza-
#   tion to create p-values. Included tests are: Kolmogorov-Smirnov, Kuiper, Cramer-
#   von Mises, Anderson-Darling, Wasserstein, and DTS. The de-
#   fault test (two_sample) is based on the DTS test statistic, as it is the most power-
#   ful, and thus most useful to most users.
# The DTS test statistic builds on the Wasserstein distance by using a weight-
#   ing scheme like that of Anderson-Darling

##  The DTS statistic is the reweighted integral of the distance
## between the two ECDFs.

# Output is a length 2 Vector with test stat and p-value in that order. That 
# vector has 3 attributes – the sample sizes of each sample, and the number of bootstraps performed for the pvalue.

# If the Wasserstein metric improves on CVM by moving it into the realm of 
# interval data, then DTS improves on AD by doing the same. Alternately – DTS 
# offers the same improvement over Wasserstein that AD offers over CVM.

# Broadly, under the Null that the two distributions are the same, we can 
# bootstrap samples with the same sizes from the joint distribution. Each 
# resampling lets us calculate the test statistic again -- generating the
# distribution of test statistics under the null. Then the observed test 
# statistic is compared to the null distribution, and the p-value is the portion
# of the null distribution that our observed test statistic is smaller than.

results <- list()
# - let's go -----------
for (metric in metrics){
  results[[metric]] = data.frame(matrix(nrow = length(dataset_name),
                                        ncol = 3, 
                                        dimnames = list(dataset_name,
                                                        c(paste0("X",c(2, 6)), "average"))))
  for (X in c(2, 6)){
    for (data in dataset_name){
      input_s <- inputs[[data]]
      tryCatch(
        expr={
          # Cada archivo se compara con las simulaciones neutras
          simul <- read.csv(paste0(input_s,"/data_simul_",metric,"_X",X,"_all_PERC",".csv"),
                            row.names = 1, check.names = F)
          # neutral <- read.csv(paste0("~/AAA/2021-06-28_my_null_model/my_neutral_model_glc_CLUS_ene23/analysis_distr/data_simul_",metric,"_X",X,".csv"),
          # row.names = 1, check.names = F)
          real <- read.csv(paste0(input_real, "/", "data_real_", metric,"_X",X,"_all_PERC",".csv"),
                           row.names = 1, check.names = F)
          # NA handling
          if (sum(is.na(real))>0){
            message(paste("A total of",sum(is.na(real)),"NAs at real dataset //",metric, X, data))
          }
          real[is.na(real)] <- 0
          
          if (sum(is.na(simul))>0){
            message(paste("A total of",sum(is.na(simul)),"NAs at simul dataset //",metric, X, data))
          }
          simul[is.na(simul)] <- 0
          
          # save result. If p-value > 0.05, ponerlo entre paréntesis.
          r <- two_sample(real[[1]], simul[[1]], nboots=2000)
          if (r[["P-Value"]] < 0.05) {
            results[[metric]][data, c(paste0("X",X))] <- round(r[["Test Stat"]],r_dig)
          } else {
            results[[metric]][data, c(paste0("X",X))] <- 0 # <- lo hago 0 porque si no es significativo quiere decir que son la misma distribución (aparentemente al menos)
            # paste0("(",round(r[["Test Stat"]],r_dig),")")
          }
        },
        error=function(cond){
          message(paste("ERROR at", metric, X, data, cond))
        }
      )
      
    }
  }
  # ALSO add mean of SIGNIFICANT && not NA values in an additional column
  results[[metric]]["average"] <-
    apply(results[[metric]], MARGIN = 1, function(x){
      mean(as.numeric(
        x [!is.na(x) & !grepl("(", x, fixed = TRUE)]
      ))
    })
  rownames(results[[metric]]) <- c("Mean percentages",
                                   "Transfer 1 start")
  
  write.csv(results[[metric]], file=paste0(output,"/",metric,"_vs_real.csv"))
}
