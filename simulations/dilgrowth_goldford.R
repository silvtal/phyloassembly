## =============================================================================
##                              About this script
## =============================================================================
## This script is a variation of silvtal/dilgrowth/scripts/dilgrowth.R, specific
## for Goldford et al.'s glucose dataset.
##  
## This old version simulates non-logistic, lineal growth for each functional
## group separately, creating separate files. It needs the initial abundances
## from one sample (a matrix column from an abundance table) and an index or list
## of OTUs (corresponding to a group; results.txt file from BacterialCore.py).
## Also simulates the changes in relative abundances over a series of
## dilution-transfer cycles.
## =============================================================================

# ===========
# mis datos [glc]
# ===========
library("untb") # simulate_timeseries uses "as.count"
library("tidyverse")
library("gsubfn") # for unpacking values
library("Rcpp")
library("parallel")
library("optparse") 
library("dplyr") # para bind_rows tras paralelización

option_list <- list(
  # input params
  make_option(c("-a", "--abuntable"), type = "character",
              default = NULL,
              help = "Abundance table"),
  make_option(c("-s", "--sample"), type = "character",
              default = NULL,
              help = "Sample name in the abundance table"),
  make_option(c("--subset"), type = "character",
              default = NULL,
              help = "Colon-separated list of OTUs to be selected from the abuntable (e. g. those belonging to a given PCG). NULL by default, meaning no selection)"),
  # simul params
  make_option(c("--dilution"), type = "character",
              default = "8 * 10 ** (-3)",
              help = "Dilution before each simulated transfer. Can be an equation (e.g. 10**(-3))"),
  make_option(c("--no_of_dil"), type = "integer",
              default = 12,
              help = "Number of dilutions/transfers"),
  make_option(c("--fixation_at"), type = "double",
              default = 1,
              help = "The simulations stop when only one bug is fixed at 100% [fixation_at=1], since it's meaningless to keep going when only one species is left. Nevertheless, we can consider an OTU is fixed when it has reached a lower relative abundance (e.g. fixation_at=0.95) and stop earlier."),
  make_option(c("--save_all"), type = "logical",
              default = FALSE,
              help = "save all intermediate states of the simulations? default FALSE"),
  make_option(c("--grow_step"), type = "integer",
              default = 1,
              help = "How many bugs are born on each iteration. 1 by default"),
  make_option(c("--fix_percentage"), type = "logical",
              default = TRUE,
              help = "If TRUE, uses a fixed percentage (relative abundance, given by --perc). If FALSE, randomly chooses a percentage in each simulation. These percentages are taken from a map given by --perc_map. Default: TRUE"),
  make_option(c("--perc"), type = "character",
              default = "1",
              help = "Desired final PCG abundance (if --fix_percentages==TRUE). 1 by default."),
  make_option(c("--perc_map"), type = "character",
              default = NULL,
              help = "Percentage map, see --fix_percentage"),
  make_option(c("--cores"), type = "integer",
              default = 1,
              help = "Number of cores to use in parallelization processes (mclapply). Default: 4.",
              metavar = "integer"),
  # output params
  make_option(c("--no_of_simulations"), type = "integer",
              default = 1,
              help = "Number of simulations"),
  make_option(c("--outputname"), type = "character",
              default = "out",
              help = "name for output .csv file"),
  make_option(c("--outdir"), type = "character",
              default = "results",
              help = "output directory")
)

# Custom functions
source("./functions_for_neutral_modelling.R")
source("./my_functions.R")

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

abuntable  <- opt$abuntable  # e.g. "./original_100percArbol/Tree/0.99/table.from_biom_0.99.txt"
s <- opt$sample      # e.g. "X2", col name from map
subset_otus <- opt$subset 

outdir <- opt$outdir # e.g. "my_neutral_model_v2_test_16_simuls_8_cores"
outputname <- opt$outputname # e.g. "X2_rep4"
cores <- opt$cores  # e.g. 16

save_all  <- opt$save_all
grow_step <- opt$grow_step

dilution <-  eval(parse(text=opt$dilution)) # e.g. 8 * 10 ** (-3)
no_of_dil <- opt$no_of_dil
fixation_at <- opt$fixation_at

no_of_simulations <- opt$no_of_simulations

fix_percentage <- opt$fix_percentage
if (fix_percentage){
  fixed_percentage <- as.numeric(opt$perc)
} else {
  perc_map_f <- opt$perc_map # e.g. "./map_final_perc.csv",
  percentage_map <- read.csv(perc_map_f, row.names = 1)
}

# leo matriz abundancias
list[exp, tax] <- get_abundance_and_tax_from_table(abuntable)
# creo output dir
system(paste("mkdir -p", outdir))

# ===========
# creo counts
# ===========
# Solo necesito la "muestra original", para ejecutar el simulate_timeseries
if (is.null(subset_otus)) {
  counts <- exp[s]
} else {
  subset_otus <- strsplit(subset_otus, ";")[[1]]
  subset_otus <- subset_otus[subset_otus %in% rownames(exp)] # we can only use those PCG OTUs that are actually present
  counts <- exp[subset_otus, s, drop=F]
}

# ==========
# simulation
# ==========

## start

## results will be stored here:
final_abund <- list()


# relative expected abundances of each PCG.
if (fix_percentage) {
  perc <- fixed_percentage
} else {
  # picked percentages from random real replicate
  perc <- sample_n(percentage_map[percentage_map$ORIG==s, ], 1)
}

# initial total abundance
total_counts <- sum(exp[s]) 
# final total abundance
abun_total   <- round(total_counts * perc)

# check first if there's anything to simulate
if (total_counts == 0) {
  message(paste0("There are no detectable OTUs in initial sample ",
                 s, ". Moving to next PCG..."))
} else if (total_counts * dilution < 3) {
  stop("EXIT: 3 or less bugs will be left after diluting! Consider changing your dilution factor.")
  
# start simulations if everything's OK
} else {
  abund_temp <- mclapply(X = 1:no_of_simulations,
                             FUN = function(iter) {
                               
                               # 1) simulation
                               
                               if (abun_total == 0) { # we can't simulate anything if it's 0
                                 start <- as.count(my_transpose(counts))
                                 start <- start[order(names(start))] # this order thing is to keep indices consistent
                                 trajectory <- matrix(0, ncol=length (start), nrow = no_of_dil+1)
                                 rownames(trajectory) <- 0:no_of_dil
                                 colnames(trajectory) <- names(start)
                                 trajectory["0",]=start
                               } else {
                                 trajectory <- simulate_timeseries(counts,
                                                                   dilution = dilution,
                                                                   no_of_dil = no_of_dil,
                                                                   fixation_at = fixation_at,
                                                                   grow_step = grow_step,
                                                                   abun_total = round(
                                                                     total_counts *
                                                                       perc),
                                                                   keep_all_timesteps = save_all)
                               }
                               
                               print(paste("Simulation", iter, "finished for", s))
                               
                               return(trajectory)
                               
                             }, mc.cores = cores)

  # save data
  message(paste0("Saving data for sample ", s, "..."))
  if (save_all == T) {
    for (timepoint in 1:(no_of_dil + 1)) {
      # pick all rows number "timepoint" from all the lists in all_abund
      # "temp" refers to a single timepoint
      temp <- lapply(abund_temp, FUN = function(traj) {
        return(traj[timepoint, , drop = F] %>% as.data.frame)}) %>%
        bind_rows %>% as_tibble # one file per timepoint
      
      write.csv(temp,
                file = paste0(outdir, "/simul_", outputname, "_t_", timepoint - 1, ".csv"))
    }
  } else {
    final_abund <- as.data.frame(bind_rows(abund_temp))
    write.csv(final_abund,
              file = paste0(outdir,"/simul_", outputname, ".csv"))
  }
}
