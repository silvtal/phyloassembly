# tr0
input_s="simulation_results/my_neutral_model_glc_PERC/"
pcgtable="../simulation_parameters/pcg_leaves.txt"
output="analysis_results/const_neutral"
m_inic=paste0("X", 1:12)

# tr1
input_s="simulation_results/my_neutral_model_glc_tr1_PERC/"
pcgtable="../simulation_parameters/pcg_leaves_tr1.txt"
output="analysis_results/const_neutral_tr1"
m_inic=c("X2", "X6")

if (!file.exists(output)) {system(paste("mkdir -p",output))}

# load real (endpoint) datasets
exp_f="../simulation_parameters/table_tr12.txt"
map_f="../simulation_parameters/maps/map_glc.csv"

limit_simul_number= F    # This is for normalizing variability. For example, if
# our input simulation files include 15 simulations
# but out samples have only 8 replicates each, this
# option selects only 8 simulations for the analysis.
# If FALSE, will just use all simulations available
################################################################################

# Cargar paquetes
library("ape")
library("tidyverse")
source("../my_functions.R")
library("gsubfn") # for unpacking values
library("stringr") # for renaming and str_detect
library("vegan")
library("microbiome")
library("PhyloMeasures")
library("untb")

# Cargar datos reales

## leo map y árbol
tn = ape::read.tree("../simulation_parameters/99_otus_nodes.tree")
map <- read.csv(map_f,sep=",")

## leo matriz abundancias
list[exp,tax]<-get_abundance_and_tax_from_table(exp_f)
exp_per_m <- list()
for (m in m_inic){
  exp_per_m[[m]] <- exp[map[map["ORIG"]==m,][,"SA"]]
}

# Cargar datos simulados
PCG_final_abund=list()
PCG_name <- c("Node35562","Node27828","others")

for (m in m_inic){
  ## por PCG
  PCG_final_abund[[m]] =list()
  for (PCG in PCG_name) {
    if (file.exists(paste0(input_s,"/simul_",PCG,"_",m,".csv"))) {
      PCG_final_abund[[m]][[PCG]] <- read.csv(paste0(input_s,"/simul_",PCG,"_",m,".csv"),
                                                        row.names = NULL, check.names = F)[,-1,drop=FALSE]
      if (limit_simul_number) {
        limit=ncol(exp_per_m[[m]])
      } else {
        limit=nrow(PCG_final_abund[[m]][[PCG]])
      }
      message(paste0("A total of ",limit," simulations will be used for the null model analysis for all PCGs together (sample ",m,")."))
      
      PCG_final_abund[[m]][[PCG]] <- PCG_final_abund[[m]][[PCG]][1:limit,,drop=FALSE] %>% na.omit() # por si algúna matriz tiene menos simulaciones que la de all.
    } else {
      message(paste0("PCG ",PCG," of sample ",m, " does not have an existing file, probably disappeared during simulations. Skipped."))
      m_inic <- m_inic[m_inic!=m]
    }
  }
}


# ==============================================================================
# ecdf: simulated vs real groups. 
#  - Diversity, evenness, richness, Faith's PD
# ==============================================================================
# For each initial sample (X1-X12), 
# - for each metric:
#    > create a ECDF distribution plot for the simulated data. See where
#      the real data falls.
#    > cuanto más a la DERECHA, más altos son los valores en general.

for (m in m_inic){
  ## ---------
  ## (*) all but with percentages
  ## ---------
  # start
  # simul
  simul_abund <- list()
  for (PCG in names(PCG_final_abund[[m]])) {
    simul_abund[[PCG]] <- PCG_final_abund[[m]][[PCG]]
  }
  # since they are already scaled (= perc have been applied), just put together all simul data
  simul_abund <- bind_cols(simul_abund)
  
  # real
  real_abund <- list()
  real_abund <- exp_per_m[[m]] %>% my_transpose
  real_abund <- real_abund[colnames(real_abund) %in% colnames(simul_abund)] 
  # there might be OTUs that are absent from real_abund. We add empty columns.
  # if absent_OTUs is empty, there are no errors either.
  absent_OTUs <- colnames(simul_abund)[colnames(simul_abund) %!in% colnames(real_abund)]
  real_abund[absent_OTUs] <- rep(0,nrow(real_abund))
  
  ## normalization so all replicates have the same sum of abundances 
  # this is done by sampling to the lowest total abund   
  bigger_real <- sum(rowSums(real_abund) > sum(simul_abund[1,]))
  if (bigger_real == nrow(real_abund)) { #if total real abund is always bigger
    for (row in 1:nrow(real_abund)){ # for each simulated sample (800)
      subsample <- isolate(real_abund[row,],
                           size=sum(simul_abund[1,]),
                           replace=FALSE)
      real_abund[row,][names(subsample)] <- subsample
    }
  } else {
    # (FIX?) to choose the target abundance, we'll choose from simul+real abunds
    a <- c(rowSums(real_abund), sum(simul_abund[1,]))
    # Normalization by sampling: use min total abund. If it's weirdly low 
    # it'll throw a warning. They should be almost identical, the only differences --> not really identical anymore because we're not using the final abundances for the simulations, but the initial ones... but whatever
    # being the removal of some OTU(s) with too low of an abundance during processing
    if (((mean(a)-min(a))/sd(a))>0.2) {
      warning(paste0("are ",m," real abundances OK?"))
      warning(paste0(rownames(real_abund),"\t"))
      warning(paste0(rowSums(real_abund),"\t"))
    }
    # normalize both simul and real !
    for (row in 1:nrow(real_abund)){ # for each real sample (8)
      if (sum(real_abund[row,]) < min(a)) {rep = TRUE} else {rep = FALSE}
      subsample <- isolate(real_abund[row,],
                           size=min(a),
                           replace=rep)
      real_abund[row,][names(subsample)] <- subsample
    } 
    for (row in 1:nrow(simul_abund)){ # for each real sample (8)
      if (sum(simul_abund[row,]) < min(a)) {rep = TRUE} else {rep = FALSE}
      subsample <- isolate(simul_abund[row,],
                           size=min(a),
                           replace=rep)
      simul_abund[row,][names(subsample)] <- subsample
    }
  }
  
  ## guardo las absolutas para poder sacar la richness
  real_abund_abs <- real_abund %>% as.data.frame()%>% my_transpose() 
  real_abund[is.na(real_abund)] <- 0 # improbable, but : change NA to 0s (for completeness/coherence with the next step)
  real_abund=real_abund_abs/colSums(real_abund_abs) # normalizo
  real_abund[is.na(real_abund)] <- 0 # improbable, but : change NA to 0s
  
  # SIMUL DATA (calculate the mean abundance of each otu in all simulations)
  simul_abund_abs <- simul_abund %>% as.data.frame() %>%my_transpose() 
  simul_abund_abs[is.na(simul_abund_abs)] <- 0 # change NA to 0s ## especialmente útil para PCG=others en por ej. X3
  simul_abund <- simul_abund_abs/colSums(simul_abund_abs) # normalizo
  simul_abund[is.na(simul_abund)] <- 0 # change NA to 0s
  
  ## calculate div, even y richness for the mean abundances for the last transfer
  message(paste0("Obtaining results for sample '",m,"', PCGs with_perc"))
  # real:
  shannon = as.data.frame(apply(as.matrix(real_abund),2,diversity,index="shannon")) %>% my_transpose()
  pielou = as.data.frame(apply(as.matrix(real_abund),2,evenness,index="pielou")) %>% my_transpose()
  obs_rich = as.data.frame(apply(as.matrix(real_abund_abs),2,richness,index="observed")) %>% my_transpose()
  faith= pd.query(tree=tn,matrix=1*(real_abund_abs!=0)%>%as.data.frame()%>%my_transpose(),standardize=FALSE,null.model = "uniform") %>% as.data.frame()
  data_real=list("shannon"=shannon,"pielou"=pielou,"obs_rich"=obs_rich, "faith"=faith)  
  # simul: 
  shannon = as.data.frame(apply(as.matrix(simul_abund),2,diversity,index="shannon")) %>% my_transpose()
  pielou = as.data.frame(apply(as.matrix(simul_abund),2,evenness,index="pielou")) %>% my_transpose()
  obs_rich = as.data.frame(apply(as.matrix(simul_abund_abs),2,richness,index="observed")) %>% my_transpose()
  faith= pd.query(tree=tn,matrix=1*(simul_abund_abs!=0)%>%as.data.frame()%>%my_transpose(),standardize=FALSE,null.model = "uniform") %>% as.data.frame()
  data_simul=list("shannon"=shannon,"pielou"=pielou,"obs_rich"=obs_rich, "faith"=faith)  
  
  for (metric in c("shannon","pielou","obs_rich", "faith")) {
    png(paste0(output,"/ECDF_",metric,"_",m,"_all_PERC.png"))
    tryCatch(expr={
      plot_ecdf_KS2(data_simul,data_real, metric, m, PCG="with_perc")
    },
    error=function(cond){
      write(append=TRUE,x=paste0("Error at ",m," ",PCG," ",metric,". ",cond),file=paste0(output,"/error_log_with_perc.txt"))
    },
    warning=function(cond){
      plot_ecdf_KS2(data_simul,data_real, metric, m, PCG="with_perc") # plot anyway !
      write(append=TRUE,x=paste0("Warning at ",m," ",PCG," ",metric,". ",cond),file=paste0(output,"/error_log_with_perc.txt"))
    })
    dev.off()
    # plot_boxplot(data_simul,data_real,metric)
    
    write.csv(data_real[[metric]],paste0(output,"/data_real_",metric,"_",m,"_all_PERC.csv"))
    write.csv(data_simul[[metric]],paste0(output,"/data_simul_",metric,"_",m,"_all_PERC.csv"))
  }
}
