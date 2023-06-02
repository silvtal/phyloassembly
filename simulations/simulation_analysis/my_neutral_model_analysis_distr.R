# tr0
input_s="simulation_results/my_neutral_model_glc/"
output="analysis_results/neutral"
m_inic <- paste0("X", c(1:10, 12))

if (!file.exists(output)) {system(paste("mkdir -p",output))}

# load real (endpoint) datasets
exp_f="../simulation_parameters/table_tr12.txt"
map_f="../simulation_parameters/maps/map_glc.csv"


# params
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
all_final_abund=list()# (1)
for (m in m_inic){
  ## (1) todos
  all_final_abund[[m]] <- read.csv(paste0(input_s,"/simul_",m,".csv"),
                                   row.names = NULL, check.names = F)[,-1,drop=FALSE]
  if (limit_simul_number) {
    limit=ncol(exp_per_m[[m]])
  } else {
    limit=nrow(all_final_abund[[m]])
  }
  message(paste0("A total of ",limit," simulations will be used for the null model analysis for all PCGs together (sample ",m,")."))

  all_final_abund[[m]] <- all_final_abund[[m]][1:limit, , drop=FALSE]
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
  simul_abund <- all_final_abund[[m]]
  real_abund <- exp_per_m[[m]] %>% my_transpose   ## traspongo porque aquí me interesan las SAMPLES no las OTUs
  real_abund <- real_abund[colnames(real_abund) %in% colnames(simul_abund)]
  
  # there might be OTUs that are absent from real_abund. We add empty columns.
  # if absent_OTUs is empty, there are no errors either.
  absent_OTUs <- colnames(simul_abund)[colnames(simul_abund) %!in% colnames(real_abund)]
  real_abund[absent_OTUs] <- rep(0,nrow(real_abund))
  # ----------------------------------------------------------------------------
  ## normalization so all replicates have the same sum of abundances 
  # this is done by sampling to the lowest total abund   
  bigger_real <- sum(rowSums(real_abund) > sum(simul_abund[1,]))# TODO CREO que todo este if no hace falta???? la línea de definición de "a" debería ser la primera??? pero not sure, "si lo hice así sería por algo"...
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
    # it'll throw a warning. They should be almost identical, the only differences   --> not identical anymore because we're not using the final abundances for the simulations, but the initial ones... so don't worry
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
  # ----------------------------------------------------------------------------
  
  # antes de normalizar a la relativa, guardo las absolutas para poder sacar la richness
  real_abund_abs <- real_abund %>% as.data.frame() %>%my_transpose() 
  real_abund=real_abund_abs/colSums(real_abund_abs) # normalizo
  real_abund[is.na(real_abund)] <- 0 # improbable, but : change NA to 0s
  # SIMUL DATA (calculate the mean abundance of each otu in all simulations)
  simul_abund_abs <- simul_abund %>% as.data.frame() %>%my_transpose() 
  simul_abund <- simul_abund_abs/colSums(simul_abund_abs) # normalizo
  simul_abund[is.na(simul_abund)] <- 0 # change NA to 0s
  
  # calculate div, even y richness for the mean abundances for the last transfer
  message(paste0("Obtaining results for sample '",m,"', all PCGs"))
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

  plots=list()
  colors=factor(c("simul"="4","real"="red"))
  for (metric in c("shannon","pielou","obs_rich", "faith")) {
    png(paste0(output,"/ECDF_",metric,"_",m,".png"))
    tryCatch(expr={
    plot_ecdf_KS2(data_simul,data_real, metric, m)
    },
    error=function(cond){
      write(append=TRUE,x=paste0("Error at ",m," ",metric,". ",cond),file=paste0(output,"/error_log.txt"))
    },
    warning=function(cond){
      plot_ecdf_KS2(data_simul,data_real, metric, m) # plot anyway !
      write(append=TRUE,x=paste0("Warning at ",m," ",metric,". ",cond),file=paste0(output,"/error_log.txt"))
    })
    dev.off()

    write.csv(data_real[[metric]],paste0(output,"/data_real_",metric,"_",m,".csv"))
    write.csv(data_simul[[metric]],paste0(output,"/data_simul_",metric,"_",m,".csv"))
  }
}