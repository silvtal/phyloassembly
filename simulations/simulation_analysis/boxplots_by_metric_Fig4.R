### Figure 4 (uno solo comparando x2 y x6)

# IDEAS : https://www.r-bloggers.com/2013/06/box-plot-with-r-tutorial/
# Plots a big plot for each metric. plots together the real, simul and const simul
# for each Xinic. Adds DTS values for simul and const simul bars.

inputs <- c("analysis_results/const_neutral/", "analysis_results/const_neutral_tr1/")
input_r <- "analysis_results/const_neutral_tr1/" # any folder does it
output  <- "analysis_results/Figures/"
dts_dir <- "analysis_results/comparison_results/"

if (!file.exists(output)) {system(paste("mkdir -p",output))}

# Cargar paquetes
library("tidyverse")
source("../my_functions.R")
# library("Rcpp")
# sourceCpp("./growth.cpp")
# source("../functions_for_neutral_modelling.R")
library("ggplot2")
library("reshape2")
library("wesanderson")
library("patchwork")
library("ggpubr") #get legend

m_inic <- c(2, 6) 
m_inic_sa <- c(8, 8) #porque algunas tienen 7, otras 8... y la 11 no existe
metric_pretty_names <- c("Faith's Evenness", "Shannon's Diversity", "Pielou's Phylogenetic Diversity", "Observed Richness")
metrics  <- c("faith", "shannon", "pielou", "obs_rich")
names(metric_pretty_names) <- metrics

all_alls <- list()

r = 0.5
rbig = 1.3 # tamaño del plot y la leyenda. invers proporc

layout <- "
    llllllll
    AABBCCDD
    AABBCCDD
    AABBCCDD
    AABBCCDD
    AABBCCDD
    AABBCCDD
    AABBCCDD
    AABBCCDD
    AABBCCDD
    AABBCCDD
  "
dts <- list()
bp <- list()

for (met in 1:length(metrics)) {
  metric <- metrics[met]
  # read the DTS tables
  dts[[metric]][["table"]] <- read.csv(paste0(dts_dir, "/", metric, "_vs_real.csv"), row.names = 1) # first row (empty rowname == all without fixed percs)
  dts[[metric]][["Experimental"]] <- rep("",3)
  names(dts[[metric]][["Experimental"]]) <- colnames(dts[[metric]][["table"]])
  dts[[metric]][["Constrained (transfer 0)"]] <- dts[[metric]][["table"]][1,]
  dts[[metric]][["Constrained (transfer 1)"]] <- dts[[metric]][["table"]][2,]
  
  ### HERE !!!
  
  all <- list()
  for (X in m_inic){
    done = 0
    tryCatch(
      expr={
        if (done==FALSE) {
          all[paste0(X,"_real")] <- list(read.csv(paste0(input_r,"/data_real_",metric,"_X",X,"_all_PERC.csv"),
                                                               row.names = 1, check.names = F)[[1]])  
          is.na(all[paste0(X,"_real")]) <- 0
          done <- 1
        }
        
        all[paste0(X,"_tr0")] <- list(read.csv(paste0(inputs[1],"/data_simul_",metric,"_X",X,"_all_PERC",".csv"),
                                                            row.names = 1, check.names = F)[[1]]) 
        
        all[paste0(X,"_tr1")] <- list(read.csv(paste0(inputs[2],"/data_simul_",metric,"_X",X,"_all_PERC",".csv"),
                                                             row.names = 1, check.names = F)[[1]]) 
        is.na(all[paste0(X,"_tr0")]) <- 0
        is.na(all[paste0(X,"_tr1")]) <- 0
        # real[[metric]][PCG] <- read.csv(paste0(input,"/data_real_",metric,"_X",X,PCG_name[PCG],".csv"),
        #                                         row.names = 1, check.names = F)[[1]]
        # all[[metric]][PCG] <- cbind(simul[[metric]][PCG], real[[metric]][PCG])
      },
      warning=function(cond){
        message(paste("ERROR, skipping",metric,"X",X,cond))
      },
      error=function(cond){
        message(paste("ERROR, skipping",metric,"X",X,cond))
      }
    )
  }
  
  melted_data <- melt(all) 
  Data <-  c() ; minic <- c()

  for (n in 1:nrow(melted_data)) {
    Data[n] <- substr(melted_data[2][n,], start = 3, stop = nchar(melted_data[2][n,]))
    Data[n] <- c("real"="Experimental",
                 "tr0"="Constrained (transfer 0)",
                 "tr1"="Constrained (transfer 1)")[Data[n]] # constrained simulations
    minic[n] <- sprintf('X%s', substr(melted_data[2][n,], start = 1, stop = 1))
  }
    Data <- factor(Data, levels=c("Experimental","Constrained (transfer 0)","Constrained (transfer 1)"))
  melted_data <- cbind(melted_data, Data, minic)

  # position for labels
  labeldat <- melted_data %>%
    group_by(Data, minic) %>%
    summarize(ypos = max(value) + max(max(value)*0.1,.75))

  for (n in 1:nrow(labeldat)) {
    l <-  dts[[metric]][[as.character(labeldat$Data[n])]][[labeldat$minic[n]]]
    # if (is.numeric(l)) {l <- round(l,3)} 
    labeldat[n, "labels"] <- l %>% as.character()
  }
  
  for (mini in unique(labeldat$minic)) {
    if (labeldat[labeldat$minic==mini & labeldat$Data=="Constrained (transfer 0)",]$labels >
             labeldat[labeldat$minic==mini & labeldat$Data=="Constrained (transfer 1)",]$labels) {
      labeldat[labeldat$minic==mini & labeldat$Data=="Constrained (transfer 0)",]$labels <- ""
      labeldat[labeldat$minic==mini & labeldat$Data=="Constrained (transfer 1)",]$labels <- "＊"
    } else {
      labeldat[labeldat$minic==mini & labeldat$Data=="Constrained (transfer 0)",]$labels <- "＊"#"⋆"
      labeldat[labeldat$minic==mini & labeldat$Data=="Constrained (transfer 1)",]$labels <- ""
    }
  }
      # A tibble: 9 × 4
      # Groups:   Data [3]
      # Data        minic  ypos labels              
      # <chr>       <chr> <dbl> <chr>               
      #   1 const simul X10    3.09 "0.0431017115157696"
      # 2 const simul X12    1.82 "0.260374805236532" 
      # 3 const simul X9     2.28 "0.104810555493115" 
      # 4 real        X10    1.87 ""                  
      # 5 real        X12    2.46 ""                  
      # 6 real        X9     1.02 ""                  
      # 7 simul       X10    6.04 "0.570877772348109" 
      # 8 simul       X12    1.93 "0.234984186430073" 
      # 9 simul       X9     2.66 "0.195510579617056" 
  
  bp[[met]] <-
    ggplot(data=melted_data, mapping=aes(y=factor(Data, levels = c("Constrained (transfer 0)", "Constrained (transfer 1)", "Experimental")),#factor(minic, levels=paste0("X", m_inic)),
                                         x=value,  # ^^^ vvv estos levels afectan a todos los subgrafos !
                                         fill=Data,
                                         color=Data), lwd=2) +
    geom_boxplot(orientation = "y") +
    geom_text(data = labeldat, aes(label = labels, x = ypos), 
              show.legend = FALSE, linewidth = 15*r) +
    facet_grid(factor(minic, levels=paste0("X", m_inic)) ~ . , switch = "y") +
    scale_fill_manual(values = c("Experimental"="#00AFBB",
                                "Constrained (transfer 1)"= "#99190b", #custom.col[1],
                                "Constrained (transfer 0)"="#FC4E07"),
                      name   = "") +#wes_palette(n=3, name="BottleRocket1")) +
    # scale_fill_brewer(palette="rocket") + #c("real"="red","tr0"="sienna","simul_perc"="royalblue2")) +
    # scale_fill_manual(values=c("real"="red","tr0"="sienna","simul_perc"="royalblue2")) +
    # scale_color_brewer(palette="viridis") +
    scale_color_manual(values = c("Experimental" = custom.col[8],
                                 "Constrained (transfer 1)"   = "#3d1f1c", #custom.col[4], #custom.col[2],
                                 "Constrained (transfer 0)"   = "#542c1c"), #custom.col[4]),
                       name   = "") +
    
    # scale_color_manual(values=c("real"="red","tr0"="sienna","simul_perc"="royalblue2"))
    theme(legend.position = "none",
          legend.text=element_text(linewidth=(55/rbig)*r),
          legend.title=element_text(linewidth=(65/rbig)*r),
          axis.text.x = element_text(linewidth=30*r),
          axis.text.y = element_blank(), #element_text(linewidth=30*r, angle = 45),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(linewidth =35*r),
          axis.title.x = element_text(linewidth =35*r),
          panel.grid = element_blank(),
          strip.text.y = element_text(linewidth=30*r, color="black"),# <- facet
          strip.text.y.left = element_text(angle = 0), 
          strip.background = element_rect(colour="black", fill="white", 
                                          linewidth=1, linetype="solid")) + 
    labs(y = NULL, x = metric_pretty_names[[metric]])
}

leg <- get_legend(bp[[3]], position = "bottom")
leg <- as_ggplot(leg) + theme(text=element_text(linewidth=80*r))

combined <- wrap_plots(#t = title,
                       A = bp[[1]], B =  bp[[2]], C = bp[[3]], D = bp[[4]], l = leg, design = layout)

png(filename=paste0(output,"Figure 4.png"),
    height=(2100/rbig)*r, width=((2970)/rbig)*r)
plot((combined + theme(legend.position = "none")))
dev.off()

