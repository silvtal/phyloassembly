# Suppl Figure 2

# IDEAS : https://www.r-bloggers.com/2013/06/box-plot-with-r-tutorial/
# Plots a big plot for each metric. plots together the real, simul and const simul
# for each X inic. Adds DTS values for simul and const simul bars.


inputs <- c("analysis_results/neutral/", "analysis_results/const_neutral/")
output  <- "analysis_results/Figures/SupplFig3/"
dts_dir <- "analysis_results/comparison_results_SupplFig3/"

if (!file.exists(output)) {system(paste("mkdir -p",output))}

# Cargar paquetes
library("tidyverse")
source("~/Apps/my_functions.R")
library("ggplot2")
library("reshape2")
library("wesanderson")
library("patchwork")
library("ggpubr") #get legend

m_inic <- c(1:10, 12) 
m_inic_sa <- c(8, 8,8,8,8,8,8,8,7,7,NULL,7) #porque algunas tienen 7, otras 8... y la 11 no existe
datasets <- c("", "_all_PERC")
metric_pretty_names <- c("Faith's Evenness", "Shannon's Diversity", "Pielou's Phylogenetic Diversity", "Observed Richness")
metrics  <- c("faith", "shannon", "pielou", "obs_rich")
names(metric_pretty_names) <- metrics

all_alls <- list()

tramos = list(c(2:5),c(6:9),c(10,11))# Ignore X1: start by 2
# tramos = as.list(1:11)

r = 0.5
rbig = 1.3 # tamaño del plot y la leyenda. invers proporc

dts <- list()
for (metric in metrics) {
  # 15feb2022 read the DTS tables
  dts[[metric]][["table"]] <- read.csv(paste0(dts_dir,"/",metric,"_vs_real.csv"))
  dts[[metric]][["real"]] <- rep("",13)
  names(dts[[metric]][["real"]]) <- names(dts[[metric]][["table"]][1,])
  dts[[metric]][["neutral"]] <- dts[[metric]][["table"]][1,]
  dts[[metric]][["const neutral"]] <- dts[[metric]][["table"]][2,]
  
  # inicializo ya el plot
  bp <- list()
  layout <- "
    tttttt
    AAAAAA
    AAAAAA
    AAAAAA
    AAAAAA    
    AAAAAA    
    AAAAAA    
    AAAAAA
    AAAAAA
    AAAAAA
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    CCCCDD
    CCCCDD
    CCCCDD
    CCCCDD
    CCCCDD
    CCCCDD
    CCCCDD
    CCCCDD
    CCCCDD
  "
  layout2 <- "
    tttttt
    AAAAAA
    AAAAAA
    AAAAAA
    AAAAAA    
    AAAAAA    
    AAAAAA    
    AAAAAA
    AAAAAA
    AAAAAA
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    BBBBBB
    CCCCCC
    CCCCCC
    CCCCCC
    CCCCCC
    CCCCCC
    CCCCCC
    CCCCCC
    CCCCCC
    CCCCCC
  "
  
  png(filename=paste0(output,"BP_",metric,".png"),
      width=(2100/rbig)*r, height=(2970/rbig)*r)
  
  for (tr in 1:length(tramos)){
    all <- list()
    for (X in m_inic[tramos[[tr]]]){
      done = 0
      for (ds in 1:length(datasets)) {
        tryCatch(
          expr={
            if (done==FALSE) {
              all[paste0(X,"_real")] <- list(read.csv(paste0(inputs[[1]],"/data_real_",metric,"_X",X,".csv"),
                                                                   row.names = 1, check.names = F)[[1]])  
              is.na(all[paste0(X,"_real")]) <- 0
              done <- 1
            }
            
            all[paste0(X,"_simul",datasets[ds])] <- list(read.csv(paste0(inputs[[ds]],"/data_simul_",metric,"_X",X,datasets[ds],".csv"),
                                                                row.names = 1, check.names = F)[[1]]) 

            is.na(all[paste0(X,"_simul",datasets[ds])]) <- 0
            # real[[metric]][ds] <- read.csv(paste0(input,"/data_real_",metric,"_X",X,datasets[ds],".csv"),
            #                                         row.names = 1, check.names = F)[[1]]
            # all[[metric]][ds] <- cbind(simul[[metric]][ds], real[[metric]][ds])
          },
          warning=function(cond){
            message(paste("ERROR, skipping",ds,metric,"X",X,cond))
          },
          error=function(cond){
            message(paste("ERROR, skipping",ds,metric,"X",X,cond))
          }
        )
      }
    }
    
    melted_data <- melt(all) 
    Data <-  c() ; minic <- c()

    for (n in 1:nrow(melted_data)) {
      Data[n] <- str_remove_all(string = melted_data[2][n,], pattern ="[_1234567890]")
      Data[n] <- c("real"="real",
                   "simul"="neutral",
                   "simulallPERC"="const neutral")[Data[n]] # constrained simulations
      minic[n] <- sprintf('X%s', str_remove_all(string = melted_data[2][n,], pattern ="[_a-zA-Z]"))
    }
    Data   <-factor(Data, levels=c("real","neutral","const neutral"))
    melted_data <- cbind(melted_data, Data, minic)

    # label position
    labeldat <- melted_data %>%
      group_by(Data, minic) %>%
      summarize(ypos = max(value) + max(max(value)*0.1,.75))

    # numbers as labels
    for (n in 1:nrow(labeldat)) {
      l <-  dts[[metric]][[as.character(labeldat$Data[n])]][[labeldat$minic[n]]]
      # if (is.numeric(l)) {l <- round(l,3)} 
      labeldat[n, "labels"] <- l %>% as.character()
    }
  
    # just asterisks (if you want numbers, comment this out)
    for (mini in unique(labeldat$minic)) {
      if (labeldat[labeldat$minic==mini & labeldat$Data=="neutral",]$labels >
          labeldat[labeldat$minic==mini & labeldat$Data=="const neutral",]$labels) {
        labeldat[labeldat$minic==mini & labeldat$Data=="neutral",]$labels <- ""
        labeldat[labeldat$minic==mini & labeldat$Data=="const neutral",]$labels <- "＊" # smaller dts == smaller difference == marked 
      } else {
        labeldat[labeldat$minic==mini & labeldat$Data=="neutral",]$labels <- "＊"
        labeldat[labeldat$minic==mini & labeldat$Data=="const neutral",]$labels <- ""
      }
    }
    
    bp[[tr]] <-
      ggplot(data=melted_data, mapping=aes(x=Data,#factor(minic, levels=paste0("X", m_inic)),
                                           y=value,  # ^^^ vvv estos levels afectan a todos los subgrafos !
                                           fill=Data,
                                           color=Data), lwd=2) +
      geom_boxplot() +
      geom_text(data = labeldat, aes(label = labels, y = ypos), 
                show.legend = FALSE, size = 10*r) +
      facet_grid(. ~ factor(minic, levels=paste0("X", m_inic))) +
      scale_fill_manual(values= c("real"="#00AFBB", "neutral"= "#E7B800", "const neutral"="#FC4E07")) +#wes_palette(n=3, name="BottleRocket1")) +
      # scale_fill_brewer(palette="rocket") + #c("real"="red","neutral"="sienna","simul_perc"="royalblue2")) +
      # scale_fill_manual(values=c("real"="red","neutral"="sienna","simul_perc"="royalblue2")) +
      # scale_color_brewer(palette="viridis") +
      scale_color_manual(values= c("real"=custom.col[8], "neutral"= custom.col[2], "const neutral"=custom.col[4])) +
      # scale_color_manual(values=c("real"="red","neutral"="sienna","simul_perc"="royalblue2"))
      theme(legend.position = "none", 
            legend.text=element_text(size=(55/rbig)*r),
            legend.title=element_text(size=(65/rbig)*r),
            axis.text.y = element_text(size=30*r),
            axis.text.x = element_blank(), #element_text(size=30*r, angle = 45),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size =35*r),
            panel.grid = element_blank(),
            strip.text.x = element_text(size=30*r, color="black"), # <- facet
            strip.background = element_rect(colour="black", fill="white", 
                                            size=1, linetype="solid")) + 
      labs(x = NULL, y = metric)
  }
  leg <- get_legend(bp[[3]], position = "right")
  leg <- as_ggplot(leg) + theme(text=element_text(size=80*r))
  
  title <- ggplot() + 
    theme_void() +
    # ggtitle(metric_pretty_names[[metric]], subtitle = NULL) +
    geom_text(aes(0,0,label=metric_pretty_names[[metric]]),size=20*r)
    # theme(plot.title = element_text(size=40*r),
    #       plot.subtitle = element_blank(),
    #       axis.title = element_blank(),
    #       legend.position = "none",
    #       axis.ticks.length = unit(0, "pt"),
    #       panel.background = element_rect(fill = 'white', color = 'white'),
    #       plot.margin = unit(c(0, 0, 0, 0), "points"))
  title
  combined <- wrap_plots(t = title, A = bp[[1]], B =  bp[[2]], C = bp[[3]], D = leg, design = layout)
  # combined <- wrap_plots(t = title, A = bp[[1]], B =  bp[[2]], C = bp[[3]], design = layout2)
  plot((combined + theme(legend.position = "none")))
  dev.off()
}
