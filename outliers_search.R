## Script que te genera tres diagramas de barras, uno para cada PCG (y noPCG)
## - En el eje X van los nombres con especie de todas las OTUs
## - En el eje Y hay dos medidas. A la izquierda, la abundancia relativa final
##   (FRA) esperada (simul), a la derecha, la FRA obtenida (real)
## - Se dibujan una barra y una flecha desde la FRA esperada hasta la real. Se 
##   colorea verde (y la flecha va hacia arriba) si la FRA ha subido; se colorea
##   gris (y la flecha va hacia abajo) si la FRA ha bajado.
## - Si ha "bajado a 0" pero es porque no se había detectado en las reads a
##   a tiempo 0, entonces sale en negro, no en gris.
## - TODO: poner significancia estadística.
# ==============================================================================
## v3
# - añado todos los puntos (real_final_dots[[m]][[PCG]][otu,])
# - añado DTS (= T-test pero mejor) 
      # <- TODO hmmm tarda mucho pero me sale casi todo ultrasignificativo. cambio de test a t test o algo?
## ------------------
## v4
# a) numeritos rojos ( # TODO )
# b) quitar others (HECHO)
# c) solo 8 y 9, y juntas (HECHO!)

## A mano mejor...
## real_final_abund[[m]][[PCG]] <- exp[map[map["ORIG"]==m & map["TIME"]==my_tr
##                                        ,][,"SA"]][pcgc %in% c(1, 2),]
## real_final_abund[[m]][[PCG]]

## X8
##          sa17  sa24 sa35   sa5  sa54  sa62 sa77 sa93
# 4399988   504 10690   14    40 10539    28 6040 9016
# 4454257 12362   227 9637 14018  2892 15486 1680 5766
# 4327501     0     0    0     0     0     0    0    0
# 4456892     0     0    0     0     0     8    0    0
# 4461938     0     0    0     0     0     0    0    0
# 4455861     0     0    0     0     0     0    0    0
# 4451011     0    15    0   258    56    62    0   24
# 4456889    36  2695 5684  1542  1280   248 4367  234

## X9
# sa101  sa49  sa51  sa67   sa8  sa80  sa96
# 4399988    37     0     0     0     0     0     0
# 4454257 13100 14174 11396 13804 12242 13704 15034
# 4327501     0     0     0     0     0     0     0
# 4456892     5     0     0     0     0     0     0
# 4461938     0     0     0     0     0     0     0
# 4455861     0     0     0     0     0     0     0
# 4451011  2674  1697  4508  1855  3621  2145   851
# 4456889    29     0     0     5    13     0     0
# ==============================================================================
# Libraries
source("./simulations/my_functions.R")
library("patchwork")
library("ggrepel")
library("ggnewscale") # source https://elartedeldato.com/blog/como-asignar-dos-escalas-diferentes-para-una-misma-estetica-en-ggplot/
library("twosamples")

# input simul files folder; we assume default naming for files inside
input_s="./simulations/simulation_analysis/simulation_results/my_neutral_model_glc_PERC/"

# Cargar datos reales (no hay que cambiarlo, es para ver los datos finales solamente)
exp_f <- "./simulations/simulation_parameters/all_glc.txt"
map_f <- "./simulations/simulation_parameters/maps/map_all_glc.csv"

# output
outfile = "Supplementary Figure 3.png"

# params
my_tr  <- "final" # studied transfer; should normally refer to the final transfer
signif <- FALSE

# ==============================================================================
# prepare real data reading
map <- read.csv(map_f,sep=",")
m_inic <- unique(map["ORIG"])[["ORIG"]]
# m_inic <- m_inic[m_inic != "X11"] # i need to remove the X11 sample because 
#                                  # i couldn't simulate those
# m_inic <- c("X2", "X6") #tr1 # TODO
m_inic <- c("X8", "X9") # v4 # Suppl Fig 4

list[exp,tax]<-get_abundance_and_tax_from_table(exp_f)

# PCG table to classify the OTUs correctly
pcg_table <- read.csv("./simulations/simulation_parameters/pcg_leaves_tr12_orig.txt", sep ="\t")[c("Core", "Leaves")]
otu_to_core <- pcg_table %>% separate_rows(Leaves, sep = ";") %>% as.data.frame()
rownames(otu_to_core) <- otu_to_core$Leaves
pcgc <- c(2, 1, 3)[otu_to_core[rownames(tax), ]$Core %>% as.factor] # accurate OTU classification
names(pcgc) <- rownames(tax)

real_final_abund <- list()

# prepare simul data reading
simul_final_abund=list()
PCG_name <- c("Node35562","Node27828") #v4,"others")

# start data loop
for (m in m_inic){
  simul_final_abund[[m]] <- list()
  real_final_abund[[m]] <- list()
  
  my_dat <- list()
  for (pcg in 1:length(PCG_name)) {
    PCG <- PCG_name[pcg]
    ## load real data (per m and per PCG)
    real_final_abund[[m]][[PCG]] <- exp[map[map["ORIG"]==m & map["TIME"]==my_tr
                                            ,][,"SA"]][pcgc == pcg,]
    "./simulations/simulation_analysis/simulation_results/my_neutral_model_glc_PERC/"
    ## load simul data (per m and per PCG)
    if (file.exists(paste0(input_s,"/simul_",PCG,"_",m,".csv"))) {
      simul_final_abund[[m]][[PCG]] <- my_transpose(read.csv(paste0(input_s,"/simul_",PCG,"_",m,".csv"),
                                                             row.names = NULL, check.names = F)[,-1,drop=FALSE])
      # TODO / DEBUG
      # This step would normally not be necessary. Basically, it excludes non-Core
      # OTUs that are included in Core simul files. The excluded ones won't be
      # added to "others" automatically.
      core.otus <- rownames(simul_final_abund[[m]][[PCG]])[pcgc[rownames(simul_final_abund[[m]][[PCG]])]==pcg] 
      non.core  <- rownames(simul_final_abund[[m]][[PCG]])[pcgc[rownames(simul_final_abund[[m]][[PCG]])]!=pcg] 
      warning(paste0(paste(non.core, collapse = ", "), " were incorrectly included in ", PCG, ". Excluded."))
      simul_final_abund[[m]][[PCG]] <- simul_final_abund[[m]][[PCG]][core.otus, ]

      limit <- ncol(simul_final_abund[[m]][[PCG]])
      
      if (limit < 800) {
          simul_final_abund[[m]][[PCG]][(limit+1):800] <-  0
          colnames(simul_final_abund[[m]][[PCG]]) <- str_replace_all(colnames(simul_final_abund[[m]][[PCG]]), pattern = "V", replacement = "")
          message(paste0("A total of ", 800, " simulations will be used for the null model analysis for all PCGs together (sample ",m,"). ", 800-limit, " empty columns were added."))
          limit <- 800
        } else {
          message(paste0("A total of ",limit," simulations will be used for the null model analysis for all PCGs together (sample ",m,")."))
        }
    } else {
      message(paste0("PCG ",PCG," of sample ",m, " does not have an existing file, probably disappeared during simulations. Skipped."))
    }
  }
}


# switch to relative abundances
for (m in m_inic){
  real_all_pcgs  <- rbind(real_final_abund[[m]][[PCG_name[1]]],
                          real_final_abund[[m]][[PCG_name[2]]],
                          real_final_abund[[m]][[PCG_name[3]]])
  simul_all_pcgs <- rbind(simul_final_abund[[m]][[PCG_name[1]]],
                          simul_final_abund[[m]][[PCG_name[2]]],
                          simul_final_abund[[m]][[PCG_name[3]]])
  for (PCG in PCG_name) {
    real_final_abund [[m]][[PCG]] <- real_final_abund [[m]][[PCG]] / colSums(real_all_pcgs)
    simul_final_abund[[m]][[PCG]] <- simul_final_abund[[m]][[PCG]] / colSums(simul_all_pcgs)
  }
}


#### v3:
real_final_dots  <- real_final_abund
simul_final_dots <- simul_final_abund

undetected <- list()
plots <- list()
my_dat <- list()

layout <- "
  AAAAAAAAsCCCCCCCC
  BBBBBBBBsDDDDDDDD
"
  
# start plotting loop
for (m in m_inic){
  undetected[[m]] <- list()
  plots[[m]] <- list()

  my_dat[[m]] <- list()
  for (pcg in 1:length(PCG_name)) {
    PCG <- PCG_name[pcg]
    #### v3:
    nonzerootus <- unique(c(rownames(real_final_abund [[m]][[PCG]][rowSums(real_final_abund [[m]][[PCG]])!=0,]),
                           rownames(simul_final_abund[[m]][[PCG]][rowSums(simul_final_abund[[m]][[PCG]])!=0,])))
    # remove upcoming NA values
    undetected[[m]][[PCG]] <- rownames(real_final_abund[[m]][[PCG]]) [rownames(real_final_abund[[m]][[PCG]]) %!in% rownames(simul_final_abund[[m]][[PCG]])]
    
    simul_final_abund[[m]][[PCG]] <- simul_final_abund[[m]][[PCG]][nonzerootus,]
    rownames(simul_final_abund[[m]][[PCG]]) <- nonzerootus
    simul_final_abund[[m]][[PCG]][is.na(simul_final_abund[[m]][[PCG]])] <- 0 # undetected
    
    real_final_abund [[m]][[PCG]] <- real_final_abund [[m]][[PCG]][nonzerootus,]
    real_final_abund [[m]][[PCG]][is.na(real_final_abund [[m]][[PCG]])] <- 0 # 0 in all m
    rownames(real_final_abund[[m]][[PCG]]) <- nonzerootus
    
    tbl_colnames <- c("real", "simul", "type")
    # create empty one (NaNs)
    my_dat[[m]][[PCG]] <- as.data.frame(matrix(nrow = length(nonzerootus),
                                          ncol = length(tbl_colnames),
                                          dimnames = list(nonzerootus, tbl_colnames)))
    
    space <- 20/length(nonzerootus) # between bars
    xmin  <-  0                # position in the plot
    for (otu in nonzerootus) {
      ### actual data
      my_dat[[m]][[PCG]][otu, "real"]    <- mean(as.numeric(real_final_abund[[m]][[PCG]][otu,]))
      # TODO OJO: varía mucho de una réplica a otra          
      #          sa16 sa25 sa36  sa4 sa53 sa61 sa78 sa94
      # 4399988 11744   28  117 2614   85 2402  605 10756
      my_dat[[m]][[PCG]][otu, "simul"]   <- mean(as.numeric(simul_final_abund[[m]][[PCG]][otu,]))
      my_dat[[m]][[PCG]][otu, "type"]    <- ifelse(test = otu %in% undetected[[m]][[PCG]],
                                              yes = "undetected",
                                              no = c("higher than predicted", "lower than predicted")[1 + as.numeric(my_dat[[m]][[PCG]][otu, "real"] <  my_dat[[m]][[PCG]][otu, "simul"])])
      #### v3:
      ### coords
      xmin <- xmin + space
      my_dat[[m]][[PCG]][otu, "xmin"] <- xmin
      my_dat[[m]][[PCG]][otu, "xmax"] <- xmin + space/1.5
      real_final_dots[[m]][[PCG]][otu, "xmin"] <- xmin
      real_final_dots[[m]][[PCG]][otu, "xmax"] <- xmin + space/1.5
      simul_final_dots[[m]][[PCG]][otu, "xmin"] <- xmin
      simul_final_dots[[m]][[PCG]][otu, "xmax"] <- xmin + space/1.5
      my_dat[[m]][[PCG]][otu, "OTU"]  <- otu
      my_dat[[m]][[PCG]][otu, "PCG"]  <- PCG
      my_dat[[m]][[PCG]][otu, "ymin"] <- min(my_dat[[m]][[PCG]][otu, "real"],
                                        my_dat[[m]][[PCG]][otu, "simul"])
      my_dat[[m]][[PCG]][otu, "ymax"] <- max(my_dat[[m]][[PCG]][otu, "real"],
                                        my_dat[[m]][[PCG]][otu, "simul"])
      my_dat[[m]][[PCG]][otu, "arrow_start"] <- ifelse(test = my_dat[[m]][[PCG]][otu, "type"] == "lower than predicted", 
                                                  yes  = my_dat[[m]][[PCG]][otu, "ymax"],
                                                  no   = my_dat[[m]][[PCG]][otu, "ymin"])
      my_dat[[m]][[PCG]][otu, "arrow_end"]   <- ifelse(test = my_dat[[m]][[PCG]][otu, "type"] == "lower than predicted", 
                                                  yes  = my_dat[[m]][[PCG]][otu, "ymin"],
                                                  no   = my_dat[[m]][[PCG]][otu, "ymax"])
      # my_dat[[m]][[PCG]][otu, "num"]   <- ifelse(test = otu %in% c("4456889", "4454257"),  # TODO
      #                                             yes  = my_dat[[m]][[PCG]][otu, "ymin"],
      #                                             no   = my_dat[[m]][[PCG]][otu, "ymax"])
      ### dts calc
      if (signif) {
        pval <- two_sample(unlist(simul_final_abund[[m]][[PCG]]),
                         unlist(real_final_abund [[m]][[PCG]]),
                         nboots=800)

        if (pval[[2]] <= 0.001) {
          my_dat[[m]][[PCG]][otu, "significance"]  <- "***"
        } else if (pval[[2]] <= 0.01) {
          my_dat[[m]][[PCG]][otu, "significance"]  <- "**"
        } else if (pval[[2]] <= 0.05) {
          my_dat[[m]][[PCG]][otu, "significance"]  <- "*"
        } else {
          my_dat[[m]][[PCG]][otu, "significance"]  <- ""
        }
      }
    }
    
    plots[[m]][[PCG]] <- ggplot() +
      ylim(0, max(real_final_abund[[m]][[PCG]], simul_final_abund[[m]][[PCG]]))
    
    #### v3:
    # Luego reshapeo el df (simul, pintadas encima, importante)
    dfseg <- NULL
    sas <- colnames(simul_final_dots[[m]][[PCG]])
    for (sa in sas[1:(length(sas)-2)]) {
      if (is.null(dfseg)){
        dfseg <- data.frame(simul_final_dots[[m]][[PCG]]["xmin"], 
                            simul_final_dots[[m]][[PCG]]["xmax"],
                            simul_final_dots[[m]][[PCG]][sa])
        colnames(dfseg) <- c("xmin","xmax", "SA")
        dfseg["Replicate"] <- "simul"
        dfseg["size"] <- 0.1 # "simul"
      } else {
        new <- data.frame(simul_final_dots[[m]][[PCG]]["xmin"], 
                          simul_final_dots[[m]][[PCG]]["xmax"],
                          simul_final_dots[[m]][[PCG]][sa])
        colnames(new) <- c("xmin","xmax", "SA")
        new["Replicate"] <- "simul"
        new["size"] <- 0.1 # "simul"
        dfseg <- rbind(dfseg, new)
      }
    } 
    
    # Reshapeo el df (real, pinto encima)
    sas <- colnames(real_final_dots[[m]][[PCG]])
    for (sa in sas[1:(length(sas)-2)]) {
      new <- data.frame(real_final_dots[[m]][[PCG]]["xmin"], 
                        real_final_dots[[m]][[PCG]]["xmax"],
                        real_final_dots[[m]][[PCG]][sa])
      colnames(new) <- c("xmin","xmax", "SA")
      new["Replicate"] <- "real"
      new["size"] <- 0.3 # "real"
      dfseg <- rbind(dfseg, new)
    }
    
    # Pinto líneas horizontales
    plots[[m]][[PCG]] <- plots[[m]][[PCG]] + 
      geom_segment(data = dfseg,
                   aes(x = xmin,
                       xend =xmax,
                       y = SA,
                       yend = SA,
                       color = Replicate),
                   size = dfseg$size
      ) +
      scale_color_manual("Replicate",
                         values = c("real"  = "red",
                                    "simul" = "black")
      )

    # Ahora las flechas
    plots[[m]][[PCG]] <- plots[[m]][[PCG]] + 
      scale_y_continuous(name="Relative abundance") +
      scale_x_continuous(name="OTUs") +
      ggtitle(label = PCG) +
      # geom_rect(data = my_dat[[m]][[PCG]],
      #           aes(xmin = xmin,
      #               xmax = xmax,
      #               ymin = ymin,
      #               ymax = ymax,
      #               fill = type)
      # ) +
      new_scale_color() +
      geom_segment(data = my_dat[[m]][[PCG]], 
                   aes(x = rowMeans(cbind(xmin, xmax)),
                       xend = rowMeans(cbind(xmin, xmax)),  
                       y = arrow_start,
                       yend = arrow_end,
                       ),
                   color = "black",
                   size = (my_dat[[m]][[PCG]]$xmax - my_dat[[m]][[PCG]]$xmin) / 1.4,
                   arrow = arrow(angle  = 25, 
                                 type   = "closed",
                                 length = unit(2, "mm"),
                   ),
                   inherit.aes = FALSE
      ) +
      geom_segment(data = my_dat[[m]][[PCG]], 
                   aes(x = rowMeans(cbind(xmin, xmax)),
                       xend = rowMeans(cbind(xmin, xmax)),  
                       y = arrow_start,
                       yend = arrow_end,
                       color = type),
                   size = (my_dat[[m]][[PCG]]$xmax - my_dat[[m]][[PCG]]$xmin) / 2,
                   arrow = arrow(angle  = 25, 
                                 type   = "closed",
                                 length = unit(2, "mm")
                                 ),
                   inherit.aes = FALSE
      ) +
      scale_color_manual('type',
                         values = c("undetected" = "black",
                                    "lower than predicted"  = "darkgrey",
                                    "higher than predicted" = "darkgreen"),  
                         guide = guide_legend(override.aes = list(alpha = 1))
      ) + 
      # Ahora los OTU names (1)
      geom_label_repel(data = my_dat[[m]][[PCG]],
                      aes(x = xmax , #- space,
                          y = ymax),
                      label = my_dat[[m]][[PCG]]$OTU,
                      size  = 3,
                      direction = "y",
                      # nudge_x = + space,
                      force = 2,
                      max.time = 1,
                      ylim = c(0, 10),
                      max.overlaps = ifelse(PCG == "others", yes = 25, no = 1000),
                      label.padding = 0.1,
                      alpha = 0.6,
                      label.size = NA,
                      seed = 1234
      ) + 
      # Ahora los OTU names (2)
      geom_label_repel(data = my_dat[[m]][[PCG]],
                       aes(x = xmax , #- space,
                           y = ymax),
                       label = my_dat[[m]][[PCG]]$OTU,
                       size  = 3,
                       direction = "y",
                       # nudge_x = + space,
                       force = 2,
                       max.time = 1,
                       ylim = c(0, 10),
                       max.overlaps = ifelse(PCG == "others", yes = 25, no = 1000),
                       label.padding = 0.1,
                       alpha = 1,
                       fill = NA,
                       label.size = NA,
                       seed = 1234
      # ) + 
      # # Ahora los numeritos rojos
      # geom_label_repel(data = my_dat[[m]][[num]],
      #                  aes(x = dfseg$xmin - 0.5,              # TODO
      #                      y = dfseg$SA),
      #                  color = "red"
      )
    
    # Ahora la significación (opcional)
    if (signif) {
      plots[[m]][[PCG]] <- plots[[m]][[PCG]] + 
        geom_text(data = my_dat[[m]][[PCG]],
                  aes(x = xmax,
                      y = ymin + (ymax+ymin)/20),
                  label = my_dat[[m]][[PCG]]$significance,
                  size  = 4
        )
    }
    ## Tema gráfico
    plots[[m]][[PCG]] <- plots[[m]][[PCG]] +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            plot.title = element_text(face = "bold", size=15),
            legend.position = "bottom",
            legend.text = element_text(size=14),
            # plot.background = element_rect("grey95"),
            panel.background = element_rect("white"),
            panel.grid = element_line(colour = "grey90")
      )
  }
}

  ## put together and save the plot
  combined <- wrap_plots(guides = "collect",
                         A = plots[["X8"]][[PCG_name[1]]], B = plots[["X8"]][[PCG_name[2]]],
                         C = plots[["X9"]][[PCG_name[1]]], D = plots[["X9"]][[PCG_name[2]]],
                         #  v4 C = plots[[PCG_name[3]]],
                         design = layout) &
    theme(legend.position = 'bottom') &
    plot_annotation(title = "\tX8\t\t\t\t\t\t\t\t\t\t\t\t\t\tX9",
                    theme = theme(plot.title = element_text(face = "bold", size=20, hjust = 0.5))) 
  
  png(file = paste0(outfile),
      width = 40, #21*2, # v4,
      height = 25, # v4 429.7,
      units = "cm",
      res = 480)
  
  plot(combined)
  dev.off()