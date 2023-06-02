source("./my_neutral_model_analysis_distr.R")
# From data generated in "my_neutral_model.R". 
# Does ecdf:
#  - simulated vs real groups. 
#  - Diversity, evenness, richness, Faith's PD x4
#  - For each initial sample (X1-X12), x4x 12
#  - for each pcg and for all x4 x12 x4
# For each metric:
#    > create a ECDF distribution plot for the simulated data. See where
#      the real data falls.
# ------------------------------------------------------------------------------
#  --> el "all" usado aquí es "all_final_abund", el csv tal cual salio del script
#      de simulación. !! Ese script puede haber usado porcentajes o no !
#  --> ademas de "all", genera "data_real_" y "data_simul_" para los PCG P, E y O

source("./my_neutral_model_analysis_distr_PERC.R") # aplic porc
# (perc == v.2)
# v.2 : Lee los PCGs separados y MULTIPLICA TODAS SUS ABUNDANCIAS POR EL 
#       PORCENTAJE ASIGNADO A CADA PCG                                              # FIX
#  --> no hace un sample, solamente multiplica. Esto significa que lo único que 
#      se verá afectado serán las abundancias relativas finales, pero no la 
#      richness/supervivencia ni la diversidad.
#! --> genera "ECDF_metric_Xm_all_PERC". 
# --> No genera para PCG por separado.
