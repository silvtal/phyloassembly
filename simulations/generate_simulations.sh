## Input datasets (PCGTABLE and TABLE) come from BacterialCore.py output
## We added "others" manually (OTUs that are not E or P but are in the endpoint (tr12) communities by Goldford et al.
PCGTABLE="simulation_parameters/glucosa_results.txt"

## Original samples (tr0) [Fig. 4, Suppl. Fig. 3]
TABLE="simulation_parameters/original_table.from_biom_0.99.txt"
SAMPLES="sa2 sa10 sa7 sa5 sa1 sa8 sa3 sa11 sa6 sa4 sa9" # ordered X1 to X12 (no X11), check map
OUTPUT="my_neutral_model_glc"

i=0
# we run the script once per sample, once per PCG
for sa in $SAMPLES
do
  i=$(($i+1)) # change "X11" to "X12"
  # Total neutrality:
  Rscript dilgrowth_goldford.R \
   --abuntable $TABLE \
   --sample $sa \
   --dilution 0.008 \
   --no_of_dil 12 \
   --fixation_at 0.90 \
   --grow_step 1 \
   --fix_percentage TRUE \
   --perc 1 \
   --cores 8 \
   --no_of_simulations 800 \
   --outputname X$i \
   --outdir "simulations_analysis/simulation_results/"$OUTPUT

  # Constrained neutrality:
  tail -3 $PCGTABLE | tr -d '"' | tr "\t" "|" | while IFS="|" read -r Core Prevalence Abundance Relative_abundances Min Max Average SD Leaves Taxonomy Leaves_number
  do
    Rscript dilgrowth_goldford.R \
     --abuntable $TABLE \
     --sample $sa \
     --subset $Leaves \
     --dilution 0.008 \
     --no_of_dil 12 \
     --fixation_at 0.90 \
     --grow_step 1 \
     --fix_percentage TRUE \
     --perc $Average \
     --cores 8 \
     --no_of_simulations 800 \
     --outputname $Core"_X"$i \
     --outdir "simulations_analysis/simulation_results/"$OUTPUT"_PERC"
  done
done



## tr1 samples [Fig. 4]
## Originally a very specific script included the PCG/simulnumber/rep's loops. It was used like this:
##  > my_neutral_model_parallel_16_tr7_X2_X6.R --m_ini 2 --cores 16 --simuln 664 --rep 6 # simuls from 664 to 680
## With the following wrapper we do the same but the loops stay outside the call to the script
TABLE="simulation_parameters/tr1_X2_X6_table.from_biom_0.99.txt" # transfer1_X2_X6100percArbol/Tree/0.99/table.from_biom_0.99.txt
OUTPUT="my_neutral_model_glc_tr1"

# we run the script once per sample, once per PCG
for i in "2" "6"
do
  if [ $i == "2" ]; then
    SAMPLES="sa1 sa3 sa6 sa9 sa10 sa11 sa13 sa15" # X2
  else
    SAMPLES="sa2 sa4 sa5 sa7 sa8 sa12 sa14 sa16"  # X6
  fi

  for sa in $SAMPLES
  do
    # Constrained neutrality:
    tail -3 $PCGTABLE | tr -d '"' | tr "\t" "|" | while IFS="|" read -r Core Prevalence Abundance Relative_abundances Min Max Average SD Leaves Taxonomy Leaves_number
    do
      Rscript dilgrowth_goldford.R \
       --abuntable $TABLE \
       --sample $sa \
       --subset $Leaves \
       --dilution 0.008 \
       --no_of_dil 12 \
       --fixation_at 0.90 \
       --grow_step 1 \
       --fix_percentage TRUE \
       --perc $Average \
       --cores 8 \
       --no_of_simulations 16 \
       --outputname $Core"_"$sa \
       --outdir "simulations_analysis/simulation_results/"$OUTPUT"_PERC"
    done
  done
done
