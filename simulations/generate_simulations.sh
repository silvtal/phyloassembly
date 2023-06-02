## Input abundance tables (TABLE) come from BacterialCore.py output
## the PCGTABLEs were generated with make_pcg_leaves.R

## Original samples (tr0) [Fig. 4, Suppl. Fig. 3]
TABLE="simulation_parameters/original_table.from_biom_0.99.txt"
PCGTABLE="simulation_parameters/pcg_leaves.txt"
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
   --outdir "simulation_analysis/simulation_results/"$OUTPUT

  # Constrained neutrality:
  tail -3 $PCGTABLE | tr -d '"' | tr "\t" "|" | while IFS="|" read -r Core Leaves Percentage
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
     --perc $Percentage \
     --cores 8 \
     --no_of_simulations 800 \
     --outputname $Core"_X"$i \
     --outdir "simulation_analysis/simulation_results/"$OUTPUT"_PERC"
  done
done



## tr1 samples [Fig. 4]
## Originally a very specific script included the PCG/simulnumber/rep's loops. It was used like this:
##  > my_neutral_model_parallel_16_tr7_X2_X6.R --m_ini 2 --cores 16 --simuln 664 --rep 6 # simuls from 664 to 680
## With the following wrapper we do the same but the loops stay outside the call to the script
TABLE="simulation_parameters/tr1_X2_X6_table.from_biom_0.99.txt" # transfer1_X2_X6100percArbol/Tree/0.99/table.from_biom_0.99.txt
PCGTABLE="simulation_parameters/pcg_leaves_tr1.txt"
OUTPUT="my_neutral_model_glc_tr1"

# we run the script once per sample, once per PCG
# we want 800 simulations and have 8 initial samples. So we get 100 simulations from each sample
# and put them together afterwards 
for i in "2" "6"
do
    if [ $i = "2" ]; then
      SAMPLES="sa1 sa3 sa6 sa9 sa10 sa11 sa13 sa15" # X2
    else
      SAMPLES="sa2 sa4 sa5 sa7 sa8 sa12 sa14 sa16"  # X6
    fi

    # Constrained neutrality:
    tail -3 $PCGTABLE | tr -d '"' | tr "\t" "|" | while IFS="|" read -r Core Leaves Percentage
    do
      for sa in $SAMPLES
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
       --perc $Percentage \
       --cores 8 \
       --no_of_simulations 100 \
       --outputname $Core"_X"$i"_"$sa \
       --outdir "simulation_analysis/simulation_results/"$OUTPUT"_PERC"
       
      # Add that replicate to the common pool. Warning: there's many replicates in X6 than don't 
      # contain "others" or Pseudo. For those, I just create a "dummy" file with one OTU column
      # that is all 0s, with length 100
      FILE="simulation_analysis/simulation_results/"$OUTPUT"_PERC/simul_"$Core"_X"$i"_"$sa".csv"
      if [ -f $FILE ]; then
        echo "simul_"$Core"_X"$i"_"$sa".csv was saved"
      else
        echo "simul_"$Core"_X"$i"_"$sa".csv could not be generated; saving dummy file"
        # save dummy file instead (empty column, not conserved in the final output)
        cp empty.csv $FILE
      fi
      done
      Rscript merge_tr1_simuls.R "simulation_analysis/simulation_results/"$OUTPUT"_PERC/simul_"$Core"_X"$i".csv" "simulation_analysis/simulation_results/"$OUTPUT"_PERC/" $Core"_X"$i
    done
done

rm "simulation_analysis/simulation_results/"$OUTPUT"_PERC/"*sa*

##create dummy file
##echo "\"\",\"227343\"" > empty.csv
##for ((i=1; i<=100; i++))
##do
##    echo "\"$i\",0" >> empty.csv
##done
