dir=$1
path=$dir/analysis/figs
cd $path
pwd
python ~/sim2/tumorEvoSim/src/plot_diff.py ../comp.csv 
python ~/sim2/tumorEvoSim/src/plot_corr.py ../comp.csv 
python ~/sim2/tumorEvoSim/src/plot_agevr.py ../comp.csv 
python ~/sim2/tumorEvoSim/src/plot_diversity.py ../comp.csv 
#sh ~/sim2/tumor_stuff/tumorEvoSim/src/animate.sh bloodvtiss