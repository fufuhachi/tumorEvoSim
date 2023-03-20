dir=$1
path=$dir/analysis/figs
cd $path
pwd
python ~/Documents/current/tumor_stuff/tumorEvoSim/src/plot_diff.py ../comp.csv 
python ~/Documents/current/tumor_stuff/tumorEvoSim/src/plot_corr.py ../comp.csv 
python ~/Documents/current/tumor_stuff/tumorEvoSim/src/plot_diversity.py ../comp.csv 
sh ~/Documents/current/tumor_stuff/tumorEvoSim/src/animate.sh bloodvtiss
#python ~/Documents/current/tumor_stuff/tumorEvoSim/src/tumor-timeplot.py ../comp.csv 
#sh ~/Documents/current/tumor_stuff/tumorEvoSim/src/animate.sh tumor-timeplot