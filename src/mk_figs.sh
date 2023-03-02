INDIR = $1
OUTDIR = $2
python ~/Documents/current/tumor_stuff/tumorEvoSim/src/plot_diff.py ../comp.csv 
python ~/Documents/current/tumor_stuff/tumorEvoSim/src/plot_corr.py ../comp.csv 
python ~/Documents/current/tumor_stuff/tumorEvoSim/src/plot_diversity.py ../comp.csv 
