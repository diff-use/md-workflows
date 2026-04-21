#Uncomment the lines that need to be run

#mkdir -p ligand
#cd ligand
#bash ../run_params_gaussian.sh
#cd -
#bash scripts/param_prot.sh 6B8X
#bash scripts/make_crystal.sh 1
bash scripts/make_waterbox.sh
bash scripts/solvate.sh
bash scripts/minimize.sh
bash scripts/equilibrate.sh
#bash scripts/resolvate.sh
