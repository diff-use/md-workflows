#How the script is called: bash scripts/minimize.sh

#"gmx grompp" line preprocesses the topology + structure + MDP settings into a single binary run-input file (md_min.tpr)
#"gmx mdrun" line runs the minimization (steepest descent or conjugate gradient, as defined in min.mdp). Uses 1 MPI rank x 26 OpenMP threads. Output is md_min.gro
gmx grompp -f min.mdp -c md_model.pdb -o md_min.tpr -p md_model.top >& grompp_min.log
gmx mdrun -ntmpi 1 -ntomp 26 -deffnm md_min -v
