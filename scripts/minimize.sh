gmx grompp -f min.mdp -c md_model.pdb -o md_min.tpr -p md_model.top >& grompp_min.log
gmx mdrun -ntmpi 1 -ntomp 26 -deffnm md_min -v
