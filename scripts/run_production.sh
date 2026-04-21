#gmx grompp -f restrained.mdp -c md_equil.gro -o md_restrained.tpr -p md_model_posre.top -r md_model.pdb -maxwarn 2 >& grompp_restrained.log
#gmx mdrun -ntmpi 1 -ntomp 26 -nb gpu -bonded gpu -pme gpu -pin on -deffnm md_restrained -v >& md_restrained.log

gmx grompp -f unrestrained.mdp -c md_restrained.gro -o md_unrestrained.tpr -p md_model_posre.top -maxwarn 2 >& grompp_unrestrained.log
gmx mdrun -ntmpi 1 -ntomp 26 -nb gpu -bonded gpu -pme gpu -pin on -deffnm md_unrestrained -v >& md_unrestrained.log

