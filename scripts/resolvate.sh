gmx solvate -cp md_equil.gro -cs waterbox/water_equil.gro -o tmp.pdb >& gmx_solvate.log
maxsol=$( grep "Number of solvent molecules:" gmx_solvate.log | awk '{print int ($NF*0.25)}' )
gmx solvate -cp md_equil.gro -cs waterbox/water_equil.gro -o md_resolv.pdb -maxsol ${maxsol} >& gmx_resolvate.log
cat >> md_model_posre.top <<EOF
WAT ${maxsol}
EOF
gmx grompp -f min.mdp -c md_resolv.pdb -o md_resolv_min.tpr -p md_model_posre.top >& grompp_resolv_min.log
gmx mdrun -ntmpi 8 -ntomp 1 -deffnm md_resolv_min -v

gmx grompp -f equil.mdp -c md_resolv_min.gro -o md_resolv_equil.tpr -p md_model_posre.top -r md_model.pdb -maxwarn 2 >& grompp_resolv_equil.log
gmx mdrun -ntmpi 8 -ntomp 1 -deffnm md_resolv_equil -v
