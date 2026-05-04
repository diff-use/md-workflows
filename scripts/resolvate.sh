#How the script is called: bash scripts/resolvate.sh

#1st gmx line conducts trial solvation to determine how many water molecules could fit in the voids that opened up during equilibration
#maxsol line takes 25% of that count as maxsol (conservative refill — avoids overpacking)
#2nd gmx line conducts actual re-solvation capped at maxsol molecules
gmx solvate -cp md_equil.gro -cs waterbox/water_equil.gro -o tmp.pdb >& gmx_solvate.log
maxsol=$( grep "Number of solvent molecules:" gmx_solvate.log | awk '{print int ($NF*0.25)}' )
gmx solvate -cp md_equil.gro -cs waterbox/water_equil.gro -o md_resolv.pdb -maxsol ${maxsol} >& gmx_resolvate.log

#Appends the new water count to the topology's "molecules" section
cat >> md_model_posre.top <<EOF
WAT ${maxsol}
EOF

#Second minimization of the re-solvated system. Note the parallelization changes to 8 MPI ranks x 1 OMP thread (domain decomposition, 
#vs. the earlier 1x26 thread-only approach).
gmx grompp -f min.mdp -c md_resolv.pdb -o md_resolv_min.tpr -p md_model_posre.top >& grompp_resolv_min.log
gmx mdrun -ntmpi 8 -ntomp 1 -deffnm md_resolv_min -v

#Second position-restrained equilibration. -maxwarn 2 suppresses warnings from the appended water topology potentially having coordinate 
#mismatches. -r md_model.pdb still uses the original coordinates as restraint reference
gmx grompp -f equil.mdp -c md_resolv_min.gro -o md_resolv_equil.tpr -p md_model_posre.top -r md_model.pdb -maxwarn 2 >& grompp_resolv_equil.log
gmx mdrun -ntmpi 8 -ntomp 1 -deffnm md_resolv_equil -v
