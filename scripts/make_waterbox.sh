mkdir -p waterbox
cd waterbox
grep CRYST1 ../xtal.pdb > cryst1_xtal.pdb
awk -v FIELDWIDTHS="6 9 9 9 7 7 7 10 4" '/^CRYST1/ {printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f\n",$2/10, $3/10, $4/10, $5, $6, $7}' ../xtal.pdb > box.pdb
gmx insert-molecules -f box.pdb -ci ../WAT.pdb -conc 58.0 -o box_solv.pdb >& insert-molecules.log
PropPDB -p box_solv.pdb -o box_solv_expand.pdb -ix 10 -iy 10 -iz 10
grep -v -e CRYST1 -e HEADER box_solv_expand.pdb > temp.pdb
cat cryst1_xtal.pdb temp.pdb > box_solv_expand.pdb
awk '/molecules/{stop=1} stop==0{print}' ../prot.top > waterbox.top
nwat=$( grep "Output configuration contains" insert-molecules.log | awk '{print $7*1000}' )
cat >> waterbox.top <<EOF
[ molecules ]
; Compound       #mols
WAT              $nwat
EOF
gmx grompp -f ../min_water.mdp -c box_solv_expand.pdb -o water_min.tpr -p waterbox.top > grompp_min.log
gmx mdrun -ntmpi 1 -ntomp 26 -deffnm water_min -v
gmx grompp -f ../equil_water.mdp -c water_min.gro -o water_equil.tpr -p waterbox.top > grompp_equil.log
gmx mdrun -ntmpi 1 -ntomp 26 -deffnm water_equil -v
cd -


