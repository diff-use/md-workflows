#How the script is called: bash scripts/make_waterbox.sh

#ank line converts CRYST1 dimensions from Angstroms to nm (divide by 10) to create GROMACS-compatible box def in box.pdb
#gmx line fills box at 58 mol/L
mkdir -p waterbox
cd waterbox
grep CRYST1 ../xtal.pdb > cryst1_xtal.pdb
awk -v FIELDWIDTHS="6 9 9 9 7 7 7 10 4" '/^CRYST1/ {printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f\n",$2/10, $3/10, $4/10, $5, $6, $7}' ../xtal.pdb > box.pdb
gmx insert-molecules -f box.pdb -ci ../WAT.pdb -conc 58.0 -o box_solv.pdb >& insert-molecules.log

#PropPDB line expands water box 10x in each direction, creating a huge reservoir so GROMACS solvation can tile it to fill any crystal geometry
#grep+cat lines replace CRYST1 in the expanded box w/the original crystal cell dimensions (in Angstroms)
#awk line copies the force-field header from prot.top (everything before [ molecules ]) into waterbox.top
#nwat line parses the water molecule count from the GROMACS log and multiplies by 1000 (field 7 is in thousands).
PropPDB -p box_solv.pdb -o box_solv_expand.pdb -ix 10 -iy 10 -iz 10
grep -v -e CRYST1 -e HEADER box_solv_expand.pdb > temp.pdb
cat cryst1_xtal.pdb temp.pdb > box_solv_expand.pdb
awk '/molecules/{stop=1} stop==0{print}' ../prot.top > waterbox.top
nwat=$( grep "Output configuration contains" insert-molecules.log | awk '{print $7*1000}' )

#Appends the [ molecules ] section with the water count.
#Minimizes the water box (min_water.mdp), then equilibrates it (equil_water.mdp) — two sequential GROMACS simulations. The result 
#water_equil.gro is a thermally equilibrated pure-water reservoir.
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


