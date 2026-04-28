#How the script is called: bash scripts/equilibrate.sh

#set up restraints

#awk line extracts the first copy of the asymmetric unit from pdb_clean.pdb. It reads lines until it hits a GOL (glycerol) residue, includes that line, then stops at the next atom. This isolates copy 1
#grep line Keeps only ATOM/HETATM/TER, removes any GOL lines -> first_copy_prot.pdb
#rm and csplit lines splits the protein at every TER record into files part00, part01, etc. (one per chain).
awk 'BEGIN{nhe=0;stop=0}{if(stop==0){if ($0 ~ /GOL/) {nhe=1;print}else{if (nhe==1){idx=$2-1;stop=1}else{print}}}}' pdb_clean.pdb |grep -v JRNL > first_copy.pdb
grep -e ATOM -e HETATM -e TER first_copy.pdb | grep -v GOL > first_copy_prot.pdb
rm -f part??
csplit -z -f part first_copy_prot.pdb '/TER/'+1

#For each chain fragment: run pdb4amber to fix naming, then gmx genrestr to generate a position restraint .itp file with force constant 209.2 kJ/mol/nm² (≈ 50 
#kcal/mol/A²) on all non-hydrogen protein atoms (Protein-H group).
for f in part??; do
    pdb4amber -i $f -o ${f}_amber.pdb
    gmx genrestr -fc 209.2 209.2 209.2 -f ${f}_amber.pdb -o posre_$f.itp<<EOF
Protein-H
q
EOF
done

#Guarded by if [ -e skip ] — effectively disabled unless a file called skip exists. Would extract AR6 ligand atoms, create 
#a custom index group for non-hydrogen ligand atoms (AR6-H), and generate restraints for them.
if [ -e skip ]; then
np=$( ls -1 part?? | wc -l)
f=$( printf "part%02d" $np )
grep -e ATOM -e HETATM first_copy.pdb | grep AR6 > $f
pdb4amber -i $f -o ${f}_amber.pdb
gmx make_ndx -f ${f}_amber.pdb -o ${f}_amber.ndx <<EOF
! a H*
name 3 AR6-H
q
EOF
gmx genrestr -fc 209.2 209.2 209.2 -f ${f}_amber.pdb -o posre_$f.itp -n ${f}_amber.ndx<<EOF
AR6-H
q
EOF
fi

#Copies md_model.top -> md_model_posre.top.
#For each chain part, inserts a conditional #ifdef POSRES_partXX / #include "posre_partXX.itp" / #endif block before the next [ moleculetype ] section. This 
#lets you selectively enable restraints per chain at compile time via -D POSRES_part00 etc.
molnum=1
cp md_model.top md_model_posre.top
for f in part??; do
    awk -v molnum=${molnum} -v partname=${f} 'BEGIN{cnt=0}{if ($0 ~ /moleculetype/) {cnt=cnt+1; if(cnt == molnum + 1) {print "#ifdef POSRES_"partname"\n#include \"posre_"partname".itp\"\n#endif\n"}}print}' md_model_posre.top > temp.top
    mv temp.top md_model_posre.top
    molnum=$(( molnum + 1 ))
done

#"gmx grompp" uses -r md_model.pdb as the reference coordinates for position restraints (atoms are pulled back toward their original positions in md_model.pdb).
#"gmx mdrun" line runs the restrained NPT equilibration defined in equil.mdp
gmx grompp -f equil.mdp -c md_min.gro -o md_equil.tpr -p md_model_posre.top -r md_model.pdb >& grompp_equil.log
gmx mdrun -ntmpi 1 -ntomp 26 -deffnm md_equil -v
