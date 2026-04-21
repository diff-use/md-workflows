#set up restraints
awk 'BEGIN{nhe=0;stop=0}{if(stop==0){if ($0 ~ /GOL/) {nhe=1;print}else{if (nhe==1){idx=$2-1;stop=1}else{print}}}}' pdb_clean.pdb |grep -v JRNL > first_copy.pdb
grep -e ATOM -e HETATM -e TER first_copy.pdb | grep -v GOL > first_copy_prot.pdb
rm -f part??
csplit -z -f part first_copy_prot.pdb '/TER/'+1
for f in part??; do
    pdb4amber -i $f -o ${f}_amber.pdb
    gmx genrestr -fc 209.2 209.2 209.2 -f ${f}_amber.pdb -o posre_$f.itp<<EOF
Protein-H
q
EOF
done
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
molnum=1
cp md_model.top md_model_posre.top
for f in part??; do
    awk -v molnum=${molnum} -v partname=${f} 'BEGIN{cnt=0}{if ($0 ~ /moleculetype/) {cnt=cnt+1; if(cnt == molnum + 1) {print "#ifdef POSRES_"partname"\n#include \"posre_"partname".itp\"\n#endif\n"}}print}' md_model_posre.top > temp.top
    mv temp.top md_model_posre.top
    molnum=$(( molnum + 1 ))
done
gmx grompp -f equil.mdp -c md_min.gro -o md_equil.tpr -p md_model_posre.top -r md_model.pdb >& grompp_equil.log
gmx mdrun -ntmpi 1 -ntomp 26 -deffnm md_equil -v
