cell_asu=$( grep CRYST1 pdb_clean.pdb | awk '{print $2","$3","$4","$5","$6","$7}' )
sg_asu=$( grep CRYST1 pdb_clean.pdb | awk '{sg=substr($0,56,10);gsub(/[ \t]/, "", sg);print sg}' )
echo $cell_asu
echo $sg_asu
if [ -e skip ];then
gmx editconf -f md_restrained.tpr -o md_restrained.pdb
grep -v -e MODEL -e ENDMDL md_restrained.pdb > temp.pdb
mv temp.pdb md_restrained.pdb
pdb4amber -i md_restrained.pdb -o md_restrained_amber.pdb
gmx trjconv -f md_restrained.xtc -s md_restrained.tpr -b 90000 -e 100000 -o md_restrained_last_10ns.xtc<<EOF
0
EOF
mpirun -np 13 python ~/packages/lunus/lunus/command_line/xtraj.py top=md_restrained_amber.pdb traj=md_restrained_last_10ns.xtc first=0 last=250 unit_cell=${cell_asu} space_group=${sg_asu}  fcalc=fcalc_md_restrained_asym.mtz icalc=icalc_md_restrained_asym.mtz diffuse=diffuse_md_restrained_asym.hkl

sel="peptide"
seln="${sel// /_}"
mpirun -np 13 python ~/packages/lunus/lunus/command_line/xtraj.py top=md_restrained_amber.pdb traj=md_restrained_last_10ns.xtc first=0 last=250 unit_cell=${cell_asu} space_group=${sg_asu}  fcalc=fcalc_md_restrained_asym_${seln}.mtz icalc=icalc_md_restrained_asym_${seln}.mtz diffuse=diffuse_md_restrained_asym_${seln}.hkl selection="${sel}"

sel="water"
seln="${sel// /_}"
mpirun -np 13 python ~/packages/lunus/lunus/command_line/xtraj.py top=md_restrained_amber.pdb traj=md_restrained_last_10ns.xtc first=0 last=250 unit_cell=${cell_asu} space_group=${sg_asu}  fcalc=fcalc_md_restrained_asym_${seln}.mtz icalc=icalc_md_restrained_asym_${seln}.mtz diffuse=diffuse_md_restrained_asym_${seln}.hkl selection="${sel}"
fi
sel="resname Na+ or resname Cl-"
seln="${sel// /_}"
mpirun -np 13 python ~/packages/lunus/lunus/command_line/xtraj.py top=md_restrained_amber.pdb traj=md_restrained_last_10ns.xtc first=0 last=250 unit_cell=${cell_asu} space_group=${sg_asu}  fcalc=fcalc_md_restrained_asym_${seln}.mtz icalc=icalc_md_restrained_asym_${seln}.mtz diffuse=diffuse_md_restrained_asym_${seln}.hkl selection="${sel}"

#mpirun -np 26 python ~/packages/lunus/lunus/command_line/xtraj.py top=md_restrained_amber.pdb traj=md_restrained_last_10ns.xtc first=0 last=250 fcalc=fcalc_md_restrained_p1.mtz icalc=icalc_md_restrained_p1.mtz diffuse=diffuse_md_restrained_p1.hkl
