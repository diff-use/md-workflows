lbl="1_ns_to_11ns"
cell_asu=$( grep CRYST1 pdb_clean.pdb | awk '{print $2","$3","$4","$5","$6","$7}' )
sg_asu=$( grep CRYST1 pdb_clean.pdb | awk '{sg=substr($0,56,10);gsub(/[ \t]/, "", sg);print sg}' )
echo $cell_asu
echo $sg_asu
gmx editconf -f md_restrained.tpr -o md_restrained.pdb
grep -v -e MODEL -e ENDMDL md_restrained.pdb > temp.pdb
mv temp.pdb md_restrained.pdb
pdb4amber -i md_restrained.pdb -o md_restrained_amber.pdb
gmx trjconv -f md_restrained.xtc -s md_restrained.tpr -b 1000 -e 11000 -o md_restrained_${lbl}.xtc<<EOF
0
EOF
mpirun -np 13 python ~/packages/lunus/lunus/command_line/xtraj.py top=md_restrained_amber.pdb traj=md_restrained_${lbl}.xtc first=0 last=250 unit_cell=${cell_asu} space_group=${sg_asu}  fcalc=fcalc_md_restrained_${lbl}_asym.mtz icalc=icalc_md_restrained_${lbl}_asym.mtz diffuse=diffuse_md_restrained_${lbl}_asym.hkl

if [ -e skip ];then
sel="peptide"
seln="${sel// /_}"
mpirun -np 13 python ~/packages/lunus/lunus/command_line/xtraj.py top=md_restrained_amber.pdb traj=md_restrained_last_10ns.xtc first=0 last=250 unit_cell=${cell_asu} space_group=${sg_asu}  fcalc=fcalc_md_restrained_asym_${seln}.mtz icalc=icalc_md_restrained_asym_${seln}.mtz diffuse=diffuse_md_restrained_asym_${seln}.hkl selection="${sel}"

sel="water"
seln="${sel// /_}"
mpirun -np 13 python ~/packages/lunus/lunus/command_line/xtraj.py top=md_restrained_amber.pdb traj=md_restrained_last_10ns.xtc first=0 last=250 unit_cell=${cell_asu} space_group=${sg_asu}  fcalc=fcalc_md_restrained_asym_${seln}.mtz icalc=icalc_md_restrained_asym_${seln}.mtz diffuse=diffuse_md_restrained_asym_${seln}.hkl selection="${sel}"
#fi
sel="resname Na+ or resname Cl-"
seln="${sel// /_}"
mpirun -np 13 python ~/packages/lunus/lunus/command_line/xtraj.py top=md_restrained_amber.pdb traj=md_restrained_last_10ns.xtc first=0 last=250 unit_cell=${cell_asu} space_group=${sg_asu}  fcalc=fcalc_md_restrained_asym_${seln}.mtz icalc=icalc_md_restrained_asym_${seln}.mtz diffuse=diffuse_md_restrained_asym_${seln}.hkl selection="${sel}"
fi
#mpirun -np 26 python ~/packages/lunus/lunus/command_line/xtraj.py top=md_restrained_amber.pdb traj=md_restrained_last_10ns.xtc first=0 last=250 fcalc=fcalc_md_restrained_p1.mtz icalc=icalc_md_restrained_p1.mtz diffuse=diffuse_md_restrained_p1.hkl
