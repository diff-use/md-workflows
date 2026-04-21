PDBID=$1
grep -v -e REMARK -e KEYWDS ${PDBID}.pdb | grep -e CRYST1 -e ATOM -e HETATM -e TER -e END > pdb_clean.pdb 
#awk '{if (substr($0, 13, 1)=="D") {print substr($0,0,12)"H"substr($0,14,64)"H"} else print $0}' | \
#awk '{if (substr($0, 13, 2)==" D") {print substr($0,0,13)"H"substr($0,15,63)"H"} else print $0}' > pdb_clean.pdb 
pdb4amber -i pdb_clean.pdb --prot -o pdb_clean_amber.pdb
#AddToBox -c pdb_clean_amber.pdb -a spce.pdb -na 1 -o pdb_clean_amber_spce.pdb -RP 3.0 -RW 3.0 -G 0.2 -V 1
#sed -i .bak 's/tip3p/spce/' leap.template.in
#sed -i .bak 's/pdb_clean_amber/pdb_clean_amber_spce/' leap.template.in
#sed -i .bak 's/set\ default\ PBradii/addions2 x Na+ 1\naddions2 x Cl- 1\nset\ default\ PBradii/' leap.template.in

#initial model building to obtain solvent .pdb files
cat > tleap_temp.in <<EOF
source leaprc.protein.ff19SB
source leaprc.DNA.OL15
source leaprc.RNA.OL3
source leaprc.water.spceb
source leaprc.gaff2
p = loadpdb pdb_clean_amber.pdb
x = combine{p}
addions2 x Cl- 0
addions2 x Na+ 0
addions2 x Na+ 1
addions2 x Cl- 1
solvateBox x SPCBOX 1.
set default PBradii mbondi3
set default nocenter on
saveAmberParm x prot.parm7 prot.rst7
quit
EOF
tleap -f tleap_temp.in
rm -f prot.top prot.pdb
cat > amber_to_gromacs.py <<EOF
import parmed as pmd
parm = pmd.load_file("prot.parm7", "prot.rst7")
parm.save("prot.top") # .top saves to GROMACS topology file format
parm.save("prot.pdb")
EOF
python amber_to_gromacs.py
grep HETATM prot.pdb | grep Na+ | head -1 > Na+.pdb
grep HETATM prot.pdb | grep Cl- | head -1 > Cl-.pdb
grep HETATM prot.pdb | grep WAT | head -3 > WAT.pdb

cat > tleap_prot.in <<EOF
source leaprc.protein.ff19SB
source leaprc.DNA.OL15
source leaprc.RNA.OL3
source leaprc.water.spceb
source leaprc.gaff2
p = loadpdb pdb_clean_amber.pdb
w = loadpdb WAT.pdb
x = combine{p w}
addions2 x Na+ 0
addions2 x Cl- 0
addions2 x Na+ 1
addions2 x Cl- 1
set default PBradii mbondi3
set default nocenter on
saveAmberParm x prot.parm7 prot.rst7
quit
EOF
tleap -f tleap_prot.in
rm -f prot.top prot.pdb
python amber_to_gromacs.py
