#How the script is called: bash scripts/make_crystal.sh 1
#You can run this with $1 $2 $3 for supercell dimensions

# pdb4amber --dry removes water molecules and saves as prot_dry.pdb
# rest of the lines ... grab CRYST1 record from the clean PDB file (to preserve unit cell) and prepends it to prot_dry.pdb
# --> removes stray ions and any wrong CRYST1 that pdb4amber may have written
export WORKDIR=${PWD}
pdb4amber -i prot.pdb -o prot_dry.pdb --dry
grep CRYST1 pdb_clean.pdb > temp.pdb
grep -v -e Na+ -e Cl- -e CRYST1 prot_dry.pdb >> temp.pdb
mv temp.pdb prot_dry.pdb

#Chimera (used in headless mode) exopands the asymetric unit into the full unit 
#cell using crystallographic symmetry operators from the PDB header. "unitcell" 
#generates all symmetry mates, combine merges them into one model, and saves as prot_dry_cell.pdb
CHIMERA_EXEC=/Applications/ChimeraX-1.10.app/Contents/bin/ChimeraX
cat > expand.cxc <<EOF
open ${WORKDIR}/prot_dry.pdb
changechains #1 A
unitcell #1
combine #2
save ${WORKDIR}/prot_dry_cell.pdb #3
quit
EOF
${CHIMERA_EXEC} --nogui expand.cxc

#Rewrites the CRYST1 spacegroup to P 1 (since symmetry has already been expanded — all molecules are 
#now explicit atoms). Prepends this to the cell PDB.
grep CRYST1 prot_dry.pdb | awk '{print substr($0,0,55)"P 1"}' > cryst1_p1.pdb
cp cryst1_p1.pdb temp.pdb
cat prot_dry_cell.pdb >> temp.pdb
mv temp.pdb prot_dry_cell.pdb

#PropPDB replicates the P1 cell into a supercell. If 3 args are given, each axis gets 
#its own count; if 1 arg, it's applied to all three. If no args, just copies (1x1x1). For run_all.sh 
#this is called with 1, producing a 1x1x1 supercell.
if [ -n "$1" ]; then
    if [ "$#" -eq 3 ]; then
	PropPDB -p prot_dry_cell.pdb -o xtal.pdb -ix $1 -iy $2 -iz $3
    else
	PropPDB -p prot_dry_cell.pdb -o xtal.pdb -ix $1 -iy $1 -iz $1
    fi
else
    cp prot_dry_cell.pdb xtal.pdb
fi
