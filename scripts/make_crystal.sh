export WORKDIR=${PWD}
pdb4amber -i prot.pdb -o prot_dry.pdb --dry
grep CRYST1 pdb_clean.pdb > temp.pdb
grep -v -e Na+ -e Cl- -e CRYST1 prot_dry.pdb >> temp.pdb
mv temp.pdb prot_dry.pdb
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
grep CRYST1 prot_dry.pdb | awk '{print substr($0,0,55)"P 1"}' > cryst1_p1.pdb
cp cryst1_p1.pdb temp.pdb
cat prot_dry_cell.pdb >> temp.pdb
mv temp.pdb prot_dry_cell.pdb
if [ -n "$1" ]; then
    if [ "$#" -eq 3 ]; then
	PropPDB -p prot_dry_cell.pdb -o xtal.pdb -ix $1 -iy $2 -iz $3
    else
	PropPDB -p prot_dry_cell.pdb -o xtal.pdb -ix $1 -iy $1 -iz $1
    fi
else
    cp prot_dry_cell.pdb xtal.pdb
fi
