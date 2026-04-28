#How the script is called: bash scripts/solvate.sh

#nats_one to ncopies lines counts atoms in the aysm unit vs the full crystal to determine how many symmetry copies exist
topcopy="system1              1"$'\n'
nats_one=$( grep -e ATOM -e HETATM prot_dry.pdb | wc -l )
nats_cell=$( grep -e ATOM -e HETATM xtal.pdb | wc -l )
ncopies=$( echo "$nats_cell/$nats_one" | bc )
echo "There are ${ncopies} copies of the asymmetric unit"

echo "${topcopy}"
topall=""
for i in $( seq 1 $ncopies )
do
    topall=${topall}${topcopy}
done

#gmx line fills crystal's void spaces w/equilibrated water molecules
#nwat line parse the resulting water count
#awk line copies the force-field params (everything above [ molecules ]) from prot.top into md_model.top)
gmx solvate -cp xtal.pdb -cs waterbox/water_equil.gro -o xtal_solv.pdb >& gmx_solvate.log
nwat=$( grep "Output configuration contains" gmx_solvate.log | awk '{print $7}' )
awk '/molecules/{stop=1} stop==0{print}' prot.top > md_model.top

#Counts existing ions per asymmetric unit and scales by ncopies for the full crystal
#Computes Na+ and Cl- counts to achieve ~0.1 M ionic strength (nwat*0.1/55) while also neutralizing the net charge
ions_pos=$( grep HETATM prot.pdb | grep Na+ | wc -l )
ions_neg=$( grep HETATM prot.pdb | grep Cl- | wc -l )
net_ion_charge=$( echo "($ions_pos - $ions_neg)*$ncopies" | bc )
echo "Positive ions, negative ions, net charge: "$ions_pos" "$ions_neg" "$net_ion_charge
if (( net_ion_charge >= 0 )); then
    ncl=$( echo "$nwat*0.1/55" | bc )
    nna=$( echo "$ncl + $net_ion_charge" | bc )
else
    nna=$( echo "$nwat*0.1/55" | bc )
    ncl=$( echo "$nna - $net_ion_charge" | bc )
fi

#gmx lines insert Cl- then Na+ ions by replacing solvent molecules at random positions
#nwat line recounts water after ion insertion (since some waters were displaced)
#cat through EOF lines writes the "molecules" section: one system1 1 entry per symmetry copy, then WAT/Cl-/Na+ counts
#cp line creates a copy of final solvate/ionized structure as md_model.pdb
gmx insert-molecules -f xtal_solv.pdb -ci Cl-.pdb -o xtal_solv_cl.pdb -replace SOL -nmol $ncl
gmx insert-molecules -f xtal_solv_cl.pdb -ci Na+.pdb -o xtal_solv_cl_na.pdb -replace SOL -nmol $nna
nwat=$( grep WAT xtal_solv_cl_na.pdb | wc | awk '{print $1/3}' )
cat >> md_model.top <<EOF
[ molecules ]
; Compound       #mols
${topall}
WAT              ${nwat}
Cl-              ${ncl}
Na+              ${nna}
EOF
cp xtal_solv_cl_na.pdb md_model.pdb
