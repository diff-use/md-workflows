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
gmx solvate -cp xtal.pdb -cs waterbox/water_equil.gro -o xtal_solv.pdb >& gmx_solvate.log
nwat=$( grep "Output configuration contains" gmx_solvate.log | awk '{print $7}' )
awk '/molecules/{stop=1} stop==0{print}' prot.top > md_model.top
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
