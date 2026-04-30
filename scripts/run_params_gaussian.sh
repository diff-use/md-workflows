#How the script is called: bash ../run_params_gaussian.sh

#resn line sets the ligand name (GOL)
#export and source lines load gaussian env
#antechamber line converts GOL.pdb to gaussian input, GOL.gau, with net charge of "-2" and multiplicity of "1"
resn=GOL
export g16root=/Users/mewall/packages
source /Users/mewall/packages/g16/bsd/g16.profile
antechamber -fi pdb -fo gcrt -i ${resn}.pdb -o ${resn}.gau -nc -2 -m 1

#awk line inserts %NprocShared=8 into Gaussian input
#sed line replaces "#HF" with "#B3LYP/6-31+G(d,p)"
awk '{if ($0 ~ /Link/) {print $0 "\n%NprocShared=8"}else{print $0}}' ${resn}.gau > tmp
sed -i -e 's/\#HF/\#B3LYP\/6\-31+G\(d,p\)/g' tmp

#NOTE: First line code isn't necessary because 1 is always equal to 1!! 
#If the log file exists, compares current input to newly generated tmp
    #unchanged --> skip
    #changed --> overwrite and rerun 
if [ 1 == 1 ]; then
    if [ ! -f ${resn}.log ]; then
        echo "Running Gaussian optimization"
        mv tmp ${resn}.gau
        OMP_NUM_THREADS=8 g16 ${resn}.gau
    else
        if cmp -s "${resn}.gau" "tmp"; then
        echo "Gaussian input file is the same. Skipping run."
        else
        echo "Gaussian input file is different. Running Gaussian optimization."
        mv tmp ${resn}.gau
        OMP_NUM_THREADS=8 g16 ${resn}.gau
        fi	
    fi


sed -e 's/opt/pop\(chelpg\,regular\)/' ${resn}.gau > tmp
sed -i -e 's/molecule/grid/' tmp


if [ ! -f ${resn}_resp.log ]; then
    echo "Running Gaussian ESP calculation."
    mv tmp ${resn}_resp.gau
    OMP_NUM_THREADS=8 g16 ${resn}_resp.gau
else
    if cmp -s "${resn}_resp.gau" "tmp"; then
	echo "Gaussian input file is the same. Skipping run."
    else
	echo "Gaussian input file is different. Running Gaussian ESP calculation."
	mv tmp ${resn}_resp.gau
	OMP_NUM_THREADS=8 g16 ${resn}_resp.gau
    fi	
fi
fi


antechamber -fi gout -i ${resn}_resp.log -cf ${resn}_resp.crg -c resp -o ${resn}_gauss.ac -fo ac -rn ${resn}
#antechamber -fi gout -i ${resn}_resp.log -c resp -o ${resn}_gauss.ac -fo ac -rn ${resn}
antechamber -fi gout -i ${resn}_resp.log -o ${resn}_gauss.pdb -fo pdb -rn ${resn}


awk 'BEGIN{anum=0}{if (FNR==NR) {if ($1 == "ATOM" || $1 == "HETATM") {a[NR]=substr($0,32,23)}} else {if ($1 == "ATOM" || $1 == "HETATM") {anum=anum+1;print substr($0,1,31)""a[anum]""substr($0,55)} else {print $0}}}' ${resn}.pdb ${resn}_gauss.ac > ${resn}_resp.ac
grep ATOM ${resn}_resp.ac |awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{c=0;mc=0}{c=c+$9;if (abs(mc)<abs($9)) {mc = $9;n=NR} }END{if(c<0){ic = int(c-0.5)}else{ic = int(c+0.5)};printf "%9.6f %d %f\n",ic-c+mc,n,ic-c}' > charge_correction.dat
newcharge=$( awk '{print $1}' charge_correction.dat )
newcharge_anum=$( awk '{print $2}' charge_correction.dat )
awk -v newcharge="$newcharge" -v newcharge_anum="$newcharge_anum" 'BEGIN{anum=0}{if ($1 == "ATOM") {anum = anum + 1; if (anum == newcharge_anum) {printf "%s%9.6f%s\n", substr($0,0,55),newcharge,substr($0,65)}else{print $0}}else {print $0}}' ${resn}_resp.ac > temp.ac
mv temp.ac ${resn}_resp.ac


antechamber -fi ac -i ${resn}_resp.ac -fo mol2 -o ${resn}_resp.mol2 -rn ${resn}
atomtype -i ${resn}_resp.ac -o ${resn}_resp_gaff.ac -p gaff
prepgen -i ${resn}_resp_gaff.ac -o ${resn}_resp_gaff.prepc -f car
parmchk2 -i ${resn}_resp_gaff.prepc -o ${resn}_resp.frcmod -f prepc


cat > tleap_lig.in <<EOF
source leaprc.protein.ff19SB
source leaprc.gaff
loadamberparams ${resn}_resp.frcmod
loadamberprep ${resn}_resp_gaff.prepc
lig = loadmol2 ${resn}_resp.mol2
savepdb lig ${resn}_resp.pdb
saveamberparm lig ${resn}_resp.parm7 ${resn}_resp.rst7
quit
EOF
tleap -f tleap_lig.in
