## calc cohesive energy from dcd file
# $1 pdb filename without drude particle
# $2 dcd trajectory file to read
# $3 pdb filename with drude particle

cp ~/des_2018/cohesive/cohesive_openMM.py ./
cp ~/des_2018/cohesive/sapt_cohesive.xml ./

~/anaconda2/bin/python cohesive_openMM.py $1 $2 $3 | tee cohesive.log

