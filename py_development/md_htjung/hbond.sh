# $1: chol/any

mole=$1
if [ $mole = 'chol' ]
then
	python ~/des/md_hbond.py -i md_nvt_pos.dcd -s start_drudes.pdb -sel ~/des/select/def_OH_Otf.hbond -dih YES -o oh_otf
fi

python ~/des/md_hbond.py -i md_nvt_pos.dcd -s start_drudes.pdb -sel ~/des/select/def_CN_Otf.hbond -dih YES -o nch_otf
python ~/des/md_hbond.py -i md_nvt_pos.dcd -s start_drudes.pdb -sel ~/des/select/def_CN_Ntf.hbond -dih YES -o nch_ntf

