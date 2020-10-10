# $1 : ch or tm

prefix=$1

## tail
if [ $prefix == 'tm' ]; then
	python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/CH3 -s2 ~/des/select/Ntf -o rdf_CH3_Ntf
	python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/CH3 -s2 ~/des/select/Otf -o rdf_CH3_Otf
	python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/CH3 -s2 ~/des/select/Ftf -o rdf_CH3_Ftf
	python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/CH3 -s2 ~/des/select/CH3 -o rdf_CH3_CH3
else
	python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/OH -s2 ~/des/select/Ntf -o rdf_OH_Ntf
        python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/OH -s2 ~/des/select/Otf -o rdf_OH_Otf
        python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/OH -s2 ~/des/select/Ftf -o rdf_OH_Ftf
        python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/OH -s2 ~/des/select/OH -o rdf_OH_OH
fi


## head
#python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/CN -s2 ~/des/select/Ntf -o rdf_CN_Ntf
#python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/CN -s2 ~/des/select/Otf -o rdf_CN_Otf
#python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/CN -s2 ~/des/select/Ftf -o rdf_CN_Ftf
#python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/CN -s2 ~/des/select/CN -o rdf_CN_CN

#if [ $prefix == 'ch' ]; then
#	python ~/des/md_rdf.py -i md_nvt_pos.dcd -s start_drudes.pdb -step 10 -rmax 2 -nbin 100 -inc except -s1 ~/des/select/CN -s2 ~/des/select/OH -o rdf_CN_OH
#fi

foldn=$prefix"-2nvt"
cp rdf_* ~/des/analysis/rdf/${foldn}

