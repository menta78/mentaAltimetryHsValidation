inPath=/eos/jeodpp/data/projects/CLIMEX/ssh-global-l3/
outPath=/eos/jeodpp/data/projects/CLIMEX/mentaAltimetryHsValidation/data/rawData/

for year in {2002..2020}
do
	for month in 01 02 03 04 05 06 07 08 09 10 11 12
	do
		echo dt_global_*phy_l3_$year$month*nc
		file2transfer=$inPath/dt_global_*phy_l3_$year$month*nc
		echo $file2transfer
		ln -sf $file2transfer $outPath/$year/$month/.
		exit
	done
done
