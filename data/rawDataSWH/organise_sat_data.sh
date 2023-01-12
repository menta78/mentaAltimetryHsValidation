sat="j1"
inPath=/eos/jeodpp/data/projects/CLIMEX/swhCMEMS/cci_obs-wave_glo_phy-swh_my_$sat-l3_PT1S/
outPath=/eos/jeodpp/data/projects/CLIMEX/mentaAltimetryHsValidation/data/rawDataSWH/

for year in {2002..2012}
do
	for month in 01 02 03 04 05 06 07 08 09 10 11 12
	do
		file2transfer=$inPath/global_vavh_l3_rep_${sat}_$year$month*nc
		ls $file2transfer && ln -sf $file2transfer $outPath/$year/$month/.
	done
done
