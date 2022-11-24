inPath=/eos/jeodpp/data/projects/CLIMEX/mentaAltimetryHsValidation/data/rawData/
outPath=/eos/jeodpp/data/projects/CLIMEX/mentaAltimetryHsValidation/data/rawDataSpecific/

for year in {2003..2009}
do
	for month in 01 02 03 04 05 06 07 08 09 10 11 12
	do
		mkdir -p $year/$month
		file2transfer=$inPath/$year/$month/dt_global_en_phy_l3_$year$month*nc
		echo $file2transfer
		ln -sf $file2transfer $outPath/$year/$month/.
		#exit	
	done
done
