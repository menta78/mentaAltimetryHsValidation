
for year in {2002..2020} 
do
	for month in 01 02 03 04 05 06 07 08 09 10 11 12
	do
		for file in $year/$month/*
		do
			ncdump -h $file || rm $file
		done
	done
done
