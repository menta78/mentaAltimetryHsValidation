while [ $# -gt 0 ] ; do
  case $1 in
	-i | --in) inFolder=${2};;
	-o | --out) outFolder=${2};;
#	-r | --root) rootFile=${2};;
  esac
  shift
done

if [ -z "${inFolder+xxx}" ]
then
	echo "Error! Undefined input folder."
	exit
fi

if [ -z "${outFolder+xxx}" ]
then
	echo "Error! Undefined output folder."
	exit
fi

#if [ -z "${rootFile+xxx}" ]
#then
#	echo "Error! Undefined root file."
#	exit
#fi

# To create the output folder structure
# for year in {1993..2022}; do for month in 01 02 03 04 05 06 07 08 09 10 11 12; do mkdir -p $year/$month; done; done 

for year in {1993..2022}
do
	for month in 01 02 03 04 05 06 07 08 09 10 11 12
	do
		res=$(printf '%s\n' $inFolder/dt_global_*_phy_l3_${year}${month}*)
		echo $res
		rsync -va $res $outFolder/$year/$month/.
		
		#ls ${inFolder}/${rootFile}_${year}${month}*
	done
done

