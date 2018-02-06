# make directory for each run_lane
run_lane="$1"
mkdir $run_lane

# extract run and lane
run=`echo $run_lane | sed -e 's/_.*//'`
lane=`echo $run_lane | sed -e 's/.*_//'`

# get the cram files
imeta qu -z seq -d id_run = $run and lane = $lane and target = 1 and type = cram \
| grep : | awk '{ print $2 }' | paste - - -d/ \
| xargs -ixxx iget -K xxx $run_lane

chmod 664 $run_lane/*

# remove phiX control
find $run_lane | grep -E '#888\.' | xargs rm

# get and format the meta info.
for cram in $(find $run_lane | grep cram$ | sed -e 's/.*\///' | sed -e 's/\.cram$//'); do
    imeta ls -d /seq/$run/$cram.cram > $run_lane/$cram.imeta
    sn=$(grep -A 1 sample_supplier_name $run_lane/$cram.imeta | tail -1 | sed 's/ /_/g')
    sample_name=${sn:7}
    echo -e "$run_lane/$cram\t$sample_name" >> $run_lane/sampleInfo.txt
done
