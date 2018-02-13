# make directory for each run_lane
run_lane="$1"

# extract run and lane
run=`echo $run_lane | sed -e 's/_.*//'`
lane=`echo $run_lane | sed -e 's/.*_//'`

# get the cram files
imeta qu -z seq -d id_run = $run and lane = $lane and target = 1 and type = cram \
| grep : | awk '{ print $2 }' | paste - - -d/ \
| xargs -ixxx iget -K xxx ./

chmod 664 *

# remove phiX control
find ./ | grep -E '#888\.' | xargs rm

# get and format the meta info.
for cram in $(find ./ | grep cram$ | sed -e 's/.*\///' | sed -e 's/\.cram$//'); do
    imeta ls -d /seq/$run/$cram.cram > $cram.imeta
    sn=$(grep -A 1 sample_supplier_name $cram.imeta | tail -1 | sed 's/ /_/g')
    sample_name=${sn:7}
    echo -e "$run_lane/$cram\t$sample_name" >> $run_$lane_sampleInfo.txt
done