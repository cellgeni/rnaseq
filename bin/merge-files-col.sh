#!/bin/bash

# Merge single columns from multiple files using standard command line
# utilities and additionally the custom 'transpose' program from the
# Github micans/reaper respository.
#
# Several checks (e.g. row name identiy) and options are available, see the -h
# output.
# Intermediate output and final output are both gzipped.
# A check option by double-transposing is available (-T option).
#
# This program uses cut(1) to isolate the columns; it writes each sequentially
# using paste(1) as a single line and appends the lines.  This leads to a
# result which is the transposed of what is needed; In terms of a standard
# application (input columns are gene counts across multiple conditions) output
# row names now correspond to samples and column names correspond to gene names.
#
# The final result is obtained by using the custom program 'transpose', which
# is a C implementation of the needed transpose operation. This combined
# approach is generally a lot faster than using R to load tables and cbind to
# combine the columns into a new table, or even using a standard non-optimised
# python solution.
#
# TODO: clean-up of extra files; perhaps work in /tmp.
#


set -euo pipefail

                        # For these options, see the documentation under -h
                        # for a description of what they do.
basename=out
test_tp=false
header=true
rowdatastart=2
rownames_check=true
persist=false
help=false
col=2
colname=
colname_expected=
stripright=
stripleft=
strippath=false
output=
skip=0
verbose=false
metafile=""

function clean_up {
   thestatus=$?
   if $help; then
      true
   elif [[ $thestatus == 0 ]]; then
      >&2 echo "[happy] Script succeeded - gzipped output in $output"
      if [[ -e "$tpfile.TP" ]]; then
         rm "$tpfile.TP"
         >&2 echo "removed file $tpfile.TP after successful test"
      fi
   else
      >&2 echo "[grumpy] Script failed"
   fi
}


   ##  register the function clean_up to be run when script exits,
   ##  after all the other code in this script.
   ##  At the moment it does not do much.
   ##
trap clean_up INT TERM EXIT


while getopts :c:i:C:E:b:o:s:x:y:VLTHPZh opt
do
    case "$opt" in
    b)
      basename=$OPTARG
      ;;
    c)
      col=$OPTARG
      ;;
    x)
      stripleft=$OPTARG
      ;;
    y)
      stripright=$OPTARG
      ;;
    C)
      colname=$OPTARG
      ;;
    Z)
      persist=true
      ;;
    H)
      header=false
      rowdatastart=1
      ;;
    E)
      colname_expected=$OPTARG
      ;;
    P)
      strippath=true
      ;;
    i)
      metafile=$OPTARG
      ;;
    o)
      output=$OPTARG
      ;;
    T)
      test_tp=true
      ;;
    L)
      rownames_check=false
      ;;
    V)
      verbose=true
      ;;
    s)
      skip=$(($OPTARG+1))
      ;;
    h)
      cat <<EOU
Usage: merge-col-files.sh [ -c <NUM> ] [ MORE OPTIONS ] FILES
Options:
-b basename    basename for output (default $basename)
-c <NUM>       column index (default 2)
-i filename    specify file that contains the file paths to read
-C <colname>   rename column to <colname> (default taken from column head)
-E <check>     require file column name to be <check> (mnemonic: expect)
                  Caveat: currently this check is only performed for the first file in the list.
-o <name>      write output to name (default: <basename>.<colname>.txt.gz)
-x pattern     for result column name, strip pattern from the left of file name (greedy match!)
-y pattern     for result column name, strip pattern from the right of file name
-P             remove leading path (this is performed before -x stripping)
-T             test use of the transpose program (re-transpose, compare with original)
-L             lax; do not test whether row names (first column) are identical between files
-V             verbosity on
-H             no header line (default assumes first line is header line)
-Z             keep going if some files fail checks
-s <num>       skip first <num> lines
EOU
      help=true
      exit
      ;;
    :) >&2 echo "Flag $opt needs argument"
        exit 1;;
    ?) >&2 echo "Flag $OPTARG unknown"
        exit 1;;
   esac
done


# Validate how inputs are specified, and set firstfile.
# Accept either trailing arguments, or metafile, but not both.
# If trailing arguments, set metafile subsesquently to /dev/null
# This is slightly hackish.
# firstfile=
# ------------------------------------------------
OPTIND=$(($OPTIND-1))
shift $OPTIND

if [[ -n $metafile ]] && (( $# > 0 )); then
  echo Use of trailing filenames with -i is currently not supported
  false
elif [[ -z $metafile ]] && (( $# == 0 )); then
   >&2 echo "Error (no files provided either on command line or with -i option)"
   false
fi
if (( $# > 0 )); then
  firstfile=$1
  metafile=/dev/null
else
  firstfile=$(head -n 1 $metafile)
fi
# ------------------------------------------------


# Manage column labels.
# Note: when skip is read from command line, it is incremented by one.
# ------------------------------------------------
if $header; then
     # get the column labels of interest, e.g. gene counts.
     #
  read -a colnames <<EOF
  $(tail -n +$skip "$firstfile" | head -n 1 | cut -f $col)
EOF

  >&2 echo "Found column name ${colnames[0]}"
  if [[ -n "$colname" ]]; then
     >&2 echo "Using supplied column name $colname"
  else
     colname=${colnames[0]}
  fi

  if [[ -n "$colname_expected" && "$colname_expected" != "$colname" ]]; then
     >&2 echo "Expected column name $colname_expected but found $colname"
     false
  fi
else
  colname=col${col}
fi
# ------------------------------------------------


tpfile=$basename.$colname.TP.gz
dstfile=$basename.$colname.txt.gz
rownamefile=$basename.$colname.rn
idx=1

   # Get the row labels of interest, e.g. gene names.
   # Save them for comparison with rownames in later files,
   # and write them to the destination file (gzipped).
   #
(tail -n +$skip "$firstfile" | cut -f 1 | tee "$rownamefile" | paste -s | gzip) > "$tpfile"
nrows=$(cat $rownamefile | wc -l)
   #
   # Note: below we append to tpfile (after gzipping). This is OK because
   # zipfiles can be concatenated.

n=0
(
for f in $@ $(cat $metafile); do
   level=$f
   if $strippath; then
      level=${level##*/}
   fi
   if [[ -n "$stripleft" ]]; then
      level=${level##$stripleft}
   fi
   if [[ -n "$stripright" ]]; then
      level=${level%$stripright}
   fi
   if $verbose; then
    >&2 printf "%5d $level\n" $idx
   fi
   idx=$(($idx+1))

   if $rownames_check && ! cmp -s "$rownamefile" <(tail -n +$skip $f | cut -f 1) ; then
      >&2 echo "rownames check failed when comparing $firstfile and $f"
      if $persist; then
        continue
      else
        false
      fi
   fi

            # Note we cannot simply add $((skip+rowdatastart)) due to UNIX tail 1-based indexing.
   (echo "$level"; cut -f $col $f | tail -n +$skip | tail -n +$rowdatastart) | transpose --quiet -i - --nozip -o -
   n=$((n+1))
   >&2 echo -n '.'
   if (( $n % 100 == 0 )); then >&2 printf " %5d\n" $n; fi
done
)  | gzip >> $tpfile


if [[ -z "$output" ]]; then
   output=$dstfile
fi

transpose -i $tpfile -o $output

if $test_tp; then
   transpose -i "$output" -o "$tpfile.TP2.gz"
   files="$tpfile and $tpfile.TP2.gz"
   start_byte=1
   if ! $header; then       # A leading tab will have been added, so skip it.
    start_byte=2
   fi
   if ! cmp -s <(zcat "$tpfile") <(zcat "$tpfile.TP2.gz" | tail -c +$start_byte); then
      >&2 echo "Files $files are not byte-identical (output in $output)"
      false
   else
      >&2 echo "Test passed - $files are byte-identical (output in $output)"
   fi
fi   



# Todo:
# - improve tpfile dstfile rownamefile; enable directories and play
#   nicely with -o option. e.g. -o output/base
# - -s option exists, perhaps also skip initial regex (e.g. initial lines
#   starting with #
