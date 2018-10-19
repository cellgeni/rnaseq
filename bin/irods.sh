#!/bin/bash

set -euo pipefail

sampleid=
studyid=
librarytype=
qc=true

while getopts :s:t:l:Qh opt
do
    case "$opt" in
    s)
      sampleid=$OPTARG
      ;;
    t)
      studyid=$OPTARG
      ;;
    l)
      librarytype="$OPTARG"
      ;;
    Q)
      qc=false
      ;;
    h)
      cat <<EOU
Usage: irods.sh -s SAMPLENAME [-t STUDYID] [-l LIBRARYTYPE] [Other-options]
Other-options:
  -Q    allow un-qc-ed data (drop 'manual_qc = 1' from imeta)
EOU
      exit
      ;;
    :) >&2 echo "Flag $opt needs argument"
        exit 1;;
    ?) >&2 echo "Flag $OPTARG unknown"
        exit 1;;
   esac
done

if [[ -z $sampleid ]]; then
  echo "I need a sample ID (-s option)"
  false
fi

part_studyid=
part_librarytype=
part_manualqc=

if $qc; then
  part_manualqc="manual_qc = 1 and"
fi
if [[ ! -z $studyid ]]; then
  part_studyid="study_id = $studyid and"
fi
if [[ ! -z $librarytype ]]; then
  ltype=$(printf "%q" "%$librarytype%")               # may contain spaces.
  part_librarytype="library_type like $ltype and"
fi

    # imeta command below:
    # (1) Eval is only slightly evil. We use it because the construction of
    # part_librarytype and part_studyid (and perhaps because librarytype can contain spaces).
    # Note that we escaped any spaces above using printf.
    #
    # (2) In imeta commands the most general (high-incident) query items are put first,
    # so do keep the ordering as it is now.
    # target = 1 "filters out files containing data such as spiked-in phiX and unconsented human",
    # This sounds very general, but testing says it's not, so we keep it rightmost in the imeta query.
    #
    # (3) In the presence of library_type clause, it seems faster if manual_qc is on the left.
    # in the absence of library_type clause, it seems faster if it is on the right.
    # We need to ask the NPG people for advice here.
    #

if [[ -z $librarytype ]]; then
  cmd="imeta qu -z seq -d $part_librarytype $part_studyid sample = $sampleid and $part_manualqc target = 1"
else
  cmd="imeta qu -z seq -d $part_librarytype $part_manualqc $part_studyid sample = $sampleid and target = 1"
fi
>&2 echo running "$cmd"

IMETA_OUTPUT=$(eval $cmd)

if [[ $(echo $IMETA_OUTPUT) == 'No rows found' ]]; then
    exit 64
else
    echo "${IMETA_OUTPUT}" \
    | perl -0777 -ne 'while (/collection:\s*(\S+)\ndataObj:\s*(\S+)/gs) { print "iget -K $1/$2\n" }' \
    | bash -e
fi
 

