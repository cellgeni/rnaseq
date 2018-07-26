#!/bin/bash

set -euo pipefail

IMETA_OUTPUT="$(imeta qu -z seq -d study_id = $1 sample = $2 and target = 1 and manual_qc = 1)"
if [ "$(echo $IMETA_OUTPUT)" != 'No rows found' ]; then
    echo "${IMETA_OUTPUT}" \
    | perl -0777 -ne 'while (/collection:\s*(\S+)\ndataObj:\s*(\S+)/gs) { print "iget -K $1/$2\n" }' \
    | bash -e
fi
