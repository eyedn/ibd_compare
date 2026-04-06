#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           create_nonPyrho_sample_lists.sh
###############################################################################

set -euo pipefail

# input
# override by prefixing: SAMPLE_SIZE=5000 ./script.sh
: "${SAMPLE_SIZE:=500}" 
POPS=(YRI CEU CHB)

# output files
BASE="${SCRATCH}/stdpopsim_runs/pop_info"
mkdir -p "${BASE}"

POP_LISTS_FILE="${BASE}/pop_sample_paths.txt"
: > "${POP_LISTS_FILE}"

# sanity check
if ! [[ "${SAMPLE_SIZE}" =~ ^[0-9]+$ ]] || [[ "${SAMPLE_SIZE}" -le 0 ]]; then
    echo "ERROR: SAMPLE_SIZE must be a positive integer; got '${SAMPLE_SIZE}'" >&2
    exit 1
fi

# create each pop file with contiguous index blocks:
# pop k gets [k*SAMPLE_SIZE, (k+1)*SAMPLE_SIZE - 1]
for idx in "${!POPS[@]}"; do
    pop="${POPS[$idx]}"
    out_file="${BASE}/${pop}_samples.txt"
    : > "${out_file}"

    start=$(( idx * SAMPLE_SIZE ))
    end=$(( (idx + 1) * SAMPLE_SIZE - 1 ))

    # write sample IDs
    for i in $(seq "${start}" "${end}"); do
    echo "tsk_${i}" >> "${out_file}"
    done

    # append path to master list
    echo "${out_file}" >> "${POP_LISTS_FILE}"

    echo "Created ${out_file} with tsk_${start}..tsk_${end}"
done

echo "Created: ${POP_LISTS_FILE}"
