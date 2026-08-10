#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           require_sim_id_contract.sh
###############################################################################

set -euo pipefail

contract_file="$1"
expected_contract="$2"

[[ -s "$contract_file" ]] || {
    echo "ERROR: missing simulation ID-contract marker: $contract_file" >&2
    exit 1
}

actual_contract="$(<"$contract_file")"
[[ "$actual_contract" == "$expected_contract" ]] || {
    echo "ERROR: simulation ID contract mismatch; expected " \
        "${expected_contract}, found ${actual_contract}" >&2
    exit 1
}
