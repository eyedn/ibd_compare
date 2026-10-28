#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           software_profiles.sh
###############################################################################

# update software parameters based on input; see cases clause for profiles
select_software_profile() {
    local software="$1"
    shift
    if [[ $# != 1 ]]; then
        echo "ERROR: expected exactly one ${software} profile argument" >&2
        return 1
    fi
    case "${software}:$1" in
        hapibd:seed2_gap1000) HAP_MIN_SEED=2.0; HAP_MAX_GAP=1000 ;;
        hapibd:seed1_gap1000) HAP_MIN_SEED=1.0; HAP_MAX_GAP=1000 ;;
        hapibd:seed2_gap2000) HAP_MIN_SEED=2.0; HAP_MAX_GAP=2000 ;;
        rapid:w100) RAP_WINDOW_SIZE=100 ;;
        rapid:w200) RAP_WINDOW_SIZE=200 ;;
        rapid:w300) RAP_WINDOW_SIZE=300 ;;
        *)
            echo "ERROR: unknown ${software} profile: $1" >&2
            return 1 ;;
    esac
    echo "profile=$1; seed=${HAP_MIN_SEED}; gap=${HAP_MAX_GAP}; window=${RAP_WINDOW_SIZE}" >&2
}
