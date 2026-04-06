#!/usr/bin/env bash

NUM_REPS=10

TMAPS=(
    "HapMapII_GRCh38"
    "DeCodeSexAveraged_GRCh38"
    "DefaultUniform"
    "Pyrho_GRCh38"
)

IMAPS=(
    "HapMapII_GRCh38"
    "DeCodeSexAveraged_GRCh38"
    "Flat1Mb1cM"
    "Pyrho_GRCh38"
)

POPS=(YRI CEU CHB)

NUM_TMAPS=${#TMAPS[@]}
NUM_IMAPS=${#IMAPS[@]}
NUM_POPS=${#POPS[@]}

decode_task () {
    local task_id="$1"

    local rep=$(( task_id % NUM_REPS + 1 ))
    local tmp=$(( task_id / NUM_REPS ))

    local tmap_idx=$(( tmp % NUM_TMAPS ))
    tmp=$(( tmp / NUM_TMAPS ))

    local imap_idx=$(( tmp % NUM_IMAPS ))
    local pop_idx=$(( tmp / NUM_IMAPS ))

    local pop=${POPS[$pop_idx]}
    local tmap=${TMAPS[$tmap_idx]}
    local imap=${IMAPS[$imap_idx]}

    echo "$task_id $pop $tmap $imap $rep"
}

FAILED_IDS=(355 356 357 394 395 396 397 434 435 436 437 439 475 476 477)

for id in "${FAILED_IDS[@]}"; do
    read task pop tmap imap rep <<<"$(decode_task "$id")"
    echo "=== replaying $id: pop=$pop tmap=$tmap imap=$imap rep=$rep ==="

    IBD_DIR="/home1/karatas/scratch/ibd_compare/rapid_sim_tmap${tmap}_imap${imap}_pop${pop}_d2_w250_r10_s2"
    IBD_OUT="${IBD_DIR}/${pop}_chr1_rep${rep}_out"
    IBD_BED="${IBD_DIR}/${pop}_chr1_rep${rep}_ibd.bed.gz"

    # clean previous attempt for this replicate
    if [[ -d "$IBD_OUT" ]]; then
        echo "Removing existing output dir: $IBD_OUT"
        rm -rf "$IBD_OUT"
    fi

    if [[ -f "$IBD_BED" ]]; then
        echo "Removing existing bed file: $IBD_BED"
        rm -f "$IBD_BED"
    fi

    bash /home1/karatas/ibd_compare/job_scripts/rapid_inference.sh \
        "$pop" "$tmap" "$imap" "$rep" \
        2 250 10 2 \
        "$IBD_DIR"
done
