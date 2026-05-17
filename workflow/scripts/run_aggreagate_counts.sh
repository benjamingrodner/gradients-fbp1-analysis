#!/bin/bash

GLOB=$1
OUT=$2
THREADS=$3

FNS=$(ls "$GLOB")

python3 workflow/scripts/aggregate_counts.py \
    --jobs "$THREADS" \
    $FNS \
    "$OUT"

