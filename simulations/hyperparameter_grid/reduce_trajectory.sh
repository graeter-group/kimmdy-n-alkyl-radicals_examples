#!/bin/bash

# Check if a parameter is provided
if [ $# -lt 2 ]; then
  echo "Usage: $0 <input_file1> <input_file2>"
  exit 1
fi

# Initial input and output filenames
input_trajectory="$1"
input_topology="$2"
wrk_dir="$3"
output_trajectory="$1"

# Iterate 6 times to extract every 10th frame
for i in {1..6}; do
    # Construct output filenames
    new_output="${output_trajectory%.xtc}_${i}.xtc"

    # Use trjconv to extract every 10th frame
    printf "0\n"$ | gmx trjconv -f "${wrk_dir}/${output_trajectory}" -s "${wrk_dir}/${input_topology}" -skip 10 -o "${wrk_dir}/${new_output}"

    # Update output_trajectory for the next iteration
    output_trajectory="${new_output}"
done