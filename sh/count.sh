#!/bin/bash

declare -A count_map

# Read input and count occurrences
while read -r item; do
    ((count_map["$item"]++))
done

# Print the counts
for item in "${!count_map[@]}"; do
    printf "%s: %d\n" "$item" "${count_map[$item]}"
done
