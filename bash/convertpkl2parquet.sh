#!/bin/bash

# Directory to search in — default is current directory if not provided
BASE_DIR="//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/"

# Temporary Python script used for conversion
PY_SCRIPT=$(mktemp)

# Create the Python conversion script
cat << 'EOF' > "$PY_SCRIPT"
import pandas as pd
import sys
import os

pkl_file = sys.argv[1]
parquet_file = sys.argv[2]

try:
    df = pd.read_pickle(pkl_file)
    df.to_parquet(parquet_file)
    print(f"Converted: {pkl_file} -> {parquet_file}")
except Exception as e:
    print(f"Failed to convert {pkl_file}: {e}", file=sys.stderr)
EOF

# Find all .pkl files and convert them
find "$BASE_DIR" -type f -name "*centsAB.pkl" | while read -r pkl_file; do
    parquet_file="${pkl_file%.pkl}.parquet"
    python "$PY_SCRIPT" "$pkl_file" "$parquet_file"
done

find "$BASE_DIR" -type f -name "*mappingsAB.pkl" | while read -r pkl_file; do
    parquet_file="${pkl_file%.pkl}.parquet"
    python "$PY_SCRIPT" "$pkl_file" "$parquet_file"
done

# Cleanup
rm "$PY_SCRIPT"