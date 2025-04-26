#!/bin/bash

# Loop through all .m files in the current directory
for file in *.m; do

    # Extract filename without extension
    filename=$(basename "$file" .m)
    
    # Submit job with filename as job name and in the MATLAB command
    module load matlab

    sbatch -J "$filename" -t "15-00" --cpus-per-task=4 --mem-per-cpu=32G --wrap="matlab -r $filename"
    
    echo "Submitted job for $filename"
done
