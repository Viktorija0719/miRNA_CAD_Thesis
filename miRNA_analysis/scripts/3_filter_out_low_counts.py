#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DESCRIPTION:
Filters out miRNAs with low expression based on median counts across groups.
Only miRNAs with a median ≥ 10 in at least one group are retained.

Also identifies miRNAs with low variability (fewer than 20 unique values).

INPUT:
- A tab-delimited count matrix file where the first column is GeneID and 
  remaining columns are samples with names starting with '0_CAD_', '1_CAD_', or '2_CAD_'.

OUTPUT:
- A file with retained miRNAs and their group-wise medians
- A file with retained miRNAs in original count matrix format
USAGE:
python filter_out_low_counts.py countData.txt kept_medians.txt kept_data.txt
"""

import os
import sys
import statistics

def main():
    if len(sys.argv) != 4:
        print("USAGE: python filter_out_low_counts.py INPUTFILENAME OUTPUT_MEDIANS OUTPUT_MATRIX")
        sys.exit(1)

    rawCountFileName = sys.argv[1]
    mediansOutputPath = sys.argv[2]
    matrixOutputPath = sys.argv[3]
    low_variability_threshold = 1  

    group_prefixes = ['0_CAD_', '1_CAD_', '2_CAD_']

    try:
        with open(rawCountFileName) as rawCounts, \
             open(mediansOutputPath, 'w') as keptWithMediansFile, \
             open(matrixOutputPath, 'w') as keptDataFile:

            header = rawCounts.readline().strip().split("\t")
            group_indices = {
                prefix: [i for i, col in enumerate(header) if col.startswith(prefix)]
                for prefix in group_prefixes
            }

            # Write headers
            medianHeader = "GeneID\tMedianGroup0\tMedianGroup1\tMedianGroup2\n"
            keptWithMediansFile.write(medianHeader)
            keptDataFile.write("\t".join(header) + '\n')

            for line in rawCounts:
                fields = line.strip().split("\t")
                geneId = fields[0]
                values = [float(x) for x in fields[1:]]

                # Skip miRNAs with low variability
                if len(set(values)) < low_variability_threshold:
                    continue

                
                # List of miRNAs to exclude
                exclude_miRNAs = {
                    "hsa-miR-1299", "hsa-miR-20a-3p", "hsa-miR-3157-5p",
                    "hsa-miR-34a-5p", "hsa-miR-4516", "hsa-miR-548ax",
                    "hsa-miR-610", "hsa-miR-6818-5p"
                }

                # Skip if the miRNA is in the exclusion list
                if geneId in exclude_miRNAs:
                    continue



                group_medians = []
                for prefix, indices in group_indices.items():
                    group_counts = [float(fields[i]) for i in indices]
                    group_median = statistics.median(group_counts)
                    group_medians.append(group_median)

                if any(median >= 10 for median in group_medians):
                    medianLine = f'{geneId}\t' + '\t'.join(map(str, group_medians)) + '\n'
                    keptWithMediansFile.write(medianLine)
                    keptDataFile.write(line)

    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(2)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass
