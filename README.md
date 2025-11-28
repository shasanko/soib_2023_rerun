This repository describes how to analyse the 2023 full-country dataset using the current SoIB codebase (that is, after the updates made by Shree Kumar).

## Required input files (the ones used to create the 2023 report)

```
00_data/
├── SoIB_mapping_2022.csv
└── spec_misid.RData

01_analyses_full/
├── dataforanalyses.RData
├── fullspecieslist.csv
├── randomgroupids.RData
└── specieslists.RData
```

## Previous codebase

1. `00_scripts/filter_data_for_species.R` processes eBird data and saves `sub_samp_locs.csv`, which contains LocalityID and groupID.

2. `00_scripts/create_random_groupids.R` processes `sub_samp_locs.csv` and saves `randomgroupids.RData`, a matrix of 1000 columns where each column contains a set of groupIDs.

3. `00_scripts/create_random_datafiles.R` processes `dataforanalyses.RData` using each column of `randomgroupids.RData` and saves 1000 filtered data files, one per random groupID column.

## Latest codebase

1. `p1s3.R` processes eBird data and saves `sub_samp_locs.csv`.

2. `p2s1a.R` processes `sub_samp_locs.csv` and saves 1000 individual randomgroupID files instead of one large matrix as in the previous codebase. Each file corresponds to one column of the matrix. GroupIDs are also converted to numeric to reduce memory usage.

3. `p2s2a.R` processes `dataforanalyses.RData` and saves a version with numeric groupIDs, enabling later filtering using the individual randomgroupIDs in `p3s1a.R`.

## Analysis of 2023 data using the new codebase

1. Run `p2s2a.R`. It processes `dataforanalyses.RData` and saves two output files in `01_analyses_full/`:
   
   - `dataforanalyses.RData-data_opt` (restructured data with numeric groupIDs)  
   
   - `dataforanalyses.RData-metadata` (metadata extracted from the original file)

2. Run `remap-rgrids.R`. It processes `randomgroupids.RData` and saves 1000 random groupID files in `01_analyses_full/`.

3. A file called `00_data/current_soib_migyears.RData` is required by the new codebase. It is derived from the raw eBird data and specifies the migratory years covered by the analysis. Because the raw eBird data used in the 2023 analysis is inaccessible, `prep_migyears_2023.R` recreates the file using the migration years listed in the 2023 report.

4. Run `config/localhost/config.R` with the following specifications:

   - `species_to_process = c("Gray Francolin", "Ashy Prinia", "Black-winged Kite")`
   
   - `my_assignment = 1:1000`
   
   The file is used to specify any required configuration.

5. Run `p3s1a.R`. This is the model fitting step and saves as many trend files as the number of assignments in `01_analyses_full/trends`.

6. Run `00_scripts/resolve_trends_only.R` to combine trend files from all assignments and resolve them. This step saves the combined file in `01_analyses_full/trends`.