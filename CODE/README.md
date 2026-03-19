# Code directory

This directory contains the code used in the paper *Using virtual patients to parametrise an age-of-infection model* by Julien Arino, Morgan Craig, Clotilde Djuikem,Kang-Ling Liao and Stéphanie Portet. The code mixes `julia` and `R`.

## To run the code

Several bash files are provided to run the computational and plotting code:

- `install-required-packages`
- `run-sensitivity-analysis`
- `run-cohort`
- `run-between-host`

Details of the code in these files is provided below.

Note that computations are memory intensive and greatly benefit from large scale parallelisation. Running a small cohort (10,000 individuals) is easily carried out on a middle range laptop, but to run the full 1 million individuals cohort simulation, you will need a serious amount of RAM (256 GB recommended) and will greatly benefit from being to call on at least 30 threads. 

We chose to save the results in a format that can be read by R (qs2) to facilitate the subsequent analysis, but this is not the most memory efficient format, since the simulation results need to be copied to R to the be saved. A potentially more memory efficient format is the Arrow format, which exists in both Julia and R. However, Arrow saves flat files, so R will then need to convert to a list of dataframes, which can be memory intensive and very long.

### 1. Install required packages

```bash
# Install Julia packages
julia install-required-packages.jl

# Install R packages
Rscript install-required-packages.R
```

### 2. Run the sensitivity analysis

```bash
# Run the sensitivity analysis
julia run-sensitivity-analysis-sims-withV0.jl

# Run sensitivity analysis PRCC
Rscript process-PRCC.R

# Plot the PRCC figures
Rscript plot-Figure-03-B1-PRCC-within-host.R
```

### 3. Run the cohort simulation

```bash
# Run the cohort simulation
julia run-virtual-cohort-sims.jl

# Process the cohort simulation
Rscript process-cohort-classify-individuals-fct-xih-xid.R

# Process the cohort simulation
Rscript process-cohort-assign-zero-transmssion.R
```

### 4. Run the between-host simulation

