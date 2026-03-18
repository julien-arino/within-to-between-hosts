# Code directory

This directory contains the code used in the paper *Using virtual patients to parametrise an age-of-infection model* by Julien Arino, Morgan Craig, Clotilde Djuikem,Kang-Ling Liao and Stéphanie Portet. The code mixes `julia` and `R`.

## To run the code

Several bash files are provided to run the code:

- `install-required-packages`
- `run-sensitivity-analysis`
- `run-cohort`
- `run-between-host`

Note however that the computations can be quite long, so details are provided below.

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
Rscript process-PRCC-withV0.R

# Plot the PRCC figures
Rscript plot-Figure-03-B1-PRCC-within-host-withV0.R
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

