This repo contains the scripts used for the paper *Species selection regime and phylogenetic tree shape*.

### To set up:

1) Clone this repo. Note the `scripts` directory contains all the scripts needed to run the simulations under the range of models explored in the paper.

2) Prepare the following directory structure. Note that the paths within each script may need to be adjusted to match your directory structure.
   - `scripts` – The `scripts` directory from this repo.
   - `outputs` – This is the subdirectory to which the simulation results are written (as .RData files) by the simulation scripts.
   - `results` – This is the directory to which the data summaries are written (as .csv files) by the data summarization scripts.
   - `final_results` – This is the directory to which the final data summaries, with parameter settings, are written (as .csv files) by the parameter summarization scripts.

### To run a simulation:

1)  Execute `setup.R`. This loads the core functions. Note that two lines in the "Helper function to simulate quasse trees with regimes" control whether or not the simulations run under a tree size constraint.

2)  Now run a simulation under a particular model. Each model has a separate simulation script that provides the parameter settings and executes the simulations. The results of each simulation are written to a list which is saved as a .RData file. The provided simulation scripts are:
    - `uniform.R`
    - `static_Gaussian_speciation.R`
    - `static_skewed_speciation.R`
    - `static_Gaussian_extinction.R`
    - `shifting_Gaussian_speciation.R`
    - `shifting_neg_skewed_speciation.R`
    - `shifting_pos_skewed_speciation.R`
    - `shifting_Gaussian_extinction.R`

### To summarize the results of a simulation:

1)  Execute `summarize.R` which summarizes the results of the simulations run under a particular model/parameter combination. The results are written to a .csv file. Note that, as written, these scripts summarize more tree attributes than those used in the paper.

2)  Execute `summarize_parameters.R` which furnishes the result summary .csv files with the appropriate parameter settings underlying each simulation. Note that you will need to select the appropriate simulation type at the top of the file and in some places ensure that the correct parameter set is chosen.

