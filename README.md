# Metapopulation_MRSA
Code Repository for "Population Structure Drives Differential Methicillin-resistant Staphylococcus aureus Colonization Dynamics"
----

The Python scripts in this repository require the 'os', 'stochpy', 'numpy', 'stats', 'argparse', 'matplotlib.pyplot, and 'random' libraries, all of which are available via pip and scipy. These scripts are single run programs performed on an array. User's slurm admission scripts (.srun) for simulation on a cluster are not included.    

The R scripts require the 'readr', 'sm', and 'vioplot' libraries.

Contents:
----
* 'SST.psc': Configuration file for the Single Staff Type (SST) model. This model assumes random mixing of all ICU staff (nurses and MD) evenly with all patients. 

* 'Nurse_MD.psc': Configuration file for a model with the Nurses and MD separated as different staff types. All nurses and the MD are still assumed to randomly mix with all patients in the ICU.

* 'MetaPop.psc': Configuration file for a model with a specific population structure of nurses and patients that represent a metapopulation structure. The ratio of nurses to patients is 3:1 with a total of 6 strict nurse-patient groups where the nurse does not interact with any patient outside of the origianlly assigned group or other staff. The MD remains to interact randomly with all patients. 

* 'Gamma.psc': Configuration file for an expanded version of the metapopulation model (MetaPop.psc) that includes a new parameter, gamma, that represents the proportion of time a nurse spends in the originally assigned group of patients. 

* 'model_comparison.py': This script runs each of the models and generates the cumulative number of MRSA acquisitions for comparison.

* 'model_analysis.R': This R script produces the violoin plot figure (Figure 4) from the manuscript.

* 'psi_fit_SST.py': One of three scripts that performs a model recalibration of the parameter psi for the SST model. This allows for a comparison of the models in a setting where their outcomes are equal and examines how each model might influence the value of an estimated parameter.

* 'psi_fit_NurseMD.py': The second of three scripts that performs a model recalibration of the parameter psi for the Nurse-MD model. 

* 'psi_fit_MetaPop.py': The third of three scripts that performs a model recalibration of the parameter psi for the metapopulation model. 

* 'psi_fit_analysis.R': This script performs an analysis for comparison between models, where the median of each parameter estimate distribution was used as the value for psi. 

* 'sensitivity_analysis.py': This script runs the sensitivity analysis for all parameters in the models. The parameters are allowed to vary +/-50% of their original values.

* 'sensitivity_analysis.R': THe R script performs a linear regression analysis for each parameter to determine the percentage change in MRSA acquisitions due to a single-percentage change in each parameter value. Figure 5 in the manuscript is also produced from this program. 

* 'gamma_model.py': This script performed a sensitivity analysis between the metapopulation model and an expanded version of the model to include a proportion of time nurses spent in their originally assigned patient groups.

* 'gamma_analysis.R': This R script produces Figure 6 from the manuscript, where a segmented Poisson regression model was fit to detect any thresholds in the value of gamma where it’s relationship to the rate of MRSA acquisitions notably changed. 
