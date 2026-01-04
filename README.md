# Weber_predator_prey_spatial_game

Code for "Unveiling the predator prey spatial game using multiple habitat selection functions"

## Updates
Newer versions or comments regarding this source code may be available, please check [here](https://github.com/abbweber/predator-prey-spatial-game-multiple-HSFs) to ensure you are using the latest version.

## Authors and Affliations

Abigail M. Weber<sup>1,2<sup>; Nicole T. Gorman^1,2^; Michael E. Egan^1,2^; Michael W. Eichholz^1,2^; Peter E. Schlichting^3^; Daniel J. Skinner^3^; Guillaume Bastille-Rousseau^1,2^

^1^Center for Wildlife Sustainability Research, Southern Illinois University, Carbondale, Illinois, USA
^2^School of Biological Sciences, Southern Illinois University, Carbondale, Illinois, USA
^3^Illinois Department of Natural Resources, Springfield, Illinois, USA

## Current Affiliations

Abigail M. Weber: Pennsylvania Cooperative Fish and Wildlife Research Unit, Pennsylvania State University, University Park, Pennsylvania, USA

Nicole T. Gorman: Department of Fish and Wildlife Conservation, Virginia Tech, Blacksburg, Virginia, USA

## Point of contact: 

Abigail M. Weber (amw8076@psu.edu) 

## Background

We provide the following R Scripts corresponding to our analyses of predator prey spatial interactions in Illinois, USA. This GitHub [repository](https://github.com/abbweber/predator-prey-spatial-game-multiple-HSFs) includes all R scripts used for analyses in the manuscript. The two data .csv's are stored on Zenodo [here](10.5281/zenodo.18147653). 

We include scripts for the analysis of data and generation of figures. GPS location data for bobcats, coyotes, and white-tailed female deer are available per reasonable request, please contact Guillaume Bastille-Rousseau (gbr@siu.edu) for more information.

## Script names and descriptions

"01_Model_fitting.R": This RScript contains the necessary steps for running the four model sets which highlight the different aspects of the predator-prey spatial game. 

Specifically, the landscape perspective model compared fawn kill site locations to 10,000 random locations (therefore a Resource Selection Function or RSF) and was fitted using generalized linear models (GLM) with a binomial distribution. The prey perspective model compared fawn kill sites to female deer locations (therefore a Latent Selection Difference Function or LSDF) and was fitted using a GLM with a clustered sandwich estimator grouped by individual which assumed that data across, but not within, an individual were independent (Latham et al. 2011). For the predator perspective model, we compared fawn kill site locations to predator locations (therefore a LSDF) which was also fitted as a GLM with a sandwich estimator clustered by individual. The predator selection model compared predator locations to 2,500 random locations resampled with replacement from the original set of 10,000 random locations (therefore a RSF). For the fourth model, we used a generalized estimating equation (GEE) to directly estimate marginal effects, representing the population-averaged response for predator selection (Koper and Manseau 2009). 

We ran ten global models across spatial buffers for the landscape perspective model, and for model selection, we focused on the three comparisons which included fawn kill sites (i.e. landscape, prey, and predator perspectives) because we were most interested in variables impacting fawn mortality risk.

All data can be found in the folder "data". 

+ Input data file: 
  + "data/data_for_HSF_models.csv"
  + "data/data_for_pred_GEE_model.csv"
+ Output data files: 
  + "output/sum_loglikelihood_fawn_models.csv"
  + "output/all_model_coefs_100m_buffer.csv"

"02_Figure_generation.R": This markdown contains the necessary steps for generating results figures. Coefficient data can be found in the folder "output" and figure 2 can be found in the folder "figures". 

+ Input data file: "output/all_model_coefs_100m_buffer.csv"
+ Output data files: "figures/Fig_2.jpg"
