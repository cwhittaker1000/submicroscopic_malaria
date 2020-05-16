# Submicroscopic Malaria Infection Systematic Review

# Contents

- [Overview](#Overview)
- [Description](#Description)
- [Repo Contents](#Repo-Contents)
- [Software Requirements](#Software-Requirements)
- [System Requirements](#System-Requirements)
- [Installation Guide and Instructions for Use](#Installation-Guide-and-Instructions-for-Use)

# Overview
This repository contains code used to generate the figures for the (submitted) publication:

**Global & Temporal Patterns of Submicroscopic Plasmodium falciparum Malaria Infection**   
*Charles Whittaker<sup>1</sup>, Hannah Slater<sup>1,2</sup>, Teun Bousema<sup>3</sup>, Chris Drakeley<sup>4</sup>, Azra Ghani<sup>1</sup> and Lucy Okell<sup>1</sup>*  
*<sup>1</sup>MRC Centre for Global Infectious Disease Analysis, Department of Infectious Disease Epidemiology, Imperial College London, London, United Kingdom*
*<sup>2</sup>PATH, 2201 Westlake Avenue, Seattle, USA*
*<sup>3</sup>Department of Medical Microbiology, Radboud University Medical Center, Nijmegen, The Netherlands*
*<sup>4</sup>Department of Immunology and Infection, Faculty of Infectious and Tropical Diseases, London School of Hygiene and Tropical Medicine, London, United Kingdom*  

Specifically, it uses a database of collated malaria cross-sectional survey data (where populations have had infection status determined by both Light Microscopy and PCR) to:

- Explore how the proportion of submicroscopically infected individuals varies across the transmission gradient **(Figure 1)**
- Assess the extent of geographical variation in the proportion of submicroscopically infected individuals **(Figure 2)** 
- Explore how historical patterns of transmission influence the extent and size of the submicroscopic reservoir **(Figure 3)**
- Examine age-related patterns of infection and see how the size of the submicroscopic reservoir varies acrosss different age groups **(Figure 4)**
- Integrate these results with previous estimates of the comparative infectivity (to mosquitoes) of microscopically detectable and submicroscopic infections to generate estimates for the contribution to onwards transmission of these two infected subgroups **(Figure 5)** 

- It also contains code required to implement and run all analyses contained within the Supplementary Information of the publication. 

# Description
A systematic literature review was carried out to compile malaria prevalence data from references where both microscopy and PCR based methods had been used to determine infection status. This represents an update to the systematic review carried out by Okell et al., in 2012 (Nature Communications, 3, p1237). Searches were carried out using PubMed and WebofScience and the search terms (“PCR” OR “Polymerase Chain Reaction”) AND “falciparum”. 

Following the screening process, a total of 60 new references were included. Together with the results of the previous systematic review, this gave a total of 387 prevalence survey pairs. For each prevalence survey pair, a variety of different relevant data were extracted, including (amongst others):

- Number sampled and number positive (for both Light Microscopy and PCR) **(from the papers themselves)**
- Location and Year of Sampling **(from the papers themselves)**
- Geographical Location **(from the papers themselves)**
- Historical and Current Patterns of Transmission **(for African surveys only, using results from the Malaria Atlas Project: https://map.ox.ac.uk/)**
- The Age Range of the Sampled Population **(from the papers themselves)**
- Various Indicators of Diagnostic Quality (Type of PCR, Number of Microscopy Fields Examined etc) **(from the papers themselves)**
- The Season in Which Sampling Was Carried Out (Either Rainy or Dry) **(from the papers themselves)**

This collated information comprises the data from which the analyses presented in the paper are derived.

# Repo Contents

- [Analyses](./Analyses): `R` code required to generate the figures featured in the paper.
- [JAGS_Model](./JAGS_Model): Code specifying the linear (on the logit scale) regression model fitted to the data using the JAGS MCMC software.   
- [Data](./Data): Containing the raw data used in the analyses presented here. This represents the data collated as part of this systematic review, as well as data from Okell et al., 2012. Available in both the raw form, as well as the version that was imported into `R` for analysis. 
- [Functions](./Functions): Containing helper-functions used to run frequently occurring analyses (such as running the JAGS model, extracting parameter value estimates etc)
- [Figures](./Figures): Containing .PDF and .ai versions of paper figures.
- [Outputs](./Outputs): Containing .rds outputs from MCMC model running in JAGS for each of the analyses conducted.


# System Requirements

Running this code requirs only a standard computer with enough RAM to support contained within this code. The machine this code has been ran and tested on has the following proprties:

Operating System: Windows 7 Entreprise   
Processor: Intel(R) Core(TM) i7-3770 CPU @ 3.40 GHz   
Installed Memory (RAM): 8.00 GB  
System Type: 64-bit Operating System  

This code has been tested only on Windows systems.

# Software Requirements

Running the code contained in this repository requires the following

- The R Programming Language (Version 3.5.2 "Eggshell Igloo" used here) (see: https://www.r-project.org/)
- The package **ssa** (Version 1.3.0 used here) (see: https://cran.r-project.org/web/packages/ssa/index.html)
- The package **binom** (Version 1.1.1 used here) (see: https://cran.r-project.org/web/packages/binom/index.html)
- The package **rjags** (Version 4.8 used here) (see: https://cran.r-project.org/web/packages/rjags/index.html)
- The software **JAGS** (version 4.3.0 used here) (see: http://mcmc-jags.sourceforge.net/)

rjags is an `R` package offering an interface to the JAGS MCMC library. JAGS is a program for analysis of Bayesian models using Markov Chain Monte Carlo (MCMC) simulation. More information and details about package itself are available here: http://mcmc-jags.sourceforge.net/ and the software itself is available for download at https://sourceforge.net/projects/mcmc-jags/ and takes on the order of minutes to install. 

# Installation Guide and Instructions for Use
The following instructions require that all packages required have been installed within `R` and that additionally, the JAGS software has also been installed. For more information about how and where to download these, see above in the Software Requirements section. 
To replicate and reproduce the analyses presented in this paper, do the following: 

1. Download `SI_Systematic_Review_Results_R_Import.csv` from the [Data](./Data) folder of this repository. 
2. Download the `R` code from  [Analyses](./Analyses) for the particular part of the analysis you are trying to reproduce. 
3. From an `R` session, import `SI_Systematic_Review_Results_R_Import.csv`.
4. Run the `R` code. The output from running this code will be a number of MCMC objects, as well as a series of plots representing the output from MCMC based fitting of the logit-linear model to the collated data. These plots form the basis of the figures presented in the publication. 

Note: The exact runtime depends on the size of the dataset being analysed. It is therefore longer for the full dataset (267 datapoint pairs, on the order of 30 minutes) compared to the age-disaggregated datasets (approx. 30-40 datapoints, on the order of 1 minute). 



