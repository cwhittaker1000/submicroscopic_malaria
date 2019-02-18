# Submicroscopic Malaria Infections

This repository contains code used to generate the figures for the (submitted) publication "Variation in the Prevalence of Submicroscopic Malaria Infections: Historical Transmission Intensity and Age as Key Determinants. Specifically, it uses a database of collated malaria cross-sectional survey data (where populations have had infection status determined by both Light Microscopy and PCR) to:

- Explore how the proportion of submicroscopically infected individuals varies across the transmission gradient (Figure 1)
- Examine age-related patterns of infection and see how the size of the submicroscopic reservoir varies acrosss different age groups (Figure 2)
- Assess the extent of geographical variation in the proportion of submicroscopically infected individuals (Figure 3) 
- Explore how historical patterns of transmission influence the extent and size of the submicroscopic reservoir (Figure 4)
- Integrate these results with previous estimates of the comparative infectivity (to mosquitoes) of microscopically detectable and submicroscopic infections to generate estimates for the contribution to onwards transmission of these two infected subgroups (Figure 5) 

## Description
A systematic literature review was carried out to compile malaria prevalence data from references where both microscopy and PCR based methods had been used to determine infection status. This represents an update to the systematic review carried out by Okell et al., in 2012 (Nature Communications, 3, p1237). Searches were carried out using PubMed and WebofScience and the search terms (“PCR” OR “Polymerase Chain Reaction”) AND “falciparum”. 

Following the screening process, a total of 55 new references were included, with a further 2 relevant references, unpublished at the time of screening, but since published, included. Together with the results of the previous systematic review, this gave atotal of 282 prevalence survey pairs. For each prevalence survey pair, a variety of different relevant data were extracted, including (amongst others):

- Number sampled and number positive (for both Light Microscopy and PCR)
- Location and Year of Sampling
- Historical and Current Patterns of Transmission (for African surveys only, using results from the Malaria Atlas Project)
- The Age Range of the Sampled Population

This collated information comprises the data from which the analyses presented in the paper are derived. This data is available in this github repo, filename: "For Submission Whittaker et al R Import.csv"

# Repo Contents

- [R](./R): `R` code required to generate the figures featured in the paper.
- [JAGS_Model](./JAGS_Model): Code specifying the linear (on the logit scale) regression model fitted to the data using the JAGS MCMC software.   
- [Data](./Data): Containing the raw data used in the analyses presented here. This represents the data collated as part of this systematic review, as well as data from Okell et al., 2012. Available in both the raw form, as well as the final version that was imported into `R` for analysis.    


# Software Requirements

Running the code contained in this repository requires the following:

- The package **ssa** (see: https://cran.r-project.org/web/packages/ssa/index.html)
- The package **binom** (see: https://cran.r-project.org/web/packages/binom/index.html)
- The package **rjags** (see: https://cran.r-project.org/web/packages/rjags/index.html)
- The software **JAGS** (see: http://mcmc-jags.sourceforge.net/)

rjags is an `R` package offering an interface to the JAGS MCMC library. JAGS is a program for analysis of Bayesian models using Markov Chain Monte Carlo (MCMC) simulation. More information and details about package itself are available here: http://mcmc-jags.sourceforge.net/ and the software itself is available for download at https://sourceforge.net/projects/mcmc-jags/. 


