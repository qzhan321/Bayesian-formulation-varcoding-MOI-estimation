# Bayesian-formulation-varcoding-MOI-estimation

## Contents
* [Overview](#Overview)
* [varcoding for MOI Estimation and Bayesian Formulation](#varcoding-for-MOI-Estimation-and-Bayesian-Formulation)
* [Applying the Bayesian Formulation to New Dataset](#Applying-the-Bayesian-Formulation-to-New-Dataset)


## Overview
We summarize the main steps of Bayesian approach for MOI (multiplicity of infection) estimation, and more details can be found in [Tiedje K.E., Zhan Q., Ruybal-Pesántez S., Tonkin-Hill G., He Q., Tan M.H., Argyropoulos D.C., Deed S.L., Ghansah A., Bangre O., Oduro A.R., Koram K.A., Pascual M., Day K.P., 2023. Measuring changes in Plasmodium falciparum var census population size and structure in response to sequential malaria control interventions. medRxiv. https://doi.org/10.1101/2023.05.18.23290210.](https://www.medrxiv.org/content/early/2023/05/19/2023.05.18.23290210.full.pdf) 

The multiplicity of infection (MOI), defined as the number of genetically distinct parasite strains co-infecting a host, is one key epidemiological parameter for measuring malaria transmission and evaluating malaria interventions. Estimating MOI remains challenging especially in high-transmission endemic settings where individuals typically carry multiple co-occurring infections, recently reviewed in [Labbé et al., 2023](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010816).   

## varcoding for MOI Estimation and Bayesian Formulation
The hyper-diversity of *var*(DBLα types) and limited repertoire overlap of the *var* multigene family encoding the major *Plasmodium falciparum* blood-stage antigen *PfEMP*1 provide a viable solution for MOI estimation. A constant repertoire size or number of DBLα types in a parasite genome can be used to convert the number of types sequenced in an isolate to its estimated MOI ([Ruybal-Pesántez et al., 2022](https://www.sciencedirect.com/science/article/pii/S0020751922000030?via%3Dihub); [Tiedje et al., 2022](https://journals.plos.org/globalpublichealth/article?id=10.1371/journal.pgph.0000285)).

Here, we extend the method to a Bayesian formulation which accounts for the measurement error in the repertoire size introduced by targeted PCR
and amplicon sequencing of *var* genes in an infection (subsampling of *var* genes) and estimate the posterior distribution for each sampled individual for the probability of different MOI values. From individual posterior distributions, we can then obtain the estimated MOI frequency distribution for the population as a whole.

- The repertoire size distribution
P(s|MOI=1): the distribution of the number of non-upsA DBLα types sequenced given MOI = 1, which is empirically available.
- Serial convolutions of this distribution: from P(s|MOI=1) to P(s|MOI=2), P(s|MOI=3), etc.
- Bayes’ rule to get P(MOI=i|s), which requires the specification of a prior distribution of MOI. We examined negative binomial distributions with a wide range of parameter value and a uniform distribution. The MOI estimation is not sensitive to the prior distribution of MOI. We use a uniform prior in this work. 
- From individual to the population level MOI distribution, we either pool the maximum a posteriori MOI estimate for each sampled individual, or use the technique called mixture distribution.

## Applying the Bayesian Formulation to New Dataset
There are two R scripts in the folder **[scripts](https://github.com/qzhan321/Bayesian-formulation-varcoding-MOI-estimation/tree/main/scripts)**.  To apply our Bayesian method to a new dataset, only the script **[MOI_estimation.R](https://github.com/qzhan321/Bayesian-formulation-varcoding-MOI-estimation/blob/main/scripts/MOI_estimation.R)** is needed. The other script **[derive-prob-s-given-MOI.R](https://github.com/qzhan321/Bayesian-formulation-varcoding-MOI-estimation/blob/main/scripts/derive-prob-s-given-MOI.R)** walks through how we derive P(s|MOI=2) and P(s|MOI=3) and so on up until P(s|MOI=20) (20 being the upper limit of MOI values, empirical determined) based on the serial convolutions of the repertoire size distribution P(s|MOI=1). But we upload the output list which stores the probability distribution of s, i.e., the number of DBLa types sequenced and typed, given any certain MOI **[s_givenMOI_list](https://github.com/qzhan321/Bayesian-formulation-varcoding-MOI-estimation/blob/main/scripts/s_givenMOI_list)**, which can be downloaded and used directly. The first element of this list corresponds to the distribution of s given MOI = 1, and the second element of this list correspondds to the distribution of s given MOI = 2, and so on so forth. 



