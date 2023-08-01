# Bayesian-formulation-varcoding-MOI-estimation

## Contents
* [Overview](#-Overview)


## Overview
We summarize the main steps of Bayesian approach for MOI (multiplicity of infection) estimation, and more details can be found in [Tiedje K.E., Zhan Q., Ruybal-Pesántez S., Tonkin-Hill G., He Q., Tan M.H., Argyropoulos D.C., Deed S.L., Ghansah A., Bangre O., Oduro A.R., Koram K.A., Pascual M., Day K.P., 2023. Measuring changes in Plasmodium falciparum var census population size and structure in response to sequential malaria control interventions. medRxiv. https://doi.org/10.1101/2023.05.18.23290210.](https://www.medrxiv.org/content/early/2023/05/19/2023.05.18.23290210.full.pdf) 

The multiplicity of infection (MOI), defined as the number of genetically distinct parasite strains co-infecting a host, is one key epidemiological parameter for measuring malaria transmission and evaluating malaria interventions. Estimating MOI remains challenging especially in high-transmission endemic settings where individuals typically carry multiple co-occurring infections, recently reviewed in [Labbé, F., He, Q., Zhan, Q., Tiedje, K.E., Argyropoulos, D.C., Tan, M.H., Ghansah, A., Day, K.P., Pascual, M., 2023. Neutral vs. non-neutral genetic footprints of Plasmodium falciparum multiclonal infections. PLOS Computational Biology 19, e1010816. https://doi.org/10.1371/journal.pcbi.1010816](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010816).   

## varcoding for MOI estimation and Bayesian formulation
The hyper-diversity of *var*(DBLα types) and limited repertoire overlap of the *var* multigene family encoding the major *Plasmodium falciparum* blood-stage antigen *PfEMP*1 provide a viable solution for MOI estimation. A constant repertoire size or number of DBLα types in a parasite genome can be used to convert the number of types sequenced in an isolate to its estimated MOI ([Ruybal-Pesántez et al., 2022](https://www.sciencedirect.com/science/article/pii/S0020751922000030?via%3Dihub); [Tiedje et al., 2022](https://journals.plos.org/globalpublichealth/article?id=10.1371/journal.pgph.0000285)).

Here, we extend the method to a Bayesian formulation which accounts for the measurement error in the repertoire size introduced by targeted PCR
and amplicon sequencing of *var* genes in an infection (subsampling of *var* genes) and estimate the posterior distribution for each sampled individual for the probability of different MOI values. From individual posterior distributions, we can then obtain the estimated MOI frequency distribution for the population as a whole.

- The repertoire size distribution
P(s|MOI=1): the distribution of the number of non-upsA DBLα types sequenced given MOI = 1, which is empirically available.
- Serial convolutions of this distribution: from P(s|MOI=1) to P(s|MOI=2), P(s|MOI=3), etc.
- Bayes’ rule to get P(MOI=i|s), which requires the specification of a prior distribution of MOI. We examined negative binomial distributions with a wide range of parameter value and a uniform distribution. The MOI estimation is not sensitive to the prior distribution of MOI. We use a uniform prior in this work. 
- From individual to the population level MOI distribution, we either pool the maximum a posteriori MOI estimate for each sampled individual, or use the technique called mixture distribution.

## Applying the Bayesian formulation to your dataset




