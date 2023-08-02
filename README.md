# Bayesian-formulation-varcoding-MOI-estimation

## Contents
* [Overview and Background](#Overview-and-Background)
* [Applying the Bayesian Formulation of MOI Estimation to New Datasets](#Applying-the-Bayesian-Formulation-of-MOI-Estimation-to-New-Datasets)


## Overview and Background
We summarize the main steps of a Bayesian approach for MOI (multiplicity of infection) estimation, proposed in [Tiedje K.E., Zhan Q., Ruybal-Pesántez S., Tonkin-Hill G., He Q., Tan M.H., Argyropoulos D.C., Deed S.L., Ghansah A., Bangre O., Oduro A.R., Koram K.A., Pascual M., Day K.P., 2023. Measuring changes in Plasmodium falciparum var census population size and structure in response to sequential malaria control interventions. medRxiv. https://doi.org/10.1101/2023.05.18.23290210.](https://www.medrxiv.org/content/early/2023/05/19/2023.05.18.23290210.full.pdf) 

The multiplicity of infection (MOI), defined as the number of genetically distinct parasite strains co-infecting a host, is one key epidemiological parameter for measuring malaria transmission and evaluating malaria interventions. Estimating MOI remains challenging especially in high-transmission endemic settings where individuals typically carry multiple co-occurring infections, recently reviewed in [Labbé et al., 2023](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010816).   

### varcoding for MOI estimation and a Bayesian formulation
The hyper-diversity of *var*(DBLα types) and limited repertoire overlap of the *var* multigene family encoding the major *Plasmodium falciparum* blood-stage antigen *PfEMP*1 provide a viable solution for MOI estimation. A constant repertoire size or number of non-upsA DBLα types in a parasite genome can be used to convert the number of types sequenced in an isolate to its estimated MOI ([Ruybal-Pesántez et al., 2022](https://www.sciencedirect.com/science/article/pii/S0020751922000030?via%3Dihub); [Tiedje et al., 2022](https://journals.plos.org/globalpublichealth/article?id=10.1371/journal.pgph.0000285)).

Here, we extend the method to a Bayesian formulation which accounts for the measurement error in the repertoire size introduced by targeted PCR
and amplicon sequencing of *var* genes in an infection (subsampling of *var* genes) and estimate the posterior distribution for each sampled individual for the probability of different MOI values. From individual posterior distributions, we can then obtain the estimated MOI frequency distribution for the population as a whole.

- **The repertoire size distribution**
P(s|MOI=1): the distribution of the number of non-upsA DBLα types sequenced given MOI = 1, which is empirically available.
- **Serial convolutions** of this distribution: from P(s|MOI=1) to P(s|MOI=2), P(s|MOI=3), etc.
- **Bayes’ rule** to get P(MOI=i|s), which requires the specification of a prior distribution of MOI. We examined negative binomial distributions with a wide range of parameter value and a uniform distribution. The MOI estimation is not sensitive to the prior distribution of MOI. We use a uniform prior in [Tiedje and Zhan et al, 2023](https://www.medrxiv.org/content/early/2023/05/19/2023.05.18.23290210.full.pdf). Nonetheless we allow users to specify their own prior and associated parameters. Details included in the following [section](#Running-the-script).
- From individual to the population level MOI distribution, we either **pool the maximum a posteriori MOI estimate for each sampled individual**, or **use the technique called mixture distribution**. More details can be found in [Tiedje and Zhan et al, 2023](https://www.medrxiv.org/content/early/2023/05/19/2023.05.18.23290210.full.pdf)

## Applying the Bayesian Formulation of MOI estimation to New Datasets
There are two R scripts in the folder **[scripts](https://github.com/qzhan321/Bayesian-formulation-varcoding-MOI-estimation/tree/main/scripts)**. To apply our Bayesian method to a new dataset, only the script **[MOI_estimation.R](https://github.com/qzhan321/Bayesian-formulation-varcoding-MOI-estimation/blob/main/scripts/MOI_estimation.R)** is needed. 

The other script **[derive-prob-s-given-MOI.R](https://github.com/qzhan321/Bayesian-formulation-varcoding-MOI-estimation/blob/main/scripts/derive-prob-s-given-MOI.R)** walks through how we derive the probability distribution of the number of non-upsA DBLα types sequenced and typed given any MOI value, i.e., P(s|MOI=2) and P(s|MOI=3) and so on up until P(s|MOI=20) (20 being the upper limit of MOI values, empirically determined) based on the serial convolutions of the repertoire size distribution P(s|MOI=1). But we upload the output list to the **[scripts](https://github.com/qzhan321/Bayesian-formulation-varcoding-MOI-estimation/tree/main/scripts)** folder, which stores the probability distribution of s given any certain MOI, i.e., **[s_givenMOI_list](https://github.com/qzhan321/Bayesian-formulation-varcoding-MOI-estimation/blob/main/scripts/s_givenMOI_list)**. It can be downloaded and used directly. The first element of this list corresponds to the distribution of s given MOI = 1, and the second element of this list corresponds to the distribution of s given MOI = 2, and so on so forth. 

### Running the script
**MOI_estimation.R** estimates MOI values for individual hosts given their number of non-upsA DBLα types. It takes a matrix as its input. The matrix has two columns, with the first one being host IDs, or identifiers of hosts, and the second column being the number of **non-upsA** DBLα types sequenced and typed in each individual corresponding host. Below is an example of the first two rows of an input matrix:

| HostID | NumDBLaTypes |
| :--: | :---------: | 
| `RS1MRS0432.MID76.76.P6.dec15` | 46 |
| `RS1MRS1967.MID88.88.P6.dec15`  | 95 |

#### Pre-processing
We recommend to perform a preprocessing step on the input matrix, removing isolates with either fewer than 10 or more than 900 non-upsA DBLα types sequenced. MOI estimates are capped at 20.

#### Command
```bash
Rscript MOI_estimation.R --input "path/to/directory/inputFile" --aggregate "pool" --util "/path/to/directory/utilFile" --output "/path/to/directory/outFile"
```
We can write console output to a text file by adding the following at the end of the command:
```bash
>consoleOutput.txt
```
Right now the output text file does not contain much information. You may modify your **MOI_estimation.R** by adding printing command to inspect input or intermediary variables.

#### Example Command 
```bash
Rscript MOI_estimation.R --input "/Users/John/Downloads/survey_1.csv" --aggregate "pool" --util "/Users/John/Downloads/s_givenMOI_list" --output "/Users/John/Downloads/survey_1_MOI.RData" >consoleOutput.txt
```

#### Command arguments
| Name | Description |
| :--: | :---------: | 
| `input` | The full path to the input matrix: both .csv and .txt file formats are acceptable |
| `aggregate`  | How to obtain the MOI distribution at the population level from individual MOI estimates, either pooling the maximum a posteriori MOI estimate for each sampled individual or using the technique called mixture distribution, "pool" vs. "mixtureDist" |
| `util`  | Local path to the downloadable object **[s_givenMOI_list](https://github.com/qzhan321/Bayesian-formulation-varcoding-MOI-estimation/blob/main/scripts/s_givenMOI_list)** |
| `output`  | Path to the directory where the output will be saved and the name of the output file (for example, in the .RData format) |

#### Output
The above example command will output a list of two objects. When set the argument **fromIndividualToPop** to be 'pool', the output list contains one matrix which records the maximum a posteriori MOI estimate for each sampled individual, and a second matrix which records the probability distribution at the population level. The matrix at the individual level looks like the example below, with a **prob** column storing the actual probability of MOI = maxAPosMOIEst (the maximum a posteriori MOI estimate). 
| HostID | NumDBLαTypes | maxAPosMOIEst | prob |
| :--: | :---------: | :--: | :--: |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 2 | 0.9870216 |
| `RS1MRS1967.MID88.88.P6.dec15`  | 95 | 3 | 0.8315176 |

The matrix at the population level:
| MOI | prob |
| :--: | :--: |
| 1 | 0.108029197 |
| 2 | 0.224817518 |
| 3 | 0.144525547 |
| 4 | 0.135766423 |
| 5 | 0.112408759 |
| 6 | 0.080291971 |
| 7 | 0.055474453 |
| 8 | 0.042335766 |
| 9 | 0.024817518 |
| 10 | 0.030656934 |
| 11 | 0.010218978 |
| 12 | 0.011678832 |
| 13 | 0.008759124 |
| 14 | 0.001459854 |
| 15 | 0.002919708 |
| 16 | 0.002919708 |
| 17 | 0 |
| 18 | 0 |
| 19 | 0 |
| 20 | 0.002919708 |

When set the argument **fromIndividualToPop** to be 'mixtureDist', the output list again contains two matrices, one recording the full probability distribution of MOI for individual hosts, and a second one recording the probability distribution at the population level.
The matrix for individual hosts (only listing one host as an example):
| HostID | NumDBLαTypes | MOI | prob |
| :--: | :---------: | :--: | :--: |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 1 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 2 | 0.9870216 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 3 | 0.01297535 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 4 | 3.035835e-06 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 5 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 6 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 7 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 8 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 9 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 10 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 11 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 12 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 13 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 14 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 15 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 16 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 17 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 18 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 19 | 0 |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 20 | 0 |

The matrix at the population level:
| MOI | prob |
| :--: | :--: |
| 1 | 0.0994607006 |
| 2 | 0.2159426383 |
| 3 | 0.1522643518 |
| 4 | 0.1259793933 |
| 5 | 0.1126470417 |
| 6 | 0.0870748021 |
| 7 | 0.0595727360 |
| 8 | 0.0420562894 |
| 9 | 0.0311597868 |
| 10 | 0.0232443687 |
| 11 | 0.0163802735 |
| 12 | 0.0119803323 |
| 13 | 0.0080522791 |
| 14 | 0.0047072134 |
| 15 | 0.0029879126 |
| 16 | 0.0020022142 |
| 17 | 0.0011265310 |
| 18 | 0.0006744842 |
| 19 | 0.0006565408 |
| 20 | 0.0020301102|

MOI distributions at the population level obtained from the two approaches differ slightly. But the difference is non-significant with the example datasets from northern Ghana as determined by the Kolmogorov-Smirnov Test. 
