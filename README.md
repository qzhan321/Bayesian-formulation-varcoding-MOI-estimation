# Bayesian-formulation-varcoding-MOI-estimation

## Contents
- [Overview and Background](#Overview-and-Background)
  - [varcoding for MOI Estimation and a Bayesian Formulation](#varcoding-for-MOI-Estimation-and-a-Bayesian-Formulation)
- [Applying the Bayesian varcoding to New Datasets](#Applying-the-Bayesian-Formulation-of-MOI-Estimation-to-New-Datasets)
  - [Running the Script](#Running-the-Script)
    - [Pre-processing](#Pre-processing)
    - [Command](#Command)
    - [Example Command](#Example-Command)
    - [Command Arguments](#Command-Arguments)
    - [Output](#Output)
    - [Large dataset](#Large-dataset)
    - [Help](#Help)
- [Contact](#Contact)

## Overview and Background
We summarize the main steps of an extended, Bayesian formulation of *var*coding for MOI (multiplicity of infection) estimation, proposed in [Tiedje and Zhan et al., *eLife*, 2023.](https://doi.org/10.7554/eLife.91411.1) 

The multiplicity of infection (MOI), defined as the number of genetically distinct parasite strains co-infecting a host, is one key epidemiological parameter for measuring malaria transmission and evaluating malaria interventions. Estimating MOI was challenging especially in high-transmission endemic regions where individuals typically carry multiple co-occurring infections, reviewed in [Labbé et al., PLoS Comput Biol, 2023](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010816).   

### varcoding for MOI Estimation and a Bayesian Formulation
The hyper-diversity of the *var* multigene family encoding the major variant surface antigen, *Plasmodium falciparum* erythrocyte membrane protein 1 (*PfEMP*1), and the limited repertoire similarity, provide a viable solution for MOI estimation. A constant repertoire size or number of *var* (non-upsA DBLα types) in a parasite genome can be used to convert the number of types sequenced in an isolate to its MOI ([Ruybal-Pesántez et al., Int. J. Parasitol., 2022](https://www.sciencedirect.com/science/article/pii/S0020751922000030?via%3Dihub); [Tiedje et al., PLOS Glob Public Health, 2022](https://journals.plos.org/globalpublichealth/article?id=10.1371/journal.pgph.0000285)).

Here, we extend the method to a Bayesian formulation which accounts for the measurement error in the repertoire size introduced by targeted PCR and amplicon sequencing of *var* genes in an infection (subsampling of *var* genes) and estimate the posterior distribution for each sampled individual for the probability of different MOI values. From individual posterior distributions, we can then obtain the estimated MOI frequency distribution for the population as a whole.

- **The repertoire size distribution**
P(s|MOI=1): repertoire size distribution, the distribution of the number of non-upsA DBLα types sequenced given MOI = 1, which is empirically available. Examples of this distribution are given in the [section](#Running-the-script).
- **Serial convolutions** of this size distribution: from P(s|MOI=1) to P(s|MOI=2), P(s|MOI=3), etc.
- **Bayes’ rule** to get P(MOI=i|s), which requires the specification of a prior distribution of MOI. The prior reflects your belief on MOI distribution of your sampled population before seeing or analyzing the sequence data, for example, centering around lower values or higher values, which depends on your rough understanding and estimation of the local transmission intensity. We examined negative binomial distributions with a wide range of parameter value and a uniform distribution. The MOI estimation is not sensitive to the prior distribution of MOI for datasets sampled from Bongo District in northern Ghana. We use a uniform prior in [Tiedje and Zhan et al., *eLife*, 2023](https://doi.org/10.7554/eLife.91411.1). Nonetheless we allow users to specify their own prior and associated parameters. Details included in the following [section](#Running-the-script).
- From individual to the population level MOI distribution, we either **pool the maximum a posteriori MOI estimate for each sampled individual**, or **use the technique called mixture distribution**. More details can be found in [Tiedje and Zhan et al., *eLife*, 2023](https://doi.org/10.7554/eLife.91411.1).

## Applying the Bayesian Formulation of MOI estimation to New Datasets
The script **[MOI_estimation.R](https://github.com/qzhan321/Bayesian-formulation-varcoding-MOI-estimation/blob/main/scripts/MOI_estimation.R)** is needed. 

### Running the Script
**MOI_estimation.R** estimates MOI values for individual hosts given their number of non-upsA DBLα types. It requires a .csv file as input. The associated dataframe has two columns, with the first one being host IDs, or any type of identifiers of hosts, and the second column being the number of **non-upsA** DBLα types sequenced and typed in each individual corresponding host. Below is an example of the first two rows of an input matrix:

| HostID | DBLa_upsBC_rep_size |
| :--: | :---------: | 
| `RS1MRS0432.MID76.76.P6.dec15` | 46 |
| `RS1MRS1967.MID88.88.P6.dec15`  | 95 |

It also requires a repertoire size distribution, which is a .csv file (see previous [section](#varcoding-for-MOI-Estimation-and-a-Bayesian-Formulation)). This repertoire size distribution could be in the raw count data format. Below is the first two rows of an example repertoire size distribution:
| DBLa_upsBC_rep_size | n |
| :---------: | :--: | 
| 9 | 2 |
| 10 | 4 |

Among all monoclonal infections in your datasets, 2 of them have 9 **non-upsA** DBLα types sequenced and typed successfully, and 4 of them have 10 **non-upsA** DBLα types sequenced and typed successfully.
Alternatively, the repertoire size distribution could be in a probability format of the raw count data, or even a smoothed version of the raw count data or its corresponing probability format.
| DBLa_upsBC_rep_size | p |
| :---------: | :---------: | 
| 9 | 0.004228330 |
| 10 | 0.002114165 |

Among all monoclonal infections in your datasets, there is a probability of 0.004228330 for having 9 **non-upsA** DBLα types sequenced and typed successfully, and there is probability of 0.002114165 for having 10 **non-upsA** DBLα types sequenced and typed successfully. 
Note that our script will check whether specific values of the number of **non-upsA** DBLα types are associated with n = 0 or p = 0, i.e., no observation of these values for the number of **non-upsA** DBLα types among your monoclonal infections. By default we impute a non-zero n or p based on the mean of the nearest two neighbors' n or p. For example, if for DBLa_upsBC_rep_size = 11, n = 0 (i.e., p = 0), we impute its n or p based on that of upsBC_rep_size = 10 and upsBC_rep_size = 12.
It is recommended to use the collection of monoclonal infection from the empirical surveys whose MOI are to be estimated. 

#### Pre-processing
Depending on the quality of your data, you may wish to perform some pre-processing steps. For example, removing isolates with very few non-upsA DBLα types sequenced. 

#### Command
```bash
Rscript MOI_estimation.R -i "path/to/directory/inputFile" -m 20 -r "path/to/directory/repertoireSizeDistributionFile" -t "count" -p "uniform" -s "medium" -v TRUE -a "pool" -o "/path/to/directory/outFile"
```

#### Example Command 
```bash
Rscript MOI_estimation.R -i "/Users/John/Downloads/survey_1.csv" -m 20 -r "/Users/John/Downloads/repertoireSizeDistribution.csv" -t "count" -p "uniform" -s "NA" -v TRUE -a "mixtureDist" -o "/Users/John/Downloads/survey_1_MOI_estimation_results.RData" 
```

#### Command Arguments
|  Name | Description |
|  :-:  | :---------: | 
|  `i`  | Input. The full path to the input matrix, which is a .csv file format with two columns: host identifiers "HostID", and the number of **non-upsA** DBLα types sequenced in each host "DBLa_upsBC_rep_size". |
|  `m`  | maxMOI. The maximum value for MOI: an estimate from your empirical dataset. All the estimated MOI values will be capped by this parameter. |
|  `r`  | repertoireSizeDist. The full path to the repertoire size distribution, i.e., the distribution of the number of **non-upsA** DBLα types sequenced in monoclonal infections, or its corresponding probablity format, or even a smoothed version of the corresponding probability format. |
|  `t`  | TypeOfRepertoireSizeDist. The type of the repertoire size distribution, which can be the raw count data, or a smoothed probability version of the raw count data. Two options: "count" vs. "probability". |
|  `p`  | The prior distribution for MOI to be estimated; two options, "uniform" or "negBinom" (short for negative binomial) |
|  `s`  | When specifying the prior to be a "negBinom", users need to specify the range for the mean of the negative binomial. There are three options: "low" or "medium" or "high", corresponding to a mean of ~1.5, ~4.3, ~6.7 respectively. If your sample is from a high- or extremely high-transmission endemic region, we recommend you to use either a uniform or a negative binomial with medium or high parameter setting. If your sample is from a low-tranmission region, we recommend you to use either a uniform or a negative binomial with a low parameter setting. For a "uniform" prior, this parameter does not apply and the user can simply specify it to be "NA". |
|  `v`  | Whether output the prior distribution of MOI or not. It is logical, "TRUE" or "FALSE". |
|  `a`  | How to obtain the MOI distribution at the population level from individual MOI estimates, either pooling the maximum a posteriori MOI estimate for each sampled individual or using the technique called mixture distribution, "pool" vs. "mixtureDist". |
|  `o`   | The full path to the directory where the output will be saved and the name of the output file. For example, "MOI_estimation_reuslts.RData". |

#### Output
The above example command will output a list of two objects (plus the prior distribution for MOI estimation when setting **verbose** to be TRUE). When set the argument **aggregate** to be 'pool', the output list contains one matrix which records the maximum a posteriori MOI estimate for each sampled individual, and a second matrix which records the probability distribution at the population level. The matrix at the individual level looks like the example below, with a **prob** column storing the actual probability of MOI = maxAPosMOIEst (the maximum a posteriori MOI estimate). 
| HostID | DBLa_upsBC_rep_size | MOI | Prob |
| :--: | :---------: | :--: | :--: |
| `RS1MRS0432.MID76.76.P6.dec15` | 46 | 2 | 0.9870216 |
| `RS1MRS1967.MID88.88.P6.dec15`  | 95 | 3 | 0.8315176 |

The matrix at the population level:
| MOI | Prob |
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

When set the argument **aggregate** to be 'mixtureDist', the output list again contains two matrices, one recording the full probability distribution of MOI for individual hosts, and a second one recording the probability distribution at the population level.
The matrix for individual hosts (only listing one host as an example):
| HostID | DBLa_upsBC_rep_size | MOI | Prob |
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
| MOI | Prob |
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

MOI distributions at the population level obtained from the two approaches differ slightly. But the difference is non-significant with the example datasets sampled from Bongo District in northern Ghana as determined by the Kolmogorov-Smirnov Test. Users can run both and compare MOI distributions at the population level to inspect the difference between the two. 

#### Large Dataset
You may embed the command line to a bash script and run it on a computational cluster. That way you can request a number of nodes and memory per node for enough computational power.

#### Help
Run the command below to print out help page.
```bash
Rscript MOI_estimation.R --help
```
Users can refer to the help page for the definition of each parameter, their default values, and all the possible options for their values.

## Contact
If you run into any issues using the method, feel free to open a new [issue](https://github.com/qzhan321/Bayesian-formulation-varcoding-MOI-estimation/issues).
