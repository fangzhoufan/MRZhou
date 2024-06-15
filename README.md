# MRZhou

This R package is designed for the convenience of medical professionals, aiming to further simplify the analysis methods of Mendelian randomization (MR). It is integrated based on the [TwosampleMR](https://github.com/MRCIEU/TwoSampleMR/), [MRpresso](https://github.com/rondolab/MR-PRESSO), and [coloc](https://github.com/chr1swallace/coloc) packages. Therefore, whether it is the online data of the [IEU GWAS database](https://gwas.mrcieu.ac.uk/) or local GWAS data (such as [FinnGen database](https://www.finngen.fi/en)), this package can be conveniently executed. This package also includes: [GSMR analysis](https://yanglab.westlake.edu.cn/software/gsmr/), [SMR analysis](https://yanglab.westlake.edu.cn/software/smr/) methods, etc. [TWAS/Fusion analysis](http://gusevlab.org/projects/fusion/) function method is under development.

When using this [MRZhou](https://github.com/fangzhoufan/MRZhou) package, it is recommended to first learn the basic syntax of [TwosampleMR](https://github.com/MRCIEU/TwoSampleMR/).

## Installation

You can install the development version of MRZhou from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fangzhoufan/MRZhou")
```

## Example

This is a basic example which shows you how to solve a common problem:

### 1.Exposure and Outcome

Let's identify the instrumental variables (IVs) for BMI and extract the BMI-associated SNPs from the CAD (coronary artery disease) outcome GWAS data.

```{r example}
library(MRZhou)
exp <- extract_instruments(outcomes = "ieu-a-2")
out <- extract_instruments(outcomes = "ieu-a-7")
```

> The above two functions are based on the [TwosampleMR](https://github.com/MRCIEU/TwoSampleMR/) package, and detailed information can be found in the tutorial.

If you would like to use the FinnGen database, you can follow: exposure:

-   Exposure: [BMI](https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_BMI_IRN.gz)

-   Outcome: [CAD](https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_IHD.gz)

```{r example}
library(MRZhou)
data("finngenR9")
data("finngenR10")
# exp<-read_finngen_exposure('BMI','finngen_R10_BMI_IRN.gz', 290820, p_threshold = 5e-08)
# out<-read_finngen_outcome(exp, 'CAD', 'finngen_R10_I9_IHD.gz', 412181)
```

### 2. IVW analysis (main analysis method)

```{r }
dat <- harmonise_data(exp, out)
result <- IVW_fix_random(dat)
test <- mr_test(dat)
presso <- MRpresso(dat)
resall <- combine_results(result, test, presso)
```

The resall dataframe is all of the MR methods results (including: MRpresso), Heterogeneity statistics and Horizontal pleiotropy results

### 3.Plot

```{r }
## Not include MRpresso method results
forest_plot <- forest_plot(result)

## Include MRpresso method results
result2 <- combine_results(result, NULL, presso)
forest_plot <- forest_plot(result2)
```

### 4.Other Visualization

-   Volcano plot

```{r }
volcano_plot <- volcano_plot(result)
```

-   Circle plot

```{r }
circle_plot <- circle_plot(result)
```

