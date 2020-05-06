mlmi implements so called Maximum Likelihood Multiple Imputation as described by von Hippel and Bartlett (2019) https://arxiv.org/abs/1210.0870v10. A number of different imputations are available, by utilising the norm, cat and mix packages. Inferences can be performed either using combination rules similar to Rubin's or using a likelihood score based approach based on theory by Wang and Robins (1998) https://doi.org/10.1093/biomet/85.4.935.

mlmi also implements a maximum likelihood MI version of reference based MNAR imputation for repeatedly measured continuous endpoints.

You can install the released version of bootImpute from CRAN with:
install.packages("mlmi")

And the development version with
install.packages("devtools")
devtools::install_github("jwb133/mlmi")
