# regPOspline
This R package includes function po_fit that can do regression analysis of arbitrarily censored and left-truncated data under a popular semiparametric proportional odds model. A new estimation approach via an EM algorithm is developed based on a novel data augmentation involving latent exponential and multinomial random variables.

To install the pacakge, please run the following code:
```r
install.packages("devtools")
library(devtools)
remotes::install_github("luwstat/regPOspline", dependencies = TRUE)
library(regPOspline)
```
To know details of the function po_fit, please type the following code in R:
```r
?po_fit

To try the function on the built-in generated data, please try the example provided in R:
```r
?generated_data
