# Author's code to 'Vulnerability-CoVaR: Investigating the Crypto-market'

This repository contains the source [R](https://www.r-project.org/) code of:

Martin Waltz, Abhay Kumar Singh, Ostap Okhrin (2022). Vulnerability-CoVaR: investigating the crypto-market, *Quantitative Finance*, 
DOI: [10.1080/14697688.2022.2063166](https://www.tandfonline.com/doi/full/10.1080/14697688.2022.2063166)

The used Community Network Data is kindly provided by [CoinMetrics](https://coinmetrics.io/). Please note that no warranties are given for the functionality and correctness of the code.

## Prerequisites

To run the code, the following packages are needed on top of base [R](https://www.r-project.org/):
[binaryLogic](https://cran.r-project.org/web/packages/binaryLogic/index.html), [copula](https://cran.r-project.org/web/packages/copula/index.html), [HAC](https://cran.r-project.org/web/packages/HAC/index.html), [nloptr](https://cran.r-project.org/web/packages/nloptr/index.html), [plotly](https://cran.r-project.org/web/packages/plotly/index.html), [rmgarch](https://cran.r-project.org/web/packages/rmgarch/index.html), [rugarch](https://cran.r-project.org/web/packages/rugarch/index.html), [timeDate](https://cran.r-project.org/web/packages/timeDate/index.html), [tseries](https://cran.r-project.org/web/packages/tseries/index.html), [WeightedPortTest](https://cran.r-project.org/web/packages/WeightedPortTest/index.html).

We recommend using version 1.0.0. of the package `copula` since we had some minor issues in version 1.0.1. Each of the above mentioned packages can be installed via, e.g.:
```
install.packages("HAC")
```

## Citation

If you use this code in one of your projects or papers, please cite it as follows.

~~~bibtex
@article{doi:10.1080/14697688.2022.2063166,
author = {Martin Waltz and Abhay Kumar Singh and Ostap Okhrin},
title = {Vulnerability-CoVaR: investigating the crypto-market},
journal = {Quantitative Finance},
year  = {2022},
publisher = {Routledge},
doi = {10.1080/14697688.2022.2063166},
URL = {https://doi.org/10.1080/14697688.2022.2063166}
}
~~~