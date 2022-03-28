# Author's code to 'Vulnerability-CoVaR: Investigating the Crypto-market'

This repository contains the source [R](https://www.r-project.org/) code of:

Martin Waltz, Abhay Kumar Singh, Ostap Okhrin (2022). *Vulnerability-CoVaR: Investigating the Crypto-market*.

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
@misc{2022VCoVaR,
  author = {Waltz, Martin and Singh, Abhay Kumar and Okhrin, Ostap},
  title = {Vulnerability-CoVaR: Investigating the Crypto-market},
  year = {2022},
  publisher = {GitHub},
  journal = {GitHub Repository},
  howpublished = {\url{https://github.com/MarWaltz/VCoVaR}}
}
~~~

Note: This citation will be updated once the paper is published.