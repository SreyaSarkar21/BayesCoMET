# CoMET: A Compressed Bayesian Mixed-Effects Model for High-Dimensional Tensors

This repositary contains:

1. **R package BayesCoMET**
    - Implements the Compressed Bayesian Mixed-Effects Model for Tensors (CoMET).
    - Install this library with:
    ```r
    devtools::install_github("SreyaSarkar21/BayesCoMET")
    library(BayesCoMET)
    ```
    
2. **Others Folder**
    - contains the codes implementing the competing methods.
    - `sampler_Oracle.R` contains the code for oracle, the Bayesian oracle benchmark of CoMET, where the true random-effects covariance structure is known.
    - `fanli2012.R` contains the code for implementing the penalized quasi-likelihood method for fixed effects selection by [1].
    - `licaili2021.R` contains the code implementing the penalized quasi-likelihood estimation and inference procedures for fixed effects selection by [2], using their published supplementary code as a reference.
    - `gee_cv.m`, `gee_equicorr_predict.m`, `gee_run_DEAM.m` are the MATLAB codes using the SparseReg MATLAB library for implementing the GEE approach [3]. 


## Acknowledgement

Sreya Sarkar was supported by the National Science Foundation grant DMS-1854667 for this work.


## Citation

If you use *BayesCoMET* in your work, please cite:
Sarkar, S., Khare, K., & Srivastava, S. (2026). **CoMET: A Compressed Bayesian Mixed-Effects Model for High-Dimensional Tensors.** *arXiv.* [https://arxiv.org/pdf/2602.19236](https://arxiv.org/pdf/2602.19236)


## Other References

[1] Fan, Y., & Li, R. (2012). **Variable Selection in Linear Mixed Effects Models.** *The Annals of Statistics*, Vol. 40, No. 4, pp. 2043–2068. DOI: [10.1214/12-AOS1028](https://doi.org/10.1214/12-AOS1028).

[2] Li, S., Cai T. T., & Li, H. (2022). **Inference for High-Dimensional Linear Mixed-Effects Models: A Quasi-Likelihood Approach.** *Journal of the American Statistical Association*, 117(540), 1835–1846. DOI: [10.1080/01621459.2021.1888740](https://doi.org/10.1080/01621459.2021.1888740).

[3] Zhang X., Li, L., Zhou, H., Shen, D., et al. (2019). **Tensor generalized estimating equations for longitudinal imaging analysis.** *Statistica Sinica*, 29(4), 1977-2005. DOI: [10.5705/ss.202017.0153](https://doi.org/10.5705/ss.202017.0153).

