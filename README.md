Nearest Neighbor GP Models in Stan
================

- [Gaussian Process](#gaussian-process)
- [NNGP approximation of a GP](#nngp-approximation-of-a-gp)
- [References](#references)

## Gaussian Process

### Inference Procedure

Let
![\mathbf{y} = (y(\mathbf{s}\_1),\ldots, y(\mathbf{s}\_n))^\prime](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D%20%3D%20%28y%28%5Cmathbf%7Bs%7D_1%29%2C%5Cldots%2C%20y%28%5Cmathbf%7Bs%7D_n%29%29%5E%5Cprime "\mathbf{y} = (y(\mathbf{s}_1),\ldots, y(\mathbf{s}_n))^\prime")
be the ![n](https://latex.codecogs.com/svg.image?n "n")–dimensional
vector of data over ![n](https://latex.codecogs.com/svg.image?n "n")
locations
![\mathbf{s}\_1,\ldots,\mathbf{s}\_n](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_1%2C%5Cldots%2C%5Cmathbf%7Bs%7D_n "\mathbf{s}_1,\ldots,\mathbf{s}_n")
in ![d](https://latex.codecogs.com/svg.image?d "d")–dimensional domain
![\mathcal{D} \in \mathbb{R}^{d}](https://latex.codecogs.com/svg.image?%5Cmathcal%7BD%7D%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bd%7D "\mathcal{D} \in \mathbb{R}^{d}").
Using marginalization with respect to
![\mathbf{z}^\prime = (z(\mathbf{s}\_1),\ldots,z(s_n))](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D%5E%5Cprime%20%3D%20%28z%28%5Cmathbf%7Bs%7D_1%29%2C%5Cldots%2Cz%28s_n%29%29 "\mathbf{z}^\prime = (z(\mathbf{s}_1),\ldots,z(s_n))"),
it can be shown that
![\mathbf{y}](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D "\mathbf{y}")
is distributed as multivariate normal with
![\mathbb{E}\[\mathbf{y}\] = \mathbf{X}\boldsymbol{\theta}](https://latex.codecogs.com/svg.image?%5Cmathbb%7BE%7D%5B%5Cmathbf%7By%7D%5D%20%3D%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D "\mathbb{E}[\mathbf{y}] = \mathbf{X}\boldsymbol{\theta}")
and
![\text{Var}\[\mathbf{y}\] = \sigma^2 \mathbf{B} + \tau^2\mathbf{I}](https://latex.codecogs.com/svg.image?%5Ctext%7BVar%7D%5B%5Cmathbf%7By%7D%5D%20%3D%20%5Csigma%5E2%20%5Cmathbf%7BB%7D%20%2B%20%5Ctau%5E2%5Cmathbf%7BI%7D "\text{Var}[\mathbf{y}] = \sigma^2 \mathbf{B} + \tau^2\mathbf{I}"),
where
![\mathbf{X}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BX%7D "\mathbf{X}")
is a
![n \times p](https://latex.codecogs.com/svg.image?n%20%5Ctimes%20p "n \times p")
design matrix based on the vector of covariates
![\mathbf{x}(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bx%7D%28%5Cmathbf%7Bs%7D_i%29 "\mathbf{x}(\mathbf{s}_i)")
and
![\mathbf{B}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BB%7D "\mathbf{B}")
is the ![n](https://latex.codecogs.com/svg.image?n "n")–dimensional
correlation matrix whose
![(i,j)](https://latex.codecogs.com/svg.image?%28i%2Cj%29 "(i,j)")th
elements is
![\mathbf{B}\_{ij} = \rho(\|\|\mathbf{s}\_i-\mathbf{s}\_j\|\|)](https://latex.codecogs.com/svg.image?%5Cmathbf%7BB%7D_%7Bij%7D%20%3D%20%5Crho%28%7C%7C%5Cmathbf%7Bs%7D_i-%5Cmathbf%7Bs%7D_j%7C%7C%29 "\mathbf{B}_{ij} = \rho(||\mathbf{s}_i-\mathbf{s}_j||)").
Under the Bayesian paradigm, model specification is complete after
assigning a prior distribution for
![\boldsymbol{\beta}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7B%5Cbeta%7D "\boldsymbol{\beta}"),
![\sigma](https://latex.codecogs.com/svg.image?%5Csigma "\sigma"),
![\ell](https://latex.codecogs.com/svg.image?%5Cell "\ell") and
![\tau](https://latex.codecogs.com/svg.image?%5Ctau "\tau"). Then,
following Bayes’ theorem, the joint posterior distribution of
![\boldsymbol{\Phi} = \\\boldsymbol{\theta}, \sigma, \ell, \tau\\](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7B%5CPhi%7D%20%3D%20%5C%7B%5Cboldsymbol%7B%5Ctheta%7D%2C%20%5Csigma%2C%20%5Cell%2C%20%5Ctau%5C%7D "\boldsymbol{\Phi} = \{\boldsymbol{\theta}, \sigma, \ell, \tau\}")
is proportional to

![\pi(\boldsymbol{\Phi} \mid \mathbf{y}) \propto \mathcal{N}\left(\mathbf{y} \mid \mathbf{X}\boldsymbol{\theta}, \mathbf{V}\right) \\ \pi(\boldsymbol{\Phi}),](https://latex.codecogs.com/svg.image?%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%20%5Cmid%20%5Cmathbf%7By%7D%29%20%5Cpropto%20%5Cmathcal%7BN%7D%5Cleft%28%5Cmathbf%7By%7D%20%5Cmid%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%2C%20%5Cmathbf%7BV%7D%5Cright%29%20%5C%3B%20%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%29%2C "\pi(\boldsymbol{\Phi} \mid \mathbf{y}) \propto \mathcal{N}\left(\mathbf{y} \mid \mathbf{X}\boldsymbol{\theta}, \mathbf{V}\right) \; \pi(\boldsymbol{\Phi}),")

where
![\mathbf{V} = \sigma^2 \mathbf{B} + \tau^2\mathbf{I}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BV%7D%20%3D%20%5Csigma%5E2%20%5Cmathbf%7BB%7D%20%2B%20%5Ctau%5E2%5Cmathbf%7BI%7D "\mathbf{V} = \sigma^2 \mathbf{B} + \tau^2\mathbf{I}")
and
![\pi(\boldsymbol{\Phi})](https://latex.codecogs.com/svg.image?%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%29 "\pi(\boldsymbol{\Phi})")
denotes the prior distribution assigned to
![\boldsymbol{\Phi}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7B%5CPhi%7D "\boldsymbol{\Phi}").
In practice, the distribution
![\pi(\boldsymbol{\Phi} \mid \mathbf{y})](https://latex.codecogs.com/svg.image?%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%20%5Cmid%20%5Cmathbf%7By%7D%29 "\pi(\boldsymbol{\Phi} \mid \mathbf{y})")
does not have a closed-form, and Markov chain Monte Carlo (MCMC)
sampling methods are commonly employed to approximate this distribution.
These methods are straightforward to implement using modern statistical
computing platforms such as , , , and . MCMC methods provide samples
from the posterior distribution, which can be used to estimate various
summary statistics. Once samples from the posterior distribution are
available, predictions to unobserved locations follow straightforwardly.

### Spatial Prediction

To predict the responses
![\mathbf{y}^{\star} = (y(\mathbf{s}\_1^\star),\ldots,y(\mathbf{s}\_{n^\star}^\star))^\prime](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D%5E%7B%5Cstar%7D%20%3D%20%28y%28%5Cmathbf%7Bs%7D_1%5E%5Cstar%29%2C%5Cldots%2Cy%28%5Cmathbf%7Bs%7D_%7Bn%5E%5Cstar%7D%5E%5Cstar%29%29%5E%5Cprime "\mathbf{y}^{\star} = (y(\mathbf{s}_1^\star),\ldots,y(\mathbf{s}_{n^\star}^\star))^\prime")
at any set of
![n^\star](https://latex.codecogs.com/svg.image?n%5E%5Cstar "n^\star")
unobserved locations
![\mathbf{s}\_1^\star,\ldots,\mathbf{s}\_{n^\star}^\star](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_1%5E%5Cstar%2C%5Cldots%2C%5Cmathbf%7Bs%7D_%7Bn%5E%5Cstar%7D%5E%5Cstar "\mathbf{s}_1^\star,\ldots,\mathbf{s}_{n^\star}^\star"),
consider the joint vector
![(\mathbf{y}^{\star\prime},\mathbf{y}^\prime)](https://latex.codecogs.com/svg.image?%28%5Cmathbf%7By%7D%5E%7B%5Cstar%5Cprime%7D%2C%5Cmathbf%7By%7D%5E%5Cprime%29 "(\mathbf{y}^{\star\prime},\mathbf{y}^\prime)"),
under a Gaussian process assumption whose distribution is
![(n^\star + n)](https://latex.codecogs.com/svg.image?%28n%5E%5Cstar%20%2B%20n%29 "(n^\star + n)")–dimensional
multivariate normal. Consequently, the conditional distribution
![\mathbf{y}^\star](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D%5E%5Cstar "\mathbf{y}^\star")
given
![\mathbf{y}](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D "\mathbf{y}")
is
![n^\star](https://latex.codecogs.com/svg.image?n%5E%5Cstar "n^\star")–dimensional
multivariate normal with conditional mean and variance, respectively,
given by

![\begin{aligned}
\mathbb{E}\[\mathbf{y}^\star \mid \mathbf{y}\] = 
\mathbf{X}^\star\boldsymbol{\theta} + \mathbf{V}^{\text{pred-to-obs}} \mathbf{V}^{-1} (\mathbf{y} - \mathbf{X}\boldsymbol{\theta})\\
\text{and} \text{Var}\[\mathbf{y}^\star \mid \mathbf{y}\] = \mathbf{V}^\star - \mathbf{V}^{\text{pred-to-obs}} \mathbf{V}^{-1} \mathbf{V}^{\text{obs-to-pred}},
\end{aligned}](https://latex.codecogs.com/svg.image?%5Cbegin%7Baligned%7D%0A%5Cmathbb%7BE%7D%5B%5Cmathbf%7By%7D%5E%5Cstar%20%5Cmid%20%5Cmathbf%7By%7D%5D%20%3D%20%0A%5Cmathbf%7BX%7D%5E%5Cstar%5Cboldsymbol%7B%5Ctheta%7D%20%2B%20%5Cmathbf%7BV%7D%5E%7B%5Ctext%7Bpred-to-obs%7D%7D%20%5Cmathbf%7BV%7D%5E%7B-1%7D%20%28%5Cmathbf%7By%7D%20-%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%29%5C%5C%0A%5Ctext%7Band%7D%20%5Ctext%7BVar%7D%5B%5Cmathbf%7By%7D%5E%5Cstar%20%5Cmid%20%5Cmathbf%7By%7D%5D%20%3D%20%5Cmathbf%7BV%7D%5E%5Cstar%20-%20%5Cmathbf%7BV%7D%5E%7B%5Ctext%7Bpred-to-obs%7D%7D%20%5Cmathbf%7BV%7D%5E%7B-1%7D%20%5Cmathbf%7BV%7D%5E%7B%5Ctext%7Bobs-to-pred%7D%7D%2C%0A%5Cend%7Baligned%7D "\begin{aligned}
\mathbb{E}[\mathbf{y}^\star \mid \mathbf{y}] = 
\mathbf{X}^\star\boldsymbol{\theta} + \mathbf{V}^{\text{pred-to-obs}} \mathbf{V}^{-1} (\mathbf{y} - \mathbf{X}\boldsymbol{\theta})\\
\text{and} \text{Var}[\mathbf{y}^\star \mid \mathbf{y}] = \mathbf{V}^\star - \mathbf{V}^{\text{pred-to-obs}} \mathbf{V}^{-1} \mathbf{V}^{\text{obs-to-pred}},
\end{aligned}")

which is used to perform the prediction, where
![\mathbf{X}^\star](https://latex.codecogs.com/svg.image?%5Cmathbf%7BX%7D%5E%5Cstar "\mathbf{X}^\star")
is the
![n^\star \times p](https://latex.codecogs.com/svg.image?n%5E%5Cstar%20%5Ctimes%20p "n^\star \times p")
design matrix of covariates at prediction locations. The covariance
matrix
![\mathbf{V}^\star](https://latex.codecogs.com/svg.image?%5Cmathbf%7BV%7D%5E%5Cstar "\mathbf{V}^\star")
is equal to
![\sigma^2 \mathbf{B}^\star + \tau^2 \mathbf{I}](https://latex.codecogs.com/svg.image?%5Csigma%5E2%20%5Cmathbf%7BB%7D%5E%5Cstar%20%2B%20%5Ctau%5E2%20%5Cmathbf%7BI%7D "\sigma^2 \mathbf{B}^\star + \tau^2 \mathbf{I}"),
where
![\mathbf{B}^\star](https://latex.codecogs.com/svg.image?%5Cmathbf%7BB%7D%5E%5Cstar "\mathbf{B}^\star")
denotes the
![n^\star](https://latex.codecogs.com/svg.image?n%5E%5Cstar "n^\star")–dimensional
spatial correlation matrix among the prediction locations. The component
![\mathbf{V}^{\text{pred-to-obs}}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BV%7D%5E%7B%5Ctext%7Bpred-to-obs%7D%7D "\mathbf{V}^{\text{pred-to-obs}}")
is equal to
![\sigma^2 \mathbf{B}^{\text{pred-to-obs}}](https://latex.codecogs.com/svg.image?%5Csigma%5E2%20%5Cmathbf%7BB%7D%5E%7B%5Ctext%7Bpred-to-obs%7D%7D "\sigma^2 \mathbf{B}^{\text{pred-to-obs}}"),
where
![\mathbf{B}^\text{pred-to-obs}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BB%7D%5E%5Ctext%7Bpred-to-obs%7D "\mathbf{B}^\text{pred-to-obs}")
denotes the
![n^\star \times n](https://latex.codecogs.com/svg.image?n%5E%5Cstar%20%5Ctimes%20n "n^\star \times n")
spatial correlation matrix between prediction and observed locations.

Note that the model above specification is referred to as the marginal
or response Gaussian model, and the inference and prediction procedures
are outlined based on it. However, this model can be represented
hierarchically as follows:

![\mathbf{y} \mid \boldsymbol{\theta}, \mathbf{z} \sim \mathcal{N} \left(\mathbf{X}\boldsymbol{\theta} + \mathbf{z}, \tau^2\mathbf{I}\right),](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D%20%5Cmid%20%5Cboldsymbol%7B%5Ctheta%7D%2C%20%5Cmathbf%7Bz%7D%20%5Csim%20%5Cmathcal%7BN%7D%20%5Cleft%28%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%20%2B%20%5Cmathbf%7Bz%7D%2C%20%5Ctau%5E2%5Cmathbf%7BI%7D%5Cright%29%2C "\mathbf{y} \mid \boldsymbol{\theta}, \mathbf{z} \sim \mathcal{N} \left(\mathbf{X}\boldsymbol{\theta} + \mathbf{z}, \tau^2\mathbf{I}\right),")

![\mathbf{z} \sim \mathcal{N}(\mathbf{0}, \sigma^2\mathbf{B}).](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cmathbf%7B0%7D%2C%20%5Csigma%5E2%5Cmathbf%7BB%7D%29. "\mathbf{z} \sim \mathcal{N}(\mathbf{0}, \sigma^2\mathbf{B}).")

In practice, the response Gaussian process model is often preferred for
efficient parameter estimation, as it circumvents the need to estimate
the latent vector
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
directly. Instead, in a Bayesian analysis, once posterior samples for
the parameters are obtained, estimates for
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
can be recovered through composition sampling techniques.

### Recovery of the Latent Component

One might be interested in the posterior distribution of the latent
spatial component
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})").
The inference using joint posterior distribution in equation~ ignores
the estimation of the latent vector
![\mathbf{z}^\prime = (z(\mathbf{s}\_1), \ldots, z(\mathbf{s}\_n))](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D%5E%5Cprime%20%3D%20%28z%28%5Cmathbf%7Bs%7D_1%29%2C%20%5Cldots%2C%20z%28%5Cmathbf%7Bs%7D_n%29%29 "\mathbf{z}^\prime = (z(\mathbf{s}_1), \ldots, z(\mathbf{s}_n))")
during model fitting. Nevertheless, we can recover the distribution of
vector
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
components via composition sampling once samples from the posterior
distribution of the parameters are available. Note that the joint
posterior distribution of
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
is

![\begin{align\*}
\pi(\mathbf{z} \mid \mathbf{y}) &= \int \pi(\boldsymbol{\Phi}, \mathbf{z} \mid \mathbf{y}) \\ \mathrm{d} \boldsymbol{\Phi}\\
&= \int \pi(\mathbf{z} \mid \boldsymbol{\Phi}, \mathbf{y}) \\ \pi(\boldsymbol{\Phi} \mid \mathbf{y}) \\ \mathrm{d} \boldsymbol{\Phi},
\end{align\*}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%2A%7D%0A%5Cpi%28%5Cmathbf%7Bz%7D%20%5Cmid%20%5Cmathbf%7By%7D%29%20%26%3D%20%5Cint%20%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%2C%20%5Cmathbf%7Bz%7D%20%5Cmid%20%5Cmathbf%7By%7D%29%20%5C%3B%20%5Cmathrm%7Bd%7D%20%5Cboldsymbol%7B%5CPhi%7D%5C%5C%0A%26%3D%20%5Cint%20%5Cpi%28%5Cmathbf%7Bz%7D%20%5Cmid%20%5Cboldsymbol%7B%5CPhi%7D%2C%20%5Cmathbf%7By%7D%29%20%5C%3B%20%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%20%5Cmid%20%5Cmathbf%7By%7D%29%20%5C%3B%20%5Cmathrm%7Bd%7D%20%5Cboldsymbol%7B%5CPhi%7D%2C%0A%5Cend%7Balign%2A%7D "\begin{align*}
\pi(\mathbf{z} \mid \mathbf{y}) &= \int \pi(\boldsymbol{\Phi}, \mathbf{z} \mid \mathbf{y}) \; \mathrm{d} \boldsymbol{\Phi}\\
&= \int \pi(\mathbf{z} \mid \boldsymbol{\Phi}, \mathbf{y}) \; \pi(\boldsymbol{\Phi} \mid \mathbf{y}) \; \mathrm{d} \boldsymbol{\Phi},
\end{align*}")

and

![\begin{align}
\pi(\mathbf{z} \mid \boldsymbol{\Phi}, \mathbf{y})
&\propto \mathcal{N}(\mathbf{z} \mid \mathbf{0}, \sigma^2 \mathbf{B}) \\ \mathcal{N}\left(\mathbf{y} \mid \mathbf{X}\boldsymbol{\theta} + \mathbf{z}, \tau^2\mathbf{I}\right)\\
&\propto \exp\left\\-\frac{1}{2\sigma^2} \mathbf{z}^\prime \mathbf{B}^{-1} \mathbf{z}\right\\ \\ \exp\left\\-\frac{1}{2\tau^2} (\mathbf{y} - \mathbf{X}\boldsymbol{\theta} -\mathbf{z})^\prime (\mathbf{y} - \mathbf{X}\boldsymbol{\theta} - \mathbf{z})\right\\ \nonumber\\
&\propto \exp\left\\-\frac{1}{2}\mathbf{z}^\prime \left(\frac{1}{\tau^2} \mathbf{I} + \frac{1}{\sigma^2}\mathbf{B}^{-1} \right) \mathbf{z} - \mathbf{z}^\prime \left(\frac{1}{\tau^2} \mathbf{I}\right) (\mathbf{y} - \mathbf{X}\boldsymbol{\theta})\right\\,
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cpi%28%5Cmathbf%7Bz%7D%20%5Cmid%20%5Cboldsymbol%7B%5CPhi%7D%2C%20%5Cmathbf%7By%7D%29%0A%26%5Cpropto%20%5Cmathcal%7BN%7D%28%5Cmathbf%7Bz%7D%20%5Cmid%20%5Cmathbf%7B0%7D%2C%20%5Csigma%5E2%20%5Cmathbf%7BB%7D%29%20%5C%3B%20%5Cmathcal%7BN%7D%5Cleft%28%5Cmathbf%7By%7D%20%5Cmid%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%20%2B%20%5Cmathbf%7Bz%7D%2C%20%5Ctau%5E2%5Cmathbf%7BI%7D%5Cright%29%5C%5C%0A%26%5Cpropto%20%5Cexp%5Cleft%5C%7B-%5Cfrac%7B1%7D%7B2%5Csigma%5E2%7D%20%5Cmathbf%7Bz%7D%5E%5Cprime%20%5Cmathbf%7BB%7D%5E%7B-1%7D%20%5Cmathbf%7Bz%7D%5Cright%5C%7D%20%5C%3B%20%5Cexp%5Cleft%5C%7B-%5Cfrac%7B1%7D%7B2%5Ctau%5E2%7D%20%28%5Cmathbf%7By%7D%20-%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%20-%5Cmathbf%7Bz%7D%29%5E%5Cprime%20%28%5Cmathbf%7By%7D%20-%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%20-%20%5Cmathbf%7Bz%7D%29%5Cright%5C%7D%20%5Cnonumber%5C%5C%0A%26%5Cpropto%20%5Cexp%5Cleft%5C%7B-%5Cfrac%7B1%7D%7B2%7D%5Cmathbf%7Bz%7D%5E%5Cprime%20%5Cleft%28%5Cfrac%7B1%7D%7B%5Ctau%5E2%7D%20%5Cmathbf%7BI%7D%20%2B%20%5Cfrac%7B1%7D%7B%5Csigma%5E2%7D%5Cmathbf%7BB%7D%5E%7B-1%7D%20%5Cright%29%20%5Cmathbf%7Bz%7D%20-%20%5Cmathbf%7Bz%7D%5E%5Cprime%20%5Cleft%28%5Cfrac%7B1%7D%7B%5Ctau%5E2%7D%20%5Cmathbf%7BI%7D%5Cright%29%20%28%5Cmathbf%7By%7D%20-%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%29%5Cright%5C%7D%2C%0A%5Cend%7Balign%7D "\begin{align}
\pi(\mathbf{z} \mid \boldsymbol{\Phi}, \mathbf{y})
&\propto \mathcal{N}(\mathbf{z} \mid \mathbf{0}, \sigma^2 \mathbf{B}) \; \mathcal{N}\left(\mathbf{y} \mid \mathbf{X}\boldsymbol{\theta} + \mathbf{z}, \tau^2\mathbf{I}\right)\\
&\propto \exp\left\{-\frac{1}{2\sigma^2} \mathbf{z}^\prime \mathbf{B}^{-1} \mathbf{z}\right\} \; \exp\left\{-\frac{1}{2\tau^2} (\mathbf{y} - \mathbf{X}\boldsymbol{\theta} -\mathbf{z})^\prime (\mathbf{y} - \mathbf{X}\boldsymbol{\theta} - \mathbf{z})\right\} \nonumber\\
&\propto \exp\left\{-\frac{1}{2}\mathbf{z}^\prime \left(\frac{1}{\tau^2} \mathbf{I} + \frac{1}{\sigma^2}\mathbf{B}^{-1} \right) \mathbf{z} - \mathbf{z}^\prime \left(\frac{1}{\tau^2} \mathbf{I}\right) (\mathbf{y} - \mathbf{X}\boldsymbol{\theta})\right\},
\end{align}")

which is the kernel of the multivariate normal distribution with mean
and covariance,

![\begin{align}
\mathbb{E}\[\mathbf{z} \mid \mathbf{y}\] &= \left(\dfrac{1}{\tau^2} \mathbf{I} + \dfrac{1}{\sigma^2}\mathbf{B}^{-1} \right)^{-1} \left(\dfrac{1}{\tau^2} \mathbf{I}\right) (\mathbf{y} - \mathbf{X}\boldsymbol{\theta}),\\
\text{and}\\ \text{Var}\[\mathbf{z} \mid \mathbf{y}\] &= \left(\dfrac{1}{\tau^2} \mathbf{I} + \dfrac{1}{\sigma^2}\mathbf{B}^{-1} \right)^{-1},
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cmathbb%7BE%7D%5B%5Cmathbf%7Bz%7D%20%5Cmid%20%5Cmathbf%7By%7D%5D%20%26%3D%20%5Cleft%28%5Cdfrac%7B1%7D%7B%5Ctau%5E2%7D%20%5Cmathbf%7BI%7D%20%2B%20%5Cdfrac%7B1%7D%7B%5Csigma%5E2%7D%5Cmathbf%7BB%7D%5E%7B-1%7D%20%5Cright%29%5E%7B-1%7D%20%5Cleft%28%5Cdfrac%7B1%7D%7B%5Ctau%5E2%7D%20%5Cmathbf%7BI%7D%5Cright%29%20%28%5Cmathbf%7By%7D%20-%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%29%2C%5C%5C%0A%5Ctext%7Band%7D%5C%3B%20%5Ctext%7BVar%7D%5B%5Cmathbf%7Bz%7D%20%5Cmid%20%5Cmathbf%7By%7D%5D%20%26%3D%20%5Cleft%28%5Cdfrac%7B1%7D%7B%5Ctau%5E2%7D%20%5Cmathbf%7BI%7D%20%2B%20%5Cdfrac%7B1%7D%7B%5Csigma%5E2%7D%5Cmathbf%7BB%7D%5E%7B-1%7D%20%5Cright%29%5E%7B-1%7D%2C%0A%5Cend%7Balign%7D "\begin{align}
\mathbb{E}[\mathbf{z} \mid \mathbf{y}] &= \left(\dfrac{1}{\tau^2} \mathbf{I} + \dfrac{1}{\sigma^2}\mathbf{B}^{-1} \right)^{-1} \left(\dfrac{1}{\tau^2} \mathbf{I}\right) (\mathbf{y} - \mathbf{X}\boldsymbol{\theta}),\\
\text{and}\; \text{Var}[\mathbf{z} \mid \mathbf{y}] &= \left(\dfrac{1}{\tau^2} \mathbf{I} + \dfrac{1}{\sigma^2}\mathbf{B}^{-1} \right)^{-1},
\end{align}")

respectively. Therefore, posterior samples for
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
can be obtained by drawing samples from
![\pi(\mathbf{z} \mid \boldsymbol{\Phi}, \mathbf{y})](https://latex.codecogs.com/svg.image?%5Cpi%28%5Cmathbf%7Bz%7D%20%5Cmid%20%5Cboldsymbol%7B%5CPhi%7D%2C%20%5Cmathbf%7By%7D%29 "\pi(\mathbf{z} \mid \boldsymbol{\Phi}, \mathbf{y})")
one-for-one for each posterior sample of
![\boldsymbol{\Phi}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7B%5CPhi%7D "\boldsymbol{\Phi}").
These are post-MCMC calculations; hence, sampling is not very expensive.
Given the posterior samples for
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
associated with observed locations and
![\boldsymbol{\Phi}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7B%5CPhi%7D "\boldsymbol{\Phi}"),
it is also possible to obtain samples of the distribution of
![n^\star](https://latex.codecogs.com/svg.image?n%5E%5Cstar "n^\star")–dimensional
vector
![\mathbf{z}^\star](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D%5E%5Cstar "\mathbf{z}^\star")
of the values of
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})")
at unobserved locations
![\mathbf{s}\_{1}^\star, \ldots, \mathbf{s}\_{n^\star}^\star](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_%7B1%7D%5E%5Cstar%2C%20%5Cldots%2C%20%5Cmathbf%7Bs%7D_%7Bn%5E%5Cstar%7D%5E%5Cstar "\mathbf{s}_{1}^\star, \ldots, \mathbf{s}_{n^\star}^\star")
via composition sampling. The procedure involves assuming joint vectors
![(\mathbf{z}^{\star\prime},\mathbf{z}^\prime)](https://latex.codecogs.com/svg.image?%28%5Cmathbf%7Bz%7D%5E%7B%5Cstar%5Cprime%7D%2C%5Cmathbf%7Bz%7D%5E%5Cprime%29 "(\mathbf{z}^{\star\prime},\mathbf{z}^\prime)")
which follows
![(n^\star + n)](https://latex.codecogs.com/svg.image?%28n%5E%5Cstar%20%2B%20n%29 "(n^\star + n)")–dimensional
multivariate normal distribution and conditional distribution of
![\mathbf{z}^\star](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D%5E%5Cstar "\mathbf{z}^\star")
given
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
is used to draw samples for
![\mathbf{z}^\star](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D%5E%5Cstar "\mathbf{z}^\star").
The conditional distribution is
![n^\star](https://latex.codecogs.com/svg.image?n%5E%5Cstar "n^\star")–dimensional
multivariate normal with mean
![\mathbf{E}\[\mathbf{z}^\star \mid \mathbf{z}\] = \mathbf{B}^{\text{pred-to-obs}} \mathbf{B}^{-1} \mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BE%7D%5B%5Cmathbf%7Bz%7D%5E%5Cstar%20%5Cmid%20%5Cmathbf%7Bz%7D%5D%20%3D%20%5Cmathbf%7BB%7D%5E%7B%5Ctext%7Bpred-to-obs%7D%7D%20%5Cmathbf%7BB%7D%5E%7B-1%7D%20%5Cmathbf%7Bz%7D "\mathbf{E}[\mathbf{z}^\star \mid \mathbf{z}] = \mathbf{B}^{\text{pred-to-obs}} \mathbf{B}^{-1} \mathbf{z}")
and variance
![\text{Var}\[\mathbf{z}^\star \mid \mathbf{z}\] = \sigma^2 (\mathbf{B}^\star - \mathbf{B}^{\text{pred-to-obs}} \mathbf{B}^{-1} \mathbf{B}^{\text{obs-to-pred}})](https://latex.codecogs.com/svg.image?%5Ctext%7BVar%7D%5B%5Cmathbf%7Bz%7D%5E%5Cstar%20%5Cmid%20%5Cmathbf%7Bz%7D%5D%20%3D%20%5Csigma%5E2%20%28%5Cmathbf%7BB%7D%5E%5Cstar%20-%20%5Cmathbf%7BB%7D%5E%7B%5Ctext%7Bpred-to-obs%7D%7D%20%5Cmathbf%7BB%7D%5E%7B-1%7D%20%5Cmathbf%7BB%7D%5E%7B%5Ctext%7Bobs-to-pred%7D%7D%29 "\text{Var}[\mathbf{z}^\star \mid \mathbf{z}] = \sigma^2 (\mathbf{B}^\star - \mathbf{B}^{\text{pred-to-obs}} \mathbf{B}^{-1} \mathbf{B}^{\text{obs-to-pred}})").

### Stan implementation of marginal model

    data {
      int<lower=0> n;
      int<lower=0> p;
      vector[n] y;
      matrix[n,p] X;
      array[n] vector[2] coords;
      
      vector<lower=0>[p] scale_beta;
      real<lower=0> scale_sigma;
      real<lower=0> scale_tau;
      
      real<lower=0> a;
      real<lower=0> b;
    }

    transformed data{
      
    }

    parameters {
      vector[p] beta_std;
      real<lower=0> phi;
      real<lower=0> sigma_std;
      real<lower=0> tau_std;
    }

    transformed parameters{
      vector[p] beta = scale_beta .* beta_std;
      real sigma = scale_sigma * sigma_std;
      real tau = scale_sigma * tau_std;
    }

    model {
      beta_std ~ std_normal();
      phi ~ inv_gamma(a,b);
      sigma_std ~ std_normal();
      tau_std ~ std_normal();
      vector[n] mu = X*beta;
      //matrix[n,n] Sigma = gp_matern32_cov(coords, sigma, phi);
      matrix[n,n] Sigma = gp_exponential_cov(coords, sigma, phi);
      matrix[n,n] L = cholesky_decompose(add_diag(Sigma, square(tau)));
      y ~ multi_normal_cholesky(mu, L);
    }

### Computational Complexity of the Inference of GP Model

## NNGP approximation of a GP

Datta et al. ([2016](#ref-datta2016hierarchical)) developed the NNGP as
a sparse approximation of to a full GP. It generalizes the idea of
Vecchia ([1988](#ref-vecchia1988estimation)) from nearest neighbor
approximation of a data likelihood to the nearest neighbor approximation
of the likelihood of realizations of the process
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})"),
where the nearest neighbors set of
![\mathbf{s}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D "\mathbf{s}")
is defined based upon the Euclidean distance ([Datta
2021](#ref-datta2021sparse)). Both the nearest neighbor approximations
builds upon the idea that the joint distribution for a random vector
![\boldsymbol{z}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7Bz%7D "\boldsymbol{z}")
can be looked upon as a directed acyclic graph (DAG) to construct sparse
models for
![\boldsymbol{z}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7Bz%7D "\boldsymbol{z}")
by limiting the size of the set of parents of each node ([Banerjee
2017](#ref-banerjee2017high); [Finley et al.
2019](#ref-finley2019efficient)).

For a GP
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})")
with mean zero and covariance function
![C](https://latex.codecogs.com/svg.image?C "C"), let
![z(\mathbf{s}\_1),\ldots,z(\mathbf{s}\_n)](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_1%29%2C%5Cldots%2Cz%28%5Cmathbf%7Bs%7D_n%29 "z(\mathbf{s}_1),\ldots,z(\mathbf{s}_n)")
denotes the ![n](https://latex.codecogs.com/svg.image?n "n")
realizations of
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})")
at
![\mathbf{s}\_1,\ldots,\mathbf{s}\_n](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_1%2C%5Cldots%2C%5Cmathbf%7Bs%7D_n "\mathbf{s}_1,\ldots,\mathbf{s}_n")
and
![\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29 "\mathcal{N}(\mathbf{s})(\mathbf{s}_i)")
denotes the set of ![m](https://latex.codecogs.com/svg.image?m "m")
nearest neighbors of them. Then the NNGP is defined as the nearest
neighbor approximation of the likelihood of
![z(\mathbf{s}\_1),\ldots,z(\mathbf{s}\_n)](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_1%29%2C%5Cldots%2Cz%28%5Cmathbf%7Bs%7D_n%29 "z(\mathbf{s}_1),\ldots,z(\mathbf{s}_n)")
is given by

![\begin{align}
\label{eq_nngp_lik}
f(z(\mathbf{s}\_1)) \prod\_{i=2}^{n} f(z(\mathbf{s}\_i) \mid \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)})
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Clabel%7Beq_nngp_lik%7D%0Af%28z%28%5Cmathbf%7Bs%7D_1%29%29%20%5Cprod_%7Bi%3D2%7D%5E%7Bn%7D%20f%28z%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%29%0A%5Cend%7Balign%7D "\begin{align}
\label{eq_nngp_lik}
f(z(\mathbf{s}_1)) \prod_{i=2}^{n} f(z(\mathbf{s}_i) \mid \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)})
\end{align}")

where
![z(\mathbf{s}\_1) \sim \mathcal{N}(0,C\_{\mathbf{s}\_1,\mathbf{s}\_1})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_1%29%20%5Csim%20%5Cmathcal%7BN%7D%280%2CC_%7B%5Cmathbf%7Bs%7D_1%2C%5Cmathbf%7Bs%7D_1%7D%29 "z(\mathbf{s}_1) \sim \mathcal{N}(0,C_{\mathbf{s}_1,\mathbf{s}_1})")
and
![\boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D "\boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}")
is the set of ![m](https://latex.codecogs.com/svg.image?m "m") nearest
neighbor of the realizations. The conditional distribution
![f(z(\mathbf{s}\_i) \mid \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)})](https://latex.codecogs.com/svg.image?f%28z%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%29 "f(z(\mathbf{s}_i) \mid \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)})")
is an univariate normal
![\mathcal{N}(z(\mathbf{s}\_i) \| \mu\_{z(\mathbf{s}\_i) \mid \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}}, \sigma^2\_{z(\mathbf{s}\_i) \mid \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}})](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%28z%28%5Cmathbf%7Bs%7D_i%29%20%7C%20%5Cmu_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%2C%20%5Csigma%5E2_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%29 "\mathcal{N}(z(\mathbf{s}_i) | \mu_{z(\mathbf{s}_i) \mid \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}}, \sigma^2_{z(\mathbf{s}_i) \mid \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}})")
and which derived from the following multivariate normal

![\begin{align}
\begin{pmatrix}
z(\mathbf{s}\_i)\\
\boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}\\
\end{pmatrix} \sim \mathcal{N}\_{m+1}\left(
\begin{bmatrix}
0\\
0
\end{bmatrix},
\begin{bmatrix}
C\_{\mathbf{s}\_i\\\mathbf{s}\_i} & \boldsymbol{C}^\top\_{\mathbf{s}\_i,\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}\\
\boldsymbol{C}\_{\mathbf{s}\_i,\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)} & \boldsymbol{C}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i),\mathcal{N}(\mathbf{s})({\mathbf{s}\_i})}
\end{bmatrix}
\right),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cbegin%7Bpmatrix%7D%0Az%28%5Cmathbf%7Bs%7D_i%29%5C%5C%0A%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%5C%5C%0A%5Cend%7Bpmatrix%7D%20%5Csim%20%5Cmathcal%7BN%7D_%7Bm%2B1%7D%5Cleft%28%0A%5Cbegin%7Bbmatrix%7D%0A0%5C%5C%0A0%0A%5Cend%7Bbmatrix%7D%2C%0A%5Cbegin%7Bbmatrix%7D%0AC_%7B%5Cmathbf%7Bs%7D_i%5C%2C%5Cmathbf%7Bs%7D_i%7D%20%26%20%5Cboldsymbol%7BC%7D%5E%5Ctop_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%5C%5C%0A%5Cboldsymbol%7BC%7D_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%20%26%20%5Cboldsymbol%7BC%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%7B%5Cmathbf%7Bs%7D_i%7D%29%7D%0A%5Cend%7Bbmatrix%7D%0A%5Cright%29%2C%0A%5Cend%7Balign%7D "\begin{align}
\begin{pmatrix}
z(\mathbf{s}_i)\\
\boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}\\
\end{pmatrix} \sim \mathcal{N}_{m+1}\left(
\begin{bmatrix}
0\\
0
\end{bmatrix},
\begin{bmatrix}
C_{\mathbf{s}_i\,\mathbf{s}_i} & \boldsymbol{C}^\top_{\mathbf{s}_i,\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}\\
\boldsymbol{C}_{\mathbf{s}_i,\mathcal{N}(\mathbf{s})(\mathbf{s}_i)} & \boldsymbol{C}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i),\mathcal{N}(\mathbf{s})({\mathbf{s}_i})}
\end{bmatrix}
\right),
\end{align}")

where

![\begin{align}
\label{eq_nngp_con_moments}
\begin{split}
\mu\_{z(\mathbf{s}\_i) \mid \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}} &= \boldsymbol{C}^\top\_{\mathbf{s}\_i,\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)} \boldsymbol{C}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i),\mathcal{N}(\mathbf{s})({\mathbf{s}\_i})}^{-1} \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}\\
\sigma^2\_{z(\mathbf{s}\_i) \mid \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}} &= C\_{\mathbf{s}\_i\\\mathbf{s}\_i} - \boldsymbol{C}^\top\_{\mathbf{s}\_i,\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)} \boldsymbol{C}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i),\mathcal{N}(\mathbf{s})({\mathbf{s}\_i})}^{-1} \boldsymbol{C}\_{\mathbf{s}\_i,\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}\\
\end{split}
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Clabel%7Beq_nngp_con_moments%7D%0A%5Cbegin%7Bsplit%7D%0A%5Cmu_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%20%26%3D%20%5Cboldsymbol%7BC%7D%5E%5Ctop_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%20%5Cboldsymbol%7BC%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%7B%5Cmathbf%7Bs%7D_i%7D%29%7D%5E%7B-1%7D%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%5C%5C%0A%5Csigma%5E2_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%20%26%3D%20C_%7B%5Cmathbf%7Bs%7D_i%5C%2C%5Cmathbf%7Bs%7D_i%7D%20-%20%5Cboldsymbol%7BC%7D%5E%5Ctop_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%20%5Cboldsymbol%7BC%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%7B%5Cmathbf%7Bs%7D_i%7D%29%7D%5E%7B-1%7D%20%5Cboldsymbol%7BC%7D_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%5C%5C%0A%5Cend%7Bsplit%7D%0A%5Cend%7Balign%7D "\begin{align}
\label{eq_nngp_con_moments}
\begin{split}
\mu_{z(\mathbf{s}_i) \mid \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}} &= \boldsymbol{C}^\top_{\mathbf{s}_i,\mathcal{N}(\mathbf{s})(\mathbf{s}_i)} \boldsymbol{C}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i),\mathcal{N}(\mathbf{s})({\mathbf{s}_i})}^{-1} \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}\\
\sigma^2_{z(\mathbf{s}_i) \mid \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}} &= C_{\mathbf{s}_i\,\mathbf{s}_i} - \boldsymbol{C}^\top_{\mathbf{s}_i,\mathcal{N}(\mathbf{s})(\mathbf{s}_i)} \boldsymbol{C}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i),\mathcal{N}(\mathbf{s})({\mathbf{s}_i})}^{-1} \boldsymbol{C}_{\mathbf{s}_i,\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}\\
\end{split}
\end{align}")

The conditional distribution
![f(z(\mathbf{s}\_i) \mid \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)})](https://latex.codecogs.com/svg.image?f%28z%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%29 "f(z(\mathbf{s}_i) \mid \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)})")
can also be viewed as the likelihood from the generative model
![z(\mathbf{s}\_i) = \boldsymbol{a}\_i \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)} + e(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_i%29%20%3D%20%5Cboldsymbol%7Ba%7D_i%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%20%2B%20e%28%5Cmathbf%7Bs%7D_i%29 "z(\mathbf{s}_i) = \boldsymbol{a}_i \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)} + e(\mathbf{s}_i)"),
where
![\boldsymbol{a}\_i = \boldsymbol{C}^\top\_{\mathbf{s}\_i,\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)} \boldsymbol{C}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i),\mathcal{N}(\mathbf{s})({\mathbf{s}\_i})}^{-1}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7Ba%7D_i%20%3D%20%5Cboldsymbol%7BC%7D%5E%5Ctop_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%20%5Cboldsymbol%7BC%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%7B%5Cmathbf%7Bs%7D_i%7D%29%7D%5E%7B-1%7D "\boldsymbol{a}_i = \boldsymbol{C}^\top_{\mathbf{s}_i,\mathcal{N}(\mathbf{s})(\mathbf{s}_i)} \boldsymbol{C}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i),\mathcal{N}(\mathbf{s})({\mathbf{s}_i})}^{-1}")
and
![e(\mathbf{s}\_i) \sim i.i.d\\ \mathcal{N}(0,d\_{i})](https://latex.codecogs.com/svg.image?e%28%5Cmathbf%7Bs%7D_i%29%20%5Csim%20i.i.d%5C%2C%20%5Cmathcal%7BN%7D%280%2Cd_%7Bi%7D%29 "e(\mathbf{s}_i) \sim i.i.d\, \mathcal{N}(0,d_{i})")
where
![d_1 = C\_{\mathbf{s}\_i,\mathbf{s}\_i}](https://latex.codecogs.com/svg.image?d_1%20%3D%20C_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathbf%7Bs%7D_i%7D "d_1 = C_{\mathbf{s}_i,\mathbf{s}_i}")
and
![d_i = \sigma^2\_{z(\mathbf{s}\_i) \mid \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}}\\ i=2,\ldots,n](https://latex.codecogs.com/svg.image?d_i%20%3D%20%5Csigma%5E2_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%5C%3B%20i%3D2%2C%5Cldots%2Cn "d_i = \sigma^2_{z(\mathbf{s}_i) \mid \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}}\; i=2,\ldots,n").

Therefore, the likelihood in can be written as the following set of
linear models,
![z(\mathbf{s}\_1) = 0 + e(\mathbf{s}\_1)](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_1%29%20%3D%200%20%2B%20e%28%5Cmathbf%7Bs%7D_1%29 "z(\mathbf{s}_1) = 0 + e(\mathbf{s}_1)")
and for
![i=2,\ldots,n](https://latex.codecogs.com/svg.image?i%3D2%2C%5Cldots%2Cn "i=2,\ldots,n"),

![\begin{aligned}
z(\mathbf{s}\_i) = a\_{i1} z(\mathbf{s}\_1) + a\_{i2} z(\mathbf{s}\_2) + \cdots + a\_{i,i-1} z(\mathbf{s}\_{i-1}) + e(\mathbf{s}\_i).
\end{aligned}](https://latex.codecogs.com/svg.image?%5Cbegin%7Baligned%7D%0Az%28%5Cmathbf%7Bs%7D_i%29%20%3D%20a_%7Bi1%7D%20z%28%5Cmathbf%7Bs%7D_1%29%20%2B%20a_%7Bi2%7D%20z%28%5Cmathbf%7Bs%7D_2%29%20%2B%20%5Ccdots%20%2B%20a_%7Bi%2Ci-1%7D%20z%28%5Cmathbf%7Bs%7D_%7Bi-1%7D%29%20%2B%20e%28%5Cmathbf%7Bs%7D_i%29.%0A%5Cend%7Baligned%7D "\begin{aligned}
z(\mathbf{s}_i) = a_{i1} z(\mathbf{s}_1) + a_{i2} z(\mathbf{s}_2) + \cdots + a_{i,i-1} z(\mathbf{s}_{i-1}) + e(\mathbf{s}_i).
\end{aligned}")

where
![a\_{ij}](https://latex.codecogs.com/svg.image?a_%7Bij%7D "a_{ij}") is
the ![k](https://latex.codecogs.com/svg.image?k "k")th element of
![\boldsymbol{a}\_i](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7Ba%7D_i "\boldsymbol{a}_i")
if
![\mathbf{s}\_j](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_j "\mathbf{s}_j")
is the ![k](https://latex.codecogs.com/svg.image?k "k")th neighbor of
![\mathbf{s}\_i](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_i "\mathbf{s}_i")
and
![a\_{ij} =0](https://latex.codecogs.com/svg.image?a_%7Bij%7D%20%3D0 "a_{ij} =0")
if
![\mathbf{s}\_j](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_j "\mathbf{s}_j")
is not a neighbor of
![\mathbf{s}\_i](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_i "\mathbf{s}_i").
In matrix form, it can be written as
![\boldsymbol{z} = \boldsymbol{A} \boldsymbol{z} + \boldsymbol{e}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7Bz%7D%20%3D%20%5Cboldsymbol%7BA%7D%20%5Cboldsymbol%7Bz%7D%20%2B%20%5Cboldsymbol%7Be%7D "\boldsymbol{z} = \boldsymbol{A} \boldsymbol{z} + \boldsymbol{e}"),
where
![\boldsymbol{A}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7BA%7D "\boldsymbol{A}")
is
![n\times n](https://latex.codecogs.com/svg.image?n%5Ctimes%20n "n\times n")
strictly lower-triangular and
![\boldsymbol{e} \sim \mathcal{N}\left(\mathbf{0}, \boldsymbol{D}\right)](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7Be%7D%20%5Csim%20%5Cmathcal%7BN%7D%5Cleft%28%5Cmathbf%7B0%7D%2C%20%5Cboldsymbol%7BD%7D%5Cright%29 "\boldsymbol{e} \sim \mathcal{N}\left(\mathbf{0}, \boldsymbol{D}\right)")
and
![\boldsymbol{D}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7BD%7D "\boldsymbol{D}")
is a diagonal with elements
![d\_{i}](https://latex.codecogs.com/svg.image?d_%7Bi%7D "d_{i}") for
![i=2,\ldots,n](https://latex.codecogs.com/svg.image?i%3D2%2C%5Cldots%2Cn "i=2,\ldots,n").
Therefore,
![\mathbf{I}-\boldsymbol{A}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BI%7D-%5Cboldsymbol%7BA%7D "\mathbf{I}-\boldsymbol{A}")
is a nonsingular matrix and
![\boldsymbol{z} \sim \mathcal{N}\left(\mathbf{0},(\mathbf{I}-\boldsymbol{A})^{-1}\boldsymbol{D}(\mathbf{I}-\boldsymbol{A})^{-\top}\right)](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7Bz%7D%20%5Csim%20%5Cmathcal%7BN%7D%5Cleft%28%5Cmathbf%7B0%7D%2C%28%5Cmathbf%7BI%7D-%5Cboldsymbol%7BA%7D%29%5E%7B-1%7D%5Cboldsymbol%7BD%7D%28%5Cmathbf%7BI%7D-%5Cboldsymbol%7BA%7D%29%5E%7B-%5Ctop%7D%5Cright%29 "\boldsymbol{z} \sim \mathcal{N}\left(\mathbf{0},(\mathbf{I}-\boldsymbol{A})^{-1}\boldsymbol{D}(\mathbf{I}-\boldsymbol{A})^{-\top}\right)")
corresponds to a generative multivariate Gaussian model
![\boldsymbol{z} \sim \mathcal{N}\left(\mathbf{0},\widetilde{\boldsymbol{C}}\right)](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7Bz%7D%20%5Csim%20%5Cmathcal%7BN%7D%5Cleft%28%5Cmathbf%7B0%7D%2C%5Cwidetilde%7B%5Cboldsymbol%7BC%7D%7D%5Cright%29 "\boldsymbol{z} \sim \mathcal{N}\left(\mathbf{0},\widetilde{\boldsymbol{C}}\right)")
for the process realizations, where
![\widetilde{\boldsymbol{C}} = (\mathbf{I}-\boldsymbol{A})^{-1}\boldsymbol{D}(\mathbf{I}-\boldsymbol{A})^{-\top}](https://latex.codecogs.com/svg.image?%5Cwidetilde%7B%5Cboldsymbol%7BC%7D%7D%20%3D%20%28%5Cmathbf%7BI%7D-%5Cboldsymbol%7BA%7D%29%5E%7B-1%7D%5Cboldsymbol%7BD%7D%28%5Cmathbf%7BI%7D-%5Cboldsymbol%7BA%7D%29%5E%7B-%5Ctop%7D "\widetilde{\boldsymbol{C}} = (\mathbf{I}-\boldsymbol{A})^{-1}\boldsymbol{D}(\mathbf{I}-\boldsymbol{A})^{-\top}")
([Datta 2021](#ref-datta2021sparse)). Datta et al.
([2016](#ref-datta2016hierarchical)) shows that this is a valid
approximation of the likelihood corresponding to the GP model
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})")
with zero mean and covariance function
![C](https://latex.codecogs.com/svg.image?C "C"). However, specification
of a valid GP over the entire domain was completed by defining
prediction distribution of
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})")
at new locations conditional on
![\boldsymbol{z}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7Bz%7D "\boldsymbol{z}").

For the prediction of the latent process
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})")
at a set of ![n_0](https://latex.codecogs.com/svg.image?n_0 "n_0") new
locations
![\\\mathbf{s}\_{01},\ldots,\mathbf{s}\_{0n_0}\\ \notin \\\mathbf{s}\_1,\ldots,\mathbf{s}\_n\\](https://latex.codecogs.com/svg.image?%5C%7B%5Cmathbf%7Bs%7D_%7B01%7D%2C%5Cldots%2C%5Cmathbf%7Bs%7D_%7B0n_0%7D%5C%7D%20%5Cnotin%20%5C%7B%5Cmathbf%7Bs%7D_1%2C%5Cldots%2C%5Cmathbf%7Bs%7D_n%5C%7D "\{\mathbf{s}_{01},\ldots,\mathbf{s}_{0n_0}\} \notin \{\mathbf{s}_1,\ldots,\mathbf{s}_n\}"),
Datta et al. ([2016](#ref-datta2016hierarchical)) specified the
conditional distribution of
![z(\mathbf{s}\_{0i}) \| \boldsymbol{z}](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%20%7C%20%5Cboldsymbol%7Bz%7D "z(\mathbf{s}_{0i}) | \boldsymbol{z}")
independently, which is given by

![\begin{aligned}
z(\mathbf{s}\_{0i}) \| \boldsymbol{z} \sim \mathcal{N}\left(\boldsymbol{C}^\top\_{\mathbf{s}\_{0i},\mathcal{N}(\mathbf{s})(\mathbf{s}\_{0i})} \boldsymbol{C}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_{0i}),\mathcal{N}(\mathbf{s})(\mathbf{s}\_{0i})}^{-1} \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_{0i})}, C\_{\mathbf{s}\_{0i}\\\mathbf{s}\_{0i}} - \boldsymbol{C}^\top\_{\mathbf{s}\_{0i},\mathcal{N}(\mathbf{s})(\mathbf{s}\_{0i})} \boldsymbol{C}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_{0i}),\mathcal{N}(\mathbf{s})({\mathbf{s}\_{0i}})}^{-1} \boldsymbol{C}\_{\mathbf{s}\_{0i},\mathcal{N}(\mathbf{s})(\mathbf{s}\_{0i})}\right),
\end{aligned}](https://latex.codecogs.com/svg.image?%5Cbegin%7Baligned%7D%0Az%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%20%7C%20%5Cboldsymbol%7Bz%7D%20%5Csim%20%5Cmathcal%7BN%7D%5Cleft%28%5Cboldsymbol%7BC%7D%5E%5Ctop_%7B%5Cmathbf%7Bs%7D_%7B0i%7D%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%7D%20%5Cboldsymbol%7BC%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%7D%5E%7B-1%7D%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%7D%2C%20C_%7B%5Cmathbf%7Bs%7D_%7B0i%7D%5C%2C%5Cmathbf%7Bs%7D_%7B0i%7D%7D%20-%20%5Cboldsymbol%7BC%7D%5E%5Ctop_%7B%5Cmathbf%7Bs%7D_%7B0i%7D%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%7D%20%5Cboldsymbol%7BC%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%7B%5Cmathbf%7Bs%7D_%7B0i%7D%7D%29%7D%5E%7B-1%7D%20%5Cboldsymbol%7BC%7D_%7B%5Cmathbf%7Bs%7D_%7B0i%7D%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%7D%5Cright%29%2C%0A%5Cend%7Baligned%7D "\begin{aligned}
z(\mathbf{s}_{0i}) | \boldsymbol{z} \sim \mathcal{N}\left(\boldsymbol{C}^\top_{\mathbf{s}_{0i},\mathcal{N}(\mathbf{s})(\mathbf{s}_{0i})} \boldsymbol{C}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_{0i}),\mathcal{N}(\mathbf{s})(\mathbf{s}_{0i})}^{-1} \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_{0i})}, C_{\mathbf{s}_{0i}\,\mathbf{s}_{0i}} - \boldsymbol{C}^\top_{\mathbf{s}_{0i},\mathcal{N}(\mathbf{s})(\mathbf{s}_{0i})} \boldsymbol{C}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_{0i}),\mathcal{N}(\mathbf{s})({\mathbf{s}_{0i}})}^{-1} \boldsymbol{C}_{\mathbf{s}_{0i},\mathcal{N}(\mathbf{s})(\mathbf{s}_{0i})}\right),
\end{aligned}")

where the
![\mathcal{N}(\mathbf{s})(\mathbf{s}\_{0i})](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29 "\mathcal{N}(\mathbf{s})(\mathbf{s}_{0i})")
is the ![m](https://latex.codecogs.com/svg.image?m "m") nearest
neighbors of
![\mathbf{s}\_{0i}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_%7B0i%7D "\mathbf{s}_{0i}")
in
![\\\mathbf{s}\_1,\ldots,\mathbf{s}\_n\\](https://latex.codecogs.com/svg.image?%5C%7B%5Cmathbf%7Bs%7D_1%2C%5Cldots%2C%5Cmathbf%7Bs%7D_n%5C%7D "\{\mathbf{s}_1,\ldots,\mathbf{s}_n\}")
instead of
![\\\mathbf{s}\_{01},\ldots,\mathbf{s}\_{0n_0}\\](https://latex.codecogs.com/svg.image?%5C%7B%5Cmathbf%7Bs%7D_%7B01%7D%2C%5Cldots%2C%5Cmathbf%7Bs%7D_%7B0n_0%7D%5C%7D "\{\mathbf{s}_{01},\ldots,\mathbf{s}_{0n_0}\}").
One of the benefits of the NNGP approximation is computational as the
matrix
![\boldsymbol{A}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7BA%7D "\boldsymbol{A}")
is sparse. Computing the likelihood in involves only computing the
conditional mean and variances in and requiring
![\mathcal{O}(nm^3)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BO%7D%28nm%5E3%29 "\mathcal{O}(nm^3)")
flops and
![\mathcal{O}(nm^2)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BO%7D%28nm%5E2%29 "\mathcal{O}(nm^2)")
storage as opposed to
![\mathcal{O}(n^3)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BO%7D%28n%5E3%29 "\mathcal{O}(n^3)")
flops and
![\mathcal{O}(n^2)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BO%7D%28n%5E2%29 "\mathcal{O}(n^2)")
storage for computing the full GP likelihood. However, Bayesian
hierarchical modeling of a large spatial data using the NNGP is still
challenging. In particular, when the latent process cannot be
marginalized in the models and treated as parameters that must be
sampled. MCMC sampler becomes inefficient due to high-dimensional
parameter space.

To a similar problem Wang and Furrer ([2019](#ref-wang2019efficient))
showed that non-centered parameterization to both location and scale
parameters of the latent NNGP improves the MCMC efficiency by
implementing it in Stan.

This implies that the model
![\mathcal{N}\left(z(\mathbf{s}\_i) \| \mu\_{z(\mathbf{s}\_i) \mid \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}}, \sigma^2\_{z(\mathbf{s}\_i) \mid \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}}\right)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%5Cleft%28z%28%5Cmathbf%7Bs%7D_i%29%20%7C%20%5Cmu_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%2C%20%5Csigma%5E2_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%5Cright%29 "\mathcal{N}\left(z(\mathbf{s}_i) | \mu_{z(\mathbf{s}_i) \mid \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}}, \sigma^2_{z(\mathbf{s}_i) \mid \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}}\right)")
in is parameterized as
![z(\mathbf{s}\_i) \|  z\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)} = \mu\_{z(\mathbf{s}\_i) \mid \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}} + \sigma\_{z(\mathbf{s}\_i) \mid \boldsymbol{z}\_{\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}} v(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_i%29%20%7C%20%20z_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%20%3D%20%5Cmu_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%20%2B%20%5Csigma_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cboldsymbol%7Bz%7D_%7B%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%20v%28%5Cmathbf%7Bs%7D_i%29 "z(\mathbf{s}_i) |  z_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)} = \mu_{z(\mathbf{s}_i) \mid \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}} + \sigma_{z(\mathbf{s}_i) \mid \boldsymbol{z}_{\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}} v(\mathbf{s}_i)")
where
![v(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?v%28%5Cmathbf%7Bs%7D_i%29 "v(\mathbf{s}_i)")
is independent Gaussian noise with mean zero and variance one.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-banerjee2017high" class="csl-entry">

Banerjee, Sudipto. 2017. “High-Dimensional Bayesian Geostatistics.”
*Bayesian Analysis* 12 (2): 583–614.
<https://doi.org/10.1214/17-BA1056R>.

</div>

<div id="ref-datta2021sparse" class="csl-entry">

Datta, Abhirup. 2021. “Sparse Cholesky Matrices in Spatial Statistics.”
*arXiv Preprint arXiv:2102.13299*.

</div>

<div id="ref-datta2016hierarchical" class="csl-entry">

Datta, Abhirup, Sudipto Banerjee, Andrew O Finley, and Alan E Gelfand.
2016. “Hierarchical Nearest-Neighbor Gaussian Process Models for Large
Geostatistical Datasets.” *Journal of the American Statistical
Association* 111 (514): 800–812.

</div>

<div id="ref-finley2019efficient" class="csl-entry">

Finley, Andrew O, Abhirup Datta, Bruce D Cook, Douglas C Morton, Hans E
Andersen, and Sudipto Banerjee. 2019. “Efficient Algorithms for Bayesian
Nearest Neighbor Gaussian Processes.” *Journal of Computational and
Graphical Statistics* 28 (2): 401–14.

</div>

<div id="ref-vecchia1988estimation" class="csl-entry">

Vecchia, Aldo V. 1988. “Estimation and Model Identification for
Continuous Spatial Processes.” *Journal of the Royal Statistical
Society: Series B (Methodological)* 50 (2): 297–312.

</div>

<div id="ref-wang2019efficient" class="csl-entry">

Wang, Craig, and Reinhard Furrer. 2019. “Efficient Inference of
Generalized Spatial Fusion Models with Flexible Specification.” *Stat* 8
(1): e216.

</div>

</div>
