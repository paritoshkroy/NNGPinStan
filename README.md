Nearest Neighbor GP Models in Stan
================

### Inference Procedure

Let $\mathbf{y} = (y(\mathbf{s}_1),\ldots, y(\mathbf{s}_n))^\prime$ be
the $n$–dimensional vector of data over $n$ locations
$\mathbf{s}_1,\ldots,\mathbf{s}_n$ in $d$–dimensional domain
$\mathcal{D} \in \mathbb{R}^{d}$. Using marginalization with respect to
$\mathbf{z}^\prime = (z(\mathbf{s}_1),\ldots,z(s_n))$, it can be shown
that $\mathbf{y}$ is distributed as multivariate normal with
$\mathbb{E}[\mathbf{y}] = \mathbf{X}\boldsymbol{\theta}$ and
$\text{Var}[\mathbf{y}] = \sigma^2 \mathbf{B} + \tau^2\mathbf{I}$, where
$\mathbf{X}$ is a $n \times p$ design matrix based on the vector of
covariates $\mathbf{x}(\mathbf{s}_i)$ and $\mathbf{B}$ is the
$n$–dimensional correlation matrix whose $(i,j)$th elements is
$\mathbf{B}_{ij} = \rho(||\mathbf{s}_i-\mathbf{s}_j||)$. Under the
Bayesian paradigm, model specification is complete after assigning a
prior distribution for $\boldsymbol{\beta}$, $\sigma$, $\ell$ and
$\tau$. Then, following Bayes’ theorem, the joint posterior distribution
of $\boldsymbol{\Phi} = \{\boldsymbol{\theta}, \sigma, \ell, \tau\}$ is
proportional to where
$\mathbf{V} = \sigma^2 \mathbf{B} + \tau^2\mathbf{I}$ and
$\pi(\boldsymbol{\Phi})$ denotes the prior distribution assigned to
$\boldsymbol{\Phi}$. In practice, the distribution
$\pi(\boldsymbol{\Phi} \mid \mathbf{y})$ does not have a closed-form,
and Markov chain Monte Carlo (MCMC) sampling methods are commonly
employed to approximate this distribution. These methods are
straightforward to implement using modern statistical computing
platforms such as , , , and . MCMC methods provide samples from the
posterior distribution, which can be used to estimate various summary
statistics. Once samples from the posterior distribution are available,
predictions to unobserved locations follow straightforwardly.

### Spatial Prediction

To predict the responses
$\mathbf{y}^{\star} = (y(\mathbf{s}_1^\star),\ldots,y(\mathbf{s}_{n^\star}^\star))^\prime$
at any set of $n^\star$ unobserved locations
$\mathbf{s}_1^\star,\ldots,\mathbf{s}_{n^\star}^\star$, consider the
joint vector $(\mathbf{y}^{\star\prime},\mathbf{y}^\prime)$, under a
Gaussian process assumption whose distribution is
$(n^\star + n)$–dimensional multivariate normal. Consequently, the
conditional distribution $\mathbf{y}^\star$ given $\mathbf{y}$ is
$n^\star$–dimensional multivariate normal with conditional mean and
variance, respectively, given by which is used to perform the
prediction, where $\mathbf{X}^\star$ is the $n^\star \times p$ design
matrix of covariates at prediction locations. The covariance matrix
$\mathbf{V}^\star$ is equal to
$\sigma^2 \mathbf{B}^\star + \tau^2 \mathbf{I}$, where
$\mathbf{B}^\star$ denotes the $n^\star$–dimensional spatial correlation
matrix among the prediction locations. The component
$\mathbf{V}^{\text{pred-to-obs}}$ is equal to
$\sigma^2 \mathbf{B}^{\text{pred-to-obs}}$, where
$\mathbf{B}^\text{pred-to-obs}$ denotes the $n^\star \times n$ spatial
correlation matrix between prediction and observed locations.

Note that the model above specification is referred to as the marginal
or response Gaussian model, and the inference and prediction procedures
are outlined based on it. However, this model can be represented
hierarchically as follows: In practice, the response Gaussian process
model is often preferred for efficient parameter estimation, as it
circumvents the need to estimate the latent vector $\mathbf{z}$
directly. Instead, in a Bayesian analysis, once posterior samples for
the parameters are obtained, estimates for $\mathbf{z}$ can be recovered
through composition sampling techniques.

### Recovery of the Latent Component

One might be interested in the posterior distribution of the latent
spatial component $z(\boldsymbol{s})$. The inference using joint
posterior distribution in equation~ ignores the estimation of the latent
vector
$\mathbf{z}^\prime = (z(\boldsymbol{s}_1), \ldots, z(\boldsymbol{s}_n))$
during model fitting. Nevertheless, we can recover the distribution of
vector $\mathbf{z}$ components via composition sampling once samples
from the posterior distribution of the parameters are available. Note
that the joint posterior distribution of $\mathbf{z}$ is and which is
the kernel of the multivariate normal distribution with mean and
covariance, respectively. Therefore, posterior samples for $\mathbf{z}$
can be obtained by drawing samples from
$\pi(\mathbf{z} \mid \boldsymbol{\Phi}, \mathbf{y})$ one-for-one for
each posterior sample of $\boldsymbol{\Phi}$. These are post-MCMC
calculations; hence, sampling is not very expensive. Given the posterior
samples for $\mathbf{z}$ associated with observed locations and
$\boldsymbol{\Phi}$, it is also possible to obtain samples of the
distribution of $n^\star$–dimensional vector $\mathbf{z}^\star$ of the
values of $z(\mathbf{s})$ at unobserved locations
$\mathbf{s}_{1}^\star, \ldots, \mathbf{s}_{n^\star}^\star$ via
composition sampling. The procedure involves assuming joint vectors
$(\mathbf{z}^{\star\prime},\mathbf{z}^\prime)$ which follows
$(n^\star + n)$–dimensional multivariate normal distribution and
conditional distribution of $\mathbf{z}^\star$ given $\mathbf{z}$ is
used to draw samples for $\mathbf{z}^\star$. The conditional
distribution is $n^\star$–dimensional multivariate normal with mean
$\mathbf{E}[\mathbf{z}^\star \mid \mathbf{z}] = \mathbf{B}^{\text{pred-to-obs}} \mathbf{B}^{-1} \mathbf{z}$
and variance
$\text{Var}[\mathbf{z}^\star \mid \mathbf{z}] = \sigma^2 (\mathbf{B}^\star - \mathbf{B}^{\text{pred-to-obs}} \mathbf{B}^{-1} \mathbf{B}^{\text{obs-to-pred}})$.

## Stan implementation of marginal model
