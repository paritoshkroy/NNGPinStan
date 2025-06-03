
# Gaussian processes and point-referenced spatial data

## Gaussian processes

A Gaussian process is a random function
![\\z(\mathbf{s}): \mathbf{s} \in \mathcal{D}\\](https://latex.codecogs.com/svg.image?%5C%7Bz%28%5Cmathbf%7Bs%7D%29%3A%20%5Cmathbf%7Bs%7D%20%5Cin%20%5Cmathcal%7BD%7D%5C%7D "\{z(\mathbf{s}): \mathbf{s} \in \mathcal{D}\}")
defined over a
![d](https://latex.codecogs.com/svg.image?d "d")–dimensional surface
(domain)
![\mathcal{D} \subset \mathbb{R}^d](https://latex.codecogs.com/svg.image?%5Cmathcal%7BD%7D%20%5Csubset%20%5Cmathbb%7BR%7D%5Ed "\mathcal{D} \subset \mathbb{R}^d"),
any finite number of which have a multivariate normal distribution.
Therefore, a Gaussian process can be completely specified by a mean
function
![\mu(\mathbf{s}) = \mathbb{E}\left\[z(\mathbf{s})\right\]](https://latex.codecogs.com/svg.image?%5Cmu%28%5Cmathbf%7Bs%7D%29%20%3D%20%5Cmathbb%7BE%7D%5Cleft%5Bz%28%5Cmathbf%7Bs%7D%29%5Cright%5D "\mu(\mathbf{s}) = \mathbb{E}\left[z(\mathbf{s})\right]")
and a covariance function
![C(\mathbf{s},\mathbf{s}') = \text{Cov}\left\[z(\mathbf{s}),z(\mathbf{s}^\prime)\right\]](https://latex.codecogs.com/svg.image?C%28%5Cmathbf%7Bs%7D%2C%5Cmathbf%7Bs%7D%27%29%20%3D%20%5Ctext%7BCov%7D%5Cleft%5Bz%28%5Cmathbf%7Bs%7D%29%2Cz%28%5Cmathbf%7Bs%7D%5E%5Cprime%29%5Cright%5D "C(\mathbf{s},\mathbf{s}') = \text{Cov}\left[z(\mathbf{s}),z(\mathbf{s}^\prime)\right]").
The covariance function is commonly specified as the product of a
variance parameter and a correlation function, a function of Euclidean
distance between locations in
![\mathcal{D}](https://latex.codecogs.com/svg.image?%5Cmathcal%7BD%7D "\mathcal{D}"),
that is,
![C(\mathbf{s},\mathbf{s}') = \sigma^2 \rho(\mathbf{s}, \mathbf{s}^\prime)](https://latex.codecogs.com/svg.image?C%28%5Cmathbf%7Bs%7D%2C%5Cmathbf%7Bs%7D%27%29%20%3D%20%5Csigma%5E2%20%5Crho%28%5Cmathbf%7Bs%7D%2C%20%5Cmathbf%7Bs%7D%5E%5Cprime%29 "C(\mathbf{s},\mathbf{s}') = \sigma^2 \rho(\mathbf{s}, \mathbf{s}^\prime)"),
where ![\sigma](https://latex.codecogs.com/svg.image?%5Csigma "\sigma")
is the marginal standard deviation and
![\rho(\mathbf{s}, \mathbf{s}^\prime)](https://latex.codecogs.com/svg.image?%5Crho%28%5Cmathbf%7Bs%7D%2C%20%5Cmathbf%7Bs%7D%5E%5Cprime%29 "\rho(\mathbf{s}, \mathbf{s}^\prime)")
is a valid correlation function that depends on the Euclidean distance
between
![\mathbf{s}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D "\mathbf{s}")
and
![\mathbf{s}^\prime](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D%5E%5Cprime "\mathbf{s}^\prime")
in
![\mathcal{D}](https://latex.codecogs.com/svg.image?%5Cmathcal%7BD%7D "\mathcal{D}").
The resulting Gaussian process is stationary and isotropic , and is also
called homogeneous Gaussian process .

A common example of a valid covariance function is the Mat'ern family of
isotropic covariance functions, which is given by

![\begin{align}
\label{Chap2:Maternfunction}
\rho(r, \nu, \ell) = \dfrac{2^{1-\nu}}{\Gamma (\nu)}\left( \sqrt{2\nu}\\ \dfrac{r}{\ell} \right)^{\nu} \\ K\_{\nu} \left(\sqrt{2\nu} \\ \dfrac{r}{\ell}\right),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Clabel%7BChap2%3AMaternfunction%7D%0A%5Crho%28r%2C%20%5Cnu%2C%20%5Cell%29%20%3D%20%5Cdfrac%7B2%5E%7B1-%5Cnu%7D%7D%7B%5CGamma%20%28%5Cnu%29%7D%5Cleft%28%20%5Csqrt%7B2%5Cnu%7D%5C%2C%20%5Cdfrac%7Br%7D%7B%5Cell%7D%20%5Cright%29%5E%7B%5Cnu%7D%20%5C%2C%20K_%7B%5Cnu%7D%20%5Cleft%28%5Csqrt%7B2%5Cnu%7D%20%5C%2C%20%5Cdfrac%7Br%7D%7B%5Cell%7D%5Cright%29%2C%0A%5Cend%7Balign%7D "\begin{align}
\label{Chap2:Maternfunction}
\rho(r, \nu, \ell) = \dfrac{2^{1-\nu}}{\Gamma (\nu)}\left( \sqrt{2\nu}\, \dfrac{r}{\ell} \right)^{\nu} \, K_{\nu} \left(\sqrt{2\nu} \, \dfrac{r}{\ell}\right),
\end{align}")

where
![r = \|\|\mathbf{s}-\mathbf{s}^\prime\|\|](https://latex.codecogs.com/svg.image?r%20%3D%20%7C%7C%5Cmathbf%7Bs%7D-%5Cmathbf%7Bs%7D%5E%5Cprime%7C%7C "r = ||\mathbf{s}-\mathbf{s}^\prime||")
is the distance between any two locations
![\mathbf{s}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D "\mathbf{s}")
and
![\mathbf{s}^\prime](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D%5E%5Cprime "\mathbf{s}^\prime")
in
![\mathcal{D}](https://latex.codecogs.com/svg.image?%5Cmathcal%7BD%7D "\mathcal{D}").
The parameter
![\nu\>0](https://latex.codecogs.com/svg.image?%5Cnu%3E0 "\nu>0")
controls the differentiability of the process, and
![\ell](https://latex.codecogs.com/svg.image?%5Cell "\ell") is called
lengthscale, which measures the distance at which the process’s
fluctuations begin to repeat. A smaller lengthscale results in a more
oscillatory function with rapid changes, capturing fine details in the
data, whereas a larger lengthscale produces a smoother function with
gradual changes, averaging out smaller fluctuations. In this sense, the
lengthscale of a Gaussian process also serves as a measure of the
smoothness or roughness of the functions it generates, impacting how it
captures patterns and variations in the data. The component
![K\_{\nu}(\cdot)](https://latex.codecogs.com/svg.image?K_%7B%5Cnu%7D%28%5Ccdot%29 "K_{\nu}(\cdot)")
in the covariance function is a modified Bessel function of the second
kind of order ![\nu](https://latex.codecogs.com/svg.image?%5Cnu "\nu").
The process becomes very rough for a small value of
![\nu](https://latex.codecogs.com/svg.image?%5Cnu "\nu") (say,
![\nu = 1/2](https://latex.codecogs.com/svg.image?%5Cnu%20%3D%201%2F2 "\nu = 1/2")),
whereas
![\nu = 3/2](https://latex.codecogs.com/svg.image?%5Cnu%20%3D%203%2F2 "\nu = 3/2")
and
![\nu = 5/2](https://latex.codecogs.com/svg.image?%5Cnu%20%3D%205%2F2 "\nu = 5/2")
are appealing cases for which the processes are once and twice
differentiable, respectively, and for
![\nu \> 7/2](https://latex.codecogs.com/svg.image?%5Cnu%20%3E%207%2F2 "\nu > 7/2")
the process is very smooth. Also, note that for
![\nu \in \\1/2, 3/2, 5/2\\](https://latex.codecogs.com/svg.image?%5Cnu%20%5Cin%20%5C%7B1%2F2%2C%203%2F2%2C%205%2F2%5C%7D "\nu \in \{1/2, 3/2, 5/2\}")
the resultant covariance function has a computationally simple form;
however, for
![\nu \notin \\1/2, 3/2, 5/2\\](https://latex.codecogs.com/svg.image?%5Cnu%20%5Cnotin%20%5C%7B1%2F2%2C%203%2F2%2C%205%2F2%5C%7D "\nu \notin \{1/2, 3/2, 5/2\}")
it is necessary to compute
![K\_{\nu}(\cdot)](https://latex.codecogs.com/svg.image?K_%7B%5Cnu%7D%28%5Ccdot%29 "K_{\nu}(\cdot)"),
which is computationally expensive. In practice, the analyst may prefer
to fix the parameter
![\nu](https://latex.codecogs.com/svg.image?%5Cnu "\nu") based on
subject knowledge.

The lengthscale parameter is closely related to the practical range in
spatial statistics, indicating the distance over which spatial
dependence between observations remains effective. Beyond this range,
observations are considered nearly independent, with correlations
diminishing to negligible levels. The practical range is crucial for
determining appropriate spatial scales for analysis and modeling,
ensuring accurate representation of spatial processes. For instance,
with
![\nu = 3/2](https://latex.codecogs.com/svg.image?%5Cnu%20%3D%203%2F2 "\nu = 3/2"),
the Mat'ern 3/2 correlation function is

![\begin{align}
C\_{3/2}(r, \ell) = (1 + \dfrac{\sqrt{3}\\ r}{\ell}) \exp\left(-\dfrac{\sqrt{3}\\r}{\ell}\right),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0AC_%7B3%2F2%7D%28r%2C%20%5Cell%29%20%3D%20%281%20%2B%20%5Cdfrac%7B%5Csqrt%7B3%7D%5C%2C%20r%7D%7B%5Cell%7D%29%20%5Cexp%5Cleft%28-%5Cdfrac%7B%5Csqrt%7B3%7D%5C%2Cr%7D%7B%5Cell%7D%5Cright%29%2C%0A%5Cend%7Balign%7D "\begin{align}
C_{3/2}(r, \ell) = (1 + \dfrac{\sqrt{3}\, r}{\ell}) \exp\left(-\dfrac{\sqrt{3}\,r}{\ell}\right),
\end{align}")

indicating how spatial correlation decreases with distance, and for
distances
![\ell/2](https://latex.codecogs.com/svg.image?%5Cell%2F2 "\ell/2"),
![\ell](https://latex.codecogs.com/svg.image?%5Cell "\ell"),
![2\ell](https://latex.codecogs.com/svg.image?2%5Cell "2\ell"),
![2.75\ell](https://latex.codecogs.com/svg.image?2.75%5Cell "2.75\ell"),
and ![4\ell](https://latex.codecogs.com/svg.image?4%5Cell "4\ell"), the
correlations are approximately 0.78, 0.48, 0.14, 0.05, and 0.008,
respectively.

By definition, for any finite set of locations
![\\\mathbf{s}\_1,\ldots,\mathbf{s}\_n\\ \in \mathcal{D}](https://latex.codecogs.com/svg.image?%5C%7B%5Cmathbf%7Bs%7D_1%2C%5Cldots%2C%5Cmathbf%7Bs%7D_n%5C%7D%20%5Cin%20%5Cmathcal%7BD%7D "\{\mathbf{s}_1,\ldots,\mathbf{s}_n\} \in \mathcal{D}")
the joint distribution of an
![n](https://latex.codecogs.com/svg.image?n "n")–dimensional vector of
possible realizations
![\mathbf{z}^\prime = (z(\mathbf{s}\_1),\ldots,z(\mathbf{s}\_n))](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D%5E%5Cprime%20%3D%20%28z%28%5Cmathbf%7Bs%7D_1%29%2C%5Cldots%2Cz%28%5Cmathbf%7Bs%7D_n%29%29 "\mathbf{z}^\prime = (z(\mathbf{s}_1),\ldots,z(\mathbf{s}_n))")
from a Gaussian process follows a multivariate normal distribution, that
is,

![\begin{align}
f(\mathbf{z} \mid \boldsymbol{\theta}) \propto \dfrac{1}{\sqrt{\|\sigma^2\mathbf{B}\|}} \exp\left\\-\dfrac{1}{2\sigma^2} (\mathbf{z} - \boldsymbol{\mu})^\prime \mathbf{B}^{-1}(\mathbf{z} - \boldsymbol{\mu})\right\\,
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0Af%28%5Cmathbf%7Bz%7D%20%5Cmid%20%5Cboldsymbol%7B%5Ctheta%7D%29%20%5Cpropto%20%5Cdfrac%7B1%7D%7B%5Csqrt%7B%7C%5Csigma%5E2%5Cmathbf%7BB%7D%7C%7D%7D%20%5Cexp%5Cleft%5C%7B-%5Cdfrac%7B1%7D%7B2%5Csigma%5E2%7D%20%28%5Cmathbf%7Bz%7D%20-%20%5Cboldsymbol%7B%5Cmu%7D%29%5E%5Cprime%20%5Cmathbf%7BB%7D%5E%7B-1%7D%28%5Cmathbf%7Bz%7D%20-%20%5Cboldsymbol%7B%5Cmu%7D%29%5Cright%5C%7D%2C%0A%5Cend%7Balign%7D "\begin{align}
f(\mathbf{z} \mid \boldsymbol{\theta}) \propto \dfrac{1}{\sqrt{|\sigma^2\mathbf{B}|}} \exp\left\{-\dfrac{1}{2\sigma^2} (\mathbf{z} - \boldsymbol{\mu})^\prime \mathbf{B}^{-1}(\mathbf{z} - \boldsymbol{\mu})\right\},
\end{align}")

where
![\boldsymbol{\mu}^\prime = (\mu(\mathbf{s}\_1),\ldots,\mu(\mathbf{s}\_n))](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7B%5Cmu%7D%5E%5Cprime%20%3D%20%28%5Cmu%28%5Cmathbf%7Bs%7D_1%29%2C%5Cldots%2C%5Cmu%28%5Cmathbf%7Bs%7D_n%29%29 "\boldsymbol{\mu}^\prime = (\mu(\mathbf{s}_1),\ldots,\mu(\mathbf{s}_n))")
is a ![n](https://latex.codecogs.com/svg.image?n "n")–dimensional vector
of means, and
![\mathbf{B}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BB%7D "\mathbf{B}")
is a ![n](https://latex.codecogs.com/svg.image?n "n")–dimensional
correlation matrix with
![(i,j)](https://latex.codecogs.com/svg.image?%28i%2Cj%29 "(i,j)")th
element is
![\rho(\mathbf{s}\_i,\mathbf{s}\_j)](https://latex.codecogs.com/svg.image?%5Crho%28%5Cmathbf%7Bs%7D_i%2C%5Cmathbf%7Bs%7D_j%29 "\rho(\mathbf{s}_i,\mathbf{s}_j)").

## Point-referenced spatial data

Point-referenced spatial data refers to observations where each data
point is associated with a precise location defined by coordinates. The
coordinates typically include latitude and longitude for global
positioning, easting and northing for local projections, or
![(x,y)](https://latex.codecogs.com/svg.image?%28x%2Cy%29 "(x,y)")–coordinates
of a surface. Analyzing point reference data aims to capture variability
and correlation in observed phenomena, predict values at unobserved
locations, and assess uncertainty. Its applications are widespread in
environmental monitoring and geophysical studies. For example, weather
stations record temperature, humidity, and air quality at some fixed
monitoring sites across a geographical area, and data analysis often
aims to obtain a predicted surface of the phenomena by estimating values
at unobserved locations.

Statistical modeling of point reference data often assumes that the
measurement variable is, theoretically, defined at locations that vary
continuously across the domain. Thus, it necessitates specifying a
random surface . A widely adopted approach involves modeling the surface
as a realization of a stochastic process. Gaussian processes provide a
practical framework for such modeling, offering a versatile tool for
representing the spatial processes that vary continuously across space.
The Gaussian processes facilitate straightforward inference and
prediction by capturing spatial correlations, interpolating data,
modeling variability, and enabling probabilistic inference. In the
following, we will explore the application of Gaussian processes in
analyzing point-referenced spatial data.

## Modeling point-referenced spatial data using GP

Within this framework, observations over a finite set of locations in a
spatial domain are assumed to be partial realizations of a spatial
Gaussian process
![\\y(\mathbf{s}): \mathbf{s} \in \mathcal{D}\\](https://latex.codecogs.com/svg.image?%5C%7By%28%5Cmathbf%7Bs%7D%29%3A%20%5Cmathbf%7Bs%7D%20%5Cin%20%5Cmathcal%7BD%7D%5C%7D "\{y(\mathbf{s}): \mathbf{s} \in \mathcal{D}\}")
that defined on spatial index
![\mathbf{s}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D "\mathbf{s}")
varying continuously throughout
![d](https://latex.codecogs.com/svg.image?d "d")–dimensional domain
![\mathcal{D} \in \mathbb{R}^{d}](https://latex.codecogs.com/svg.image?%5Cmathcal%7BD%7D%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bd%7D "\mathcal{D} \in \mathbb{R}^{d}").
Therefore, the joint distribution of measurements at any finite set of
locations is assumed to be multivariate normal, and the properties of
the multivariate normal distribution ensure closed-form marginal and
conditional distributions, leading to straightforward computation for
model fitting and prediction.

It assumes that the measurement
![y(\mathbf{s})](https://latex.codecogs.com/svg.image?y%28%5Cmathbf%7Bs%7D%29 "y(\mathbf{s})"),
at location
![\mathbf{s}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D "\mathbf{s}"),
be a generated as

![\begin{align}
y(\mathbf{s}) = \mathbf{x}(\mathbf{s})^\prime \boldsymbol{\theta} + z(\mathbf{s}) + \epsilon(\mathbf{s}),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0Ay%28%5Cmathbf%7Bs%7D%29%20%3D%20%5Cmathbf%7Bx%7D%28%5Cmathbf%7Bs%7D%29%5E%5Cprime%20%5Cboldsymbol%7B%5Ctheta%7D%20%2B%20z%28%5Cmathbf%7Bs%7D%29%20%2B%20%5Cepsilon%28%5Cmathbf%7Bs%7D%29%2C%0A%5Cend%7Balign%7D "\begin{align}
y(\mathbf{s}) = \mathbf{x}(\mathbf{s})^\prime \boldsymbol{\theta} + z(\mathbf{s}) + \epsilon(\mathbf{s}),
\end{align}")

where
![\boldsymbol{\theta}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7B%5Ctheta%7D "\boldsymbol{\theta}")
is a ![p](https://latex.codecogs.com/svg.image?p "p")–dimensional vector
of coefficient associated with
![p](https://latex.codecogs.com/svg.image?p "p")–dimensional vector of
covariates,
![\mathbf{x}(\mathbf{s})](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bx%7D%28%5Cmathbf%7Bs%7D%29 "\mathbf{x}(\mathbf{s})"),
including intercept. The component
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})")
is assumed to be a zero mean isotropic Gaussian process with marginal
variance
![\sigma^2](https://latex.codecogs.com/svg.image?%5Csigma%5E2 "\sigma^2")
and correlation function
![\rho(\cdot)](https://latex.codecogs.com/svg.image?%5Crho%28%5Ccdot%29 "\rho(\cdot)")
capturing the spatial correlation and ensures that
![y(\mathbf{s})](https://latex.codecogs.com/svg.image?y%28%5Cmathbf%7Bs%7D%29 "y(\mathbf{s})")
is defined at every location
![\mathbf{s} \in \mathcal{D}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D%20%5Cin%20%5Cmathcal%7BD%7D "\mathbf{s} \in \mathcal{D}"),
and
![\epsilon(\mathbf{s})](https://latex.codecogs.com/svg.image?%5Cepsilon%28%5Cmathbf%7Bs%7D%29 "\epsilon(\mathbf{s})")
is assumed to be an independent random measurement error which follows a
normal distribution with mean zero and variance
![\tau^2](https://latex.codecogs.com/svg.image?%5Ctau%5E2 "\tau^2").

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

![\begin{align}
\mathbb{E}\[\mathbf{y}\] = \mathbf{X}\boldsymbol{\theta} \qquad \text{and} \qquad \text{Var}\[\mathbf{y}\] = \sigma^2 \mathbf{B} + \tau^2\mathbf{I},
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cmathbb%7BE%7D%5B%5Cmathbf%7By%7D%5D%20%3D%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%20%5Cqquad%20%5Ctext%7Band%7D%20%5Cqquad%20%5Ctext%7BVar%7D%5B%5Cmathbf%7By%7D%5D%20%3D%20%5Csigma%5E2%20%5Cmathbf%7BB%7D%20%2B%20%5Ctau%5E2%5Cmathbf%7BI%7D%2C%0A%5Cend%7Balign%7D "\begin{align}
\mathbb{E}[\mathbf{y}] = \mathbf{X}\boldsymbol{\theta} \qquad \text{and} \qquad \text{Var}[\mathbf{y}] = \sigma^2 \mathbf{B} + \tau^2\mathbf{I},
\end{align}")

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

![\begin{align}
\pi(\boldsymbol{\Phi} \mid \mathbf{y}) \propto \mathcal{N}\left(\mathbf{y} \mid \mathbf{X}\boldsymbol{\theta}, \mathbf{V}\right) \\ \pi(\boldsymbol{\Phi}),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%20%5Cmid%20%5Cmathbf%7By%7D%29%20%5Cpropto%20%5Cmathcal%7BN%7D%5Cleft%28%5Cmathbf%7By%7D%20%5Cmid%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%2C%20%5Cmathbf%7BV%7D%5Cright%29%20%5C%3B%20%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%29%2C%0A%5Cend%7Balign%7D "\begin{align}
\pi(\boldsymbol{\Phi} \mid \mathbf{y}) \propto \mathcal{N}\left(\mathbf{y} \mid \mathbf{X}\boldsymbol{\theta}, \mathbf{V}\right) \; \pi(\boldsymbol{\Phi}),
\end{align}")

where
![\mathbf{V} = \sigma^2 \mathbf{B} + \tau^2\mathbf{I}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BV%7D%20%3D%20%5Csigma%5E2%20%5Cmathbf%7BB%7D%20%2B%20%5Ctau%5E2%5Cmathbf%7BI%7D "\mathbf{V} = \sigma^2 \mathbf{B} + \tau^2\mathbf{I}")
and
![\pi(\boldsymbol{\Phi})](https://latex.codecogs.com/svg.image?%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%29 "\pi(\boldsymbol{\Phi})")
denotes the prior distribution assigned to
![\boldsymbol{\Phi}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7B%5CPhi%7D "\boldsymbol{\Phi}").
In general, the distribution
![\pi(\boldsymbol{\Phi} \mid \mathbf{y})](https://latex.codecogs.com/svg.image?%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%20%5Cmid%20%5Cmathbf%7By%7D%29 "\pi(\boldsymbol{\Phi} \mid \mathbf{y})")
does not have a closed-form, and Markov chain Monte Carlo (MCMC)
sampling methods are commonly employed to approximate this distribution.
These methods are straightforward to implement using modern statistical
computing platforms such as `BUGS`, `JAGS`, `NIMBLE`, and `Stan`. MCMC
methods provide samples from the posterior distribution, which can be
used to estimate various summary statistics. Once samples from the
posterior distribution are available, predictions to unobserved
locations follow straightforwardly.

The described inference procedure utilizes a model that is marginalized
with respect to the latent GP by integrating out the latent variable
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}").
This approach allows the model to directly predict the observed
responses, which is why it is referred to as a response GP or marginal
GP.

## Response GP in Stan

    data {
      int<lower=0> n;
      int<lower=0> p;
      vector[n] y;
      matrix[n,p] X;
      array[n] vector[2] coords;
      
      vector<lower=0>[p] scale_theta;
      real<lower=0> scale_sigma;
      real<lower=0> scale_tau;
      
      real<lower = 0> a; // Shape parameters in the prior for ell
      real<lower = 0> b; // Scale parameters in the prior for ell
      
    }

    transformed data{
      
    }

    parameters {
      vector[p] theta_std;
      real<lower=0> ell;
      real<lower=0> sigma_std;
      real<lower=0> tau_std;
    }

    transformed parameters{
      vector[p] theta = scale_theta .* theta_std;
      real sigma = scale_sigma * sigma_std;
      real tau = scale_sigma * tau_std;
    }

    model {
      theta_std ~ std_normal();
      ell ~ inv_gamma(a,b);
      sigma_std ~ std_normal();
      tau_std ~ std_normal();
      vector[n] mu = X*theta;
      matrix[n,n] Sigma = gp_matern32_cov(coords, sigma, ell);
      matrix[n,n] L = cholesky_decompose(add_diag(Sigma, square(tau)));
      y ~ multi_normal_cholesky(mu, L);
    }

### Spatial Interpolation

One main interest in point-referenced spatial data analysis is obtaining
a predicted surface for the process through pointwise prediction. Let
![\\\mathbf{s}\_1^\star, \ldots, \mathbf{s}\_{n^\star}^\star\\](https://latex.codecogs.com/svg.image?%5C%7B%5Cmathbf%7Bs%7D_1%5E%5Cstar%2C%20%5Cldots%2C%20%5Cmathbf%7Bs%7D_%7Bn%5E%5Cstar%7D%5E%5Cstar%5C%7D "\{\mathbf{s}_1^\star, \ldots, \mathbf{s}_{n^\star}^\star\}")
be a set of
![n^\star](https://latex.codecogs.com/svg.image?n%5E%5Cstar "n^\star")
high resoluted grid locations covering
![\mathcal{D}](https://latex.codecogs.com/svg.image?%5Cmathcal%7BD%7D "\mathcal{D}")
and suppose that vector of
![p](https://latex.codecogs.com/svg.image?p "p") covariates values
![\mathbf{x}(\mathbf{s}^\star)](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bx%7D%28%5Cmathbf%7Bs%7D%5E%5Cstar%29 "\mathbf{x}(\mathbf{s}^\star)")
at each site is available. To obtain predicted surface along with
reporting uncertainty measures, we need posterior predicted distribution
of
![\mathbf{y}^\star = (y(\mathbf{s}\_1^\star),\ldots, y(\mathbf{s}\_{n^\star}^\star))^\prime](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D%5E%5Cstar%20%3D%20%28y%28%5Cmathbf%7Bs%7D_1%5E%5Cstar%29%2C%5Cldots%2C%20y%28%5Cmathbf%7Bs%7D_%7Bn%5E%5Cstar%7D%5E%5Cstar%29%29%5E%5Cprime "\mathbf{y}^\star = (y(\mathbf{s}_1^\star),\ldots, y(\mathbf{s}_{n^\star}^\star))^\prime"),
conditional on the observed data
![\mathbf{y}](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D "\mathbf{y}").
For this purpose, consider the joint vector
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

![\begin{align}
\mathbb{E}\[\mathbf{y}^\star \mid \mathbf{y}\] = 
\mathbf{X}^\star\boldsymbol{\theta} + \mathbf{V}^{\text{pred-to-obs}} \mathbf{V}^{-1} (\mathbf{y} - \mathbf{X}\boldsymbol{\theta})\\
\text{and} \text{Var}\[\mathbf{y}^\star \mid \mathbf{y}\] = \mathbf{V}^\star - \mathbf{V}^{\text{pred-to-obs}} \mathbf{V}^{-1} \mathbf{V}^{\text{obs-to-pred}},
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cmathbb%7BE%7D%5B%5Cmathbf%7By%7D%5E%5Cstar%20%5Cmid%20%5Cmathbf%7By%7D%5D%20%3D%20%0A%5Cmathbf%7BX%7D%5E%5Cstar%5Cboldsymbol%7B%5Ctheta%7D%20%2B%20%5Cmathbf%7BV%7D%5E%7B%5Ctext%7Bpred-to-obs%7D%7D%20%5Cmathbf%7BV%7D%5E%7B-1%7D%20%28%5Cmathbf%7By%7D%20-%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%29%5C%5C%0A%5Ctext%7Band%7D%20%5Ctext%7BVar%7D%5B%5Cmathbf%7By%7D%5E%5Cstar%20%5Cmid%20%5Cmathbf%7By%7D%5D%20%3D%20%5Cmathbf%7BV%7D%5E%5Cstar%20-%20%5Cmathbf%7BV%7D%5E%7B%5Ctext%7Bpred-to-obs%7D%7D%20%5Cmathbf%7BV%7D%5E%7B-1%7D%20%5Cmathbf%7BV%7D%5E%7B%5Ctext%7Bobs-to-pred%7D%7D%2C%0A%5Cend%7Balign%7D "\begin{align}
\mathbb{E}[\mathbf{y}^\star \mid \mathbf{y}] = 
\mathbf{X}^\star\boldsymbol{\theta} + \mathbf{V}^{\text{pred-to-obs}} \mathbf{V}^{-1} (\mathbf{y} - \mathbf{X}\boldsymbol{\theta})\\
\text{and} \text{Var}[\mathbf{y}^\star \mid \mathbf{y}] = \mathbf{V}^\star - \mathbf{V}^{\text{pred-to-obs}} \mathbf{V}^{-1} \mathbf{V}^{\text{obs-to-pred}},
\end{align}")

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

However, the above joint prediction is computationally expensive as the
conditional distribution of a multivariate normal distribution of
dimension
![(n^\star+n)](https://latex.codecogs.com/svg.image?%28n%5E%5Cstar%2Bn%29 "(n^\star+n)"),
computing conditional mean
![\mu\_{y(\mathbf{s}^\star) \mid \mathbf{y}}](https://latex.codecogs.com/svg.image?%5Cmu_%7By%28%5Cmathbf%7Bs%7D%5E%5Cstar%29%20%5Cmid%20%5Cmathbf%7By%7D%7D "\mu_{y(\mathbf{s}^\star) \mid \mathbf{y}}")
and variance
![\sigma^2\_{y(\mathbf{s}^\star) \mid \mathbf{y}}](https://latex.codecogs.com/svg.image?%5Csigma%5E2_%7By%28%5Cmathbf%7Bs%7D%5E%5Cstar%29%20%5Cmid%20%5Cmathbf%7By%7D%7D "\sigma^2_{y(\mathbf{s}^\star) \mid \mathbf{y}}")
involves expensive matrix calculations. In practice, this can be avoided
by performing predictions for each unobserved location separately. In
that case, at a generic prediction location
![\mathbf{s}^\star \in \mathcal{D}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D%5E%5Cstar%20%5Cin%20%5Cmathcal%7BD%7D "\mathbf{s}^\star \in \mathcal{D}"),
the posterior predictive distribution of
![y(\mathbf{s})](https://latex.codecogs.com/svg.image?y%28%5Cmathbf%7Bs%7D%29 "y(\mathbf{s})")
at
![\mathbf{s}^\star](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D%5E%5Cstar "\mathbf{s}^\star")
is given by

![\begin{align}
\pi(y(\mathbf{s}^\star) \mid \mathbf{y}) = \int\_{\boldsymbol{\Phi}} \mathcal{N}(y(\mathbf{s}^\star) \mid \mathbf{X}\boldsymbol{\theta}, \mathbf{V}) \pi(\boldsymbol{\Phi} \mid \mathbf{y}) \mathrm{d}\boldsymbol{\Phi}
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cpi%28y%28%5Cmathbf%7Bs%7D%5E%5Cstar%29%20%5Cmid%20%5Cmathbf%7By%7D%29%20%3D%20%5Cint_%7B%5Cboldsymbol%7B%5CPhi%7D%7D%20%5Cmathcal%7BN%7D%28y%28%5Cmathbf%7Bs%7D%5E%5Cstar%29%20%5Cmid%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%2C%20%5Cmathbf%7BV%7D%29%20%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%20%5Cmid%20%5Cmathbf%7By%7D%29%20%5Cmathrm%7Bd%7D%5Cboldsymbol%7B%5CPhi%7D%0A%5Cend%7Balign%7D "\begin{align}
\pi(y(\mathbf{s}^\star) \mid \mathbf{y}) = \int_{\boldsymbol{\Phi}} \mathcal{N}(y(\mathbf{s}^\star) \mid \mathbf{X}\boldsymbol{\theta}, \mathbf{V}) \pi(\boldsymbol{\Phi} \mid \mathbf{y}) \mathrm{d}\boldsymbol{\Phi}
\end{align}")

This procedure is known as a univariate prediction, each step of which
involves calculating the matrix inversion of order
![n](https://latex.codecogs.com/svg.image?n "n") and is still expensive
if ![n](https://latex.codecogs.com/svg.image?n "n") is large.

### Recovery of the Latent Component

One might be interested in the posterior distribution of the latent
spatial component
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})").
The inference using joint posterior distribution in equation ignores the
estimation of the latent vector
![\mathbf{z}^\prime = (z(\mathbf{s}\_1), \ldots, z(\mathbf{s}\_n))](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D%5E%5Cprime%20%3D%20%28z%28%5Cmathbf%7Bs%7D_1%29%2C%20%5Cldots%2C%20z%28%5Cmathbf%7Bs%7D_n%29%29 "\mathbf{z}^\prime = (z(\mathbf{s}_1), \ldots, z(\mathbf{s}_n))")
during model fitting. Nevertheless, we can recover the distribution of
vector
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
components via composition sampling once samples from the posterior
distribution of the parameters are available. Note that the joint
posterior distribution of
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
is

![\begin{align}
\pi(\mathbf{z} \mid \mathbf{y}) &= \int \pi(\boldsymbol{\Phi}, \mathbf{z} \mid \mathbf{y}) \\ \mathrm{d} \boldsymbol{\Phi}\\
&= \int \pi(\mathbf{z} \mid \boldsymbol{\Phi}, \mathbf{y}) \\ \pi(\boldsymbol{\Phi} \mid \mathbf{y}) \\ \mathrm{d} \boldsymbol{\Phi},
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cpi%28%5Cmathbf%7Bz%7D%20%5Cmid%20%5Cmathbf%7By%7D%29%20%26%3D%20%5Cint%20%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%2C%20%5Cmathbf%7Bz%7D%20%5Cmid%20%5Cmathbf%7By%7D%29%20%5C%3B%20%5Cmathrm%7Bd%7D%20%5Cboldsymbol%7B%5CPhi%7D%5C%5C%0A%26%3D%20%5Cint%20%5Cpi%28%5Cmathbf%7Bz%7D%20%5Cmid%20%5Cboldsymbol%7B%5CPhi%7D%2C%20%5Cmathbf%7By%7D%29%20%5C%3B%20%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%20%5Cmid%20%5Cmathbf%7By%7D%29%20%5C%3B%20%5Cmathrm%7Bd%7D%20%5Cboldsymbol%7B%5CPhi%7D%2C%0A%5Cend%7Balign%7D "\begin{align}
\pi(\mathbf{z} \mid \mathbf{y}) &= \int \pi(\boldsymbol{\Phi}, \mathbf{z} \mid \mathbf{y}) \; \mathrm{d} \boldsymbol{\Phi}\\
&= \int \pi(\mathbf{z} \mid \boldsymbol{\Phi}, \mathbf{y}) \; \pi(\boldsymbol{\Phi} \mid \mathbf{y}) \; \mathrm{d} \boldsymbol{\Phi},
\end{align}")

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
![(\mathbf{z}^{\star\prime}, \mathbf{z}^\prime)](https://latex.codecogs.com/svg.image?%28%5Cmathbf%7Bz%7D%5E%7B%5Cstar%5Cprime%7D%2C%20%5Cmathbf%7Bz%7D%5E%5Cprime%29 "(\mathbf{z}^{\star\prime}, \mathbf{z}^\prime)")
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

## Hierarchical representation of the above model

Note that the model above specification is referred to as the marginal
or response Gaussian model, and the inference and prediction procedures
are outlined based on it. However, this model can be represented
hierarchically as follows:

![\begin{align}
\mathbf{y} \mid \boldsymbol{\theta}, \mathbf{z} \sim \mathcal{N} \left(\mathbf{X}\boldsymbol{\theta} + \mathbf{z}, \tau^2\mathbf{I}\right),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cmathbf%7By%7D%20%5Cmid%20%5Cboldsymbol%7B%5Ctheta%7D%2C%20%5Cmathbf%7Bz%7D%20%5Csim%20%5Cmathcal%7BN%7D%20%5Cleft%28%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%20%2B%20%5Cmathbf%7Bz%7D%2C%20%5Ctau%5E2%5Cmathbf%7BI%7D%5Cright%29%2C%0A%5Cend%7Balign%7D "\begin{align}
\mathbf{y} \mid \boldsymbol{\theta}, \mathbf{z} \sim \mathcal{N} \left(\mathbf{X}\boldsymbol{\theta} + \mathbf{z}, \tau^2\mathbf{I}\right),
\end{align}")

![\mathbf{z} \sim \mathcal{N}(\mathbf{0}, \sigma^2\mathbf{B}).](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cmathbf%7B0%7D%2C%20%5Csigma%5E2%5Cmathbf%7BB%7D%29. "\mathbf{z} \sim \mathcal{N}(\mathbf{0}, \sigma^2\mathbf{B}).")

In practice, the response Gaussian process model is often preferred for
efficient parameter estimation, as it circumvents the need to estimate
the latent vector
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
directly. Instead, in a Bayesian analysis, once posterior samples for
the parameters are obtained, estimates for
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
can be recovered through composition sampling techniques.

## Latent GP in Stan

    data {
      int<lower=0> n;
      int<lower=0> p;
      vector[n] y;
      matrix[n,p] X;
      array[n] vector[2] coords;
      
      vector<lower=0>[p] scale_theta;
      real<lower=0> scale_sigma;
      real<lower=0> scale_tau;
      
      real<lower=0> a;
      real<lower=0> b;
    }

    transformed data{
      
    }

    parameters {
      vector[p] theta_std;
      real<lower=0> ell;
      real<lower=0> sigma_std;
      real<lower=0> tau_std;
      vector[n] noise;
    }

    transformed parameters{
      vector[p] theta = scale_theta .* theta_std;
      real sigma = scale_sigma * sigma_std;
      real tau = scale_sigma * tau_std;
      //vector[n] z = cholesky_decompose(add_diag(gp_matern32_cov(coords, sigma, ell), 1e-8)) * noise;
      vector[n] z = cholesky_decompose(add_diag(gp_exponential_cov(coords, sigma, ell), 1e-8)) * noise;
      
      //matrix[n,n] Sigma = gp_exponential_cov(coords, sigma, phi);
      //matrix[n,n] V = add_diag(V, 1e-8);
      //matrix[n,n] L = cholesky_decompose(V);
      //vector[n] z = L * noise;
    }

    model {
      theta_std ~ std_normal();
      ell ~ inv_gamma(a,b);
      sigma_std ~ std_normal();
      tau_std ~ std_normal();
      noise ~ std_normal();
      vector[n] mu = X*theta;
      
      y ~ normal(mu + z, tau);
    }

    generated quantities{
      
    }

## Computational complexity in analysing large datasets

The evaluation of the likelihood of a GP-based model requires the
inversion of an
![n \times n](https://latex.codecogs.com/svg.image?n%20%5Ctimes%20n "n \times n")
covariance matrix, which has a computational complexity of
![\mathcal{O}(n^3)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BO%7D%28n%5E3%29 "\mathcal{O}(n^3)"),
making it impractical for analyzing large datasets.

## Vecchia’s approximation and NNGP

Vecchia’s approximation is a method that efficiently approximates the
likelihood of a GP-based model where exact computation becomes
computationally infeasible. Vecchia’s approximation reduces this
computational burden by approximating the joint probability distribution
of the spatial data as a product of conditional distributions, each
depending only on a small subset of the data, typically nearby locations
called neighbors or neighboring sets.

Specifically, as the exact joint density
![\mathcal{N}(\mathbf{y} \mid \mathbf{X}\boldsymbol{\theta}, \mathbf{V})](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%28%5Cmathbf%7By%7D%20%5Cmid%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%2C%20%5Cmathbf%7BV%7D%29 "\mathcal{N}(\mathbf{y} \mid \mathbf{X}\boldsymbol{\theta}, \mathbf{V})")
can be written as

![\begin{align}
\mathcal{N}(\mathbf{y} \mid \mathbf{X}\boldsymbol{\theta}, \mathbf{V})  = f(y(\mathbf{s}\_1))\\ f(y(\mathbf{s}\_2) \mid y(\mathbf{s}\_1))\\ f(y(\mathbf{s}\_3) \mid y(\mathbf{s}\_1), y(\mathbf{s}\_2)), \cdots, f(y(\mathbf{s}\_n) \mid y(\mathbf{s}\_1), \ldots, y(\mathbf{s}\_{n-1})),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cmathcal%7BN%7D%28%5Cmathbf%7By%7D%20%5Cmid%20%5Cmathbf%7BX%7D%5Cboldsymbol%7B%5Ctheta%7D%2C%20%5Cmathbf%7BV%7D%29%20%20%3D%20f%28y%28%5Cmathbf%7Bs%7D_1%29%29%5C%3B%20f%28y%28%5Cmathbf%7Bs%7D_2%29%20%5Cmid%20y%28%5Cmathbf%7Bs%7D_1%29%29%5C%3B%20f%28y%28%5Cmathbf%7Bs%7D_3%29%20%5Cmid%20y%28%5Cmathbf%7Bs%7D_1%29%2C%20y%28%5Cmathbf%7Bs%7D_2%29%29%2C%20%5Ccdots%2C%20f%28y%28%5Cmathbf%7Bs%7D_n%29%20%5Cmid%20y%28%5Cmathbf%7Bs%7D_1%29%2C%20%5Cldots%2C%20y%28%5Cmathbf%7Bs%7D_%7Bn-1%7D%29%29%2C%0A%5Cend%7Balign%7D "\begin{align}
\mathcal{N}(\mathbf{y} \mid \mathbf{X}\boldsymbol{\theta}, \mathbf{V})  = f(y(\mathbf{s}_1))\; f(y(\mathbf{s}_2) \mid y(\mathbf{s}_1))\; f(y(\mathbf{s}_3) \mid y(\mathbf{s}_1), y(\mathbf{s}_2)), \cdots, f(y(\mathbf{s}_n) \mid y(\mathbf{s}_1), \ldots, y(\mathbf{s}_{n-1})),
\end{align}")

based on some fixed ordering of observation locations
![\\\mathbf{s}\_1,\ldots,\mathbf{s}\_n\\](https://latex.codecogs.com/svg.image?%5C%7B%5Cmathbf%7Bs%7D_1%2C%5Cldots%2C%5Cmathbf%7Bs%7D_n%5C%7D "\{\mathbf{s}_1,\ldots,\mathbf{s}_n\}")
in
![\mathcal{D}](https://latex.codecogs.com/svg.image?%5Cmathcal%7BD%7D "\mathcal{D}"),
the Vecchia method approximate the joint density as

![\begin{align}
\mathcal{N}\left(\mathbf{y} \\\|\\\mathbf{X} \boldsymbol{\theta}, \mathbf{V}\right) = \mathcal{N}\bigl(y(\mathbf{s}\_1) \\\|\\\mathbf{X}(\mathbf{s}\_1)\boldsymbol{\theta}, V(\mathbf{s}\_1)\bigr) \\ \times \\\prod\_{i=2}^{n} \pi\bigl(y(\mathbf{s}\_i) \\\|\\ \mathbf{y}\_{\mathbb{N}(\mathbf{s}\_i)}\bigr),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cmathcal%7BN%7D%5Cleft%28%5Cmathbf%7By%7D%20%5C%3B%7C%5C%3B%5Cmathbf%7BX%7D%20%5Cboldsymbol%7B%5Ctheta%7D%2C%20%5Cmathbf%7BV%7D%5Cright%29%20%3D%20%5Cmathcal%7BN%7D%5Cbigl%28y%28%5Cmathbf%7Bs%7D_1%29%20%5C%3B%7C%5C%3B%5Cmathbf%7BX%7D%28%5Cmathbf%7Bs%7D_1%29%5Cboldsymbol%7B%5Ctheta%7D%2C%20V%28%5Cmathbf%7Bs%7D_1%29%5Cbigr%29%20%5C%2C%20%5Ctimes%20%5C%2C%5Cprod_%7Bi%3D2%7D%5E%7Bn%7D%20%5Cpi%5Cbigl%28y%28%5Cmathbf%7Bs%7D_i%29%20%5C%3B%7C%5C%3B%20%5Cmathbf%7By%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%5Cbigr%29%2C%0A%5Cend%7Balign%7D "\begin{align}
\mathcal{N}\left(\mathbf{y} \;|\;\mathbf{X} \boldsymbol{\theta}, \mathbf{V}\right) = \mathcal{N}\bigl(y(\mathbf{s}_1) \;|\;\mathbf{X}(\mathbf{s}_1)\boldsymbol{\theta}, V(\mathbf{s}_1)\bigr) \, \times \,\prod_{i=2}^{n} \pi\bigl(y(\mathbf{s}_i) \;|\; \mathbf{y}_{\mathbb{N}(\mathbf{s}_i)}\bigr),
\end{align}")

where
![\mathbb{N}(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29 "\mathbb{N}(\mathbf{s}_i)")
is a set of up to ![m](https://latex.codecogs.com/svg.image?m "m")
locations among the preceding observations in the ordering, called the
nearest neighbors conditioning set, and
![\mathbf{y}\_{\mathbb{N}(\mathbf{s}\_i)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D "\mathbf{y}_{\mathbb{N}(\mathbf{s}_i)}")
is at most ![m](https://latex.codecogs.com/svg.image?m "m")–dimensional
vector formed by stacking
![y(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?y%28%5Cmathbf%7Bs%7D_i%29 "y(\mathbf{s}_i)")
for
![\mathbf{s}\_i \in \mathbb{N}(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_i%20%5Cin%20%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29 "\mathbf{s}_i \in \mathbb{N}(\mathbf{s}_i)").
Usually, the nearest neighbors set of
![\mathbf{s}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D "\mathbf{s}")
is defined based upon the Euclidean distance. With
![m \times p](https://latex.codecogs.com/svg.image?m%20%5Ctimes%20p "m \times p")
design matrix
![\mathbf{X}\_{\mathbb{N}(\mathbf{s}\_i)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BX%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D "\mathbf{X}_{\mathbb{N}(\mathbf{s}_i)}")
corresponding to the
![\mathbf{y}\_{\mathbb{N}(\mathbf{s}\_i)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D "\mathbf{y}_{\mathbb{N}(\mathbf{s}_i)}"),
the conditional distribution
![\pi(y(\mathbf{s}\_i) \\\|\\\mathbf{y}\_{\mathbb{N}(\mathbf{s}\_i)})](https://latex.codecogs.com/svg.image?%5Cpi%28y%28%5Cmathbf%7Bs%7D_i%29%20%5C%3B%7C%5C%3B%5Cmathbf%7By%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%29 "\pi(y(\mathbf{s}_i) \;|\;\mathbf{y}_{\mathbb{N}(\mathbf{s}_i)})")
is a univariate normal distribution which is obtained from the joint
distribution given by

![\begin{align}
\begin{pmatrix}
y(\mathbf{s}\_i)\\
\mathbf{y}\_{\mathbb{N}(\mathbf{s}\_i)}
\end{pmatrix} 
\sim \mathcal{N}\left(
\begin{bmatrix}
\mathbf{x}(\mathbf{s}\_i)\boldsymbol{\theta}\\
\mathbf{X}\_{\mathbb{N}(\mathbf{s}\_i)}\boldsymbol{\theta}
\end{bmatrix}, 
\begin{bmatrix}
V(\mathbf{s}\_i) & \mathbf{V}\_{\mathbf{s}\_i,\mathbb{N}(\mathbf{s}\_i)}\\
\mathbf{V}\_{\mathbf{s}\_i,\mathbb{N}(\mathbf{s}\_i)} & \mathbf{V}\_{\mathbb{N}(\mathbf{s}\_i)}
\end{bmatrix}\right),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cbegin%7Bpmatrix%7D%0Ay%28%5Cmathbf%7Bs%7D_i%29%5C%5C%0A%5Cmathbf%7By%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%0A%5Cend%7Bpmatrix%7D%20%0A%5Csim%20%5Cmathcal%7BN%7D%5Cleft%28%0A%5Cbegin%7Bbmatrix%7D%0A%5Cmathbf%7Bx%7D%28%5Cmathbf%7Bs%7D_i%29%5Cboldsymbol%7B%5Ctheta%7D%5C%5C%0A%5Cmathbf%7BX%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%5Cboldsymbol%7B%5Ctheta%7D%0A%5Cend%7Bbmatrix%7D%2C%20%0A%5Cbegin%7Bbmatrix%7D%0AV%28%5Cmathbf%7Bs%7D_i%29%20%26%20%5Cmathbf%7BV%7D_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%5C%5C%0A%5Cmathbf%7BV%7D_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%20%26%20%5Cmathbf%7BV%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%0A%5Cend%7Bbmatrix%7D%5Cright%29%2C%0A%5Cend%7Balign%7D "\begin{align}
\begin{pmatrix}
y(\mathbf{s}_i)\\
\mathbf{y}_{\mathbb{N}(\mathbf{s}_i)}
\end{pmatrix} 
\sim \mathcal{N}\left(
\begin{bmatrix}
\mathbf{x}(\mathbf{s}_i)\boldsymbol{\theta}\\
\mathbf{X}_{\mathbb{N}(\mathbf{s}_i)}\boldsymbol{\theta}
\end{bmatrix}, 
\begin{bmatrix}
V(\mathbf{s}_i) & \mathbf{V}_{\mathbf{s}_i,\mathbb{N}(\mathbf{s}_i)}\\
\mathbf{V}_{\mathbf{s}_i,\mathbb{N}(\mathbf{s}_i)} & \mathbf{V}_{\mathbb{N}(\mathbf{s}_i)}
\end{bmatrix}\right),
\end{align}")

where
![V(\mathbf{s}\_i) = \sigma^2 + \tau^2](https://latex.codecogs.com/svg.image?V%28%5Cmathbf%7Bs%7D_i%29%20%3D%20%5Csigma%5E2%20%2B%20%5Ctau%5E2 "V(\mathbf{s}_i) = \sigma^2 + \tau^2"),
![\mathbf{V}\_{\mathbb{N}(\mathbf{s}\_i)} = \sigma^2 \mathbf{B}\_{\mathbb{N}(\mathbf{s}\_i)} + \tau^2 \mathbf{I}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BV%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%20%3D%20%5Csigma%5E2%20%5Cmathbf%7BB%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%20%2B%20%5Ctau%5E2%20%5Cmathbf%7BI%7D "\mathbf{V}_{\mathbb{N}(\mathbf{s}_i)} = \sigma^2 \mathbf{B}_{\mathbb{N}(\mathbf{s}_i)} + \tau^2 \mathbf{I}")
is the covariance among
![\mathbf{y}\_{\mathbb{N}(\mathbf{s}\_i)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D "\mathbf{y}_{\mathbb{N}(\mathbf{s}_i)}"),
and
![\mathbf{V}\_{\mathbf{s}\_i, \mathbb{N}(\mathbf{s}\_i)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BV%7D_%7B%5Cmathbf%7Bs%7D_i%2C%20%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D "\mathbf{V}_{\mathbf{s}_i, \mathbb{N}(\mathbf{s}_i)}")
is the covariance between
![y(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?y%28%5Cmathbf%7Bs%7D_i%29 "y(\mathbf{s}_i)")
and
![\mathbf{y}\_{\mathbf{N}(\mathbf{s}\_i)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D_%7B%5Cmathbf%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D "\mathbf{y}_{\mathbf{N}(\mathbf{s}_i)}").
Therefore, using the conditional means

![\begin{align}
\mu\_{\mathbf{y}(\mathbf{s}\_i) \\\|\\\mathbf{y}\_{\mathbb{N}(\mathbf{s}\_i)}} = \mathbf{x}(\mathbf{s}\_i)\boldsymbol{\theta} + \mathbf{V}\_{\mathbf{s}\_i,\mathbb{N}(\mathbf{s}\_i)} \mathbf{V}\_{\mathbb{N}(\mathbf{s}\_i)}^{-1} (\mathbf{y}\_{\mathbb{N}(\mathbf{s}\_i)} - \mathbf{X}\_{\mathbb{N}(\mathbf{s}\_i)}\boldsymbol{\theta})
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cmu_%7B%5Cmathbf%7By%7D%28%5Cmathbf%7Bs%7D_i%29%20%5C%3B%7C%5C%3B%5Cmathbf%7By%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%7D%20%3D%20%5Cmathbf%7Bx%7D%28%5Cmathbf%7Bs%7D_i%29%5Cboldsymbol%7B%5Ctheta%7D%20%2B%20%5Cmathbf%7BV%7D_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%20%5Cmathbf%7BV%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%5E%7B-1%7D%20%28%5Cmathbf%7By%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%20-%20%5Cmathbf%7BX%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%5Cboldsymbol%7B%5Ctheta%7D%29%0A%5Cend%7Balign%7D "\begin{align}
\mu_{\mathbf{y}(\mathbf{s}_i) \;|\;\mathbf{y}_{\mathbb{N}(\mathbf{s}_i)}} = \mathbf{x}(\mathbf{s}_i)\boldsymbol{\theta} + \mathbf{V}_{\mathbf{s}_i,\mathbb{N}(\mathbf{s}_i)} \mathbf{V}_{\mathbb{N}(\mathbf{s}_i)}^{-1} (\mathbf{y}_{\mathbb{N}(\mathbf{s}_i)} - \mathbf{X}_{\mathbb{N}(\mathbf{s}_i)}\boldsymbol{\theta})
\end{align}")

and conditional variances

![\begin{align}
\boldsymbol{\Sigma}\_{\mathbf{y}(\mathbf{s}\_i) \\\|\\\mathbf{y}\_{\mathbb{N}(\mathbf{s}\_i)}} = V(\mathbf{s}\_i) - \mathbf{V}\_{\mathbf{s}\_i,\mathbb{N}(\mathbf{s}\_i)} \mathbf{V}\_{\mathbb{N}(\mathbf{s}\_i)}^{-1} \mathbf{V}\_{\mathbf{s}\_i,\mathbb{N}(\mathbf{s}\_i)}^\prime
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cboldsymbol%7B%5CSigma%7D_%7B%5Cmathbf%7By%7D%28%5Cmathbf%7Bs%7D_i%29%20%5C%3B%7C%5C%3B%5Cmathbf%7By%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%7D%20%3D%20V%28%5Cmathbf%7Bs%7D_i%29%20-%20%5Cmathbf%7BV%7D_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%20%5Cmathbf%7BV%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%5E%7B-1%7D%20%5Cmathbf%7BV%7D_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_i%29%7D%5E%5Cprime%0A%5Cend%7Balign%7D "\begin{align}
\boldsymbol{\Sigma}_{\mathbf{y}(\mathbf{s}_i) \;|\;\mathbf{y}_{\mathbb{N}(\mathbf{s}_i)}} = V(\mathbf{s}_i) - \mathbf{V}_{\mathbf{s}_i,\mathbb{N}(\mathbf{s}_i)} \mathbf{V}_{\mathbb{N}(\mathbf{s}_i)}^{-1} \mathbf{V}_{\mathbf{s}_i,\mathbb{N}(\mathbf{s}_i)}^\prime
\end{align}")

for
![i = 1, \ldots, n,](https://latex.codecogs.com/svg.image?i%20%3D%201%2C%20%5Cldots%2C%20n%2C "i = 1, \ldots, n,")
the joint distribution
![\mathcal{N}\left(\mathbf{y} \\\|\\\mathbf{X} \boldsymbol{\theta}, \mathbf{V}\right)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%5Cleft%28%5Cmathbf%7By%7D%20%5C%3B%7C%5C%3B%5Cmathbf%7BX%7D%20%5Cboldsymbol%7B%5Ctheta%7D%2C%20%5Cmathbf%7BV%7D%5Cright%29 "\mathcal{N}\left(\mathbf{y} \;|\;\mathbf{X} \boldsymbol{\theta}, \mathbf{V}\right)")
can be approximated by evaluating
![q](https://latex.codecogs.com/svg.image?q "q")–dimensional normal
distributions for ![n](https://latex.codecogs.com/svg.image?n "n")
times.

This approach facilitates computationally efficient operations while
minimizing storage requirements. Specifically, Vecchia’s approximation
for computing the likelihood involves calculating only the conditional
mean and variances, which requires
![\mathcal{O}(nm^3)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BO%7D%28nm%5E3%29 "\mathcal{O}(nm^3)")
floating point operations (flops) and
![\mathcal{O}(nm^2)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BO%7D%28nm%5E2%29 "\mathcal{O}(nm^2)")
storage. In contrast, computing the exact GP likelihood would require
![\mathcal{O}(n^3)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BO%7D%28n%5E3%29 "\mathcal{O}(n^3)")
flops and
![\mathcal{O}(n^2)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BO%7D%28n%5E2%29 "\mathcal{O}(n^2)")
storage.

**Spatial interpolation using Vecchia’s approximation**

Following the Vecchia’s approximation, the posterior predictive
distribution of
![y(\mathbf{s}^\star)](https://latex.codecogs.com/svg.image?y%28%5Cmathbf%7Bs%7D%5E%5Cstar%29 "y(\mathbf{s}^\star)")
at any unobserved location
![\mathbf{s}^\star \in \mathcal{D}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D%5E%5Cstar%20%5Cin%20%5Cmathcal%7BD%7D "\mathbf{s}^\star \in \mathcal{D}")
is given by

![\begin{align}
\pi(y(\mathbf{s}^\star) \mid \mathbf{y}) &= \int \mathcal{N}\left(y(\mathbf{s}^\star) \mid \mu\_{y(\mathbf{s}^\star)\mid \mathbf{y}}, \sigma^2\_{y(\mathbf{s}^\star) \mid \mathbf{y}} \right)\\ \pi(\boldsymbol{\Phi}, \mid \mathbf{y}) \\ d\boldsymbol{\Phi},
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cpi%28y%28%5Cmathbf%7Bs%7D%5E%5Cstar%29%20%5Cmid%20%5Cmathbf%7By%7D%29%20%26%3D%20%5Cint%20%5Cmathcal%7BN%7D%5Cleft%28y%28%5Cmathbf%7Bs%7D%5E%5Cstar%29%20%5Cmid%20%5Cmu_%7By%28%5Cmathbf%7Bs%7D%5E%5Cstar%29%5Cmid%20%5Cmathbf%7By%7D%7D%2C%20%5Csigma%5E2_%7By%28%5Cmathbf%7Bs%7D%5E%5Cstar%29%20%5Cmid%20%5Cmathbf%7By%7D%7D%20%5Cright%29%5C%3B%20%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%2C%20%5Cmid%20%5Cmathbf%7By%7D%29%20%5C%3B%20d%5Cboldsymbol%7B%5CPhi%7D%2C%0A%5Cend%7Balign%7D "\begin{align}
\pi(y(\mathbf{s}^\star) \mid \mathbf{y}) &= \int \mathcal{N}\left(y(\mathbf{s}^\star) \mid \mu_{y(\mathbf{s}^\star)\mid \mathbf{y}}, \sigma^2_{y(\mathbf{s}^\star) \mid \mathbf{y}} \right)\; \pi(\boldsymbol{\Phi}, \mid \mathbf{y}) \; d\boldsymbol{\Phi},
\end{align}")

which also does not have a closed form and composition sampling is used
to obtain samples from
![\pi(y(\mathbf{s}^\star) \mid \mathbf{y})](https://latex.codecogs.com/svg.image?%5Cpi%28y%28%5Cmathbf%7Bs%7D%5E%5Cstar%29%20%5Cmid%20%5Cmathbf%7By%7D%29 "\pi(y(\mathbf{s}^\star) \mid \mathbf{y})").
It requires sampling from
![\mathcal{N}\bigl(y(\mathbf{s}^\star) \mid \mu\_{y(\mathbf{s}^\star)\mid \mathbf{y}}, \sigma^2\_{y(\mathbf{s}^\star) \mid \mathbf{y}}\bigr)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%5Cbigl%28y%28%5Cmathbf%7Bs%7D%5E%5Cstar%29%20%5Cmid%20%5Cmu_%7By%28%5Cmathbf%7Bs%7D%5E%5Cstar%29%5Cmid%20%5Cmathbf%7By%7D%7D%2C%20%5Csigma%5E2_%7By%28%5Cmathbf%7Bs%7D%5E%5Cstar%29%20%5Cmid%20%5Cmathbf%7By%7D%7D%5Cbigr%29 "\mathcal{N}\bigl(y(\mathbf{s}^\star) \mid \mu_{y(\mathbf{s}^\star)\mid \mathbf{y}}, \sigma^2_{y(\mathbf{s}^\star) \mid \mathbf{y}}\bigr)")
for each posterior samples of
![\boldsymbol{\Phi}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7B%5CPhi%7D "\boldsymbol{\Phi}").

However, when one is interested on the latent process model the
Vecchia’s approximation does not work. Considering this situation, Datta
et al. ([2016](#ref-datta2016hierarchical)) developed the NNGP as a
sparse approximation of to a full GP. It generalizes the idea of Vecchia
([1988](#ref-vecchia1988estimation)) from nearest neighbor approximation
of a data likelihood to the nearest neighbor approximation of the
likelihood of realizations of the process
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})").

## Response NNGP in Stan

    functions {

    real responseNNGP_matern32_lpdf(vector y, vector mu, real sigma, real tau, real lscale, matrix site2neiDist, matrix neiDistMat, array[,] int neiID, int n, int m) {
        
        vector[n] V;
        vector[n] resid = y - mu;
        vector[n] U = resid;
        real tausq = square(tau);
        real sigmasq = square(sigma);
        
        real variance_ratio_plus_1 = tausq * inv(sigmasq) + 1; // variance ratio plus 1
        int h;
        for (i in 2:n) {
          int dim = (i < (m + 1))? (i - 1) : m;
          matrix[dim, dim] neiCorMat;
          matrix[dim, dim] neiCorChol;
          vector[dim] site2neiCor;
          vector[dim] v;
          row_vector[dim] v2;
          
          if(dim == 1){
            neiCorMat[1, 1] = variance_ratio_plus_1;
            } else {
              h = 0;
              for (j in 1:(dim - 1)){
                for (k in (j + 1):dim){
                  h = h + 1;
                  neiCorMat[j, k] = (1 + sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale)) * exp(-sqrt(3) * neiDistMat[(i - 1), h] * inv(lscale));
                  neiCorMat[k, j] = neiCorMat[j, k];
                  }
                }
                for(j in 1:dim){
                  neiCorMat[j, j] = variance_ratio_plus_1;
                }
            }

            neiCorChol = cholesky_decompose(add_diag(neiCorMat, 1e-7));
            site2neiCor = to_vector((1 + sqrt(3) * site2neiDist[(i - 1), 1: dim] * inv(lscale)) .* exp(-sqrt(3) * site2neiDist[(i - 1), 1: dim] * inv(lscale)));
            v = mdivide_left_tri_low(neiCorChol, site2neiCor);
            V[i] = variance_ratio_plus_1 - dot_self(v);
            v2 = mdivide_right_tri_low(v', neiCorChol);
            U[i] = U[i] - v2 * resid[neiID[(i - 1), 1:dim]];
            }
            V[1] = variance_ratio_plus_1;
            return - 0.5 * ( 1 / sigmasq * dot_product(U, (U ./ V)) + sum(log(V)) + n * log(sigmasq));
          }
     
    }


    data{
      int n;
      int m;
      int p;
      vector[n] y;
      matrix[n,p] X;
      
      vector<lower=0>[p] scale_theta;
      real<lower = 0> scale_tau;
      real<lower = 0> scale_sigma;
      
      real<lower = 0> a; // Shape parameters in the prior for ell
      real<lower = 0> b; // Scale parameters in the prior for ell
      
      array[n-1, m] int neiID;
      matrix[n-1, m] site2neiDist;
      matrix[n-1, (m*(m-1))%/%2] neiDistMat;
      array[n-1] int nNei;

    }

    transformed data {
    }

    parameters {
      vector[p] theta_std;
      real<lower=0> tau_std;
      real<lower=0> ell;
      real<lower=0> sigma_std;
    }

    transformed parameters{
      vector[p] theta = scale_theta .* theta_std;
      real tau = scale_tau * tau_std;
      real sigma = scale_sigma * sigma_std;
    }


    model {
      theta_std ~ std_normal();
      tau_std ~ std_normal();
      sigma_std ~ std_normal();
      ell ~ inv_gamma(a,b);
      
      vector[n] mu = X*theta;
      target += responseNNGP_matern32_lpdf(y | mu, sigma, tau, ell, site2neiDist, neiDistMat, neiID, n, m);
      
    }

    generated quantities{
      
    }

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-datta2016hierarchical" class="csl-entry">

Datta, Abhirup, Sudipto Banerjee, Andrew O Finley, and Alan E Gelfand.
2016. “Hierarchical Nearest-Neighbor Gaussian Process Models for Large
Geostatistical Datasets.” *Journal of the American Statistical
Association* 111 (514): 800–812.

</div>

<div id="ref-vecchia1988estimation" class="csl-entry">

Vecchia, Aldo V. 1988. “Estimation and Model Identification for
Continuous Spatial Processes.” *Journal of the Royal Statistical
Society: Series B (Methodological)* 50 (2): 297–312.

</div>

</div>
