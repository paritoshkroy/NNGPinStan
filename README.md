
# 1 Gaussian Processes and Geostatistical Data

## 1.1 Gaussian Processes

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
C\_{3/2}(r, \ell) = (1 + \sqrt{3}\\ r/\ell) \exp\left(-\sqrt{3}\\r/\ell\right),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0AC_%7B3%2F2%7D%28r%2C%20%5Cell%29%20%3D%20%281%20%2B%20%5Csqrt%7B3%7D%5C%2C%20r%2F%5Cell%29%20%5Cexp%5Cleft%28-%5Csqrt%7B3%7D%5C%2Cr%2F%5Cell%5Cright%29%2C%0A%5Cend%7Balign%7D "\begin{align}
C_{3/2}(r, \ell) = (1 + \sqrt{3}\, r/\ell) \exp\left(-\sqrt{3}\,r/\ell\right),
\end{align}")

indicating how spatial correlation decreases with distance. For
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
![{i,j}](https://latex.codecogs.com/svg.image?%7Bi%2Cj%7D "{i,j}")th
element is
![\rho(\mathbf{s}\_i,\mathbf{s}\_j)](https://latex.codecogs.com/svg.image?%5Crho%28%5Cmathbf%7Bs%7D_i%2C%5Cmathbf%7Bs%7D_j%29 "\rho(\mathbf{s}_i,\mathbf{s}_j)").

## 1.2 Geostatistical Data

Geostatistical data refers to information explicitly linked with
locations across any surface or geographical area. This thesis delves
into a specific data type where each observation or data point is
associated with a precise location defined by coordinates, known as
point reference spatial data or geostatistical data. The coordinates
typically include latitude and longitude for global positioning, easting
and northing for local projections, or
![(x,y)](https://latex.codecogs.com/svg.image?%28x%2Cy%29 "(x,y)")–coordinates
of a surface. Analyzing point reference data aims to capture variability
and correlation in observed phenomena and predict values at unobserved
locations while assessing uncertainty. Its applications are widespread
in environmental monitoring and geophysical studies. For example,
weather stations record temperature, humidity, and air quality at some
fixed monitoring sites across a geographical area, and data analysis
often aims to obtain a predicted surface of the phenomena by estimating
values at unobserved locations. Spatial data provide valuable insights
into the geographical distribution of phenomena. However, spatial data
analysis presents significant challenges, underscoring the critical need
for a deep understanding of statistical methods and computational tools.
Proficiency in navigating these intricacies is vital to unlocking the
full potential of these data types.

Statistical modeling of point reference data necessitates specifying a
random surface . One method to achieve this is through basis function
representations, including splines, wavelets, and radial basis
functions. An alternative, widely adopted approach involves modeling the
surface as a realization of a stochastic process. Gaussian processes
provide a practical framework for such modeling, offering a versatile
tool for spatial processes. They facilitate straightforward inference
and prediction by capturing spatial correlations, interpolating data,
modeling variability, and enabling probabilistic inference. Therefore,
in the following, we will outline Gaussian processes and explore their
application in modeling spatial processes. We will then review the
literature on non-Gaussian spatial and spatio-temporal modeling.
Finally, we will outline the objectives of this thesis.

## 1.3 Modeling Geostatistical Data

Analysis of spatial data observed over a finite set of locations over a
fixed domain often assumes that the measurement variable is,
theoretically, defined at locations varying continuously across the
domain. Consequently, a Gaussian process is thus a valuable tool for
modeling spatial data. For this topic, several textbooks are available,
including~, while illustrates how a Gaussian process is the most
valuable tool in the analysis of spatial data. Within this framework,
observations over a finite set of locations in a spatial domain are
assumed to be partial realizations of a spatial Gaussian process
![\\y(\mathbf{s}): \mathbf{s} \in \mathcal{D}\\](https://latex.codecogs.com/svg.image?%5C%7By%28%5Cmathbf%7Bs%7D%29%3A%20%5Cmathbf%7Bs%7D%20%5Cin%20%5Cmathcal%7BD%7D%5C%7D "\{y(\mathbf{s}): \mathbf{s} \in \mathcal{D}\}")
that defined on spatial index
![\mathbf{s}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D "\mathbf{s}")
varying continuously throughout
![d](https://latex.codecogs.com/svg.image?d "d")–dimensional domain
![\mathcal{D} \in \mathbb{R}^{d}](https://latex.codecogs.com/svg.image?%5Cmathcal%7BD%7D%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bd%7D "\mathcal{D} \in \mathbb{R}^{d}").
Therefore, the joint distribution of measurements at any finite set of
locations is assumed to be multivariate normal, and modeling only
requires specifying the mean and a valid covariance function. The
properties of multivariate normal distribution ensure the closed-form
marginal and conditionals, leading to straightforward computation for
model fitting and prediction.

To illustrate the above, consider that the measurement
![y(\mathbf{s})](https://latex.codecogs.com/svg.image?y%28%5Cmathbf%7Bs%7D%29 "y(\mathbf{s})"),
at location
![\mathbf{s}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D "\mathbf{s}"),
is generated as

![\begin{align}
\label{chap2:gp_marginal_model}
y(\mathbf{s}) = \mathbf{x}(\mathbf{s})^\prime \boldsymbol{\theta} + z(\mathbf{s}) + \epsilon(\mathbf{s}),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Clabel%7Bchap2%3Agp_marginal_model%7D%0Ay%28%5Cmathbf%7Bs%7D%29%20%3D%20%5Cmathbf%7Bx%7D%28%5Cmathbf%7Bs%7D%29%5E%5Cprime%20%5Cboldsymbol%7B%5Ctheta%7D%20%2B%20z%28%5Cmathbf%7Bs%7D%29%20%2B%20%5Cepsilon%28%5Cmathbf%7Bs%7D%29%2C%0A%5Cend%7Balign%7D "\begin{align}
\label{chap2:gp_marginal_model}
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

In practice, a single partial realization of a spatial Gaussian process
is available to infer the parameters and prediction.

### 1.3.1 Inference Procedure

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
In practice, the distribution
![\pi(\boldsymbol{\Phi} \mid \mathbf{y})](https://latex.codecogs.com/svg.image?%5Cpi%28%5Cboldsymbol%7B%5CPhi%7D%20%5Cmid%20%5Cmathbf%7By%7D%29 "\pi(\boldsymbol{\Phi} \mid \mathbf{y})")
does not have a closed-form, and Markov chain Monte Carlo (MCMC)
sampling methods are commonly employed to approximate this distribution.
These methods are straightforward to implement using modern statistical
computing platforms such as , , , and . MCMC methods provide samples
from the posterior distribution, which can be used to estimate various
summary statistics. Once samples from the posterior distribution are
available, predictions to unobserved locations follow straightforwardly.

## 1.4 Latent GP in Stan

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

### 1.4.1 Marginalization of Latent Process

### 1.4.2 Spatial Interpolation

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

### 1.4.3 Recovery of the Latent Component

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

## 1.5 Response GP in Stan

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

### 1.5.1 Computational complexity in analysing large datasets

## 1.6 Vecchia’s approximation and NNGP

Datta et al. ([2016](#ref-datta2016hierarchical)) developed the NNGP as
a sparse approximation of to a full GP. It generalizes the idea of
Vecchia ([1988](#ref-vecchia1988estimation)) from nearest neighbor
approximation of a data likelihood to the nearest neighbor approximation
of the likelihood of realizations of the process
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})"),
where the nearest neighbors set of
![\mathbf{s}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D "\mathbf{s}")
is defined based upon the Euclidean distance ([Datta
2022](#ref-datta2022nearest)). Both the nearest neighbor approximations
builds upon the idea that the joint distribution for a random vector
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
can be looked upon as a directed acyclic graph (DAG) to construct sparse
models for
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}")
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
![\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29 "\mathbb{N}(\mathbf{s})(\mathbf{s}_i)")
denotes the set of ![m](https://latex.codecogs.com/svg.image?m "m")
nearest neighbors of them. Then the NNGP is defined as the nearest
neighbor approximation of the likelihood of
![z(\mathbf{s}\_1),\ldots,z(\mathbf{s}\_n)](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_1%29%2C%5Cldots%2Cz%28%5Cmathbf%7Bs%7D_n%29 "z(\mathbf{s}_1),\ldots,z(\mathbf{s}_n)")
is given by

![\begin{align}
\label{eq_nngp_lik}
f(z(\mathbf{s}\_1)) \prod\_{i=2}^{n} f(z(\mathbf{s}\_i) \mid \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)})
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Clabel%7Beq_nngp_lik%7D%0Af%28z%28%5Cmathbf%7Bs%7D_1%29%29%20%5Cprod_%7Bi%3D2%7D%5E%7Bn%7D%20f%28z%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%29%0A%5Cend%7Balign%7D "\begin{align}
\label{eq_nngp_lik}
f(z(\mathbf{s}_1)) \prod_{i=2}^{n} f(z(\mathbf{s}_i) \mid \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)})
\end{align}")

where
![z(\mathbf{s}\_1) \sim \mathcal{N}(0,C\_{\mathbf{s}\_1,\mathbf{s}\_1})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_1%29%20%5Csim%20%5Cmathcal%7BN%7D%280%2CC_%7B%5Cmathbf%7Bs%7D_1%2C%5Cmathbf%7Bs%7D_1%7D%29 "z(\mathbf{s}_1) \sim \mathcal{N}(0,C_{\mathbf{s}_1,\mathbf{s}_1})")
and
![\mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D "\mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}")
is the set of ![m](https://latex.codecogs.com/svg.image?m "m") nearest
neighbor of the realizations. The conditional distribution
![f(z(\mathbf{s}\_i) \mid \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)})](https://latex.codecogs.com/svg.image?f%28z%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%29 "f(z(\mathbf{s}_i) \mid \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)})")
is an univariate normal
![\mathcal{N}(z(\mathbf{s}\_i) \| \mu\_{z(\mathbf{s}\_i) \mid \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}}, \sigma^2\_{z(\mathbf{s}\_i) \mid \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}})](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%28z%28%5Cmathbf%7Bs%7D_i%29%20%7C%20%5Cmu_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%2C%20%5Csigma%5E2_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%29 "\mathcal{N}(z(\mathbf{s}_i) | \mu_{z(\mathbf{s}_i) \mid \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}}, \sigma^2_{z(\mathbf{s}_i) \mid \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}})")
and which derived from the following multivariate normal

![\begin{align}
\begin{pmatrix}
z(\mathbf{s}\_i)\\
\mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}\\
\end{pmatrix} \sim \mathcal{N}\_{m+1}\left(
\begin{bmatrix}
0\\
0
\end{bmatrix},
\begin{bmatrix}
C\_{\mathbf{s}\_i\\\mathbf{s}\_i} & \mathbf{C}^\top\_{\mathbf{s}\_i,\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}\\
\mathbf{C}\_{\mathbf{s}\_i,\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)} & \mathbf{C}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i),\mathbb{N}(\mathbf{s})({\mathbf{s}\_i})}
\end{bmatrix}
\right),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cbegin%7Bpmatrix%7D%0Az%28%5Cmathbf%7Bs%7D_i%29%5C%5C%0A%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%5C%5C%0A%5Cend%7Bpmatrix%7D%20%5Csim%20%5Cmathcal%7BN%7D_%7Bm%2B1%7D%5Cleft%28%0A%5Cbegin%7Bbmatrix%7D%0A0%5C%5C%0A0%0A%5Cend%7Bbmatrix%7D%2C%0A%5Cbegin%7Bbmatrix%7D%0AC_%7B%5Cmathbf%7Bs%7D_i%5C%2C%5Cmathbf%7Bs%7D_i%7D%20%26%20%5Cmathbf%7BC%7D%5E%5Ctop_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%5C%5C%0A%5Cmathbf%7BC%7D_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%20%26%20%5Cmathbf%7BC%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%7B%5Cmathbf%7Bs%7D_i%7D%29%7D%0A%5Cend%7Bbmatrix%7D%0A%5Cright%29%2C%0A%5Cend%7Balign%7D "\begin{align}
\begin{pmatrix}
z(\mathbf{s}_i)\\
\mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}\\
\end{pmatrix} \sim \mathcal{N}_{m+1}\left(
\begin{bmatrix}
0\\
0
\end{bmatrix},
\begin{bmatrix}
C_{\mathbf{s}_i\,\mathbf{s}_i} & \mathbf{C}^\top_{\mathbf{s}_i,\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}\\
\mathbf{C}_{\mathbf{s}_i,\mathbb{N}(\mathbf{s})(\mathbf{s}_i)} & \mathbf{C}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i),\mathbb{N}(\mathbf{s})({\mathbf{s}_i})}
\end{bmatrix}
\right),
\end{align}")

where

![\begin{align}
\label{eq_nngp_con_moments}
\begin{split}
\mu\_{z(\mathbf{s}\_i) \mid \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}} &= \mathbf{C}^\top\_{\mathbf{s}\_i,\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)} \mathbf{C}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i),\mathcal{N}(\mathbf{s})({\mathbf{s}\_i})}^{-1} \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}\\
\sigma^2\_{z(\mathbf{s}\_i) \mid \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}} &= C\_{\mathbf{s}\_i\\\mathbf{s}\_i} - \mathbf{C}^\top\_{\mathbf{s}\_i,\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)} \mathbf{C}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i),\mathcal{N}(\mathbf{s})({\mathbf{s}\_i})}^{-1} \mathbf{C}\_{\mathbf{s}\_i,\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)}\\
\end{split}
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Clabel%7Beq_nngp_con_moments%7D%0A%5Cbegin%7Bsplit%7D%0A%5Cmu_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%20%26%3D%20%5Cmathbf%7BC%7D%5E%5Ctop_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%20%5Cmathbf%7BC%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%7B%5Cmathbf%7Bs%7D_i%7D%29%7D%5E%7B-1%7D%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%5C%5C%0A%5Csigma%5E2_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%20%26%3D%20C_%7B%5Cmathbf%7Bs%7D_i%5C%2C%5Cmathbf%7Bs%7D_i%7D%20-%20%5Cmathbf%7BC%7D%5E%5Ctop_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%20%5Cmathbf%7BC%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%7B%5Cmathbf%7Bs%7D_i%7D%29%7D%5E%7B-1%7D%20%5Cmathbf%7BC%7D_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%5C%5C%0A%5Cend%7Bsplit%7D%0A%5Cend%7Balign%7D "\begin{align}
\label{eq_nngp_con_moments}
\begin{split}
\mu_{z(\mathbf{s}_i) \mid \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}} &= \mathbf{C}^\top_{\mathbf{s}_i,\mathcal{N}(\mathbf{s})(\mathbf{s}_i)} \mathbf{C}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i),\mathcal{N}(\mathbf{s})({\mathbf{s}_i})}^{-1} \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}\\
\sigma^2_{z(\mathbf{s}_i) \mid \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}} &= C_{\mathbf{s}_i\,\mathbf{s}_i} - \mathbf{C}^\top_{\mathbf{s}_i,\mathcal{N}(\mathbf{s})(\mathbf{s}_i)} \mathbf{C}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i),\mathcal{N}(\mathbf{s})({\mathbf{s}_i})}^{-1} \mathbf{C}_{\mathbf{s}_i,\mathcal{N}(\mathbf{s})(\mathbf{s}_i)}\\
\end{split}
\end{align}")

The conditional distribution
![f(z(\mathbf{s}\_i) \mid \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)})](https://latex.codecogs.com/svg.image?f%28z%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%29 "f(z(\mathbf{s}_i) \mid \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)})")
can also be viewed as the likelihood from the generative model
![z(\mathbf{s}\_i) = \boldsymbol{a}\_i \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)} + e(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_i%29%20%3D%20%5Cboldsymbol%7Ba%7D_i%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%20%2B%20e%28%5Cmathbf%7Bs%7D_i%29 "z(\mathbf{s}_i) = \boldsymbol{a}_i \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)} + e(\mathbf{s}_i)"),
where
![\boldsymbol{a}\_i = \mathbf{C}^\top\_{\mathbf{s}\_i,\mathcal{N}(\mathbf{s})(\mathbf{s}\_i)} \mathbf{C}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i),\mathcal{N}(\mathbf{s})({\mathbf{s}\_i})}^{-1}](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7Ba%7D_i%20%3D%20%5Cmathbf%7BC%7D%5E%5Ctop_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%20%5Cmathbf%7BC%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%7B%5Cmathbf%7Bs%7D_i%7D%29%7D%5E%7B-1%7D "\boldsymbol{a}_i = \mathbf{C}^\top_{\mathbf{s}_i,\mathcal{N}(\mathbf{s})(\mathbf{s}_i)} \mathbf{C}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i),\mathcal{N}(\mathbf{s})({\mathbf{s}_i})}^{-1}")
and
![e(\mathbf{s}\_i) \sim i.i.d\\ \mathcal{N}(0,d\_{i})](https://latex.codecogs.com/svg.image?e%28%5Cmathbf%7Bs%7D_i%29%20%5Csim%20i.i.d%5C%2C%20%5Cmathcal%7BN%7D%280%2Cd_%7Bi%7D%29 "e(\mathbf{s}_i) \sim i.i.d\, \mathcal{N}(0,d_{i})")
where
![d_1 = C\_{\mathbf{s}\_i,\mathbf{s}\_i}](https://latex.codecogs.com/svg.image?d_1%20%3D%20C_%7B%5Cmathbf%7Bs%7D_i%2C%5Cmathbf%7Bs%7D_i%7D "d_1 = C_{\mathbf{s}_i,\mathbf{s}_i}")
and
![d_i = \sigma^2\_{z(\mathbf{s}\_i) \mid \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}}\\ i=2,\ldots,n](https://latex.codecogs.com/svg.image?d_i%20%3D%20%5Csigma%5E2_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%5C%3B%20i%3D2%2C%5Cldots%2Cn "d_i = \sigma^2_{z(\mathbf{s}_i) \mid \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}}\; i=2,\ldots,n").

Therefore, the likelihood in can be written as the following set of
linear models,
![z(\mathbf{s}\_1) = 0 + e(\mathbf{s}\_1)](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_1%29%20%3D%200%20%2B%20e%28%5Cmathbf%7Bs%7D_1%29 "z(\mathbf{s}_1) = 0 + e(\mathbf{s}_1)")
and for
![i=2,\ldots,n](https://latex.codecogs.com/svg.image?i%3D2%2C%5Cldots%2Cn "i=2,\ldots,n"),

![\begin{align}
z(\mathbf{s}\_i) = a\_{i1} z(\mathbf{s}\_1) + a\_{i2} z(\mathbf{s}\_2) + \cdots + a\_{i,i-1} z(\mathbf{s}\_{i-1}) + e(\mathbf{s}\_i).
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0Az%28%5Cmathbf%7Bs%7D_i%29%20%3D%20a_%7Bi1%7D%20z%28%5Cmathbf%7Bs%7D_1%29%20%2B%20a_%7Bi2%7D%20z%28%5Cmathbf%7Bs%7D_2%29%20%2B%20%5Ccdots%20%2B%20a_%7Bi%2Ci-1%7D%20z%28%5Cmathbf%7Bs%7D_%7Bi-1%7D%29%20%2B%20e%28%5Cmathbf%7Bs%7D_i%29.%0A%5Cend%7Balign%7D "\begin{align}
z(\mathbf{s}_i) = a_{i1} z(\mathbf{s}_1) + a_{i2} z(\mathbf{s}_2) + \cdots + a_{i,i-1} z(\mathbf{s}_{i-1}) + e(\mathbf{s}_i).
\end{align}")

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
![\mathbf{z} = \boldsymbol{A} \mathbf{z} + \boldsymbol{e}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D%20%3D%20%5Cboldsymbol%7BA%7D%20%5Cmathbf%7Bz%7D%20%2B%20%5Cboldsymbol%7Be%7D "\mathbf{z} = \boldsymbol{A} \mathbf{z} + \boldsymbol{e}"),
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
![\mathbf{z} \sim \mathcal{N}\left(\mathbf{0},(\mathbf{I}-\boldsymbol{A})^{-1}\boldsymbol{D}(\mathbf{I}-\boldsymbol{A})^{-\top}\right)](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D%20%5Csim%20%5Cmathcal%7BN%7D%5Cleft%28%5Cmathbf%7B0%7D%2C%28%5Cmathbf%7BI%7D-%5Cboldsymbol%7BA%7D%29%5E%7B-1%7D%5Cboldsymbol%7BD%7D%28%5Cmathbf%7BI%7D-%5Cboldsymbol%7BA%7D%29%5E%7B-%5Ctop%7D%5Cright%29 "\mathbf{z} \sim \mathcal{N}\left(\mathbf{0},(\mathbf{I}-\boldsymbol{A})^{-1}\boldsymbol{D}(\mathbf{I}-\boldsymbol{A})^{-\top}\right)")
corresponds to a generative multivariate Gaussian model
![\mathbf{z} \sim \mathcal{N}\left(\mathbf{0},\widetilde{\mathbf{C}}\right)](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D%20%5Csim%20%5Cmathcal%7BN%7D%5Cleft%28%5Cmathbf%7B0%7D%2C%5Cwidetilde%7B%5Cmathbf%7BC%7D%7D%5Cright%29 "\mathbf{z} \sim \mathcal{N}\left(\mathbf{0},\widetilde{\mathbf{C}}\right)")
for the process realizations, where
![\widetilde{\mathbf{C}} = (\mathbf{I}-\boldsymbol{A})^{-1}\boldsymbol{D}(\mathbf{I}-\boldsymbol{A})^{-\top}](https://latex.codecogs.com/svg.image?%5Cwidetilde%7B%5Cmathbf%7BC%7D%7D%20%3D%20%28%5Cmathbf%7BI%7D-%5Cboldsymbol%7BA%7D%29%5E%7B-1%7D%5Cboldsymbol%7BD%7D%28%5Cmathbf%7BI%7D-%5Cboldsymbol%7BA%7D%29%5E%7B-%5Ctop%7D "\widetilde{\mathbf{C}} = (\mathbf{I}-\boldsymbol{A})^{-1}\boldsymbol{D}(\mathbf{I}-\boldsymbol{A})^{-\top}")
([Datta 2022](#ref-datta2022nearest)). Datta et al.
([2016](#ref-datta2016hierarchical)) shows that this is a valid
approximation of the likelihood corresponding to the GP model
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})")
with zero mean and covariance function
![C](https://latex.codecogs.com/svg.image?C "C"). However, specification
of a valid GP over the entire domain was completed by defining
prediction distribution of
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})")
at new locations conditional on
![\mathbf{z}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D "\mathbf{z}").

For the prediction of the latent process
![z(\mathbf{s})](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D%29 "z(\mathbf{s})")
at a set of ![n_0](https://latex.codecogs.com/svg.image?n_0 "n_0") new
locations
![\\\mathbf{s}\_{01},\ldots,\mathbf{s}\_{0n_0}\\ \notin \\\mathbf{s}\_1,\ldots,\mathbf{s}\_n\\](https://latex.codecogs.com/svg.image?%5C%7B%5Cmathbf%7Bs%7D_%7B01%7D%2C%5Cldots%2C%5Cmathbf%7Bs%7D_%7B0n_0%7D%5C%7D%20%5Cnotin%20%5C%7B%5Cmathbf%7Bs%7D_1%2C%5Cldots%2C%5Cmathbf%7Bs%7D_n%5C%7D "\{\mathbf{s}_{01},\ldots,\mathbf{s}_{0n_0}\} \notin \{\mathbf{s}_1,\ldots,\mathbf{s}_n\}"),
Datta et al. ([2016](#ref-datta2016hierarchical)) specified the
conditional distribution of
![z(\mathbf{s}\_{0i}) \| \mathbf{z}](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%20%7C%20%5Cmathbf%7Bz%7D "z(\mathbf{s}_{0i}) | \mathbf{z}")
independently, which is given by
![z(\mathbf{s}\_{0i}) \| \mathbf{z} \sim \mathcal{N}\left(\mu\_{z(\mathbf{s}\_{0i}) \| \mathbf{z}}, \sigma^2\_{z(\mathbf{s}\_{0i}) \| \mathbf{z}}\right)](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%20%7C%20%5Cmathbf%7Bz%7D%20%5Csim%20%5Cmathcal%7BN%7D%5Cleft%28%5Cmu_%7Bz%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%20%7C%20%5Cmathbf%7Bz%7D%7D%2C%20%5Csigma%5E2_%7Bz%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%20%7C%20%5Cmathbf%7Bz%7D%7D%5Cright%29 "z(\mathbf{s}_{0i}) | \mathbf{z} \sim \mathcal{N}\left(\mu_{z(\mathbf{s}_{0i}) | \mathbf{z}}, \sigma^2_{z(\mathbf{s}_{0i}) | \mathbf{z}}\right)"),
where

![\begin{align}
\sigma^2\_{z(\mathbf{s}\_{0i}) \| \mathbf{z}} = \mathbf{C}^\top\_{\mathbf{s}\_{0i},\mathcal{N}(\mathbf{s})(\mathbf{s}\_{0i})} \mathbf{C}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_{0i})}
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Csigma%5E2_%7Bz%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%20%7C%20%5Cmathbf%7Bz%7D%7D%20%3D%20%5Cmathbf%7BC%7D%5E%5Ctop_%7B%5Cmathbf%7Bs%7D_%7B0i%7D%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%7D%20%5Cmathbf%7BC%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%7D%0A%5Cend%7Balign%7D "\begin{align}
\sigma^2_{z(\mathbf{s}_{0i}) | \mathbf{z}} = \mathbf{C}^\top_{\mathbf{s}_{0i},\mathcal{N}(\mathbf{s})(\mathbf{s}_{0i})} \mathbf{C}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_{0i})}
\end{align}")

![\begin{align}
\sigma^2\_{z(\mathbf{s}\_{0i}) \| \mathbf{z}} = \mathcal{N}(\mathbf{s})(\mathbf{s}\_{0i})}^{-1} \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_{0i})}, C\_{\mathbf{s}\_{0i}\\\mathbf{s}\_{0i}} - \mathbf{C}^\top\_{\mathbf{s}\_{0i},\mathcal{N}(\mathbf{s})(\mathbf{s}\_{0i})} \mathbf{C}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_{0i}),\mathcal{N}(\mathbf{s})({\mathbf{s}\_{0i}})}^{-1} \mathbf{C}\_{\mathbf{s}\_{0i},\mathcal{N}(\mathbf{s})(\mathbf{s}\_{0i})}\right),
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Csigma%5E2_%7Bz%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%20%7C%20%5Cmathbf%7Bz%7D%7D%20%3D%20%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%7D%5E%7B-1%7D%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%7D%2C%20C_%7B%5Cmathbf%7Bs%7D_%7B0i%7D%5C%2C%5Cmathbf%7Bs%7D_%7B0i%7D%7D%20-%20%5Cmathbf%7BC%7D%5E%5Ctop_%7B%5Cmathbf%7Bs%7D_%7B0i%7D%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%7D%20%5Cmathbf%7BC%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%7B%5Cmathbf%7Bs%7D_%7B0i%7D%7D%29%7D%5E%7B-1%7D%20%5Cmathbf%7BC%7D_%7B%5Cmathbf%7Bs%7D_%7B0i%7D%2C%5Cmathcal%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_%7B0i%7D%29%7D%5Cright%29%2C%0A%5Cend%7Balign%7D "\begin{align}
\sigma^2_{z(\mathbf{s}_{0i}) | \mathbf{z}} = \mathcal{N}(\mathbf{s})(\mathbf{s}_{0i})}^{-1} \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_{0i})}, C_{\mathbf{s}_{0i}\,\mathbf{s}_{0i}} - \mathbf{C}^\top_{\mathbf{s}_{0i},\mathcal{N}(\mathbf{s})(\mathbf{s}_{0i})} \mathbf{C}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_{0i}),\mathcal{N}(\mathbf{s})({\mathbf{s}_{0i}})}^{-1} \mathbf{C}_{\mathbf{s}_{0i},\mathcal{N}(\mathbf{s})(\mathbf{s}_{0i})}\right),
\end{align}")

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
![\mathcal{N}\left(z(\mathbf{s}\_i) \| \mu\_{z(\mathbf{s}\_i) \mid \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}}, \sigma^2\_{z(\mathbf{s}\_i) \mid \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}}\right)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%5Cleft%28z%28%5Cmathbf%7Bs%7D_i%29%20%7C%20%5Cmu_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%2C%20%5Csigma%5E2_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%5Cright%29 "\mathcal{N}\left(z(\mathbf{s}_i) | \mu_{z(\mathbf{s}_i) \mid \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}}, \sigma^2_{z(\mathbf{s}_i) \mid \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}}\right)")
in is parameterized as
![z(\mathbf{s}\_i) \|  z\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)} = \mu\_{z(\mathbf{s}\_i) \mid \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}} + \sigma\_{z(\mathbf{s}\_i) \mid \mathbf{z}\_{\mathbb{N}(\mathbf{s})(\mathbf{s}\_i)}} v(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?z%28%5Cmathbf%7Bs%7D_i%29%20%7C%20%20z_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%20%3D%20%5Cmu_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%20%2B%20%5Csigma_%7Bz%28%5Cmathbf%7Bs%7D_i%29%20%5Cmid%20%5Cmathbf%7Bz%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D%29%28%5Cmathbf%7Bs%7D_i%29%7D%7D%20v%28%5Cmathbf%7Bs%7D_i%29 "z(\mathbf{s}_i) |  z_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)} = \mu_{z(\mathbf{s}_i) \mid \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}} + \sigma_{z(\mathbf{s}_i) \mid \mathbf{z}_{\mathbb{N}(\mathbf{s})(\mathbf{s}_i)}} v(\mathbf{s}_i)")
where
![v(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?v%28%5Cmathbf%7Bs%7D_i%29 "v(\mathbf{s}_i)")
is independent Gaussian noise with mean zero and variance one.

One main interest in spatial data analysis is obtaining an estimated
process surface through pointwise prediction. Suppose that a vector of
![p](https://latex.codecogs.com/svg.image?p "p") covariates values
![\mathbf{x}(\mathbf{s}\_0)](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bx%7D%28%5Cmathbf%7Bs%7D_0%29 "\mathbf{x}(\mathbf{s}_0)")
at a generic prediction location
![\mathbf{s}\_0 \in \mathcal{D}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_0%20%5Cin%20%5Cmathcal%7BD%7D "\mathbf{s}_0 \in \mathcal{D}")
is available. Following model, the posterior predictive distribution of
![y(\mathbf{s})](https://latex.codecogs.com/svg.image?y%28%5Cmathbf%7Bs%7D%29 "y(\mathbf{s})")
at
![\mathbf{s}\_0](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_0 "\mathbf{s}_0")
is given by

![\begin{align}
\pi(y(\mathbf{s}\_0) \mid \mathbf{y}) = \int \mathcal{N}\left(y(\mathbf{s}\_0) \mid \mu\_{y(\mathbf{s}\_0) \mid \mathbf{y}}, \tau^2 \right) \\ \pi(z_1(\mathbf{s}\_0) \mid \mathbf{z}\_1, \sigma_1, \ell_1) \\ \pi(z_2(\mathbf{s}\_0) \mid \mathbf{z}\_2, \sigma_2, \ell_2) \\
\pi(\boldsymbol{\Theta}, \mathbf{z}\_1, \mathbf{z}\_2 \mid \mathbf{y}) \\ dz_1(\mathbf{s}\_0)\\ dz_2(\mathbf{s}\_0)\\ d\mathbf{z}\_1 \\ d\mathbf{z}\_2 \\ d\boldsymbol{\Theta},
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cpi%28y%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7By%7D%29%20%3D%20%5Cint%20%5Cmathcal%7BN%7D%5Cleft%28y%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmu_%7By%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7By%7D%7D%2C%20%5Ctau%5E2%20%5Cright%29%20%5C%3B%20%5Cpi%28z_1%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7Bz%7D_1%2C%20%5Csigma_1%2C%20%5Cell_1%29%20%5C%3B%20%5Cpi%28z_2%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7Bz%7D_2%2C%20%5Csigma_2%2C%20%5Cell_2%29%20%5C%5C%0A%5Cpi%28%5Cboldsymbol%7B%5CTheta%7D%2C%20%5Cmathbf%7Bz%7D_1%2C%20%5Cmathbf%7Bz%7D_2%20%5Cmid%20%5Cmathbf%7By%7D%29%20%5C%3B%20dz_1%28%5Cmathbf%7Bs%7D_0%29%5C%3B%20dz_2%28%5Cmathbf%7Bs%7D_0%29%5C%3B%20d%5Cmathbf%7Bz%7D_1%20%5C%3B%20d%5Cmathbf%7Bz%7D_2%20%5C%3B%20d%5Cboldsymbol%7B%5CTheta%7D%2C%0A%5Cend%7Balign%7D "\begin{align}
\pi(y(\mathbf{s}_0) \mid \mathbf{y}) = \int \mathcal{N}\left(y(\mathbf{s}_0) \mid \mu_{y(\mathbf{s}_0) \mid \mathbf{y}}, \tau^2 \right) \; \pi(z_1(\mathbf{s}_0) \mid \mathbf{z}_1, \sigma_1, \ell_1) \; \pi(z_2(\mathbf{s}_0) \mid \mathbf{z}_2, \sigma_2, \ell_2) \\
\pi(\boldsymbol{\Theta}, \mathbf{z}_1, \mathbf{z}_2 \mid \mathbf{y}) \; dz_1(\mathbf{s}_0)\; dz_2(\mathbf{s}_0)\; d\mathbf{z}_1 \; d\mathbf{z}_2 \; d\boldsymbol{\Theta},
\end{align}")

where
![\mu\_{y(\mathbf{s}\_0)\mid\mathbf{y}} = \mathbf{x}'(\mathbf{s}\_0)\boldsymbol{\theta} + \gamma \exp\\z_1(\mathbf{s}\_0)\\ + z_2(\mathbf{s}\_0)](https://latex.codecogs.com/svg.image?%5Cmu_%7By%28%5Cmathbf%7Bs%7D_0%29%5Cmid%5Cmathbf%7By%7D%7D%20%3D%20%5Cmathbf%7Bx%7D%27%28%5Cmathbf%7Bs%7D_0%29%5Cboldsymbol%7B%5Ctheta%7D%20%2B%20%5Cgamma%20%5Cexp%5C%7Bz_1%28%5Cmathbf%7Bs%7D_0%29%5C%7D%20%2B%20z_2%28%5Cmathbf%7Bs%7D_0%29 "\mu_{y(\mathbf{s}_0)\mid\mathbf{y}} = \mathbf{x}'(\mathbf{s}_0)\boldsymbol{\theta} + \gamma \exp\{z_1(\mathbf{s}_0)\} + z_2(\mathbf{s}_0)")
and
![\boldsymbol{\Theta} = \\\boldsymbol{\theta}, \gamma, \sigma_1, \ell_1, \sigma_2, \ell_2, \tau\\](https://latex.codecogs.com/svg.image?%5Cboldsymbol%7B%5CTheta%7D%20%3D%20%5C%7B%5Cboldsymbol%7B%5Ctheta%7D%2C%20%5Cgamma%2C%20%5Csigma_1%2C%20%5Cell_1%2C%20%5Csigma_2%2C%20%5Cell_2%2C%20%5Ctau%5C%7D "\boldsymbol{\Theta} = \{\boldsymbol{\theta}, \gamma, \sigma_1, \ell_1, \sigma_2, \ell_2, \tau\}").
The integral in does not have a closed form; therefore, realizations
from the posterior predictive distribution of
![y(\mathbf{s}\_0)](https://latex.codecogs.com/svg.image?y%28%5Cmathbf%7Bs%7D_0%29 "y(\mathbf{s}_0)")
are obtained using composition sampling. For each set of posterior
samples of
![\\\boldsymbol{\Theta}, \mathbf{z}\_1, \mathbf{z}\_2\\](https://latex.codecogs.com/svg.image?%5C%7B%5Cboldsymbol%7B%5CTheta%7D%2C%20%5Cmathbf%7Bz%7D_1%2C%20%5Cmathbf%7Bz%7D_2%5C%7D "\{\boldsymbol{\Theta}, \mathbf{z}_1, \mathbf{z}_2\}"),
a sample from the posterior predictive can be obtained by sampling from
![\pi(z_k(\mathbf{s}\_0) \mid \mathbf{z}\_k, \sigma_k, \ell_k)](https://latex.codecogs.com/svg.image?%5Cpi%28z_k%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7Bz%7D_k%2C%20%5Csigma_k%2C%20%5Cell_k%29 "\pi(z_k(\mathbf{s}_0) \mid \mathbf{z}_k, \sigma_k, \ell_k)")
for
![k = 1, 2](https://latex.codecogs.com/svg.image?k%20%3D%201%2C%202 "k = 1, 2"),
and
![\mathcal{N}\left(y(\mathbf{s}\_0) \mid \mu\_{y(\mathbf{s}\_0)\mid\mathbf{y}}, \tau^2 \right)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%5Cleft%28y%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmu_%7By%28%5Cmathbf%7Bs%7D_0%29%5Cmid%5Cmathbf%7By%7D%7D%2C%20%5Ctau%5E2%20%5Cright%29 "\mathcal{N}\left(y(\mathbf{s}_0) \mid \mu_{y(\mathbf{s}_0)\mid\mathbf{y}}, \tau^2 \right)").
Note that once the value of
![z_1(\mathbf{s}\_0)](https://latex.codecogs.com/svg.image?z_1%28%5Cmathbf%7Bs%7D_0%29 "z_1(\mathbf{s}_0)")
and
![z_2(\mathbf{s}\_0)](https://latex.codecogs.com/svg.image?z_2%28%5Cmathbf%7Bs%7D_0%29 "z_2(\mathbf{s}_0)")
are available, the predicted value of
![y(\mathbf{s}\_0)](https://latex.codecogs.com/svg.image?y%28%5Cmathbf%7Bs%7D_0%29 "y(\mathbf{s}_0)")
is easily sampled from the univariate distribution
![\mathcal{N}\left(y(\mathbf{s}\_0) \mid \mu\_{y(\mathbf{s}\_0)\mid\mathbf{y}}, \tau^2 \right)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%5Cleft%28y%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmu_%7By%28%5Cmathbf%7Bs%7D_0%29%5Cmid%5Cmathbf%7By%7D%7D%2C%20%5Ctau%5E2%20%5Cright%29 "\mathcal{N}\left(y(\mathbf{s}_0) \mid \mu_{y(\mathbf{s}_0)\mid\mathbf{y}}, \tau^2 \right)").

While following model, the posterior predictive distribution of
![y(\mathbf{s}\_0)](https://latex.codecogs.com/svg.image?y%28%5Cmathbf%7Bs%7D_0%29 "y(\mathbf{s}_0)")
at any unobserved location
![\mathbf{s}\_0 \in \mathcal{D}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_0%20%5Cin%20%5Cmathcal%7BD%7D "\mathbf{s}_0 \in \mathcal{D}")
is given by

![\begin{align}
\pi(y(\mathbf{s}\_0) \mid \mathbf{y}) &= \int \mathcal{N}\left(y(\mathbf{s}\_0) \mid \mu\_{y(\mathbf{s}\_0)\mid \mathbf{y}}, \sigma^2\_{y(\mathbf{s}\_0) \mid \mathbf{y}} \right)\\ \pi(z_1(\mathbf{s}\_0) \mid \mathbf{z}\_1, \sigma_1, \ell_1)\\
& \pi(\boldsymbol{\Theta}, \mathbf{z}\_1 \mid \mathbf{y}) \\ dz_1(\mathbf{s}\_0)\\ d\mathbf{z}\_1 \\ d\boldsymbol{\Theta},
\end{align}](https://latex.codecogs.com/svg.image?%5Cbegin%7Balign%7D%0A%5Cpi%28y%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7By%7D%29%20%26%3D%20%5Cint%20%5Cmathcal%7BN%7D%5Cleft%28y%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmu_%7By%28%5Cmathbf%7Bs%7D_0%29%5Cmid%20%5Cmathbf%7By%7D%7D%2C%20%5Csigma%5E2_%7By%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7By%7D%7D%20%5Cright%29%5C%3B%20%5Cpi%28z_1%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7Bz%7D_1%2C%20%5Csigma_1%2C%20%5Cell_1%29%5C%5C%0A%26%20%5Cpi%28%5Cboldsymbol%7B%5CTheta%7D%2C%20%5Cmathbf%7Bz%7D_1%20%5Cmid%20%5Cmathbf%7By%7D%29%20%5C%3B%20dz_1%28%5Cmathbf%7Bs%7D_0%29%5C%3B%20d%5Cmathbf%7Bz%7D_1%20%5C%3B%20d%5Cboldsymbol%7B%5CTheta%7D%2C%0A%5Cend%7Balign%7D "\begin{align}
\pi(y(\mathbf{s}_0) \mid \mathbf{y}) &= \int \mathcal{N}\left(y(\mathbf{s}_0) \mid \mu_{y(\mathbf{s}_0)\mid \mathbf{y}}, \sigma^2_{y(\mathbf{s}_0) \mid \mathbf{y}} \right)\; \pi(z_1(\mathbf{s}_0) \mid \mathbf{z}_1, \sigma_1, \ell_1)\\
& \pi(\boldsymbol{\Theta}, \mathbf{z}_1 \mid \mathbf{y}) \; dz_1(\mathbf{s}_0)\; d\mathbf{z}_1 \; d\boldsymbol{\Theta},
\end{align}")

which also does not have a closed form and composition sampling is used
to obtain samples from
![\pi(y(\mathbf{s}\_0) \mid \mathbf{y})](https://latex.codecogs.com/svg.image?%5Cpi%28y%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7By%7D%29 "\pi(y(\mathbf{s}_0) \mid \mathbf{y})").
It requires sampling from
![\pi(z_1(\mathbf{s}\_0) \mid \mathbf{z}\_1, \sigma_1, \ell_1)](https://latex.codecogs.com/svg.image?%5Cpi%28z_1%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7Bz%7D_1%2C%20%5Csigma_1%2C%20%5Cell_1%29 "\pi(z_1(\mathbf{s}_0) \mid \mathbf{z}_1, \sigma_1, \ell_1)")
at first and then from
![\mathcal{N}\bigl(y(\mathbf{s}\_0) \mid \mu\_{y(\mathbf{s}\_0)\mid \mathbf{y}}, \sigma^2\_{y(\mathbf{s}\_0) \mid \mathbf{y}}\bigr)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%5Cbigl%28y%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmu_%7By%28%5Cmathbf%7Bs%7D_0%29%5Cmid%20%5Cmathbf%7By%7D%7D%2C%20%5Csigma%5E2_%7By%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7By%7D%7D%5Cbigr%29 "\mathcal{N}\bigl(y(\mathbf{s}_0) \mid \mu_{y(\mathbf{s}_0)\mid \mathbf{y}}, \sigma^2_{y(\mathbf{s}_0) \mid \mathbf{y}}\bigr)")
for each posterior samples of
![\\\boldsymbol{\Theta}, \mathbf{z}\_1\\](https://latex.codecogs.com/svg.image?%5C%7B%5Cboldsymbol%7B%5CTheta%7D%2C%20%5Cmathbf%7Bz%7D_1%5C%7D "\{\boldsymbol{\Theta}, \mathbf{z}_1\}").
However, with known value of
![z_1(\mathbf{s}\_0)](https://latex.codecogs.com/svg.image?z_1%28%5Cmathbf%7Bs%7D_0%29 "z_1(\mathbf{s}_0)"),
sampling from the distribution
![\mathcal{N}\bigl(y(\mathbf{s}\_0) \mid \mu\_{y(\mathbf{s}\_0)\mid \mathbf{y}}, \sigma^2\_{y(\mathbf{s}\_0) \mid \mathbf{y}}\bigr)](https://latex.codecogs.com/svg.image?%5Cmathcal%7BN%7D%5Cbigl%28y%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmu_%7By%28%5Cmathbf%7Bs%7D_0%29%5Cmid%20%5Cmathbf%7By%7D%7D%2C%20%5Csigma%5E2_%7By%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7By%7D%7D%5Cbigr%29 "\mathcal{N}\bigl(y(\mathbf{s}_0) \mid \mu_{y(\mathbf{s}_0)\mid \mathbf{y}}, \sigma^2_{y(\mathbf{s}_0) \mid \mathbf{y}}\bigr)")
is computationally expensive. As this distribution is the univariate
conditional distribution of a multivariate normal distribution of
dimension
![(n+1)](https://latex.codecogs.com/svg.image?%28n%2B1%29 "(n+1)"),
computing conditional mean
![\mu\_{y(\mathbf{s}\_0) \mid \mathbf{y}}](https://latex.codecogs.com/svg.image?%5Cmu_%7By%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7By%7D%7D "\mu_{y(\mathbf{s}_0) \mid \mathbf{y}}")
and variance
![\sigma^2\_{y(\mathbf{s}\_0) \mid \mathbf{y}}](https://latex.codecogs.com/svg.image?%5Csigma%5E2_%7By%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7By%7D%7D "\sigma^2_{y(\mathbf{s}_0) \mid \mathbf{y}}")
involves expensive matrix calculations. This can be avoided using
Vecchia’s NN method for approximating the density of
![\mathbf{y}](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D "\mathbf{y}"),
in which the conditional mean and variance are approximated,
respectively, as
![\mu^{\text{NN}}\_{y(\mathbf{s}\_0)\mid \mathbf{y}} = \mathbf{x}(\mathbf{s}\_0)'\boldsymbol{\theta} + \gamma \exp\left\\z_1(\mathbf{s}\_0)\right\\ + \mathbf{c}\_{2,\mathbf{s}\_0,\mathbb{N}(\mathbf{s}\_0)} (\mathbf{C}\_{2,\mathbb{N}(\mathbf{s}\_0)} + (\tau^2/\sigma_2^2)\mathbf{I})^{-1}(\mathbf{y}\_{\mathbb{N}(\mathbf{s}\_0)} - \mathbf{X}\_{\mathbb{N}(\mathbf{s}\_0)}\boldsymbol{\theta} - \gamma\\\exp\\\mathbf{z}\_{1,\mathbb{N}(\mathbf{s}\_0)}\\)](https://latex.codecogs.com/svg.image?%5Cmu%5E%7B%5Ctext%7BNN%7D%7D_%7By%28%5Cmathbf%7Bs%7D_0%29%5Cmid%20%5Cmathbf%7By%7D%7D%20%3D%20%5Cmathbf%7Bx%7D%28%5Cmathbf%7Bs%7D_0%29%27%5Cboldsymbol%7B%5Ctheta%7D%20%2B%20%5Cgamma%20%5Cexp%5Cleft%5C%7Bz_1%28%5Cmathbf%7Bs%7D_0%29%5Cright%5C%7D%20%2B%20%5Cmathbf%7Bc%7D_%7B2%2C%5Cmathbf%7Bs%7D_0%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D%20%28%5Cmathbf%7BC%7D_%7B2%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D%20%2B%20%28%5Ctau%5E2%2F%5Csigma_2%5E2%29%5Cmathbf%7BI%7D%29%5E%7B-1%7D%28%5Cmathbf%7By%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D%20-%20%5Cmathbf%7BX%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D%5Cboldsymbol%7B%5Ctheta%7D%20-%20%5Cgamma%5C%2C%5Cexp%5C%7B%5Cmathbf%7Bz%7D_%7B1%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D%5C%7D%29 "\mu^{\text{NN}}_{y(\mathbf{s}_0)\mid \mathbf{y}} = \mathbf{x}(\mathbf{s}_0)'\boldsymbol{\theta} + \gamma \exp\left\{z_1(\mathbf{s}_0)\right\} + \mathbf{c}_{2,\mathbf{s}_0,\mathbb{N}(\mathbf{s}_0)} (\mathbf{C}_{2,\mathbb{N}(\mathbf{s}_0)} + (\tau^2/\sigma_2^2)\mathbf{I})^{-1}(\mathbf{y}_{\mathbb{N}(\mathbf{s}_0)} - \mathbf{X}_{\mathbb{N}(\mathbf{s}_0)}\boldsymbol{\theta} - \gamma\,\exp\{\mathbf{z}_{1,\mathbb{N}(\mathbf{s}_0)}\})")
and
![\sigma^{2, \text{NN}}\_{y(\mathbf{s}\_0) \mid \mathbf{y}} = \sigma_2^2\[1 + (\tau^2/\sigma_2^2) - \mathbf{c}\_{2,\mathbf{s}\_0,\mathbb{N}(\mathbf{s}\_0)}' (\mathbf{C}\_{2,\mathbb{N}(\mathbf{s}\_0)}^{-1} + (\tau^2/\sigma_2^2)\mathbf{I})^{-1}\mathbf{c}\_{2,\mathbf{s}\_0,\mathbb{N}(\mathbf{s}\_0)}\]](https://latex.codecogs.com/svg.image?%5Csigma%5E%7B2%2C%20%5Ctext%7BNN%7D%7D_%7By%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7By%7D%7D%20%3D%20%5Csigma_2%5E2%5B1%20%2B%20%28%5Ctau%5E2%2F%5Csigma_2%5E2%29%20-%20%5Cmathbf%7Bc%7D_%7B2%2C%5Cmathbf%7Bs%7D_0%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D%27%20%28%5Cmathbf%7BC%7D_%7B2%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D%5E%7B-1%7D%20%2B%20%28%5Ctau%5E2%2F%5Csigma_2%5E2%29%5Cmathbf%7BI%7D%29%5E%7B-1%7D%5Cmathbf%7Bc%7D_%7B2%2C%5Cmathbf%7Bs%7D_0%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D%5D "\sigma^{2, \text{NN}}_{y(\mathbf{s}_0) \mid \mathbf{y}} = \sigma_2^2[1 + (\tau^2/\sigma_2^2) - \mathbf{c}_{2,\mathbf{s}_0,\mathbb{N}(\mathbf{s}_0)}' (\mathbf{C}_{2,\mathbb{N}(\mathbf{s}_0)}^{-1} + (\tau^2/\sigma_2^2)\mathbf{I})^{-1}\mathbf{c}_{2,\mathbf{s}_0,\mathbb{N}(\mathbf{s}_0)}]"),
where
![\mathbb{N}(\mathbf{s}\_{0})](https://latex.codecogs.com/svg.image?%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_%7B0%7D%29 "\mathbb{N}(\mathbf{s}_{0})")
denotes the ![m](https://latex.codecogs.com/svg.image?m "m") nearest
neighbors of
![\mathbf{s}\_0](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_0 "\mathbf{s}_0")
in the observed locations
![\\\mathbf{s}\_1,\ldots,\mathbf{s}\_n\\](https://latex.codecogs.com/svg.image?%5C%7B%5Cmathbf%7Bs%7D_1%2C%5Cldots%2C%5Cmathbf%7Bs%7D_n%5C%7D "\{\mathbf{s}_1,\ldots,\mathbf{s}_n\}"),
![\mathbf{y}\_{\mathbb{N}(\mathbf{s}\_0)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D "\mathbf{y}_{\mathbb{N}(\mathbf{s}_0)}")
is a vector formed by stacking
![y(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?y%28%5Cmathbf%7Bs%7D_i%29 "y(\mathbf{s}_i)")
for
![\mathbf{s}\_i \in \mathbb{N}(\mathbf{s}\_0)](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_i%20%5Cin%20%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29 "\mathbf{s}_i \in \mathbb{N}(\mathbf{s}_0)"),
![\mathbf{X}\_{\mathbb{N}(\mathbf{s}\_0)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BX%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D "\mathbf{X}_{\mathbb{N}(\mathbf{s}_0)}")
is the design matrix corresponding to
![\mathbf{y}\_{\mathbb{N}(\mathbf{s}\_0)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7By%7D_%7B%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D "\mathbf{y}_{\mathbb{N}(\mathbf{s}_0)}").
With
![k = 1, 2](https://latex.codecogs.com/svg.image?k%20%3D%201%2C%202 "k = 1, 2"),
![\mathbf{z}\_{k,\mathbb{N}(\mathbf{s}\_0)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bz%7D_%7Bk%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D "\mathbf{z}_{k,\mathbb{N}(\mathbf{s}_0)}")
is a vector formed by stacking
![z_k(\mathbf{s}\_i)](https://latex.codecogs.com/svg.image?z_k%28%5Cmathbf%7Bs%7D_i%29 "z_k(\mathbf{s}_i)")
where
![\mathbf{s}\_i \in \mathbb{N}(\mathbf{s}\_0)](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_i%20%5Cin%20%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29 "\mathbf{s}_i \in \mathbb{N}(\mathbf{s}_0)"),
![\mathbf{C}\_{k,\mathbb{N}(\mathbf{s}\_0)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7BC%7D_%7Bk%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D "\mathbf{C}_{k,\mathbb{N}(\mathbf{s}_0)}")
is the correlation matrix obtained by evaluating the correlation
function
![C(\cdot, \ell_k)](https://latex.codecogs.com/svg.image?C%28%5Ccdot%2C%20%5Cell_k%29 "C(\cdot, \ell_k)")
at the distance matrix between locations in
![\mathbb{N}(\mathbf{s}\_0)](https://latex.codecogs.com/svg.image?%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29 "\mathbb{N}(\mathbf{s}_0)"),
and
![\mathbf{c}\_{k,\mathbf{s}\_0,\mathbb{N}(\mathbf{s}\_0)}](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bc%7D_%7Bk%2C%5Cmathbf%7Bs%7D_0%2C%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29%7D "\mathbf{c}_{k,\mathbf{s}_0,\mathbb{N}(\mathbf{s}_0)}")
is the vector of correlations with elements evaluated
![C(\cdot, \ell_k)](https://latex.codecogs.com/svg.image?C%28%5Ccdot%2C%20%5Cell_k%29 "C(\cdot, \ell_k)")
at the vector of distances from
![\mathbf{s}\_0](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_0 "\mathbf{s}_0")
to
![\mathbf{s}\_i \in \mathbb{N}(\mathbf{s}\_0)](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_i%20%5Cin%20%5Cmathbb%7BN%7D%28%5Cmathbf%7Bs%7D_0%29 "\mathbf{s}_i \in \mathbb{N}(\mathbf{s}_0)").

Sampling from
![\pi(z_k(\mathbf{s}\_0) \mid \mathbf{z}\_k, \sigma_k, \ell_k)](https://latex.codecogs.com/svg.image?%5Cpi%28z_k%28%5Cmathbf%7Bs%7D_0%29%20%5Cmid%20%5Cmathbf%7Bz%7D_k%2C%20%5Csigma_k%2C%20%5Cell_k%29 "\pi(z_k(\mathbf{s}_0) \mid \mathbf{z}_k, \sigma_k, \ell_k)")
for posterior predictive value for the latent processes
![z_k(\mathbf{s})](https://latex.codecogs.com/svg.image?z_k%28%5Cmathbf%7Bs%7D%29 "z_k(\mathbf{s})")
at prediction location
![\mathbf{s}\_0](https://latex.codecogs.com/svg.image?%5Cmathbf%7Bs%7D_0 "\mathbf{s}_0")
depend on how the latent processes are approximated. In what follows, we
describe the procedure under each approximating method.

## 1.7 Response NNGP in Stan

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

## 1.8 References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-banerjee2017high" class="csl-entry">

Banerjee, Sudipto. 2017. “High-Dimensional Bayesian Geostatistics.”
*Bayesian Analysis* 12 (2): 583–614.
<https://doi.org/10.1214/17-BA1056R>.

</div>

<div id="ref-datta2022nearest" class="csl-entry">

Datta, Abhirup. 2022. “Nearest-Neighbor Sparse Cholesky Matrices in
Spatial Statistics.” *Wiley Interdisciplinary Reviews: Computational
Statistics* 14 (5): e1574.

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
