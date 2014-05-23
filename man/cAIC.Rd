\name{cAIC}
\alias{cAIC}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Conditional Akaike information for lmer & glmer
}

\description{
Estimates the conditional Akaike information for models that were fitted in 
lme4. This is possible for all distributions, i.e. \code{\link[stats]{family}} 
arguments, based on parametric conditional bootstrap. For the Gaussian 
distribution (from a \code{\link[lme4]{lmer}} call) and the Poisson distribution 
analytical estimators for the degrees of freedom are available, based on Stein 
type formulas. Also the conditional Akaike information for generalized additive 
models based on a fit via the \code{gamm4}-package can be estimated. 
}

\usage{
cAIC(object, method = NULL, B = NULL)
}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{
An object of class merMod either fitted by \code{\link[lme4]{lmer}} or 
\code{\link[lme4]{glmer}} of the lme4-package. Also objects returned form a 
\code{gamm4} call are possible.
}
  \item{method}{
Either \code{"conditionalBootstrap"} for the estimation of the degrees of 
freedom with the help of conditional Bootstrap or \code{"steinian"} for 
analytical representations based on Stein type formulas. The default is 
\code{NULL}. In this case the method is choosen automatically based on the 
\code{family} argument of the \code{(g)lmer}-object. For \code{"gaussian"} and 
\code{"poisson"} this is the Steinian type estimator, for all others it is the 
conditional Bootstrap. 
}
  \item{B}{
Number of Bootstrap replications. The default is \code{NULL}. Then B is the 
minimum of 100 and the length of the response vector. 
}
}

\details{
For \code{method = "steinian"} and an object of class \code{merMod} computed the 
analytic representation of the corrected conditional AIC in Greven and Kneib 
(2010). This is based on a the Stein formula and uses implicit differentiation 
to calculate the derivative of the random effects covariance parameters w.r.t. 
the data. The code is adapted form the one provided in the supplementary 
material of the paper by Greven and Kneib (2010). Although the inverse of the 
matrix B from Theorem 3 in that paper is not explicitly calculated. Instead a 
Cholesky decomposition of the (actually positive definite) matrix B is used, in 
order to increase stability and speed. Furthermore the model needs to be checked 
if a random effects covariance parameter has an optimum on the boundary. And if 
so the model needs to be refitted without the these parameters. This is also 
done by the function and the refitted model is also returned. Notice that the 
\code{boundary.tol} argument in \code{\link[lme4]{lmerControl}} has an impact on 
whether a parameter is estimated to lie on the boundary of the parameter space. 

If the object is of class \code{\link[lme4]{merMod}} and has \code{family = 
"poisson"} there is also an analytic representation of the conditional AIC based 
on the Chen-Stein formula, see for instance Saefken et. al (2014). For the 
calculation the model needs to be refitted for each observed response variable 
minus the number of response variables that are exactly zero. The calculation 
therefore takes longer then for models with Gaussian responses. Due to the speed 
and stability of lme4 this is still possible, also for larger datasets. 

If the model has Bernoulli distributed responses and \code{method = "steinian"}, 
\code{\link{cAIC}} calculates the degrees of freedom based on a proposed 
estimator by Efron (2004). This estimator is asymptotically unbiased if the 
estimated conditional mean is consistent. The calculation needs as many model 
refits as there are data points. 

Another more general method for the estimation of the degrees of freedom is the 
conditional bootstrap. This is proposed in Efron (2004). For the B boostrap 
samples the degrees of freedom are estimated by 
\deqn{\frac{1}{B - 1}\sum_{i=1}^n\theta_i(z_i)(z_i-\bar{z}),}
where \eqn{\theta_i(z_i)} is the i-th element of the estimated natural 
parameter.
}

\value{
A list consisting of: 1. the conditional log likelihood, i.e. the log likelihood 
with the random effects as penalized parameters; 2. the estimated degrees of 
freedom; 3. a list element that is either \code{NULL} if no new model was fitted 
otherwise the new (reduced) model, see details; 4. a boolean variable indicating 
whether a new model was fitted or not; 5. the estimator of the conditional 
Akaike information, i.e. minus twice the log likelihood plus twice the degrees 
of freedom.
}

\references{
Saefken, B., Kneib T., van Waveren C.-S. and Greven, S. (2014) A unifying 
approach to the estimation of the conditional Akaike information in generalized 
linear mixed models. Electronic Journal Statistics Vol. 8, 201-225.

Greven, S. and Kneib T. (2010) On the behaviour of marginal and conditional AIC 
in linear mixed models. Biometrika 97(4), 773-789.

Efron , B. (2004) The estimation of prediction error. J. Amer. Statist. Ass. 
99(467), 619-632. 
}

\author{
Benjamin Saefken \email{bsaefke@uni-goettingen.de}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\section{WARNINGS }{
Currently the cAIC can only be estimated for \code{family} equal to 
\code{"gaussian"}, \code{"poisson"} and \code{"binomial"}. Neither negative 
binomial nor gamma distributed responses are available. 

Weighted Gaussian models are not yet implemented.
}

\seealso{
\code{\link[lme4]{lme4-package}}, \code{\link[lme4]{lmer}}, 
\code{\link[lme4]{glmer}} 
}

\examples{
b <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

cAIC(b)

b2 <- lmer(Reaction ~ (1 | Days) + (1 | Subject), sleepstudy)

cAIC(b2)

b2ML <- lmer(Reaction ~ (1 + Days | Subject), sleepstudy, REML = FALSE)

cAIC(b2ML)
}

\keyword{regression}
