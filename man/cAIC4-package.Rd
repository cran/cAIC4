\name{cAIC4-package}
\alias{cAIC4-package}
\docType{package}

\title{
Conditional Akaike information criterion for lme4
}

\description{
Provides functions for the estimation of the conditional Akaike 
information in Generalized Mixed-effects models fitted with (g)lmer 
form lme4.
}

\details{
\tabular{ll}{
Package: \tab cAIC\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-05-23\cr
License: \tab GPL (>=2)\cr
}
}

\author{
Benjamin Saefken, David Ruegamer

based on the work of Thomas Kneib and Sonja Greven.

Maintainer: Benjamin Saefken <bsaefke@uni-goettingen.de>
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


\keyword{ package }

\examples{
b <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

cAIC(b)
}
