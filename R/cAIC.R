cAIC <-
function(object, method = NULL, B = NULL) {
  # A function that calls the bias correction functions.
  #
  # Args: 
  #   object = Object of class lmerMod or glmerMod. Obtained by lmer() or glmer()
  #   method = How the bias correction should be evaluated. If NULL than method 
  #            is chosen by family, i.e. analytical if family is Poisson or 
  #            Gaussian and with parametric bootstrap for other. Method may also
  #            be specified before, either "steinian" or "conditionalBootstrap".
  #            "steinian" only available for Gaussian, Poisson and Bernoulli.
  #   B      = Number of Bootstrap replications. Default is NULL then it is 
  #            chosen as maximum of the number of observations and 100.
  #
  # Returns:
  #   list   = The list contains the conditional log-likelihood; the estimated 
  #            conditional prediction error (degrees of freedom); If a new model
  #            was fitted, the new model and an boolean indicating weather a new
  #            model was fitted; the conditional Akaike information, caic.
  #
  if (any(names(object) == "mer")) {
    warning("Treat with care: gamm4 models are not sufficiently tested!") 
    object <- object$mer
  }
  
  if (!inherits(object, c("lmerMod", "glmerMod"))) {
    stop("Class of object is not known")
  }
  
  if (family(object)$family == "binomial" && length(unique(getME(object, "y"))) > 2) {
    warning("Method not yet supplied for binomial data with n larger 2. 
            Therefore the conditional parametric bootstrap is returned")
    method <- "conditionalBootstrap"
  }
  
  dfList   <- bcMer(object , method = method, B = B)
  if (mode(dfList) == "list") {
    bc       <- dfList$bc
    newModel <- dfList$newModel
    new      <- dfList$new
  } else {
    bc       <- dfList
    newModel <- NULL
    new      <- FALSE
  }
  
  cll  <- getcondLL(object)
  caic <- - 2 * cll + 2 * bc
  return(list(loglikelihood = cll, 
              df            = bc, 
              newModel      = newModel, 
              new           = new, 
              caic          = caic))
}
