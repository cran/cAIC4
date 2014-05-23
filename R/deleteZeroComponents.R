deleteZeroComponents <-
function(m) {
  # A function that deletes all random effects terms if corresponding variance
  # parameter is estimated to zero.
  #
  # Args: 
  #   m     = Object of class lmerMod. Obtained by lmer()
  #
  # Returns:
  #   m/newMod = A model without zero estimated variance component
  #
  theta <- getME(m, "theta")
  thetazero  <- which(theta == 0)
  
  if (length(thetazero) == 0) {  # every thing's fine
    return(m)
  }
    
  if (length(theta) == length(thetazero)) {  #  only lm left
    warning("Model has no random effects variance components larger than zero.")
    return(lm(nobars(formula(m)), model.frame(m)))
  }

  varBlockMatrices <- getME(m, "ST")
  cnms <- m@cnms
  
  for(i in 1:length(varBlockMatrices)){
    cnms[[i]] <- cnms[[i]][which(diag(varBlockMatrices[[i]]) != 0)]
  }
  
  reFormula  <- cnms2formula(cnms)
  rhs        <- c(attr(terms(nobars(formula(m))), "term.labels"), reFormula)
  lhs        <- formula(m)[[2]]  # left hand side of the formlua
  newFormula <- reformulate(rhs, lhs)  # merge both sides           
  newMod     <- update(m, formula. = newFormula, evaluate = TRUE)
  
  return(deleteZeroComponents(newMod))
}
