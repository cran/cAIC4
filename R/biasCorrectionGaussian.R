biasCorrectionGaussian <-
function(m) {
  zeroLessModel <- deleteZeroComponents(m)
  if (inherits(zeroLessModel, "lm")) {
    return(ncol(model.matrix(zeroLessModel)))
  }
  model <- getModelComponents(zeroLessModel)
  if (identical(m, zeroLessModel)) {
    bc       <- calculateGaussianBc(model)
    newModel <- NULL
    new      <- FALSE
  } else {
    bc       <- calculateGaussianBc(model)
    newModel <- zeroLessModel
    new      <- FALSE
  }
  return(list(bc = bc, newModel = newModel, new = new))
}
