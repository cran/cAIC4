calculateGaussianBc <-
function(model) {
  # A function that calculates the analytic representation of the bias 
  # corrections in linear mixed models, see Greven & Kneib (2010).
  #
  # Args: 
  #   model = From getAllModelComponents()
  #
  # Returns:
  #   df = Bias correction (i.e. degrees of freedom) for a linear mixed model.
  #           
  C     <- model$C
  B     <- model$B
  e     <- model$e
  A     <- model$A
  tye   <- model$tye
  V0inv <- model$V0inv
    for (j in 1:length(model$theta)) {
        Wj     <- model$Wlist[[j]]
        eWje   <- model$eWelist[[j]]
        C[j, ] <- as.vector((e %*% Wj) %*% A - eWje * e/(2 * tye))
        for (k in j:length(model$theta)) {
            Wk <- model$Wlist[[k]]
            eWke   <- model$eWelist[[k]]
            if (!model$isREML) {
              B[j, k] <- B[k, j] <-  - tye * 
                sum(t(Wk %*% V0inv) * (Wj %*% V0inv))/(2 * model$n) - 
                eWje * eWke/(2 * tye) + 
                as.numeric(e %*% Wk %*% (A %*% (Wj %*% e)))
            } else {
              B[j, k] <- B[k, j] <- - tye * 
                sum(t(Wk %*% A) * (Wj %*% A))/(2*(model$n - ncol(model$X))) - 
                eWje * eWke/(2 * tye) + 
                as.numeric(e %*% Wk %*% (model$A %*% (Wj %*% e)))
            }
        }
    }
# Naive solver:
#  Lambday <- solve(B, C)

#if (solver != "chol") { # SVD solver:
#  USVt    <- svd(B)
#  UtC     <- crossprod(USVt$u, C)
#  LambdaV <- backsolve(diag(USVt$d), UtC)
#  Lambday <- crossprod(USVt$v, LambdaV)
#} else {  # Cholesky solver:
  Rchol   <- chol(B)
# Rchol   <- chol(B, pivot = TRUE)
# Rchol   <- Rchol[, order(attr(Rchol, "pivot"))]
# Rchol   <- qr.R(qr(Rchol, tol = 1e-20))
  L1      <- backsolve(Rchol, C, transpose = TRUE)
  Lambday <- backsolve(Rchol, L1)


    df <- model$n - sum(diag(A))
    for (j in 1:length(model$theta)) {
        df <- df + sum(Lambday[j,] %*% (A %*% (model$Wlist[[j]] %*% e)))  
    }
return(df + 1)
}
