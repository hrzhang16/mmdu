#' @title Model-Based Multidimensional Unfolding
#'
#' @description Model-Based Multidimensional Unfolding
#'
#' @param data the data matrix. Entries are binary. Missing entries should be \code{NA}.
#' @param K the dimension of space one wants to project the ideal points onto. Default value is 2.
#' @param delta the pre-specified positive intercept parameter, which determines
#' the probability when the distance is 0. Default value is 1e-1.
#' @param M the positive constant which constrains the norm of estimated ideal points. Default value is 3.
#' @param tol tolerance for termination. When the increment for loglikelihood between two
#' iterations is less than tol, the iteration is terminated. Default value is 1.
#' @param TA_init staring ideal points for optimization. If missing, the results from classical
#' multidimensional unfolding will be used as starting points.
#' @return The function returns a list with the following components:
#' \describe{
#'   \item{A}{The estimated item ideal points.}
#'   \item{Theta}{The estimated individual ideal points.}
#' }
#'
#' @export
#'
#' @import Rcpp RcppArmadillo smacof mvtnorm stats
#' @useDynLib mmdu, .registration = TRUE
#'
#' @examples
#' data <- vote108$data
#' tmp <- munfolding(data)
#' Theta <- tmp$Theta
#' A <- tmp$A
#'
#' \dontrun{
#' plot(A[, 1], A[, 2], pch = 2, xlim = c(-4, 4), ylim = c(-3, 3), ann = FALSE)
#' points(Theta[, 1], Theta[, 2], col = 6, pch = 1)
#' legend("bottomleft", c("senator", "roll call"), col = c(6, 1),pch = c(1, 2), inset = 0)
#' }
#'
#' @references Chen, Y., Ying, Z., & Zhang, H. (2020). Unfolding-Model-Based Visualization: Theory, Method and Applications.
#' \url{https://arxiv.org/abs/2009.01551}.
munfolding <- function(data,
                       K,
                       delta,
                       M,
                       tol,
                       TA_init) {

  # checking input data
  if(missing(data)) stop("\"data\" is missing")
  N <- nrow(data)
  J <- ncol(data)

  data.na <- is.na(data)
  Omega <- 1 - data.na

  uni <- unique(c(data))
  uni <- uni[which(is.na(uni) == FALSE)]
  uni <- sort(uni, decreasing = FALSE)

  num.cat <- length(uni)
  if (num.cat <= 1) stop("only 1 category in \"data\"")
  if (num.cat >= 3) stop("more than 2 categories in \"data\"")

  Y <- data
  Y[Omega == 0] <- 0
  Y[data == uni[1]] <- 0
  Y[data == uni[2]] <- 1

  # checking input K
  if (missing(K)) {
    if (missing(TA_init) == FALSE) {
      K <- ncol(TA_init)
    } else {
      K <- 2
    }
  }

  # checking input delta
  if (missing(delta)) {
    delta <- 0.1
  } else if (delta < 0) stop("\"delta\" is negative")

  # checking input M
  if (missing(M)) {
    M <- 3
  } else if (M <= 0.5) {
    warning("\"M\" is less than 0.5")
  } else if (M >= 5) warning("\"M\" is larger than 5")

  # checking input tol
  if (missing(tol)) {
    tol <- 1
  } else if (tol <= 0) {
    stop("\"tol\" is negative")
  }

  # initial values for ideal points
  Theta_init <- matrix(0, N, K)
  A_init <- matrix(0, J, K)

  if (missing(TA_init) == FALSE) {

    if(ncol(TA_init) != K) {
      stop("dimension of \"TA_init\" does not match the given \"K\"")
    }
    Theta_init <- TA_init[1 : N, ]
    A_init <- TA_init[(N + 1) : (N + J), ]

  } else {

    Y_dist <- 1 - Y
    tmp <- unfolding(Y_dist, ndim = K, type = "ordinal", weightmat = Omega)
    Theta_init <- tmp$conf.row
    A_init <- tmp$conf.col

  }

  res <- opti_fn(Theta_init, A_init, delta, Y, Omega, M, tol)
  res

}
