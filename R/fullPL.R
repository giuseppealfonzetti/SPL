FUN_logPL <- function(THETA, X, Y, LINK, DICT1, DICT2) {
  cpp_llikFullPool2D(
    Y = Y,
    X = X,
    DICT1 = DICT1,
    DICT2 = DICT2,
    THETA = THETA,
    LINK = LINK,
    NCAT = length(unique(Y))
  )
}

FUN_grlogPL <- function(THETA, X, Y, LINK, DICT1, DICT2) {
  cpp_grllFullPool2D(
    Y = Y,
    X = X,
    DICT1 = DICT1,
    DICT2 = DICT2,
    THETA = THETA,
    LINK = LINK,
    NCAT = length(unique(Y))
  )
}

#' Compute full Pairwise likelihood
#'
#' @param Y Response vector
#' @param X Design matrix
#' @param F1 Row factor vector
#' @param F2 Column factor vector
#' @param MODEL Choose among "probit", "logit" and "ordprobit"
#' @param START Starting vector for model parameters
#'
#' @return fitted model
#'
#' @export
fullPL <- function(
  Y,
  X,
  F1,
  F2,
  MODEL = c("probit", "logit", "ordprobit"),
  START
) {
  MODEL <- match.arg(MODEL)
  list_dict <- get_list_ind_dict(F1, F2)
  fit <- optimx::optimx(
    START,
    FUN_logPL,
    gr = FUN_grlogPL,
    X = X,
    Y = Y * 1.0,
    DICT1 = as.matrix(list_dict$dict1),
    DICT2 = as.matrix(list_dict$dict2),
    LINK = MODEL,
    control = list(maximize = TRUE),
    method = "nlm"
  )

  out <- list()
  out$fit <- fit
  return(out)
}
