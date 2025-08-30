#' Fit via stochastic pairwise likelihood
#'
#' @param Y Response vector
#' @param X Design matrix
#' @param F1 Vector with row memberships
#' @param F2 Vector with column memberships
#' @param MODEL Choose among "probit", "logit" and "ordprobit"
#' @param START Starting vector for model parameters
#' @param CONTROL control list to be passed to the optimiser
#' @param VERBOSE Verbose output
#'
#' @export
saPL <- function(
  Y,
  X,
  F1,
  F2,
  MODEL = c("probit", "logit", "ordprobit"),
  START = NULL,
  CONTROL = list(),
  VERBOSE = 0
) {
  stopifnot(nrow(X) == length(Y))
  stopifnot(length(F1) == length(Y))
  stopifnot(length(F2) == length(Y))
  if (is.null(START)) {
    thr <- qnorm(
      cumsum(table(Y)) / sum(table(Y))
    )[1:4]
    fe <- rep(0, ncol(X))

    if (MODEL == "ordprobit") {
      START <- c(thr, fe, 1 / 3, 1 / 3)
    } else {
      START <- c(fe, 1 / 3, 1 / 3)
    }
  }

  MODEL <- match.arg(MODEL)
  list_dict <- get_list_ind_dict(F1, F2)
  if (VERBOSE > 0) {
    cat(
      "Pairs F1:",
      nrow(list_dict$dict1),
      ", Pairs F2:",
      nrow(list_dict$dict2),
      "\n"
    )
  }

  ctrl_args <- check_sa_args(
    CONTROL,
    N = length(Y),
    R = length(unique(F1)),
    C = length(unique(F2))
  )
  args <- c(
    # CONTROL,
    ctrl_args,
    list(
      START = START,
      Y = Y,
      X = X,
      LINK = MODEL,
      DICT1 = as.matrix(list_dict$dict1),
      DICT2 = as.matrix(list_dict$dict2),
      NCAT = length(unique(Y)),
      VERBOSE = VERBOSE
    )
  )
  fit <- do.call(cpp_SA2, args)

  fit$types <- c(rep("fixed", ncol(X)), rep("var", 2))
  fit$names <- c(colnames(X), "vF1", "vF2")
  if (MODEL == "ordprobit") {
    ncat <- length(unique(Y))
    th_names <- paste0(1:(ncat - 1), "|", 2:ncat)
    fit$names <- c(th_names, fit$names)
    fit$types <- c(rep("threshold", ncat - 1), fit$types)
  }
  return(fit)
}
