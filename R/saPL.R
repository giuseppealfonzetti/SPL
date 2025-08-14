#' @export
saPL <- function(
  Y,
  X,
  F1,
  F2,
  MODEL = c("probit", "logit", "ordprobit"),
  START,
  CONTROL,
  VERBOSE = 0
) {
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

  args <- c(
    CONTROL,
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
  return(fit)
}
