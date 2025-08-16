#' Creat factor dictionaries
#'
#' @param F1 Row factor
#' @param F2 Column factor
#'
#' @return List with one dictionary per factor
#' @export
get_list_ind_dict <- function(F1, F2) {
  list1 <- split(1:length(F1), F1)
  len1 <- unlist(lapply(list1, length))
  which.list1 <- which(len1 > 1)
  list1 <- list1[which.list1]
  list2 <- split(1:length(F2), F2)
  len2 <- unlist(lapply(list2, length))
  which.list2 <- which(len2 > 1)
  list2 <- list2[which.list2]
  np1 <- sum(len1 * (len1 - 1) / 2)
  np2 <- sum(len2 * (len2 - 1) / 2)

  return(list(
    dict1 = cpp_get_dict(list1, np1),
    dict2 = cpp_get_dict(list2, np2)
  ))
}

#' From marginal to conditional parameterisation
#' @param PAR Marginal parameter vector
#' @export
marg2cond <- function(PAR) {
  p <- length(PAR) - 2
  gamma <- PAR[1:p]
  sigma2A <- PAR[p + 1] / (1 - PAR[p + 1] - PAR[p + 2])
  sigma2B <- PAR[p + 2] / (1 - PAR[p + 1] - PAR[p + 2])
  betaAB <- gamma * sqrt(1 + sigma2A + sigma2B)
  return(c(betaAB, sqrt(sigma2A), sqrt(sigma2B)))
}
