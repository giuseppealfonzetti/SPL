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

# FUN_logPL <- function(THETA, X, Y, LINK, DICT1, DICT2) {
#   cpp_llikFullPool2D(
#     Y = Y,
#     X = X,
#     DICT1 = DICT1,
#     DICT2 = DICT2,
#     THETA = THETA,
#     LINK = LINK,
#     NCAT = length(unique(Y))
#   )
# }

# FUN_grlogPL <- function(THETA, X, Y, LINK, DICT1, DICT2) {
#   cpp_grllFullPool2D(
#     Y = Y,
#     X = X,
#     DICT1 = DICT1,
#     DICT2 = DICT2,
#     THETA = THETA,
#     LINK = LINK,
#     NCAT = length(unique(Y))
#   )
# }
