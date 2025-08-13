objFUN <- function(
  PAR,
  Y,
  X,
  DICT1,
  DICT2,
  LINK
) {
  FUN_logPL(
    Y = Y,
    X = X,
    DICT1 = DICT1,
    DICT2 = DICT2,
    LINK = LINK,
    THETA = PAR
  )
}

grFUN <- function(
  PAR,
  Y,
  X,
  DICT1,
  DICT2,
  LINK
) {
  FUN_grlogPL(
    Y = Y,
    X = X,
    DICT1 = DICT1,
    DICT2 = DICT2,
    LINK = LINK,
    THETA = PAR
  )
}


test_that("link probit", {
  skip_if_not_installed("glmm")
  skip_if_not_installed("numDeriv")

  data(salamander, package = "glmm")

  link <- "probit"
  mod.glm <- glm(Mate ~ Cross - 1, family = binomial(link), data = salamander)
  list_dict <- get_list_ind_dict(salamander$Male, salamander$Female)
  init <- c(mod.glm$coefficients, 0.5, 0.5)

  expect_equal(
    numDeriv::grad(
      func = objFUN,
      x = init,
      Y = salamander$Mate,
      X = model.matrix(mod.glm),
      DICT1 = as.matrix(list_dict$dict1),
      DICT2 = as.matrix(list_dict$dict2),
      LINK = link
    ),
    grFUN(
      PAR = init,
      Y = salamander$Mate,
      X = model.matrix(mod.glm),
      DICT1 = as.matrix(list_dict$dict1),
      DICT2 = as.matrix(list_dict$dict2),
      LINK = link
    ),
    tolerance = 1e-3
  )
})

test_that("link logit", {
  skip_if_not_installed("glmm")
  skip_if_not_installed("numDeriv")

  data(salamander, package = "glmm")
  link <- "logit"
  mod.glm <- glm(Mate ~ Cross - 1, family = binomial(link), data = salamander)
  list_dict <- get_list_ind_dict(salamander$Male, salamander$Female)
  init <- c(mod.glm$coefficients, 0.5, 0.5)

  expect_equal(
    numDeriv::grad(
      func = objFUN,
      x = init,
      Y = salamander$Mate,
      X = model.matrix(mod.glm),
      DICT1 = as.matrix(list_dict$dict1),
      DICT2 = as.matrix(list_dict$dict2),
      LINK = link
    ),
    grFUN(
      PAR = init,
      Y = salamander$Mate,
      X = model.matrix(mod.glm),
      DICT1 = as.matrix(list_dict$dict1),
      DICT2 = as.matrix(list_dict$dict2),
      LINK = link
    ),
    tolerance = 1e-3
  )
})


test_that("link ordprobit", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("numDeriv")

  data(InstEval, package = "lme4")
  InstRid <- InstEval[1:500, ]
  InstRid$s <- droplevels(InstRid$s)
  InstRid$d <- droplevels(InstRid$d)
  link <- "ordprobit"
  mod.polr <- MASS::polr(
    factor(y) ~ studage + lectage + service + dept,
    method = "probit",
    data = InstRid
  )

  list_dict <- get_list_ind_dict(InstRid$s, InstRid$d)
  init <- c(mod.polr$zeta, mod.polr$coefficients, 0.25, 0.25) ###polychoric correlations of 0.5 were too off

  expect_equal(
    numDeriv::grad(
      func = objFUN,
      x = init,
      Y = InstRid$y * 1.0,
      X = model.matrix(mod.polr)[, -1],
      DICT1 = as.matrix(list_dict$dict1),
      DICT2 = as.matrix(list_dict$dict2),
      LINK = link
    ),
    grFUN(
      PAR = init,
      Y = InstRid$y * 1.0,
      X = model.matrix(mod.polr)[, -1],
      DICT1 = as.matrix(list_dict$dict1),
      DICT2 = as.matrix(list_dict$dict2),
      LINK = link
    ),
    tolerance = 1e-3
  )
})

# marg2cond <- function(theta) {
#   p <- length(theta) - 2
#   gamma <- theta[1:p]
#   sigma2A <- theta[p + 1] / (1 - theta[p + 1] - theta[p + 2])
#   sigma2B <- theta[p + 2] / (1 - theta[p + 1] - theta[p + 2])
#   betaAB <- gamma * sqrt(1 + sigma2A + sigma2B)
#   return(c(betaAB, sqrt(sigma2A), sqrt(sigma2B)))
# }

# init <- c(mod.glm$coefficients, 0.5, 0.5)
# ple <- optimx::optimx(
#   init,
#   FUN_logPL,
#   gr = FUN_grlogPL,
#   X = x,
#   Y = y,
#   DICT1 = as.matrix(list_dict$dict1),
#   DICT2 = as.matrix(list_dict$dict2),
#   LINK = link,
#   control = list(maximize = TRUE),
#   method = "nlm"
# )
# marg2cond(unlist(ple))
# fit.arc <- arcProbit::arcbin.fit(
#   x,
#   y,
#   f1 = salamander$Male,
#   f2 = salamander$Female,
#   obj_glm = mod.glm
# )

# ple <- fullPL(
#   Y = y,
#   X = x,
#   F1 = salamander$Male,
#   F2 = salamander$Female,
#   START = init
# )
# ple
