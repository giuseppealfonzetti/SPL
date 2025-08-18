check_sa_args <- function(LIST, N, R, C) {
  out <- list()

  # pair per dimension per iteration
  if (is.null(LIST$PPI)) {
    LIST$PPI <- 16
  }
  stopifnot(is.numeric(LIST$PPI))
  stopifnot(LIST$PPI > 0)
  out$PAIRS_PER_ITERATION <- as.integer(LIST$PPI)

  # LR SCHEDULE
  if (is.null(LIST$SCHEDULE)) {
    LIST$SCHEDULE <- "PJ"
  }
  stopifnot(LIST$SCHEDULE %in% c("PJ", "Xu"))
  out$SCHEDULE <- switch(LIST$SCHEDULE, "PJ" = 1, "Xu" = 2)

  # UPDATE TYPE
  if (is.null(LIST$UPDATE)) {
    LIST$UPDATE <- "AdaptBurn"
  }
  stopifnot(LIST$UPDATE %in% c("AdaptBurn", "Plain"))
  out$UPDATE <- switch(LIST$UPDATE, "AdaptBurn" = 1, "Plain" = 0)

  # UPDATE SWITCH
  if (is.null(LIST$SWITCH)) {
    LIST$SWITCH <- TRUE
  }
  stopifnot(is.logical(LIST$SWITCH))
  out$SWITCH <- as.numeric(LIST$SWITCH)

  # BASE TOL
  if (is.null(LIST$TOL)) {
    LIST$TOL <- 1e-2
  }
  stopifnot(is.numeric(LIST$TOL))
  stopifnot(LIST$TOL > 0)
  out$TOL <- LIST$TOL

  # BASE LR
  if (is.null(LIST$STEP0)) {
    LIST$STEP0 <- 1e-2
  }
  stopifnot(is.numeric(LIST$STEP0))
  stopifnot(LIST$STEP0 > 0)
  out$STEP0 <- LIST$STEP0

  # Total number of cycles
  if (is.null(LIST$MAXE)) {
    LIST$MAXE <- 10
  }
  stopifnot(is.numeric(LIST$MAXE))
  stopifnot(LIST$MAXE > 0)
  stopifnot(LIST$MAXE > LIST$BURNE)
  out$MAXE <- as.integer(LIST$MAXE)

  # Number of cycles to burn before averaging
  if (is.null(LIST$BURNE)) {
    LIST$BURNE <- LIST$MAXE - 1
  }

  stopifnot(is.numeric(LIST$BURNE))
  stopifnot(LIST$BURNE > 0)
  out$BURNE <- as.integer(LIST$BURNE)

  # Update per Cycle
  if (is.null(LIST$UPE)) {
    LIST$UPE <- N / (out$MAXE)
  }
  stopifnot(is.numeric(LIST$UPE))
  stopifnot(LIST$UPE > 0)
  out$UPE <- as.integer(LIST$UPE)

  # Sampling SEED
  if (is.null(LIST$SEED)) {
    LIST$SEED <- 123
  }

  stopifnot(is.numeric(LIST$SEED))
  stopifnot(LIST$SEED > 0)
  out$SEED <- as.integer(LIST$SEED)

  # sampler
  if (is.null(LIST$REPLACEMENT)) {
    LIST$REPLACEMENT <- TRUE
  }
  stopifnot(is.logical(LIST$REPLACEMENT))

  out$REPLACEMENT <- LIST$REPLACEMENT

  # sampler
  if (is.null(LIST$CHECK)) {
    LIST$CHECK <- TRUE
  }
  stopifnot(is.logical(LIST$CHECK))
  out$CHECK <- LIST$CHECK

  # to be deprecated
  out$STEP1 <- 1
  out$STEP2 <- 1e-5
  out$STEP3 <- 3 / 4

  #
  return(out)
}
