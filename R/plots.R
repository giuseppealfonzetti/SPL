#' Plot stochastic parameters trajectories
#'
#' @param FIT Fitted object from [saPL]
#' @param MARG TRUE for marginal parameterisation. FALSE for conditional
#'
#' @importFrom dplyr tibble as_tibble mutate left_join filter
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes theme theme_minimal labs element_text
#' @importFrom ggplot2 geom_line geom_point facet_wrap
#' @importFrom ggplot2 scale_fill_gradient2 guide_colorbar guides unit
#' @importFrom patchwork plot_layout
#' @importFrom rlang .data
#'
#' @export
plot_sa_traj <- function(FIT, MARG = FALSE) {
  if (MARG) {
    saEstTraj <- Reduce(
      rbind,
      FIT$path_avtheta
    )
  } else {
    saEstTraj <- Reduce(
      rbind,
      lapply(FIT$path_avtheta, function(par) marg2cond(unlist(par)))
    )
  }

  colnames(saEstTraj) <- FIT$names

  dictEst <- tibble(
    par = FIT$names,
    type = FIT$types
  )

  gg1 <- saEstTraj |>
    as_tibble() |>
    mutate(iter = FIT$path_iters) |>
    pivot_longer(-.data$iter, names_to = "par", values_to = "est") |>
    left_join(dictEst, by = "par") |>
    filter(.data$type == "threshold") |>
    ggplot(aes(
      x = .data$iter,
      y = .data$est,
      group = .data$par,
      col = .data$par
    )) +
    theme_minimal() +
    geom_line() +
    geom_point(size = 1) +
    labs(
      title = "Thresholds",
      y = "Stochastic estimates",
      x = "Iteration",
      col = "Category"
    ) +
    theme(
      # legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  gg2 <- saEstTraj |>
    as_tibble() |>
    mutate(iter = FIT$path_iters) |>
    pivot_longer(-.data$iter, names_to = "par", values_to = "est") |>
    left_join(dictEst, by = "par") |>
    filter(.data$type == "fixed") |>
    ggplot(aes(x = .data$iter, y = .data$est, group = .data$par)) +
    theme_minimal() +
    geom_line() +
    geom_point(size = 1) +
    labs(
      title = "Fixed effects",
      y = "Stochastic estimates",
      x = "Iteration"
    ) +
    facet_wrap(~ .data$par, scales = "free") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  gg3 <- saEstTraj |>
    as_tibble() |>
    mutate(iter = FIT$path_iters) |>
    pivot_longer(-.data$iter, names_to = "par", values_to = "est") |>
    left_join(dictEst, by = "par") |>
    filter(.data$type == "var") |>
    ggplot(aes(
      x = .data$iter,
      y = .data$est,
      group = .data$par,
      col = .data$par
    )) +
    theme_minimal() +
    geom_line() +
    geom_point(size = 1) +
    labs(
      title = "Variance components",
      y = "Stochastic estimates",
      x = "Iteration",
      col = "Group"
    ) +
    theme(
      # legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  (gg2 | (gg1 / gg3)) + plot_layout(ncol = 2, widths = c(4, 1))
}
