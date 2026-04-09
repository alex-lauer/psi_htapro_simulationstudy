required_events <- function(
    HR,
    power = 0.80,
    alpha = 0.05,
    allocation = 0.5,
    sided = 2
) {
  # Z-values
  z_alpha <- if (sided == 2) {
    qnorm(1 - alpha / 2)
  } else {
    qnorm(1 - alpha)
  }
  z_beta <- qnorm(power)

  # Schoenfeld formula
  D <- (z_alpha + z_beta)^2 /
    ((log(HR))^2 * allocation * (1 - allocation))

  ceiling(D)
}

required_events(
  HR = 0.3,
  power = 0.8,
  alpha = 0.05,
  allocation = 0.5
)

# 22 total events are required to achieve 80% power
# (for a very strong effect like HR=0.3)
