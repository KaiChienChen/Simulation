library(ggplot2)
library(gridExtra)
library(dplyr)
library(mcmcse)

initialize_lattice <- function(L) {
  matrix(sample(c(-1, 1), L * L, replace = TRUE), L, L)
}

compute_neighbors_sum <- function(lattice, i, j) {
  L <- nrow(lattice)
  lattice[(i %% L) + 1, j] +
    lattice[(i - 2) %% L + 1, j] +
    lattice[i, (j %% L) + 1] +
    lattice[i, (j - 2) %% L + 1]
}

gibbs_step <- function(lattice, beta) {
  L <- nrow(lattice)
  for (i in 1:L) {
    for (j in 1:L) {
      neighbors_sum <- compute_neighbors_sum(lattice, i, j)
      prob <- 1 / (1 + exp(-2 * beta * neighbors_sum))
      lattice[i, j] <- ifelse(runif(1) < prob, 1, -1)
    }
  }
  return(lattice)
}

metropolis_step <- function(lattice, beta) {
  L <- nrow(lattice)
  for (k in 1:(L * L)) {
    i <- sample(1:L, 1)
    j <- sample(1:L, 1)
    spin <- lattice[i, j]
    neighbors_sum <- compute_neighbors_sum(lattice, i, j)
    delta_E <- 2 * spin * neighbors_sum
    if (delta_E <= 0 || runif(1) < exp(-beta * delta_E)) {
      lattice[i, j] <- -spin
    }
  }
  return(lattice)
}

monitor_convergence <- function(magnetizations) {
  mcse <- mcse(magnetizations)$est
  mean_mag <- mean(magnetizations)
  data.frame(Mean_Magnetization = mean_mag, MCSE = mcse)
}

L <- 32
steps <- 100
temperatures <- c(1.5, 2.5, 3.5)
results <- list()

par(mfrow = c(length(temperatures), 2))
for (T in temperatures) {
  beta <- 1 / T
  lattice_gibbs <- initialize_lattice(L)
  lattice_mh <- initialize_lattice(L)

  for (step in 1:steps) {
    lattice_gibbs <- gibbs_step(lattice_gibbs, beta)
    lattice_mh <- metropolis_step(lattice_mh, beta)
  }

  image(t(apply(lattice_gibbs, 2, rev)), col = c("blue", "red"), main = paste("Gibbs Sampling: T =", T))
  image(t(apply(lattice_mh, 2, rev)), col = c("blue", "red"), main = paste("Metropolis-Hastings: T =", T))
}

for (T in temperatures) {
  lattice <- initialize_lattice(L)
  magnetizations <- numeric(steps)
  
  for (step in 1:steps) {
    lattice <- gibbs_step(lattice, 1 / T)
    magnetizations[step] <- mean(lattice)
  }
  
  conv_metrics <- monitor_convergence(magnetizations)
  print(conv_metrics)
}
