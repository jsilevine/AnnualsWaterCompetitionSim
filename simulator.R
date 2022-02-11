## this is a simulator

library(rootSolve)
library(ggplot2)

## define constants
nu <- 1.5       ## allometric constant
k <- 0.05       ## xylem conductivity
a_m <- 0.3      ## maximum rate of carbon accumulation
Gamma <- 50     ##
gamma <- 0.5    ## root shoot ratio
rr <- 0.1       ## root respiration rate
rl <- 0.2       ## leaf respiration rate
Vmax <- 1.0545  ## light availability
C_a <- 400      ## atmospheric carbon concentration
omega <- 450    ## empirical constant
D_s <- 1.5       ## vapor pressure deficit
D_0 <- 1.5       ## empirical constant
m <- 5.6        ## empirical constant

## calculate soil water content, W, from soil water potential (see App. Eq. 7)
calc_W <- function(psi_s, Wmax = 0.6, Wmin = 0.01, lambda = 4, psi_max = -0.05) {

  Wscale <- Wmax - Wmin

  W <- ( ( ( psi_max / psi_s ) ^ (1/lambda) ) * Wscale ) + Wmin

  return(W)

}

## calculate soil water potential, psi, from soil water content, W (see App. Eq. 7)
calc_psi <- function(W, Wmax = 0.6, Wmin = 0.01, lambda = 4, psi_max = -0.05) {

  Wscale <- Wmax - Wmin

  psi <- psi_max * (Wscale / (W - Wmin))^lambda

  return(psi)

}


## calculate phi (collection of constants, see equation 2 from man.) from parameters
calc_phi <- function() {

  phi <- m / ((C_a - Gamma) * (1 + (D_s / D_0)))
  return(phi)

}

## calculate a_m, the maximum carbon accumulation rate (minus leaf respiration) per unit leaf area
calc_a_m <- function(Vmax = 1.054545) {

  phi <- calc_phi()

  c_i <- C_a - (1 / phi)

  a_m <- ((Vmax * c_i) / (omega + c_i)) - rl

  return(a_m)

}

## function to calculate a
calc_a <- function(psi_s, psi_0, Vmax = 1.054545) {

  afun <- function(a) {

    ((a * (Vmax - rl - a)) / (((gamma * k) / D_s) * ((C_a * (Vmax - rl - a)) -
                             (omega * (a + rl))))) - (psi_s - psi_0)

  }

  ## solve for a
  a <- min(uniroot.all(afun, c(0, 1), tol = 1e-5))

  return(a)

}


## function to calculate G
calc_G <- function(a, delta = 1) {

  G <- ((nu - 1)/nu) * ((a - gamma * rr)/(delta))
  return(G)

}

## calculate psi_star, the value of psi_soil at which psi_leaf first reaches the critical value (see App. Eq. 5)
calc_psi_star <- function(psi_0 = -1.5, Vmax = 1.054545) {

  phi <- calc_phi()

  a_m <- calc_a_m(Vmax = Vmax)

  psi_star <- (( phi * a_m * D_s) / (gamma * k)) + psi_0

  return(psi_star)

}

calc_dtau <- function(t, dt) {

  ((t + dt)^(nu/(nu-1))) - (t^(nu/(nu-1)))

}

## function to simulate single year
simulate_year <- function(W0, dt = 0.25, psi_0, Vmax, N, evap = 0) {

  ## calculate initial soil water potential
  psi_init <- calc_psi(W0)

  Bdata <- data.frame(t = 0, psi_s = psi_init, W = W0)
  B <- as.data.frame(t(as.matrix(rep(0, times = length(psi_0)), nrow = 1)))
  colnames(B) <- as.character(seq(1, length(B)))
  Bdata <- cbind(Bdata, B)

  ## generate data frame of species characteristics
  Sdata <- data.frame(spp = as.character(seq(1, length(N))),
                      psi_0 = psi_0,
                      N = N,
                      Vmax = Vmax,
                      a_m = sapply(Vmax, FUN = calc_a_m),
                      psi_star = mapply(psi_0, Vmax, FUN = calc_psi_star))

  ## indicator variable
  i <- 1

  ## continue looping until latest phenology species is out of water
  while (Bdata[i, "psi_s"] > min(Sdata$psi_0)) {

    a <- rep(0, times = nrow(Sdata))

    if (any(Bdata[i, "psi_s"] < Sdata$psi_0)) {
      a[Bdata[i, "psi_s"] < Sdata$psi_0] <- 0
    }

    ## calculate vector of carbon accumulations
    if (any(Bdata[i, "psi_s"] < Sdata$psi_star & Bdata[i, "psi_s"] > Sdata$psi_0)) {
      a[Bdata[i, "psi_s"] < Sdata$psi_star & Bdata[i, "psi_s"] > Sdata$psi_0] <- mapply(Sdata[Bdata[i, "psi_s"] < Sdata$psi_star & Bdata[i, "psi_s"] > Sdata$psi_0, "psi_0"],
                                                                                        Sdata[Bdata[i, "psi_s"] < Sdata$psi_star & Bdata[i, "psi_s"] > Sdata$psi_0, "Vmax"],
                                                                                        FUN = calc_a, psi_s = Bdata[i, "psi_s"])
    }

    ## if psi_s > psi_star, a = a_m
    if (any(Bdata[i, "psi_s"] > Sdata$psi_star)) {
      a[Bdata[i, "psi_s"] > Sdata$psi_star] <- Sdata[Bdata[i, "psi_s"] > Sdata$psi_star, "a_m"]
    }

    ## truncate at 0 to prevent negative growth
    a[a <= 0] <- 0

    ## calculate growth rates
    G <- sapply(a, FUN = calc_G)
    G[G <= 0] <- 0

    ## impose some cutoff for to reduce simulation time, makes almost no difference.
    ## This is biologically reasonable given Beta shutoff operator
    if (all(G <= 1e-2)) break

    G_B <- G^(nu/(nu-1))
    G_L <- G^(1/(nu-1))

    dtau <- (Bdata[i, "t"] + dt)^(nu/(nu-1)) - Bdata[i, "t"]^(nu/(nu-1))

    nB <- Bdata[i, as.character(seq(1, nrow(Sdata)))] + (G_B * dtau)
    colnames(nB) <- as.character(seq(1, nrow(Sdata)))

    ## transpire water
    E <- 1e-2 * ((nu-1)/nu) * sum((Sdata[, "N"] * G_L)) * dtau

    nW <- Bdata[i, "W"] - E

    nW <- nW - (evap * nW)

    npsi_s <- calc_psi(nW)

    nBd <- data.frame(t = Bdata[nrow(Bdata), "t"] + dt,
                         W = nW,
                         psi_s = npsi_s)

    nBd <- cbind(nBd, nB)

    Bdata <- rbind(Bdata, nBd)
    rownames(Bdata) <- NULL

   i <- i + 1

  }

  return(list(Bdata, Bdata[i,]))

}

simulate_dynamics <- function(W0, dt = 0.01, Tmax = 100, psi_0, Vmax, N0, F, evap = 0) {

  ## set up data
  data <- data.frame(T = 0)
  data <- cbind(data, as.data.frame(matrix(N0, ncol = length(N0))))
  colnames(data) <- c("T", as.character(seq(1, length(N0))))

  equil <- FALSE
  i <- 1
  while (equil == FALSE) {

    out <- simulate_year(W0 = W0, dt = dt, psi_0 = psi_0,
                         Vmax = Vmax, N = as.numeric(data[i, 2:(length(N0)+1)]),
                         evap = evap)[[2]]

    ndata <- as.data.frame(c(i, out[,4:(3+length(N0))] * as.numeric(data[i, 2:(length(N0)+1)]) * F))
    colnames(ndata) <- c("T", as.character(seq(1, length(N0))))

    data <- rbind(data, ndata)

    deltas <- abs((data[i+1, 2:(length(N0)+1)] - data[i, 2:(length(N0)+1)]) / data[i, 2:(length(N0)+1)])

    if (all(deltas < 0.0001) | i > Tmax) equil <- TRUE
    print(i)
    print(deltas)

    i <- i + 1

  }

  return(data)

}
