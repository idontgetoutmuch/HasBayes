library(moments)

nrep <- 105000
tau  <- array(,nrep)
mu   <- array(,nrep)

## Read the simulated data into a data frame.
fn <- read.table("example1.data", header=FALSE)

x      <- fn[,1]
n      <- length(x)
xBar   <- mean(x)
x2     <- sum(x^2)
lambda <- 1.0 / (0.5 * n * var(x))

tTau <- rgamma(1, shape = n / 2, rate = lambda)

for (i in 1:nrep) {
    tMu  <- rnorm(1, mean = xBar, sd = 1 / sqrt(n * tTau))
    scale  <- 0.5 * (x2  + n * tMu^2 - 2 * n * tMu * xBar)
    tTau <- rgamma(1, shape = n / 2, scale)
    tau[i] <- tTau
    mu[i]  <- tMu
}

nb  <- 5000
nb1 <- nb + 1

tauBurnt <- tau[nb1:nrep]
muBurnt  <- mu[nb1:nrep]

summary(muBurnt)

