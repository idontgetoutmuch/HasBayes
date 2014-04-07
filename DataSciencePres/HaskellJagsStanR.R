## Import the library that allows R to inter-work with jags.
library(rjags)

## Read the simulated data into a data frame.
fn <- read.table("example1.data", header=FALSE)

jags <- jags.model('example1.bug',
                   data = list('x' = fn[,1], 'N' = 20),
                   n.chains = 4,
                   n.adapt = 100)

## Burnin for 10000 samples
update(jags, 10000);

mcmc_samples <- coda.samples(jags, variable.names=c("mu", "sigma"), n.iter=20000)
plot(mcmc_samples)


