set.seed(1729)
y <- rnorm(n = 20, mean = 10, sd = 5)
mean(y)
sd(y)

write.table(y,
            file = 'example1.data',
            row.names = FALSE,
            col.names = FALSE)

jags <- jags.model('example1.bug',
                   data = "example1.data",
                   n.chains = 4,
                   n.adapt = 100)

library(rjags)

fn <- read.table("example1.data", header=FALSE)

jags <- jags.model('example1.bug',
                   data = list('x' = fn[,1], 'N' = 20),
                   n.chains = 4,
                   n.adapt = 100)

# The model specification
model_string <- "model{
  for(i in 1:length(y)) {
    y[i] ~ dnorm(mu, tau)
  }
  mu ~ dnorm(0, 0.0001)
  sigma ~ dlnorm(0, 0.0625)
  tau <- 1 / pow(sigma, 2)
}"

# Running the model
model <- jags.model(textConnection(model_string), data = list(y = y), n.chains = 3, n.adapt= 10000)
update(model, 10000); # Burnin for 10000 samples
mcmc_samples <- coda.samples(model, variable.names=c("mu", "sigma"), n.iter=20000)
