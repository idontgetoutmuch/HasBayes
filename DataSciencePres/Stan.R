library(rstan)

## Read the simulated data into a data frame.
fn <- read.table("example1.data", header=FALSE)

## Running the model
fit1 <- stan(file = 'Stan.stan',
             data = list('x' = fn[,1], 'N' = 20),
             pars=c("mu", "sigma"),
             chains=3,
             iter=30000,
             warmup=10000)

png(file="diagrams/stan.png",width=400,height=350)
plot(fit1)
dev.off()
