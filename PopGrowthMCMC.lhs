% Bayesian Analysis: A Conjugate Prior and Markov Chain Monte Carlo
% Dominic Steinitz
% 9th March 2014

---
bibliography: Bayes.bib
---

Introduction
============




Preamble
--------

> {-# OPTIONS_GHC -Wall                      #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults    #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods  #-}
> {-# OPTIONS_GHC -fno-warn-orphans          #-}

> {-# LANGUAGE NoMonomorphismRestriction     #-}

> module PopGrowthMCMC where
>
> import qualified Data.Vector.Unboxed as V
> import Data.Random.Source.PureMT
> import Data.Random
> import Control.Monad.State
> import qualified Control.Monad.Writer as W
> import Control.Monad.Loops
> -- import Data.Histogram ( asList )
> import qualified Data.Histogram as H
> import Data.Histogram.Fill
> import Data.Histogram.Generic ( Histogram )
> import Data.List
> import Control.Parallel.Strategies
>
> import Diagrams.Backend.Cairo.CmdLine
>
> import Diagrams.Backend.CmdLine
> import Diagrams.Prelude hiding ( sample, render )


A Simple Example
================

$$
\begin{aligned}
\dot{p} & =  rp\Big(1 - \frac{p}{k}\Big)
\end{aligned}
$$

$$
p = \frac{kp_0\exp rt}{k + p_0(\exp rt - 1)}
$$


Let us see if we can estimate the parameter for population growth
using MCMC in the
[example](https://idontgetoutmuch.wordpress.com/2014/09/09/fun-with-extended-kalman-filters-4/)
in which we used Kalman filtering.

Our data is given by

$$
\begin{aligned}
p_i &= \frac{kp_0\exp r\Delta T i}{k + p_0(\exp r\Delta T i - 1)} \\
y_i &= p_i + \epsilon_i
\end{aligned}
$$

In other words $y_i \sim \cal{N}(p_i, \sigma^2)$, where $\sigma$ is
known so the likelihood is

$$
p(y\,|\,r) \propto \prod_{i=1}^n \exp{\bigg( -\frac{(y_i - p_i)^2}{2\sigma^2}\bigg)} =
\exp{\bigg( -\sum_{i=1}^n \frac{(y_i - p_i)^2}{2\sigma^2}\bigg)}
$$

Let us assume a prior of $r \sim {\cal{N}}(\mu_0,\sigma_0^2)$ then the posterior becomes

$$
p(r\,|\,y) \propto \exp{\bigg( -\frac{(r - \mu_0)^2}{2\sigma_0^2} \bigg)} \exp{\bigg( -\sum_{i=1}^n \frac{(y_i - p_i)^2}{2\sigma^2}\bigg)}
$$

Let us assume a growth rate

> mu0 :: Double
> mu0 = 10.0

The known variance for the samples

> sigma0 :: Double
> sigma0 = 1e1

> sigma :: Double
> sigma = 1e-2

> prior' :: Double -> Double
> prior' r = exp (-(r - mu0)**2 / (2 * sigma0**2))
>
> likelihood' :: Double -> [Double] -> Double
> likelihood' r ys = exp (-sum (zipWith (\y mu -> (y - mu)**2 / (2 * sigma**2)) ys mus))
>   where
>     mus :: [Double]
>     mus = map (logit p0 k . (* (r * deltaT))) (map fromIntegral [0..])
>
> posterior' :: Double -> [Double] -> Double
> posterior' r ys = likelihood' r ys * prior' r

Here's the implementation of the logistic function

> logit :: Double -> Double -> Double -> Double
> logit p0 k x = k * p0 * (exp x) / (k + p0 * (exp x - 1))

We assume most of the parameters are known with the exception of the
the growth rate $r$. We fix this also in order to generate test data.

> k, p0 :: Double
> k = 1.0
> p0 = 0.1

> r, deltaT :: Double
> r = 10.0
> deltaT = 0.0005

> singleSample :: Double -> RVarT (W.Writer [Double]) Double
> singleSample p0 = do
>   epsilon <- rvarT (Normal 0.0 sigma)
>   let p1 = logit p0 k (r * deltaT)
>   lift $ W.tell [p1 + epsilon]
>   return p1

> nObs :: Int
> nObs = 300

> streamSample :: RVarT (W.Writer [Double]) Double
> streamSample = iterateM_ singleSample p0

> samples :: [Double]
> samples = take nObs $ snd $
>           W.runWriter (evalStateT (sample streamSample) (pureMT 3))

The [Metropolis
algorithm](http://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm)
tells us that we always jump to a better place but only sometimes jump
to a worse place. We count the number of acceptances as we go.

> acceptanceProb' :: Double -> Double -> [Double] -> Double
> acceptanceProb' r r' ys = min 1.0 ((posterior' r' ys) / (posterior' r ys))

> oneStep' :: (Double, Int) -> (Double, Double) -> (Double, Int)
> oneStep' (mu, nAccs) (proposedJump, acceptOrReject) =
>   if acceptOrReject < acceptanceProb' r (r + proposedJump) samples
>   then (mu + proposedJump, nAccs + 1)
>   else (mu, nAccs)



and express our uncertainty about it with a largish prior variance

And also arbitrarily let us pick the known variance for the samples as

Via Markov Chain Monte Carlo
----------------------------

> normalisedProposals :: Int -> Double -> Int -> [Double]
> normalisedProposals seed sigma nIters =
>   evalState (replicateM nIters (sample (Normal 0.0 sigma)))
>   (pureMT $ fromIntegral seed)

We also need samples from the uniform distribution

> acceptOrRejects :: Int -> Int -> [Double]
> acceptOrRejects seed nIters =
>   evalState (replicateM nIters (sample stdUniform))
>   (pureMT $ fromIntegral seed)

And now we can calculate the (un-normalised) prior, likelihood and posterior

The [Metropolis
algorithm](http://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm)
tells us that we always jump to a better place but only sometimes jump
to a worse place. We count the number of acceptances as we go.

Now we can actually run our simulation. We set the number of jumps and
a burn in but do not do any thinning.

> nIters, burnIn :: Int
> nIters = 3000
> burnIn = nIters `div` 10

Let us start our chain at

> startMu :: Double
> startMu = 10.0

and set the variance of the jumps to

> jumpVar :: Double
> jumpVar = 0.01

> test :: Int -> [(Double, Int)]
> test seed =
>   drop burnIn $
>   scanl oneStep' (startMu, 0) $
>   zip (normalisedProposals seed jumpVar nIters)
>       (acceptOrRejects seed nIters)

We put the data into a histogram

> numBins :: Int
> numBins = 400

> hb :: HBuilder Double (Data.Histogram.Generic.Histogram V.Vector BinD Double)
> hb = forceDouble -<< mkSimple (binD lower numBins upper)
>   where
>     lower = startMu - 1.5*sigma0
>     upper = startMu + 1.5*sigma0
>
> hist :: Int -> Histogram V.Vector BinD Double
> hist seed = fillBuilder hb (map fst $ test seed)

```{.dia width='800'}
dia = image "diagrams/HistMCMC.png" 1.0 1.0
````

Not bad but a bit lumpy. Let's try a few runs and see if we can smooth
things out.

> hists :: [Histogram V.Vector BinD Double]
> hists = parMap rpar hist [3,4..102]

> emptyHist :: Histogram V.Vector BinD Double
> emptyHist = fillBuilder hb (replicate numBins 0)
>
> smoothHist :: Histogram V.Vector BinD Double
> smoothHist = foldl' (H.zip (+)) emptyHist hists

```{.dia width='800'}
dia = image "diagrams/SmoothHistMCMC.png" 1.0 1.0
````

Quite nice and had my machine running at 750% with +RTS -N8.

