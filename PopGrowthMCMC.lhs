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

Our data is IID normal, $x_i \sim \cal{N}(\mu, \sigma)$, where
$\sigma$ is known, so the likelihood is

$$
p(x\,|\,\mu, \sigma) \propto \prod_{i=1}^n \exp{\bigg( -\frac{(x_i - \mu)^2}{2\sigma^2}\bigg)}
$$

where

$$
\mu = 
$$

$$
\eta = \log{\frac{\theta}{1 - \theta}}
$$

> mu0 :: Double
> mu0 = 11.0

and express our uncertainty about it with a largish prior variance

> sigma_0 :: Double
> sigma_0 = 2.0

And also arbitrarily let us pick the know variance for the samples as

> sigma :: Double
> sigma = 1.0

> hierarchicalSample :: MonadRandom m => m Double
> hierarchicalSample = do
>   mu <- sample (Normal mu0 sigma_0)
>   x  <- sample (Normal mu sigma)
>   return x

> nSamples :: Int
> nSamples = 10

> arbSeed :: Int
> arbSeed = 2

> simpleXs :: [Double]
> simpleXs =
>   evalState (replicateM nSamples hierarchicalSample)
>             (pureMT $ fromIntegral arbSeed)

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

> prior :: Double -> Double
> prior mu = exp (-(mu - mu0)**2 / (2 * sigma_0**2))
>
> likelihood :: Double -> [Double] -> Double
> likelihood mu xs = exp (-sum (map (\x -> (x - mu)**2 / (2 * sigma**2)) xs))
>
> posterior :: Double -> [Double] -> Double
> posterior mu xs = likelihood mu xs * prior mu

The [Metropolis
algorithm](http://en.wikipedia.org/wiki/Metropolis–Hastings_algorithm)
tells us that we always jump to a better place but only sometimes jump
to a worse place. We count the number of acceptances as we go.

> acceptanceProb :: Double -> Double -> [Double] -> Double
> acceptanceProb mu mu' xs = min 1.0 ((posterior mu' xs) / (posterior mu xs))

> oneStep :: (Double, Int) -> (Double, Double) -> (Double, Int)
> oneStep (mu, nAccs) (proposedJump, acceptOrReject) =
>   if acceptOrReject < acceptanceProb mu (mu + proposedJump) simpleXs
>   then (mu + proposedJump, nAccs + 1)
>   else (mu, nAccs)

Now we can actually run our simulation. We set the number of jumps and
a burn in but do not do any thinning.

> nIters, burnIn :: Int
> nIters = 300000
> burnIn = nIters `div` 10

Let us start our chain at

> startMu :: Double
> startMu = 10.0

and set the variance of the jumps to

> jumpVar :: Double
> jumpVar = 0.4

> test :: Int -> [(Double, Int)]
> test seed =
>   drop burnIn $
>   scanl oneStep (startMu, 0) $
>   zip (normalisedProposals seed jumpVar nIters)
>       (acceptOrRejects seed nIters)

We put the data into a histogram

> numBins :: Int
> numBins = 400

> hb :: HBuilder Double (Data.Histogram.Generic.Histogram V.Vector BinD Double)
> hb = forceDouble -<< mkSimple (binD lower numBins upper)
>   where
>     lower = startMu - 1.5*sigma_0
>     upper = startMu + 1.5*sigma_0
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

PostAmble
=========

Normally with BlogLiteratelyD, we can generate diagrams on the
fly. However, here we want to run the simulations in parallel so we
need to actually compile something.

~~~~ { .shell }
ghc -O2 ConjMCMCSimple.lhs -main-is ConjMCMCSimple -threaded -fforce-recomp
~~~~

> displayHeader :: FilePath -> Diagram B R2 -> IO ()
> displayHeader fn =
>   mainRender ( DiagramOpts (Just 900) (Just 700) fn
>              , DiagramLoopOpts False Nothing 0
>              )

> main :: IO ()
> main = undefined


