% Gibbs Sampling in R, Haskell, Jags and Stan
% Dominic Steinitz
% 4th April 2014

---
bibliography: Bayes.bib
---

Introduction
============

It's possible to Gibbs sampling in most languages and since I am doing
some work in R and some work in Haskell, I thought I'd present a
simple example in both languages: estimating the mean from a normal
distribution with unknown mean and variance. Although one can do Gibbs
sampling directly in R, it is more common to use a specialised
language such as JAGS or STAN to do the actual sampling and do
pre-processing and post-processing in R. This blog post presents
implementations in native R, JAGS and STAN.

Preamble
--------

> {-# OPTIONS_GHC -Wall                      #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults    #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods  #-}
> {-# OPTIONS_GHC -fno-warn-orphans          #-}

> {-# LANGUAGE NoMonomorphismRestriction     #-}

> module ConjMCMCSimple where
>
> import qualified Data.Vector.Unboxed as V
> import Data.Random.Source.PureMT
> import Data.Random
> import Control.Monad.State
> import Data.Histogram ( asList )
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
>
> import LinRegAux

Normal Distribution with Unknown Mean and Variance
==================================================

It is fairly standard to use an improper prior

$$
\begin{aligned}
\pi(\mu, \tau) \propto \frac{1}{\tau} & & -\infty < \mu < \infty\, \textrm{and}\, 0 < \tau < \infty
\end{aligned}
$$

 is $\mu \sim \cal{N}(\mu_0, \sigma_0)$, that is

$$
\pi(\mu) \propto \exp{\bigg( -\frac{(\mu - \mu_0)^2}{2\sigma_0^2}\bigg)}
$$

Our data is IID normal, $x_i \sim \cal{N}(\mu, \sigma)$, where
$\sigma$ is known, so the likelihood is

$$
p(x\,|\,\mu, \sigma) \propto \prod_{i=1}^n \exp{\bigg( -\frac{(x_i - \mu)^2}{2\sigma^2}\bigg)}
$$

The assumption that $\sigma$ is known is unlikely but the point of
this post is to demonstrate MCMC matching an analytic formula.

This gives a posterior of

$$
\begin{aligned}
p(\mu\,|\, \boldsymbol{x}) &\propto \exp{\bigg(
-\frac{(\mu - \mu_0)^2}{2\sigma_0^2}
- \frac{\sum_{i=1}^n(x_i - \mu)^2}{2\sigma^2}\bigg)} \\
&\propto \exp{\bigg[-\frac{1}{2}\bigg(\frac{\mu^2 \sigma^2 -2\sigma^2\mu\mu_0 - 2\sigma_0^2n\bar{x}\mu + \sigma_0^2 n\mu^2}{\sigma^2\sigma_0^2}\bigg)\bigg]} \\
&= \exp{\bigg[-\frac{1}{2}\bigg(\frac{ (n\sigma_0^2 + \sigma^2)\mu^2 - 2(\sigma^2\mu_0 - \sigma_0^2n\bar{x})\mu}{\sigma^2\sigma_0^2}\bigg)\bigg]} \\
&= \exp{\Bigg[-\frac{1}{2}\Bigg(\frac{ \mu^2 - 2\mu\frac{(\sigma^2\mu_0 - \sigma_0^2n\bar{x})}{(n\sigma_0^2 + \sigma^2)}}{\frac{\sigma^2\sigma_0^2}{(n\sigma_0^2 + \sigma^2)}}\Bigg)\Bigg]} \\
&\propto \exp{\Bigg[-\frac{1}{2}\Bigg(\frac{\big(\mu - \frac{(\sigma^2\mu_0 - \sigma_0^2n\bar{x})}{(n\sigma_0^2 + \sigma^2)}\big)^2}{\frac{\sigma^2\sigma_0^2}{(n\sigma_0^2 + \sigma^2)}}\Bigg)\Bigg]}
\end{aligned}
$$

In other words

$$
\mu\,|\, \boldsymbol{x} \sim \cal{N}\bigg(\frac{\sigma^2\mu_0 + n\sigma_0^2\bar{x}}{n\sigma_0^2 + \sigma^2}, \frac{\sigma^2\sigma_0^2}{n\sigma_0^2 + \sigma^2} \bigg)
$$

Writing

$$
\sigma_n^2 = \frac{\sigma^2\sigma_0^2}{n\sigma_n^2 + \sigma^2}
$$

we get

$$
\frac{1}{\sigma_n^2} = \frac{n}{\sigma^2} + \frac{1}{\sigma_0^2}
$$

Thus the precision (the inverse of the variance) of the posterior is
the precision of the prior plus the precision of the data scaled by
the number of observations. This gives a nice illustration of how
Bayesian statistics improves our beliefs.

Writing

$$
\mu_n = \frac{\sigma^2\mu_0 + n\sigma_0^2\bar{x}}{n\sigma_0^2 + \sigma^2}
$$

and

$$
\lambda = 1 / \sigma^2, \, \lambda_0 = 1 / \sigma_0^2, \, \lambda_n = 1 / \sigma_n^2
$$

we see that

$$
\mu_n = \frac{n\bar{x}\lambda + \mu_0\lambda_0}{\lambda_n}
$$

Thus the mean of the posterior is a weight sum of the mean of the
prior and the sample mean scaled by preciscion of the prior and the
precision of the data itself scaled by the number of observations.

Rather arbitrarily let us pick a prior mean of

> mu0 :: Double
> mu0 = 11.0

and express our uncertainty about it with a largish prior variance

> sigma_0 :: Double
> sigma_0 = 2.0

And also arbitrarily let us pick the know variance for the samples as

> sigma :: Double
> sigma = 1.0

```{.dia height='600'}
import ConjMCMCSimple
import LinRegAux

dia = diagNormals [(mu0, sigma_0, blue, "Prior")]
````

We can sample from this in way that looks very similar to
[STAN](http://mc-stan.org) and
[JAGS](http://mcmc-jags.sourceforge.net):

> hierarchicalSample :: MonadRandom m => m Double
> hierarchicalSample = do
>   mu <- sample (Normal mu0 sigma_0)
>   x  <- sample (Normal mu sigma)
>   return x

and we didn't need to write a new language for this.

Again arbitrarily let us take

> nSamples :: Int
> nSamples = 10

and use

> arbSeed :: Int
> arbSeed = 2

And then actually generate the samples.

> simpleXs :: [Double]
> simpleXs =
>   evalState (replicateM nSamples hierarchicalSample)
>             (pureMT $ fromIntegral arbSeed)

Using the formulae we did above we can calculate the posterior

> mu_1, sigma1, simpleNumerator :: Double
> simpleNumerator = fromIntegral nSamples * sigma_0**2 + sigma**2
> mu_1 = (sigma**2 * mu0 + sigma_0**2 * sum simpleXs) / simpleNumerator
> sigma1 = sigma**2 * sigma_0**2 / simpleNumerator

and then compare it against the prior

```{.dia height='600'}
import ConjMCMCSimple
import LinRegAux

dia = diagNormals [(mu0, sigma_0, blue, "Prior"), (mu_1, sigma1, red, "Posterior")]
````

The red posterior shows we are a lot more certain now we have some evidence.

Via Markov Chain Monte Carlo
----------------------------

The theory behinde MCMC is described in a [previous
post](http://idontgetoutmuch.wordpress.com/2013/12/07/haskell-ising-markov-metropolis/). We
need to generate some proposed steps for the chain. We sample from the
normal distribution but we could have used e.g. the gamma.

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

Comparison
----------

Let's create the same histogram but from the posterior created analytically.

> analPosterior :: [Double]
> analPosterior =
>   evalState (replicateM (nIters - burnIn) (sample (Normal mu_1 (sqrt sigma1))))
>   (pureMT $ fromIntegral 5)
>
> histAnal :: Histogram V.Vector BinD Double
> histAnal = fillBuilder hb analPosterior

And then compare them. Because they overlap so well, we show the MCMC, both and the analytic on separate charts.

```{.dia width='800'}
dia = image "diagrams/HistMCMC.png" 1.0 1.0
      ===
      image "diagrams/HistMCMCAnal.png" 1.0 1.0
      ===
      image "diagrams/HistAnal.png" 1.0 1.0
````

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
> main = do
>   displayHeader "diagrams/HistMCMC.png"
>     (barDiag MCMC
>      (zip (map fst $ asList (hist 2)) (map snd $ asList (hist 2)))
>      (zip (map fst $ asList histAnal) (map snd $ asList histAnal)))
>
>   displayHeader "diagrams/HistMCMCAnal.png"
>     (barDiag MCMCAnal
>      (zip (map fst $ asList (hist 2)) (map snd $ asList (hist 2)))
>      (zip (map fst $ asList histAnal) (map snd $ asList histAnal)))
>
>   displayHeader "diagrams/HistAnal.png"
>     (barDiag Anal
>      (zip (map fst $ asList (hist 2)) (map snd $ asList (hist 2)))
>      (zip (map fst $ asList histAnal) (map snd $ asList histAnal)))
>
>   displayHeader "diagrams/SmoothHistMCMC.png"
>     (barDiag MCMC
>      (zip (map fst $ asList smoothHist) (map snd $ asList smoothHist))
>      (zip (map fst $ asList histAnal) (map snd $ asList histAnal)))

