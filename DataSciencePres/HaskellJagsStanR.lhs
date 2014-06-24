% Gibbs Sampling in R, Haskell, Jags and Stan
% Dominic Steinitz
% 4th April 2014

---
bibliography: ../Bayes.bib
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
implementations in native R, JAGS and STAN as well as Haskell.

Preamble
--------

> {-# OPTIONS_GHC -Wall                      #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults    #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods  #-}
> {-# OPTIONS_GHC -fno-warn-orphans          #-}

> {-# LANGUAGE NoMonomorphismRestriction     #-}

> module Gibbs (
>     main
>   , m
>   , Moments(..)
>   ) where
>
> import qualified Data.Vector.Unboxed as V
> import qualified Control.Monad.Loops as ML
> import Data.Random.Source.PureMT
> import Data.Random
> import Control.Monad.State
> import Data.Histogram ( asList )
> import Data.Histogram.Fill
> import Data.Histogram.Generic ( Histogram )
> import Data.List
> import qualified Control.Foldl as L
>
> import Diagrams.Backend.Cairo.CmdLine
>
> import LinRegAux
>
> import Diagrams.Backend.CmdLine
> import Diagrams.Prelude hiding ( sample, render )

The length of our chain and the burn-in.

> nrep, nb :: Int
> nb   = 5000
> nrep = 105000

Data generated from ${\cal{N}}(10.0, 5.0)$.

> xs :: [Double]
> xs = [
>     11.0765808082301
>   , 10.918739177542
>   , 15.4302462747137
>   , 10.1435649220266
>   , 15.2112705014697
>   , 10.441327659703
>   , 2.95784054883142
>   , 10.2761068139607
>   , 9.64347295100318
>   , 11.8043359297675
>   , 10.9419989262713
>   , 7.21905367667346
>   , 10.4339807638017
>   , 6.79485294803006
>   , 11.817248658832
>   , 6.6126710570584
>   , 12.6640920214508
>   , 8.36604701073303
>   , 12.6048485320333
>   , 8.43143879537592
>   ]


A Bit of Theory
===============

Gibbs Sampling
--------------

For a multi-parameter situation, Gibbs sampling is a special case of
Metropolis-Hastings in which the proposal distributions are the
posterior conditional distributions.

Referring back to the explanation of the [metropolis
algorithm](http://idontgetoutmuch.wordpress.com/2013/12/07/haskell-ising-markov-metropolis/),
let us describe the state by its parameters $i \triangleq
\boldsymbol{\theta}^{(i)} \triangleq (\theta^{(i)}_1,\ldots,
\theta^{(i)}_n)$ and the conditional posteriors by
$\pi\big({\theta}_{k}^{(j)} \,\big|\, {\boldsymbol{\theta}}_{-k}^{(i)}\big)$ where ${\boldsymbol{\theta}}^{(i)}_{-k} = \big(\theta_1^{(i)},\ldots,\theta_{k-1}^{(i)},\theta_{k+1}^{(i)}\ldots\theta_n^{(i)}\big)$
then

$$
\begin{aligned}
\frac{\pi\big(\boldsymbol{\theta}^{(j)}\big)q\big(\boldsymbol{\theta}^{(j)}, \boldsymbol{\theta}^{(i)}\big)}
{\pi(\boldsymbol{\theta}^{(i)})q(\boldsymbol{\theta}^{(i)}, \boldsymbol{\theta}^{(j)})}
&=
\frac{
\pi\big({\theta}_{k}^{(j)} \,\big|\, {\boldsymbol{\theta}}_{-k}^{(j)}\big)\pi\big({\boldsymbol{\theta}}_{-k}^{(j)}\big)\pi\big({\theta}_{k}^{(i)} \,\big|\, {\boldsymbol{\theta}}_{-k}^{(j)}\big)
}
{
\pi\big({\theta}_{k}^{(i)} \,\big|\, {\boldsymbol{\theta}}_{-k}^{(i)}\big)\pi\big({\boldsymbol{\theta}}_{-k}^{(i)}\big)\pi\big({\theta}_{k}^{(j)} \,\big|\, {\boldsymbol{\theta}}_{-k}^{(i)}\big)
} \\
&=
\frac{
\pi\big({\theta}_{k}^{(j)} \,\big|\, {\boldsymbol{\theta}}_{-k}^{(j)}\big)\pi\big({\boldsymbol{\theta}}_{-k}^{(j)}\big)\pi\big({\theta}_{k}^{(i)} \,\big|\, {\boldsymbol{\theta}}_{-k}^{(j)}\big)
}
{
\pi\big({\theta}_{k}^{(i)} \,\big|\, {\boldsymbol{\theta}}_{-k}^{(j)}\big)\pi\big({\boldsymbol{\theta}}_{-k}^{(j)}\big)\pi\big({\theta}_{k}^{(j)} \,\big|\, {\boldsymbol{\theta}}_{-k}^{(j)}\big)
} \\
&= 1
\end{aligned}
$$

where we have used the rules of conditional probability and the fact that $\boldsymbol{\theta}_i^{(-k)} = \boldsymbol{\theta}_j^{(-k)}$

Thus we always accept the proposed jump. Note that the chain is not in
general reversible as the order in which the updates are done matters.

Normal Distribution with Unknown Mean and Variance
--------------------------------------------------

It is fairly standard to use an improper prior

$$
\begin{aligned}
\pi(\mu, \tau) \propto \frac{1}{\tau} & & -\infty < \mu < \infty\, \textrm{and}\, 0 < \tau < \infty
\end{aligned}
$$

The likelihood is

$$
p(\boldsymbol{x}\,|\,\mu, \sigma) = \prod_{i=1}^n \bigg(\frac{1}{\sigma\sqrt{2\pi}}\bigg)\exp{\bigg( -\frac{(x_i - \mu)^2}{2\sigma^2}\bigg)}
$$

re-writing in terms of precision

$$
p(\boldsymbol{x}\,|\,\mu, \tau) \propto \prod_{i=1}^n \sqrt{\tau}\exp{\bigg( -\frac{\tau}{2}{(x_i - \mu)^2}\bigg)} = \tau^{n/2}\exp{\bigg( -\frac{\tau}{2}\sum_{i=1}^n{(x_i - \mu)^2}\bigg)}
$$

Thus the posterior is

$$
p(\mu, \tau \,|\, \boldsymbol{x}) \propto \tau^{n/2 - 1}\exp{\bigg( -\frac{\tau}{2}\sum_{i=1}^n{(x_i - \mu)^2}\bigg)}
$$

We can re-write the sum in terms of the sample mean $\bar{x} =
\frac{1}{n}\sum_{i=1}^n x_i$ and variance $s^2 =
\frac{1}{n-1}\sum_{i=1}^n (x_i - \bar{x})^2$ using

$$
\begin{aligned}
\sum_{i=1}^n (x_i - \mu)^2 &= \sum_{i=1}^n (x_i - \bar{x} + \bar{x} - \mu)^2 \\
&= \sum_{i=1}^n (x_i - \bar{x})^2 - 2\sum_{i=1}^n (x_i - \bar{x})(\bar{x} - \mu) + \sum_{i=1}^n (\bar{x} - \mu)^2 \\
&= \sum_{i=1}^n (x_i - \bar{x})^2 - 2(\bar{x} - \mu)\sum_{i=1}^n (x_i - \bar{x}) + \sum_{i=1}^n (\bar{x} - \mu)^2 \\
&= (n - 1)s^2 + n(\bar{x} - \mu)^2
\end{aligned}
$$

Thus the conditional posterior for $\mu$ is

$$
\begin{aligned}
p(\mu \,|\, \tau, \boldsymbol{x}) &\propto \exp{\bigg( -\frac{\tau}{2}\bigg(\nu s^2 + \sum_{i=1}^n{(\mu - \bar{x})^2}\bigg)\bigg)} \\
&\propto \exp{\bigg( -\frac{n\tau}{2}{(\mu - \bar{x})^2}\bigg)} \\
\end{aligned}
$$

which we recognise as a normal distribution with mean of $\bar{x}$ and
a variance of $(n\tau)^{-1}$.

The conditional posterior for $\tau$ is

$$
\begin{aligned}
p(\tau \,|\, , \mu, \boldsymbol{x}) &\propto \tau^{n/2 -1}\exp\bigg(-\tau\frac{1}{2}\sum_{i=1}^n{(x_i - \mu)^2}\bigg)
\end{aligned}
$$

which we recognise as a gamma distribution with a shape of $n/2$ and a scale of $\frac{1}{2}\sum_{i=1}^n{(x_i - \mu)^2}$

In this particular case, we can calculate the marginal posterior of
$\mu$ analytically. Writing $z = \frac{\tau}{2}\sum_{i=1}^n{(x_i -
\mu)^2}$ we have

$$
\begin{aligned}
p(\mu \,|\, \boldsymbol{x}) &= \int_0^\infty p(\mu, \tau \,|\, \boldsymbol{x}) \textrm{d}\tau \\
&\propto \int_0^\infty \tau^{n/2 - 1}\exp{\bigg( -\frac{\tau}{2}\sum_{i=1}^n{(x_i - \mu)^2}\bigg)} \textrm{d}\tau \\
&\propto \bigg( \sum_{i=1}^n{(x_i - \mu)^2} \bigg)^{-n/2} \int_0^\infty z^{n/2 - 1}\exp{-z}\textrm{d}\tau \\
&\propto \bigg( \sum_{i=1}^n{(x_i - \mu)^2} \bigg)^{-n/2} \\
\end{aligned}
$$

Finally we can calculate

$$
\begin{aligned}
p(\mu \,|\, \boldsymbol{x}) &\propto \bigg( (n - 1)s^2 + n(\bar{x} - \mu)^2 \bigg)^{-n/2} \\
&\propto \bigg( 1 + \frac{n(\mu - \bar{x})^2}{(n - 1)s^2} \bigg)^{-n/2} \\
\end{aligned}
$$

This is the
[non-standardized Student's t-distribution](http://en.wikipedia.org/wiki/Student%27s_t-distribution#Non-standardized_Student.%27s_t-distribution) $t_{n-1}(\bar{x}, s^2/n)$.

Alternatively the marginal posterior of $\mu$ is

$$
\frac{\mu - \bar{x}}{s/\sqrt{n}}\bigg|\, x \sim t_{n-1}
$$

where $t_{n-1}$ is the standard t distribution with $n - 1$ degrees of freedom.

The Model in Haskell
====================

Following up on a
[comment](http://idontgetoutmuch.wordpress.com/2014/04/02/students-t-and-space-leaks/#comments)
from a previous blog post, let us try using the
[foldl](http://hackage.haskell.org/package/foldl) package to calculate
the length, the sum and the sum of squares traversing the list only
once. An improvement on creating your own strict record and using
*foldl'* but maybe it is not suitable for some methods
e.g. calculating the skewness and kurtosis incrementally, see below.

> x2Sum, xSum, n :: Double
> (x2Sum, xSum, n) = L.fold stats xs
>   where
>     stats = (,,) <$>
>             (L.premap (\x -> x * x) L.sum) <*>
>             L.sum <*>
>             L.genericLength

And re-writing the sample variance $s^2 = \frac{1}{n-1}\sum_{i=1}^n
(x_i - \bar{x})^2$ using

$$
\begin{aligned}
\sum_{i=1}^n (x_i - \bar{x})^2 &= \sum_{i=1}^n (x_i^2 - 2x_i\bar{x} + \bar{x}^2) \\
&= \sum_{i=1}^n x_i^2 - 2\bar{x}\sum_{i=1}^n x_i + \sum_{i=1}^n \bar{x}^2 \\
&= \sum_{i=1}^n x_i^2 - 2n\bar{x}^2 + n\bar{x}^2 \\
&= \sum_{i=1}^n x_i^2 - n\bar{x}^2 \\
\end{aligned}
$$

we can then calculate the sample mean and variance using the sums we
have just calculated.

> xBar, varX :: Double
> xBar = xSum / n
> varX = n * (m2Xs - xBar * xBar) / (n - 1)
>   where m2Xs = x2Sum / n

In random-fu, the Gamma distribution is specified by the rate paratmeter, $\beta$.

> beta, initTau :: Double
> beta = 0.5 * n * varX
> initTau = evalState (sample (Gamma (n / 2) beta)) (pureMT 1)

Our sampler takes an old value of $\tau$ and creates new values of
$\mu$ and $\tau$.

> gibbsSampler :: MonadRandom m => Double -> m (Maybe ((Double, Double), Double))
> gibbsSampler oldTau = do
>   newMu <- sample (Normal xBar (recip (sqrt (n * oldTau))))
>   let shape = 0.5 * n
>       scale = 0.5 * (x2Sum + n * newMu^2 - 2 * n * newMu * xBar)
>   newTau <- sample (Gamma shape (recip scale))
>   return $ Just ((newMu, newTau), newTau)

From which we can create an infinite stream of samples.

> gibbsSamples :: [(Double, Double)]
> gibbsSamples = evalState (ML.unfoldrM gibbsSampler initTau) (pureMT 1)

As our chains might be very long, we calculate the mean, variance,
skewness and kurtosis using an [incremental
method](http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance).

> data Moments = Moments { mN :: !Double
>                        , m1 :: !Double
>                        , m2 :: !Double
>                        , m3 :: !Double
>                        , m4 :: !Double
>                        }
>   deriving Show

> moments :: [Double] -> Moments
> moments xs = foldl' f (Moments 0.0 0.0 0.0 0.0 0.0) xs
>   where
>     f :: Moments -> Double -> Moments
>     f m x = Moments n' m1' m2' m3' m4'
>       where
>         n = mN m
>         n'  = n + 1
>         delta = x - (m1 m)
>         delta_n = delta / n'
>         delta_n2 = delta_n * delta_n
>         term1 = delta * delta_n * n
>         m1' = m1 m + delta_n
>         m4' = m4 m +
>               term1 * delta_n2 * (n'*n' - 3*n' + 3) +
>               6 * delta_n2 * m2 m - 4 * delta_n * m3 m
>         m3' = m3 m + term1 * delta_n * (n' - 2) - 3 * delta_n * m2 m
>         m2' = m2 m + term1

In order to examine the posterior, we create a histogram.

> numBins :: Int
> numBins = 400

> hb :: HBuilder Double (Data.Histogram.Generic.Histogram V.Vector BinD Double)
> hb = forceDouble -<< mkSimple (binD lower numBins upper)
>   where
>     lower = xBar - 2.0 * sqrt varX
>     upper = xBar + 2.0 * sqrt varX

And fill it with the specified number of samples preceeded by a burn-in.

> hist :: Histogram V.Vector BinD Double
> hist = fillBuilder hb (take (nrep - nb) $ drop nb $ map fst gibbsSamples)

Now we can plot this.

```{.dia width='800'}
dia = image "diagrams/DataScienceHaskPost.png" 1.0 1.0
````

And calculate the skewness and kurtosis.

> m :: Moments
> m = moments (take (nrep - nb) $ drop nb $ map fst gibbsSamples)

    [ghci]
    import Gibbs
    putStrLn $ show $ (sqrt (mN m)) * (m3 m) / (m2 m)**1.5
    putStrLn $ show $ (mN m) * (m4 m) / (m2 m)**2

We expect a skewness of 0 and a kurtosis of $3 + 6 / \nu - 4 = 3.4$
for $\nu = 19$. Not too bad.

The Model in JAGS
=================

[JAGS](http://mcmc-jags.sourceforge.net) is a mature, declarative,
domain specific language for building Bayesian statistical models
using Gibbs sampling.

Here is our model as expressed in JAGS. Somewhat terse.

~~~~ {.r include="example1.bug"}
~~~~

To run it and examine its results, we wrap it up in some R

~~~~{.r include="HaskellJagsStanR.R"}
~~~~

And now we can look at the posterior for $\mu$.

```{.dia width='800'}
dia = image "diagrams/jags.png" 1.0 1.0
````

The Model in STAN
=================

[STAN](http://mc-stan.org) is a domain specific language for building
Bayesian statistical models similar to JAGS but newer and which allows
variables to be re-assigned and so cannot really be described as
declarative.

~~~~ {.r include="Stan.stan"}
~~~~

Just as with JAGS, to run it and examine its results, we wrap it up in
some R.

~~~~{.r include="Stan.R"}
~~~~

Again we can look at the posterior although we only seem to get
medians and 80% intervals.

```{.dia width='800'}
dia = image "diagrams/stan.png" 1.0 1.0
````

PostAmble
=========

Write the histogram produced by the Haskell code to a file.

> displayHeader :: FilePath -> Diagram B R2 -> IO ()
> displayHeader fn =
>   mainRender ( DiagramOpts (Just 900) (Just 700) fn
>              , DiagramLoopOpts False Nothing 0
>              )

> main :: IO ()
> main = do
>   displayHeader "diagrams/DataScienceHaskPost.png"
>     (barDiag
>      (zip (map fst $ asList hist) (map snd $ asList hist)))

The code can be downloaded from [github](https://github.com/idontgetoutmuch/HasBayes).