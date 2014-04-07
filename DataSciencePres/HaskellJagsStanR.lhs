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

> module Gibbs where
>
> import qualified Data.Vector.Unboxed as V
> import qualified Control.Monad.Loops as ML
> import Data.Random.Source.PureMT
> import Data.Random
> import Control.Monad.State
> import Data.Histogram ( asList )
> import qualified Data.Histogram as H
> import Data.Histogram.Fill
> import Data.Histogram.Generic ( Histogram )
> import Data.List
> import qualified Control.Foldl as L
> import Control.Parallel.Strategies
>
> import Diagrams.Backend.Cairo.CmdLine
>
> import LinRegAux
>
> import Diagrams.Backend.CmdLine
> import Diagrams.Prelude hiding ( sample, render )

Normal Distribution with Unknown Mean and Variance
==================================================

It is fairly standard to use an improper prior

$$
\begin{aligned}
\pi(\mu, \tau) \propto \frac{1}{\tau} & & -\infty < \mu < \infty\, \textrm{and}\, 0 < \tau < \infty
\end{aligned}
$$

The likelihood is

$$
p(x\,|\,\mu, \sigma) = \prod_{i=1}^n \bigg(\frac{1}{\sigma\sqrt{2\pi}}\bigg)\exp{\bigg( -\frac{(x_i - \mu)^2}{2\sigma^2}\bigg)}
$$

re-writing in terms of precision

$$
p(x\,|\,\mu, \tau) \propto \prod_{i=1}^n \sqrt{\tau}\exp{\bigg( -\frac{\tau}{2}{(x_i - \mu)^2}\bigg)} = \tau^{n/2}\exp{\bigg( -\frac{\tau}{2}\sum_{i=1}^n{(x_i - \mu)^2}\bigg)}
$$

Thus the posterior is

$$
p(\mu, \tau \,|\, x) \propto \tau^{n/2 - 1}\exp{\bigg( -\frac{\tau}{2}\sum_{i=1}^n{(x_i - \mu)^2}\bigg)}
$$

The marginal posterior for $\mu$ is

$$
\begin{aligned}
p(\mu \,|\, , \tau, x) &\propto \exp{\bigg( -\frac{\tau}{2}\bigg(\nu s^2 + \sum_{i=1}^n{(\mu - \bar{x})^2}\bigg)\bigg)} \\
&\propto \exp{\bigg( -\frac{n\tau}{2}{(\mu - \bar{x})^2}\bigg)} \\
\end{aligned}
$$

which we recognise as a normal distribution with mean of $\bar{x}$ and
a variance of $(n\tau)^{-1}$.

The marginal posterior for $\tau$ is

$$
\begin{aligned}
p(\tau \,|\, , \mu, x) &\propto \tau^{n/2 -1}\exp\bigg(-\tau\frac{1}{2}\sum_{i=1}^n{(x_i - \mu)^2}\bigg)
\end{aligned}
$$

which we recognise as a gamma distribution with a shape of $n/2$ and a scale of $\frac{1}{2}\sum_{i=1}^n{(x_i - \mu)^2}$

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

Calculate the length, the sum and the sum of squares traversing the
list only once using the
[foldl](http://hackage.haskell.org/package/foldl) package. Much better
than creating your own strict record and using *foldl'*.

> x2Sum, xSum, n :: Double
> (x2Sum, xSum, n) = L.fold stats xs
>   where
>     stats = (,,) <$>
>             (L.premap (\x -> x * x) L.sum) <*>
>             L.sum <*>
>             L.genericLength
>
> xBar :: Double
> xBar = xSum / n

> gibbsSampler :: MonadRandom m => Double -> m (Maybe ((Double, Double), Double))
> gibbsSampler oldTau = do
>   newMu <- sample (Normal xBar (recip (n * oldTau)))
>   let shape = 0.5 * n
>       scale = 0.5 * (x2Sum + n * newMu^2 - 2 * n * newMu * xBar)
>   newTau <- sample (Gamma shape scale)
>   return $ Just ((newMu, newTau), newTau)


> gibbsSamples :: [(Double, Double)]
> gibbsSamples = evalState (ML.unfoldrM gibbsSampler 1.0) (pureMT 1)

> g4Sum, g3Sum, g2Sum, gSum, m :: Double
> (g4Sum, g3Sum, g2Sum, gSum, m) = L.fold stats (take 10000 $ map fst gibbsSamples)
>   where
>     stats = (,,,,) <$>
>             (L.premap (\x -> x * x * x * x) L.sum) <*>
>             (L.premap (\x -> x * x * x) L.sum) <*>
>             (L.premap (\x -> x * x) L.sum) <*>
>             L.sum <*>
>             L.genericLength

> moments xs = foldl' f (0.0, 0.0, 0.0, 0.0, 0.0) xs
>   where
>     f (m4, m3, m2, m1, n) x = (m4', m3', m2', m1', n')
>       where
>         n'  = n + 1
>         delta = x - m1
>         delta_n = delta / n'
>         delta_n2 = delta_n * delta_n
>         term1 = delta * delta_n * n
>         m1' = m1 + delta_n
>         m4' = m4 + term1 * delta_n2 * (n'*n' - 3*n' + 3) + 6 * delta_n2 * m2 - 4 * delta_n * m3
>         m3' = m3 + term1 * delta_n * (n' - 2) - 3 * delta_n * m2
>         m2' = m2 + term1

> norms :: [Double]
> norms = evalState (replicateM 10000 (sample (Normal 0.0 1.0))) (pureMT 1)


> nSumSumSqr :: Fractional a => L.Fold a (a, a, a)
> nSumSumSqr = (,,) <$>
>              L.genericLength <*>
>              L.sum <*>
>              (L.premap (\x -> x * x) L.sum)
>
> 

> numBins :: Int
> numBins = 400

> hb :: HBuilder Double (Data.Histogram.Generic.Histogram V.Vector BinD Double)
> hb = forceDouble -<< mkSimple (binD lower numBins upper)
>   where
>     lower = xBar - 1.0*0.001
>     upper = xBar + 1.0*0.001
>
> hist :: Histogram V.Vector BinD Double
> hist = fillBuilder hb (take 10000000 $ map fst gibbsSamples)

```{.dia width='800'}
dia = image "diagrams/DataScienceHaskPost.png" 1.0 1.0
````

~~~~ {.r include="example1.bug"}
~~~~

> displayHeader :: FilePath -> Diagram B R2 -> IO ()
> displayHeader fn =
>   mainRender ( DiagramOpts (Just 900) (Just 700) fn
>              , DiagramLoopOpts False Nothing 0
>              )

> main :: IO ()
> main =
>   displayHeader "diagrams/DataScienceHaskPost.png"
>     (barDiag
>      (zip (map fst $ asList hist) (map snd $ asList hist)))
