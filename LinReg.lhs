% Bayesian Regession and Gibbs Sampler in Haskell
% Dominic Steinitz
% 9th March 2014

---
bibliography: Bayes.bib
---

Introduction
============

Suppose you have some data to which you wish to fit a straight
line. For simplicity, let us consider a line which always goes through
origin. Clearly, this is not realistic but will simplify the
calculations and prepare us for doing multivariate regression.

Preamble
--------

> {-# OPTIONS_GHC -Wall                      #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults    #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods  #-}
> {-# OPTIONS_GHC -fno-warn-orphans          #-}

> {-# LANGUAGE NoMonomorphismRestriction #-}

> module LinReg where
>
> import qualified Data.Vector.Unboxed as V
> import Data.Random.Source.PureMT
> import Data.Random
> import Control.Monad.State
>
> import Diagrams.Backend.Cairo.CmdLine
>
> import Diagrams.Backend.CmdLine
> import Diagrams.Prelude hiding ( sample, render )
>
> import LinRegAux

Conjugate Prior
===============

Ultimately we are going to use a Monte Carlo Markov Chain (MCMC)
method to calculate the posterior distribution. Before we do and so
that we have something to test against let us first consider a
conjugate prior. A conjugate prior is a member of a family of
distributions in which the posterior is a member of the same family.

> testXs :: Int -> V.Vector Double
> testXs m =
>   V.fromList $
>   evalState (replicateM m (sample StdUniform))
>   (pureMT 2)

    [ghci]
    import LinReg
    testXs 5

The likelihood

$$
f(x | \theta) \propto  h^{1/2} \exp \Bigg[-\frac{h}{2}(\beta - \hat{\beta})^2\sum_{i = 1}^N x_i^2\Bigg] h^{\nu / 2}\exp{\Bigg[-\frac{h\nu}{2s^{-2}}\Bigg]}
$$

The prior

$$
\pi(\theta) \propto h^{1/2} \exp\Bigg[{-\frac{h}{2V}(\beta -d)^2}\Bigg]h^{a - 1}\exp\Bigg[{-bh}\Bigg]
$$

The posterior

$$
f(\theta | x) \propto h^{a - 1 + \nu / 2 + 1 / 2 + 1 / 2} \exp\Bigg[{-bh}{-\frac{h}{2V}(\beta -d)^2}{-\frac{h\nu}{2s^{-2}}}{-\frac{h}{2}(\beta - \hat{\beta})^2\sum_{i = 1}^N x_i^2}\Bigg]
$$

Equating powers of $h$ we get

$$
a' = a + \frac{\nu}{2} + \frac{1}{2} = a + \frac{N}{2}
$$

Now let us examine the factors inside the exponential

$$
-b'h - \frac{h}{2V'}(\beta - d')^2 =
-bh -  \frac{h}{2V}(\beta - d)^2
-\frac{h\nu}{2s^{-2}} - \frac{h}{2}(\beta - \hat{\beta})^2\sum_{i = 1}^N x_i^2
$$

Equating terms in $\beta^2$

$$
\frac{1}{V'} = \frac{1}{V} + \sum_{i = 1}^N x_i^2
$$

and so

$$
V' = \bigg(V^{-1} + \sum_{i = 1}^N x_i^2\bigg)^{-1}
$$

Equating terms in $\beta$

$$
\frac{d'}{V'} = \frac{d}{V} + \hat{\beta}\sum_{i = 1}^N x_i^2
$$

and so

$$
d' = V'\bigg(dV^{-1} + \hat{\beta}\sum_{i = 1}^N x_i^2\bigg)
$$

Equating constant terms

$$
b' + \frac{d'^2}{2V'} =
b + \frac{d^2}{2V} + \frac{\nu}{2s^{-2}} +
\frac{\hat{\beta}^2\sum_{i = 1}^N x_i^2}{2}
$$

and so

$$
b' =
b + \frac{1}{2}
\bigg(\frac{d^2}{V} - \frac{d'^2}{V'} +
\frac{\nu}{s^{-2}} +
\hat{\beta}^2\sum_{i = 1}^N x_i^2\bigg)
$$

But we know that

$$
\hat{\beta} = \frac{\sum x_i y_i}{\sum x^2_i}\, \text{and}\, s^2 = \frac{\sum (y_i - \hat{\beta}x_i)^2}{\nu}
$$

Thus we can rewrite

$$
\begin{aligned}
\frac{\nu}{s^{-2}} +
\hat{\beta}^2\sum_{i = 1}^N x_i^2 &=
\sum (y_i - \hat{\beta}x_i)^2 + \hat{\beta}^2\sum x_i^2 \\
&= \sum y_i^2 - 2\hat{\beta}\sum x_i y_i + \hat{\beta}^2\sum x_i^2 + \hat{\beta}^2\sum x_i^2 \\
&= \sum y_i^2
\end{aligned}
$$

giving

$$
b' =
b + \frac{1}{2}
\bigg(\frac{d^2}{V} - \frac{d'^2}{V'} +
\sum y_i^2\bigg)
$$

>
> nSamples :: Int
> nSamples = 10

> a :: Double
> a = 1.0
>
> a' :: Double
> a' = a + (fromIntegral nSamples) / 2.0

    [ghci]
    import LinReg
    a'

> v :: Double
> v = 1.0

> xs :: V.Vector Double
> xs = testXs nSamples

> xs2 :: Double
> xs2 = V.sum $ V.map (**2) xs

> v' :: Double
> v' = recip (recip v + xs2)

    [ghci]
    v'

> d :: Double
> d = 2.0
>
> betaHat :: Double
> betaHat = (V.sum $ V.zipWith (*) xs ys) / xs2
>
> d' :: Double
> d' = v' * (d * recip v + betaHat * xs2)

    [ghci]
    d'

> b :: Double
> b = v / (2.0 * recip v**2)
>
> ys :: V.Vector Double
> ys = testYs nSamples
>
> ys2 :: Double
> ys2 = V.sum $ V.map (**2) ys
>
> b' :: Double
> b' = b + 0.5 * (d**2 / v - d'**2 / v' + ys2)

    [ghci]
    b'

We can plot the prior and posterior gamma distributions. Even if we
only sample 10 observations we can see that the posterior gives quite
a tight estimate. In the example below we take the prior gamma to have
shape 4.0 and rate 1.0 but this does not seem to have much influence
of the posterior,

```{.dia height='600'}
import LinReg
import LinRegAux

dia = diagGamma 4.0 1.0 a' b'
````

We can also plot the prior and posterior normal distributions assuming
that $h = 1.0$. Note that after 10 observations we do not get a very
good estimate but we get a much better one after 100 observations.

```{.dia height='600'}
import LinReg
import LinRegAux

dia = diagNormal d (sqrt v) d' (sqrt v') 1.879657598238415 (sqrt 3.146906057947862e-2)
````

```{.dia height='600'}
import LinReg
import LinRegAux

dia = (diag red prices # scaleX 0.4 # scaleY 0.4 # translate (r2 (0.1, 0.0)))
      |||
      (diag blue prices' # scaleX 0.4 # scaleY 0.4 # translate (r2 (0.2, 0.0)))
````

> testXs' :: Int ->V.Vector Double
> testXs' m =
>   V.fromList $
>   evalState (replicateM m (sample StdUniform))
>   (pureMT 3)
>
> testEpsilons :: Int ->V.Vector Double
> testEpsilons m =
>   V.fromList $
>   evalState (replicateM m (sample StdNormal))
>   (pureMT 2)
>
> testEpsilons' :: Int ->V.Vector Double
> testEpsilons' m =
>   V.fromList $
>   evalState (replicateM m (sample StdNormal))
>   (pureMT 3)
>
> testYs :: Int ->V.Vector Double
> testYs m = V.zipWith (\x e -> d * x + e)
>            (testXs m) (testEpsilons m)
>
> testYs' :: Int ->V.Vector Double
> testYs' m = V.zipWith (\x e -> d * x + e)
>            (testXs' m) (testEpsilons' m)
>
> prices :: [(Double,Double)]
> prices = V.toList $ V.zip (testXs nSamples) (testYs nSamples)
>
> prices' :: [(Double,Double)]
> prices' = V.toList $ V.zip (testXs' nSamples) (testYs' nSamples)
>
> displayHeader :: FilePath -> Diagram B R2 -> IO ()
> displayHeader fn =
>   mainRender ( DiagramOpts (Just 900) (Just 600) fn
>              , DiagramLoopOpts False Nothing 0
>              )
>



A Gibbs Sampler
===============

> main :: IO ()
> main = do
>   displayHeader "Gamma.png" (diagGamma 4.0 1.0 a' b')
>   displayHeader "TestInteractive.png"
>     ((diag red prices # scaleX 0.4 # scaleY 0.4 # translate (r2 (0.1, 0.0)))
>      ||| (diag blue prices' # scaleX 0.4 # scaleY 0.4 # translate (r2 (0.2, 0.0))))
>   putStrLn "Hello"
