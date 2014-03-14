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

Test Data
---------

Let's create some noisy test data.

> testXs :: Int -> Int -> V.Vector Double
> testXs seed nSamples =
>   V.fromList $
>   evalState (replicateM nSamples (sample StdUniform))
>   (pureMT $ fromIntegral seed)

> testEpsilons :: Int -> Int ->V.Vector Double
> testEpsilons seed nSamples =
>   V.fromList $
>   evalState (replicateM nSamples (sample StdNormal))
>   (pureMT $ fromIntegral seed)
>
> testYs :: V.Vector Double -> V.Vector Double -> V.Vector Double
> testYs xs es =
>   V.zipWith (\x e -> d * x + e) xs es
>
> testVs :: Int -> Int -> V.Vector (Double, Double)
> testVs seed nSamples = V.zip xs ys
>   where
>     xs = testXs seed nSamples
>     es = testEpsilons (seed + 1) nSamples
>     ys = testYs xs es
>
> testData :: Int -> Int -> [(Double,Double)]
> testData seed nSamples = V.toList $ testVs seed nSamples

We can look at this

    [ghci]
    import LinReg
    testXs 2 5
    testEpsilons 2 5

but a picture paints a thousand words, and given we have taken the
variance of the added noise to be pretty big, it is really only after
about 1,000 samples that it starts to look like there is a linear
relationship.

```{.dia height='600'}
import LinReg
import LinRegAux

dia = ((diag red (testData 2 10) # scaleX 0.3 # scaleY 0.3)
       |||
       (diag blue (testData 3 100) # scaleX 0.3 # scaleY 0.3))
      ===
      ((diag blue (testData 3 1000) # scaleX 0.3 # scaleY 0.3)
       |||
       (diag green (testData 6 10000) # scaleX 0.3 # scaleY 0.3))
````

Conjugate Prior
===============

Ultimately we are going to use a Monte Carlo Markov Chain (MCMC)
method to calculate the posterior distribution. Before we do and so
that we have something to test against let us first consider a
conjugate prior. A conjugate prior is a member of a family of
distributions in which the posterior is a member of the same family.

Let us take a prior from the normal-Gamma distribution $NG(\mu,
\lambda | \mu_0, \kappa_0, \alpha_0, \beta_0)$

$$
\pi(\theta) \propto \lambda^{1/2} \exp\Bigg[{-\frac{\lambda}{2V}(\beta -\mu_0)^2}\Bigg]\lambda^{a - 1}\exp\Bigg[{-b\lambda}\Bigg]
$$


The likelihood

$$
f(x | \theta) \propto  \lambda^{1/2} \exp \Bigg[-\frac{\lambda}{2}(\beta - \hat{\beta})^2\sum_{i = 1}^N x_i^2\Bigg] \lambda^{\nu / 2}\exp{\Bigg[-\frac{\lambda\nu}{2s^{-2}}\Bigg]}
$$

The posterior

$$
f(\theta | x) \propto \lambda^{a - 1 + \nu / 2 + 1 / 2 + 1 / 2} \exp\Bigg[{-b\lambda}{-\frac{\lambda}{2V}(\beta -\mu_0)^2}{-\frac{\lambda\nu}{2s^{-2}}}{-\frac{\lambda}{2}(\beta - \hat{\beta})^2\sum_{i = 1}^N x_i^2}\Bigg]
$$

Equating powers of $\lambda$ we get

$$
a' = a + \frac{\nu}{2} + \frac{1}{2} = a + \frac{N}{2}
$$

Now let us examine the factors inside the exponential

$$
-b'\lambda - \frac{\lambda}{2V'}(\beta - \mu_0')^2 =
-b\lambda -  \frac{\lambda}{2V}(\beta - \mu_0)^2
-\frac{\lambda\nu}{2s^{-2}} - \frac{\lambda}{2}(\beta - \hat{\beta})^2\sum_{i = 1}^N x_i^2
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
\frac{\mu_0'}{V'} = \frac{\mu_0}{V} + \hat{\beta}\sum_{i = 1}^N x_i^2
$$

and so

$$
\mu_0' = V'\bigg(\mu_0V^{-1} + \hat{\beta}\sum_{i = 1}^N x_i^2\bigg)
$$

Equating constant terms

$$
b' + \frac{\mu_0'^2}{2V'} =
b + \frac{\mu_0^2}{2V} + \frac{\nu}{2s^{-2}} +
\frac{\hat{\beta}^2\sum_{i = 1}^N x_i^2}{2}
$$

and so

$$
b' =
b + \frac{1}{2}
\bigg(\frac{\mu_0^2}{V} - \frac{\mu_0'^2}{V'} +
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
\bigg(\frac{\mu_0^2}{V} - \frac{\mu_0'^2}{V'} +
\sum y_i^2\bigg)
$$

>
> nSamples, seed :: Int
> nSamples = 10
> seed = 2

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

> zs :: V.Vector (Double, Double)
> zs = testVs seed nSamples
>
> xs :: V.Vector Double
> xs = V.map fst zs

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
> ys = V.map snd zs
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


> displayHeader :: FilePath -> Diagram B R2 -> IO ()
> displayHeader fn =
>   mainRender ( DiagramOpts (Just 900) (Just 600) fn
>              , DiagramLoopOpts False Nothing 0
>              )

A Gibbs Sampler
===============

> main :: IO ()
> main = do
>   displayHeader "Gamma.png" (diagGamma 4.0 1.0 a' b')
>   displayHeader "TestInteractive.png"
>     (((diag red (testData 2 10) # scaleX 0.4 # scaleY 0.4 # translate (r2 (0.1, 0.0)))
>       === (diag blue (testData 4 100) # scaleX 0.4 # scaleY 0.4 # translate (r2 (0.2, 0.0))))
>     |||
>     ((diag red (testData 2 10) # scaleX 0.4 # scaleY 0.4 # translate (r2 (0.1, 0.0)))
>       === (diag blue (testData 4 100) # scaleX 0.4 # scaleY 0.4 # translate (r2 (0.2, 0.0)))))
>   putStrLn "Hello"
