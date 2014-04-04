% Student's T and Space Leaks
% Dominic Steinitz
% 2nd April 2014

---
bibliography: Bayes.bib
---

Introduction
============

The other speaker at the Machine Learning Meetup at which I gave my
talk on automatic differentiation gave a very interesting talk on [A/B
testing](http://en.wikipedia.org/wiki/A/B_testing). Apparently this is
big business these days as attested by the fact I got 3 ads above the
wikipedia entry when I googled for it.

It seems that people tend to test with small sample sizes and to do so
very often, resulting in spurious results. Of course readers of
[XKCD](http://xkcd.com/) will be well aware of some of the
[pitfalls](http://xkcd.com/882/).

I thought a [Bayesian
approach](http://www.indiana.edu/~kruschke/BEST/BEST.pdf) might
circumvent some of the problems and set out to write a blog article
only to discover that there was no Haskell library for sampling from
[Student's t](http://en.wikipedia.org/wiki/Student%27s_t-distribution). Actually
there was one but is currently an unreleased part of
[random-fu](http://hackage.haskell.org/package/random-fu). So I set
about fixing this shortfall.

I thought I had better run a few tests so I calculated the sampled
mean, variance, [skewness](http://en.wikipedia.org/wiki/Skewness) and
[kurtosis](http://en.wikipedia.org/wiki/Kurtosis).

I wasn't really giving this my full attention and as a result ran into
a few problems with space. I thought these were worth sharing and that
is what this blog post is about. Hopefully, I will have time soon to
actually blog about the Bayesian equivalent of A/B testing.

Preamble
--------

> {-# OPTIONS_GHC -Wall                      #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults    #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods  #-}
> {-# OPTIONS_GHC -fno-warn-orphans          #-}
>
> {-# LANGUAGE NoMonomorphismRestriction     #-}
>
> module StudentTest (
>     main
>   ) where
>
> import qualified Data.Vector.Unboxed as V
> import Data.Random.Source.PureMT
> import Data.Random
> import Data.Random.Distribution.T
> import Control.Monad.State
> import Data.Histogram.Fill
> import Data.Histogram.Generic ( Histogram )
> import Data.List

Space Analysis
==============

Let's create a reasonable number of samples as the higher moments
converge quite slowly.

> nSamples :: Int
> nSamples = 1000000

An arbitrary seed for creating the samples.

> arbSeed :: Int
> arbSeed = 8

Student's t only has one parameter, the number of degrees of freedom.

> nu :: Integer
> nu = 6

Now we can do our tests by calculating the sampled values.

> ts :: [Double]
> ts =
>   evalState (replicateM nSamples (sample (T nu)))
>             (pureMT $ fromIntegral arbSeed)

> mean, variance, skewness, kurtosis :: Double
> mean = (sum ts) / fromIntegral nSamples
> variance = (sum (map (**2) ts)) / fromIntegral nSamples
> skewness = (sum (map (**3) ts) / fromIntegral nSamples) / variance**1.5
> kurtosis = (sum (map (**4) ts) / fromIntegral nSamples) / variance**2

This works fine for small sample sizes but not for the number we have chosen.

~~~~ { .shell }
./StudentTest +RTS -hc
Stack space overflow: current size 8388608 bytes.
Use `+RTS -Ksize -RTS' to increase it.
~~~~

It seems a shame that the function in the Prelude has this behaviour
but never mind let us ensure that we consume values strictly (they are
being produced lazily).

> mean' = (foldl' (+) 0 ts) / fromIntegral nSamples
> variance' = (foldl' (+) 0 (map (**2) ts)) / fromIntegral nSamples
> skewness' = (foldl' (+) 0 (map (**3) ts) / fromIntegral nSamples) / variance'**1.5
> kurtosis' = (foldl' (+) 0 (map (**4) ts) / fromIntegral nSamples) / variance'**2

We now have a space leak on the heap as using the ghc profiler below
shows. What went wrong?

```{.dia width='800'}
dia = image "diagrams/StudentTestFoldl.png" 1.0 1.0
````

If we only calculate the mean using foldl then all is well. Instead of
35M we only use 45K.

```{.dia width='800'}
dia = image "diagrams/StudentTestFoldlSingle.png" 1.0 1.0
````

Well that gives us a clue. The garbage collector cannot reclaim the
samples as they are needed for other calculations. What we need to do
is calculate the moments strictly altogether.

Let's create a strict record to do this.

> data Moments = Moments { m1 :: !Double
>                        , m2 :: !Double
>                        , m3 :: !Double
>                        , m4 :: !Double
>                        }
>   deriving Show

And calculate the results strictly.

>
> m  = foldl' (\m x -> Moments { m1 = m1 m + x
>                              , m2 = m2 m + x**2
>                              , m3 = m3 m + x**3
>                              , m4 = m4 m + x**4
>                              }) (Moments 0.0 0.0 0.0 0.0) ts
>
> mean''       = m1 m / fromIntegral nSamples
> variance''   = m2 m / fromIntegral nSamples
> skewness''   = (m3 m / fromIntegral nSamples) / variance''**1.5
> kurtosis'' = (m4 m / fromIntegral nSamples) / variance''**2

Now we have what we want; the program runs in small constant space.

```{.dia width='800'}
dia = image "diagrams/StudentTestFoldlFinal.png" 1.0 1.0
````

> main :: IO ()
> main = do
>   putStrLn $ show mean''
>   putStrLn $ show variance''
>   putStrLn $ show skewness''
>   putStrLn $ show kurtosis''

Oh and the moments give the expected answers.

    [ghci]
    mean''
    variance''
    skewness''
    kurtosis''

Running the Code
================

To run this you will need my
[version](https://github.com/idontgetoutmuch/random-fu) of
random-fu. The code for this article is
[here](https://github.com/idontgetoutmuch/HasBayes). You will need to
compile everything with profiling, something like

~~~~ { .shell }
ghc -O2 -main-is StudentTest StudentTest.lhs -prof
    -package-db=.cabal-sandbox/x86_64-osx-ghc-7.6.2-packages.conf.d
~~~~

Since you need all the packages to be built with profiling, you will
probably want to build using a sandbox as above. The only slightly
tricky aspect is building random-fu so it is in your sandbox.

~~~~ { .shell }
runghc Setup.lhs configure --enable-library-profiling
--package-db=<whatever>/HasBayes/.cabal-sandbox/x86_64-osx-ghc-7.6.2-packages.conf.d
--libdir=<whatever>/HasBayes/.cabal-sandbox/lib
~~~~