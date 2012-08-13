\documentclass[12pt]{article}
%include polycode.fmt
\usepackage[pdftex,pagebackref,letterpaper=true,colorlinks=true,pdfpagemode=none,urlcolor=blue,linkcolor=blue,citecolor=blue,pdfstartview=FitH]{hyperref}

\usepackage{amsmath,amsfonts}
\usepackage{graphicx}
\usepackage{color}

\setlength{\oddsidemargin}{0pt}
\setlength{\evensidemargin}{0pt}
\setlength{\textwidth}{6.0in}
\setlength{\topmargin}{0in}
\setlength{\textheight}{8.5in}

\setlength{\parindent}{0in}
\setlength{\parskip}{5px}


\begin{document}

\section{A Simple Example}
In Section 7.1 of \href{http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/}{\lq\lq Doing Bayesian Data Analysis\rq\rq}, John Kruschke gives a simple example of the Metropolis algorithm, in which we generate samples from a distribution without knowing the distribution itself. Of course the example is contrived as we really do know the distribtion. In the particular example, ${\cal P}(X = i) = i/k$ for $i = 1 \ldots n$ where $n$ is some fixed natural number and $k$ is a normalising constant.

Here's the algorithm in Haskell. We use the \href{http://hackage.haskell.org/package/random-fu-0.2.3.0}{random-fu} package and the \href{http://hackage.haskell.org/package/rvar-0.2.0.1}{rvar} package for random variables and the \href{http://hackage.haskell.org/package/random-1.0.1.1}{random} package to supply random numbers.

\begin{code}
{-# LANGUAGE ScopedTypeVariables, NoMonomorphismRestriction #-}

import Data.Random
import Data.RVar
import System.Random
import Control.Monad.State
import Data.List
\end{code}

We pick an explicit seed and set $n = 7$ (in Kruschke's example this is the number of islands that a politician visits).

\begin{code}
seed :: Int
seed = 2
numIslands :: Int
numIslands = 7
\end{code}

And we pick an arbitrary sample size.

\begin{code}
n = 11112
\end{code}

We generate proposed jumps which are either one step up or one step down.

\begin{code}
proposedJumps :: Int -> [Int]
proposedJumps seed =
  map f $ fst $
  runState (replicateM n $ sampleRVar $ uniform False True) (mkStdGen seed)
    where f False = negate 1
          f True  = 1
\end{code}

And we generate samples from the uniform distribution on $[0, 1]$ which will allow us to determine whether to accept or reject a proposed jump.

\begin{code}
acceptOrRejects :: Int -> [Double]
acceptOrRejects seed =
  fst $ runState (replicateM n $ sampleRVar $ uniform 0 1) (mkStdGen seed)
\end{code}

We pretend we only know a measure of how often we pick a given number but not the normalising constant to make this measure a probability measure.

\begin{code}
p n | n >= 1 && n <= numIslands = n
    | otherwise                 = 0
\end{code}

We define a function which defines one step of the Metropolis algorithm.

\begin{code}
f currentPosition (proposedJump, acceptOrReject) =
  if acceptOrReject < probAccept
    then
      currentPosition + proposedJump
    else
      currentPosition
  where
    probAccept = min 1 (pNew / pCur)
    pNew = fromIntegral $ currentPosition + proposedJump
    pCur = fromIntegral currentPosition
\end{code}

Let's try this out with a somewhat arbitrary burn in period starting at position 3.

\begin{code}
runMC seed = map (/ total) numVisits
  where
    total     = sum numVisits
    numVisits = map (\j -> fromIntegral $ length $ filter (==j) $ xs) [1 .. numIslands]
    xs        = drop (n `div` 10) $ scanl f 3 (zip (proposedJumps seed) (acceptOrRejects seed))
\end{code}

We can then compute the root mean square error for this particular sample size.

\begin{code}
actual = map (/s) xs
           where
             xs = map fromIntegral [1 .. numIslands]
             s  = sum xs

rmsError n = sqrt $ (/(fromIntegral numIslands)) $
             sum $ map (^2) $ zipWith (-) actual (runMC n)
\end{code}

\section{A Bernoulli Example}

We can also code the Metropolis algorithm for the example in which samples are drawn from a Bernoulli distribution.
First we can define the likelihood for drawing drawing a specific sequence of 0's and 1's from a binomial distribution with parrameter $\theta$.
We also define the example sample of data given in \lq\lq Doing Bayesian Data Analysis\rq\rq.

\begin{code}

likelihood :: Int -> Int -> Double -> Double
likelihood z n theta = theta^z * (1 - theta)^(n - z)

myData = [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
\end{code}

We define a function which defines one step of the Metropolis algorithm.

\begin{code}
oneStep currentPosition (proposedJump, acceptOrReject) =
  if acceptOrReject < probAccept
    then
      currentPosition + proposedJump
    else
      currentPosition
  where
    probAccept = min 1 (p (currentPosition + proposedJump) / p currentPosition)
    p x | x < 0 = 0
        | x > 1 = 0
        | otherwise = pAux myData x
    pAux :: [Double] -> Double -> Double
    pAux xs position = likelihood z n position
      where
        n = length xs
        z = length $ filter (== 1) xs
\end{code}

Finally we need some proposed jumps; following the example we generate these from a normal distribution ${\cal N}(0, 0.1)$.

\begin{code}
normals :: [Double]
normals =  fst $ runState (replicateM n (sampleRVar (normal 0 0.1))) (mkStdGen seed)
\end{code}

And now we can run the Metropolis algorithm and, for example, find the mean of the posterior.

\begin{code}
accepteds = drop (n `div` 10) $ scanl oneStep 0.5 (zip normals (acceptOrRejects seed))

mean = total / count
         where
           (total, count) = foldl' f (0.0, 0) accepteds
           f (s, c) x = (s+x, c+1)
\end{code}

\end{document}
