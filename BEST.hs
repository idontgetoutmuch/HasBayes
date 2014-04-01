{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

{-# LANGUAGE NoMonomorphismRestriction     #-}

module ConjMCMCSimple where

import qualified Data.Vector.Unboxed as V
import Data.Random.Source.PureMT
import Data.Random
import Control.Monad.State
import Data.Histogram ( asList )
import qualified Data.Histogram as H
import Data.Histogram.Fill
import Data.Histogram.Generic ( Histogram )
import Data.List
import Control.Parallel.Strategies

import Diagrams.Backend.Cairo.CmdLine

import Diagrams.Backend.CmdLine
import Diagrams.Prelude hiding ( sample, render )

import LinRegAux

mu0 :: Double
mu0 = 11.0

sigma_0 :: Double
sigma_0 = 2.0

sigma :: Double
sigma = 1.0

hierarchicalSample :: MonadRandom m => m Double
hierarchicalSample = do
  mu <- sample (Normal mu0 sigma_0)
  x  <- sample (Normal mu sigma)
  return x
