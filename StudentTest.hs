{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

{-# LANGUAGE NoMonomorphismRestriction     #-}

module StudentTest where

import qualified Data.Vector.Unboxed as V
import Data.Random.Source.PureMT
import Data.Random
-- import Data.Random.Distribution.ChiSquare
import Data.Random.Distribution.T
import Control.Monad.State
import Data.Histogram ( asList )
import qualified Data.Histogram as H
import Data.Histogram.Fill
import Data.Histogram.Generic ( Histogram )
-- import Data.List
-- import Control.Parallel.Strategies

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

nSamples :: Int
nSamples = 100000

arbSeed :: Int
arbSeed = 2

sampleT :: MonadRandom m => m Double
sampleT = do
  x <- sample (T 2)
  return x

sampleNorm :: MonadRandom m => m Double
sampleNorm = do
  x <- sample StdNormal
  return x

simpleXs :: [Double]
simpleXs =
  evalState (replicateM nSamples sampleT)
            (pureMT $ fromIntegral arbSeed)

simpleYs :: [Double]
simpleYs =
  evalState (replicateM nSamples sampleNorm)
            (pureMT $ fromIntegral arbSeed)

numBins :: Int
numBins = 400

startMu :: Double
startMu = 0.0

hb :: HBuilder Double (Data.Histogram.Generic.Histogram V.Vector BinD Double)
hb = forceDouble -<< mkSimple (binD (-20.0) numBins 20.0)

hist :: Histogram V.Vector BinD Double
hist = fillBuilder hb simpleXs

histNorm :: Histogram V.Vector BinD Double
histNorm = fillBuilder hb simpleYs

displayHeader :: FilePath -> Diagram B R2 -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

main :: IO ()
main = do
  displayHeader "diagrams/StudentHist.png"
    (barDiag MCMCAnal
     (zip (map fst $ asList hist)     (map snd $ asList hist))
     (zip (map fst $ asList histNorm) (map snd $ asList histNorm)))
  displayHeader "diagrams/wikiBarChart.png" (wikiDiag True)

