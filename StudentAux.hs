{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

{-# LANGUAGE NoMonomorphismRestriction     #-}

module StudentTest (
    main
  , Moments (..)
  ) where

import qualified Data.Vector.Unboxed as V
import Data.Random.Source.PureMT
import Data.Random
import Data.Random.Distribution.T
import Control.Monad.State
import Data.Histogram ( asList )
import Data.Histogram.Fill
import Data.Histogram.Generic ( Histogram )

-- import Diagrams.Backend.Cairo.CmdLine

-- import Diagrams.Backend.CmdLine
-- import Diagrams.Prelude hiding ( sample, render )

import Data.List

-- import LinRegAux

nSamples :: Int
nSamples = 1000000

arbSeed :: Int
arbSeed = 8

nu :: Integer
nu = 6

sampleT :: MonadRandom m => m Double
sampleT = do
  x <- sample (T nu)
  return x

simpleXs :: [Double]
simpleXs =
  evalState (replicateM nSamples sampleT)
            (pureMT $ fromIntegral arbSeed)

mean = (sum simpleXs) / fromIntegral nSamples
variance = (sum (map (**2) simpleXs)) / fromIntegral nSamples
skewness = (sum (map (**3) simpleXs) / fromIntegral nSamples) / variance**1.5
exKurtosis = (sum (map (**4) simpleXs) / fromIntegral nSamples) / variance**2


data Moments = Moments { m1 :: !Double
                       , m2 :: !Double
                       , m3 :: !Double
                       , m4 :: !Double
                       }
  deriving Show

-- m  = foldl' (\m x -> Moments { m1 = m1 m + x
--                               , m2 = m2 m + x**2
--                               , m3 = m3 m + x**3
--                               , m4 = m4 m + x**4
--                               }) (Moments 0.0 0.0 0.0 0.0) simpleXs

-- mean       = m1 m / fromIntegral nSamples
-- variance   = m2 m / fromIntegral nSamples
-- skewness   = (m3 m / fromIntegral nSamples) / variance**1.5
-- exKurtosis = (m4 m / fromIntegral nSamples) / variance**2

-- mean = (foldl' (+) 0 simpleXs) / fromIntegral nSamples
-- variance = (foldl' (+) 0 (map (**2) simpleXs)) / fromIntegral nSamples
-- skewness = (foldl' (+) 0 (map (**3) simpleXs) / fromIntegral nSamples) / variance**1.5
-- exKurtosis = (foldl' (+) 0 (map (**4) simpleXs) / fromIntegral nSamples) / variance**2


sampleNorm :: MonadRandom m => m Double
sampleNorm = do
  x <- sample StdNormal
  return x

simpleYs :: [Double]
simpleYs =
  evalState (replicateM nSamples sampleNorm)
            (pureMT $ fromIntegral arbSeed)

numBins :: Int
numBins = 400

hb :: HBuilder Double (Data.Histogram.Generic.Histogram V.Vector BinD Double)
hb = forceDouble -<< mkSimple (binD (-20.0) numBins 20.0)

hist :: Histogram V.Vector BinD Double
hist = fillBuilder hb simpleXs

histNorm :: Histogram V.Vector BinD Double
histNorm = fillBuilder hb simpleYs

-- displayHeader :: FilePath -> Diagram B R2 -> IO ()
-- displayHeader fn =
--   mainRender ( DiagramOpts (Just 900) (Just 700) fn
--              , DiagramLoopOpts False Nothing 0
--              )

main :: IO ()
main = do
  -- putStrLn $ show foo
  putStrLn $ show mean
  putStrLn $ show variance
  putStrLn $ show skewness
  putStrLn $ show exKurtosis
  -- displayHeader "diagrams/StudentHist.png"
  --   (barDiag MCMCAnal
  --    (zip (map fst $ asList hist)     (map snd $ asList hist))
  --    (zip (map fst $ asList histNorm) (map snd $ asList histNorm)))
  -- displayHeader "diagrams/wikiBarChart.png" (wikiDiag True)

