{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

{-# LANGUAGE NoMonomorphismRestriction #-}

module Main where

import qualified Data.Vector.Unboxed as V
import Data.Random.Source.PureMT
import Data.Random
import Control.Monad.State

import Data.Colour
import Control.Lens hiding ( (#) )
import Graphics.Rendering.Chart hiding ( translate )
import Diagrams.Backend.Cairo.CmdLine

import Diagrams.Backend.CmdLine
import Diagrams.Prelude hiding ( sample, render )

import Data.Default.Class

import Graphics.Rendering.Chart.Backend.Diagrams

import System.IO.Unsafe

nSamples :: Int
nSamples = 100

betaHyperPrior :: Double
betaHyperPrior = 2.0

nuHyperPrior :: Double
nuHyperPrior = 1.0

vHyperPrior :: Double
vHyperPrior = 1.0

vHyperPriorInv :: Double
vHyperPriorInv = recip vHyperPrior

s2HyperPrior :: Double
s2HyperPrior = 1.0

testXs :: Int -> V.Vector Double
testXs m =
  V.fromList $
  evalState (replicateM m (sample StdUniform))
  (pureMT 2)

testXs' :: Int -> V.Vector Double
testXs' m =
  V.fromList $
  evalState (replicateM m (sample StdUniform))
  (pureMT 3)

testEpsilons :: Int -> V.Vector Double
testEpsilons m =
  V.fromList $
  evalState (replicateM m (sample StdNormal))
  (pureMT 2)

testEpsilons' :: Int -> V.Vector Double
testEpsilons' m =
  V.fromList $
  evalState (replicateM m (sample StdNormal))
  (pureMT 3)

testYs :: Int -> V.Vector Double
testYs m = V.zipWith (\x e -> betaHyperPrior * x + e)
           (testXs m) (testEpsilons m)

testYs' :: Int -> V.Vector Double
testYs' m = V.zipWith (\x e -> betaHyperPrior * x + e)
           (testXs' m) (testEpsilons' m)

posterior :: Int -> (Double, Double)
posterior m = (nuPosterior, vPosteriorInv)
  where
    xs       = testXs m
    xsSquare = V.sum $ V.zipWith (*) xs xs

    nuPosterior = nuHyperPrior + fromIntegral m

    vPosteriorInv = vHyperPriorInv + xsSquare

prices :: [(Double,Double)]
prices = V.toList $ V.zip (testXs nSamples) (testYs nSamples)

prices' :: [(Double,Double)]
prices' = V.toList $ V.zip (testXs' nSamples) (testYs' nSamples)

chart :: Colour Double -> [(Double, Double)] -> Graphics.Rendering.Chart.Renderable ()
chart c prices = toRenderable layout
  where

    price1 = plot_points_style . point_color .~ opaque c
           $ plot_points_values .~ prices
           $ plot_points_title .~ "price 1"
           $ def

    layout = layoutlr_title .~"Price History"
           $ layoutlr_left_axis . laxis_override .~ axisGridHide
           $ layoutlr_right_axis . laxis_override .~ axisGridHide
           $ layoutlr_x_axis . laxis_override .~ axisGridHide
           $ layoutlr_plots .~ [Left (toPlot price1),
                                Right (toPlot price1)]
           $ layoutlr_grid_last .~ False
           $ def

displayHeader :: FilePath -> Diagram B R2 -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 600) fn
             , DiagramLoopOpts False Nothing 0
             )

denv :: DEnv
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 500 500

diag :: Colour Double -> [(Double, Double)] -> QDiagram Cairo R2 Any
diag c prices = fst $ runBackend denv (render (chart c prices) (500, 500))

main :: IO ()
main = do
  displayHeader "TestInteractive.png"
    ((diag red prices # scaleX 0.4 # scaleY 0.4 # translate (r2 (0.1, 0.0)))
     <> (diag blue prices' # scaleX 0.6 # scaleY 0.6 # translate (r2 (0.2, 0.0))))
  putStrLn "Hello"
