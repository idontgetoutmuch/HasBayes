{-# LANGUAGE NoMonomorphismRestriction #-}

    module Main (
  main
  ) where

import qualified Data.Vector.Unboxed as V
import Data.Random.Source.PureMT
import Data.Random
import Control.Monad.State

import System.Environment(getArgs)
import Data.Colour.Names
import Data.Colour
import Control.Lens
import Data.Default.Class
import Data.Time.LocalTime
import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Cairo

nSamples :: Int
nSamples = 10000

beta :: Double
beta = 2.0

testXs :: Int -> V.Vector Double
testXs m =
  V.fromList $
  evalState (replicateM m (sample StdUniform))
  (pureMT 2)

testEpsilons :: Int -> V.Vector Double
testEpsilons m =
  V.fromList $
  evalState (replicateM m (sample StdNormal))
  (pureMT 2)

testYs m = V.zipWith (\x e -> beta * x + e)
           (testXs m) (testEpsilons m)

prices' :: [(Double,Double)]
prices' = V.toList $ V.zip (testXs nSamples) (testYs nSamples)

chart = toRenderable layout
  where

    price1 = plot_points_style . point_color .~ opaque red
           $ plot_points_values .~ prices'
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

main = renderableToFile def chart "example2_big.png"
