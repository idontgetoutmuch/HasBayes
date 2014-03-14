{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

module LinRegAux (
    denv
  , diagGamma
  , diagNormal
  , diag
  ) where

import System.IO.Unsafe

import Graphics.Rendering.Chart.Backend.Diagrams
import Graphics.Rendering.Chart

import Numeric.SpecFunctions

import Text.Printf

import Data.Colour
import Control.Lens hiding ( (#) )
import Diagrams.Backend.Cairo.CmdLine

import Diagrams.Prelude hiding ( sample, render )

import Data.Default.Class


denv :: DEnv
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 500 500

logGammaPdf :: Double -> Double -> Double -> Double
logGammaPdf alpha beta x = unNorm - logGamma alpha
  where
    unNorm = alpha * (log beta) + (alpha - 1) * log x - beta * x

gammaPlot :: Double -> Double -> Double -> Double -> Graphics.Rendering.Chart.Renderable ()
gammaPlot a b a' b' = toRenderable layout
  where
    am :: Double -> Double
    am x = exp (logGammaPdf a b x)

    am' :: Double -> Double
    am' x = exp (logGammaPdf a' b' x)

    gammaPlot1 = plot_lines_values .~ [[ (x,(am x)) | x <- [0,(0.05)..20]]]
                 $ plot_lines_style  . line_color .~ opaque blue
                 $ plot_lines_title .~ "prior shape = " ++ printf "%3.3f" (a :: Double) ++
                                             " rate = " ++ printf "%3.3f" (b :: Double)
                 $ def

    gammaPlot2 = plot_lines_values .~ [[ (x,(am' x)) | x <- [0,(0.05)..20]]]
                 $ plot_lines_style   . line_color .~ opaque red
                 $ plot_lines_title   .~ "posterior shape = " ++ printf "%3.3f" a' ++
                                                   " rate = " ++ printf "%3.3f" b'
                 $ def

    layout = layout_title .~ "Gamma Prior and Posterior (10 Observations)"
             $ layout_plots .~ [toPlot gammaPlot1,
                                toPlot gammaPlot2]
             $ def

diagGamma :: Double -> Double -> Double -> Double -> QDiagram Cairo R2 Any
diagGamma a b a' b' = fst $ runBackend denv (render (gammaPlot a b a' b') (500, 500))

normalPdf :: Double -> Double -> Double -> Double
normalPdf mu sigma x = recip (sigma * sqrt (2.0 * pi)) * exp(-(x - mu)**2 / (2 * sigma**2))

normalPlot :: Double -> Double ->
              Double -> Double ->
              Double -> Double ->
              Graphics.Rendering.Chart.Renderable ()
normalPlot a b a' b' a'' b'' = toRenderable layout
  where
    am :: Double -> Double
    am x =normalPdf a b x

    am' :: Double -> Double
    am' x = normalPdf a' b' x

    am'' :: Double -> Double
    am'' x = normalPdf a'' b'' x

    normalPlot1 = plot_lines_values .~ [[ (x,(am x)) | x <- [-2,(-1.95)..5]]]
                  $ plot_lines_style  . line_color .~ opaque blue
                  $ plot_lines_title .~ "prior mean = " ++ printf "%3.3f" (a :: Double) ++
                                              " var = " ++ printf "%3.3f" (b :: Double)
                 $ def

    normalPlot2 = plot_lines_values .~ [[ (x,(am' x)) | x <- [-2,(-1.95)..5]]]
                 $ plot_lines_style   . line_color .~ opaque red
                 $ plot_lines_title   .~ "post mean = " ++ printf "%3.3f" a' ++
                                              " var = " ++ printf "%3.3f" b'
                 $ def

    normalPlot3 = plot_lines_values .~ [[ (x,(am'' x)) | x <- [-2,(-1.95)..5]]]
                 $ plot_lines_style   . line_color .~ opaque green
                 $ plot_lines_title   .~ "post mean = " ++ printf "%3.3f" a'' ++
                                              " var = " ++ printf "%3.3f" b''
                 $ def

    layout = layout_title .~ "Normal Prior and Posterior (10 & 100 Observations)"
             $ layout_plots .~ [toPlot normalPlot1,
                                toPlot normalPlot2,
                                toPlot normalPlot3]
             $ def

diagNormal :: Double -> Double ->
              Double -> Double ->
              Double -> Double ->
              QDiagram Cairo R2 Any
diagNormal a b a' b' a'' b'' = fst $ runBackend denv
                               (render (normalPlot a b a' b' a'' b'') (500, 500))


chart :: Colour Double ->[(Double, Double)] -> Graphics.Rendering.Chart.Renderable ()
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

diag :: Colour Double -> [(Double, Double)] -> QDiagram Cairo R2 Any
diag c prices = fst $ runBackend denv (render (chart c prices) (500, 500))

