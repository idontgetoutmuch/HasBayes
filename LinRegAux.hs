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
  , diagNormals
  , diag
  , barDiag
  , MCMCAnal (..)
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

normalPlots :: [(Double, Double, Colour Double, String)] ->
               Graphics.Rendering.Chart.Renderable ()
normalPlots abs = toRenderable layout
  where

    lower a b = a - 5*b
    upper a b = a + 5*b
    gap     b = 10*b / 1000.0

    am :: Double -> [Double]
    am x = map (\(a, b, _, _) -> normalPdf a b x) abs

    normalPlots = zipWith (
      \(a, b, c, l) n ->
      plot_lines_values .~ [[ (x, (am x)!!n) |
                              x <- [lower a b,lower a b + gap b .. upper a b]]]
      $ plot_lines_style  . line_color .~ opaque c
      $ plot_lines_title .~ l ++ " mean = " ++ printf "%3.3f" a ++
                                  " var = " ++ printf "%3.3f" b
                   $ def
      ) abs [0..]

    layout =   layout_title .~ "Normal Prior and Posterior"
             $ layout_plots .~ (map toPlot normalPlots)
             $ def

diagNormals :: [(Double, Double, Colour Double, String)] ->
               QDiagram Cairo R2 Any
diagNormals abs = fst $ runBackend denv
                  (render (normalPlots abs) (500, 500))

data MCMCAnal = MCMC | Anal | MCMCAnal

barChart :: MCMCAnal ->
            [(Double, Double)] ->
            [(Double, Double)] ->
            Graphics.Rendering.Chart.Renderable ()
barChart pt bvs bvs' = toRenderable layout
  where
    layout =
      layout_title .~ title pt
      $ layout_y_axis . laxis_title .~ "Frequency"
      $ layout_plots .~ (map plotBars $ plots pt)
      $ def

    title MCMC     = "Posterior via MCMC"
    title Anal     = "Analytic Posterior"
    title MCMCAnal = "MCMC and Analytic Posteriors Overlaid"

    plots MCMC     = [ bars1 ]
    plots Anal     = [ bars2 ]
    plots MCMCAnal = [ bars1, bars2 ]

    bars1 =
      plot_bars_titles .~ ["MCMC"]
      $ plot_bars_values .~ addIndexes (map return $ map snd bvs)
      $ plot_bars_style .~ BarsClustered
      $ plot_bars_item_styles .~ [(solidFillStyle (blue `withOpacity` 0.25), Nothing)]
      $ def

    bars2 =
      plot_bars_titles .~ ["Analytic"]
      $ plot_bars_values .~ addIndexes (map return $ map snd bvs')
      $ plot_bars_style .~ BarsClustered
      $ plot_bars_item_styles .~ [(solidFillStyle (red `withOpacity` 0.25), Nothing)]
      $ def

barDiag :: MCMCAnal ->
           [(Double, Double)] ->
           [(Double, Double)] ->
           QDiagram Cairo R2 Any
barDiag pt bvs bvs' = fst $ runBackend denv (render (barChart pt bvs bvs') (500, 500))

