{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

module PopGrowthMCMCChart (
    diagEstsL
  , diagEsts3
  ) where

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Diagrams
import Diagrams.Backend.Cairo.CmdLine
import Diagrams.Prelude hiding ( render, Renderable )
import Data.Default.Class

import System.IO.Unsafe

denv :: DEnv
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 500 500

diagEstsL :: String -> String -> String ->
            [(Double, Double)] -> [(Double, Double)] -> QDiagram Cairo R2 Any
diagEstsL t l1 l2 ls ps =
  fst $ runBackend denv (render (chartEstsL t l1 l2 ls ps) (500, 500))

chartEstsL :: String -> String -> String ->
             [(Double, Double)] -> [(Double, Double)] -> Renderable ()
chartEstsL t l1 l2 lineVals pointVals = toRenderable layout
  where

    actuals = plot_lines_values .~ [lineVals]
              $ plot_lines_style  . line_color .~ opaque blue
              $ plot_lines_title .~ l1
              $ def

    estimas = plot_lines_values .~ [pointVals]
              $ plot_lines_style  . line_color .~ opaque green
              $ plot_lines_title .~ l2
              $ def

    layout = layout_title .~ t
           $ layout_plots .~ [toPlot actuals,
                              toPlot estimas]
           $ def

diagEsts3 :: String -> String -> String -> String ->
            [(Double, Double)] ->
            [(Double, Double)] ->
            [(Double, Double)] ->
            QDiagram Cairo R2 Any
diagEsts3 t l1 l2 l3 ls ps qs =
  fst $ runBackend denv (render (chartEsts3 t l1 l2 l3 ls ps qs) (500, 500))

chartEsts3 :: String -> String -> String -> String ->
             [(Double, Double)] ->
             [(Double, Double)] ->
             [(Double, Double)] ->
             Renderable ()
chartEsts3 t l1 l2 l3 lineVals pointVals qVals = toRenderable layout
  where

    actuals = plot_lines_values .~ [lineVals]
              $ plot_lines_style  . line_color .~ opaque blue
              $ plot_lines_title .~ l1
              $ def

    estimas = plot_lines_values .~ [pointVals]
              $ plot_lines_style  . line_color .~ opaque green
              $ plot_lines_title .~ l2
              $ def

    tbdimas = plot_lines_values .~ [qVals]
              $ plot_lines_style  . line_color .~ opaque red
              $ plot_lines_title .~ l3
              $ def

    layout = layout_title .~ t
           $ layout_plots .~ [toPlot actuals,
                              toPlot estimas,
                              toPlot tbdimas
                             ]
           $ def
