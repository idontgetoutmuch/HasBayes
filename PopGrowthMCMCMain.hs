{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

{-# LANGUAGE ViewPatterns                 #-}

module Main (
    main
    ) where

import Diagrams.Prelude
import Diagrams.Backend.CmdLine
import Diagrams.Backend.Cairo.CmdLine

import qualified Data.Vector as V
import qualified Data.ByteString.Lazy as BL
import Data.Csv
import Data.Time
import System.Locale
import Data.Char
import qualified Data.ByteString as B
import Control.Monad

import PopGrowthMCMC
import PopGrowthMCMCChart

displayHeader :: FilePath -> Diagram B R2 -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

main :: IO ()
main = do
  displayHeader "diagrams/AutoregressionVary1.png"
                (diag "Predicted Flow at Kingston Bridge (Varying Parameters)"
                      (zip [0..] (V.toList $ V.map exp $ V.take 800 flows))
                      (zip [0..] (V.toList $ V.map exp $ V.take 800 preds)))
