Conjugate Prior
===============

The likelihood

$$
f(x | \theta) \propto  h^{1/2} \exp \Bigg[-\frac{h}{2}(\beta - \hat{\beta})^2\sum_{i = 1}^N x_i^2\Bigg] h^{\nu / 2}\exp{\Bigg[-\frac{h\nu}{2s^{-2}}\Bigg]}
$$

The prior

$$
\pi(\theta) \propto h^{1/2} \exp\Bigg[{-\frac{h}{2V}(\beta -d)^2}\Bigg]h^{a - 1}\exp\Bigg[{-bh}\Bigg]
$$

The posterior

$$
f(\theta | x) \propto h^{a - 1 + \nu / 2 + 1 / 2 + 1 / 2} \exp\Bigg[{-bh}{-\frac{h}{2V}(\beta -d)^2}{-\frac{h\nu}{2s^{-2}}}{-\frac{h}{2}(\beta - \hat{\beta})^2\sum_{i = 1}^N x_i^2}\Bigg]
$$

Equating powers of $h$

$$
a' = a + \frac{\nu}{2} + \frac{1}{2} = a + \frac{N}{2}
$$

Now let us examine the factors inside the exponential

$$
-b'h - \frac{h}{2V'}(\beta - d')^2 =
-bh -  \frac{h}{2V}(\beta - d)^2
-\frac{h\nu}{2s^{-2}} - \frac{h}{2}(\beta - \hat{\beta})^2\sum_{i = 1}^N x_i^2
$$

Equating terms in $\beta^2$

$$
\frac{1}{V'} = \frac{1}{V} + \sum_{i = 1}^N x_i^2
$$

and so

$$
V' = \bigg(V^{-1} + \sum_{i = 1}^N x_i^2\bigg)^{-1}
$$

Equating terms in $\beta$

$$
\frac{d'}{V'} = \frac{d}{V} + \hat{\beta}\sum_{i = 1}^N x_i^2
$$

and so

$$
d' = V'\bigg(dV^{-1} + \hat{\beta}\sum_{i = 1}^N x_i^2\bigg)
$$

Equating constant terms

$$
b' + \frac{d'^2}{2V'} =
b + \frac{d^2}{2V} + \frac{\nu}{2s^{-2}} +
\frac{\hat{\beta}^2\sum_{i = 1}^N x_i^2}{2}
$$

and so

$$
b' =
b + \frac{1}{2}
\bigg(\frac{d^2}{V} - \frac{d'^2}{V'} +
\frac{\nu}{s^{-2}} +
\hat{\beta}^2\sum_{i = 1}^N x_i^2\bigg)
$$

But we know that

$$
\hat{\beta} = \frac{\sum x_i y_i}{\sum x^2_i}\, \text{and}\, s^2 = \frac{\sum (y_i - \hat{\beta}x_i)^2}{\nu}
$$

Thus we can rewrite

$$
\begin{aligned}
\frac{\nu}{s^{-2}} +
\hat{\beta}^2\sum_{i = 1}^N x_i^2 &=
\sum (y_i - \hat{\beta}x_i)^2 + \hat{\beta}^2\sum x_i^2 \\
&= \sum y_i^2 - 2\hat{\beta}\sum x_i y_i + \hat{\beta}^2\sum x_i^2 + \hat{\beta}^2\sum x_i^2 \\
&= \sum y_i^2
\end{aligned}
$$

giving

$$
b' =
b + \frac{1}{2}
\bigg(\frac{d^2}{V} - \frac{d'^2}{V'} +
\sum y_i^2\bigg)
$$

> {-# OPTIONS_GHC -Wall                      #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults    #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods  #-}
> {-# OPTIONS_GHC -fno-warn-orphans          #-}

> {-# LANGUAGE NoMonomorphismRestriction #-}

> module Main where
>
> import qualified Data.Vector.Unboxed as V
> import Data.Random.Source.PureMT
> import Data.Random
> import Control.Monad.State
>
> import Data.Colour
> import Control.Lens hiding ( (#) )
> import Graphics.Rendering.Chart hiding ( translate )
> import Diagrams.Backend.Cairo.CmdLine
>
> import Diagrams.Backend.CmdLine
> import Diagrams.Prelude hiding ( sample, render )
>
> import Data.Default.Class
>
> import Graphics.Rendering.Chart.Backend.Diagrams
>
> import System.IO.Unsafe
>
> nSamples :: Int
> nSamples = 10
>
> d :: Double
> d = 2.0
>
> s2InvPrior :: Double
> s2InvPrior = 1.0

> nu :: Double
> nu = 1.0
>
> v :: Double
> v = 1.0
>
> a :: Double
> a = v / 2.0
>
> b :: Double
> b = v / (2.0 * recip v**2)
>
> testXs :: Int ->V.Vector Double
> testXs m =
>   V.fromList $
>   evalState (replicateM m (sample StdUniform))
>   (pureMT 2)
>
> testXs' :: Int ->V.Vector Double
> testXs' m =
>   V.fromList $
>   evalState (replicateM m (sample StdUniform))
>   (pureMT 3)
>
> testEpsilons :: Int ->V.Vector Double
> testEpsilons m =
>   V.fromList $
>   evalState (replicateM m (sample StdNormal))
>   (pureMT 2)
>
> testEpsilons' :: Int ->V.Vector Double
> testEpsilons' m =
>   V.fromList $
>   evalState (replicateM m (sample StdNormal))
>   (pureMT 3)
>
> testYs :: Int ->V.Vector Double
> testYs m = V.zipWith (\x e -> d * x + e)
>            (testXs m) (testEpsilons m)
>
> testYs' :: Int ->V.Vector Double
> testYs' m = V.zipWith (\x e -> d * x + e)
>            (testXs' m) (testEpsilons' m)
>
> prices :: [(Double,Double)]
> prices = V.toList $ V.zip (testXs nSamples) (testYs nSamples)
>
> prices' :: [(Double,Double)]
> prices' = V.toList $ V.zip (testXs' nSamples) (testYs' nSamples)
>
> chart :: Colour Double ->[(Double, Double)] ->Graphics.Rendering.Chart.Renderable ()
> chart c prices = toRenderable layout
>   where
>
>     price1 = plot_points_style . point_color .~ opaque c
>            $ plot_points_values .~ prices
>            $ plot_points_title .~ "price 1"
>            $ def
>
>     layout = layoutlr_title .~"Price History"
>            $ layoutlr_left_axis . laxis_override .~ axisGridHide
>            $ layoutlr_right_axis . laxis_override .~ axisGridHide
>            $ layoutlr_x_axis . laxis_override .~ axisGridHide
>            $ layoutlr_plots .~ [Left (toPlot price1),
>                                 Right (toPlot price1)]
>            $ layoutlr_grid_last .~ False
>            $ def
>
> displayHeader :: FilePath -> Diagram B R2 -> IO ()
> displayHeader fn =
>   mainRender ( DiagramOpts (Just 900) (Just 600) fn
>              , DiagramLoopOpts False Nothing 0
>              )
>
> denv :: DEnv
> denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 500 500
>
> diag :: Colour Double -> [(Double, Double)] -> QDiagram Cairo R2 Any
> diag c prices = fst $ runBackend denv (render (chart c prices) (500, 500))
>
> a' = a + (fromIntegral nSamples) / 2.0
>
> xs = testXs nSamples
> ys = testYs nSamples
>
> xs2 = V.sum $ V.map (**2) xs
> ys2 = V.sum $ V.map (**2) ys
>
> v' = recip (recip v + xs2)
>
> betaHat = (V.sum $ V.zipWith (*) xs ys) / xs2
>
> d' = v' * (d * recip v + betaHat * xs2)
>
> b' = b + 0.5 * (d**2 / v - d'**2 / v' + ys2)
>
> main :: IO ()
> main = do
>   displayHeader "TestInteractive.png"
>     ((diag red prices # scaleX 0.4 # scaleY 0.4 # translate (r2 (0.1, 0.0)))
>      <> (diag blue prices' # scaleX 0.6 # scaleY 0.6 # translate (r2 (0.2, 0.0))))
>   putStrLn "Hello"
