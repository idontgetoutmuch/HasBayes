{-# LANGUAGE ScopedTypeVariables, NoMonomorphismRestriction #-}

import Data.Random
import Data.RVar
import Control.Monad.Trans as MTL
import System.Random
import Control.Monad.State
import System.Random.MWC (create)

import Debug.Trace

logNormal :: Double -> Double -> RVar Double
logNormal mu sigmaSq = do
    x <- normal mu sigmaSq
    return (exp x)

main = do
    mwc <- create
    y <- sampleFrom mwc (logNormal 5 1)
    print y

rwalkState :: RVarT (State Double) Double
rwalkState = do
  prev <- MTL.lift get
  change  <- rvarT StdNormal

  let new = prev + change
  MTL.lift (put new)
  return new

rwalk :: Int -> Double -> StdGen -> ([Double], StdGen)
rwalk count start gen =
  flip evalState start .
  flip runStateT gen .
  sampleRVarTWith MTL.lift $
  replicateM count rwalkState

n = 11112
seed = 0
m = 7

proposedJumps :: [Int]
proposedJumps = map f $ fst $ runState (replicateM n (sampleRVar (uniform False True))) (mkStdGen seed)
  where f False = negate 1
        f True  = 1

acceptOrRejects :: [Double]
acceptOrRejects = fst $ runState (replicateM n (sampleRVar (uniform 0 1))) (mkStdGen seed)

f currentPosition (proposedJump, acceptOrReject) =
  if acceptOrReject < probAccept
    then
      currentPosition + proposedJump
    else
      currentPosition
  where probAccept = min 1 (fromIntegral (p $ currentPosition + proposedJump) / fromIntegral (p currentPosition))

p n | n >=1 && n <= m = n
    | otherwise       = 0

main' = map (\j -> length $ filter (==j) $ drop (n `div` 10) $ scanl f 3 (zip proposedJumps acceptOrRejects)) [1..m]
-- [323,655,915,1309,1627,1901,2271]

likelihood :: Int -> Int -> Double -> Double
likelihood z n theta = theta^z * (1 - theta)^(n - z)

normals :: [Double]
normals =  fst $ runState (replicateM n (sampleRVar (normal 0 0.1))) (mkStdGen seed)

myData = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]

q :: Double -> [Double] -> Double
q position xs = likelihood z n position
  where
    n = length xs
    z = length $ filter (== 1) xs

g currentPosition (proposedJump, acceptOrReject) =
  if acceptOrReject < probAccept
    then
      currentPosition + proposedJump
    else
      currentPosition
  where
    probAccept = min 1 (p (currentPosition + proposedJump) / p currentPosition)
    p x | x < 0 = 0
        | x > 1 = 0
        | otherwise = f x
    f = flip q myData

main'' = drop (n `div` 10) $ scanl g 0.5 (zip normals acceptOrRejects)