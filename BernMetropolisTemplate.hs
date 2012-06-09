{-# LANGUAGE ScopedTypeVariables, NoMonomorphismRestriction #-}

import System.Random.MWC.Monad
import System.Random.MWC.Distributions.Monad
import Control.Monad.Primitive
import Control.Applicative
import Debug.Trace

myData = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]

likelihood theta xs = theta^z * (1-theta)^(n-z)
  where
    z = sum $ filter (== 1) xs
    n = length xs

-- Define the prior density function. For purposes of computing p(D),
-- at the end of this program, we want this prior to be a proper density.
prior theta = runWithCreate $ uniformR (0,1)
prior' theta = uniformR (0,1)

-- TBD: Use a bimodal prior.

-- Define the relative probability of the target distribution,
-- as a function of vector theta. For our application, this
-- target distribution is the unnormalized posterior distribution.
targetRelProb theta xs = do
  y <- prior theta
  trace (show y) $ return $ y * likelihood theta xs

targetRelProb' theta xs = do
  y <- prior' theta
  trace (show y) $ return $ y * likelihood theta xs

-- Specify the length of the trajectory, i.e., the number of jumps to try:
trajLength = 10 -- 11112 -- arbitrary large number

-- Specify the burn-in period:
burnIn = ceiling ( 0.1 * fromIntegral trajLength ) -- arbitrary number, less than trajLength

-- Use the proposal distribution (normal for this example) to generate a proposed jump.
proposedJumps :: PrimMonad m => Rand m [Double]
proposedJumps = sequence $ take trajLength $ repeat $ normal 0.0 0.1

acceptOrRejects :: PrimMonad m => Rand m [Double]
acceptOrRejects = sequence $ take trajLength $ repeat $ uniformR (0,1)

newTrajectory :: IO [Double]
newTrajectory = do
  xs :: [Double] <- runWithCreate $ proposedJumps
  ys :: [Double] <- runWithCreate $ acceptOrRejects
  scanlM f 0.5 (zip xs ys)

-- newTrajectory' :: IO [Double]
newTrajectory' = do
  xs <- proposedJumps
  ys <- acceptOrRejects
  scanlM f' 0.5 (zip xs ys)

scanlM :: Monad m => (a -> b -> m a) -> a -> [b] -> m [a]
scanlM f q ls = do
  ms <- case ls of
          []   -> return []
          x:xs -> do y <- f q x
                     scanlM f y xs
  return $ q : ms

f :: Double -> (Double, Double) -> IO Double
f currentPosition (proposedJump, acceptOrReject) = do
  -- Compute the probability of accepting the proposed jump.
  x <- targetRelProb (currentPosition + proposedJump) myData
  y <- targetRelProb currentPosition myData
  let probAccept = min 1 (x / y)
  if acceptOrReject < probAccept
    then
      return $ currentPosition + proposedJump
    else
      return currentPosition

f' currentPosition (proposedJump, acceptOrReject) = do
  -- Compute the probability of accepting the proposed jump.
  x <- targetRelProb' (currentPosition + proposedJump) myData
  y <- targetRelProb' currentPosition myData
  let probAccept = min 1 (x / y)
  if acceptOrReject < probAccept
    then
      return $ currentPosition + proposedJump
    else
      return currentPosition


acceptedTraj = (drop burnIn) <$> newTrajectory
