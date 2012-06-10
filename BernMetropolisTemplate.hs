{-# LANGUAGE ScopedTypeVariables, NoMonomorphismRestriction #-}

import System.Random.MWC.Monad
import System.Random.MWC.Distributions.Monad
import Control.Applicative
import Control.Monad.State
import Control.Monad.ST

myData = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]

-- Define the prior density function. For purposes of computing p(D),
-- at the end of this program, we want this prior to be a proper density.
prior theta = uniformR (0,1)

-- TBD: Use a bimodal prior.

likelihood theta xs = theta^z * (1-theta)^(n-z)
  where
    z = sum $ filter (== 1) xs
    n = length xs

-- Define the relative probability of the target distribution,
-- as a function of vector theta. For our application, this
-- target distribution is the unnormalized posterior distribution.
targetRelProb theta xs = do
  y <- prior theta
  return $ y * likelihood theta xs

-- Specify the length of the trajectory, i.e., the number of jumps to try:
trajLength = 111112 -- arbitrary large number

-- Specify the burn-in period:
burnIn = ceiling ( 0.1 * fromIntegral trajLength ) -- arbitrary number, less than trajLength

-- Use the proposal distribution (normal for this example) to generate a proposed jump.
proposedJumps = sequence $ take trajLength $ repeat $ normal 0.0 0.1

acceptOrRejects = sequence $ take trajLength $ repeat $ uniformR (0,1)

stepOnce (proposedJump, acceptOrReject) = do
  currentPosition : oldPositions <- get
  x <- lift $ targetRelProb (currentPosition + proposedJump) myData
  y <- lift $ targetRelProb currentPosition myData
  let probAccept = min 1 (x / y)
  if acceptOrReject < probAccept
    then
      put $ (currentPosition + proposedJump) : currentPosition : oldPositions
    else
      put $ currentPosition : currentPosition : oldPositions
  return ()

stepOnce' = do
  currentPosition : oldPositions <- get
  proposedJump <- lift $ normal 0.0 0.1
  acceptOrReject <- lift $ uniformR (0,1)
  x <- lift $ targetRelProb (currentPosition + proposedJump) myData
  y <- lift $ targetRelProb currentPosition myData
  let probAccept = min 1 (x / y)
  if acceptOrReject < probAccept
    then
      put $ (currentPosition + proposedJump) : currentPosition : oldPositions
    else
      put $ currentPosition : currentPosition : oldPositions
  return ()

newTrajectory = do
  xs <- lift proposedJumps
  ys <- lift acceptOrRejects
  mapM_ stepOnce (zip xs ys)

-- newTrajectory' = do
--   mapM_ stepOnce'

acceptedTraj = (drop burnIn . snd) <$> runStateT newTrajectory [0.5]

foo = runWithCreate $ asRandST acceptedTraj
-- main = do
--   ys <- runWithCreate acceptedTraj
--   putStrLn $ show $ ys