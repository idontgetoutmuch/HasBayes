{-# LANGUAGE ScopedTypeVariables, NoMonomorphismRestriction #-}

import System.Random.MWC.Monad
import System.Random.MWC.Distributions.Monad
import Control.Applicative
import Control.Monad.State
import Debug.Trace

myData = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]

-- Define the prior density function. For purposes of computing p(D),
-- at the end of this program, we want this prior to be a proper density.
-- prior theta = uniformR (0,1)

-- TBD: Use a bimodal prior.

-- Define the relative probability of the target distribution,
-- as a function of vector theta. For our application, this
-- target distribution is the unnormalized posterior distribution.
-- targetRelProb theta xs = do
--   y <- prior theta
--   trace (show y) $ return $ y * likelihood theta xs

-- Specify the length of the trajectory, i.e., the number of jumps to try:
trajLength = 3 -- 11112 -- arbitrary large number

-- Specify the burn-in period:
burnIn = ceiling ( 0.1 * fromIntegral trajLength ) -- arbitrary number, less than trajLength

-- Use the proposal distribution (normal for this example) to generate a proposed jump.
proposedJumps = sequence $ take trajLength $ repeat $ normal 0.0 0.1

acceptOrRejects = sequence $ take trajLength $ repeat $ uniformR (0,1)

newTrajectory = do
  xs <- proposedJumps
  ys <- acceptOrRejects
  scanlM f 0.5 (zip xs ys)

scanlM :: Monad m => (a -> b -> m a) -> a -> [b] -> m [a]
scanlM f q ls = do
  ms <- case ls of
          []   -> return []
          x:xs -> do y <- f q x
                     scanlM f y xs
  return $ q : ms

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

g (proposedJump, acceptOrReject) = do
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

newTrajectory' = do
  xs <- lift proposedJumps
  ys <- lift acceptOrRejects
  mapM_ g (zip xs ys)

prior theta = uniformR (0,1)

targetRelProb theta xs = do
  y <- prior theta
  trace (show y) $ return $ y * likelihood theta xs

likelihood theta xs = theta^z * (1-theta)^(n-z)
  where
    z = sum $ filter (== 1) xs
    n = length xs

acceptedTraj = (drop burnIn) <$> newTrajectory

acceptedTraj' = (drop burnIn . snd) <$> runStateT newTrajectory' [0.5]

main = do
  xs <- runWithCreate acceptedTraj
  putStrLn $ show $ sum xs
  ys <- runWithCreate acceptedTraj'
  putStrLn $ show $ sum ys