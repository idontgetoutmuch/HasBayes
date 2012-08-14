## An arbitrary number of jumps to try in Metropolis algorithm.
trajLength = 11112

## The data represented as a vector.
myData = c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 )

## An arbitrary seed presumably used in generating random values from
## the various distributions.
set.seed(47405)

## We need some proposed jumps; following the example we generate
## these from a normal distribution.
normals = rep (0, trajLength)
for (j in 1:trajLength) { normals[j] = rnorm(1, mean = 0, sd = 0.1) }

## We generate proposed jumps which are either one step up or one step
## down.
proposedJumps = rep (0, trajLength)
for (j in 1:trajLength) { proposedJumps[j] = if ( runif( 1 ) < 0.5 ) { -1 } else { 1 }}

## We generate samples from the uniform distribution on [0, 1]
## which will allow us to determine whether to accept or reject a
## proposed jump.
acceptOrRejects = rep (0, trajLength)
for (j in 1:trajLength) { acceptOrRejects[j] = runif( 1 ) }

likelihood = function( z, n, theta ) {
  x = theta^z * (1 - theta)^(n - z)
  x[ theta > 1 | theta < 0 ] = 0
  return ( x )
}

prior = function( position, xs ) {
  n = length( xs )
  z = sum ( xs == 1)
  return ( likelihood ( z, n, position ) )
}

## We define a function which defines one step of the Metropolis algorithm.
oneStep = function ( currentPosition, proposedJumpAcceptOrReject ) {
  proposedJump   = proposedJumpAcceptOrReject [1]
  acceptOrReject = proposedJumpAcceptOrReject [2]
  probAccept = min( 1,
    prior( currentPosition + proposedJump , myData )
    / prior( currentPosition , myData ) )
  if ( acceptOrReject < probAccept ) {
    trajectory = currentPosition + proposedJump
  } else {
    trajectory = currentPosition
  }
  return ( trajectory )
}

## Pair together the proposed jumps and the probability whether a
## given proposal will be accepted or rejected.
nsAorRs <- list ()
for (i in 1:trajLength) nsAorRs[[i]] <- c(normals[i], acceptOrRejects[i])

## Fold (well really scanl) over the pairs returning all the
## intermediate results.
trajectory = Reduce(function(a,b) oneStep (a, b), nsAorRs, accumulate=T,init= 0.5 )

## Drop the values for the burn in period.
burnIn = ceiling( .1 * trajLength )
acceptedTraj = trajectory[burnIn:trajLength]

## Finally return the mean of the posterior.
result = mean ( acceptedTraj )
