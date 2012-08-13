trajLength = 11112

myData = c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 )

set.seed(47405)

normals = rep (0, trajLength)
for (j in 1:trajLength) { normals[j] = rnorm(1, mean = 0, sd = 0.1) }

proposedJumps = rep (0, trajLength)
for (j in 1:trajLength) { proposedJumps[j] = if ( runif( 1 ) < 0.5 ) { -1 } else { 1 }}

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

nsAorBs <- list ()
for (i in 1:trajLength) nsAorBs[[i]] <- c(normals[i], acceptOrRejects[i])

trajectory = Reduce(function(a,b) oneStep (a, b), nsAorBs, accumulate=T,init= 0.5 )

burnIn = ceiling( .1 * trajLength )

acceptedTraj = trajectory[burnIn:trajLength]

result = mean ( acceptedTraj )
