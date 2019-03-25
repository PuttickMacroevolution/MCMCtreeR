# cauchyMCMCtree - internal function

cauchyMCMCtree <- function(xRange=c(0, 2), tL=1, p=0.5, c=0.2, minProb=0.025, maxProb=0.975) {

	A <- 0.5 + (1/pi) * atan(p/c)
	omega <- (maxProb/minProb) * (1 / (pi * A * c * (1 + (p/c)^2)))
	xx <- seq(xRange[1], xRange[2], by=0.001)
	rc <- 1
	numero <- c()
	for(x in xx) {
		t <- x
		if(t > tL) numero[rc] <- maxProb * (1 / (A * pi * c * tL * (1 + ((t - tL * (1 + p)) / (c * tL)) ^ 2)))
		if(t <= tL) numero[rc] <- minProb * (omega/tL) * (t/tL) ^ (omega-1)
		rc <- rc + 1
		}
	time <- xx ; density <- numero
	return(cbind(time, density))
}


