#library(ggplot2)


#http://www.ms314u.me/2018/03/02/monitoring-signals-deviations-from-the-normal-behavior/

getCoeffs <- function(y){
  N <- length(y)
  A <- rep(0, (floor(N/2)+1))
  B <- rep(0, (floor(N/2)+1))
  for(p in 1:(floor(N/2)+1)){
    for(k in 1:N){
      A[p] = A[p] + (2/N)*y[k]*cos(2*pi*(p-1)*k/N)
      B[p] = B[p] + (2/N)*y[k]*sin(2*pi*(p-1)*k/N)
    }
  }
  A[floor(N/2)+1] <- A[floor(N/2)+1]/2
  return(list(A,B))
}

predictSignal<-function(origSig, maxComponents, horizon){
  fCoeffs = getCoeffs(origSig)
  N = length(origSig)
  A = fCoeffs[[1]]
  B = fCoeffs[[2]]
  fourierEstim <- rep(0,N+horizon) + A[1]/2 
  t_moments = seq(1, N + horizon)
  for(p in 2:maxComponents){
    fourierEstim = fourierEstim + A[p]*cos(2*pi*(p-1)*t_moments/N) + B[p]*sin(2*pi*(p-1)*t_moments/N)
  }

  return(fourierEstim)
}

n_components = 5#5
dat <- sin(seq(1,100)* 2*pi/30) + rnorm(100)*0.1
predictions <- predictSignal(dat, n_components, 30)
dat <- c(dat, rep(NA, 30))

plot(predictions,type='l', cex=.35, xlab = "t", ylab = "Y", col="red")
lines(dat, col="black")

#t = seq(1, length(dat))
#pdat = data.frame(x=t, y=dat, y2 = predictions)
#p <- ggplot(data = pdat, aes(x=x, y=y)) + theme_bw()
#p <- p + xlab("Time") + ylab("Value")
#p <- p + geom_line()
#p <- p + geom_line(aes(y=predictions), colour = "red")
#p

#lines(dat - predictions)

# Reconstruct Signal - pmax is number of frequencies used in increasing order
#pmax <- 30
#fourierEstim <- rep(0,N)
#for(n in 1:N){
#  fourierEstim[n] <- A[1]/2
#  for(p in 2:pmax){
#    fourierEstim[n] <- fourierEstim[n] + A[p]*cos(2*pi*(p-1)*n/N) + B[p]*sin(2*pi*(p-1)*n/N)
#  }
#}
#plot(dat,type='l', cex=.35, xlab = "t", ylab = "Y")
#lines(fourierEstim, col="red")
#lines(dat - fourierEstim)

#for(n in 1:N){
#  fourierEstim[n] <- A[1]/2
#  for(p in 2:maxComponents){
#    fourierEstim[n] <- fourierEstim[n] + A[p]*cos(2*pi*(p-1)*n/N) + B[p]*sin(2*pi*(p-1)*n/N)
#  }
#}