#' ---
#' title: Lab 6 solutions
#' author: 
#' date: 19-20th February 2018
#' ---
#' 
#' ## Exercise 1
#'
Env = read.table("http://www.massey.ac.nz/~jcmarsha/rcourse/EnvData.txt",header=T)
summary(Env)
#' A plot of the DANT station can be done using
TempDant <- Env[Env$Station == "DANT" ,]
plot(TempDant$MyTime,TempDant$T,type="l", xlab="Date", ylab="Temperature", main="DANT")
#' We can loop through them all using
par(mfrow=c(2,2))
AllStations = unique(Env$Station)
N = length(AllStations)
for (i in 1:N) {
  Station.i = as.character(AllStations[i])
  print(Station.i)
  TPi = Env[Env$Station == Station.i,]
  plot(TPi$MyTime, TPi$T, type="l", xlab="Date", ylab="Temperature", main=Station.i)
}
#'
#' ## Exercise 2
#' 
#' The following function counts the number of NA fields in a data frame.
NAcount = function(DF) {
  D1 = is.na(DF)
  x = colSums(D1)
  return(x)
}
NAcount(Env)
#'
#' ## Exercise 3
#'
#' This function plots a circle.
circ = function(...) {
  r = seq(0, 2*pi, length=1000)
  xCoords = sin(r)
  yCoords = cos(r)
  plot(xCoords, yCoords, type="l", ...)
}
par(mfrow=c(1,2))
circ()
circ(col="red", asp=1)
#' We can alter it to change the radius and coordinates by adding additional parameters
circ = function(radius=1, centerx=0, centery=0, ...) {
  r = seq(0, 2*pi, length=1000)
  xCoords = centerx + radius*sin(r)
  yCoords = centery + radius*cos(r)
  plot(xCoords, yCoords, type="l", ...)
}
circ(0.5, 1, 2, col="green", asp=1)
