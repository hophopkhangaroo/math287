library(tuneR)
library(readr)
library(ProDenICA)
library(whitening)
library(gridBase)
library(gridExtra)
library(ggplot2)
install.packages("remotes")
remotes::install_github("DesiQuintans/desiderata")
library(desiderata)
 
setwd("H:/My Drive/MATH287A")

set.seed(123456)

# Read audio files
whale = readWave("whale5s.wav")
ship = readWave("ship5s.wav")
freq = readWave("freq5s.wav")

# The files are stereo but we only need 1 vector for each sound source
whale_mat = whale@left
ship_mat = ship@left
freq_mat = freq@left

# Our source matrix
S = cbind(whale_mat, ship_mat, freq_mat)

# Our mixing matrix. The rows need not add up to 1
A = matrix(c(0.1, 0.5, 0.3,
             0.3, 0.6, 0.9,
             1, 1, 1), byrow=T, nrow=3)

X = S%*%A

# Example of mixed audio
play(Wave(X[,1], samp.rate = 44100, bit = 16))

# Centralize our observation vector X and make cov(X)=I
X_white = whiten(X, center=TRUE, method="PCA")

# Perform Projection Density ICA
results = ProDenICA(X_white, k=3, Gfunc=GPois,trace=TRUE, density=TRUE)

#FastICA
fastica = ProDenICA(X_white, k=3, trace=TRUE,Gfunc=G1)

#Normalizing wave files to perceptible levels 
result1 = results$s[,1]*32767/max(abs(results$s[,1]))
result2 = results$s[,2]*32767/max(abs(results$s[,2]))/10 #lowering vol bc I know this is the freqency noise
result3 = results$s[,3]*32767/max(abs(results$s[,3]))

fastica1 = fastica$s[,1]*32767/max(abs(fastica$s[,1]))
fastica2 = fastica$s[,2]*32767/max(abs(fastica$s[,2])) / 10
fastica3 = fastica$s[,3]*32767/max(abs(fastica$s[,3]))

# Formatting for playback
sound1 = Wave(result1, samp.rate = 44100, bit = 16)
sound2 = Wave(result2, samp.rate = 44100, bit = 16)
sound3 = Wave(result3, samp.rate = 44100, bit = 16)
sound4 = Wave(fastica1, samp.rate = 44100, bit = 16)
sound5 = Wave(fastica2, samp.rate = 44100, bit = 16)
sound6 = Wave(fastica3, samp.rate = 44100, bit = 16)



play(sound1) # whale 
play(sound2) # frequency
play(sound3) # ship
play(sound4) # whale
play(sound5) # freq
play(sound6) # ship

# Amari Distance
amari(solve(results$W), ICAorthW(solve(A)))
amari(solve(fastica$W), ICAorthW(solve(A)))


plot_arrange(plot(whale_mat,type="l", main="Whale"),
plot(ship_mat,type="l", main="Ship"),
plot(freq_mat,type="l", main="Frequency"),
plot(X[,1], type="l", main="Mixed Signal 1"),
plot(X[,2], type="l", main="Mixed Signal 2"),
plot(X[,3], type="l", main="Mixed Signal 3"),
plot(result1, type="l", main="Extracted Whale Sound"),
plot(result2, type="l", main = "Extracted Freqency Sound"),
plot(result3, type="l", main = "Extracted Ship Sound"), nrow=3, ncol=3)
