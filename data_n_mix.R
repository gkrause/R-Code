################################################################
##  FUNCTION TO GENERATE DATA WITH ADDIDIVE NOISE SEQUENCE    ##
##                 AS MIXED GAUSSIAN DISTRIBUTION             ## 
################################################################

#############Description of Input-Variables######################
# numsamp = length of time series; numpar = number of parameters#
# a,b = vector of amplitudes; om = vector of frequencies        #  
# stdev1 = Standard Deviation of gaussian distrib Nr. 1         #
# stdev2 = Standard Deviation of gaussian distrib Nr. 2         #
#################################################################


data_n_mix <- function(numsamp,numpar,a,b,om, stdev1,stdev2){
  
  y = c(rep(0,numsamp))
  x = c(rep(0,numsamp))
  #Generation of empty vectors that be filled in for-loops  
  
  for(j in 1:numsamp){
    
    for(i in 1:numpar){
      
      x[j] =x[j]+(a[i]*cos(om[i]*j)+b[i]*sin(om[i]*j))    
      #vector of data without noise
    }
    bin = rbinom(1,1,0.6)
    y[j] = x[j] + ((bin*rnorm(1,0,stdev1)) +  (1-bin)*rnorm(1,0,stdev2))
    #data with additive noise sequence
  }
  return(y) 
}


