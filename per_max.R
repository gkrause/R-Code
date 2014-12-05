####   PERIODOGRAM FUNCTION  ##
###############################

#####################################################################################################
# Used to find starting values for optimization in non-linear least square method                  ##
# (see: "Sequential estimation of the sum of sinusoidal model parameters"; Prasad, Kundu & Mitra ) ##
#####################################################################################################

#############Description of Input-Variables######################
## yadj =Daten - Zeitreihe; numsamp = Laenge der Zeitreihe      ##
## gridnum = Anzahl der grid points                            ##
#################################################################

per_max <- function(yadj, numsamp, gridnum){
  
  start_om = 0
  end_om= pi
   
  
  
  for(k in 1:gridnum+1){
    om_grid[k] = start_om+((end_om-start_om)/gridnum)*(k-1)
     
    
      for(kk in 1:numsamp){
        per[k] = per[k] + yadj[kk]*exp(complex(1,0,-1)*kk*om_grid[k])
         
      }
    
    per[k] = (per[k]/numsamp)
    per[k] = (abs(per[k])^2)
     #Berechnung des Periodograms
  }
  I = which.max(per)    
  om = om_grid[I]
   #Maxima des Periodograms
  
          for( kk in 1:numsamp){
           a_om[kk,1] = cos(om*kk)
           a_om[kk,2] = sin(om*kk)
           yvec[kk] = yadj[kk]
            #Definiere Startwerte
           
          }
  
lin = (ginv((t(a_om)%*%a_om)))%*%t(a_om)%*%yvec
 #Berechnung der Startwerte mittel Pseudoinversen

inits[1] = lin[1]
inits[2] = lin[2]
inits[3] = om

return(inits)

}
