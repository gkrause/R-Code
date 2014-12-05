################################################################################
##MONTE-CARLO SIMULATION TO GET LEAST SQUARE ESTIMATES, CALCULATE LIKELIHOODS ##
##AND DIFFERENT INFORMATION CRITERIAS ON M-ESTIMATOR BASIS FOR SINUSOIDAL     ##
##SIGNAL MODEL.                                                               ##
################################################################################
###R-Packages###
################################
#install.packages("MASS")
#install.packages("pracma")
require(MASS)
#Fuer per_max Funktion

require(pracma)
#Fuer Optimierung
################################
setwd("F://R//Thesis//Kapitel 3//R Code")
#Work Dir. auf Kapitel 3 - R Code

source("Definition_of_Variables.R")


source("per_max.R")
source("Obj_Andrew.R")
source("Obj_Huber.R")
source("Obj_Hampel.R")
source("obj_Least.R")
source("Obj_Ramsey.R")
source("data_n.R")
source("data_n_mix.R")
source("data_t.R")

#Begin der Monte-Carlo Simulation
for(sim in 1: 100){
  
  #Nullsetzen aller Variablen zur Berechnung
  yp = c(rep(0,numsamp))
  obh = c(rep(0,numsamp))
  seq_a = c(rep(0,para))
  seq_b = c(rep(0,para))
  seq_om = c(rep(0,para))
  a_mh = c(rep(0,para))
  b_mh = c(rep(0,para))
  om_mh = c(rep(0,para))
  a = c(rep(0,para))
  b = c(rep(0,para))
  om = c(rep(0,para))
  om_grid = c(rep(0,(gridnum+1)))
  per = c(rep(0,(gridnum+1)))
  a_om = matrix(0,numsamp,2)
  yvec = c(rep(0,numsamp))
  inits = c(0,0,0)
  a_m = c(rep(0,para))
  b_m = c(rep(0,para))
  om_m = c(rep(0,para))
  diff_full_model = c(rep(0,numsamp))
  diff_med_dev = c(rep(0,numsamp))
  
  
#Generierung der Daten

#y = data_n(numsamp,numpar,a_o,b_o,om_o,stdev)
 #Zeitreihe mit Gaussian Normal Fehler

#y = data_n_mix(numsamp, numpar,a_o,b_o,om_o, stdev1,stdev2)
#y = data_n_mix(numsamp, numpar,a_o,b_o,om_o, stdev3,stdev4)
#y = data_n_mix(numsamp, numpar,a_o,b_o,om_o, stdev5,stdev6)
y = data_n_mix(numsamp, numpar,a_o,b_o,om_o, stdev7,stdev8)
 #Zeitreihe mit mixture Gaussian normal Fehler

#y = data_t(numsamp,numpar,a_o,b_o,om_o,dft)
 #Zeitreihe mit Student´s t Fehler

########################################################################
#Speichern der Zeitreihen in Matrix

#data_gaussian_normal_m[sim,] = y
#data_gaussian_normal_mix_01_001_m[sim,] = y
#data_gaussian_normal_mix_06_001_m[sim,] = y
#data_gaussian_normal_mix_2_01_m[sim,] = y
data_gaussian_normal_mix_3_01_m[sim,] = y
#data_gaussian_student_t_m[sim,] = y


 #Start der Schaetzung von robustem Sigma

yadj = y
#Definition von yadj als Daten, da diese im Laufe der for-Schleife veraendert werden, die Zeitreihe aber immer gleich bleiben sollte


for(ks in 1:para){
    inits = per_max(yadj,numsamp,gridnum)
      #Startwerte fuer Optimierung aus Periodogram 
  
    seq_a[ks] = inits[1]
	  seq_b[ks] = inits[2]
	  seq_om[ks] = inits[3]
  
	for(ns in 1:numsamp){
    
		yadj[ns] = (yadj[ns] - inits[1]*cos(inits[3]*ns)-inits[2]*sin(inits[3]*ns))
		 #Die Startwerte fuer Alpha, Beta und Omega werden ueber die Maxima des Periodograms ermittelt. Da es hier
     #mehrere Maxima gibt wird ueber Differenzieren versucht diese Step-by-Step fuer alle Komponenten = para zu
     #errechnen.
     
	}
}

		vect <- c(seq_a,seq_b,seq_om)
		#Vektor der Startwerte, notwenig fuer fminsearch-Befehl
		
		mm_esta_sig = fminsearch(obj_hub, vect, comp = para)
		#Optimierung/Minimierung von Huber´s Funktion zur Berchnung des robusten Sigmas. Es handelt sich hier um eine
    # M-Schaetzung im klassischen Sinne. Der fminsearch-Befehl vom Paket pracma ist eigentlich ein Nelder-Mead Algorithmus.

		m_est_sig = mm_esta_sig$xval
		#Ergebnis der M-Schaetzung --> Schaetzer
		
			a_m[1:para] = m_est_sig[1:para]
			b_m[1:para] = m_est_sig[(para+1):(2*para)]
			om_m[1:para]= m_est_sig[(2*para+1):(3*para)]
			#Definition der Schaetzer nach Logik des Models
			
				for (k in 1:numsamp){
  
					for(kp in 1:para){
						yp[k] = yp[k] + (a_m[kp]*cos(om_m[kp]*k)+b_m[kp]*sin(om_m[kp]*k))
						#Erstellung des Models mit den ermittelten Schaetzern
     
					}
				diff_full_model[k] = y[k] - yp[k] 
				}
		med_diff = median(diff_full_model)

						for (k in 1:numsamp){
       
							diff_med_dev[k] = abs(diff_full_model[k] - med_diff)
						}
		#Berechnung des robusten Sigmas nach Mitra & Mitra (15)

		sig_est = 1.483*median(diff_med_dev)
		#Robustes Sigma
##################################################################################

#Anwendung der Huber-methode

for (comp in 1:para){
	#Beginn der for-Schleife fuer Berechnung aller Informations Krieterien fuer p = (1,...,6); Schleife laeuft
  #ueber comp mit para = 6
  
		nummax = comp
		yadj = y
		
		 if(comp > 1){
       #sobald die Schleife einmal gelaufen ist, sollen alle Schaetzer geloescht werden um keinen Konflikt mit neuen
       # Schaetzern von Alpha, Beta und Omega fuer p>1 zu erzielen.
		 rm(m_est)}
		 
		yp = c(rep(0,numsamp))
		obh = c(rep(0,numsamp))
		
		seq_a = c(rep(0,comp))
		seq_b = c(rep(0,comp))
		seq_om = c(rep(0,comp))
		 #Nullsetzen der Vektoren in jeder Iteration fuer Komponenten -> sonst kommt es zu Rechenfehlern, 
     #da R alte und neue Werte vermischt. 
    
		for(ks in 1:nummax){
		  #Startwerte fuer M-Schaetzung/Minimierung ueber Periodogram 
		  
		  inits = per_max(yadj,numsamp,gridnum)
			seq_a[ks] = inits[1]
			seq_b[ks] = inits[2]
			seq_om[ks] = inits[3]
		  
		  for(ns in 1:numsamp){
		    yadj[ns] = (yadj[ns] - inits[1]*cos(inits[3]*ns)-inits[2]*sin(inits[3]*ns))
        # Da mehrere Maxima im Periodogram existieren muessen diese ueber Differeinzierung 
        # berechnet werden.
		  }
		}
		
		
		
		vect <- c(seq_a,seq_b,seq_om)
		#Startwerte fuer fminsearch - Befehl		
				
					mm_esta = fminsearch(obj_hub, vect, comp = comp)
          # M-Schaetzung; Anwendung des Nelder-Mead-Algorithmus ueber fminsearch-Befehl aus Paket pracma.
          # Ein wichtiger Grund fuer die Verwendung von 'fminsearch' und nicht von z.B. 'optim' oder einer
          # anderen Optimierung in R ist die Moeglichkeit ueber den Befehl Variablen der Objektiven Funktion 
          # (hier Huber´s Funktion) zu veraendern. Das ist wichtig fuer comp, da comp mit jeder Iteration um
          # +1 waechst. Andere Optimierungsbefehle in R benoetigen eine fixe Objektive Funktion.
    
    
					m_est = mm_esta$xval
					#Schaetzer
    
			if (comp == 1){
				a_hub_1 =m_est[1:comp]
				b_hub_1 =m_est[(comp+1):(2*comp)]
				om_hub_1 =m_est[(2*comp+1):(3*comp)]}
				if (comp ==2){
					a_hub_2 =m_est[1:comp]
					b_hub_2 =m_est[(comp+1):(2*comp)]
					om_hub_2 =m_est[(2*comp+1):(3*comp)]}
					if (comp ==3){
						a_hub_3 =m_est[1:comp]
						b_hub_3 =m_est[(comp+1):(2*comp)]
						om_hub_3 =m_est[(2*comp+1):(3*comp)]}
						if (comp ==4){
							a_hub_4 =m_est[1:comp]
							b_hub_4 =m_est[(comp+1):(2*comp)]
							om_hub_4=m_est[(2*comp+1):(3*comp)]}
							if (comp ==5){
								a_hub_5 =m_est[1:comp]
								b_hub_5 =m_est[(comp+1):(2*comp)]
								om_hub_5 =m_est[(2*comp+1):(3*comp)]}
								if (comp ==6){
									a_hub_6 =m_est[1:comp]
									b_hub_6 =m_est[(comp+1):(2*comp)]
									om_hub_6 =m_est[(2*comp+1):(3*comp)]}
    #Da fuer die Analysen (Konfidenzintervalle, Standardabweichung etc.) die jeweiligen Schaetzer benoetigt werden muessen
    #diese gesondert gespeichert werden. Natuerlich existieren andere Schaetzer je nach Wert von der Komponenten p = comp.
    
###
ob_hub = c(rep(0,numsamp))
differ = c(rep(0,numsamp))

		for(k in 1:numsamp){
      
			yp[k] = 0
			for(kp in 1:comp){
					if(comp == 1){
						yp[k] = yp[k]+(a_hub_1[kp]*cos(om_hub_1[kp]*k)+b_hub_1[kp]*sin(om_hub_1[kp]*k))
					}
						if(comp == 2){
							yp[k] = yp[k]+(a_hub_2[kp]*cos(om_hub_2[kp]*k)+b_hub_2[kp]*sin(om_hub_2[kp]*k))
						}
							if(comp == 3){
								yp[k] = yp[k]+(a_hub_3[kp]*cos(om_hub_3[kp]*k)+b_hub_3[kp]*sin(om_hub_3[kp]*k))
							}
								if(comp == 4){
									yp[k] = yp[k]+(a_hub_4[kp]*cos(om_hub_4[kp]*k)+b_hub_4[kp]*sin(om_hub_4[kp]*k))
								}
									if(comp == 5){
										yp[k] = yp[k]+(a_hub_5[kp]*cos(om_hub_5[kp]*k)+b_hub_5[kp]*sin(om_hub_5[kp]*k))
									}
										if(comp == 6){
											yp[k] = yp[k]+(a_hub_6[kp]*cos(om_hub_6[kp]*k)+b_hub_6[kp]*sin(om_hub_6[kp]*k))
										}
          #Je nach Wert von p = comp (Komponente) muss die geschaetzte Funktion/Model gebildet werden.
				}
			differ[k] = (y[k] - yp[k])/sig_est
				if(abs(differ[k])<=ht){
					ob_hub[k] = ob_hub[k] + ((differ[k])^2)/2
				}
				else{
         			ob_hub[k] = ob_hub[k]+abs(differ[k])*ht - (ht^2)/2
				}
			}
			 ob_hub_sum = sum(ob_hub)
       # Berechnung der Huber Funktion Likelihood fuer Informationskriterien

###
#Berechnung aller robusten Informationskriterien mit Huber´s Methode:

obj_hub_bic[comp] = ob_hub_sum +0.5*(log(numsamp))*(3*comp)
obj_hub_aic[comp] = ob_hub_sum+(3*comp)
obj_hub_aicc[comp] = ob_hub_sum + (3*comp*numsamp)/(numsamp-3*comp-1)
ca = (3*comp*numsamp)/(numsamp-3*comp-1)
cb = 0.5*(log(numsamp))*(3*comp)
obj_hub_wic[comp] = ((ca)/(ca+cb))*obj_hub_aicc[comp]+((cb)/(ca+cb))*obj_hub_bic[comp]
obj_hub_hq[comp] = ob_hub_sum + (3*comp)*log(log(numsamp))	

####
#Ramsey´s, Hampel´s, Andrew´s und die Nicht-Robusten Schaetzungen werden alle nach dem gleichen Prinzip
#erstellt wie Huber´s Methode hier. Deswegen wird die Schaetzung mithilfe von fminsearch sowie die Aufteilung 
#der Schaetzer nicht mehr zusaetzlich kommentiert.
#Der Unterschied liegt in der Objektiven Funktion ueber welche der M-Schaetzer ermittelt wird, sowie in der Berechnung
#der Likelihood nach Logik der jeweils angewendeten Methode.


##RAMSEY´s METHODE

		
		nummax = comp
		yadj = y
		
		 if(comp > 1){
		 rm(m_est)}
			yp = c(rep(0,numsamp))
			differ = c(rep(0,numsamp))
		 
			
			mm_esta = fminsearch(obj_ram, vect, comp = comp)
			m_est = mm_esta$xval
			if (comp == 1){
				a_ram_1 =m_est[1:comp]
				b_ram_1 =m_est[(comp+1):(2*comp)]
				om_ram_1 =m_est[(2*comp+1):(3*comp)]}
				if (comp ==2){
					a_ram_2 =m_est[1:comp]
					b_ram_2 =m_est[(comp+1):(2*comp)]
					om_ram_2 =m_est[(2*comp+1):(3*comp)]}
					if (comp ==3){
						a_ram_3 =m_est[1:comp]
						b_ram_3 =m_est[(comp+1):(2*comp)]
						om_ram_3 =m_est[(2*comp+1):(3*comp)]}
						if (comp ==4){
							a_ram_4 =m_est[1:comp]
							b_ram_4 =m_est[(comp+1):(2*comp)]
							om_ram_4=m_est[(2*comp+1):(3*comp)]}
							if (comp ==5){
								a_ram_5 =m_est[1:comp]
								b_ram_5 =m_est[(comp+1):(2*comp)]
								om_ram_5 =m_est[(2*comp+1):(3*comp)]}
								if (comp ==6){
									a_ram_6 =m_est[1:comp]
									b_ram_6 =m_est[(comp+1):(2*comp)]
									om_ram_6 =m_est[(2*comp+1):(3*comp)]}
###									
		ob_ram = c(rep(0,numsamp))
		differ = c(rep(0,numsamp))
		
		for(k in 1:numsamp){
      
			yp[k] = 0
			for(kp in 1:comp){
					if(comp == 1){
						yp[k] = yp[k]+(a_ram_1[kp]*cos(om_ram_1[kp]*k)+b_ram_1[kp]*sin(om_ram_1[kp]*k))
					}
						if(comp == 2){
							yp[k] = yp[k]+(a_ram_2[kp]*cos(om_ram_2[kp]*k)+b_ram_2[kp]*sin(om_ram_2[kp]*k))
						}
							if(comp == 3){
								yp[k] = yp[k]+(a_ram_3[kp]*cos(om_ram_3[kp]*k)+b_ram_3[kp]*sin(om_ram_3[kp]*k))
							}
								if(comp == 4){
									yp[k] = yp[k]+(a_ram_4[kp]*cos(om_ram_4[kp]*k)+b_ram_4[kp]*sin(om_ram_4[kp]*k))
								}
									if(comp == 5){
										yp[k] = yp[k]+(a_ram_5[kp]*cos(om_ram_5[kp]*k)+b_ram_5[kp]*sin(om_ram_5[kp]*k))
									}
										if(comp == 6){
											yp[k] = yp[k]+(a_ram_6[kp]*cos(om_ram_6[kp]*k)+b_ram_6[kp]*sin(om_ram_6[kp]*k))
										}
				}
		differ[k] = (y[k] - yp[k])/sig_est
		ob_ram[k] = ob_ram[k] + (ra^(-2))*(1-exp(-ra*abs(differ[k]))*(1+ra*abs(differ[k])))
}
		ob_ram_sum =sum(ob_ram)
###
#ROBUSTE INFORMATIONSKRITERIEN MITTELS RAMSEY´S METHODE
obj_ram_bic[comp] =  ob_ram_sum+0.5*(log(numsamp))*(3*comp)
obj_ram_aic[comp] = ob_ram_sum+(3*comp)
obj_ram_aicc[comp] = ob_ram_sum + (3*comp*numsamp)/(numsamp-3*comp-1)
ca = (3*comp*numsamp)/(numsamp-3*comp-1)
cb = 0.5*(log(numsamp))*(3*comp)
obj_ram_wic[comp] = ((ca)/(ca+cb))*obj_ram_aicc[comp]+((cb)/(ca+cb))*obj_ram_bic[comp]
obj_ram_hq[comp] = ob_ram_sum + (3*comp)*log(log(numsamp))

####ANDREW´S METHODE	

		
		nummax = comp
		yadj = y
		
		 if(comp > 1){
		 rm(m_est)}
		 
			yp = c(rep(0,numsamp))
			differ = c(rep(0,numsamp))
		
			
			mm_esta = fminsearch(obj_and, vect, comp = comp)
			m_est = mm_esta$xval
			if (comp == 1){
				a_and_1 =m_est[1:comp]
				b_and_1 =m_est[(comp+1):(2*comp)]
				om_and_1 =m_est[(2*comp+1):(3*comp)]}
				if (comp ==2){
					a_and_2 =m_est[1:comp]
					b_and_2 =m_est[(comp+1):(2*comp)]
					om_and_2 =m_est[(2*comp+1):(3*comp)]}
					if (comp ==3){
						a_and_3 =m_est[1:comp]
						b_and_3 =m_est[(comp+1):(2*comp)]
						om_and_3 =m_est[(2*comp+1):(3*comp)]}
						if (comp ==4){
							a_and_4 =m_est[1:comp]
							b_and_4 =m_est[(comp+1):(2*comp)]
							om_and_4=m_est[(2*comp+1):(3*comp)]}
							if (comp ==5){
								a_and_5 =m_est[1:comp]
								b_and_5 =m_est[(comp+1):(2*comp)]
								om_and_5 =m_est[(2*comp+1):(3*comp)]}
								if (comp ==6){
									a_and_6 =m_est[1:comp]
									b_and_6 =m_est[(comp+1):(2*comp)]
									om_and_6 =m_est[(2*comp+1):(3*comp)]}
###				
        ob_and = c(rep(0,numsamp))
		differ = c(rep(0,numsamp))
		
		for(k in 1:numsamp){
      
			yp[k] = 0
			for(kp in 1:comp){
					if(comp == 1){
						yp[k] = yp[k]+(a_and_1[kp]*cos(om_and_1[kp]*k)+b_and_1[kp]*sin(om_and_1[kp]*k))
					}
						if(comp == 2){
							yp[k] = yp[k]+(a_and_2[kp]*cos(om_and_2[kp]*k)+b_and_2[kp]*sin(om_and_2[kp]*k))
						}
							if(comp == 3){
								yp[k] = yp[k]+(a_and_3[kp]*cos(om_and_3[kp]*k)+b_and_3[kp]*sin(om_and_3[kp]*k))
							}
								if(comp == 4){
									yp[k] = yp[k]+(a_and_4[kp]*cos(om_and_4[kp]*k)+b_and_4[kp]*sin(om_and_4[kp]*k))
								}
									if(comp == 5){
										yp[k] = yp[k]+(a_and_5[kp]*cos(om_and_5[kp]*k)+b_and_5[kp]*sin(om_and_5[kp]*k))
									}
										if(comp == 6){
											yp[k] = yp[k]+(a_and_6[kp]*cos(om_and_6[kp]*k)+b_and_6[kp]*sin(om_and_6[kp]*k))
										}
				}
		differ[k] = (y[k] - yp[k])/sig_est	
			
		if(abs(differ[k])<=aa * pi){
        ob_and[k] = ob_and[k] + aa*(1-cos(differ[k]/aa))
      }
      
		else{
        ob_and[k] = ob_and[k] + 2*aa
		}
	}
	ob_and_sum = sum(ob_and)
###				
#ROBUSTE INFORMATIONSKRITERIEN MITTELS ANDREW´S METHODE		
obj_and_bic[comp] = ob_and_sum +0.5*(log(numsamp))*(3*comp)
obj_and_aic[comp] = ob_and_sum+(3*comp)
obj_and_aicc[comp] = ob_and_sum + (3*comp*numsamp)/(numsamp-3*comp-1)
ca = (3*comp*numsamp)/(numsamp-3*comp-1)
cb = 0.5*(log(numsamp))*(3*comp)
obj_and_wic[comp] = ((ca)/(ca+cb))*obj_and_aicc[comp]+((cb)/(ca+cb))*obj_and_bic[comp]
obj_and_hq[comp] = ob_and_sum + (3*comp)*log(log(numsamp))		

###HAMPEL`S METHODS	

		
		nummax = comp
		yadj = y
		
		 if(comp > 1){
		 rm(m_est)}
		 
			yp = c(rep(0,numsamp))
			differ = c(rep(0,numsamp))
			
			mm_esta = fminsearch(obj_ham, vect, comp = comp)
			m_est = mm_esta$xval
			
			if (comp == 1){
				a_ham_1 =m_est[1:comp]
				b_ham_1 =m_est[(comp+1):(2*comp)]
				om_ham_1 =m_est[(2*comp+1):(3*comp)]}
				if (comp ==2){
					a_ham_2 =m_est[1:comp]
					b_ham_2 =m_est[(comp+1):(2*comp)]
					om_ham_2 =m_est[(2*comp+1):(3*comp)]}
					if (comp ==3){
						a_ham_3 =m_est[1:comp]
						b_ham_3 =m_est[(comp+1):(2*comp)]
						om_ham_3 =m_est[(2*comp+1):(3*comp)]}
						if (comp ==4){
							a_ham_4 =m_est[1:comp]
							b_ham_4 =m_est[(comp+1):(2*comp)]
							om_ham_4=m_est[(2*comp+1):(3*comp)]}
							if (comp ==5){
								a_ham_5 =m_est[1:comp]
								b_ham_5 =m_est[(comp+1):(2*comp)]
								om_ham_5 =m_est[(2*comp+1):(3*comp)]}
								if (comp ==6){
									a_ham_6 =m_est[1:comp]
									b_ham_6 =m_est[(comp+1):(2*comp)]
									om_ham_6 =m_est[(2*comp+1):(3*comp)]}
###
      ob_ham = c(rep(0,numsamp))
	  differ = c(rep(0,numsamp))
	  
	  for(k in 1:numsamp){
      
			yp[k] = 0
			for(kp in 1:comp){
					if(comp == 1){
						yp[k] = yp[k]+(a_ham_1[kp]*cos(om_ham_1[kp]*k)+b_ham_1[kp]*sin(om_ham_1[kp]*k))
					}
						if(comp == 2){
							yp[k] = yp[k]+(a_ham_2[kp]*cos(om_ham_2[kp]*k)+b_ham_2[kp]*sin(om_ham_2[kp]*k))
						}
							if(comp == 3){
								yp[k] = yp[k]+(a_ham_3[kp]*cos(om_ham_3[kp]*k)+b_ham_3[kp]*sin(om_ham_3[kp]*k))
							}
								if(comp == 4){
									yp[k] = yp[k]+(a_ham_4[kp]*cos(om_ham_4[kp]*k)+b_ham_4[kp]*sin(om_ham_4[kp]*k))
								}
									if(comp == 5){
										yp[k] = yp[k]+(a_ham_5[kp]*cos(om_ham_5[kp]*k)+b_ham_5[kp]*sin(om_ham_5[kp]*k))
									}
										if(comp == 6){
											yp[k] = yp[k]+(a_ham_6[kp]*cos(om_ham_6[kp]*k)+b_ham_6[kp]*sin(om_ham_6[kp]*k))
										}
				}
	  	differ[k] = (y[k] - yp[k])/sig_est
	  	
	  	if(abs(differ[k])<=ha){
           ob_ham[k] = ob_ham[k] + (differ[k]^2)/2
      }
      
		if(ha < abs(differ[k]) & abs(differ[k])<=hb){
			ob_ham[k] = ob_ham[k] + abs(differ[k])*ha - (ha^2)/2
		}
      
			if(hb < abs(differ[k]) & abs(differ[k])<= hc){
				ob_ham[k] = ob_ham[k] + ha*(hc*abs(differ[k])) - ((abs(differ[k])^2/2)/(hc - hb)) - (7/6)*(ha^2)
			}
      
				else{
					ob_ham[k] = ob_ham[k] + ha*(hb+hc-ha)
				}
    }
    ob_ham_sum = sum(ob_ham)

#ROBUSTE INFORMATIONSKRITERIEN MITTELS HAMPEL´S METHODE
obj_ham_bic[comp] = ob_ham_sum +0.5*(log(numsamp))*(3*comp)
obj_ham_aic[comp] = ob_ham_sum+(3*comp)
obj_ham_aicc[comp] = ob_ham_sum + (3*comp*numsamp)/(numsamp-3*comp-1)
ca = (3*comp*numsamp)/(numsamp-3*comp-1)
cb = 0.5*(log(numsamp))*(3*comp)
obj_ham_wic[comp] = ((ca)/(ca+cb))*obj_ham_aicc[comp]+((cb)/(ca+cb))*obj_ham_bic[comp]
obj_ham_hq[comp] = ob_ham_sum + (3*comp)*log(log(numsamp))
	
###NICHT ROBUSTE METHODEN

		
		nummax = comp
		yadj = y
		
		 if(comp > 1){
		 rm(m_est)}
		 
			yp = c(rep(0,numsamp))
			differ = c(rep(0,numsamp))
			
			mm_esta = fminsearch(obj_ll, vect, comp = comp)
			m_est = mm_esta$xval
			
			if (comp == 1){
				a_1 =m_est[1:comp]
				b_1 =m_est[(comp+1):(2*comp)]
				om_1 =m_est[(2*comp+1):(3*comp)]}
				if (comp ==2){
					a_2 =m_est[1:comp]
					b_2 =m_est[(comp+1):(2*comp)]
					om_2 =m_est[(2*comp+1):(3*comp)]}
					if (comp ==3){
						a_3 =m_est[1:comp]
						b_3 =m_est[(comp+1):(2*comp)]
						om_3 =m_est[(2*comp+1):(3*comp)]}
						if (comp ==4){
							a_4 =m_est[1:comp]
							b_4 =m_est[(comp+1):(2*comp)]
							om_4=m_est[(2*comp+1):(3*comp)]}
							if (comp ==5){
								a_5 =m_est[1:comp]
								b_5 =m_est[(comp+1):(2*comp)]
								om_5 =m_est[(2*comp+1):(3*comp)]}
								if (comp ==6){
									a_6 =m_est[1:comp]
									b_6 =m_est[(comp+1):(2*comp)]
									om_6 =m_est[(2*comp+1):(3*comp)]}
				
###
        ob_non = c(rep(0,numsamp))
		differ = c(rep(0,numsamp))
		
		for(k in 1:numsamp){
      
			yp[k] = 0
			for(kp in 1:comp){
					if(comp == 1){
						yp[k] = yp[k]+(a_1[kp]*cos(om_1[kp]*k)+b_1[kp]*sin(om_1[kp]*k))
					}
						if(comp == 2){
							yp[k] = yp[k]+(a_2[kp]*cos(om_2[kp]*k)+b_2[kp]*sin(om_2[kp]*k))
						}
							if(comp == 3){
								yp[k] = yp[k]+(a_3[kp]*cos(om_3[kp]*k)+b_3[kp]*sin(om_3[kp]*k))
							}
								if(comp == 4){
									yp[k] = yp[k]+(a_4[kp]*cos(om_4[kp]*k)+b_4[kp]*sin(om_4[kp]*k))
								}
									if(comp == 5){
										yp[k] = yp[k]+(a_5[kp]*cos(om_5[kp]*k)+b_5[kp]*sin(om_5[kp]*k))
									}
										if(comp == 6){
											yp[k] = yp[k]+(a_6[kp]*cos(om_6[kp]*k)+b_6[kp]*sin(om_6[kp]*k))
										}
			} 
		differ[k] = (y[k] - yp[k])
        ob_non[k] = ob_non[k]+differ[k]^2
	   }
	     ob_non_sum =sum(ob_non)
#ROBUSTE INFORMATIONSKRITERIEN MITTELS NICHT ROBUSTEN METHODEN --> SIMULATIONEN FUER KAPITEL 2 
obj_bic[comp] = numsamp*log(ob_non_sum)+(log(numsamp))*(3*comp)
bic_cor[comp]=numsamp*log(ob_non_sum)+(log(numsamp))*(5*comp);
obj_aic[comp]=numsamp*log(ob_non_sum)+2*(3*comp)
obj_aicc[comp]=numsamp*log(ob_non_sum)+2*(3*comp*numsamp)/(numsamp-3*comp-1)
ca=(3*comp*numsamp)/(numsamp-3*comp-1) 
cb=0.5*(log(numsamp))*(3*comp)
obj_wic[comp]=((ca)/(ca+cb))*obj_aicc[comp]+((cb)/(ca+cb))*bic_cor[comp]
obj_hq[comp]=numsamp*log(ob_non_sum)+2*(3*comp)*log(log(numsamp))   
}
#Ende der for-Schleife fuer die Komponenten 'comp'

#####################################################################################################################
#####################################################################################################################
obj_hub_aic_m[sim,] = obj_hub_aic
obj_hub_bic_m[sim,] = obj_hub_bic
obj_hub_aicc_m[sim,] = obj_hub_aicc
obj_hub_hq_m[sim,] = obj_hub_hq
obj_hub_wic_m[sim,] = obj_hub_wic
###
obj_ram_aic_m[sim,] = obj_ram_aic
obj_ram_bic_m[sim,] = obj_ram_bic
obj_ram_aicc_m[sim,] = obj_ram_aicc
obj_ram_hq_m[sim,] = obj_ram_hq
obj_ram_wic_m[sim,] = obj_ram_wic
###
obj_and_aic_m[sim,] = obj_and_aic
obj_and_bic_m[sim,] = obj_and_bic
obj_and_aicc_m[sim,] = obj_and_aicc
obj_and_hq_m[sim,] = obj_and_hq
obj_and_wic_m[sim,] = obj_and_wic
###
obj_ham_aic_m[sim,] = obj_ham_aic
obj_ham_bic_m[sim,] = obj_ham_bic
obj_ham_aicc_m[sim,] = obj_ham_aicc
obj_ham_hq_m[sim,] = obj_ham_hq
obj_ham_wic_m[sim,] = obj_ham_wic
###
obj_aic_m[sim,] = obj_aic
obj_bic_m[sim,] = obj_bic
obj_aicc_m[sim,] = obj_aicc
obj_hq_m[sim,] = obj_hq
obj_wic_m[sim,] = obj_wic
#
#Fuer Jede Iteration der Monte-Carlo Simulation werden die Werte der Informationskriterien in Matritzen gespeichert,
#um spaeter das Minimum zu ermitteln.

########################################################################
a_hub_1_m[sim,] = a_hub_1
b_hub_1_m[sim,]= b_hub_1 
om_hub_1_m[sim,] = om_hub_1

a_hub_2_m[sim,] =a_hub_2
b_hub_2_m[sim,]= b_hub_2
om_hub_2_m[sim,] =om_hub_2

a_hub_3_m[sim,] =a_hub_3
b_hub_3_m[sim,]= b_hub_3
om_hub_3_m[sim,] =om_hub_3

a_hub_4_m[sim,] = a_hub_4
b_hub_4_m[sim,] = b_hub_4
om_hub_4_m[sim,] =om_hub_4

a_hub_5_m[sim,] = a_hub_5
b_hub_5_m[sim,]= b_hub_5
om_hub_5_m[sim,] =om_hub_5

a_hub_6_m[sim,] = a_hub_6
b_hub_6_m[sim,]= b_hub_6
om_hub_6_m[sim,] =om_hub_6

a_hub_ges_m = cbind(a_hub_1_m,a_hub_2_m,a_hub_3_m,a_hub_4_m,a_hub_5_m,a_hub_6_m)
b_hub_ges_m = cbind(b_hub_1_m,b_hub_2_m,b_hub_3_m,b_hub_4_m,b_hub_5_m,b_hub_6_m)
om_hub_ges_m = cbind(om_hub_1_m,om_hub_2_m,om_hub_3_m,om_hub_4_m,om_hub_5_m,om_hub_6_m)

#Speichern aller Schaetzer der Huber Methode in Matrizen fuer alle Iterationen der MC-Simulation.
#Die Matritzen mit ges_m verbinden alle Schaetzer, dass ist wichtig fuer spaetere Funktionen.

###
a_ram_1_m[sim,] = a_ram_1
b_ram_1_m[sim,]= b_ram_1 
om_ram_1_m[sim,] = om_ram_1

a_ram_2_m[sim,] =a_ram_2
b_ram_2_m[sim,]= b_ram_2
om_ram_2_m[sim,] =om_ram_2

a_ram_3_m[sim,] =a_ram_3
b_ram_3_m[sim,]= b_ram_3
om_ram_3_m[sim,] =om_ram_3

a_ram_4_m[sim,] = a_ram_4
b_ram_4_m[sim,] = b_ram_4
om_ram_4_m[sim,] =om_ram_4

a_ram_5_m[sim,] = a_ram_5
b_ram_5_m[sim,]= b_ram_5
om_ram_5_m[sim,] =om_ram_5

a_ram_6_m[sim,] = a_ram_6
b_ram_6_m[sim,]= b_ram_6
om_ram_6_m[sim,] =om_ram_6

a_ram_ges_m = cbind(a_ram_1_m,a_ram_2_m,a_ram_3_m,a_ram_4_m,a_ram_5_m,a_ram_6_m)
b_ram_ges_m = cbind(b_ram_1_m,b_ram_2_m,b_ram_3_m,b_ram_4_m,b_ram_5_m,b_ram_6_m)
om_ram_ges_m = cbind(om_ram_1_m,om_ram_2_m,om_ram_3_m,om_ram_4_m,om_ram_5_m,om_ram_6_m)

#Speichern aller Schaetzer der Ramsey Methode in Matrizen fuer alle Iterationen der MC-Simulation.
###
a_and_1_m[sim,] = a_and_1
b_and_1_m[sim,]= b_and_1 
om_and_1_m[sim,] = om_and_1

a_and_2_m[sim,] =a_and_2
b_and_2_m[sim,]= b_and_2
om_and_2_m[sim,] =om_and_2

a_and_3_m[sim,] =a_and_3
b_and_3_m[sim,]= b_and_3
om_and_3_m[sim,] =om_and_3

a_and_4_m[sim,] = a_and_4
b_and_4_m[sim,] = b_and_4
om_and_4_m[sim,] =om_and_4

a_and_5_m[sim,] = a_and_5
b_and_5_m[sim,]= b_and_5
om_and_5_m[sim,] =om_and_5

a_and_6_m[sim,] = a_and_6
b_and_6_m[sim,]= b_and_6
om_and_6_m[sim,] =om_and_6

a_and_ges_m = cbind(a_and_1_m,a_and_2_m,a_and_3_m,a_and_4_m,a_and_5_m,a_and_6_m)
b_and_ges_m = cbind(b_and_1_m,b_and_2_m,b_and_3_m,b_and_4_m,b_and_5_m,b_and_6_m)
om_and_ges_m = cbind(om_and_1_m,om_and_2_m,om_and_3_m,om_and_4_m,om_and_5_m,om_and_6_m)

#Speichern aller Schaetzer der Andrew Methode in Matrizen fuer alle Iterationen der MC-Simulation.
###
a_ham_1_m[sim,] = a_ham_1
b_ham_1_m[sim,]= b_ham_1 
om_ham_1_m[sim,] = om_ham_1

a_ham_2_m[sim,] =a_ham_2
b_ham_2_m[sim,]= b_ham_2
om_ham_2_m[sim,] =om_ham_2

a_ham_3_m[sim,] =a_ham_3
b_ham_3_m[sim,]= b_ham_3
om_ham_3_m[sim,] =om_ham_3

a_ham_4_m[sim,] = a_ham_4
b_ham_4_m[sim,] = b_ham_4
om_ham_4_m[sim,] =om_ham_4

a_ham_5_m[sim,] = a_ham_5
b_ham_5_m[sim,]= b_ham_5
om_ham_5_m[sim,] =om_ham_5

a_ham_6_m[sim,] = a_ham_6
b_ham_6_m[sim,]= b_ham_6
om_ham_6_m[sim,] =om_ham_6

a_ham_ges_m = cbind(a_ham_1_m,a_ham_2_m,a_ham_3_m,a_ham_4_m,a_ham_5_m,a_ham_6_m)
b_ham_ges_m = cbind(b_ham_1_m,b_ham_2_m,b_ham_3_m,b_ham_4_m,b_ham_5_m,b_ham_6_m)
om_ham_ges_m = cbind(om_ham_1_m,om_ham_2_m,om_ham_3_m,om_ham_4_m,om_ham_5_m,om_ham_6_m)

#Speichern aller Schaetzer der Hampel Methode in Matrizen fuer alle Iterationen der MC-Simulation.
###
a_1_m[sim,] = a_1
b_1_m[sim,]= b_1 
om_1_m[sim,] = om_1

a_2_m[sim,] =a_2
b_2_m[sim,]= b_2
om_2_m[sim,] =om_2

a_3_m[sim,] =a_3
b_3_m[sim,]= b_3
om_3_m[sim,] =om_3

a_4_m[sim,] = a_4
b_4_m[sim,] = b_4
om_4_m[sim,] =om_4

a_5_m[sim,] = a_5
b_5_m[sim,]= b_5
om_5_m[sim,] =om_5

a_6_m[sim,] = a_6
b_6_m[sim,]= b_6
om_6_m[sim,] =om_6

a_ges_m = cbind(a_1_m,a_2_m,a_3_m,a_4_m,a_5_m,a_6_m)
b_ges_m = cbind(b_1_m,b_2_m,b_3_m,b_4_m,b_5_m,b_6_m)
om_ges_m = cbind(om_1_m,om_2_m,om_3_m,om_4_m,om_5_m,om_6_m)

#Speichern aller Schaetzer nicht robuster Methoden in Matrizen fuer alle Iterationen der MC-Simulation.
########################################################################
min_hub_aic[sim] = which.min(obj_hub_aic_m[sim,])
min_hub_bic[sim] = which.min(obj_hub_bic_m[sim,])
min_hub_aicc[sim] = which.min(obj_hub_aicc_m[sim,])
min_hub_hq[sim] = which.min(obj_hub_hq_m[sim,])
min_hub_wic[sim] = which.min(obj_hub_wic_m[sim,])
###
min_ram_aic[sim] = which.min(obj_ram_aic_m[sim,])
min_ram_bic[sim] = which.min(obj_ram_bic_m[sim,])
min_ram_aicc[sim] = which.min(obj_ram_aicc_m[sim,])
min_ram_hq[sim] = which.min(obj_ram_hq_m[sim,])
min_ram_wic[sim] = which.min(obj_ram_wic_m[sim,])
###
min_and_aic[sim] = which.min(obj_and_aic_m[sim,])
min_and_bic[sim] = which.min(obj_and_bic_m[sim,])
min_and_aicc[sim] = which.min(obj_and_aicc_m[sim,])
min_and_hq[sim] = which.min(obj_and_hq_m[sim,])
min_and_wic[sim] = which.min(obj_and_wic_m[sim,])
###
min_ham_aic[sim] = which.min(obj_ham_aic_m[sim,])
min_ham_bic[sim] = which.min(obj_ham_bic_m[sim,])
min_ham_aicc[sim] = which.min(obj_ham_aicc_m[sim,])
min_ham_hq[sim] = which.min(obj_ham_hq_m[sim,])
min_ham_wic[sim] = which.min(obj_ham_wic_m[sim,])
###
min_aic[sim] = which.min(obj_aic_m[sim,])
min_bic[sim] = which.min(obj_bic_m[sim,])
min_aicc[sim] = which.min(obj_aicc_m[sim,])
min_hq[sim] = which.min(obj_hq_m[sim,])
min_wic[sim] = which.min(obj_wic_m[sim,])
#Berechnung des Minimums der Informationskriterien

print(sim)
}
#Ende der Monte-Calo Simulation
########################################################################
