###########################################################
## TO START THE MONTE-CARLO SIMULATION WE HAVE TO DEFINE ##
## ALL VARIABLES IN USE.                                 ##
###########################################################



#########################################################
para = 6
#Number of max. Parameter p
numsamp = 50
#length of time series
gridnum = 1000
#number of grid points for periodogram
a_o = c(5,3)
#original value of amplitude vector a 
b_o = c(4,2)
#original value of amplitude vector b
om_o = c(0.4,0.6)
##original value of frequency vector omega
isim =  100
#number of iteration for Monte-Carlo Simulation
numpar= length(a_o)
#number of parameters
##########################################################
#Variables for Distributions
stdev=sqrt(0.01) 
stdev1=sqrt(0.1)
stdev2=sqrt(0.01)
stdev3=sqrt(0.6)
stdev4=sqrt(0.01)
stdev5=sqrt(2)
stdev6=sqrt(0.1)
stdev7=sqrt(3)
stdev8=sqrt(0.1)
dft=2
##########################################################
#per_max function variables
om_grid = c(rep(0,(gridnum+1)))
per = c(rep(0,(gridnum+1)))
a_om = matrix(0,numsamp,2)
yvec = c(rep(0,numsamp))
inits = c(0,0,0)
#calc variables
yp = c(rep(0,numsamp))
differ = c(rep(0,numsamp))
#strarting values for periodogram
seq_a = c(rep(0,para))
seq_b = c(rep(0,para))
seq_om = c(rep(0,para))
###
a = c(rep(0,para))
b = c(rep(0,para))
om = c(rep(0,para))
#robust var
a_m = c(rep(0,para))
b_m = c(rep(0,para))
om_m = c(rep(0,para))
#Huber
diff_full_model = c(rep(0,numsamp))
diff_med_dev = c(rep(0,numsamp))
med_dev= c(rep(0,numsamp))
###
ht = 1.345#Huber
ra=0.3 #Ramsy
aa=1.339 #Andrew
ha=1.7; hb=3.4; hc=8.5; #Hampel
########################################################################
########################################################################

####Define all variables#####
obj_hub_aic  = c(rep(0,para))
obj_hub_bic = c(rep(0,para))
obj_hub_aicc = c(rep(0,para))
obj_hub_hq = c(rep(0,para))
obj_hub_wic = c(rep(0,para))
###
obj_hub_aic_m = matrix(0,isim,para)
obj_hub_bic_m = matrix(0,isim,para)
obj_hub_aicc_m = matrix(0,isim,para)
obj_hub_hq_m = matrix(0,isim,para)
obj_hub_wic_m = matrix(0,isim,para)
###
obj_ram_aic  = c(rep(0,para))
obj_ram_bic = c(rep(0,para))
obj_ram_aicc = c(rep(0,para))
obj_ram_hq = c(rep(0,para))
obj_ram_wic = c(rep(0,para))
###
obj_ram_aic_m = matrix(0,isim,para)
obj_ram_bic_m = matrix(0,isim,para)
obj_ram_aicc_m = matrix(0,isim,para)
obj_ram_hq_m = matrix(0,isim,para)
obj_ram_wic_m = matrix(0,isim,para)
###
obj_and_aic  = c(rep(0,para))
obj_and_bic = c(rep(0,para))
obj_and_aicc = c(rep(0,para))
obj_and_hq = c(rep(0,para))
obj_and_wic = c(rep(0,para))
###
obj_and_aic_m = matrix(0,isim,para)
obj_and_bic_m = matrix(0,isim,para)
obj_and_aicc_m = matrix(0,isim,para)
obj_and_hq_m = matrix(0,isim,para)
obj_and_wic_m = matrix(0,isim,para)
###
obj_ham_aic  = c(rep(0,para))
obj_ham_bic = c(rep(0,para))
obj_ham_aicc = c(rep(0,para))
obj_ham_hq = c(rep(0,para))
obj_ham_wic = c(rep(0,para))
###
obj_ham_aic_m = matrix(0,isim,para)
obj_ham_bic_m = matrix(0,isim,para)
obj_ham_aicc_m = matrix(0,isim,para)
obj_ham_hq_m = matrix(0,isim,para)
obj_ham_wic_m = matrix(0,isim,para)
###
obj_aic  = c(rep(0,para))
obj_bic = c(rep(0,para))
obj_aicc = c(rep(0,para))
obj_hq = c(rep(0,para))
obj_wic = c(rep(0,para))
bic_cor = c(rep(0,para))
###
obj_aic_m = matrix(0,isim,para)
obj_bic_m = matrix(0,isim,para)
obj_aicc_m = matrix(0,isim,para)
obj_hq_m = matrix(0,isim,para)
obj_wic_m = matrix(0,isim,para)
################################
a_hub_1_m = matrix(0,isim,1)
b_hub_1_m= matrix(0,isim,1)
om_hub_1_m =matrix(0,isim,1)
###
a_hub_2_m = matrix(0,isim,2)
b_hub_2_m= matrix(0,isim,2)
om_hub_2_m =matrix(0,isim,2)
###
a_hub_3_m = matrix(0,isim,3)
b_hub_3_m= matrix(0,isim,3)
om_hub_3_m =matrix(0,isim,3)
###
a_hub_4_m = matrix(0,isim,4)
b_hub_4_m= matrix(0,isim,4)
om_hub_4_m =matrix(0,isim,4)
###
a_hub_5_m = matrix(0,isim,5)
b_hub_5_m= matrix(0,isim,5)
om_hub_5_m =matrix(0,isim,5)
###
a_hub_6_m = matrix(0,isim,6)
b_hub_6_m= matrix(0,isim,6)
om_hub_6_m =matrix(0,isim,6)
###########
a_ram_1_m = matrix(0,isim,1)
b_ram_1_m= matrix(0,isim,1)
om_ram_1_m =matrix(0,isim,1)
###
a_ram_2_m = matrix(0,isim,2)
b_ram_2_m= matrix(0,isim,2)
om_ram_2_m =matrix(0,isim,2)
###
a_ram_3_m = matrix(0,isim,3)
b_ram_3_m= matrix(0,isim,3)
om_ram_3_m =matrix(0,isim,3)
###
a_ram_4_m = matrix(0,isim,4)
b_ram_4_m= matrix(0,isim,4)
om_ram_4_m =matrix(0,isim,4)
###
a_ram_5_m = matrix(0,isim,5)
b_ram_5_m= matrix(0,isim,5)
om_ram_5_m =matrix(0,isim,5)
###
a_ram_6_m = matrix(0,isim,6)
b_ram_6_m= matrix(0,isim,6)
om_ram_6_m =matrix(0,isim,6)
###########
a_and_1_m = matrix(0,isim,1)
b_and_1_m= matrix(0,isim,1)
om_and_1_m =matrix(0,isim,1)
###
a_and_2_m = matrix(0,isim,2)
b_and_2_m= matrix(0,isim,2)
om_and_2_m =matrix(0,isim,2)
###
a_and_3_m = matrix(0,isim,3)
b_and_3_m= matrix(0,isim,3)
om_and_3_m =matrix(0,isim,3)
###
a_and_4_m = matrix(0,isim,4)
b_and_4_m= matrix(0,isim,4)
om_and_4_m =matrix(0,isim,4)
###
a_and_5_m = matrix(0,isim,5)
b_and_5_m= matrix(0,isim,5)
om_and_5_m =matrix(0,isim,5)
###
a_and_6_m = matrix(0,isim,6)
b_and_6_m= matrix(0,isim,6)
om_and_6_m =matrix(0,isim,6)
###########
a_ham_1_m = matrix(0,isim,1)
b_ham_1_m= matrix(0,isim,1)
om_ham_1_m =matrix(0,isim,1)
###
a_ham_2_m = matrix(0,isim,2)
b_ham_2_m= matrix(0,isim,2)
om_ham_2_m =matrix(0,isim,2)
###
a_ham_3_m = matrix(0,isim,3)
b_ham_3_m= matrix(0,isim,3)
om_ham_3_m =matrix(0,isim,3)
###
a_ham_4_m = matrix(0,isim,4)
b_ham_4_m= matrix(0,isim,4)
om_ham_4_m =matrix(0,isim,4)
###
a_ham_5_m = matrix(0,isim,5)
b_ham_5_m= matrix(0,isim,5)
om_ham_5_m =matrix(0,isim,5)
###
a_ham_6_m = matrix(0,isim,6)
b_ham_6_m= matrix(0,isim,6)
om_ham_6_m =matrix(0,isim,6)
##########
a_1_m = matrix(0,isim,1)
b_1_m= matrix(0,isim,1)
om_1_m =matrix(0,isim,1)
###
a_2_m = matrix(0,isim,2)
b_2_m= matrix(0,isim,2)
om_2_m =matrix(0,isim,2)
###
a_3_m = matrix(0,isim,3)
b_3_m= matrix(0,isim,3)
om_3_m =matrix(0,isim,3)
###
a_4_m = matrix(0,isim,4)
b_4_m= matrix(0,isim,4)
om_4_m =matrix(0,isim,4)
###
a_5_m = matrix(0,isim,5)
b_5_m= matrix(0,isim,5)
om_5_m =matrix(0,isim,5)
###
a_6_m = matrix(0,isim,6)
b_6_m= matrix(0,isim,6)
om_6_m =matrix(0,isim,6)
############################
a_hub_1 = c(rep(0,1))
b_hub_1= c(rep(0,1))
om_hub_1 =c(rep(0,1))
###
a_hub_2 = c(rep(0,2))
b_hub_2 = c(rep(0,2))
om_hub_2 =c(rep(0,2))
###
a_hub_3 = c(rep(0,3))
b_hub_3 = c(rep(0,3))
om_hub_3 =c(rep(0,3))
###
a_hub_4 =c(rep(0,4))
b_hub_4 = c(rep(0,4))
om_hub_4 =c(rep(0,4))
###
a_hub_5 = c(rep(0,5))
b_hub_5 = c(rep(0,5))
om_hub_5 =c(rep(0,5))
###
a_hub_6 = c(rep(0,6))
b_hub_6= c(rep(0,6))
om_hub_6 =c(rep(0,6))
###########
a_ram_1 = c(rep(0,1))
b_ram_1= c(rep(0,1))
om_ram_1 =c(rep(0,1))
###
a_ram_2 = c(rep(0,2))
b_ram_2 = c(rep(0,2))
om_ram_2 =c(rep(0,2))
###
a_ram_3 = c(rep(0,3))
b_ram_3 = c(rep(0,3))
om_ram_3 =c(rep(0,3))
###
a_ram_4 =c(rep(0,4))
b_ram_4 = c(rep(0,4))
om_ram_4 =c(rep(0,4))
###
a_ram_5 = c(rep(0,5))
b_ram_5 = c(rep(0,5))
om_ram_5 =c(rep(0,5))
###
a_ram_6 = c(rep(0,6))
b_ram_6= c(rep(0,6))
om_ram_6 =c(rep(0,6))
#########
a_and_1 = c(rep(0,1))
b_and_1= c(rep(0,1))
om_and_1 =c(rep(0,1))
###
a_and_2 = c(rep(0,2))
b_and_2 = c(rep(0,2))
om_and_2 =c(rep(0,2))
###
a_and_3 = c(rep(0,3))
b_and_3 = c(rep(0,3))
om_and_3 =c(rep(0,3))
###
a_and_4 =c(rep(0,4))
b_and_4 = c(rep(0,4))
om_and_4 =c(rep(0,4))
###
a_and_5 = c(rep(0,5))
b_and_5 = c(rep(0,5))
om_and_5 =c(rep(0,5))
###
a_and_6 = c(rep(0,6))
b_and_6= c(rep(0,6))
om_and_6 =c(rep(0,6))
##########
a_ham_1 = c(rep(0,1))
b_ham_1= c(rep(0,1))
om_ham_1 =c(rep(0,1))
###
a_ham_2 = c(rep(0,2))
b_ham_2 = c(rep(0,2))
om_ham_2 =c(rep(0,2))
###
a_ham_3 = c(rep(0,3))
b_ham_3 = c(rep(0,3))
om_ham_3 =c(rep(0,3))
###
a_ham_4 =c(rep(0,4))
b_ham_4 = c(rep(0,4))
om_ham_4 =c(rep(0,4))
###
a_ham_5 = c(rep(0,5))
b_ham_5 = c(rep(0,5))
om_ham_5 =c(rep(0,5))
###
a_ham_6 = c(rep(0,6))
b_ham_6= c(rep(0,6))
om_ham_6 =c(rep(0,6))
#########
a_1 = c(rep(0,1))
b_1= c(rep(0,1))
om_1 =c(rep(0,1))
###
a_2 = c(rep(0,2))
b_2 = c(rep(0,2))
om_2 =c(rep(0,2))
###
a_3 = c(rep(0,3))
b_3 = c(rep(0,3))
om_3 =c(rep(0,3))
###
a_4 =c(rep(0,4))
b_4 = c(rep(0,4))
om_4 =c(rep(0,4))
###
a_5 = c(rep(0,5))
b_5 = c(rep(0,5))
om_5 =c(rep(0,5))
###
a_6 = c(rep(0,6))
b_6= c(rep(0,6))
om_6 =c(rep(0,6))
########################################################################
min_hub_aic = c(rep(0,isim))
min_hub_bic = c(rep(0,isim))
min_hub_aicc = c(rep(0,isim))
min_hub_hq = c(rep(0,isim))
min_hub_wic = c(rep(0,isim))
###
min_ram_aic = c(rep(0,isim))
min_ram_bic = c(rep(0,isim))
min_ram_aicc = c(rep(0,isim))
min_ram_hq = c(rep(0,isim))
min_ram_wic = c(rep(0,isim))
###
min_and_aic = c(rep(0,isim))
min_and_bic = c(rep(0,isim))
min_and_aicc = c(rep(0,isim))
min_and_hq = c(rep(0,isim))
min_and_wic = c(rep(0,isim))
###
min_ham_aic = c(rep(0,isim))
min_ham_bic = c(rep(0,isim))
min_ham_aicc = c(rep(0,isim))
min_ham_hq = c(rep(0,isim))
min_ham_wic = c(rep(0,isim))
###
min_aic = c(rep(0,isim))
min_bic = c(rep(0,isim))
min_aicc = c(rep(0,isim))
min_hq = c(rep(0,isim))
min_wic = c(rep(0,isim))
#########################################################################
data_gaussian_normal_m = matrix(0,isim,numsamp)
data_gaussian_normal_mix_001_01_m = matrix(0,isim,numsamp)
data_gaussian_normal_mix_01_001_m = matrix(0,isim,numsamp)
data_gaussian_normal_mix_06_001_m = matrix(0,isim,numsamp)
data_gaussian_normal_mix_2_01_m = matrix(0,isim,numsamp)
data_gaussian_normal_mix_3_01_m = matrix(0,isim,numsamp)
data_gaussian_student_t_m = matrix(0,isim,numsamp)