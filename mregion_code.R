library(readr)
#define functions we will need 
loghconst=function(delta, Omega){
  delta/2*as.numeric(determinant(Omega, logarithm = TRUE)$modulus)-0.5*nrow(Omega)*delta*log(2)- logGamma_p(0.5*delta, p=nrow(Omega))
}

logGamma_p=function(a, p){
  p*(p-1)/4*log(pi)+sum(lgamma(a+(1-1:p)/2))
}
#read in new 3D dataset
new_dat <- read_csv("csv_files_3D/118.csv")
Y_new<-new_dat[,c("ADC","KTRANS","KEP","AUGC")]
Y_new$ADC<-Y_new$ADC/100
Y_new[,2:4]<-log(Y_new[,2:4])
Y_new<-as.matrix(Y_new)
zones<-new_dat$zone
m<-4

#load trained quantities to use in calculation
trained_quantities <- readRDS("~/Documents/dissertation/paper3/trained_quantities.rds")
n11<-trained_quantities$ns[1]
n01<-trained_quantities$ns[2]
n12<-trained_quantities$ns[3]
n02<-trained_quantities$ns[4]

y_tilde_11<-trained_quantities$y_tildes[[1]]
y_tilde_01<-trained_quantities$y_tildes[[2]]
y_tilde_12<-trained_quantities$y_tildes[[3]]
y_tilde_02<-trained_quantities$y_tildes[[4]]

term1_11<-trained_quantities$term1s[[1]]
term1_01<-trained_quantities$term1s[[2]]
term1_12<-trained_quantities$term1s[[3]]
term1_02<-trained_quantities$term1s[[4]]

#other terms based on trained quantities
S_tilde_11<-term1_11-n11*(y_tilde_11)%*%t(y_tilde_11)
S_tilde_01<-term1_01-n01*(y_tilde_01)%*%t(y_tilde_01)
S_tilde_12<-term1_12-n12*(y_tilde_12)%*%t(y_tilde_12)
S_tilde_02<-term1_02-n02*(y_tilde_02)%*%t(y_tilde_02)

#region specific cancer probabilities based off training data 
phat_1=n11/(n11+n01)
phat_2=n12/(n12+n02)



#calculate log posterior predictive densities for each voxel
log_f0<-NA
log_f1<-NA
for(i in 1:length(zones)){
  zone<-zones[i]

  if(zone=="CG"){
   
    n0<-n01
    n1<-n11
    y_tilde_0=(y_tilde_01*n0+Y_new[i,])/(n0+1)
    y_tilde_1=(y_tilde_11*n1+Y_new[i,])/(n1+1)
    term1_0<-term1_01+Y_new[i,]%*%t(Y_new[i,])
    term1_1<-term1_11+Y_new[i,]%*%t(Y_new[i,])
    omega_0=term1_0-(n0+1)*y_tilde_0%*%t(y_tilde_0)
    omega_1=term1_1-(n1+1)*y_tilde_1%*%t(y_tilde_1)
    S_tilde_0<-S_tilde_01
    S_tilde_1<-S_tilde_11
  }else{
    
    n0<-n02
    n1<-n12
    y_tilde_0=(y_tilde_02*n0+Y_new[i,])/(n0+1)
    y_tilde_1=(y_tilde_12*n1+Y_new[i,])/(n1+1)
    term1_0<-term1_02+Y_new[i,]%*%t(Y_new[i,])
    term1_1<-term1_12+Y_new[i,]%*%t(Y_new[i,])
    omega_0=term1_0-(n0+1)*y_tilde_0%*%t(y_tilde_0)
    omega_1=term1_1-(n1+1)*y_tilde_1%*%t(y_tilde_1)
    S_tilde_0<-S_tilde_02
    S_tilde_1<-S_tilde_12
  }
  
  log_f0[i]<--m/2*log(2*pi)-m/2*(log(n0+1)-log(n0))+loghconst(delta=n0+m-1,Omega = S_tilde_0+omega_0)-loghconst(delta=n0+m,Omega = omega_0+omega_0)
  log_f1[i]<--m/2*log(2*pi)-m/2*(log(n1+1)-log(n1))+loghconst(delta=n1+m-1,Omega = S_tilde_1+omega_1)-loghconst(delta=n1+m,Omega = omega_1+omega_1)

}



phat<-rep(phat_1,length(zones))
phat[zones=="PZ"]<-phat_2

#calculate posterior predictive probability 
p=phat*exp(log_f1)/(phat*exp(log_f1)+(1-phat)*exp(log_f0))

new_dat$mregion<-p
