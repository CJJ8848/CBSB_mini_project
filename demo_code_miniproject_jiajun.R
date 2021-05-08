# fraction binding
#install.packages("rootSolve")
library(rootSolve)
library(ggplot2)
#model construction
#model 1 phenotype=dimer expression level
#model 2 phenotype = dimer percent
model_1<- function(a,deltadelta_G1) { # a1 can be seen as a mutationa that affects deltaG of dimerize
  pr1=1 # protein expression level
  deltag1_wt = -3 # WT folding energy, kcal/mol
  R= 1.98*10^(-3) # gas constant
  Temp= 310.15 # abs temperature for 37 degrees, Kelvin
  deltaG1= deltag1_wt + deltadelta_G1  # Mut folding energy
  k = exp(-deltaG1/(R*Temp))  # dimerize kinetic
  dslnex <- function(x) {
    y <- numeric(2) # x[1]= dimerized, x[2]= unfolded. 
    y[1] <- 2*x[1]+x[2]- pr1*a  # equation for total protein amount
    y[2] <- x[2]*x[2]*k-x[1] # equation for folding-unfolding kinetics
    y
  } 
  xstart <- c(0.9*pr1*a,0.1*pr1*a)  # give a starting position to search for the solution
  components= multiroot( dslnex, start=xstart, positive = TRUE,rtol=1e-6,atol= 1e-10, ctol= 1e-10)$root
  return(components[1]/pr1) 
}
model_2<- function(a,deltadelta_G1) { # a1 can be seen as a mutationa that affects deltaG of dimerize
  pr1=1 # protein expression level
  deltag1_wt = -3 # WT folding energy, kcal/mol
  R= 1.98*10^(-3) # gas constant
  Temp= 310.15 # abs temperature for 37 degrees, Kelvin
  deltaG1= deltag1_wt + deltadelta_G1  # Mut folding energy
  k = exp(-deltaG1/(R*Temp))  # dimerize kinetic
  dslnex <- function(x) {
    y <- numeric(2) # x[1]= dimerized, x[2]= unfolded. 
    y[1] <- 2*x[1]+x[2]- pr1*a  # equation for total protein amount
    y[2] <- x[2]*x[2]*k-x[1] # equation for folding-unfolding kinetics
    y
  } 
  xstart <- c(0.9*pr1*a,0.1*pr1*a)  # give a starting position to search for the solution
  components= multiroot( dslnex, start=xstart, positive = TRUE,rtol=1e-6,atol= 1e-10, ctol= 1e-10)$root
  return(2*components[1]/(pr1*a)) 
}
#q1 for dimerization mutation
a1=1
deltadelta_G1=seq(-4,9,by=0.1)
dimerized_percent={}
for (i in 1:length(deltadelta_G1)){
  dimerized_percent=c(dimerized_percent,model_1(a1,deltadelta_G1[i]))
}
df=data.frame(deltadelta_G1,dimerized_percent)

p1<-ggplot(df,aes(x=deltadelta_G1,y=dimerized_percent))+xlab('ΔΔG (kcal/mol)')+ylab('Dimer expression')+geom_line(colour='red')+geom_vline(xintercept=0,colour='black',linetype="dotted")
p1
#q2 for expression mut
a2=seq(0,10,by=0.01)
a3=seq(0,10,by=0.01)
dimerized_percent1={}
dimerized_percent2={}
deltadelta_G2=0
for (i in 1:length(a2)){
  dimerized_percent1=c(dimerized_percent1,model_1(a2[i],deltadelta_G2))
}
df2=data.frame(a2,dimerized_percent1)
for (i in 1:length(a3)){
  dimerized_percent2=c(dimerized_percent2,model_2(a3[i],deltadelta_G2))
}
df3=data.frame(a3,dimerized_percent2)
p2<-ggplot(df2,aes(x=a2,y=dimerized_percent1))+xlab('expression mut')+ylab('Dimer expression')+geom_line(colour='red')+geom_vline(xintercept=1,colour='black',linetype="dotted")+geom_hline(yintercept = 0.5, colour = "red",linetype="dotted") 
p2
p3<-ggplot(df3,aes(x=a3,y=dimerized_percent2))+xlab('expression mut')+ylab('Dimer percent conversion')+geom_line(colour='red')+geom_vline(xintercept=1,colour='black',linetype="dotted")+geom_hline(yintercept = 0.5, colour = "red",linetype="dotted") 
p3
#q2 a : double mut for simple model
# combine expression mut and deltadeltaG
a_double=seq(0,2,by=0.01)
deltadelta_G1=seq(-4,9,by=0.1)
dexp={}
dfq2_double=as.data.frame(matrix(nrow=0,ncol=3))
new_col = c("expression","deltadelta_G","dexp")
colnames(dfq2_double) <- new_col

for (i in 1:length(a_double)){
  for (j in 1:length(deltadelta_G1)){
    #dimerized_percentfd=c(dimerized_percentfd,modol_extended(a3,deltadelta_GF[i],deltadelta_GD[j]))
    dexp=model_1(a_double[i],deltadelta_G1[j])
    dfq2_double[nrow(dfq2_double)+1,]<-c(a_double[i],deltadelta_G1[j],dexp)
  }}


p_double= ggplot(data= dfq2_double)+geom_tile(aes(x=expression, y=deltadelta_G , fill= log2(dexp)))+
  theme_classic() + labs(x= 'expression mut', y='ΔΔG (kcal/mol)' ) + scale_fill_gradient2(low="green", mid="yellow", high="red", midpoint=mean(log2(dexp)))  +
  xlim(0, 2) + ylim(-4, 9) +
  geom_contour(aes(x= expression, y= deltadelta_G,  z =log2(dexp),colour = ..level..),colour = "black",alpha=0.5 ,binwidth=2 ) +
  geom_point(aes(x= 1, y= 0), shape=4, col='black')
p_double







#q3 b extended mdoel
modol_extended<- function(a1, deltadelta_GF, deltadelta_GD) { # a1 can be seen as a mutationa that affects expression of the protein
  pr1=1 # protein expression level
  deltaGF_wt = -2 # WT folding energy, kcal/mol
  R= 1.98*10^(-3) # gas constant
  Temp= 310.15 # abs temperature for 37 degrees, Kelvin
  deltaGF= deltaGF_wt + deltadelta_GF  # Mut folding energy
  
  k1 = exp(-deltaGF/(R*Temp))  # Folding kinetic
  
  deltaGD_wt = -3 # WT dimerizing energy
  deltaGD = deltaGD_wt + deltadelta_GD # Mut binding energy
  k2= exp(-deltaGD/(R*Temp))
  
  require(rootSolve)
  dslnex <- function(x) {
    y <- numeric(3) # x[1]= dimerized, x[2]= folded,x[3]= unfolded. 
    y[1] <- 2*x[1]+ x[2] + x[3]- pr1*a1  # equation for total protein amount
    y[2] <- x[2]*x[2]*k2-x[1] # equation for dimering-undimering kinetics
    y[3] <- x[3]*k1-x[2] # equation for folding-unfolding kinetics
    y
  } 
  xstart <- c(0.9*pr1*a1,0.08*pr1*a1,0.02*pr1*a1)  # give a starting position to search for the solution
  components= multiroot( dslnex, start=xstart, positive = TRUE,rtol=1e-6,atol= 1e-10, ctol= 1e-10)$root
  return(components[1]/pr1) 
}
modol_extended2<- function(a1, deltadelta_GF, deltadelta_GD) { # a1 can be seen as a mutationa that affects expression of the protein
  pr1=1 # protein expression level
  deltaGF_wt = -2 # WT folding energy, kcal/mol
  R= 1.98*10^(-3) # gas constant
  Temp= 310.15 # abs temperature for 37 degrees, Kelvin
  deltaGF= deltaGF_wt + deltadelta_GF  # Mut folding energy
  
  k1 = exp(-deltaGF/(R*Temp))  # Folding kinetic
  
  deltaGD_wt = -3 # WT dimerizing energy
  deltaGD = deltaGD_wt + deltadelta_GD # Mut binding energy
  k2= exp(-deltaGD/(R*Temp))
  
  require(rootSolve)
  dslnex <- function(x) {
    y <- numeric(3) # x[1]= dimerized, x[2]= folded,x[3]= unfolded. 
    y[1] <- 2*x[1]+ x[2] + x[3]- pr1*a1  # equation for total protein amount
    y[2] <- x[2]*x[2]*k2-x[1] # equation for dimering-undimering kinetics
    y[3] <- x[3]*k1-x[2] # equation for folding-unfolding kinetics
    y
  } 
  xstart <- c(0.9*pr1*a1,0.08*pr1*a1,0.02*pr1*a1)  # give a starting position to search for the solution
  components= multiroot( dslnex, start=xstart, positive = TRUE,rtol=1e-6,atol= 1e-10, ctol= 1e-10)$root
  return(components[1]/(pr1*a1)) 
}

#q3 extend the model # k=k2*k1*k1
# 1 GF change (folding mut)
GF_change<-function(){
  a3=1
  deltadelta_GD=0
  deltadelta_GF=seq(-4,9,by=0.1)
  dimerized_percentf={}
  for (i in 1:length(deltadelta_GF)){
    dimerized_percentf=c(dimerized_percentf,modol_extended(a3,deltadelta_GF[i],deltadelta_GD))
  }
  dfq3f=data.frame(deltadelta_GF,dimerized_percentf)
  
  p4<-ggplot(dfq3f,aes(x=deltadelta_GF,y=dimerized_percentf))+xlab('ΔΔGf (kcal/mol)')+ylab('Dimer expression')+geom_line(colour='red')+geom_vline(xintercept=0,colour='black',linetype="dotted")
  p4
}

# 2 GD change
GD_change<-function(){
  a3=1
  deltadelta_GD=seq(-4,9,by=0.1)
  deltadelta_GF=0
  dimerized_percentd={}
  for (i in 1:length(deltadelta_GD)){
    dimerized_percentd=c(dimerized_percentd,modol_extended(a3,deltadelta_GF,deltadelta_GD[i]))
  }
  dfq3d=data.frame(deltadelta_GD,dimerized_percentd)
  
  p5<-ggplot(dfq3d,aes(x=deltadelta_GD,y=dimerized_percentd))+xlab('ΔΔGd (kcal/mol)')+ylab('Dimer expression')+geom_line(colour='red')+geom_vline(xintercept=0,colour='black',linetype="dotted")
  p5
  
}

# 3 expression mut 
expression1<-function(){
  a3=seq(0,10,by=0.01)
  deltadelta_GD=0
  deltadelta_GF=0
  dimerized_percenta={}
  for (i in 1:length(a3)){
    dimerized_percenta=c(dimerized_percenta,modol_extended(a3[i],deltadelta_GF,deltadelta_GD))
  }
  dfq3a=data.frame(a3,dimerized_percenta)
  
  p6<-ggplot(dfq3a,aes(x=a3,y=dimerized_percenta))+xlab('expression mut')+ylab('Dimer expression')+geom_line(colour='red')+geom_vline(xintercept=1,colour='black',linetype="dotted")
  p6
  
}

# 3 expression mut y<-propotion of dimers 
expression2<-function(){
  a3=seq(0,10,by=0.01)
  deltadelta_GD=0
  deltadelta_GF=0
  dimerized_percenta={}
  for (i in 1:length(a3)){
    dimerized_percenta=c(dimerized_percenta,modol_extended2(a3[i],deltadelta_GF,deltadelta_GD))
  }
  dfq3ab=data.frame(a3,dimerized_percenta)
  
  p7<-ggplot(dfq3ab,aes(x=a3,y=dimerized_percenta))+xlab('expression mut')+ylab('Dimer percent conversion')+geom_line(colour='red')+geom_vline(xintercept=1,colour='black',linetype="dotted")
  p7
  
}

# 4 combine Gd and Gf change
double_gdgf<- function(){
  a3=1
  deltadelta_GD=seq(-4,9,by=0.1)
  deltadelta_GF=seq(-4,9,by=0.1)
  dexp={}
  dfq3fd=as.data.frame(matrix(nrow=0,ncol=3))
  new_col = c("deltadelta_GF","deltadelta_GD","dexp")
  colnames(dfq3fd) <- new_col
  
  for (i in 1:length(deltadelta_GF)){
    for (j in 1:length(deltadelta_GD)){
      #dimerized_percentfd=c(dimerized_percentfd,modol_extended(a3,deltadelta_GF[i],deltadelta_GD[j]))
      dexp=modol_extended(a3,deltadelta_GF[i],deltadelta_GD[j])
      dfq3fd[nrow(dfq3fd)+1,]<-c(deltadelta_GF[i],deltadelta_GD[j],dexp)
    }}
  p8=ggplot(data= dfq3fd)+geom_tile(aes(x=deltadelta_GF, y=deltadelta_GD , fill= log2(dexp)))+ 
    theme_classic() + labs(x= 'ΔΔGf (kcal/mol)', y='ΔΔGd (kcal/mol)')+ scale_fill_gradient2(low="green", mid="yellow", high="red", midpoint=mean(log2(dexp)))  + 
    xlim(-4, 9) + ylim(-4, 9) +  
    geom_contour(aes(x= deltadelta_GF, y= deltadelta_GD,  z =log2(dexp),colour = ..level..),colour = "black",alpha=0.5 ,binwidth=2 ) + 
    geom_point(aes(x= 0, y= 0), shape=4, col='black')
  p8
}
# 5 combine Gd and expression change
double_gd_exp<-function(){
  a_double=seq(0,2,by=0.01)
  deltadelta_GD=seq(-4,9,by=0.1)
  deltadelta_GF=0
  dexp={}
  dfq3d_exp=as.data.frame(matrix(nrow=0,ncol=3))
  new_col = c("expression1","deltadelta_GD","dexp")
  colnames(dfq3d_exp) <- new_col
  
  for (i in 1:length(a_double)){
    for (j in 1:length(deltadelta_GD)){
      #dimerized_percentd_exp=c(dimerized_percentd_exp,modol_extended(a3,deltadelta_GF[i],deltadelta_GD[j]))
      dexp=modol_extended(a_double[i],deltadelta_GF,deltadelta_GD[j])
      dfq3d_exp[nrow(dfq3d_exp)+1,]<-c(a_double[i],deltadelta_GD[j],dexp)
    }}
  p9=ggplot(data= dfq3d_exp)+geom_tile(aes(x=expression1, y=deltadelta_GD , fill= log2(dexp)))+ 
    theme_classic() + labs(x= 'expression mut', y='ΔΔGd (kcal/mol)' ) + scale_fill_gradient2(low="green", mid="yellow", high="red", midpoint=mean(log2(dexp)))  + 
    xlim(0,2) + ylim(-4, 9) +  
    geom_contour(aes(x= expression1, y= deltadelta_GD,  z =log2(dexp),colour = ..level..),colour = "black",alpha=0.5 ,binwidth=2 ) + 
    geom_point(aes(x= 1, y= 0), shape=4, col='black')
  p9
}
# 6 combine Gf and expression change
double_gf_exp<-function(){
  a_double=seq(0,2,by=0.01)
  deltadelta_GF=seq(-4,9,by=0.1)
  deltadelta_GD=0
  dexp={}
  dfq3f_exp=as.data.frame(matrix(nrow=0,ncol=3))
  new_col = c("expression","deltadelta_GF","dexp")
  colnames(dfq3f_exp) <- new_col
  
  for (i in 1:length(a_double)){
    for (j in 1:length(deltadelta_GF)){
      #dimerized_percentd_exp=c(dimerized_percentd_exp,modol_extended(a3,deltadelta_GF[i],deltadelta_GD[j]))
      dexp=modol_extended(a_double[i],deltadelta_GF[j],deltadelta_GD)
      dfq3f_exp[nrow(dfq3f_exp)+1,]<-c(a_double[i],deltadelta_GF[j],dexp)
    }}
  p10=ggplot(data= dfq3f_exp)+geom_tile(aes(x=expression, y=deltadelta_GF , fill= log2(dexp)))+ 
    theme_classic() + labs(x= 'expression mut', y='ΔΔGf (kcal/mol)' ) + scale_fill_gradient2(low="green", mid="yellow", high="red", midpoint=mean(log2(dexp)))  + 
    xlim(0,2) + ylim(-4, 9) +  
    geom_contour(aes(x= expression, y= deltadelta_GF,  z =log2(dexp),colour = ..level..),colour = "black",alpha=0.5 ,binwidth=2 ) + 
    geom_point(aes(x= 1, y= 0), shape=4, col='black')
  p10 
}
 GF_change()
 GD_change()
 expression1()
 expression2()
 double_gdgf()
 double_gd_exp()
 double_gf_exp()