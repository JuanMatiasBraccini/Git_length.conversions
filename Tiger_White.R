# Estimation of Tiger and White shark size after release
handl_OneDrive=function(x)paste('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias',x,sep='/')
setwd(handl_OneDrive("Analyses/Length conversions"))
#Growth Parameters

  #Tiger
#source: Holmes et al. 2015
Linf.c_T=376    #(cm, TL)
K.c_T=0.1
Lo.T= 78 

#source: Stevens & McLaughlin 1991
a.T=0.88
b.T=-15.71


  #White
#source: O'Connor 2011
Linf.f_W=719   #(cm, TL)
K.f_W=0.056
Linf.m_W=798.94
K.m_W=0.047
Lo.W=140    #(cm, TL)

#source: Kohler et al 1996
a.W=0.944
b.W=-5.744


#Growth function
Age.to.Len.fn=function(Lo,Linf,k,age) Lo+(Linf-Lo)*(1-exp(-k*age))
Len.to.Age.fn=function(Lo,Linf,k,L) log(1-((L-Lo)/(Linf-Lo)))/-k     
TL.to.FL=function(a,b,TL) (TL*a)+b
FL.to.TL=function(a,b,FL) (FL-b)/a

#Procedure

Species=c("Tiger","White")

  #put parameters in list
FL_range=list(Tiger=c(100,300),White=c(140,400))
FL_TL_pars=list(Tiger=c(a.T,b.T),White=c(a.W,b.W))
Grw_par_combined=list(Tiger="YES",White="NO")

  #years at liberty
Yrs=0:5
Yr_mnths=seq(Yrs[1],Yrs[length(Yrs)],by=(1/12))
Yr_mnths=Yr_mnths[-1]
Names.Months=1:(Yrs[length(Yrs)]*12)

for(i in 1:length(Species))
{
  #1. Get TL from FL of tagged sharks
  MinFL=FL_range[[i]][1]
  MaxFL=FL_range[[i]][2]
  FL.tag=seq(MinFL,MaxFL,by=5)
  
  a=FL_TL_pars[[i]][1]
  b=FL_TL_pars[[i]][2]
  TL.tag=FL.to.TL(a,b,FL.tag)
  
  #2. Get age of tagged shark
  if(Grw_par_combined[[i]]=="NO")
  {
    Age.tag_F=Len.to.Age.fn(Lo.W,Linf.f_W,K.f_W,TL.tag)
    Age.tag_M=Len.to.Age.fn(Lo.W,Linf.m_W,K.m_W,TL.tag)
  }
  if(Grw_par_combined[[i]]=="YES") Age.tag=Len.to.Age.fn(Lo.T,Linf.c_T,K.c_T,TL.tag)
  
  
  #2. Get age after X time at liberty
  Mat=as.data.frame(t(Yr_mnths))
  names(Mat)=Names.Months
  if(Grw_par_combined[[i]]=="NO")
  {
    Mat_F=cbind(Age.tag_F,Mat)
    Mat_M=cbind(Age.tag_M,Mat)
    
    for(p in 2:ncol(Mat_F)) Mat_F[,p]=Mat_F[,p]+Mat_F[,1]
    for(p in 2:ncol(Mat_M)) Mat_M[,p]=Mat_M[,p]+Mat_M[,1]
    
    Mat_F=Mat_F[,-1]
    Mat_M=Mat_M[,-1]
   
  }
  
  if(Grw_par_combined[[i]]=="YES")
    {
      MaT=cbind(Age.tag,Mat)
      for(p in 2:ncol(MaT)) MaT[,p]=MaT[,p]+MaT[,1]
      MaT=MaT[,-1]
    }

  
  #3 Get size after X time at liberty
  if(Grw_par_combined[[i]]=="NO")
  {
    Mat_F_L=apply(Mat_F,1:2,function(x) Age.to.Len.fn(Lo=Lo.W,Linf=Linf.f_W,k=K.f_W,age=x))
    Mat_M_L=apply(Mat_M,1:2,function(x) Age.to.Len.fn(Lo=Lo.W,Linf=Linf.m_W,k=K.m_W,age=x))
    
    #add length at tagging
    Mat_F_L=cbind(FL.tag,TL.tag,Mat_F_L)
    Mat_M_L=cbind(FL.tag,TL.tag,Mat_M_L)
  }
  
  if(Grw_par_combined[[i]]=="YES")
  {
    Mat_L=apply(MaT,1:2,function(x) Age.to.Len.fn(Lo=Lo.T,Linf=Linf.c_T,k=K.c_T,age=x))
    
    #add length at tagging
    Mat_L=cbind(FL.tag,TL.tag,Mat_L)
    
  }
  
}


Mat_F_L=round(Mat_F_L,1)
Mat_M_L=round(Mat_M_L,1)
Mat_L=round(Mat_L,1)

write.csv(Mat_F_L,"White_female.csv",row.names=F)
write.csv(Mat_M_L,"White_male.csv",row.names=F)
write.csv(Mat_L,"Tiger.csv",row.names=F)






