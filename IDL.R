
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Source_Shark_bio.R"))
Comm.disc.sp=read.csv(handl_OneDrive("Analyses/Ecosystem indices and multivariate/Shark-bycatch/SPECIES+PCS+FATE.csv"),stringsAsFactors = F)


library(ggplot2)
library(tidyverse)
library(MASS)
library(ggpubr)
library(ggrepel)
library(Hmisc)
library(ggpmisc)


#---------Build length-length relationships------------     
#use robust regression to deal with outliers   MISSING: plot recorder in different color and fit model with recorder as factor ACA
Min.obs=5

setwd(handl_OneDrive('Analyses/Length conversions/Interdorsal length conversions'))
IDL.data=DATA%>%
  rename(IDL=TrunkL)%>%
  filter(CAAB_code<50000)%>%
  filter((!is.na(TL) & !is.na(FL)) | (!is.na(TL) & !is.na(IDL)) | (!is.na(FL) & !is.na(IDL)))%>%
  dplyr::select(SHEET_NO,LINE_NO,RECORDER,Mid.Lat,Mid.Long,Lat.round,Long.round,zone,date,Day,Month,year,
                SPECIES,COMMON_NAME,SCIENTIFIC_NAME,Taxa,RetainedFlag,SEX,Number,TL,FL,PL,IDL)
  
names(IDL.data)=tolower(names(IDL.data))

idl.species=table(IDL.data$species,1000*round(IDL.data$idl/1000))
idl.species=rownames(idl.species)[which(idl.species>Min.obs)]

r2ww <- function(x)
{
  SSe <- sum(x$w*(x$resid)^2)
  observed <- x$resid+x$fitted
  SSt <- sum(x$w*(observed-weighted.mean(observed,x$w))^2)
  value <- 1-SSe/SSt;
  return(value);
}

ggplotRegression <- function(dat, xvar, yvar,tolerance=2,SP)
{
  dat=dat%>%
    filter(!is.na((!!sym(xvar))) & !is.na((!!sym(yvar))))
  if(xvar=='tl' | yvar=='tl') dat=dat%>%filter(tl<=400)
  if(xvar=='fl' | yvar=='fl') dat=dat%>%filter(fl<=350)
  if(xvar=='idl' | yvar=='idl') dat=dat%>%filter(idl>10)
  
  p=OUT=NULL
  if(nrow(dat)>=Min.obs)
  {
    if(SP=="GM")   #systematic mis measurment by these recorders
    {
      dat=dat%>%
        filter(!recorder%in%c("AT","NB"))%>%
        filter(!sheet_no%in%c('PA0039'))%>%
        mutate(drop=case_when(idl<45 & tl>160~1,
                              TRUE~0))%>%
        filter(!drop==1)
    }
    if(SP=="PJ")   
    {
      dat=dat%>%
        filter(!recorder%in%c("AT","AS","NB"))%>%
        filter(!sheet_no%in%c('PA0122'))%>%
        mutate(drop=case_when(idl<30 & fl>80~1,
                              idl>35 & tl>150~1,
                              TRUE~0))%>%
        filter(!drop==1)
    }
    if(SP=="TK")   
    {
      dat=dat%>%
        filter(!recorder%in%c("AT"))%>%
        filter(!sheet_no%in%c('PA0017'))%>%
        mutate(drop=case_when(idl>47 & fl<85~1,
                              idl<44 & fl>95~1,
                              idl>42 & fl<75~1,
                              TRUE~0))%>%
        filter(!drop==1)
      
    }
    if(SP=="WH")   
    {
      dat=dat%>%
        mutate(drop=case_when(idl<35 & fl>100~1,
                              idl>40 & fl<60~1,
                              TRUE~0))%>%
        filter(!drop==1)%>%
        filter(idl>30)
    }
    if(SP=="SO")   
    {
      dat=dat%>%
        mutate(drop=case_when(idl>70 & fl<90~1,
                              TRUE~0))%>%
        filter(!drop==1)
    }
    if(SP=="MI")   
    {
      dat=dat%>%
        mutate(drop=case_when(idl>36 & fl<74~1,
                              TRUE~0))%>%
        filter(!drop==1)
    }
    fml <- as.formula(paste(yvar, "~", xvar))
    fit <- rlm(fml, dat)   #robust regression
    d1=dat%>%
      mutate(pred=predict(fit),
             delta=abs(pred-(!!sym(yvar))))
    d1=subset(d1,delta>tolerance*sd(d1$delta))
    p=dat%>%
      ggplot(aes_string(xvar,yvar,label = 'sheet_no'))+geom_point(na.rm= TRUE,colour='black')+
    #  geom_text_repel(data=  d1,size=3,na.rm= TRUE,segment.alpha=.4)+
      stat_smooth(method = "lm", col = "red")+
      labs(title = paste("Adj.R2= ",round(r2ww(fit),2),
                         "Inter.=",signif(fit$coef[[1]],5),
                         " Slope=",signif(fit$coef[[2]], 5)))+
      theme(plot.title = element_text(size=10))
    
    #get coefficients and SE
    out <- summary(fit)
    OUT=data.frame(out$coefficients[,1:2])
    OUT$r2=round(r2ww(fit),3)
    OUT$Species=unique(dat$common_name)
    
  }
  return(list(p=p,OUT=OUT))
  
}
fn.idl=function(SP)
{
  d=IDL.data%>%filter(species==SP)
  TL.FL=ggplotRegression(dat=d, xvar="tl", yvar="fl",SP=SP)
  FL.TL=ggplotRegression(dat=d, xvar="fl", yvar="tl",SP=SP)
  idl.TL=ggplotRegression(dat=d, xvar="idl", yvar="tl",SP=SP)
  idl.FL=ggplotRegression(dat=d, xvar="idl", yvar="fl",SP=SP)
  La.lista=list(TL.FL$p, idl.TL$p, FL.TL$p, idl.FL$p)
  La.lista=La.lista%>% discard(is.null)
  ggarrange(plotlist=La.lista)
  W=8
  H=8
  if(length(La.lista)==2)
  {
    W=8
    H=4
  }
  ggsave(paste(unique(d$common_name),".tiff",sep=''), width = W,height = H, dpi = 300,
         compression = "lzw")
  
  return(list(TL.FL=TL.FL$OUT,FL.TL=FL.TL$OUT,idl.TL=idl.TL$OUT,idl.FL=idl.FL$OUT))
}
Store=vector('list',length(idl.species))

for(i in 1:length(idl.species))  Store[[i]]=fn.idl(SP=idl.species[i])

#export coefficients
for(i in 1:length(idl.species))
{
  write.csv(do.call(rbind,Store[[i]]),paste("Coefficients_",idl.species[i],".csv",sep=""),row.names=T)
}


# Size comp for Maddie --------------------------------------------------------------
bwidth <- 10   # Set binwidth

Maddie=DATA%>%
  filter(Method=="GN" & Mid.Lat<(-26) & CAAB_code<50000 &
           !BOAT%in%c("NAT","HAM","HOU","RV BREAKSEA","RV GANNET","RV SNIPE 2","FLIN") )%>%
  mutate(Whaler=case_when(grepl('Carcharhinus',SCIENTIFIC_NAME)~"Whaler",
                          !grepl('Carcharhinus',SCIENTIFIC_NAME)~"Non-whaler"))%>%
  left_join(Comm.disc.sp%>%dplyr::select(FATE,SPECIES),by="SPECIES")%>%
  filter(FATE=='C')

Most.common=rev(sort(table(Maddie$COMMON_NAME)))
Most.common=c('Sawsharks','Thresher sharks','Shortfin mako','Spotted wobbegong',
              'Bronze whaler','Tiger shark','Wobbegongs',
              'Banded wobbegong','Common sawshark','Pencil shark','Eastern school shark','Spinner shark',
              'Smooth hammerhead','Whiskery shark','Gummy shark','Sandbar shark','Dusky shark')

Maddie=Maddie%>%
  filter(COMMON_NAME%in%Most.common & FL<1000)%>%
  mutate(Whaler.sp=paste(Whaler,COMMON_NAME,sep=' - '))

#set IDL equations
Interdorsal=data.frame(COMMON_NAME=c('Dusky shark','Gummy shark','Smooth hammerhead',
                                     'Pencil shark','Sandbar shark','Spinner shark',
                                     'Whiskery shark','Tiger shark'),
                       a=c(2.3085,1.7064,1.7452,
                           1.1377,2.0587,2.4193,
                           1.6214,1.6723),
                       b=c(-1.1192,16.985,10.657,
                           41.211,-1.7763,-5.4286,
                           25.54,25.732))%>%
            mutate(IDL=70,
                   FL=IDL*a+b,
                   xvar=bwidth*round(FL/bwidth))
  
Maddie%>%
  left_join(Interdorsal%>%dplyr::select(COMMON_NAME,xvar),by='COMMON_NAME')%>%
  ggplot( aes(x=FL, color=Whaler, fill=Whaler)) +
  geom_histogram(alpha=0.6, binwidth = bwidth)+
  facet_wrap(~Whaler.sp,scales="free_y")+
  xlab("Fork length (cm)")+ylab("Count")+
  theme(legend.position="none",
        strip.text.x = element_text(size = 11),
        axis.text=element_text(size=12),
        axis.title=element_text(size=16))+
  geom_segment(aes(x=xvar, xend=xvar,
                   y=1 + 1.5, yend=1),color='black',
               arrow=arrow(length=unit(4, "mm")))
  
ggsave("Maddie_Size.comp.TDGDLF.tiff",width = 15,height = 10,compression = "lzw")


