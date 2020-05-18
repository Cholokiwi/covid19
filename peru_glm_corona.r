#################################################
# Scripts de R para reproducir los analisis de
# Costilla y Rojas (2020)
# Predicción de corto plazo del número de fallecimientos por COVID19 en Perú:
# enfoque usando modelos lineales generalizados
#################################################

# Setting up
rm(list=ls())
set.seed(12345678)
require(MASS)
require(tidyverse)
require(readxl)
require(ggplot2)
require(glm.predict)

# Reading data directly from repository
reportes_minsa.xlsx =tempfile()
download.file("https://github.com/jincio/COVID_19_PERU/blob/master/docs/reportes_minsa.xlsx?raw=true", "reportes_minsa.xlsx", mode="wb")
alldata=read_excel("reportes_minsa.xlsx", sheet = "Sheet2")
str(alldata)
alldata$Dia[nrow(alldata)]
dataperu=data.frame(alldata[,
              c("Dia","Fallecidos","Hospitalizados","Pruebas_dia","TasaPositivos")])
colnames(dataperu)=c("date","ycum","hosp","pruebasd","tpositivos")

# Plotting Parameters
myvar="y"
# Uncomment line below for models with moving averages
#myvar="y_3dma"
maxpred=ifelse(myvar=="y",160,125)

# subsetting data from first deaths (19 March)
dataperu=dataperu[!is.na(dataperu$ycum),]
dataperu$y=dataperu$ycum-lag(dataperu$ycum, 1)
dataperu$y[1]=dataperu$ycum[1]
# until 16 May
dataperu <- dataperu %>% filter(date<"2020-05-17")

# set prediction date (7 days in future)
npred=dim(dataperu)[1]+7

# 3d ma
dataperu$y_3dma=round((dataperu$y+lag(dataperu$y,1)+lag(dataperu$y,2))/3,0)

# temporal variables
dataperu$t=1:nrow(dataperu)
dataperu$t2=(1:nrow(dataperu))^2
dataperu$t3=(1:nrow(dataperu))^3
dataperu$t4=(1:nrow(dataperu))^4

# weekdays effect for Poisson
dataperu$date=as.Date(dataperu$date)
dataperu$weekdays=as.factor(weekdays(dataperu$date, abbr=T))
dataperu$dia = recode_factor(dataperu$weekdays, 
                           Fri = "Viernes",
                           Sat = "Sabado",
                           Sun = "Domingo",
                           Mon = "Lunes",
                           Tue = "Martes",
                           Wed = "Miercoles",
                           Thu = "Jueves",
                           .default = levels(dataperu$weekdays))

#####   Poisson
mymodel=glm(get(myvar) ~ 1+t+t2+dia,
            data=dataperu, family=poisson)  
summary(mymodel)


# out-of-sample prediction
myt=seq(1,npred, 1)
newdates=as.Date(seq(dataperu$date[1], dataperu$date[1]+npred-1, by=1))
newweekend=ifelse(weekdays(newdates, abbr=T)%in%c("Sun","Mon"),1,0)
newweekdays=as.factor(weekdays(newdates, abbr=T))
newdata=data.frame(1,t=myt, t2=myt^2, t3=myt^3, t4=myt^4, 
                   weekend=newweekend, weekdays=newweekdays)
newdata$dia = recode_factor(newdata$weekdays, 
                             Fri = "Viernes",
                             Sat = "Sabado",
                             Sun = "Domingo",
                             Mon = "Lunes",
                             Tue = "Martes",
                             Wed = "Miercoles",
                             Thu = "Jueves",
                             .default = levels(newdata$weekdays))
newdata$date=dataperu$date[1]+(newdata$t-1)

# Prediction
mypred=predict.glm(mymodel,type="response", se.fit = T,
                   newdata = newdata)
newdata$mypred=mypred$fit
newdata$mypred.min=ifelse(mypred$fit-1.96*mypred$se.fit>0,mypred$fit-1.96*mypred$se.fit,0)
newdata$mypred.max=mypred$fit+1.96*mypred$se.fit

# max point
yhatmax=which(mypred$fit==max(mypred$fit[!is.na(mypred$fit)]))
as.Date(dataperu$date[1])+as.numeric(yhatmax)

# plot out-of-sample
plt <- ggplot(newdata,aes(x = date, y = mypred)) +
  geom_line(color="red", size=1.5)+
  geom_ribbon(aes(ymin = mypred.min, ymax = mypred.max),
              alpha = 0.15, size=2)+
  labs(x = 'Fecha', y = myvar)+
  ylim(0,maxpred) + xlim(as.Date(c("2020-03-15",max(newdata$date))))+
  theme_bw()+ 
  geom_point(data=dataperu, aes(x=date, y=get(myvar)), size=1.25)
plt

pdf(paste0("poisson_peru_",
           myvar,'_',paste(attr(mymodel$terms,"term.labels"), collapse = '_'),
           "_n",nrow(dataperu),"_",npred-nrow(dataperu),"days.pdf"), 
    height = 3, width = 4.5)
plt
dev.off()

# GOF competing models
mycovars=c("1","t","t2","t3","t4")
mygof=matrix(NA,nrow=length(mycovars),ncol=4, dimnames=list(mycovars, c("npar","DF","AIC","BIC")))
thiscovars=NULL
for (mycovar in mycovars){
  if (mycovar=="1") thiscovars=paste(mycovar,"dia",sep="+") else
    thiscovars=paste(thiscovars, mycovar, sep = "+")
  myformula=as.formula(paste(myvar, " ~ ", thiscovars))
  message("Model: ",myformula, '\n')
  mymodel=glm(myformula, 
                    data=dataperu,family=poisson)
  npar=nrow(dataperu)-mymodel$df.residual
  mygof[mycovar,]=c(npar,mymodel$df.residual,AIC(mymodel),BIC(mymodel))
  rownames(mygof)[which(rownames(mygof)==mycovar)]=thiscovars
}
round(mygof,1)
write.csv(mygof,paste0("poi_gof_",myvar,"_",thiscovars,".csv"), quote=F)


#############  NB
mymodel.nb=glm.nb(get(myvar)~1+t+t2+dia,
                data=dataperu)
rbind(summary(mymodel.nb)$coef, 
      theta=c(summary(mymodel.nb)$theta,summary(mymodel.nb)$SE.theta,NA,NA))
write.csv(summary(mymodel.nb)$coef,paste0("nb_coef_",myvar,".csv"), quote=F)


# out-of-sample prediction
myt=seq(1,npred, 1)
newdates=as.Date(seq(dataperu$date[1], dataperu$date[1]+npred-1, by=1))
newweekend=ifelse(weekdays(newdates, abbr=T)%in%c("Sun","Mon"),1,0)
newweekdays=as.factor(weekdays(newdates, abbr=T))
newdata=data.frame(1,t=myt, t2=myt^2, t3=myt^3, t4=myt^4, 
                   weekend=newweekend, weekdays=newweekdays)
newdata$dia = recode_factor(newdata$weekdays, 
                            Fri = "Viernes",
                            Sat = "Sabado",
                            Sun = "Domingo",
                            Mon = "Lunes",
                            Tue = "Martes",
                            Wed = "Miercoles",
                            Thu = "Jueves",
                            .default = levels(newdata$weekdays))

newdata$date=dataperu$date[1]+(newdata$t-1)

# Prediction
mypred=predict.glm(mymodel.nb,type="response", se.fit = T,
                   newdata = newdata)
# prediction
newdata$mypred=mypred$fit
newdata$mypred.min=ifelse(mypred$fit-1.96*mypred$se.fit>0,mypred$fit-1.96*mypred$se.fit,0)
newdata$mypred.max=mypred$fit+1.96*mypred$se.fit
# max point
yhatmax=which(mypred$fit==max(mypred$fit[!is.na(mypred$fit)]))
newdata[yhatmax,]
as.Date(dataperu$date[1])+as.numeric(yhatmax)

plt <- ggplot(newdata,aes(x = date, y = mypred)) +
  geom_line(color="red", size=1.5)+
  geom_ribbon(aes(ymin = mypred.min, ymax = mypred.max),
              alpha = 0.15)+
  labs(x = 'Fecha', y = myvar)+
  ylim(0,maxpred) + xlim(as.Date(c("2020-03-15",max(newdata$date))))+
  theme_bw()+ 
  geom_point(data=dataperu, aes(x=date, y=get(myvar)), size=1.25)
plt

pdf(paste0("nb_peru_",myvar,'_',
           paste(attr(mymodel.nb$terms,"term.labels"), collapse = '_'),
           "_n",nrow(dataperu),"_",npred-nrow(dataperu),"days.pdf"), 
    height = 3, width = 4.5)
plt
dev.off()

# GOF NB
mycovars=c("1","t","t2","t3","t4")
mygof.nb=matrix(NA,nrow=length(mycovars),ncol=4, dimnames=list(mycovars, c("npar","DF","AIC","BIC")))
thiscovars=NULL
for (mycovar in mycovars){
#      if (myvar=="1") mycovars=myvar else
      if (mycovar=="1") thiscovars=paste(mycovar,"dia",sep="+") else
      thiscovars=paste(thiscovars, mycovar, sep = "+")
  myformula=as.formula(paste(myvar, " ~ ", thiscovars))
  message("Model: ",myformula, '\n')
  mymodel.nb=glm.nb(myformula, 
                    data=dataperu)
  theta=c(summary(mymodel.nb)$theta, summary(mymodel.nb)$SE.theta, NA,NA)
  npar=nrow(dataperu)-mymodel.nb$df.residual+1
  mygof.nb[mycovar,]=c(npar,mymodel.nb$df.residual,AIC(mymodel.nb),BIC(mymodel.nb))
  rownames(mygof.nb)[which(rownames(mygof.nb)==mycovar)]=thiscovars
}
write.csv(mygof.nb,paste0("nb_gof_",myvar,"_",thiscovars,".csv"), quote=F)

cat("GOF Poisson \n")
round(mygof,1)

cat("GOF Negative Binomial \n")
round(mygof.nb,1)
