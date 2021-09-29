## Cleaning environment:
rm(list = ls())

#########################################################
### Multivariable multinomial logit with random sampling:
#########################################################

## Dataset used:
# Students:  information about students and households (id: mrun)
# opciones:  data set containing all schools for all households (1211x80772)

## For loop:
j <- 1000 # iteraciones
k <- 1    # otras X escuelas no aleatorias
l <- 10   # choice set por alumno
sf <- c(paste0('000',1:9),paste0('00',10:99),paste0('0',100:999),1000)

## Set memoria limits R:
memory.limit(size = 1000000000)
options(scipen=999)

## Packages:
{
  library(data.table)
  library(dplyr)
  library(mlogit)
  library(survey)
}  

## Formula:
{
  f_MNL_M1 <- mFormula(choice ~ distancia + SIMCE_mate + SNI + valor_copago +
                                porc_prioritarios + AXS + orientacion + regimen - 1)
  f_MNL_M2 <- mFormula(choice ~ distancia:prioritario + SIMCE_mate:prioritario + SNI:prioritario + 
                                valor_copago:prioritario + porc_prioritarios:prioritario + 
                                AXS:prioritario + orientacion:prioritario + regimen:prioritario - 1)
  f_MNL_M3 <- mFormula(choice ~ distancia:GSE + SIMCE_mate:GSE + SNI:GSE + valor_copago:GSE +
                                porc_prioritarios:GSE + AXS:GSE + orientacion:GSE + regimen:GSE - 1)
  alt.names <- c(paste0('Pref',1:k),paste0('Rand',(k+1):l))
}

## Modifying dataset:
{
  students <- fread('03 DATABASES FOR MODELS/STUDENTS.csv')
  ## GSE of household's neighnorhood
  students[soc_dim < 0.5, GSE := 50]
  students[soc_dim >= 0.5 & soc_dim < 0.7, GSE := 70]
  students[soc_dim > 0.7, GSE := 100]
  students[,GSE := as.factor(GSE)]
  ## Priority student
  students[,prioritario := as.factor(prioritario)]
  
  opciones <- fread('03 DATABASES FOR MODELS/MOST_PREFERRED_SCHOOL.csv')
  opciones[,choice:=ifelse(choice==1,1,0)]
  names(opciones)[names(opciones) == 'Promedio alumnos por curso'] <- 'AXS'
  choices <- opciones[choice == T]
  choices <- cbind(choices,alt=c(paste0('Pref',1:k)))
  choices <- cbind(choices,order=1:k)
  opciones <- opciones[choice != 1]
}

## Loop:

for(i in 1:j){
  ## Bases
  RBDS <- NULL
  COEFICIENTES_M1 <- NULL
  COEFICIENTES_M2 <- NULL
  COEFICIENTES_M3 <- NULL
  INTERVALS_M1 <- NULL
  INTERVALS_M2 <- NULL
  INTERVALS_M3 <- NULL
  LGL_M1 <- NULL
  LGL_M2 <- NULL
  LGL_M3 <- NULL
  TABLA <- NULL
  
  print(paste0('Iteration N°',i))
  base <- opciones[,.SD[sample(.N, min(l-k,.N))],by = mrun]
  base <- cbind(base,alt=paste0('Rand',(k+1):l))
  base <- cbind(base,order=(k+1):l)
  base <- rbind(choices,base)
  base[,valor_copago := valor_copago/1000][,SNI := SNI*100]
  base <- merge(base,students[,.(mrun,prioritario,GSE)],by='mrun',all.x=T)
  base <- base[order(mrun,order)]
  
  RBDS <- rbind(RBDS,data.table(iteracion=i,base[,.(mrun,rbd,alt)]))
  
  base <- data.frame(base)
  base <- mlogit.data(data = base,choice = 'choice',shape = 'long',id.var = 'mrun',
                      alt.levels = alt.names)
  
  print(paste0('Iteration N°',i,', models...'))
  reg_MNL_M1 <- mlogit(formula = f_MNL_M1,data = base)
  reg_MNL_M2 <- mlogit(formula = f_MNL_M2,data = base)
  reg_MNL_M3 <- mlogit(formula = f_MNL_M3,data = base)
  
  print(paste0('Iteration N°',i,', Wald Test...'))
  aM2 <- seq(1,16,2)
  bM2 <- seq(2,16,2)
  aM3 <- seq(1,24,3)
  bM3 <- seq(2,24,3)
  cM3 <- seq(3,24,3)
  for(z in 1:8){
    coefs <- names(reg_MNL_M1$coefficients)
    pM2 <- names(reg_MNL_M2$coefficients)
    pM3 <- names(reg_MNL_M3$coefficients)
    lM2 <- linearHypothesis(reg_MNL_M2,paste0(pM2[aM2[z]],' = ',pM2[bM2[z]]),test="Chisq")
    lM3 <- linearHypothesis(reg_MNL_M3,c(paste0(pM3[aM3[z]],' = ',pM3[bM3[z]]),
                                         paste0(pM3[aM3[z]],' = ',pM3[cM3[z]])),test="Chisq")
    TABLA <- rbind(TABLA,data.table(coefficients = coefs[z],
                                    `P-value M2` = lM2$`Pr(>Chisq)`[2],
                                    `P-value M3` = lM3$`Pr(>Chisq)`[2]))
  }
  
  print(paste0('Iteration N°',i,', intervals...'))
  e <- confint(reg_MNL_M1,level = 0.999)
  f <- confint(reg_MNL_M1,level = 0.99)
  g <- confint(reg_MNL_M1,level = 0.95)
  m <- confint(reg_MNL_M2,level = 0.999)
  n <- confint(reg_MNL_M2,level = 0.99)
  o <- confint(reg_MNL_M2,level = 0.95)
  p <- confint(reg_MNL_M3,level = 0.999)
  q <- confint(reg_MNL_M3,level = 0.99)
  r <- confint(reg_MNL_M3,level = 0.95)
  
  print(paste0('Iteration N°',i,', coefficients...'))
  s <- data.table(iteracion=i,variable=row.names(data.frame(summary(reg_MNL_M1)$CoefTable)),
                  data.frame(summary(reg_MNL_M1)$CoefTable))
  t <- data.table(iteracion=i,variable=row.names(data.frame(summary(reg_MNL_M1)$CoefTable)),
                  cbind(data.table(e),data.table(f),data.table(g)))
  u <- data.table(iteracion=i,variable=row.names(data.frame(summary(reg_MNL_M2)$CoefTable)),
                  data.frame(summary(reg_MNL_M2)$CoefTable))
  v <- data.table(iteracion=i,variable=row.names(data.frame(summary(reg_MNL_M2)$CoefTable)),
                  cbind(data.table(m),data.table(n),data.table(o)))
  w <- data.table(iteracion=i,variable=row.names(data.frame(summary(reg_MNL_M3)$CoefTable)),
                  data.frame(summary(reg_MNL_M3)$CoefTable))
  x <- data.table(iteracion=i,variable=row.names(data.frame(summary(reg_MNL_M3)$CoefTable)),
                  cbind(data.table(p),data.table(q),data.table(r)))
  
  print(paste0('Iteration N°',i,', Log-likelihood...'))
  LGL_M1 <- data.table(iteracion=i,`Log_Likelihood`=reg_MNL_M1$logLik)
  LGL_M2 <- data.table(iteracion=i,`Log_Likelihood`=reg_MNL_M2$logLik)
  LGL_M3 <- data.table(iteracion=i,`Log_Likelihood`=reg_MNL_M3$logLik)
  
  print(paste0('Iteration N°',i,', saving data...'))
  fwrite(RBDS,paste0('01 Data/04 Loop/RBDs/RBDs ',sf[i],'.csv'))
  fwrite(s,paste0('01 Data/04 Loop/Coeficientes/Coefs M1 ',sf[i],'.csv'))
  fwrite(u,paste0('01 Data/04 Loop/Coeficientes/Coefs M2 ',sf[i],'.csv'))
  fwrite(w,paste0('01 Data/04 Loop/Coeficientes/Coefs M3 ',sf[i],'.csv'))
  fwrite(t,paste0('01 Data/04 Loop/Intervalos/Intervalos M1 ',sf[i],'.csv'))
  fwrite(v,paste0('01 Data/04 Loop/Intervalos/Intervalos M2 ',sf[i],'.csv'))
  fwrite(x,paste0('01 Data/04 Loop/Intervalos/Intervalos M3 ',sf[i],'.csv'))
  fwrite(LGL_M1,paste0('01 Data/04 Loop/Log-likelihood/Log-likelihood M1 ',sf[i],'.csv'))
  fwrite(LGL_M2,paste0('01 Data/04 Loop/Log-likelihood/Log-likelihood M2 ',sf[i],'.csv'))
  fwrite(LGL_M3,paste0('01 Data/04 Loop/Log-likelihood/Log-likelihood M3 ',sf[i],'.csv'))
  fwrite(TABLA,paste0('01 Data/04 Loop/Wald Test/Wald Test ',sf[i],'.csv'))
  
  rm(base,reg_MNL_M1,reg_MNL_M2,reg_MNL_M3,e,f,g,m,n,o,p,q,r,s,t,u,v,w,x,coefs,pM2,pM3,lM2,lM3)
}
