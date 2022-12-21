
# Resources:
enviro<-read.table("enviromental.txt",header=TRUE, row.names = 1)

land<-(enviro[c(2:8)])
local<-(enviro[c(9:25)])


# Run GLMs to analyze the relationship between land-use and environmental conditions
# (consider regressing environmental metrics against land-use metrics, 
# and show this in the Supplementary Material)
library(lme4)	


options (na.action="na.fail")
fullmodel= glm(width~(native_vegetation  +
                        pasture +
                        agriculture +
                        waterbodies),
                  data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.1=get.models(D_AIC, subset=1)[[1]]
summary(best.1)


visreg(best.1, "native_vegetation", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
visreg(best.1, "native_vegetation", by = "width", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
(deviance(best.1) / df.residual(best.1)) 
Anova(best.1)
vif(best.1)
simulationOutput <- simulateResiduals(fittedModel = best.1, plot = TRUE)
check_overdispersion(best.1)
testZeroInflation(best.1) 



fullmodel= glm(depth~(native_vegetation  +
                         pasture +
                         agriculture +
                         waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.2=get.models(D_AIC, subset=1)[[1]]
summary(best.2)
visreg(best.2, "pasture", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
visreg(best.2, "pasture", by = "depth", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
(deviance(best.2) / df.residual(best.2)) 
Anova(best.2)
vif(best.2)
simulationOutput <- simulateResiduals(fittedModel = best.2, plot = TRUE)
check_overdispersion(best.2)
testZeroInflation(best.2) 


fullmodel= glm(oxygen~(native_vegetation  +
                         pasture +
                         agriculture +
                         waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.3=get.models(D_AIC, subset=1)[[1]]            # NO TIENE PREDICTORES
summary(best.3)
visreg(best.3)
(deviance(best.3) / df.residual(best.3)) 
Anova(best.3)
vif(best.3)
simulationOutput <- simulateResiduals(fittedModel = best.3, plot = TRUE)
check_overdispersion(best.3)
testZeroInflation(best.3) 




fullmodel= glm(temperature~(native_vegetation  +
                              pasture +
                              agriculture +
                              waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.4=get.models(D_AIC, subset=1)[[1]]
summary(best.4)
visreg(best.4, "agriculture", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
visreg(best.4, "agriculture", by = "temperature", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
visreg(best.4, "native_vegetation", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
visreg(best.4, "native_vegetation", by = "temperature", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
(deviance(best.4) / df.residual(best.4))
Anova(best.4)
vif(best.4)
simulationOutput <- simulateResiduals(fittedModel = best.4, plot = TRUE)
check_overdispersion(best.4)
testZeroInflation(best.4) 


fullmodel= glm(conductivity~(native_vegetation  +
                               pasture +
                               agriculture +
                               waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.5=get.models(D_AIC, subset=1)[[1]]            # NO TIENE PREDICTORES
summary(best.5)
visreg(best.5)


fullmodel= glm(turbidity~(native_vegetation  +
                            pasture +
                            agriculture +
                            waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.6=get.models(D_AIC, subset=1)[[1]]       # NO TIENE PREDICTORES     
summary(best.6)
visreg(best.6)



fullmodel= glm(STS~(native_vegetation  +
                      pasture +
                      agriculture +
                      waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.7=get.models(D_AIC, subset=1)[[1]]      # NO TIENE PREDICTORES    
summary(best.7)
visreg(best.7)


fullmodel= glm(pH~(native_vegetation  +
                     pasture +
                     agriculture +
                     waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.8=get.models(D_AIC, subset=1)[[1]]         
summary(best.8)
visreg(best.8, "agriculture", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
visreg(best.8, "agriculture", by = "pH", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
(deviance(best.8) / df.residual(best.8)) 
Anova(best.8)
vif(best.8)
simulationOutput <- simulateResiduals(fittedModel = best.8, plot = TRUE)
check_overdispersion(best.8)
testZeroInflation(best.8) 


fullmodel= glm(chlorophyll~(native_vegetation  +
                              pasture +
                              agriculture +
                              waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.9=get.models(D_AIC, subset=1)[[1]]         
summary(best.9)                              # NO TIENE PREDICTORES 
visreg(best.9)


fullmodel= glm(CPOM~(native_vegetation  +
                       pasture +
                       agriculture +
                       waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.10=get.models(D_AIC, subset=1)[[1]]         
summary(best.10)                           # No fue significativo
visreg(best.10)
(deviance(best.10) / df.residual(best.10)) 
Anova(best.10)
vif(best.10)
simulationOutput <- simulateResiduals(fittedModel = best.10, plot = TRUE)
check_overdispersion(best.10)
testZeroInflation(best.10)


fullmodel= glm(macrophytes~(native_vegetation  +
                              pasture +
                              agriculture +
                              waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.11=get.models(D_AIC, subset=1)[[1]]        # NO TIENE PREDICTORES     
summary(best.11)
visreg(best.11)


fullmodel= glm(submerged_riparian~(native_vegetation  +
                                     pasture +
                                     agriculture +
                                     waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.12=get.models(D_AIC, subset=1)[[1]]        # NO TIENE PREDICTORES     
summary(best.12)
visreg(best.12)




fullmodel= glm(clay~(native_vegetation  +
                       pasture +
                       agriculture +
                       waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.13=get.models(D_AIC, subset=1)[[1]]        # NO TIENE PREDICTORES     
summary(best.13)
visreg(best.13)


fullmodel= glm(sand~(native_vegetation  +
                       pasture +
                       agriculture +
                       waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.14=get.models(D_AIC, subset=1)[[1]]             
summary(best.14)
visreg(best.14, "agriculture", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
visreg(best.14, "agriculture", by = "sand", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
(deviance(best.14) / df.residual(best.14)) 
Anova(best.14)
vif(best.14)
simulationOutput <- simulateResiduals(fittedModel = best.14, plot = TRUE)
check_overdispersion(best.14)
testZeroInflation(best.14)



fullmodel= glm(gravel~(native_vegetation  +
                         pasture +
                         agriculture +
                         waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.15=get.models(D_AIC, subset=1)[[1]]             
summary(best.15)                               # No es significativo 
visreg(best.15)
(deviance(best.15) / df.residual(best.15)) 
Anova(best.15)
vif(best.15)
simulationOutput <- simulateResiduals(fittedModel = best.15, plot = TRUE)
check_overdispersion(best.15)
testZeroInflation(best.15)



fullmodel= glm(pebble~(native_vegetation  +
                         pasture +
                         agriculture +
                         waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.16=get.models(D_AIC, subset=1)[[1]]             
summary(best.16)
visreg(best.16, "waterbodies", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
visreg(best.16, "waterbodies", by = "pebble", line=list(col="black",cex=0.5), points=list(cex=1, pch=19, col="blue"))
(deviance(best.16) / df.residual(best.16)) 
Anova(best.16)
vif(best.16)
simulationOutput <- simulateResiduals(fittedModel = best.16, plot = TRUE)
check_overdispersion(best.16)
testZeroInflation(best.16)



fullmodel= glm(boulder~(native_vegetation  +
                          pasture +
                          agriculture +
                          waterbodies),
               data = enviro, family = "gaussian")

D1 = dredge(fullmodel,rank = "AIC", m.max=10)
D_AIC<-D1[1:5,]
head(D_AIC)
best.17=get.models(D_AIC, subset=1)[[1]]             
summary(best.17)                     # NO TIENE PREDICTORES 
visreg(best.17)



#-----Solo las variables que el RDA indica que son afectadas por el uso de la tierra:

mod_1= glm(width~(native_vegetation  +
                        pasture),
               data = enviro, family = "gaussian")

summary(mod_1)
Anova(mod_1)
vif(mod_1)
simulationOutput <- simulateResiduals(fittedModel = mod_1, plot = TRUE)


mod_2= glm(temperature~(native_vegetation  +
                              pasture),
               data = enviro, family = "gaussian")

summary(mod_2)
Anova(mod_2)
vif(mod_2)
simulationOutput <- simulateResiduals(fittedModel = mod_2, plot = TRUE)

mod_3= glm(STS~(native_vegetation  +
                          pasture),
           data = enviro, family = "gaussian")

summary(mod_3)
Anova(mod_3)
vif(mod_3)
simulationOutput <- simulateResiduals(fittedModel = mod_3, plot = TRUE)


mod_4= glm(conductivity~(waterbodies  +
                  pasture),
           data = enviro, family = "gaussian")

summary(mod_4)
Anova(mod_4)
vif(mod_4)
simulationOutput <- simulateResiduals(fittedModel = mod_4, plot = TRUE)


mod_5= glm(CPOM~(waterbodies  +
                           pasture),
           data = enviro, family = "gaussian")

summary(mod_5)
Anova(mod_5)
vif(mod_5)
simulationOutput <- simulateResiduals(fittedModel = mod_5, plot = TRUE)


mod_6= glm(sand~(waterbodies  +
                   pasture),
           data = enviro, family = "gaussian")

summary(mod_6)
Anova(mod_6)
vif(mod_6)
simulationOutput <- simulateResiduals(fittedModel = mod_6, plot = TRUE)


ggplot(enviro, aes(pasture, sand)) +
  geom_point(pch = 21, size = 2, alpha = 0.9,
             fill = "grey") +
  geom_smooth(method = "glm") +
  labs(x = "Pasture",
       y = "Sand")+
  theme_bw(base_size = 12)



mod_7= glm(pebble~(waterbodies  +
                     pasture),
           data = enviro, family = "gaussian")

summary(mod_7)
Anova(mod_7)
vif(mod_7)
simulationOutput <- simulateResiduals(fittedModel = mod_7, plot = TRUE)

ggplot(enviro, aes(waterbodies, pebble)) +
  geom_point(pch = 21, size = 2, alpha = 0.9,
             fill = "grey") +
  geom_smooth(method = "glm") +
  labs(x = "Waterbodies",
       y = "Pebble")+
  theme_bw(base_size = 12)



mod_8= glm(boulder~(waterbodies  +
                     pasture),
           data = enviro, family = "gaussian")

summary(mod_8)
Anova(mod_8)
vif(mod_8)
simulationOutput <- simulateResiduals(fittedModel = mod_8, plot = TRUE)

















