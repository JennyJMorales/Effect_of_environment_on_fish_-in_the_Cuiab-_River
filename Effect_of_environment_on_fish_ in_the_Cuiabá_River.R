
# Land use effects fish assemblages in tropical savanna streams

# Jenny J. Morales1, Luiza Peluso1, Lúcia Mateus2 and Jerry Penha2 

# 1 Programa de Pós-Graduação em Ecologia e Conservação da Biodiversidade, Instituto de Biociências, Universidade Federal de Mato Grosso, Av. Fernando Correia da Costa 2367, 78060-900 Cuiabá, MT, Brasil. 
# 2 Laboratório de Ecologia e Manejo de Recursos Pesqueiros, Instituto de Biociências, Universidade Federal de Mato Grosso, Av. Fernando Correia da Costa 2367, 78060-900 Cuiabá, MT, Brasil. 
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(adespatial)
library(MASS)
library(vegan)        #RDA
library(car)          #funcion: ANOVA
library(visreg)
library(ggplot2)
library(lme4)	       #para GLM
library(MuMIn)
library(performance) #funcion sobredispersion
library(DHARMa)


##-------REPOSITOTY: Effect_of_the_environment_on_fish_in_the_Cuiabá_River

#####----RESOURCES:

#Species abundance of the 85 species of fish captured:
spe<-read.table("comun_85.txt",header=TRUE, row.names = 1)
spe.H <- decostand (spe, "hellinger")     #with Hellinger transformation

#Biota data such as richness, abundance, number of common species and number of dominant species by site
biota<-read.table("biota.txt",header=TRUE, row.names = 1)

#Environmental variables (land_use: 7 + local: 17) measured in 33 streams:
enviro <-read.table("enviromental.txt",header=TRUE, row.names = 1)
envir <-(enviro [c(-15,-11)])             #Remove pH and water temperature since they have normal distribution
env.log <-log(envir+1)                    #Transform the other variables with log+1
pH <-(enviro[15])                         #Call the variable pH
Tem <-(enviro[11])                        #Call the variable temperature
env<- data.frame(env.log, pH, Tem)        #Unify the transformed environmental variables with pH and temperature
env<-(env[c(-2,-5,-6,-16)])               #Eliminate variables with zero values >90% of the sites, i.e., reforestation, urban, mining and macrophytes.
env_catchment<-(env[c(1:4)])              #Land-use variables
env_site<-(env[c(5:20)])                  #Local variables

 
####-------To test the hypothesis that in-catchment land use structures fish assemblages 
#----------by influencing local stream conditions, we conducted three redundancy analyses (RDA):-------------------------------------------------------------------------------------------------------------------------

# 1) land use variables in a predictors matrix of in-stream environmental variables; 

rda.1 <-rda(env_site~.,env_catchment)      #Redundancy Analysis (RDA)
summary(rda.1)
anova(rda.1,permutations=how(nperm=999))   #Validated by Analysis of Variance (ANOVA)        
(R2<- RsquareAdj(rda.1)$r.squared)         #Coefficient of determination       
(r2adj <- RsquareAdj(rda.1)$adj.r.squared) #Adjusted determination coefficient    
vif.cca(rda.1)                             #To evaluate multicollinearity
site_catchment.for<-forward.sel(env_site,env_catchment,adjR2thresh=r2adj, nperm=9999) #to select the predictors that best explain the variation in the data

plot(rda.1, scaling= 1, main = "Habitat variables ~ catchment variables",
     xlab="RDA 1 = 51.35 %", 
     ylab = "RDA 2 = 26.55 %",
     xlim=c(-25,10), ylim=c(-10,10),
     font = 1,bty = "L", type = "p", family = "serif")

# 2) land use variables as the predictor matrix and fish abundances as the response matrix:

rda.2 <-rda(spe.H~.,env_catchment)         #Redundancy Analysis (RDA)
summary(rda.2)
anova(rda.2,permutations=how(nperm=999))   #Validated by Analysis of Variance (ANOVA)      
(R2<- RsquareAdj(rda.2)$r.squared)         #Coefficient of determination      
(r2adj <- RsquareAdj(rda.2)$adj.r.squared) #Adjusted determination coefficient      
vif.cca(rda.2)                             #To evaluate multicollinearity      

plot(rda.2, scaling=1, main = "Fish abundance ~ catchment variables",
     xlab="", 
     ylab = "",
     xlim=c(-0,10), ylim=c(-10,10),
     font = 1,bty = "L", type = "p", family = "serif")

# 3) habitat variables that are not influenced by land use (previously identified in the first RDA model)
# as predictors of a fish abundance matrix

env_site<-(env_site[c(-1,-4,-6,-8,-11,-13,-14,-16)])  #Eliminate habitat variables influenced by land use,
                                                      #i.e., channel width, conductivity, total suspended solids, 
                                                      #coarse particulate organic matter (CPOM), sand, pebble, boulder and temperature

rda.3 <-rda(spe.H~.,env_site)
summary(rda.3)
anova(rda.3,permutations=how(nperm=999))        
(R2<- RsquareAdj(rda.3)$r.squared)              
(r2adj <- RsquareAdj(rda.3)$adj.r.squared)       
vif.cca(rda.3)                                  
env_site.for<-forward.sel(spe.H,env_site,adjR2thresh=r2adj, nperm=9999)

plot(rda.3, scaling=1, main = "Fish abundance ~ Habitat variables", #Corresponds to Figure 3 of the manuscript  
     xlab="RDA 1 = 14.90 %", 
     ylab = "RDA 2 = 7.89 %",
     xlim=c(-10,15), ylim=c(-10,10),
     font = 1,bty = "L", type = "p", family = "serif",
     cex.lab=1.3, cex.axis=1.3, cex.sub=1.3)


###--------To assess whether anthropogenic land-use activities in the catchment and their associated local variables affected fish species richness 
#---------in these headwater streams, we performed a generalized linear model (GLM):-------------------------------------------------------------------------

matrix<-cbind(biota,enviro)   #unify biota data with environmental variables

options (na.action="na.fail")
fullmodel= glm.nb(Richness~(native_vegetation  +
                       pasture + 
                       agriculture +
                       oxygen +
                       STS + 
                       turbidity +
                       CPOM +
                       submerged_riparian +
                       clay +
                       temperature),data = matrix)


D1 = dredge(fullmodel,rank = "AIC", m.max=10) #automatically selects the best models according to AIC's criteria (more parsimonious)
D_AIC<-D1[1:5,]
head(D_AIC)
best=get.models(D_AIC, subset=1)[[1]]
summary(best)
Anova(best)
vif(best)
simulationOutput <- simulateResiduals(fittedModel = best, plot = TRUE)
check_overdispersion(best)
testZeroInflation(best) 
visreg(best)  #for displaying the results of a fitted model
visreg(best, "agriculture", line=list(col="black"), points=list(cex=1, pch=1))
visreg(best, "native_vegetation", line=list(col="black"), points=list(cex=1, pch=1))
visreg(best, "oxygen", line=list(col="black"), points=list(cex=1, pch=1))
visreg(best, "STS", line=list(col="black"), points=list(cex=1, pch=1))
visreg(best, "submerged_riparian", line=list(col="black"), points=list(cex=1, pch=1))
visreg(best, "turbidity", line=list(col="black"), points=list(cex=1, pch=1))





















