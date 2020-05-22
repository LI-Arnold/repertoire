
source("~/TER/Scripts/initialiserSansAnomalies.r")

### Charger les données du questionnaire
Questionnaire<-read.csv("questionnaire_participants.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)

###-------------------------------- Analyse de variance entre PM2.5 Eet (activity,event)----------------------------------------

### table avec les mesures de PM2.5 pour chaque activité et événement
 reqDf<-" SELECT \"PM2.5\",activity,event
  		FROM Questionnaire,dfSansAnomalies 
  		WHERE Questionnaire.participant_virtual_id = dfSansAnomalies.participant_virtual_id";
dfVariancePM2.5<-sqldf(reqDf)

### Convertir en factor les variables qualitatives
dfVariancePM2.5$activity <- as.factor(dfVariancePM2.5$activity)
dfVariancePM2.5$event <- as.factor(dfVariancePM2.5$event)

### Appliquer modele Anova à deux facteurs (event et activity)
modelePM2.5<- aov( PM2.5 ~ activity + event + activity:event,data=dfVariancePM2.5)

### Affichage les resultats de l'anova
summary(modelePM2.5)

###-------------------------------- Analyse de variance entre PM2.5 Eet (activity,event)---------------------------------------
 reqDf<-" SELECT NO2,activity,event
  		FROM Questionnaire,dfSansAnomalies 
  		WHERE Questionnaire.participant_virtual_id = dfSansAnomalies.participant_virtual_id";
dfVarianceNO2<-sqldf(reqDf)
### Convertir en factor les variables qualitatives
dfVarianceNO2$activity <- as.factor(dfVarianceNO2$activity)
dfVarianceNO2$event <- as.factor(dfVarianceNO2$event)
modelePM2.5<- aov( NO2 ~ activity + event + activity:event,data=dfVarianceNO2)
summary(modeleNO2)