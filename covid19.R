library(ggplot2)
library(nlme)
library(GGally)
library("RColorBrewer")
library(lmtest)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 100))

covid<-read.table("cases",head=T,sep="\t")

#WHO BCG immunization rates, assuming NAs=0
bcg<-read.table("BCGrates.csv",head=T,sep="\t",quote="")
bcg$BCG[is.na(bcg$BCG)]<-0

#WHO Rubella immunization rates, assuming NAs=0
rubella<-read.table("Rubellarates.csv",head=T,sep="\t",quote="")
rubella$Rubella[is.na(rubella$Rubella)]<-0

#From CIA data via Wikipedia
age<-read.table("age",head=T,sep="\t")

#From CIA data via Wikipedia
population<-read.table("population",head=T,sep="\t")

#Original policy data from medarxiv publication data table https://www.researchgate.net/publication/340263333_Correlation_between_universal_BCG_vaccination_policy_and_reduced_morbidity_and_mortality_for_COVID-19_an_epidemiological_study
publication<-read.table("BCG.csv",sep="\t",head=T,quote="",fill=T)
publication$X<-NULL #Cleanup empy column

BMI<-read.table("BMI",head=T,sep="\t",quote="")
BMI$BMI<-as.character(BMI$BMI)
BMI$BMI[BMI$BMI=="-"]<-0
BMI$BMI<-as.numeric(BMI$BMI)

#Merge
agecovid<-merge(covid,age,by.x="X",by.y="Country")
#agecovid<-merge(agecovid,BMI,by.x="X",by.y="Country")
agecovidpop<-merge(agecovid,population,by.x="X",by.y="Country")
agecovidpop<-merge(agecovidpop,rubella,by.x="X",by.y="Country",all.x=T)
agecovidpop<-merge(agecovidpop,bcg,by.x="X",by.y="Country",all.x=T)
agecovidpop$BCG[is.na(agecovidpop$BCG)]<-0
agecovidpop$Rubella[is.na(agecovidpop$Rubella)]<-0

agecovidpop$Population2018<-as.numeric(as.character(gsub(",","",agecovidpop$Population2018)))
agecovidpop$Cases.b.<-as.numeric(as.character(gsub(",","",agecovidpop$Cases.b.)))
agecovidpop$CasesPerM<-(agecovidpop$Cases.b.*1000000)/agecovidpop$Population2018
agecovidpop$Deaths.c.<-as.numeric(as.character(gsub(",","",agecovidpop$Deaths.c.)))
agecovidpop$Recov..d.<-as.numeric(as.character(gsub(",","",agecovidpop$Recov..d.)))
#Recoverd-NA warning- these are legitimate 0s marked with -
#Only use countries wiht >1M people
agecovidpopM<-agecovidpop[agecovidpop$Population2018>1000000,]

#Get BMI data
covid.stage1<-merge(agecovidpopM,BMI,by.x="X",by.y="Country")
covid.stage1$DeathPerM<-(covid.stage1$Deaths.c.*1000000)/covid.stage1$Population2018

#Correlation between cases and BMI per country, add pseudocount to account for 0s
cor(log(covid.stage1$CasesPerM),covid.stage1$BMI)
cor(log(covid.stage1$Deaths.c.+0.00000001),covid.stage1$BMI)

#Add publication categories
covid.stage2<-merge(covid.stage1,publication,by.x="X",by.y="Country")
covid.stage2$Policy<-as.integer(as.character(covid.stage2$Policy))
covid.stage2$IncomeLevel<-as.factor(as.character(covid.stage2$IncomeLevel))
covid.stage2$rubella50<-"Low"
covid.stage2$rubella50[covid.stage2$Rubella>50]<-"High"
covid.stage2$rubella50<-as.factor(covid.stage2$rubella50)

#Correlations between infection rates and Median Age/BCG policy
cor(log(covid.stage2$CasesPerM),covid.stage2$Median)
cor(log(covid.stage2$CasesPerM),covid.stage2$Policy)

#Convert to factor the policy
covid.stage2$Policy<-as.factor(as.character(covid.stage2$Policy))
covid.stage2$CasesPerMLog<-log(covid.stage2$CasesPerM)

fullm<-lm(log(CasesPerM) ~  Median + Policy   + IncomeLevel , data= covid.stage2)
summary(fullm)
policylm<-lm(log(CasesPerM) ~  Policy + IncomeLevel, data= covid.stage2)
summary(policylm)
medianlm<-lm(log(CasesPerM) ~  Median + IncomeLevel , data= covid.stage2)

#Log likelyhood for the full model and reduced model without policy
lrtest(fullm,medianlm)

#Same for BCG rates
bcgfull<-lm(log(CasesPerM) ~  Median + BCG   + IncomeLevel , data= covid.stage2)
summary(bcgfull)
lrtest(bcgfull,medianlm)

#Plots
covid.stage2$CasesPerMLog<-log(covid.stage2$CasesPerM)

g1a<-ggplot(covid.stage2,aes(x=Median,y=CasesPerMLog,col=Policy))+geom_point(size=2) +
  geom_smooth(method='lm') + xlab("Median Age")+ ylab("Cases per million, log transformed")  + 
  ggtitle("Correlation between normalized infections per country and the median age(years)") +
  theme(text = element_text(size=10), axis.text.x = element_text( hjust=1)) 
ggsave("MedianAgePolicyPoints.png",plot=g1a,width=8,height = 6)

g1b<-ggplot(covid.stage2,aes(x=Median,y=CasesPerMLog,col=BCG))+geom_point(size=2) + sc +
  geom_smooth(method='lm') + xlab("Median Age")+ ylab("Cases per million, log transformed")  + 
  ggtitle("Correlation between normalized infections per country and the median age(years)") +
  theme(text = element_text(size=10), axis.text.x = element_text( hjust=1)) 
ggsave("MedianAgeBCGRatesPoints.png",plot=g1b,width=8,height = 6)

g1c<-ggplot(covid.stage2,aes(x=Median,y=CasesPerMLog,col=Rubella))+geom_point(size=2) + sc +
  geom_smooth(method='lm') + xlab("Median Age")+ ylab("Cases per million, log transformed")  + 
  ggtitle("Correlation between normalized infections per country and the median age(years)") +
  theme(text = element_text(size=10), axis.text.x = element_text( hjust=1)) 
ggsave("MedianAgeRubellaRatesPoints.png",plot=g1c,width=8,height = 6)

#Rubella associations
anova(lm(log(CasesPerM) ~  rubella50 , data= covid.stage2))
rubellalm<-lm(log(CasesPerM) ~  Median + rubella50   + IncomeLevel , data= covid.stage2)
summary(rubellalm)

#Compare BCG policy to median age
anova(lm(log(CasesPerM) ~  Policy + IncomeLevel, data= covid.stage2))
anova(lm(log(CasesPerM) ~  Median + IncomeLevel, data= covid.stage2))
anova(lm(log(CasesPerM) ~  Policy , data= covid.stage2))
anova(lm(log(CasesPerM) ~  Median, data= covid.stage2))

#IncomeLevel as a random factor
randModelPol<-lme(log(CasesPerM) ~  Median + Policy, random= ~1|IncomeLevel , data= covid.stage2)
summary(randModelPol)
randModelBBCG<-lme(log(CasesPerM) ~  Median + BCG, random= ~1|IncomeLevel , data= covid.stage2)
summary(randModelBBCG)


#Suppl plots
g5<-ggplot(covid.stage2,aes(x=Policy,y=Median,col=IncomeLevel))+geom_boxplot() + ylab("Median Age (years)")+ xlab("BCG policy")  + 
  ggtitle("Median age per country according to BCG policy and income level")
ggsave("IncomeMedianAgePolicyBxpl.png",plot=g5,width=6,height = 6)


#Deaths
covid.exploreD<-covid.stage2
covid.exploreD<-covid.exploreD[covid.exploreD$DeathPerM>0,]
covid.exploreD$Policy<-as.factor(as.character(covid.exploreD$Policy))

#We should look at policy as a continuous var
covid.exploreD$Policy<-as.integer(as.character(covid.exploreD$Policy))
g6<-ggplot(covid.exploreD,aes(x=Median,y=log(DeathPerM),col=Policy))+geom_point() +
  geom_smooth(method='lm') + xlab("Median Age")+ ylab("Deaths per million, log transformed")  + 
  ggtitle("Correlation between normalized deaths per country and the median age(years)")
ggsave("DeathsPerPolicyMedianAgePoints.png",plot=g6,width=9,height = 6)
cor(covid.exploreD$Median,log(covid.exploreD$DeathPerM))
covid.exploreD$StartYear<-as.numeric(covid.exploreD$StartYear)
cor(covid.exploreD$StartYear,log(covid.exploreD$DeathPerM))

#BMI-deaths association
covid.exploreD$Policy<-as.factor(as.character(covid.exploreD$Policy))
covid.exploreD$BMIcat<-"High"
covid.exploreD$BMIcat[covid.exploreD$BMI<25]<-"Normal"
covid.exploreD$BMIcat<-as.factor(covid.exploreD$BMIcat)
g7<-ggplot(covid.exploreD,aes(x=BMIcat,y=log(DeathPerM),col=Policy)) + geom_boxplot() + xlab("BMI category")+ylab("Deaths per M, log transformed")
ggsave("BMIDeathsBox.png",plot=g7,width=6,height = 6)
