library( car )
library( ellipse )
library( faraway )
library( leaps )
library(MASS)
library( GGally)
library( rms )
library(arm)
library(ResourceSelection)
library(pROC)
library(PRROC)
options(rgl.debug=TRUE)
library(rgl)
library(BAS)

dataset <- read.csv("C:/Users/claca/Desktop/OneDrive - Politecnico di Milano/Inferenza statistica/R/classification/diabetes.csv")

# guardiamo le caratteristiche generali del dataset
head(dataset)
View(dataset)
dim(dataset)
summary(dataset)

# togliamo eventuali NA
newdata <- na.omit(dataset) 
dim(newdata)

# confrontando la dimensione del dataset e di newdata osserviamo che non ci sono dati mancanti;
# osserviamo però che ci sono dei valori nulli in colonne dove non ha senso ci siano: 
# per esempio la pressione del sangue non può essere zero. Di conseguenza consideriamo il valore nullo
# come un dato mancante e teniamo soltanto i valori significativi (quelli >0).
newdata = newdata[which(newdata$BMI>0),]
newdata = newdata[which(newdata$SkinThickness>0),]
newdata = newdata[which(newdata$BloodPressure>0),]
newdata = newdata[which(newdata$Glucose>0),]

# notiamo inoltre che se togliessimo tutti i valori nulli alla colonna insulina rimarrebbero la metà delle 
# osservazioni; decidiamo quindi di non considerare tale colonna 
newdata = newdata[,-5]
dim(newdata)
# restano circa il 70% delle osservazioni di partenza 

# estraiamo sample casuale
indici<-sample(1:nrow(newdata),0.9*nrow(newdata))
testdata<-newdata[-indici,]
newdata<-newdata[indici,]

attach(newdata)

# visualizziamo i dati senza la variabile Outcome 
x11()
ggcorr(newdata[ , c('BMI', 'Pregnancies','Glucose','SkinThickness','BloodPressure', 'DiabetesPedigreeFunction','Age')], label=TRUE)
x11()
ggpairs(newdata[ , c('BMI','Pregnancies','Glucose','SkinThickness','BloodPressure', 'DiabetesPedigreeFunction','Age')], pch = 16)

# procediamo con la costruzione del modello lineare, come y prendiamo BMI
g1=lm(BMI ~ Pregnancies + Glucose + SkinThickness + BloodPressure + DiabetesPedigreeFunction + Age)
summary(g1)
# R-adj^2=45%


#IC al 95% per pregnancy
n = dim(newdata)[1] 
p=g1$rank
alpha = 0.05
t_alpha2 = qt( 1-alpha/2, n-p )
beta_hat_glucose = g1$coefficients[2]
se_beta_hat_glucose = summary( g1 )[[4]][2,2]
IC_glucose = c( beta_hat_glucose - t_alpha2 * se_beta_hat_glucose,
                beta_hat_glucose + t_alpha2 * se_beta_hat_glucose )
IC_glucose

#IC al 95% per SkinThickness
beta_hat_skinthickness= g1$coefficients[4]
se_beta_hat_skinthickness = summary( g1 )[[4]][4,2]
IC_skinthickness = c( beta_hat_skinthickness - t_alpha2 * se_beta_hat_skinthickness,beta_hat_skinthickness + t_alpha2 * se_beta_hat_skinthickness )
IC_skinthickness

#ellisse con skinthickness e glucosio
x11()
plot( ellipse( g1, c( 2, 4 ) ), type = "l", xlim = c( -0.4,0.1),ylim = c( -0.001,0.45),main='Ellipse' )
points( 0, 0 )
points( g1$coef[ 2 ] , g1$coef[ 4 ] , pch = 18 )
abline( h = c( IC_skinthickness[1], IC_skinthickness[2] ), lty = 2 )
abline( v = c( IC_glucose[1], IC_glucose[2] ), lty = 2 )
# lo zero non è all'interno dell'ellisse, quindi vuol dire che almeno una covariata è significativa

# pregnancies e age non sono significative; applichiamo il metodo di step
step( g1, direction = "backward" , trace = T)
step( g1, direction = "backward", k = log(478), trace = T)
# dal modello iniziale tutti e due i metodi toglierebbero solo pregnacies

# decidiamo di fare un modello parsimonioso ma significativo: come covariate teniamo solo 
# SkinThickness e BloodPressure
g1=lm(formula = BMI ~ SkinThickness + BloodPressure)
summary(g1)
# R_adj^2 si è un po' abbassato (45%); ora tutte le covariate sono significative 
vif(g1)
AIC(g1)
BIC(g1)
# diagnostica dei residui:
# omoschedasticità
x11()
plot( g1$fit, g1$res, xlab = "Fitted", ylab = "Residuals", 
      main = "Residuals vs Fitted Values", pch = 16 )
abline( h = 0, lwd = 2, lty = 2, col = 'red' )  
plot(g1, which=1)
# residui non bruttissimi

# normalità
x11()
qqnorm( g1$res, ylab = "Raw Residuals", pch = 16 )
qqline( g1$res )
# code un po' discostate
shapiro.test( g1$res ) 
# p-value molto basso (dell'ordine di 10^-8)
x11()
hist( g1$res, 10, probability = TRUE, col = 'lavender', main = 'Hist residuals'  )
x11()
boxplot( g1$res, main = "Boxplot of residuals", pch = 16, col = 'lavender' )
# boxplot simmetrico ma si vedono alcuni outliers

# analisi punti influenti: 
# punti di leva
X = model.matrix( g1 )
lev = hat( X )
p = g1$rank
p
# p=3
n = dim(newdata)[1] 
n
# n=478

# plot dei punti di leva 
x11()
plot( g1$fitted.values, lev, ylab = "Leverages", main = "Plot of Leverages", 
      pch = 16, col = 'black' )
abline( h = 2 * p/n, lty = 2, col = 'red' ) 
watchout_points_lev = lev[ which( lev > 2 * p/n  ) ]
watchout_ids_lev = seq_along( lev )[ which( lev > 2 * p/n ) ] ## identify the rows relative to leverage points
points( g1$fitted.values[ watchout_ids_lev ], watchout_points_lev, col = 'red', pch = 16 )

# residui standardizzati 
gs = summary(g1)
res_std = g1$res/gs$sigma
watchout_ids_rstd = which( abs( res_std ) > 2 )
watchout_rstd = res_std[ watchout_ids_rstd ]


# plot dei residui standardizzati
x11()
plot( g1$fitted.values, res_std, ylab = "Standardized Residuals", main = "Standardized Residuals" )
abline( h = c(-2,2), lty = 2, col = 'orange' )
points( g1$fitted.values[watchout_ids_rstd], 
        res_std[watchout_ids_rstd], col = 'red', pch = 16 )
points( g1$fitted.values[watchout_ids_lev], 
        res_std[watchout_ids_lev], col = 'orange', pch = 16 )
legend('topright', col = c('red','orange'), 
       c('Standardized Residuals', 'Leverages'), pch = rep( 16, 2 ), bty = 'n' )

# residui studentizzati 
stud = rstandard( g1 )
watchout_ids_stud = which( abs( stud ) > 2 )
watchout_stud = stud[ watchout_ids_stud ]

# plot dei residui studentizzati
x11()
plot( g1$fitted.values, stud, ylab = "Studentized Residuals", main = "Studentized Residuals", pch = 16 )
points( g1$fitted.values[watchout_ids_stud], 
        stud[watchout_ids_stud], col = 'pink', pch = 16 )
points( g1$fitted.values[watchout_ids_lev], 
        stud[watchout_ids_lev], col = 'orange', pch = 16 )
abline( h = c(-2,2), lty = 2, col = 'orange' )
legend('topright', col = c('pink','orange'), 
       c('Studentized Residual', 'Leverages'), pch = rep( 16, 3 ), bty = 'n' )

# distanza di cook 
Cdist = cooks.distance( g1 )
watchout_ids_Cdist = which( Cdist > 4/(n-p) ) 
watchout_Cdist = Cdist[ watchout_ids_Cdist ]

# plot dei punti la cui distanza di cook supera la regola base (in verde)
x11()
plot( g1$fitted.values, Cdist, pch = 16, xlab = 'Fitted values', 
      ylab = 'Cooks Distance', main = 'Cooks Distance' )
points( g1$fitted.values[ watchout_ids_Cdist ], Cdist[ watchout_ids_Cdist ], 
        col = 'green', pch = 16 )

# togliamo sia i punti di leva sia i punti a distanza di cook elevata dal modello 
id_to_keep = !( 1:n %in% watchout_ids_Cdist )
g1=lm(formula = BMI ~  SkinThickness + BloodPressure, subset = ( lev < 2*p/n), newdata[ id_to_keep, ] )
summary(g1)
vif(g1)
AIC(g1)
BIC(g1)
# le covariate solo comunque significative, l'R_edj^2 aumenta a circa il 50%

# diagnostica dei residui: 
# omoschedasticità
x11()
plot( g1$fit, g1$res, xlab = "Fitted", ylab = "Residuals", 
      main = "Residuals vs Fitted Values", pch = 16 )
abline( h = 0, lwd = 2, lty = 2, col = 'red' )  
plot(g1, which=1)
# residui non bruttissimi

# normalità
x11()
qqnorm( g1$res, ylab = "Raw Residuals", pch = 16 )
qqline( g1$res )
# le code sono un po' migliorate 
shapiro.test( g1$res )
# il p-value è salito a 0.20, molto più alto di prima e ben sopra le soglie classiche (0.1, 0.05)
x11()
hist( g1$res, 10, probability = TRUE, col = 'lavender', main = 'Hist residuals'  )
# istogramma molto migliore di prima 
x11()
boxplot( g1$res, main = "Boxplot of residuals", pch = 16, col = 'lavender' )
# boxplot simmetrico, ora si vedono solo due puntini sopra 1.5*IQR

# previsione
# controllo il range delle covariate 
range(newdata$SkinThickness)
range(testdata$SkinThickness)
# ok va bene 
range(newdata$BloodPressure)
range(testdata$BloodPressure)
# ok va bene 

y.pred = predict( g1, testdata, interval = "confidence", se = T )
# valori predetti 
y = y.pred$fit[ ,1 ] 
# lower bound dell'intervallo di confidenza
y.pred$fit[ ,2 ]
# upper bound dell'intervallo di confidenza
y.pred$fit[ ,3 ] 
# errori in valore assoluto
err = abs(y-testdata$BMI) 
max(err)          
# errore in valore assoluto medio
MAD = sum(err)/dim(testdata)[1]
MAD
# mean squared error
MSE=sum((y-testdata$BMI)^2)/dim(testdata)[1] 
MSE
# root mean squared error
RMSE = sqrt(MSE) 
RMSE
# Percentuale di y vere che si trovano negli intervalli di confidenza
tabella = as.numeric(testdata$BMI>y.pred$fit[,2])*as.numeric(testdata$BMI<y.pred$fit[,3])
perc= sum(tabella)/dim(testdata)[1]*100 
perc


# regressione logistica: analisi con la variabile Outcome
# visualizziamo i dati 
x11()
par( mfrow = c( 2, 4 ) )
plot( newdata$Pregnancies, newdata$Outcome, pch = ifelse( newdata$Outcome == 1, 3, 4 ),
      col = ifelse( newdata$Outcome == 1, 'blue', 'deeppink' ),
      xlab = 'Pregnancies', ylab = 'Outcome', main = 'Outcome vs. Pregnancies', lwd = 2, cex = 1.5 )
x = seq(min(newdata$Pregnancies), max(newdata$Pregnancies), length=9)
mid = c((x[2:9]+x[1:8])/2)
mid
GR = cut(newdata$Pregnancies , breaks = x, include.lowest = TRUE, right = FALSE )
y = tapply (newdata$Outcome, GR, mean)
y
points( mid, y, col = "blue", pch = 16 )


plot( newdata$Glucose, newdata$Outcome, pch = ifelse( newdata$Outcome == 1, 3, 4 ),
      col = ifelse( newdata$Outcome == 1,'blue', 'deeppink'),
      xlab = 'Glucose', ylab = 'Outcome', main = 'Outcome vs. Glucose', lwd = 2, cex = 1.5 )
x = seq(min(newdata$Glucose), max(newdata$Glucose), length=9)
mid = c((x[2:9]+x[1:8])/2)
mid
GR = cut(newdata$Glucose , breaks = x, include.lowest = TRUE, right = FALSE )
y = tapply (newdata$Outcome, GR, mean)
y
points( mid, y, col = "blue", pch = 16 )



plot( newdata$BloodPressure, newdata$Outcome, pch = ifelse( newdata$Outcome == 1, 3, 4 ),
      col = ifelse( newdata$Outcome == 1, 'blue', 'deeppink' ),
      xlab = 'BloodPressure', ylab = 'Outcome', main = 'Outcome vs. Blood Pressure', lwd = 2, cex = 1.5 )
x = seq(min(newdata$BloodPressure), max(newdata$BloodPressure), length=9)
mid = c((x[2:9]+x[1:8])/2)
mid
GR = cut(newdata$BloodPressure , breaks = x, include.lowest = TRUE, right = FALSE )
y = tapply (newdata$Outcome, GR, mean)
y
points( mid, y, col = "blue", pch = 16 )



plot( newdata$SkinThickness, newdata$Outcome, pch = ifelse( newdata$Outcome == 1, 3, 4 ),
      col = ifelse( newdata$Outcome == 1, 'blue', 'deeppink' ),
      xlab = 'SkinThickness', ylab = 'Outcome', main = 'Outcome vs. Skin Thickness', lwd = 2, cex = 1.5 )
x = seq(min(newdata$SkinThickness), max(newdata$SkinThickness), length=9)
mid = c((x[2:9]+x[1:8])/2)
mid
GR = cut(newdata$SkinThickness , breaks = x, include.lowest = TRUE, right = FALSE )
y = tapply (newdata$Outcome, GR, mean)
y
points( mid, y, col = "blue", pch = 16 )


plot( newdata$BMI, newdata$Outcome, pch = ifelse( newdata$Outcome == 1, 3, 4 ),
      col = ifelse( newdata$Outcome == 1, 'blue', 'deeppink' ),
      xlab = 'BMI', ylab = 'Outcome', main = 'Outcome vs. BMI', lwd = 2, cex = 1.5 )
x = seq(min(newdata$BMI), max(newdata$BMI), length=9)
mid = c((x[2:9]+x[1:8])/2)
mid
GR = cut(newdata$BMI , breaks = x, include.lowest = TRUE, right = FALSE )
y = tapply (newdata$Outcome, GR, mean)
y
points( mid, y, col = "blue", pch = 16 )


plot( newdata$DiabetesPedigreeFunction, newdata$Outcome, pch = ifelse( newdata$Outcome == 1, 3, 4 ),
      col = ifelse( newdata$Outcome == 1, 'blue', 'deeppink' ),
      xlab = 'Diabetes Pedigree Function', ylab = 'Outcome', main = 'Outcome vs. Diabetes Pedigree Function', lwd = 2, cex = 1.5 )
x = seq(min(newdata$DiabetesPedigreeFunction), max(newdata$DiabetesPedigreeFunction), length=9)
mid = c((x[2:9]+x[1:8])/2)
mid
GR = cut(newdata$DiabetesPedigreeFunction , breaks = x, include.lowest = TRUE, right = FALSE )
y = tapply (newdata$Outcome, GR, mean)
y
points( mid, y, col = "blue", pch = 16 )


plot( newdata$Age, newdata$Outcome, pch = ifelse( newdata$Outcome == 1, 3, 4 ),
      col = ifelse( newdata$Outcome == 1, 'blue', 'deeppink' ),
      xlab = 'Age', ylab = 'Outcome', main = 'Outcome vs. Age', lwd = 2, cex = 1.5 )
x = seq(min(newdata$Age), max(newdata$Age), length=9)
mid = c((x[2:9]+x[1:8])/2)
mid
GR = cut(newdata$Age , breaks = x, include.lowest = TRUE, right = FALSE )
y = tapply (newdata$Outcome, GR, mean)
y
points( mid, y, col = "blue", pch = 16 )

#percentuale di uni
sum(Outcome)/478
# creiamo il modello di regressione logistica 
g=glm(Outcome ~ Pregnancies+Glucose+BloodPressure+SkinThickness+BMI+DiabetesPedigreeFunction+Age, family=binomial(link=logit))
summary(g)
# BloodPressure, SkinThickness e age non sono significative 

# metodo di step
step( g, direction = "backward" , trace = T)
# toglierebbe BloodPressure e skinthickness

AIC(g)
# 437
BIC(g)
# 471

step( g, direction = "backward", k = log(478), trace = T) 
# terrebbe pregnancies, glucose e BMI e DiabetesPedigreeFunction

# modello più parsimonioso
g=glm(Outcome ~ Pregnancies + Glucose + BMI + DiabetesPedigreeFunction, family=binomial(link=logit))
summary(g)
# tutte le covariate sono significative; residual deviance=425, AIC=435

# predittori lineari
g$linear.predictors

# calcoliamo i valori stimati per la probabilità di avere il diabete 
exp(g$linear.predictors)/(1+exp(g$linear.predictors))

# visualizziamo i residui
x11()
binnedplot(expected, rstandard(g))
# fit buono del modello, i residui sono attorno allo zero 

# test di GOF
hoslem.test( g$y, fitted( g ), g = 6 ) 
# p-value 0.28 che è alto

# tabella di (mis) classificazione
soglia = 0.5
valori.reali  = Outcome
valori.predetti = as.numeric( g$fitted.values > soglia ) 
# 1 se > soglia, 0 se < = soglia
valori.predetti
tab = table( valori.reali, valori.predetti )
tab

# accurancy
# % di casi classificati correttamente:
round( sum( diag( tab ) ) / sum( tab ), 2 )
# 79%

# % di casi misclassificati (complementare a quello sopra)
round( ( tab [ 1, 2 ] + tab [ 2, 1 ] ) / sum( tab ), 2 )
#21%

# sensitivity (indice di maggiore interesse perchè misura quanto il modello è buono a classificare gli uni)
sensitivita =  tab [ 2, 2 ] /( tab [ 2, 1 ] + tab [ 2, 2 ] ) 
sensitivita
#60%

# specificity 
specificita = tab[ 1, 1 ] /( tab [ 1, 2 ] + tab [ 1, 1 ] )
specificita
#88%

# costuzione della curva ROC
fit2 = g$fitted
soglia_roc  = seq( 0, 1, length.out = 2e2 )
lens = length( soglia_roc )-1
ascissa_roc  = rep( NA, lens )
ordinata_roc = rep( NA, lens )
for ( k in 1 : lens )
{
  soglia = soglia_roc [ k ]
  
  classification = as.numeric( sapply( fit2, function( x ) ifelse( x < soglia, 0, 1 ) ) )
  
  ordinata_roc[ k ] = sum( classification[ which( Outcome == 1 ) ] == 1 ) /
    length( which( Outcome == 1 ) )
  
  ascissa_roc[ k ] = sum( classification[ which( Outcome == 0 ) ] == 1 ) /
    length( which( Outcome == 0 ) )
}

# plot curva ROC.
x11()
plot( ascissa_roc, ordinata_roc, type = "l", xlab = "1 - Specificity", ylab = "Sensitivity",
      main = "Curva ROC", lwd = 2, col = 'darkblue', ylim = c( 0, 1 ), xlim = c( 0, 1 ) )
abline( h = c( 0, 1 ), v = c( 0, 1 ), lwd = 1, lty = 2, col = 'red' )
abline( a = 0, b = 1, lty = 2, col = 'black' )
abline( v = 1 - specificita,  h = sensitivita, lty = 3, col = 'blue' )
points( 1 - specificita, sensitivita, pch = 4, lwd = 3, cex = 1.5, col = 'blue' )

# plot alternatio della curva ROC 
PRROC_obj <- roc.curve(scores.class0 = g$fitted.values, weights.class0=as.numeric(paste(Outcome)),
                       curve=TRUE)
x11()
plot(PRROC_obj)

# proviamo ora a cambiare soglia per vedere se riusciamo a migliorare le prestazioni del modello 

# tabella di (mis) classificazione
soglia = 0.3
valori.reali  = Outcome
valori.predetti = as.numeric( g$fitted.values > soglia ) 
# 1 se > soglia, 0 se < = soglia
valori.predetti
tab = table( valori.reali, valori.predetti )
tab

# accurancy
# % di casi classificati correttamente:
round( sum( diag( tab ) ) / sum( tab ), 2 )

# % di casi misclassificati (complementare a quello sopra)
round( ( tab [ 1, 2 ] + tab [ 2, 1 ] ) / sum( tab ), 2 )

# sensitivity
sensitivita =  tab [ 2, 2 ] /( tab [ 2, 1 ] + tab [ 2, 2 ] ) 
sensitivita

# specificity 
specificita = tab[ 1, 1 ] /( tab [ 1, 2 ] + tab [ 1, 1 ] )
specificita


# previsione
summary(testdata)
soglia = 0.3
fitted.values= predict(g, testdata, type = 'response')
valori.reali  = testdata$Outcome
valori.predetti = as.numeric( fitted.values > soglia )
# 1 se > soglia, 0 se < = soglia
valori.predetti
tab.test = table( valori.reali, valori.predetti )
tab.test

# accuracy
# % di casi classificati correttamente:
round( sum( diag( tab.test ) ) / sum( tab.test ), 2 )
# % di casi misclassificati:
round( ( tab.test [ 1, 2 ] + tab.test [ 2, 1 ] ) / sum( tab.test ), 2 )

# sensitivity
sensitivita =  tab.test [ 2, 2 ] /( tab.test [ 2, 1 ] + tab.test [ 2, 2 ] ) 
sensitivita

# specificity 
specificita = tab.test[ 1, 1 ] /( tab.test [ 1, 2 ] + tab.test [ 1, 1 ] )
specificita

detach(newdata)

