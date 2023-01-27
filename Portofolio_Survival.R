#Survival Package Exploration
#Laurensius Steven
#linkedin.com/in/laurensius-steven24/

library(survival) #survival library
library(KMsurv) #survival data library
library(stats)
library(flexsurv)
library(matlib)
library(dplyr)

data('larynx')#using larynx cancer data within KMsurv
attach(larynx)
View(larynx) #view data
help(larynx) #description

#Estimating Survival Curves
#kaplan-meier estimate
ts <- Surv(time, delta)
fs <- survfit(ts~1, data = larynx)
summary(fs)
plot(fs, main = 'Survival Curve')

#cumulative hazard curve
plot(fs$time, fs$cumhaz, main = 'Cumulative Hazard Curve')

#survival curve for each stadium
fs_stad <- survfit(ts~stage, data = larynx)
plot(fs_stad, main = 'Survival Curve per Stadium', 
     col = c('red', 'yellow', 'green', 'blue'))
legend(9.5, 1, legend = c('Stad 1', 'Stad 2', 'Stad 3', 'Stad 4'), lty = 1,
       col = c('red', 'yellow', 'green', 'blue'), cex = 0.5)


#survival curve fitted to a parametric curve
#exponential
fs_exp <- flexsurvreg(ts~1, data = larynx, dist = 'exponential')
summary(fs_exp)

#weibull
fs_wei <- flexsurvreg(ts~1, data = larynx, dist = 'weibull')
summary(fs_wei)

#gamma
fs_gam <- flexsurvreg(ts~1, data = larynx, dist = 'gamma')
summary(fs_gam)

plot(fs, main = 'Parametric Survival Curve', conf.int = FALSE)
lines(fs_exp, col = 'red', lwd = 0.5)
lines(fs_wei, col = 'green', lwd = 0.5)
lines(fs_gam, col = 'blue', lwd = 0.5)
legend(9, 1, 
       legend = c('Kaplan Meier', 'Exponential', 'Weibul', 'Gamma'), lty = 1,
       col = c('black', 'red', 'green', 'blue'), cex = 0.5)


#Cox Proportional Hazard Model
mod <- coxph(ts~ as.factor(stage) + age + diagyr, data = larynx)
summary(mod)

#interpretation:
#stadium 2 patients have 1.16 times the risk of death compared to stadium 1 patients, all things equal
#stadium 3 patients have 1.90 times the risk of death compared to stadium 1 patients, all things equal
#stadium 4 patients have 5.65 times the risk of death compared to stadium 1 patients, all things equal
#1 year older patient have 1.01 times the risk of death compared to others, all things equal
#earlier detection reduces the risk of death by a factor of 0.98, all things equal

#nested model testing
#Q: is age a useful variable to consider?
mod_a <- coxph(ts~1, data = larynx, ties = 'breslow')
mod_b <- coxph(ts~age, data = larynx, ties = 'breslow')

anova(mod_a, mod_b)
#p value: 0.10, not enough to break significance of 5%
#model containing age doesn't justify the added variable

#Q: is the stage of cancer a useful variable to consider?
mod_c <- coxph(ts~as.factor(stage), data = larynx, ties = 'breslow')

anova(mod_a, mod_c)
#p value: 0.001, more than enough to break significance of 5%
#model containing stage of cancer justifies the added complexity for added accuracy


#Difference of Survival Experience
#Q: is there a significance difference of the survival experience between 
#   patients of different stage of cancer?

survdiff(ts~as.factor(stage))
#p value: 5e-5, way lower than significance of 5%
#There is a significant difference between the stages of cancer

#trend test
#Q: is there an increasing risk of survival between patients of 
#   different stage of cancer?

zz <- c(1,2,3)
test.num <- zz %*% coef(mod_c)
test.var <- zz %*% mod_c$var %*% zz
ztren<-test.num/sqrt(test.var)

pnorm(ztren,lower.tail = FALSE)
#p value: 4e-4, way lower than significance of 5%
#There is a significant increasing risk between the stages of cancer

#stratified test

#Q: is the survival experience old between old aged and young aged patient
#   identical, across all stage of cancer?

df<-mutate(larynx, age_bin = ntile(age, n=2)) #binning age group
df<-df[-c(3,4)]
df1<-df[which(df$age_bin==1),]
df2<-df[which(df$age_bin==2),]

mod1<-coxph(Surv(time,delta)~factor(stage),data=df1)
mod2<-coxph(Surv(time,delta)~factor(stage),data=df2)

z<-mod1$coefficients + mod2$coefficients
sigma<-mod1$var + mod2$var
q<-z %*% inv(sigma) %*% z

pchisq(q, 1, lower.tail = FALSE)
#p value: 6e-5, enough significance 
#survival experience in some stages of cancer is different,
#patient from the same age group but different stage experience different risk of death

#Q: is the survival experience old between patients of different stage
#   identical, across all age group?

df1<-df[which(df$stage==1),]
df2<-df[which(df$stage==2),]
df3<-df[which(df$stage==3),]
df4<-df[which(df$stage==4),]

mod1<-coxph(Surv(time,delta)~age_bin,data=df1)
mod2<-coxph(Surv(time,delta)~age_bin,data=df2)
mod3<-coxph(Surv(time,delta)~age_bin,data=df3)
mod4<-coxph(Surv(time,delta)~age_bin,data=df4)

z<-mod1$coefficients + mod2$coefficients + mod3$coefficients + mod4$coefficients
sigma<-mod1$var + mod2$var + mod3$var + mod4$var
q<-z %*% 1/sigma %*% z

pchisq(q, 3, lower.tail = FALSE)
#p value: 0.904, not significant
#among patients of the same stage, the risk of death is identical
# regardless of the age group
