library(survminer)
library(dplyr)
library(survival)
library(gtsummary)
library(pec)
library(tidyr)
library(broom)
library(tibble)




#####start working on T26
###get T26 survival data with all predictors
getwd()
setwd("~/Desktop/Meningioma_2/several_data_revision")

T26_survival = read.csv(file = "meta_survival_current.csv", header = T)
T26_survival$gender = factor(T26_survival$gender, levels = c("M", "F"))
T26_survival$setting = factor(T26_survival$setting, levels = c("Primary", "Recurrence"))
T26_survival$grading2021_new = factor(T26_survival$grading2021_new, levels = c("1", "2", "3"))
T26_survival$risk_new = factor(T26_survival$risk_new, levels = c("low", "intermediate",
                                                         "high"))
T26_survival$invasionhisto = factor(T26_survival$invasionhisto, levels = c("no", "yes"))
T26_survival$MCconsensus = factor(T26_survival$MCconsensus, levels = c("Merlinintact", "Immuneenriched",
                                                                       "hypermetabolic","proliferative"))
T26_survival$cluster1000_new = factor(T26_survival$cluster1000_new, levels = c("1", "2"))
T26_survival$EOR = factor(T26_survival$EOR, levels = c("GTR", "STR"))
T26_survival$H3K27 = factor(T26_survival$H3K27, levels = c("wildtype", "mutated"))
T26_survival$priorRT = factor(T26_survival$priorRT, levels = c("no", "yes"))
T26_survival$adjuvantRT = factor(T26_survival$adjuvantRT, levels = c("no", "yes"))
str(T26_survival)


####adjuvant RT is not at baseline, use tmerge to generate a counting process format

data_all <- tmerge(data1=T26_survival[, 1:28], data2=T26_survival, id=ID, tstop=futime)
data_all <- tmerge(data_all, T26_survival, id=ID, adjuvantRT = tdc(adRTtime))
str(data_all)


#first check assumption of Cox proportional hazards (i.e. Schoenfeld residuals are independent from time)
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new + setting + grading2021_new + EOR +gender + adjuvantRT + age + H3K27 +invasionhisto, data=data_all)
fit

test.ph = cox.zph(fit)
test.ph
ggcoxzph = ggcoxzph(test.ph, font.main=7)
ggcoxzph[1]



####univariate analysis of all covariates
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~age, data=data_all)
summary(fit)
c_index_age = concordance(fit)$concordance
c_index_age

fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~setting, data=data_all)
summary(fit)
c_index_setting = concordance(fit)$concordance
c_index_setting

fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~grading2021_new, data=data_all)
summary(fit)
c_index_grade = concordance(fit)$concordance
c_index_grade


fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~EOR, data=data_all)
summary(fit)
c_index_EOR = concordance(fit)$concordance
c_index_EOR

fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~invasionhisto, data=data_all)
summary(fit)
c_index_invasion = concordance(fit)$concordance
c_index_invasion


str(T26_survival)
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~H3K27, data=data_all)
summary(fit)
c_index_H3K27 = concordance(fit)$concordance
c_index_H3K27

str(T26_survival)
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new, data=data_all)
summary(fit)
c_index_cluster = concordance(fit)$concordance
c_index_cluster


fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~adjuvantRT, data=data_all)
summary(fit)
c_index_RT = concordance(fit)$concordance
c_index_RT


fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~gender, data=data_all)
summary(fit)
c_index_gender = concordance(fit)$concordance
c_index_gender





######perform multivariate analysis for the entire cohort with potentially relevant covariates
#####sex and tumor setting are accommodated using stratification

fit.all <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new + grading2021_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_all)
summary(fit.all)

c_index_full_model_cluster = concordance(fit.all)$concordance
c_index_full_model_cluster


#check performance with molecular group instead of methylation cluster
fit.group <- coxph(Surv(time=tstart, time2 = tstop, event = status)~MCconsensus + grading2021_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_all)
summary(fit.group)

c_index_full_model_group = concordance(fit.group)$concordance
c_index_full_model_group


#make graph for c-indices 
c_indices = rbind(c_index_age,c_index_gender,
                  c_index_setting,c_index_grade,
                  c_index_EOR,c_index_invasion,
                  c_index_RT,c_index_H3K27,c_index_cluster,
                  c_index_full_model_group,
                  c_index_full_model_cluster)
c_indices = as.data.frame(c_indices)
c_indices$model = row.names(c_indices)



ggplot(c_indices, aes(x=V1, y=reorder(model,-V1))) + 
  geom_bar(stat = "identity", color="dodgerblue", fill="dodgerblue")+
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "red", size=1.5)+
  theme_classic()





####check Kaplan Meier curves for selected covariates

#cluster all
str(T26_survival)
sfit_cluster <- survfit(Surv(futime, status)~cluster1000_new, data=T26_survival)
ggsurvplot(sfit_cluster, risk.table=TRUE, palette=c("cyan4","violetred4"), 
           risk.table.height=.35,conf.int=TRUE,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_cluster, risk.table=FALSE, palette=c("cyan4","violetred4"),conf.int=TRUE,surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48, size=1.5)




#check survival of all METHlow, separated by molecular groups, and same for METHhigh, excluding conditions with less than 5 survival points
#do multivariate analysis accounting for potential confounding factors

#cluster 1
T26_survival_clust1 = read.csv(file = "T26_survival_cluster1.csv", header = T)
str(T26_survival_clust1)
T26_survival_clust1$MCconsensus = factor(T26_survival_clust1$MCconsensus, levels = c("Merlinintact", "Immuneenriched",
                                                                                     "hypermetabolic", "proliferative"))
T26_survival_clust1$grading2021_new = factor(T26_survival_clust1$grading2021_new, levels = c("1", "2", "3"))
T26_survival_clust1$EOR = factor(T26_survival_clust1$EOR, levels = c("GTR", "STR"))
T26_survival_clust1$adjuvantRT = factor(T26_survival_clust1$adjuvantRT, levels = c("no", "yes"))

sfit_consensus_clust1 <- survfit(Surv(futime, status)~MCconsensus, data=T26_survival_clust1)
ggsurvplot(sfit_consensus_clust1, risk.table=TRUE, palette=c("royalblue2","red3","forestgreen","darkorange2"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_consensus_clust1, risk.table=FALSE, palette=c("royalblue2","red3","forestgreen","darkorange2"), 
           xlim=c(1,200), break.x.by=48, size=1.5, conf.int = TRUE,conf.int.alpha=.2,surv.median.line = c("hv"),)

conf.int.alpha=.2


data_cluster1 <- tmerge(data1=T26_survival_clust1, data2=T26_survival_clust1, id=ID, tstop=futime)
data_cluster1 <- tmerge(data_cluster1, T26_survival_clust1, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~MCconsensus +grading2021_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_cluster1)
summary(fit)




#cluster 2
T26_survival_clust2 = read.csv(file = "T26_survival_cluster2.csv", header = T)
T26_survival_clust2$MCconsensus = factor(T26_survival_clust2$MCconsensus, levels = c("hypermetabolic","proliferative"))
T26_survival_clust2$grading2021_new = factor(T26_survival_clust2$grading2021_new, levels = c("2", "3"))
T26_survival_clust2$EOR = factor(T26_survival_clust2$EOR, levels = c("GTR", "STR"))
T26_survival_clust2$adjuvantRT = factor(T26_survival_clust2$adjuvantRT, levels = c("no", "yes"))

sfit_consensus_clust2 <- survfit(Surv(futime, status)~MCconsensus, data=T26_survival_clust2)
ggsurvplot(sfit_consensus_clust2, risk.table=TRUE, palette=c("forestgreen","darkorange2"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_consensus_clust2, risk.table=FALSE, palette=c("forestgreen","darkorange2"),xlim=c(1,200),
           surv.median.line = c("hv"),break.x.by=48, size=1.5, conf.int = TRUE,conf.int.alpha=.2)


data_cluster2 <- tmerge(data1=T26_survival_clust2, data2=T26_survival_clust2, id=ID, tstop=futime)
data_cluster2 <- tmerge(data_cluster2, T26_survival_clust2, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~MCconsensus +grading2021_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_cluster2)
summary(fit)


#show only Proliferative meningiomas stratified by cluster for reviewing purposes
T26_proliferative = T26_survival[T26_survival$MCconsensus == "proliferative",]
T26_proliferative$MCconsensus = factor(T26_proliferative$MCconsensus, levels = c("1", "2"))


sfit_proliferative <- survfit(Surv(futime, status)~cluster1000_new, data=T26_proliferative)

ggsurvplot(sfit_proliferative, risk.table=TRUE, palette=c("forestgreen","darkorange2"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)

ggsurvplot(sfit_proliferative, risk.table=FALSE, palette=c("darkorange2","darkorange2"),xlim=c(1,200),
           linetype = c("dashed", "solid"),break.x.by=48, size=1.5, risk.tabel=TRUE)


save.image()


#split grades by cluster and test

#grade1
T26_survival_grade1 = read.csv(file = "T26_survival_grade1.csv", header = T)
T26_survival_grade1$cluster1000_new = factor(T26_survival_grade1$cluster1000_new, levels = c("1","2"))
T26_survival_grade1$EOR = factor(T26_survival_grade1$EOR, levels = c("GTR", "STR"))
T26_survival_grade1$adjuvantRT = factor(T26_survival_grade1$adjuvantRT, levels = c("no", "yes"))

sfit_grade1 <- survfit(Surv(futime, status)~cluster1000_new, data=T26_survival_grade1)
ggsurvplot(sfit_grade1, risk.table=TRUE, palette=c("cyan4","violetred4"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_grade1, risk.table=FALSE, palette=c("cyan4","violetred4"),
           surv.median.line = c("hv"), xlim=c(1,200), break.x.by=48, size=1.5)

#correct for confounding factors
data_grade1 <- tmerge(data1=T26_survival_grade1, data2=T26_survival_grade1, id=ID, tstop=futime)
data_grade1 <- tmerge(data_grade1, T26_survival_grade1, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_grade1)
summary(fit)


#grade2
T26_survival_grade2 = read.csv(file = "T26_survival_grade2.csv", header = T)
T26_survival_grade2$cluster1000_new = factor(T26_survival_grade2$cluster1000_new, levels = c("1","2"))
T26_survival_grade2$EOR = factor(T26_survival_grade2$EOR, levels = c("GTR", "STR"))
T26_survival_grade2$adjuvantRT = factor(T26_survival_grade2$adjuvantRT, levels = c("no", "yes"))

sfit_grade2 <- survfit(Surv(futime, status)~cluster1000_new, data=T26_survival_grade2)
ggsurvplot(sfit_grade2, risk.table=TRUE, palette=c("cyan4","violetred4"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_grade2, risk.table=FALSE, palette=c("cyan4","violetred4"),surv.median.line = c("hv"), xlim=c(1,200), break.x.by=48, size=1.5)

#correct for confounding factors
data_grade2 <- tmerge(data1=T26_survival_grade2, data2=T26_survival_grade2, id=ID, tstop=futime)
data_grade2 <- tmerge(data_grade2, T26_survival_grade2, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_grade2)
summary(fit)



#grade3
T26_survival_grade3 = read.csv(file = "T26_survival_grade3.csv", header = T)
T26_survival_grade3$cluster1000_new = factor(T26_survival_grade3$cluster1000_new, levels = c("1","2"))
T26_survival_grade3$EOR = factor(T26_survival_grade3$EOR, levels = c("GTR", "STR"))
T26_survival_grade3$adjuvantRT = factor(T26_survival_grade3$adjuvantRT, levels = c("no", "yes"))

sfit_grade3 <- survfit(Surv(futime, status)~cluster1000_new, data=T26_survival_grade3)
ggsurvplot(sfit_grade3, risk.table=TRUE, palette=c("cyan4","violetred4"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_grade3, risk.table=FALSE, palette=c("cyan4","violetred4"),surv.median.line = c("hv"), xlim=c(1,200), break.x.by=48, size=1.5)

#correct for confounding factors
data_grade3 <- tmerge(data1=T26_survival_grade3, data2=T26_survival_grade3, id=ID, tstop=futime)
data_grade3 <- tmerge(data_grade3, T26_survival_grade2, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_grade3)
summary(fit)



#check superior separation of patients by cluster over grade, MC groups, and risk score

#get data for grade2
grade2_cluster = read.csv("T26_survival_grade2.csv", header = T)
str(grade2_cluster)
grade2_cluster$grade_cluster = factor(grade2_cluster$grade_cluster, levels = c("grade2_1", "grade2_2"))
grade2_cluster$setting = factor(grade2_cluster$setting, levels = c("Primary", "Recurrence"))
grade2_cluster$EOR = factor(grade2_cluster$EOR, levels = c("GTR", "STR"))
grade2_cluster$adjuvantRT = factor(grade2_cluster$adjuvantRT, levels = c("no", "yes"))


sfit_grade_cluster <- survfit(Surv(futime, status)~grade_cluster, data=grade2_cluster)
ggsurvplot(sfit_grade_cluster, risk.table=TRUE, palette=c("cyan4","violetred4"),conf.int=TRUE, 
           risk.table.height=.35, linetype = c("dashed", "solid"), surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_grade_cluster, risk.table=FALSE, palette=c("cyan4","violetred4"),conf.int=TRUE,
           size = 1.5,surv.median.line = c("hv"),conf.int.alpha=.2,xlim=c(1,200), break.x.by=48)

#correct for confounding factors
data_grade2 <- tmerge(data1=grade2_cluster, data2=grade2_cluster, id=ID, tstop=futime)
data_grade2 <- tmerge(data_grade2, grade2_cluster, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_grade2)
summary(fit)


#get data for hypermetabolic
hypermetabolic_cluster = read.csv("T26_survival_hypermetabolic_cluster.csv", header = T)
str(hypermetabolic_cluster)
hypermetabolic_cluster$consensus_cluster = factor(hypermetabolic_cluster$consensus_cluster, 
                                                  levels = c("hypermetabolic_1", "hypermetabolic_2"))
hypermetabolic_cluster$setting = factor(hypermetabolic_cluster$setting, levels = c("Primary", "Recurrence"))
hypermetabolic_cluster$EOR = factor(hypermetabolic_cluster$EOR, levels = c("GTR", "STR"))
hypermetabolic_cluster$adjuvantRT = factor(hypermetabolic_cluster$adjuvantRT, levels = c("no", "yes"))
hypermetabolic_cluster$grading2021_new = factor(hypermetabolic_cluster$grading2021_new, levels = c("1", "2","3"))


sfit_hypermetabolic_cluster <- survfit(Surv(futime, status)~consensus_cluster, data=hypermetabolic_cluster)
ggsurvplot(sfit_hypermetabolic_cluster, risk.table=TRUE, palette=c("cyan4","violetred4"),conf.int=TRUE, 
           risk.table.height=.35, linetype = c("dashed", "solid"), surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_hypermetabolic_cluster, risk.table=FALSE, palette=c("cyan4","violetred4"),conf.int=TRUE,
           size=1.5,surv.median.line = c("hv"),conf.int.alpha=.2,xlim=c(1,200), break.x.by=48)

#correct for confounding factors
data_hyper <- tmerge(data1=hypermetabolic_cluster, data2=hypermetabolic_cluster, id=ID, tstop=futime)
data_hyper <- tmerge(data_hyper, hypermetabolic_cluster, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new +grading2021_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_hyper)
summary(fit)



#get data for risk

risk_cluster = read.csv("T26_survival_intermediate_cluster.csv", header = T)
str(risk_cluster)
risk_cluster$risk_cluster = factor(risk_cluster$risk_cluster, 
                                   levels = c("intermediate_1", "intermediate_2"))
risk_cluster$setting = factor(risk_cluster$setting, levels = c("Primary", "Recurrence"))
risk_cluster$EOR = factor(risk_cluster$EOR, levels = c("GTR", "STR"))
risk_cluster$adjuvantRT = factor(risk_cluster$adjuvantRT, levels = c("no", "yes"))
risk_cluster$grading2021_new = factor(risk_cluster$grading2021_new, levels = c("1", "2","3"))

sfit_risk_cluster <- survfit(Surv(futime, status)~risk_cluster, data=risk_cluster)
ggsurvplot(sfit_risk_cluster, risk.table=TRUE, palette=c("cyan4","violetred4"),conf.int=TRUE, 
           risk.table.height=.35, linetype = c("dashed", "solid"), surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_risk_cluster, risk.table=FALSE, palette=c("cyan4","violetred4"),conf.int=TRUE,
           size=1.5,surv.median.line = c("hv"),conf.int.alpha=.2,xlim=c(1,200), break.x.by=48)

#correct for confounding factors
data_intermediate <- tmerge(data1=risk_cluster, data2=risk_cluster, id=ID, tstop=futime)
data_intermediate <- tmerge(data_intermediate, risk_cluster, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new +grading2021_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_intermediate)
summary(fit)




###calculate Brier scores
#total cohort
# Time points at which to evaluate Brier scores
times <- seq(0, 214, by = 3)

fit_grade <- coxph(Surv(futime, status)~grading2021_new, data=T26_survival, x = TRUE, y = TRUE)
fit_cluster <- coxph(Surv(futime, status)~cluster1000_new, data=T26_survival, x = TRUE, y = TRUE)
fit_group <- coxph(Surv(futime,status)~MCconsensus, data=T26_survival, x = TRUE, y = TRUE)

# Compute Brier score for Cox model

brier_model_grade <- pec::pec(
  object = list("Cox" = fit_grade),
  formula = Surv(futime, status) ~ 1,
  data = T26_survival,
  times = times,
  exact = FALSE
)
ibs_grade <- crps(brier_model_grade)
ibs_grade
#integrated is 0.186

brier_model_cluster <- pec::pec(
  object = list("Cox" = fit_cluster),
  formula = Surv(futime, status) ~ 1,
  data = T26_survival,
  times = times,
  exact = FALSE
)
ibs_cluster <- crps(brier_model_cluster)
ibs_cluster
#integrated is 0.178

brier_model_group <- pec::pec(
  object = list("Cox" = fit_group),
  formula = Surv(futime, status) ~ 1,
  data = T26_survival,
  times = times,
  exact = FALSE
)
ibs_group <- crps(brier_model_group)
ibs_group
#integrated is 0.174


# Extract Brier scores for Cox model
Brier_grade <- brier_model_grade$AppErr$Cox
Brier_cluster <- brier_model_cluster$AppErr$Cox
Brier_group <- brier_model_group$AppErr$Cox


# Define S3 class and method for dummy survival model
simulate_random_model <- function(data, times) {
  n <- nrow(data)
  m <- length(times)
  
  # Generate random survival probabilities between 0 and 1
  probs <- matrix(runif(n * m), nrow = n, ncol = m)
  colnames(probs) <- as.character(times)
  
  dummy_model <- list(
    predictSurvProb = function(object, newdata, times) {
      t_index <- match(as.character(times), colnames(probs))
      probs[, t_index, drop = FALSE]
    }
  )
  
  class(dummy_model) <- "randomSurvModel"
  return(dummy_model)
}

# Register predictSurvProb method for custom class
predictSurvProb.randomSurvModel <- function(object, newdata, times, ...) {
  object$predictSurvProb(object, newdata, times)
}

lung = survival::lung
lung <- na.omit(lung)
lung$status <- lung$status - 1

# Run simulations
set.seed(123)
n_sim <- 1000
brier_random <- matrix(NA, nrow = n_sim, ncol = length(times))

for (i in 1:n_sim) {
  dummy_model <- simulate_random_model(T26_survival, times)
  
  sim_result <- tryCatch({
    b <- pec::pec(
      object = list("Random" = dummy_model),
      formula = Surv(time, status) ~ 1,
      data = lung,
      times = times,
      exact = FALSE
    )
    b$AppErr$Random
  }, error = function(e) {
    warning(paste("Simulation", i, "failed:", e$message))
    rep(NA, length(times))
  })
  
  # Check that result is assignable
  if (length(sim_result) != length(times)) {
    warning(paste("Simulation", i, "returned wrong length; setting NA"))
    sim_result <- rep(NA, length(times))
  }
  
  brier_random[i, ] <- sim_result
}
sim_result

#bring together
Brier_scores = cbind(Brier_grade, Brier_cluster, Brier_group, sim_result)
Brier_scores = as.data.frame(Brier_scores)
Brier_scores$times = seq(from = 3, by = 3, length.out = 72)
colnames(Brier_scores) = c("grade", "cluster", "group", "random", "time")
write.table(Brier_scores, file="Brier_scores_full_cohort.txt", sep = "\t")



#same for only CNS WHO grade 2
T26_survival_grade2 = T26_survival[T26_survival$grading2021 == "2",]

# Time points at which to evaluate Brier scores
times <- seq(0, 214, by = 3)

fit_cluster_grade2 <- coxph(Surv(futime, status)~cluster1000_new, data=T26_survival_grade2, x = TRUE, y = TRUE)
fit_group_grade2 <- coxph(Surv(futime,status)~MCconsensus, data=T26_survival_grade2, x = TRUE, y = TRUE)

# Compute Brier score for Cox model
brier_model_cluster <- pec::pec(
  object = list("Cox" = fit_cluster_grade2),
  formula = Surv(futime, status) ~ 1,
  data = T26_survival_grade2,
  times = times,
  exact = FALSE
)
ibs_cluster_grade <- crps(brier_model_cluster)
ibs_cluster_grade
#integrated is 0.16

brier_model_group <- pec::pec(
  object = list("Cox" = fit_group_grade2),
  formula = Surv(futime, status) ~ 1,
  data = T26_survival_grade2,
  times = times,
  exact = FALSE
)
ibs_group_grade <- crps(brier_model_group)
ibs_group_grade
#integrated is 0.163

# Extract Brier scores for Cox model

Brier_cluster_grade2<- brier_model_cluster$AppErr$Cox
Brier_group_grade2 <- brier_model_group$AppErr$Cox


# Define S3 class and method for dummy survival model
simulate_random_model <- function(data, times) {
  n <- nrow(data)
  m <- length(times)
  
  # Generate random survival probabilities between 0 and 1
  probs <- matrix(runif(n * m), nrow = n, ncol = m)
  colnames(probs) <- as.character(times)
  
  dummy_model <- list(
    predictSurvProb = function(object, newdata, times) {
      t_index <- match(as.character(times), colnames(probs))
      probs[, t_index, drop = FALSE]
    }
  )
  
  class(dummy_model) <- "randomSurvModel"
  return(dummy_model)
}

# Register predictSurvProb method for custom class
predictSurvProb.randomSurvModel <- function(object, newdata, times, ...) {
  object$predictSurvProb(object, newdata, times)
}

# Run simulations
set.seed(123)
n_sim <- 1000
brier_random <- matrix(NA, nrow = n_sim, ncol = length(times))

for (i in 1:n_sim) {
  dummy_model <- simulate_random_model(T26_survival_grade2, times)
  
  sim_result <- tryCatch({
    b <- pec::pec(
      object = list("Random" = dummy_model),
      formula = Surv(time, status) ~ 1,
      data = lung,
      times = times,
      exact = FALSE
    )
    b$AppErr$Random
  }, error = function(e) {
    warning(paste("Simulation", i, "failed:", e$message))
    rep(NA, length(times))
  })
  
  # Check that result is assignable
  if (length(sim_result) != length(times)) {
    warning(paste("Simulation", i, "returned wrong length; setting NA"))
    sim_result <- rep(NA, length(times))
  }
  
  brier_random[i, ] <- sim_result
}
sim_result

#bring together
Brier_scores_grade2 = cbind(Brier_cluster_grade2, Brier_group_grade2, sim_result)
Brier_scores_grade2 = as.data.frame(Brier_scores_grade2)
Brier_scores_grade2$times = seq(from = 3, by = 3, length.out = 72)
colnames(Brier_scores_grade2) = c("cluster", "group","random", "time")
write.table(Brier_scores_grade2, file="Brier_scores_WHO_grade2.txt", sep = "\t")






####1p deletions are significantly associated with genome instability and CNVs
###check effect of 1p deletions on progression in METHlow and METHhigh (1p intact in MEHThigh exlcuded, only 2 events)
#get data
chr1p_cluster = read.csv("chr1p_cluster_data.csv", header = T)
str(chr1p_cluster)
chr1p_cluster$setting = factor(chr1p_cluster$setting, levels = c("Primary", "Recurrence"))
chr1p_cluster$cluster1000_new = factor(chr1p_cluster$cluster1000_new, levels = c("1", "2"))
chr1p_cluster$EOR = factor(chr1p_cluster$EOR, levels = c("GTR", "STR"))
chr1p_cluster$adjuvantRT = factor(chr1p_cluster$adjuvantRT, levels = c("no", "yes"))
chr1p_cluster$grading2021_new = factor(chr1p_cluster$grading2021_new, levels = c("1", "2","3"))
chr1p_cluster$chr1p_cluster = factor(chr1p_cluster$chr1p_cluster, 
                                   levels = c("intact_1", "deletion_1","deletion_2"))

sfit_chr1p_cluster <- survfit(Surv(futime, status)~chr1p_cluster, data=chr1p_cluster)
ggsurvplot(sfit_chr1p_cluster, risk.table=TRUE, palette=c("#00D7D1","#008080","violetred4"),
           conf.int=TRUE, linetype = c("dashed", "solid","solid"),
           risk.table.height=.35, surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_chr1p_cluster, risk.table=FALSE, palette=c("#00D7D1","#008080","violetred4"),conf.int=TRUE,
           linetype = c("solid", "dashed","dashed"),
           size=1.8,surv.median.line = c("hv"),conf.int.alpha=.2,xlim=c(1,200), break.x.by=48)

#correct for confounding factors
#with reference chr1p_intact
data_1p <- tmerge(data1=chr1p_cluster, data2=chr1p_cluster, id=ID, tstop=futime)
data_1p <- tmerge(data_1p, chr1p_cluster, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~chr1p_cluster + cluster1000_new +grading2021_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_1p)
summary(fit)


#comparison 1p deleted cases only
chr1p_deleted = chr1p_cluster[!chr1p_cluster$chr1p_cluster == "intact_1",]
chr1p_deleted$chr1p_cluster = factor(chr1p_deleted$chr1p_cluster, levels = c("deletion_1", "deletion_2")) 
str(chr1p_deleted)
data_1p_deleted <- tmerge(data1=chr1p_deleted, data2=chr1p_deleted, id=ID, tstop=futime)
data_1p_deleted <- tmerge(data_1p_deleted, chr1p_deleted, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~chr1p_cluster +grading2021_new + EOR +cluster1000_new +strata(setting) + adjuvantRT + strata(gender), data=data_1p_deleted)
summary(fit)


###check effect of genome instability and unfavorable CNVs on progression in the entire TUE cohort
#get data
T26_survival$CNV_intervals = factor(T26_survival$CNV_intervals, levels = c("zero", "up2",
                                                                           "more2"))
T26_survival$stability_intervals = factor(T26_survival$stability_intervals, levels = c("less1", "less10",
                                                                                       "less20", "more20"))

#check instability as a factor
sfit_T26_instability <- survfit(Surv(futime, status)~stability_intervals, data=T26_survival)
ggsurvplot(sfit_T26_instability, risk.table=TRUE, palette=c("grey35","indianred1","darkorchid3", "firebrick3"),
           conf.int=TRUE, 
           risk.table.height=.35, surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_T26_instability, risk.table=FALSE, palette=c("grey35","indianred1","darkorchid3", "firebrick3"),conf.int=TRUE,
           size=1.8,surv.median.line = c("hv"),conf.int.alpha=.2,xlim=c(1,200), break.x.by=48)

#correct for confounding factors
data_all <- tmerge(data1=T26_survival[, 1:28], data2=T26_survival, id=ID, tstop=futime)
data_all <- tmerge(data_all, T26_survival, id=ID, adjuvantRT = tdc(adRTtime))

fit.all <- coxph(Surv(time=tstart, time2 = tstop, event = status)~stability_intervals +cluster1000_new + grading2021_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_all)
summary(fit.all)


#check number of unfavorable CNVs as a factor
sfit_T26_CNVs <- survfit(Surv(futime, status)~CNV_intervals, data=T26_survival)
ggsurvplot(sfit_T26_CNVs, risk.table=TRUE, palette=c("grey35","indianred1", "firebrick3"),
           conf.int=TRUE, 
           risk.table.height=.35, surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_T26_CNVs, risk.table=FALSE, palette=c("grey35","indianred1", "firebrick3"),conf.int=TRUE,
           size=1.8,surv.median.line = c("hv"),conf.int.alpha=.2,xlim=c(1,200), break.x.by=48)

#correct for confounding factors
fit.all <- coxph(Surv(time=tstart, time2 = tstop, event = status)~CNV_intervals +cluster1000_new + grading2021_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_all)
summary(fit.all)




####make survival analyses for CNVs (DEL and AMP separately) for CNVs that are present in at least 10% of either METHlow or METHhigh
#modify CNV dataframe
data_CNV = read.csv(file="CNVs_DEL_AMP.csv", header = T)
str(data_CNV)

data_CNV <- data_CNV %>%
  mutate(across(-1, ~ factor(ifelse(. %in% c("DEL", "AMP"), "yes", "no"))))

write.csv(data_CNV, file="CNVs_as_yes_no.csv")


#get survival data
T26_survival_CNV = read.csv(file = "CNV_AMP_or_DEL_survival.csv", header = T)
str(T26_survival_CNV)
T26_survival_CNV$gender = factor(T26_survival_CNV$gender, levels = c("F", "M"))
T26_survival_CNV$setting = factor(T26_survival_CNV$setting, levels = c("Primary", "Recurrence"))
T26_survival_CNV$grading2021 = factor(T26_survival_CNV$grading2021, levels = c("1", "2", "3"))
T26_survival_CNV$EOR = factor(T26_survival_CNV$EOR, levels = c("GTR", "STR"))
T26_survival_CNV$RT = factor(T26_survival_CNV$RT, levels = c("no", "yes"))
T26_survival_CNV$cluster_2000 = factor(T26_survival_CNV$cluster_2000, levels = c("1", "2"))

library(broom)       
library(dplyr)
library(purrr) 

T26_survival_CNV <- T26_survival_CNV %>% 
  mutate(across(2:33, ~ factor(., levels = c("no", "yes"))))
str(T26_survival_CNV)

#start performing univariate analyses to decimate potential CNVs
#univariate analysis of all CNVs
fit <- coxph(Surv(time, status)~CNV_22q_loss, data=T26_survival_CNV)
summary(fit)


#get pval from univariate analysis and make padj

CNV_uni_pval = read.csv(file="CNV_univariate_for_padj.csv", header = F)
CNV_uni_padj = p.adjust(CNV_uni_pval$V1, method = "BH")

#make subset of sign CNVs
T26_survival_CNV_sub = T26_survival_CNV_sub[,-30]




#loop for each CNV of interest for multivariate analysis
cnv_vars   <- grep("^CNV", names(T26_survival_CNV), value = TRUE)   # all CNV columns start with "CNV"

covariates <- c( "grading2021","EOR", "RT", "cluster_2000", "gender")

stopifnot(all(c("time", "status", cnv_vars, covariates) %in% names(T26_survival_CNV)))


library(stringr)   # str_c(), str_starts()

fit_one_cnv <- function(cnv) {
  frm <- as.formula(
    str_c("Surv(time, status) ~ ", str_c(c(cnv, covariates), collapse = " + "))
  )
  
  fit <- coxph(frm, data = T26_survival_CNV)
  
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term == str_c(cnv, "yes")) |>
    mutate(CNV = cnv) |>
    select(CNV, HR = estimate, CI_lo = conf.low, CI_hi = conf.high, p_val = p.value)
}

#Apply over every CNV and row‑bind the outputs
results <- map_dfr(cnv_vars, fit_one_cnv)
results <- results %>%
  mutate(q_val = p.adjust(p_val, method = "BH"))
write.csv(results, file="Multivariate_per_CNV.csv")



#check PFS for 3 clusters solution of discovery cohort
T26_survival = read.csv(file = "meta_survival_current.csv", header = T)
T26_survival$gender = factor(T26_survival$gender, levels = c("M", "F"))
T26_survival$setting = factor(T26_survival$setting, levels = c("Primary", "Recurrence"))
T26_survival$grading2021_new = factor(T26_survival$grading2021_new, levels = c("1", "2", "3"))
T26_survival$EOR = factor(T26_survival$EOR, levels = c("GTR", "STR"))
T26_survival$adjuvantRT = factor(T26_survival$adjuvantRT, levels = c("no", "yes"))
T26_survival$cluster1000_new_3clusters= factor(T26_survival$cluster1000_new_3clusters, levels = c("METHlow-low", "METHlow-high", "METHhigh"))
str(T26_survival)

sfit_3clusters <- survfit(Surv(futime, status)~cluster1000_new_3clusters_2, data=T26_survival)
ggsurvplot(sfit_3clusters, risk.table=TRUE, palette=c("cyan2","cyan4","violetred4"), 
           risk.table.height=.35,conf.int=F,xlim=c(1,200), break.x.by=48,linetype = c("solid", "dashed","solid"),size=1.5,conf.int.alpha=.2)
ggsurvplot(sfit_3clusters, risk.table=FALSE, palette=c("cyan2","cyan4","violetred4"),conf.int=F,
           ,linetype = c("solid", "dashed","solid"),conf.int.alpha=.2,surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48,size=1.5)


#correct for confounding factors for METHlow
data_3clusters <- tmerge(data1=T26_survival, data2=T26_survival, id=ID, tstop=futime)
data_3clusters <- tmerge(data_3clusters, T26_survival, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new_3clusters_2 +adjuvantRT +grading2021_new + EOR + strata(setting) + strata(gender), data=data_3clusters)
summary(fit)




##separate by prior RT yes or no, entire cohort and clusters
#check entire cohort fo priorRT effect
T26_survival$priorRT = factor(T26_survival$priorRT, levels = c("no", "yes"))
sfit_prior <- survfit(Surv(futime, status)~priorRT, data=T26_survival)
ggsurvplot(sfit_prior, risk.table=TRUE, palette=c("blue","red"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48,size=1.8)
ggsurvplot(sfit_prior, risk.table=FALSE, palette=c("blue","red"),xlim=c(1,200), break.x.by=48,size=1.8)

#correct for confounding factors 
data_prior_all <- tmerge(data1=T26_survival, data2=T26_survival, id=ID, tstop=futime)
data_prior_all <- tmerge(data_prior_all, T26_survival, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~priorRT +cluster1000_new +adjuvantRT +grading2021_new + EOR + strata(setting) + strata(gender), data=data_prior_all)
summary(fit)




##split by cluster
#METHlow
T26_survival_clust1 = T26_survival[T26_survival$cluster1000_new == "1",]

sfit_prior_clust1 <- survfit(Surv(futime, status)~priorRT, data=T26_survival_clust1)
ggsurvplot(sfit_prior_clust1, risk.table=TRUE, palette=c("blue","red"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48,size=1.8)
ggsurvplot(sfit_prior_clust1, risk.table=FALSE, palette=c("blue","red"),surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48,size=1.8)

#correct for confounding factors for METHlow
data_prior_clust1 <- tmerge(data1=T26_survival_clust1, data2=T26_survival_clust1, id=ID, tstop=futime)
data_prior_clust1 <- tmerge(data_prior_clust1, T26_survival_clust1, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~priorRT +adjuvantRT +grading2021_new + EOR + strata(setting) + strata(gender), data=data_prior_clust1)
summary(fit)


#METHhigh
T26_survival_clust2 = T26_survival[T26_survival$cluster1000_new == "2",]

sfit_prior_clust2 <- survfit(Surv(futime, status)~priorRT, data=T26_survival_clust2)
ggsurvplot(sfit_prior_clust2, risk.table=TRUE, palette=c("blue","red"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48,size=1.8)
ggsurvplot(sfit_prior_clust2, risk.table=FALSE, palette=c("blue","red"),surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48,size=1.8)

#correct for confounding factors for METHlow
data_prior_clust2 <- tmerge(data1=T26_survival_clust2, data2=T26_survival_clust2, id=ID, tstop=futime)
data_prior_clust2 <- tmerge(data_prior_clust2, T26_survival_clust2, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~priorRT +adjuvantRT +grading2021_new + EOR + strata(setting) + strata(gender), data=data_prior_clust2)
summary(fit)


#since all priorRT cases are recurrences, also compare priorRT only in recurrent meningioma
#recurrent cases
T26_survival_recurrent = T26_survival[T26_survival$setting == "Recurrence",]

sfit_prior_rec <- survfit(Surv(futime, status)~priorRT, data=T26_survival_recurrent)
ggsurvplot(sfit_prior_rec, risk.table=TRUE, palette=c("blue","red"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48,size=1.8)
ggsurvplot(sfit_prior_rec, risk.table=FALSE, palette=c("blue","red"),surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48,size=1.8)

#correct for confounding factors
data_prior_rec <- tmerge(data1=T26_survival_recurrent, data2=T26_survival_recurrent, id=ID, tstop=futime)
data_prior_rec <- tmerge(data_prior_rec, T26_survival_recurrent, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~priorRT +adjuvantRT +grading2021_new + EOR + strata(gender), data=data_prior_rec)
summary(fit)










########next question: does adjuvant RT benefit METHhigh as compared to METHlow
#start with entire cohort

T26_survival = read.csv(file = "meta_survival_current.csv", header = T)
T26_survival$gender = factor(T26_survival$gender, levels = c("M", "F"))
T26_survival$setting = factor(T26_survival$setting, levels = c("Primary", "Recurrence"))
T26_survival$grading2021_new = factor(T26_survival$grading2021_new, levels = c("1", "2", "3"))
T26_survival$MCconsensus = factor(T26_survival$MCconsensus, levels = c("Merlinintact", "Immuneenriched",
                                                                       "hypermetabolic","proliferative"))
T26_survival$cluster1000_new = factor(T26_survival$cluster1000_new, levels = c("1", "2"))
T26_survival$EOR = factor(T26_survival$EOR, levels = c("GTR", "STR"))
T26_survival$adjuvantRT = factor(T26_survival$adjuvantRT, levels = c("no", "yes"))
T26_survival$adRT_cluster = factor(T26_survival$adRT_cluster, levels = c("no_1", "yes_1", "no_2", "yes_2"))

sfit_adjuvantRT_cluster <- survfit(Surv(futime, status)~adRT_cluster, data=T26_survival)
ggsurvplot(sfit_adjuvantRT_cluster, risk.table=TRUE, palette=c("cyan4","cyan2","violetred4", "violetred3"), 
           risk.table.height=.35,conf.int=TRUE,xlim=c(1,200), break.x.by=48,linetype = c("solid", "solid","dashed", "dashed"),size=1.5,conf.int.alpha=.2)
ggsurvplot(sfit_adjuvantRT_cluster, risk.table=FALSE, palette=c("cyan4","cyan2","violetred4", "violetred3"),conf.int=F,
           ,linetype = c("solid", "dashed","solid", "dashed"),conf.int.alpha=.2,surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48,size=1.5)


#correct for confounding factors for METHlow
T26_survival_clust1 = T26_survival[T26_survival$cluster1000_new == "1",]
str(T26_survival_clust1)
T26_survival_clust1$adjuvantRT = factor(T26_survival_clust1$adjuvantRT, levels = c("no", "yes"))

data_clusterRT_low <- tmerge(data1=T26_survival_clust1, data2=T26_survival_clust1, id=ID, tstop=futime)
data_clusterRT_low <- tmerge(data_clusterRT_low, T26_survival_clust1, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~adjuvantRT +grading2021_new + EOR + strata(setting) + strata(gender), data=data_clusterRT_low)
summary(fit)

#correct for confounding factors for METHhigh
T26_survival_clust2 = T26_survival[T26_survival$cluster1000_new == "2",]
str(T26_survival_clust1)
T26_survival_clust2$adjuvantRT = factor(T26_survival_clust2$adjuvantRT, levels = c("no", "yes"))

data_clusterRT_high <- tmerge(data1=T26_survival_clust2, data2=T26_survival_clust2, id=ID, tstop=futime)
data_clusterRT_high <- tmerge(data_clusterRT_high, T26_survival_clust2, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~adjuvantRT +grading2021_new + EOR + strata(setting) + strata(gender), data=data_clusterRT_high)
summary(fit)


####additionally, add an interaction term to test whether there is an interaction of cluster and adjuvantRT in the entire cohort
data_clusterRT_all <- tmerge(data1=T26_survival, data2=T26_survival, id=ID, tstop=futime)
data_clusterRT_all<- tmerge(data_clusterRT_all, T26_survival, id=ID, adjuvantRT = tdc(adRTtime))
fit_no_interaction <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new + adjuvantRT +grading2021_new +EOR + strata(setting) + strata(gender), data=data_clusterRT_all)
summary(fit_no_interaction)

fit_interaction <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new*adjuvantRT +grading2021_new +EOR + strata(setting) + strata(gender), data=data_clusterRT_all)
summary(fit_interaction)

anova(fit_no_interaction, fit_interaction, test = "LRT")
#while adjuvant RT only improves significantly in METHlow, interaction term analysis does not provide evidence for significant modulation of RT sensitivity in between clusters



#######split to clusters, with information on extend of resection and RT

#cluster 1
sfit_grade_clust1 <- survfit(Surv(futime, status)~EOR_adRT, data=T26_survival_clust1)
ggsurvplot(sfit_grade_clust1, risk.table=TRUE, palette=c("blue2","dodgerblue2","red2","deeppink2"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48,linetype = c("solid", "dashed","solid", "dashed"),size=1.8)
ggsurvplot(sfit_grade_clust1, risk.table=FALSE, palette=c("blue2","dodgerblue2","red2","deeppink2"),xlim=c(1,200), break.x.by=48,
           linetype = c("solid", "dashed","solid", "dashed"),size=1.8,surv.median.line = c("hv"))

#correct for confounding factors for cluster1/GTR
T26_survival_clust1_GTR = T26_survival_clust1[T26_survival_clust1$EOR == "GTR",]
str(T26_survival_clust1_GTR)
T26_survival_clust1_GTR$adjuvantRT = factor(T26_survival_clust1_GTR$adjuvantRT, levels = c("no", "yes"))

data_cluster1_GTR <- tmerge(data1=T26_survival_clust1_GTR, data2=T26_survival_clust1_GTR, id=ID, tstop=futime)
data_cluster1_GTR <- tmerge(data_cluster1_GTR, T26_survival_clust1_GTR, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~adjuvantRT +grading2021_new + strata(setting) + strata(gender), data=data_cluster1_GTR)
summary(fit)

#correct for confounding factors for cluster1/STR
T26_survival_clust1_STR = T26_survival_clust1[T26_survival_clust1$EOR == "STR",]
str(T26_survival_clust1_STR)
T26_survival_clust1_STR$adjuvantRT = factor(T26_survival_clust1_STR$adjuvantRT, levels = c("no", "yes"))

data_cluster1_STR <- tmerge(data1=T26_survival_clust1_STR, data2=T26_survival_clust1_STR, id=ID, tstop=futime)
data_cluster1_STR <- tmerge(data_cluster1_STR, T26_survival_clust1_STR, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~adjuvantRT +grading2021_new + strata(setting) + strata(gender), data=data_cluster1_STR)
summary(fit)


#cluster 2
sfit_grade_clust2 <- survfit(Surv(futime, status)~EOR_adRT, data=T26_survival_clust2)
ggsurvplot(sfit_grade_clust2, risk.table=TRUE, palette=c("blue2","red2","deeppink2"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48,linetype = c("solid", "dashed","solid", "dashed"),size=1.8)
ggsurvplot(sfit_grade_clust2, risk.table=FALSE, palette=c("blue2","red2","deeppink2"),xlim=c(1,200), break.x.by=48,
           linetype = c("solid","solid", "dashed"),size=1.8,surv.median.line = c("hv"))

#correct for confounding factors for cluster1/STR
T26_survival_clust2_STR = T26_survival_clust2[T26_survival_clust2$EOR == "STR",]
str(T26_survival_clust2_STR)
T26_survival_clust2_STR$adjuvantRT = factor(T26_survival_clust2_STR$adjuvantRT, levels = c("no", "yes"))

data_cluster2_STR <- tmerge(data1=T26_survival_clust2_STR, data2=T26_survival_clust2_STR, id=ID, tstop=futime)
data_cluster2_STR <- tmerge(data_cluster2_STR, T26_survival_clust2_STR, id=ID, adjuvantRT = tdc(adRTtime))
fit <- coxph(Surv(time=tstart, time2 = tstop, event = status)~adjuvantRT +grading2021_new + strata(setting) + strata(gender), data=data_cluster2_STR)
summary(fit)



###for text
#calculate median time to progression across entire cohort
km_fit <- survfit(Surv(futime, status) ~ 1, data = T26_survival)
summary(km_fit)$table["median"]
quantile(km_fit, probs = c(0.25, 0.5, 0.75))





save.image()





