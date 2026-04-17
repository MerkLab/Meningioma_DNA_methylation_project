library(survminer)
library(dplyr)
library(survival)
library(gtsummary)
library(pec)
library(tidyr)
library(broom)
library(tibble)
library(survIDINRI)
library(purrr) 
library(stringr)
library(rms)
library(splines)
library(scales)




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

data_all <- tmerge(data1=T26_survival, data2=T26_survival, id=ID, tstop=futime)
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

fit <- coxph(
  Surv(time = tstart, time2 = tstop, event = status) ~ top1000_avg,
  data = data_all
)
summary(fit)
c_index_cluster_nonlinear = concordance(fit)$concordance
c_index_cluster_nonlinear

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
#binary
fit.all <- coxph(Surv(time=tstart, time2 = tstop, event = status)~cluster1000_new + grading2021_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_all)
summary(fit.all)

c_index_full_model_cluster = concordance(fit.all)$concordance
c_index_full_model_cluster

#continuous
fit.all.cont <- coxph(Surv(time=tstart, time2 = tstop, event = status)~top1000_avg + grading2021_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_all)
summary(fit.all.cont)

c_index_full_model_cluster_cont = concordance(fit.all.cont)$concordance
c_index_full_model_cluster_cont


#check performance with molecular group instead of methylation cluster
fit.group <- coxph(Surv(time=tstart, time2 = tstop, event = status)~MCconsensus + grading2021_new + EOR + adjuvantRT + strata(setting) + strata(gender), data=data_all)
summary(fit.group)

c_index_full_model_group = concordance(fit.group)$concordance
c_index_full_model_group


#make graph for c-indices 
c_indices = rbind(c_index_age,c_index_gender,
                  c_index_setting,c_index_grade,
                  c_index_EOR,c_index_invasion,
                  c_index_RT,c_index_H3K27,c_index_cluster,c_index_cluster_nonlinear,
                  c_index_full_model_group,
                  c_index_full_model_cluster, c_index_full_model_cluster_cont)
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



###use restricted cubic splines to investigate an association between average marker methylation (continous) and outcome

dd = datadist(T26_survival)
options(datadist = "dd")

#fit spline Cox model
cox_spline <- cph(
  Surv(futime, status) ~ rcs(top1000_avg, 3),
  data = T26_survival,
  x = TRUE,
  y = TRUE,
  surv = TRUE
)

#calculate Hazard ratio predictions
pred <- Predict(
  cox_spline,
  top1000_avg,
  fun = exp,
  ref.zero = TRUE   # HR = 1 at median methylation
)

#convert to dataframe
pred_df <- as.data.frame(pred)

#change condition assignment
T26_survival <- T26_survival %>%
  mutate(WHO_grade = factor(grading2021_new, levels = c(1,2,3), labels = c("WHO I", "WHO II", "WHO III")))

# Convert numeric variables to factors with unique labels
T26_survival <- T26_survival %>%
  mutate(
    cluster1000_new_f = factor(cluster1000_new, levels = c(1,2), labels = c("Cluster 1", "Cluster 2")),
    MCconsensus_f     = factor(MCconsensus),
    risk_f = factor(risk_new)
  )

str(T26_survival)

## --- Define automatic rug geometry ---
ymin   <- min(pred_df$lower, na.rm = TRUE)
ymax   <- max(pred_df$upper, na.rm = TRUE)
yrange <- ymax - ymin

buffer <- yrange * 0.06
rug_band <- yrange * 0.34
rug_height <- (rug_band - buffer) / 4
gap <- rug_height * 0.12

y1 <- ymin - rug_band + gap
y2 <- y1 + rug_height + gap
y3 <- y2 + rug_height + gap
y4 <- y3 + rug_height + gap

ggplot(pred_df, aes(x = top1000_avg, y = yhat)) +
  geom_ribbon(
    aes(ymin = pmax(lower, ymin + yrange * 0.02), ymax = upper),
    fill = "#9ecae1", alpha = 0.4
  ) +
  geom_line(color = "#084594", linewidth = 2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_segment(
    data = T26_survival,
    aes(x = top1000_avg, xend = top1000_avg,
        y = y1, yend = y1 + rug_height,
        color = cluster1000_new_f),
    inherit.aes = FALSE,
    linewidth = 1.5,
    alpha = 0.9
  ) +
  geom_segment(
    data = T26_survival,
    aes(x = top1000_avg, xend = top1000_avg,
        y = y2, yend = y2 + rug_height,
        color = risk_f),
    inherit.aes = FALSE,
    linewidth = 1.5,
    alpha = 0.9
  ) +
  geom_segment(
    data = T26_survival,
    aes(x = top1000_avg, xend = top1000_avg,
        y = y3, yend = y3 + rug_height,
        color = MCconsensus_f),
    inherit.aes = FALSE,
    linewidth = 1.5,
    alpha = 0.9
  ) +
  geom_segment(
    data = T26_survival,
    aes(x = top1000_avg, xend = top1000_avg,
        y = y4, yend = y4 + rug_height,
        color = WHO_grade),
    inherit.aes = FALSE,
    linewidth = 1.5,
    alpha = 0.9
  ) +
  coord_cartesian(ylim = c(ymin - rug_band, ymax), clip = "off") +
  scale_color_manual(values = c(
    "Cluster 1" = "cyan4",
    "Cluster 2" = "violetred4",
    "WHO I" = "#9EBCDA",
    "WHO II" = "#8C6BB1",
    "WHO III" = "#810F7C",
    "Merlinintact" = "royalblue2",
    "Immuneenriched" = "red3",
    "hypermetabolic" = "forestgreen",
    "proliferative" = "darkorange2",
    "low" = "dodgerblue2",
    "intermediate" = "purple1",
    "high" = "firebrick1"
  )) +
  labs(
    x = "DNA methylation level",
    y = "Hazard ratio for progression",
    title = "Continuous Association Between Methylation and Progression Risk"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.line = element_line(color = "black"),
    plot.margin = margin(10, 20, 30, 10)
  )












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
fit_cluster_spline = coxph(
  Surv(futime, status) ~ ns(top1000_avg, df = 3),
  data = T26_survival,
  x = TRUE, y = TRUE
)
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

brier_model_cluster_spline <- pec(
  object = list("Spline signature" = fit_cluster_spline),
  formula = Surv(futime, status) ~ 1,
  data = T26_survival,
  times = times,
  exact = FALSE
)
ibs_cluster_spline <- crps(brier_model_cluster_spline)
ibs_cluster_spline
#integrated is 0.166


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
Brier_cluster_spline <- brier_model_cluster_spline$AppErr$`Spline signature`


#simulate a random survival model for grade 2 cohort
simulate_random_model <- function(data, times) {
  n <- nrow(data)
  
  # random subject-specific risk
  lp <- rnorm(n, mean = 0, sd = 0.75)
  
  # simple baseline hazard chosen to give reasonable decay
  base_hazard <- 0.005
  
  # survival probabilities: S(t) = exp(-H0(t) * exp(lp))
  probs <- sapply(times, function(t) {
    exp(-base_hazard * t * exp(lp))
  })
  
  probs <- as.matrix(probs)
  colnames(probs) <- as.character(times)
  
  dummy_model <- list(probs = probs)
  class(dummy_model) <- "randomSurvModel"
  dummy_model
}

# S3 method required by pec
predictSurvProb.randomSurvModel <- function(object, newdata, times, ...) {
  t_index <- match(as.character(times), colnames(object$probs))
  object$probs[, t_index, drop = FALSE]
}

#dataset is used for number of patients and times
dat <- T26_survival

# Choose your time points
times <- seq(0, 214, by = 3)

#fit random model and extract Brier scores
set.seed(123)

dummy_model <- simulate_random_model(dat, times)

brier_random_model <- pec(
  object = list("Random" = dummy_model),
  formula = Surv(futime, status) ~ 1,
  data = dat,
  times = times,
  exact = FALSE
)

# Vector of time-specific Brier scores
Brier_random <- brier_random_model$AppErr$Random
Brier_random


#bring together
Brier_scores = cbind(Brier_grade, Brier_cluster, Brier_cluster_spline, Brier_group, Brier_random)
Brier_scores = as.data.frame(Brier_scores)
Brier_scores$times = seq(from = 3, by = 3, length.out = 72)
colnames(Brier_scores) = c("grade", "cluster","cluster_nonlinear", "group", "random", "time")
write.table(Brier_scores, file="Brier_scores_full_cohort.txt", sep = "\t")

###make plots for Brier scores
#full cohort
df_long <- Brier_scores %>%
  pivot_longer(
    cols = -time,
    names_to = "model",
    values_to = "value"
  )
str(df_long)
df_long$model = factor(df_long$model, levels = c("random","grade", "group", "cluster", "cluster_nonlinear"))

cairo_pdf(filename = "Brier_scores_full_cohort.pdf", width = 9, height = 5)
ggplot(df_long, aes(x = time, y = value, color = model)) +
  geom_line(linewidth = 1.8) +
  scale_y_log10(
    breaks = c(0.05, 0.1, 0.2, 0.3),
    labels = c("0.05", "0.1", "0.2", "0.3")
  ) +
  scale_color_manual(values = c(
    "grade" = "dodgerblue1",
    "cluster" = "magenta1",
    "cluster_nonlinear" = "violetred4",
    "group" = "darkorange",
    "random" = "grey50"
  )) +
  labs(
    x = "Time",
    y = "Brier score",
    color = "Model"
  ) +
  scale_x_continuous(breaks = c(0, 48, 96, 144, 192))+
  theme_classic()
dev.off()

#test difference
#continous vs group
wilcox.test(Brier_scores$cluster_nonlinear,
            Brier_scores$group,
            paired = TRUE)

#binary vs group
wilcox.test(Brier_scores$cluster,
            Brier_scores$group,
            paired = TRUE)





#same for only CNS WHO grade 2
T26_survival_grade2 = T26_survival[T26_survival$grading2021 == "2",]

# Time points at which to evaluate Brier scores
times <- seq(0, 214, by = 3)

T26_survival_grade2 = subset(T26_survival, grade_bin == "grade2")

fit_cluster_grade2 <- coxph(Surv(futime, status)~cluster1000_new, data=T26_survival_grade2, x = TRUE, y = TRUE)
fit_cluster_spline_grade2 = coxph(
  Surv(futime, status) ~ ns(top1000_avg, df = 3),
  data = T26_survival_grade2,
  x = TRUE, y = TRUE
)
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
#integrated is 0.168

brier_model_cluster_spline <- pec(
  object = list("Spline signature" = fit_cluster_spline_grade2),
  formula = Surv(futime, status) ~ 1,
  data = T26_survival_grade2,
  times = times,
  exact = FALSE
)
ibs_cluster_spline <- crps(brier_model_cluster_spline)
ibs_cluster_spline
#integrated is 0.166


brier_model_group <- pec::pec(
  object = list("Cox" = fit_group_grade2),
  formula = Surv(futime, status) ~ 1,
  data = T26_survival_grade2,
  times = times,
  exact = FALSE
)
ibs_group_group <- crps(brier_model_group)
ibs_group_group
#integrated is 0.170

# Extract Brier scores for Cox model

Brier_cluster_grade2<- brier_model_cluster$AppErr$Cox
Brier_cluster_spline_grade2 <- brier_model_cluster_spline$AppErr$`Spline signature`
Brier_group_grade2 <- brier_model_group$AppErr$Cox




#simulate a random survival model for grade 2 cohort
simulate_random_model <- function(data, times) {
  n <- nrow(data)
  
  # random subject-specific risk
  lp <- rnorm(n, mean = 0, sd = 0.75)
  
  # simple baseline hazard chosen to give reasonable decay
  base_hazard <- 0.005
  
  # survival probabilities: S(t) = exp(-H0(t) * exp(lp))
  probs <- sapply(times, function(t) {
    exp(-base_hazard * t * exp(lp))
  })
  
  probs <- as.matrix(probs)
  colnames(probs) <- as.character(times)
  
  dummy_model <- list(probs = probs)
  class(dummy_model) <- "randomSurvModel"
  dummy_model
}

# S3 method required by pec
predictSurvProb.randomSurvModel <- function(object, newdata, times, ...) {
  t_index <- match(as.character(times), colnames(object$probs))
  object$probs[, t_index, drop = FALSE]
}

#dataset is used for number of patients and times
dat <- T26_survival_grade2

# Choose your time points
times <- seq(0, 214, by = 3)

#fit random model and extract Brier scores
set.seed(123)

dummy_model <- simulate_random_model(dat, times)

brier_random_model <- pec(
  object = list("Random" = dummy_model),
  formula = Surv(futime, status) ~ 1,
  data = dat,
  times = times,
  exact = FALSE
)

# Vector of time-specific Brier scores
Brier_random_grade2 <- brier_random_model$AppErr$Random
Brier_random_grade2



#bring together
Brier_scores_grade2 = cbind(Brier_cluster_grade2, Brier_cluster_spline_grade2, Brier_group_grade2, Brier_random_grade2)
Brier_scores_grade2 = as.data.frame(Brier_scores_grade2)
Brier_scores_grade2$times = seq(from = 3, by = 3, length.out = 72)
colnames(Brier_scores_grade2) = c("cluster", "cluster_nonlinear", "group","random", "time")
write.table(Brier_scores_grade2, file="Brier_scores_WHO_grade2.txt", sep = "\t")

#make statistical test for comparison of Brier scores within WHO grade 2 cases
shapiro.test(Brier_scores$cluster_nonlinear - Brier_scores$group)

#not normal distributed, use 
wilcox.test(Brier_scores_grade2$cluster_nonlinear,
            Brier_scores_grade2$group,
            paired = TRUE)
wilcox.test(Brier_scores_grade2$cluster,
            Brier_scores_grade2$group,
            paired = TRUE)


diff <- Brier_cluster_grade2 - Brier_group_grade2
summary(diff)
median(diff)
mean(diff)
#group has higher Brier score, so better predictive performance for cluster


###make plots for Brier scores
#grade2
df_long <- Brier_scores_grade2 %>%
  pivot_longer(
    cols = -time,
    names_to = "model",
    values_to = "value"
  )
str(df_long)
df_long$model = factor(df_long$model, levels = c("random", "group", "cluster", "cluster_nonlinear"))

cairo_pdf(filename = "Brier_scores_grade2.pdf", width = 9, height = 5)
ggplot(df_long, aes(x = time, y = value, color = model)) +
  geom_line(linewidth = 1.8) +
  scale_y_log10(
    breaks = c(0.05, 0.1, 0.2, 0.3),
    labels = c("0.05", "0.1", "0.2", "0.3")
  ) +
  scale_color_manual(values = c(
    "grade" = "dodgerblue1",
    "cluster" = "magenta1",
    "cluster_nonlinear" = "violetred4",
    "group" = "darkorange",
    "random" = "grey50"
  )) +
  labs(
    x = "Time",
    y = "Brier score",
    color = "Model"
  ) +
  scale_x_continuous(breaks = c(0, 48, 96, 144, 192))+
  theme_classic()
dev.off()






###perform integrated discrimination improvement analyses (IDI) and net reclassification improvement (NRI)
###IDI is performed for intermediate-risk subsets defined by CNS WHO grade 2, integrated intermediate-risk, or hypermetabolic

T26_survival$mol_bin <- ifelse(T26_survival$MCconsensus %in% c("Merlinintact", "Immuneenriched", "proliferative"),
                               "well-defined","hypermetabolic")

T26_survival$grade_bin <- ifelse(T26_survival$grading2021_new %in% c("1", "3"),
                               "well-defined","grade2")

T26_survival$risk_bin <- ifelse(T26_survival$risk_new %in% c("low", "high"),
                                 "well-defined","intermediate")


#based on molecular groups
#full cohort
#outcome matrix: first column = time, second = event
indata <- as.matrix(T26_survival[, c("futime", "status")])

#baseline model design matrix
covs0 <- model.matrix(
  ~ MCconsensus + grading2021_new + EOR + adjuvantRT,
  data = T26_survival
)[, -1, drop = FALSE]   # remove intercept

#new model design matrix
covs1 <- model.matrix(
  ~ MCconsensus + cluster1000_new + grading2021_new + EOR + adjuvantRT,
  data = T26_survival
)[, -1, drop = FALSE]

#choose time horizon
t0 <- 36   # or whatever unit matches futime

set.seed(123)
#run IDI / continuous NRI
x <- IDI.INF(indata, covs0, covs1, t0, npert = 2000)

#print results
IDI.INF.OUT(x)


#hypermetabolic cases
T26_intermediate <- subset(T26_survival, mol_bin == "hypermetabolic")

indata_int <- as.matrix(T26_intermediate[, c("futime", "status")])

covs0_int <- model.matrix(
  ~ grading2021_new + EOR + adjuvantRT,
  data = T26_intermediate
)[, -1, drop = FALSE]

covs1_int <- model.matrix(
  ~ grading2021_new + EOR + adjuvantRT + cluster1000_new,
  data = T26_intermediate
)[, -1, drop = FALSE]

x_int <- IDI.INF(indata_int, covs0_int, covs1_int, t0, npert = 300)
IDI.INF.OUT(x_int)



#based on grade
#full cohort
#outcome matrix: first column = time, second = event
indata <- as.matrix(T26_survival[, c("futime", "status")])

#baseline model design matrix
covs0 <- model.matrix(
  ~ grading2021_new + EOR + adjuvantRT,
  data = T26_survival
)[, -1, drop = FALSE]   # remove intercept

#new model design matrix
covs1 <- model.matrix(
  ~ grading2021_new + cluster1000_new + EOR + adjuvantRT,
  data = T26_survival
)[, -1, drop = FALSE]

#choose time horizon
t0 <- 36   # or whatever unit matches futime

set.seed(123)
#run IDI / continuous NRI
x <- IDI.INF(indata, covs0, covs1, t0, npert = 2000)

#print results
IDI.INF.OUT(x)


#WHO CNS grade 2 cases
T26_intermediate <- subset(T26_survival, grade_bin == "grade2")

indata_int <- as.matrix(T26_intermediate[, c("futime", "status")])

covs0_int <- model.matrix(
  ~ EOR + adjuvantRT,
  data = T26_intermediate
)[, -1, drop = FALSE]

covs1_int <- model.matrix(
  ~ cluster1000_new + EOR + adjuvantRT,
  data = T26_intermediate
)[, -1, drop = FALSE]

set.seed(123)
x_int <- IDI.INF(indata_int, covs0_int, covs1_int, t0, npert = 2000)
IDI.INF.OUT(x_int)




#based on integrated risk
#full cohort
#outcome matrix: first column = time, second = event
indata <- as.matrix(T26_survival[, c("futime", "status")])

#baseline model design matrix
covs0 <- model.matrix(
  ~ risk_bin + grading2021_new + EOR + adjuvantRT,
  data = T26_survival
)[, -1, drop = FALSE]   # remove intercept

#new model design matrix
covs1 <- model.matrix(
  ~ cluster1000_new  +grading2021_new + risk_bin + EOR + adjuvantRT,
  data = T26_survival
)[, -1, drop = FALSE]

#choose time horizon
t0 <- 36   # or whatever unit matches futime

set.seed(123)
#run IDI / continuous NRI
x <- IDI.INF(indata, covs0, covs1, t0, npert = 2000)

#print results
IDI.INF.OUT(x)


#intermediate integrated risk
T26_intermediate <- subset(T26_survival, risk_bin == "intermediate")

indata_int <- as.matrix(T26_intermediate[, c("futime", "status")])

covs0_int <- model.matrix(
  ~ EOR + adjuvantRT + grading2021_new,
  data = T26_intermediate
)[, -1, drop = FALSE]

covs1_int <- model.matrix(
  ~ cluster1000_new + EOR + adjuvantRT + grading2021_new,
  data = T26_intermediate
)[, -1, drop = FALSE]

set.seed(123)
x_int <- IDI.INF(indata_int, covs0_int, covs1_int, t0, npert = 2000)
IDI.INF.OUT(x_int)



#visualize IDI and NRI findings
#grade
plot_df <- data.frame(
  group = c(
    "Full cohort", "Full cohort",
    "grade2 subgroup", "grade2 subgroup"
  ),
  metric = c("IDI", "NRI", "IDI", "NRI"),
  estimate = c(0.051, 0.224, 0.065, 0.315),
  lower = c(0.005, -0.009, -0.003, 0.038),
  upper = c(0.114, 0.46, 0.18, 0.49),
  pvalue = c(0.025, 0.062, 0.066, 0.042)
)

plot_df$label <- sprintf(
  "%.3f (%.3f to %.3f), p = %s",
  plot_df$estimate,
  plot_df$lower,
  plot_df$upper,
  ifelse(plot_df$pvalue < 0.001, "<0.001", sprintf("%.3f", plot_df$pvalue))
)

plot_df$row_label <- paste(plot_df$group, "—", plot_df$metric)

plot_df$row_label <- factor(
  plot_df$row_label,
  levels = rev(plot_df$row_label)
)

cairo_pdf(filename = "IDI_NRI_grade.pdf", width = 5, height = 8)
ggplot(plot_df, aes(x = estimate, y = row_label)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 2) +
  geom_text(aes(label = label), hjust = 0.5, size = 3.3) +
  labs(
    x = "Change in discrimination / reclassification",
    y = NULL
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(min(plot_df$lower) - 0.05, max(plot_df$upper) + 0.18))
dev.off()

cairo_pdf(filename = "IDI_NRI_grade_clean.pdf", width = 5, height = 8)
ggplot(plot_df, aes(x = estimate, y = row_label)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 4) +
  labs(
    x = "Change in discrimination / reclassification",
    y = NULL
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  coord_cartesian(xlim = c(min(plot_df$lower) - 0.05, max(plot_df$upper) + 0.18))
dev.off()



#molecular group
plot_df <- data.frame(
  group = c(
    "Full cohort", "Full cohort",
    "hypermetabolic subgroup", "hypermetabolic subgroup"
  ),
  metric = c("IDI", "NRI", "IDI", "NRI"),
  estimate = c(0.008, 0.119, 0.021, 0.161),
  lower = c(-0.002, -0.206, -0.005, -0.18),
  upper = c(0.034, 0.386, 0.133, 0.415),
  pvalue = c(0.18, 0.373, 0.233, 0.379)
)

plot_df$label <- sprintf(
  "%.3f (%.3f to %.3f), p = %s",
  plot_df$estimate,
  plot_df$lower,
  plot_df$upper,
  ifelse(plot_df$pvalue < 0.001, "<0.001", sprintf("%.3f", plot_df$pvalue))
)

plot_df$row_label <- paste(plot_df$group, "—", plot_df$metric)

plot_df$row_label <- factor(
  plot_df$row_label,
  levels = rev(plot_df$row_label)
)

cairo_pdf(filename = "IDI_NRI_group.pdf", width = 5, height = 8)
ggplot(plot_df, aes(x = estimate, y = row_label)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 2) +
  geom_text(aes(label = label), hjust = 0.5, size = 3.3) +
  labs(
    x = "Change in discrimination / reclassification",
    y = NULL
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(min(plot_df$lower) - 0.05, max(plot_df$upper) + 0.18))
dev.off()

cairo_pdf(filename = "IDI_NRI_group_clean.pdf", width = 5, height = 8)
ggplot(plot_df, aes(x = estimate, y = row_label)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 4) +
  labs(
    x = "Change in discrimination / reclassification",
    y = NULL
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  coord_cartesian(xlim = c(min(plot_df$lower) - 0.05, max(plot_df$upper) + 0.18))
dev.off()



#integrated risk score
plot_df <- data.frame(
  group = c(
    "Full cohort", "Full cohort",
    "intermediate subgroup", "intermediate subgroup"
  ),
  metric = c("IDI", "NRI", "IDI", "NRI"),
  estimate = c(0.027, 0.12, 0.005, 0.141),
  lower = c(-0.004, -0.087, -0.019, -0.259),
  upper = c(0.08, 0.365, 0.076, 0.353),
  pvalue = c(0.112, 0.15, 0.754, 0.439)
)

plot_df$label <- sprintf(
  "%.3f (%.3f to %.3f), p = %s",
  plot_df$estimate,
  plot_df$lower,
  plot_df$upper,
  ifelse(plot_df$pvalue < 0.001, "<0.001", sprintf("%.3f", plot_df$pvalue))
)

plot_df$row_label <- paste(plot_df$group, "—", plot_df$metric)

plot_df$row_label <- factor(
  plot_df$row_label,
  levels = rev(plot_df$row_label)
)

cairo_pdf(filename = "IDI_NRI_risk.pdf", width = 5, height = 8)
ggplot(plot_df, aes(x = estimate, y = row_label)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 2) +
  geom_text(aes(label = label), hjust = 0.5, size = 3.3) +
  labs(
    x = "Change in discrimination / reclassification",
    y = NULL
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(min(plot_df$lower) - 0.05, max(plot_df$upper) + 0.18))
dev.off()

cairo_pdf(filename = "IDI_NRI_risk_clean.pdf", width = 3, height = 8)
ggplot(plot_df, aes(x = estimate, y = row_label)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 4) +
  labs(
    x = "Change in discrimination / reclassification",
    y = NULL
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  coord_cartesian(xlim = c(min(plot_df$lower) - 0.05, max(plot_df$upper) + 0.18))
dev.off()





###check IDI/NRI with cluster signature as a nonlinear variable

#molecular group
# full cohort
indata <- as.matrix(T26_survival[, c("futime", "status")])

# Model A: mol groups + clinical covariates
covs0_mol <- model.matrix(
  ~ MCconsensus + grading2021_new + EOR + adjuvantRT,
  data = T26_survival
)[, -1, drop = FALSE]


# spline basis for continuous signature
sig_ns <- ns(T26_survival$top1000_avg, df = 3)
sig_ns <- as.matrix(sig_ns)
colnames(sig_ns) <- paste0("top1000_ns", seq_len(ncol(sig_ns)))

# baseline
covs_clin <- model.matrix(
  ~ MCconsensus + grading2021_new + EOR + adjuvantRT,
  data = T26_survival
)[, -1, drop = FALSE]

# combine baseline + spline signature
covs1_sig <- cbind(covs_clin, sig_ns)

t0 <- 36

set.seed(123)
x_compare <- IDI.INF(
  indata,
  covs0_mol,
  covs1_sig,
  t0,
  npert = 2000
)

# Print results
IDI.INF.OUT(x_compare)



#hypermetabolic cases
T26_intermediate <- subset(T26_survival, mol_bin == "hypermetabolic")
indata <- as.matrix(T26_intermediate[, c("futime", "status")])

# baseline model
covs0 <- model.matrix(
  ~ grading2021_new + EOR + adjuvantRT,
  data = T26_intermediate
)[, -1, drop = FALSE]

# spline basis for continuous signature
sig_ns <- ns(T26_intermediate$top1000_avg, df = 3)

# convert to matrix and name columns
sig_ns <- as.matrix(sig_ns)
colnames(sig_ns) <- paste0("top1000_ns", seq_len(ncol(sig_ns)))

# baseline 
base_covs <- model.matrix(
  ~ grading2021_new + EOR + adjuvantRT,
  data = T26_intermediate
)[, -1, drop = FALSE]

# new model = baseline + spline basis
covs1_spline <- cbind(base_covs, sig_ns)

t0 <- 36

set.seed(123)
x_spline <- IDI.INF(indata, covs0, covs1_spline, t0, npert = 2000)

IDI.INF.OUT(x_spline)




#grade
# full cohort
indata <- as.matrix(T26_survival[, c("futime", "status")])

# baseline model
covs0 <- model.matrix(
  ~ grading2021_new + EOR + adjuvantRT,
  data = T26_survival
)[, -1, drop = FALSE]

# spline basis for continuous signature
sig_ns <- ns(T26_survival$top1000_avg, df = 3)

# convert to matrix and name columns
sig_ns <- as.matrix(sig_ns)
colnames(sig_ns) <- paste0("top1000_ns", seq_len(ncol(sig_ns)))

# baseline covariates as matrix
base_covs <- model.matrix(
  ~ grading2021_new + EOR + adjuvantRT,
  data = T26_survival
)[, -1, drop = FALSE]

# new model = baseline + spline basis
covs1_spline <- cbind(base_covs, sig_ns)

t0 <- 36

set.seed(123)
x_spline <- IDI.INF(indata, covs0, covs1_spline, t0, npert = 2000)

IDI.INF.OUT(x_spline)


#grade2 cases
T26_intermediate <- subset(T26_survival, grade_bin == "grade2")
indata <- as.matrix(T26_intermediate[, c("futime", "status")])

# baseline model
covs0 <- model.matrix(
  ~ EOR + adjuvantRT,
  data = T26_intermediate
)[, -1, drop = FALSE]

# spline basis for continuous signature
sig_ns <- ns(T26_intermediate$top1000_avg, df = 3)

# convert to matrix and name columns
sig_ns <- as.matrix(sig_ns)
colnames(sig_ns) <- paste0("top1000_ns", seq_len(ncol(sig_ns)))

# baseline covariates as matrix
base_covs <- model.matrix(
  ~ EOR + adjuvantRT,
  data = T26_intermediate
)[, -1, drop = FALSE]

# new model = baseline + spline basis
covs1_spline <- cbind(base_covs, sig_ns)

t0 <- 36

set.seed(123)
x_spline <- IDI.INF(indata, covs0, covs1_spline, t0, npert = 2000)

IDI.INF.OUT(x_spline)




#visualize IDI and NRI findings
#grade
plot_df <- data.frame(
  group = c(
    "Full cohort", "Full cohort",
    "grade2 subgroup", "grade2 subgroup"
  ),
  metric = c("IDI", "NRI", "IDI", "NRI"),
  estimate = c(0.074, 0.338, 0.147, 0.415),
  lower = c(0.022, 0.08, 0.053, 0.185),
  upper = c(0.141, 0.505, 0.273, 0.609),
  pvalue = c(0.006, 0.013, 0.001, 0.003)
)

plot_df$label <- sprintf(
  "%.3f (%.3f to %.3f), p = %s",
  plot_df$estimate,
  plot_df$lower,
  plot_df$upper,
  ifelse(plot_df$pvalue < 0.001, "<0.001", sprintf("%.3f", plot_df$pvalue))
)

plot_df$row_label <- paste(plot_df$group, "—", plot_df$metric)

plot_df$row_label <- factor(
  plot_df$row_label,
  levels = rev(plot_df$row_label)
)

cairo_pdf(filename = "IDI_NRI_grade_nonlinear.pdf", width = 5, height = 8)
ggplot(plot_df, aes(x = estimate, y = row_label)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 2) +
  geom_text(aes(label = label), hjust = 0.5, size = 3.3) +
  labs(
    x = "Change in discrimination / reclassification",
    y = NULL
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(min(plot_df$lower) - 0.05, max(plot_df$upper) + 0.18))
dev.off()

cairo_pdf(filename = "IDI_NRI_grade_nonlinear_clean.pdf", width = 5, height = 8)
ggplot(plot_df, aes(x = estimate, y = row_label)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 4) +
  labs(
    x = "Change in discrimination / reclassification",
    y = NULL
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  coord_cartesian(xlim = c(min(plot_df$lower) - 0.05, max(plot_df$upper) + 0.18))
dev.off()



#molecular group
plot_df <- data.frame(
  group = c(
    "Full cohort", "Full cohort",
    "hyper subgroup", "hyper subgroup"
  ),
  metric = c("IDI", "NRI", "IDI", "NRI"),
  estimate = c(0.035, 0.287, 0.092, 0.437),
  lower = c(0.004, 0.006, -0.024, -0.113),
  upper = c(0.078, 0.437, 0.321, 0.652),
  pvalue = c(0.025, 0.043, 0.1, 0.115)
)

plot_df$label <- sprintf(
  "%.3f (%.3f to %.3f), p = %s",
  plot_df$estimate,
  plot_df$lower,
  plot_df$upper,
  ifelse(plot_df$pvalue < 0.001, "<0.001", sprintf("%.3f", plot_df$pvalue))
)

plot_df$row_label <- paste(plot_df$group, "—", plot_df$metric)

plot_df$row_label <- factor(
  plot_df$row_label,
  levels = rev(plot_df$row_label)
)

cairo_pdf(filename = "IDI_NRI_group_nonlinear.pdf", width = 5, height = 8)
ggplot(plot_df, aes(x = estimate, y = row_label)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 2) +
  geom_text(aes(label = label), hjust = 0.5, size = 3.3) +
  labs(
    x = "Change in discrimination / reclassification",
    y = NULL
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(min(plot_df$lower) - 0.05, max(plot_df$upper) + 0.18))
dev.off()

cairo_pdf(filename = "IDI_NRI_group_nonlinear_clean.pdf", width = 5, height = 8)
ggplot(plot_df, aes(x = estimate, y = row_label)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 4) +
  labs(
    x = "Change in discrimination / reclassification",
    y = NULL
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  coord_cartesian(xlim = c(min(plot_df$lower) - 0.05, max(plot_df$upper) + 0.18))
dev.off()




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


####investigate mol_group stratification of 1p deleted cases
###only one proliferation case with intact chr1p, so excluded
#get data 
chr1p_group = read.csv("chr1p_group_data.csv", header = T)
chr1p_group$setting = factor(chr1p_group$setting, levels = c("Primary", "Recurrence"))
chr1p_group$cluster1000_new = factor(chr1p_group$cluster1000_new, levels = c("1", "2"))
chr1p_group$EOR = factor(chr1p_group$EOR, levels = c("GTR", "STR"))
chr1p_group$adjuvantRT = factor(chr1p_group$adjuvantRT, levels = c("no", "yes"))
chr1p_group$grading2021_new = factor(chr1p_group$grading2021_new, levels = c("1", "2","3"))
chr1p_group$chr1p_mol = factor(chr1p_group$chr1p_mol, 
                                     levels = c("intact_Merlinintact","deletion_Merlinintact","intact_Immuneenriched",
                                                "deletion_Immuneenriched","intact_hypermetabolic","deletion_hypermetabolic"))
str(chr1p_group)
chr1p_group$chr1p_mol


sfit_chr1p_group <- survfit(Surv(futime, status)~chr1p_mol, data=chr1p_group)
ggsurvplot(sfit_chr1p_group, risk.table=TRUE, palette=c("royalblue1","royalblue2", "red2", "red3","springgreen4","forestgreen"),
           conf.int=F, linetype = c("solid", "dashed","solid", "dashed", "solid", "dashed"),
           risk.table.height=.35, surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_chr1p_group, risk.table=FALSE, palette=c("royalblue1","royalblue2", "red2", "red3","springgreen4","forestgreen"),conf.int = F,
           linetype = c("solid", "dashed","solid", "dashed", "solid", "dashed"),
           size=1.8,xlim=c(1,200), break.x.by=48)

#pairwise comparisons

NF2 <- chr1p_group[chr1p_group$chr1p_mol %in% c("intact_Merlinintact", "deletion_Merlinintact"), ]
NF2$chr1p_mol = factor(NF2$chr1p_mol, levels = c("intact_Merlinintact", "deletion_Merlinintact")) 
Immune <- chr1p_group[chr1p_group$chr1p_mol %in% c("intact_Immuneenriched", "deletion_Immuneenriched"), ]
Immune$chr1p_mol = factor(Immune$chr1p_mol, levels = c("intact_Immuneenriched", "deletion_Immuneenriched")) 
Hyper <- chr1p_group[chr1p_group$chr1p_mol %in% c("intact_hypermetabolic", "deletion_hypermetabolic"), ]
Hyper$chr1p_mol = factor(Hyper$chr1p_mol, levels = c("intact_hypermetabolic", "deletion_hypermetabolic")) 


data_1p_NF2 <- tmerge(data1=NF2, data2=NF2, id=ID, tstop=futime)
data_1p_NF2 <- tmerge(data_1p_NF2, NF2, id=ID, adjuvantRT = tdc(adRTtime))

data_1p_Immune <- tmerge(data1=Immune, data2=Immune, id=ID, tstop=futime)
data_1p_Immune <- tmerge(data_1p_Immune, Immune, id=ID, adjuvantRT = tdc(adRTtime))

data_1p_Hyper <- tmerge(data1=Hyper, data2=Hyper, id=ID, tstop=futime)
data_1p_Hyper <- tmerge(data_1p_Hyper, Hyper, id=ID, adjuvantRT = tdc(adRTtime))

fit_NF2 <- coxph(Surv(time=tstart, time2 = tstop, event = status)~chr1p_mol +grading2021_new + EOR +cluster1000_new +strata(setting) + adjuvantRT + strata(gender), data=data_1p_NF2)
summary(fit_NF2)

fit_Immune <- coxph(Surv(time=tstart, time2 = tstop, event = status)~chr1p_mol +grading2021_new + EOR +cluster1000_new +strata(setting) + adjuvantRT + strata(gender), data=data_1p_Immune)
summary(fit_Immune)

fit_Hyper <- coxph(Surv(time=tstart, time2 = tstop, event = status)~chr1p_mol +grading2021_new + EOR +cluster1000_new +strata(setting) + adjuvantRT + strata(gender), data=data_1p_Hyper)
summary(fit_Hyper)

###comparisons within 1p-deleted cases with respect to mol_group
chr1p_all = read.csv("chr1p_group_data_all.csv", header = T)
chr1p_all$chr1p_mol = factor(chr1p_all$chr1p_mol, 
                               levels = c("deletion_Merlinintact","intact_Merlinintact","intact_Immuneenriched",
                                          "deletion_Immuneenriched","intact_hypermetabolic","deletion_hypermetabolic", "deletion_proliferative"))
str(chr1p_all)

sfit_chr1p_all <- survfit(Surv(futime, status)~chr1p_mol, data=chr1p_all)
ggsurvplot(sfit_chr1p_all, risk.table=TRUE, palette=c("royalblue1","royalblue2", "red2", "red3","springgreen4","forestgreen", "darkorange"),
           conf.int=F, linetype = c("dashed", "solid","solid", "dashed", "solid", "dashed", "dashed"),
           risk.table.height=.35, surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_chr1p_all, risk.table=FALSE, palette=c("royalblue1","royalblue2", "red2", "red3","springgreen4","forestgreen", "darkorange"),conf.int = F,
           linetype = c("dashed", "solid","solid", "dashed", "solid", "dashed", "dashed"),
           size=1.8,xlim=c(1,200), break.x.by=48)


data_1p_all <- tmerge(data1=chr1p_all, data2=chr1p_all, id=ID, tstop=futime)
data_1p_all <- tmerge(data_1p_all, chr1p_all, id=ID, adjuvantRT = tdc(adRTtime))

fit_all <- coxph(Surv(time=tstart, time2 = tstop, event = status)~chr1p_mol +grading2021_new + EOR +cluster1000_new +strata(setting) + adjuvantRT + strata(gender), data=data_1p_all)
summary(fit_all)

















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





