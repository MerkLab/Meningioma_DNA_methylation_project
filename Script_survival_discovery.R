library(dplyr)
library(survival)
library(survminer)

#get survival data
T26_survival = read.csv(file = "T26_survival.csv", header = T)
str(T26_survival)
T26_survival$gender = factor(T26_survival$gender, levels = c("M", "F"))
T26_survival$RT = factor(T26_survival$RT, levels = c("no", "yes"))
T26_survival$setting = factor(T26_survival$setting, levels = c("Primary", "Recurrence"))
T26_survival$MCconsensus = factor(T26_survival$MCconsensus, levels = c("Merlinintact", "Immuneenriched",
                                                                       "hypermetabolic","proliferative"))
T26_survival$MCsubtype = factor(T26_survival$MCsubtype, levels = c("Benign", "Intermediate",
                                                                   "Malignant"))
T26_survival$grading2016 = factor(T26_survival$grading2016, levels = c("1", "2", "3"))
T26_survival$grading2021 = factor(T26_survival$grading2021, levels = c("1", "2", "3"))
T26_survival$cluster_2000 = factor(T26_survival$cluster_2000, levels = c("1", "2"))
T26_survival$H3K27 = factor(T26_survival$H3K27, levels = c("wildtype", "mutated"))
T26_survival$risk = factor(T26_survival$risk, levels = c("low", "intermediate",
                                                         "high"))
T26_survival$invasionhisto = factor(T26_survival$invasionhisto, levels = c("no", "yes"))

#first check assumption of Cox proportional hazards (i.e. Schoenfeld residuals are independent from time)
fit <- coxph(Surv(time, status)~cluster_2000 + setting + grading2021 + EOR +gender + RT + age + H3K27 +invasionhisto, data=T26_survival)
fit

test.ph = cox.zph(fit)
test.ph
ggcoxzph = ggcoxzph(test.ph, font.main=9)
ggcoxzph

#baseline validation in univariate analysis
#check survival pheno for grading, setting, MC groups, and risk score
sfit_consensus <- survfit(Surv(time, status)~MCconsensus, data=T26_survival)
ggsurvplot(sfit_consensus, risk.table=TRUE, palette=c("royalblue2","red3","forestgreen","darkorange2"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_consensus, risk.table=FALSE, palette=c("royalblue2","red3","forestgreen","darkorange2"),xlim=c(1,200), break.x.by=48)

fit <- coxph(Surv(time, status)~MCconsensus, data=T26_survival)
summary(fit)
fit

sfit_subtype <- survfit(Surv(time, status)~MCsubtype, data=T26_survival)
ggsurvplot(sfit_subtype, risk.table=TRUE, palette=c("#386CB0","#F0027F","#BF5B17"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_subtype, risk.table=FALSE, palette=c("#386CB0","#F0027F","#BF5B17"),xlim=c(1,200), break.x.by=48)

fit <- coxph(Surv(time, status)~MCsubtype, data=T26_survival)
fit

sfit_grade21 <- survfit(Surv(time, status)~grading2021, data=T26_survival)
ggsurvplot(sfit_grade21, risk.table=TRUE, palette=c("#9EBCDA","#8C6BB1","#810F7C"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_grade21, risk.table=FALSE, palette=c("#9EBCDA","#8C6BB1","#810F7C"),xlim=c(1,200), break.x.by=48)
fit <- coxph(Surv(time, status)~grading2021, data=T26_survival)
summary(fit)


sfit_setting <- survfit(Surv(time, status)~setting, data=T26_survival)
ggsurvplot(sfit_setting, risk.table=TRUE, palette=c("#1B9E77","#D95F02"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_setting, risk.table=FALSE, palette=c("#1B9E77","#D95F02"),xlim=c(1,200), break.x.by=48)
fit <- coxph(Surv(time, status)~setting, data=T26_survival)
summary(fit)

sfit_risk <- survfit(Surv(time, status)~risk, data=T26_survival)
ggsurvplot(sfit_risk, risk.table=TRUE, palette=c("dodgerblue2","purple1","firebrick1"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_risk, risk.table=FALSE, palette=c("dodgerblue2","purple1","firebrick1"),xlim=c(1,200), break.x.by=48)
fit <- coxph(Surv(time, status)~risk, data=T26_survival)
summary(fit)


#univariate analysis of all features and associated c-indices
fit <- coxph(Surv(time, status)~age, data=T26_survival)
summary(fit)
c_index_age = concordance(fit)$concordance
c_index_age

fit <- coxph(Surv(time, status)~setting, data=T26_survival)
summary(fit)
c_index_setting = concordance(fit)$concordance
c_index_setting

fit <- coxph(Surv(time, status)~grading2021, data=T26_survival)
summary(fit)
c_index_grade = concordance(fit)$concordance
c_index_grade


fit <- coxph(Surv(time, status)~EOR, data=T26_survival)
summary(fit)
c_index_EOR = concordance(fit)$concordance
c_index_EOR

fit <- coxph(Surv(time, status)~invasionhisto, data=T26_survival)
summary(fit)
c_index_invasion = concordance(fit)$concordance
c_index_invasion


str(T26_survival)
fit <- coxph(Surv(time, status)~H3K27, data=T26_survival)
summary(fit)
c_index_H3K27 = concordance(fit)$concordance
c_index_H3K27

str(T26_survival)
fit <- coxph(Surv(time, status)~cluster_2000, data=T26_survival)
summary(fit)
c_index_cluster = concordance(fit)$concordance
c_index_cluster

fit <- coxph(Surv(time, status)~MCconsensus, data=T26_survival)
summary(fit)
c_index_group = concordance(fit)$concordance
c_index_group

fit <- coxph(Surv(time, status)~gender, data=T26_survival)
summary(fit)
c_index_gender = concordance(fit)$concordance
c_index_gender

fit <- coxph(Surv(time, status)~RT, data=T26_survival)
summary(fit)
c_index_RT = concordance(fit)$concordance
c_index_RT


#split dataset by time for analysis of Gender and RT as Schoenfeld residuals are different
T26_split = survSplit(Surv(time, status)~.,
                      data = T26_survival, cut = c(36,72),
                      episode = "tgroup", id="id")
T26_split$RT = relevel(T26_split$RT, "yes")

fit <- coxph(Surv(tstart,time,status)~ RT:strata(tgroup), data=T26_split)
summary(fit)
c_index_RT = concordance(fit)$concordance
c_index_RT

T26_split$gender = relevel(T26_split$gender, "F")

fit <- coxph(Surv(tstart,time,status)~ gender:strata(tgroup), data=T26_split)
summary(fit)
c_index_gender = concordance(fit)$concordance
c_index_gender

######perform multivariate analysis for the entire cohort with potentially relevant covariates
#####gender and RT violate the PH assumption
####use step function to split gender and RT by time

T26_split = survSplit(Surv(time, status)~.,
                      data = T26_survival, cut = c(36,72),
                      episode = "tgroup", id="id")
head(T26_split)
str(T26_split)
T26_split$RT = relevel(T26_split$RT, "no")
T26_split$cluster_2000 = factor(T26_split$cluster_2000, levels = c("1", "2"))

fit.split <- coxph(Surv(tstart,time,status)~gender:strata(tgroup) + RT:strata(tgroup) +cluster_2000 + setting + grading2021 + EOR, data=T26_split)
summary(fit.split)
cox.zph(fit.split)
c_index_full_model_cluster = concordance(fit.split)$concordance
c_index_full_model_cluster

#sensitivity check with molecular groups
fit.split <- coxph(Surv(tstart,time,status)~gender:strata(tgroup) + RT:strata(tgroup) +cluster_2000 + MCconsensus + setting + grading2021 + EOR, data=T26_split)
summary(fit.split)
c_index_full_model_cluster_group = concordance(fit.split)$concordance
c_index_full_model_cluster_group

#check c-index with only group in multivariate
fit.split <- coxph(Surv(tstart,time,status)~gender:strata(tgroup) + RT:strata(tgroup) + MCconsensus + setting + grading2021 + EOR, data=T26_split)
summary(fit.split)
c_index_full_model_group = concordance(fit.split)$concordance
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


#check survival of all METHlow, separated by grades or groups, and same for METHhigh, excluding conditions with less than 10 survival points
#do multivariate analysis accounting for potential confounding factors

#METH-low
T26_survival_clust1 = read.csv(file = "T26_survival_cluster1.csv", header = T)
str(T26_survival_clust1)
T26_survival_clust1$MCconsensus = factor(T26_survival_clust1$MCconsensus, levels = c("Merlinintact", "Immuneenriched",
                                                                                     "hypermetabolic"))
T26_survival_clust1$grading2021 = factor(T26_survival_clust1$grading2021, levels = c("1", "2"))
T26_survival_clust1$EOR = factor(T26_survival_clust1$EOR, levels = c("GTR", "STR"))
T26_survival_clust1$RT = factor(T26_survival_clust1$RT, levels = c("no", "yes"))
T26_survival_clust1_split = survSplit(Surv(time, status)~.,
                                      data = T26_survival_clust1, cut = c(36,72),
                                      episode = "tgroup", id="id")

sfit_grade_clust1 <- survfit(Surv(time, status)~grading2021, data=T26_survival_clust1)
ggsurvplot(sfit_grade_clust1, risk.table=TRUE, palette=c("#9EBCDA","#8C6BB1"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_grade_clust1, risk.table=FALSE, palette=c("#9EBCDA","#8C6BB1"),xlim=c(1,200), break.x.by=48)

fit <- coxph(Surv(tstart,time,status)~grading2021 + gender:strata(tgroup) + setting + EOR + RT:strata(tgroup), data=T26_survival_clust1_split)
summary(fit)

sfit_consensus_clust1 <- survfit(Surv(time, status)~MCconsensus, data=T26_survival_clust1)
ggsurvplot(sfit_consensus_clust1, risk.table=TRUE, palette=c("royalblue2","red3","forestgreen"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_consensus_clust1, risk.table=FALSE, palette=c("royalblue2","red3","forestgreen"), 
           xlim=c(1,200), break.x.by=48)

fit <- coxph(Surv(tstart,time,status)~MCconsensus + grading2021 + gender:strata(tgroup) + setting + EOR + RT:strata(tgroup), data=T26_survival_clust1_split)
fit


#METH-high
T26_survival_clust2 = read.csv(file = "T26_survival_cluster2.csv", header = T)
T26_survival_clust2$MCconsensus = factor(T26_survival_clust2$MCconsensus, levels = c("hypermetabolic","proliferative"))
T26_survival_clust2$grading2021 = factor(T26_survival_clust2$grading2021, levels = c("2", "3"))
T26_survival_clust2$EOR = factor(T26_survival_clust2$EOR, levels = c("GTR", "STR"))
T26_survival_clust2$RT = factor(T26_survival_clust2$RT, levels = c("no", "yes"))
T26_survival_clust2_split = survSplit(Surv(time, status)~.,
                                      data = T26_survival_clust2, cut = c(36,72),
                                      episode = "tgroup", id="id")
str(T26_survival_clust2_split)

sfit_grade_clust2 <- survfit(Surv(time, status)~grading2021, data=T26_survival_clust2)
ggsurvplot(sfit_grade_clust2, risk.table=TRUE, palette=c("#8C6BB1","#810F7C"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_grade_clust2, risk.table=FALSE, palette=c("#8C6BB1","#810F7C"),conf.int=TRUE,surv.median.line = c("hv"), xlim=c(1,200), break.x.by=48)

fit <- coxph(Surv(tstart,time,status)~grading2021 + setting + EOR + gender:strata(tgroup) + RT:strata(tgroup), data=T26_survival_clust2_split)
summary(fit)


sfit_consensus_clust2 <- survfit(Surv(time, status)~MCconsensus, data=T26_survival_clust2)
ggsurvplot(sfit_consensus_clust2, risk.table=TRUE, palette=c("forestgreen","darkorange2"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_consensus_clust2, risk.table=FALSE, palette=c("forestgreen","darkorange2"),xlim=c(1,200), break.x.by=48)

fit <- coxph(Surv(tstart,time,status)~MCconsensus + grading2021 +EOR + RT +gender:strata(tgroup) + RT:strata(tgroup), data=T26_survival_clust2_split)
fit

#check superior separation of patients by cluster over grade, MC groups, and risk score

#get data for grade2
getwd()
setwd("/Users/lab/Desktop/Meningioma/T26_discovery_survival/data")
grade2_cluster = read.csv("Grade2_cluster_data.csv", header = T)
str(grade2_cluster)
grade2_cluster$grade_cluster = factor(grade2_cluster$grade_cluster, levels = c("grade2_1", "grade2_2"))
grade2_cluster$setting = factor(grade2_cluster$setting, levels = c("Primary", "Recurrence"))
grade2_cluster$EOR = factor(grade2_cluster$EOR, levels = c("GTR", "STR"))
grade2_cluster$RT = factor(grade2_cluster$RT, levels = c("no", "yes"))
grade2_cluster_split = survSplit(Surv(time, status)~.,
                                 data = grade2_cluster, cut = c(36,72),
                                 episode = "tgroup", id="id")

sfit_grade_cluster <- survfit(Surv(time, status)~grade_cluster, data=grade2_cluster)
ggsurvplot(sfit_grade_cluster, risk.table=TRUE, palette=c("#8C6BB1","#8C6BB1"),conf.int=TRUE, 
           risk.table.height=.35, linetype = c("dashed", "solid"), surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_grade_cluster, risk.table=FALSE, palette=c("cyan4","violetred4"),conf.int=TRUE,
           size = 1.8,surv.median.line = c("hv"),conf.int.alpha=.2,xlim=c(1,200), break.x.by=48)

#correct for confounding factors
fit <- coxph(Surv(tstart,time,status)~grade_cluster + gender:strata(tgroup) + setting + EOR + RT:strata(tgroup), data=grade2_cluster_split)
fit

#get data for hypermetabolic
hypermetabolic_cluster = read.csv("Hypermetabolic_cluster_data.csv", header = T)
str(hypermetabolic_cluster)
hypermetabolic_cluster$consensus_cluster = factor(hypermetabolic_cluster$consensus_cluster, 
                                                  levels = c("hypermetabolic_1", "hypermetabolic_2"))
hypermetabolic_cluster$setting = factor(hypermetabolic_cluster$setting, levels = c("Primary", "Recurrence"))
hypermetabolic_cluster$EOR = factor(hypermetabolic_cluster$EOR, levels = c("GTR", "STR"))
hypermetabolic_cluster$RT = factor(hypermetabolic_cluster$RT, levels = c("no", "yes"))
hypermetabolic_cluster$grading2016 = factor(hypermetabolic_cluster$grading2016, levels = c("1", "2"))
hypermetabolic_cluster_split = survSplit(Surv(time, status)~.,
                                         data = hypermetabolic_cluster, cut = c(36,72),
                                         episode = "tgroup", id="id")

sfit_hypermetabolic_cluster <- survfit(Surv(time, status)~consensus_cluster, data=hypermetabolic_cluster)
ggsurvplot(sfit_hypermetabolic_cluster, risk.table=TRUE, palette=c("forestgreen","forestgreen"),conf.int=TRUE, 
           risk.table.height=.35, linetype = c("dashed", "solid"), surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_hypermetabolic_cluster, risk.table=FALSE, palette=c("cyan4","violetred4"),conf.int=TRUE,
           size=1.8,surv.median.line = c("hv"),conf.int.alpha=.2,xlim=c(1,200), break.x.by=48)

fit <- coxph(Surv(time, status)~consensus_cluster, data=hypermetabolic_cluster)
fit
#correct for confounding factors
fit <- coxph(Surv(tstart,time,status)~consensus_cluster + gender:strata(tgroup) + setting + EOR + RT:strata(tgroup) + grading2016, data=hypermetabolic_cluster_split)
fit

#get data for risk
getwd()
setwd("/Users/lab/Desktop/Meningioma/T26_discovery_survival/data")
risk_cluster = read.csv("Intermediate_risk_cluster_data.csv", header = T)
str(risk_cluster)
risk_cluster$risk_cluster = factor(risk_cluster$risk_cluster, 
                                   levels = c("intermediate_1", "intermediate_2"))
risk_cluster$setting = factor(risk_cluster$setting, levels = c("Primary", "Recurrence"))
risk_cluster$EOR = factor(risk_cluster$EOR, levels = c("GTR", "STR"))
risk_cluster$RT = factor(risk_cluster$RT, levels = c("no", "yes"))
risk_cluster$grading2021 = factor(risk_cluster$grading2021, levels = c("1", "2"))
risk_cluster_split = survSplit(Surv(time, status)~.,
                               data = risk_cluster, cut = c(36,72),
                               episode = "tgroup", id="id")

sfit_risk_cluster <- survfit(Surv(time, status)~risk_cluster, data=risk_cluster)
ggsurvplot(sfit_risk_cluster, risk.table=TRUE, palette=c("purple1","purple1"),conf.int=TRUE, 
           risk.table.height=.35, linetype = c("dashed", "solid"), surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_risk_cluster, risk.table=FALSE, palette=c("cyan4","violetred4"),conf.int=TRUE,
           size=1.8,surv.median.line = c("hv"),conf.int.alpha=.2,xlim=c(1,200), break.x.by=48)

fit <- coxph(Surv(time, status)~risk_cluster, data=risk_cluster)
fit
#correct for confounding factors
fit <- coxph(Surv(tstart,time,status)~risk_cluster + gender:strata(tgroup) + setting + EOR + RT:strata(tgroup) + grading2021, data=risk_cluster_split)
summary(fit)


####1p deletions are significantly associated with genome instability and CNVs
###check effect of 1p deletions on progression in METHlow and METHhigh
#get data
chr1p_cluster = read.csv("1p_cluster_data.csv", header = T)
str(chr1p_cluster)
chr1p_cluster$setting = factor(chr1p_cluster$setting, levels = c("Primary", "Recurrence"))
chr1p_cluster$EOR = factor(chr1p_cluster$EOR, levels = c("GTR", "STR"))
chr1p_cluster$RT = factor(chr1p_cluster$RT, levels = c("no", "yes"))
chr1p_cluster$grading2021 = factor(chr1p_cluster$grading2021, levels = c("1", "2","3"))
chr1p_cluster$X1p_cluster = factor(chr1p_cluster$X1p_cluster, 
                                   levels = c("deletion_1", "deletion_2","intact_1"))

chr1p_cluster_split = survSplit(Surv(time, status)~.,
                                data = chr1p_cluster, cut = c(36,72),
                                episode = "tgroup", id="id")

sfit_chr1p_cluster <- survfit(Surv(time, status)~X1p_cluster, data=chr1p_cluster)
ggsurvplot(sfit_chr1p_cluster, risk.table=TRUE, palette=c("#00D7D1","#008080","violetred4"),
           conf.int=TRUE, linetype = c("dashed", "solid","solid"),
           risk.table.height=.35, surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_chr1p_cluster, risk.table=FALSE, palette=c("#00D7D1","#008080","violetred4"),conf.int=TRUE,
           linetype = c("solid", "dashed","dashed"),
           size=1.8,surv.median.line = c("hv"),conf.int.alpha=.2,xlim=c(1,200), break.x.by=48)

#correct for confounding factors
fit <- coxph(Surv(tstart,time,status)~X1p_cluster + gender:strata(tgroup) + setting + EOR + RT:strata(tgroup) + grading2021, data=chr1p_cluster_split)
summary(fit)




###check effect of genome instability and unfavorable CNVs on progression in the entire TUE cohort
#get data
setwd("/Users/lab/Desktop/Meningioma/T26_discovery_survival/data")
T26_survival = read.csv(file = "T26_survival.csv", header = T)
str(T26_survival)
T26_survival$gender = factor(T26_survival$gender, levels = c("M", "F"))
T26_survival$RT = factor(T26_survival$RT, levels = c("no", "yes"))
T26_survival$setting = factor(T26_survival$setting, levels = c("Primary", "Recurrence"))
T26_survival$MCconsensus = factor(T26_survival$MCconsensus, levels = c("Merlinintact", "Immuneenriched",
                                                                       "hypermetabolic","proliferative"))
T26_survival$MCsubtype = factor(T26_survival$MCsubtype, levels = c("Benign", "Intermediate",
                                                                   "Malignant"))
T26_survival$grading2021 = factor(T26_survival$grading2021, levels = c("1", "2", "3"))
T26_survival$cluster_2000 = factor(T26_survival$cluster_2000, levels = c("1", "2"))
T26_survival$CNV_intervals = factor(T26_survival$CNV_intervals, levels = c("zero", "up2",
                                                                           "more2"))
T26_survival$stability_intervals = factor(T26_survival$stability_intervals, levels = c("less1", "less10",
                                                                                       "less20", "more20"))

T26_split = survSplit(Surv(time, status)~.,
                      data = T26_survival, cut = c(36,72),
                      episode = "tgroup", id="id")

#check instability as a factor
sfit_T26_instability <- survfit(Surv(time, status)~stability_intervals, data=T26_survival)
ggsurvplot(sfit_T26_instability, risk.table=TRUE, palette=c("grey35","indianred1","darkorchid3", "firebrick3"),
           conf.int=TRUE, 
           risk.table.height=.35, surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_T26_instability, risk.table=FALSE, palette=c("grey35","indianred1","darkorchid3", "firebrick3"),conf.int=TRUE,
           size=1.8,surv.median.line = c("hv"),conf.int.alpha=.2,xlim=c(1,200), break.x.by=48)

#correct for confounding factors
fit <- coxph(Surv(tstart,time,status)~stability_intervals + cluster_2000 +gender:strata(tgroup) + setting + EOR + RT:strata(tgroup) + grading2021, data=T26_split)
summary(fit)


#check number of unfavorable CNVs as a factor
sfit_T26_CNVs <- survfit(Surv(time, status)~CNV_intervals, data=T26_survival)
ggsurvplot(sfit_T26_CNVs, risk.table=TRUE, palette=c("grey35","indianred1", "firebrick3"),
           conf.int=TRUE, 
           risk.table.height=.35, surv.median.line = c("hv"),xlim=c(1,200), break.x.by=48)
ggsurvplot(sfit_T26_CNVs, risk.table=FALSE, palette=c("grey35","indianred1", "firebrick3"),conf.int=TRUE,
           size=1.8,surv.median.line = c("hv"),conf.int.alpha=.2,xlim=c(1,200), break.x.by=48)

#correct for confounding factors
fit <- coxph(Surv(tstart,time,status)~CNV_intervals + cluster_2000 +gender:strata(tgroup) + setting + EOR + RT:strata(tgroup) + grading2021, data=T26_split)
summary(fit)




######test Brier prediction with random reference using pec package
# Time points at which to evaluate Brier scores
times <- seq(0, 214, by = 3)

fit_grade <- coxph(Surv(time, status)~grading2021, data=T26_survival, x = TRUE, y = TRUE)
fit_cluster <- coxph(Surv(time, status)~cluster_2000, data=T26_survival, x = TRUE, y = TRUE)
fit_group <- coxph(Surv(time,status)~ MCconsensus, data=T26_survival, x = TRUE, y = TRUE)


# Compute Brier score for Cox model
brier_model_grade <- pec::pec(
  object = list("Cox" = fit_grade),
  formula = Surv(time, status) ~ 1,
  data = T26_survival,
  times = times,
  exact = FALSE
)
ibs_grade <- crps(brier_model_grade)
ibs_grade
#integrated is 0.187

brier_model_cluster <- pec::pec(
  object = list("Cox" = fit_cluster),
  formula = Surv(time, status) ~ 1,
  data = T26_survival,
  times = times,
  exact = FALSE
)
ibs_cluster <- crps(brier_model_cluster)
ibs_cluster
#integrated is 0.175

brier_model_group <- pec::pec(
  object = list("Cox" = fit_group),
  formula = Surv(time, status) ~ 1,
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
Brier_scores = cbind(Brier_grade, Brier_cluster, Brier_group,Brier_risk, sim_result)
Brier_scores = as.data.frame(Brier_scores)
Brier_scores$times = rownames(Brier_scores)
colnames(Brier_scores) = c("grade", "cluster", "group","risk", "random", "time")
write.csv(Brier_scores, file="Brier_scores_full_cohort.csv")


####make the same again, but only for WHO grade 2 cases in T26_survival

T26_survival_grade2 = T26_survival[T26_survival$grading2021 == "2",]

# Time points at which to evaluate Brier scores
times <- seq(0, 214, by = 3)

fit_cluster_grade2 <- coxph(Surv(time, status)~cluster_2000, data=T26_survival_grade2, x = TRUE, y = TRUE)
fit_group_grade2 <- coxph(Surv(time,status)~ MCconsensus, data=T26_survival_grade2, x = TRUE, y = TRUE)
fit_risk_grade2 <- coxph(Surv(time, status)~risk, data=T26_survival_grade2, x = TRUE, y = TRUE)

# Compute Brier score for Cox model
brier_model_cluster <- pec::pec(
  object = list("Cox" = fit_cluster_grade2),
  formula = Surv(time, status) ~ 1,
  data = T26_survival_grade2,
  times = times,
  exact = FALSE
)
ibs_cluster_grade <- crps(brier_model_cluster)
ibs_cluster_grade
#integrated is 0.155

brier_model_group <- pec::pec(
  object = list("Cox" = fit_group_grade2),
  formula = Surv(time, status) ~ 1,
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
Brier_risk_grade2 <- brier_model_risk$AppErr$Cox


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
Brier_scores_grade2 = cbind(Brier_cluster_grade2, Brier_group_grade2, Brier_risk_grade2, sim_result)
Brier_scores_grade2 = as.data.frame(Brier_scores_grade2)
Brier_scores_grade2$times = rownames(Brier_scores_grade2)
colnames(Brier_scores_grade2) = c("cluster", "group","risk", "random", "time")
write.csv(Brier_scores_grade2, file="Brier_scores_WHO_grade2.csv")



























