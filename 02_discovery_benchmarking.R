library(sesame)
library(ggplot2)
library(survival)
library(survminer)

#get processed beta values
betas <- readRDS("Betas_discovery_preprocessed_corrected.rds")
meta = read.csv(file="meta_discovery_current.csv")


#estimate leukocyte infiltration to validate immuneenriched group and in general low/high grade meningiomas
#subset to probes in EPIC
dim(betas)
betas_EPIC = mLiftOver(betas, "EPIC")
dim(betas_EPIC)
betas_EPIC = betas_EPIC[complete.cases(betas_EPIC),]

# save corrected beta values matched to EPIC array to .RDS file for later use
saveRDS(betas_EPIC, file = "Betas_discovery_preprocessed_corrected_EPIC_match.RDS")


#estimate immune cell infiltration
Leuko.infiltration <- apply(betas_EPIC, 2, estimateLeukocyte)
infiltration = as.data.frame(Leuko.infiltration)
infiltration$ID = rownames(infiltration)


all(rownames(infiltration) == meta$ID)
infiltration$grade = meta$grading2021_new
infiltration$grade = factor(infiltration$grade)
infiltration$class = meta$MCsubtype
infiltration$class = factor(infiltration$class, levels = c("Benign", "Intermediate", "Malignant"))
infiltration$group = meta$mol_group
infiltration$group = factor(infiltration$group, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))



cairo_pdf(filename = "T26_discovery_infiltration_by_grade.pdf", width = 2.8, height = 9)
ggplot(infiltration, 
       aes(grade, Leuko.infiltration,fill = grade, alpha = 0.98))+
  scale_fill_manual(values=c("#9EBCDA","#8C6BB1","#810F7C"))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(), alpha=0.55, aes(color=grade))+
  scale_color_manual(values=c("#9EBCDA","#8C6BB1","#810F7C"))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


cairo_pdf(filename = "T26_discovery_infiltration_by_group.pdf", width = 3.5, height = 9)
ggplot(infiltration, 
       aes(group, Leuko.infiltration,fill = group, alpha = 0.98))+
  scale_fill_manual(values=c("royalblue2","red3","forestgreen","darkorange2"))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(), alpha=0.55, aes(color=group))+
  scale_color_manual(values=c("royalblue2","red3","forestgreen","darkorange2"))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()







#baseline validation in univariate analysis
#check survival pheno for grading, setting, MC groups, and risk score

#get PFS data
T26_survival = read.csv(file = "T26_survival.csv", header = T)
str(T26_survival)
T26_survival$grading2021_new = factor(T26_survival$grading2021_new, levels = c("1", "2", "3"))
T26_survival$setting = factor(T26_survival$setting, levels = c("Primary", "Recurrence"))
T26_survival$MCconsensus = factor(T26_survival$MCconsensus, levels = c("Merlinintact", "Immuneenriched",
                                                                       "hypermetabolic","proliferative"))
T26_survival$risk_new = factor(T26_survival$risk_new, levels = c("low", "intermediate",
                                                         "high"))
str(T26_survival)


#first check assumption of Cox proportional hazards (i.e. Schoenfeld residuals are independent from time)
fit <- coxph(Surv(time, status)~grading2021_new + setting + MCconsensus + risk, data=T26_survival)
fit
test.ph = cox.zph(fit)
test.ph
ggcoxzph = ggcoxzph(test.ph, font.main=9)
ggcoxzph
#all factors fullfill the Cox requirement of proportinal hazards



#Kaplan-Meier curves and univariate analysis Cox regression of survival
sfit_grade21 <- survfit(Surv(time, status)~grading2021_new, data=T26_survival)
ggsurvplot(sfit_grade21, risk.table=TRUE, palette=c("#9EBCDA","#8C6BB1","#810F7C"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)

cairo_pdf(filename = "Surv_discovery_grade_univariate.pdf", width = 9, height = 7)
ggsurvplot(sfit_grade21, risk.table=FALSE, palette=c("#9EBCDA","#8C6BB1","#810F7C"),xlim=c(1,200), break.x.by=48)
dev.off()

fit <- coxph(Surv(time, status)~grading2021_new, data=T26_survival)
summary(fit)

sfit_setting <- survfit(Surv(time, status)~setting, data=T26_survival)
ggsurvplot(sfit_setting, risk.table=TRUE, palette=c("#1B9E77","#D95F02"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)

cairo_pdf(filename = "Surv_discovery_setting_univariate.pdf", width = 9, height = 7)
ggsurvplot(sfit_setting, risk.table=FALSE, palette=c("#1B9E77","#D95F02"),xlim=c(1,200), break.x.by=48)
dev.off()

fit <- coxph(Surv(time, status)~setting, data=T26_survival)
summary(fit)

sfit_consensus <- survfit(Surv(time, status)~MCconsensus, data=T26_survival)
ggsurvplot(sfit_consensus, risk.table=TRUE, palette=c("royalblue2","red3","forestgreen","darkorange2"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)

cairo_pdf(filename = "Surv_discovery_group_univariate.pdf", width = 9, height = 7)
ggsurvplot(sfit_consensus, risk.table=FALSE, palette=c("royalblue2","red3","forestgreen","darkorange2"),xlim=c(1,200), break.x.by=48)
dev.off()

fit <- coxph(Surv(time, status)~MCconsensus, data=T26_survival)
summary(fit)

sfit_risk <- survfit(Surv(time, status)~risk_new, data=T26_survival)
ggsurvplot(sfit_risk, risk.table=TRUE, palette=c("dodgerblue2","purple1","firebrick1"), 
           risk.table.height=.35,xlim=c(1,200), break.x.by=48)

cairo_pdf(filename = "Surv_discovery_risk_univariate.pdf", width = 9, height = 7)
ggsurvplot(sfit_risk, risk.table=FALSE, palette=c("dodgerblue2","purple1","firebrick1"),xlim=c(1,200), break.x.by=48)
dev.off()

fit <- coxph(Surv(time, status)~risk_new, data=T26_survival)
summary(fit)



