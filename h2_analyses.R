library(readxl)
library(tidyverse)
library(car)
library(lme4)
library(lmerTest)
library(ggpubr)
library(merTools)

#######################-
## Days until Diapause ----

#######################-

### Import data of diapause experiment ----
y.n = read_xlsx("G:\\My Drive\\GRAD SCHOOL\\RESEARCH\\D. carinulata\\Diapause analyses and manuscript\\Diapause heritability.2019\\data entry\\dipause_yn_8.9.19.xlsx",
                sheet = "import to R (2)")

### Data Clean up

y.n <- y.n %>%
  filter(is.na(diapause) == F ) %>% # filter out non-trials
  select(1:7) # remove columns for each day
y.n$Sire <- as.factor(y.n$Sire)
str(y.n)

y.n %>% group_by(Treatment, Sire) %>% summarise(
  noff = n()) %>% summarise(
    range = range(noff),
    nSire = n(),
    mean = mean(noff)
  )

short <- y.n %>% # create dataset for short treatment
  filter(Treatment == "S" & is.na(dead) == TRUE) # filter out dead beetles 
hist(short$time, breaks = 50)

long <- y.n %>% # create dataset for long treatment
  filter(Treatment == "L" & is.na(dead) == TRUE) # filter out dead beetles 
hist(long$time, breaks = 50)
str(long)

### Calculate heritability with homemade function ----
h2.short = heritability.phs(response = short$time, sire = short$Sire)
h2.long = heritability.phs(response = long$time, sire = long$Sire)

### Plot results - old ----
# baseR
fam.means = as.data.frame(cbind(h2.long$z, h2.short$z))
fam.means.sort = fam.means[order(fam.means$V1),]

plot(fam.means.sort$V1, pch = 19, ylab = "", xlab = "", ylim = c(0,40), 
     col = rgb(0/255, 160/255, 215/255), xaxt = "n")
points(fam.means.sort$V2, pch = 1, col = rgb(0/255, 160/255, 215/255))
legend(x = "topleft", y = 35, legend = c("Long days (home)", "Short days (away)"), bty = "n",
       pch = c(19,1), col = rgb(0/255, 160/255, 215/255))
title(ylab="days until diapause", mgp=c(2.3,1,0), cex.lab=1.5)
title(xlab="rank order of families", mgp=c(.8,1,0), cex.lab=1.5)

### Plot - correlation
plot(fam.means$V1, fam.means$V2, ylab = "family mean short", xlab = "family mean long")
cor.test(fam.means$V1, fam.means$V2)
cor.model = lm(fam.means$V1 ~ fam.means$V2)
summary(cor.model)

### Plot - h2 in each enviro
{plot(c(1,2), c(h2.long$h2, h2.short$h2), pch = c(19,1), cex = 2.2,
      col = rgb(0/255, 160/255, 215/255),
      ylim = c(-.75,1.5), ylab = "Heritability (95% CI)", cex.lab = 1.5,
      xlim = c(0.5,2.5), xaxt = 'n', xlab = "Environment",
      panel.first = c(abline(h=seq(0,1,.2), col = "grey"),
                      arrows(2, (h2.short$CI.95.SE[1]), 2, (h2.short$CI.95.SE[2]), length = 0.1, angle = 90, code = 3),
                      arrows(1, (h2.long$CI.95.SE[1]), 1, (h2.long$CI.95.SE[2]), length = 0.1, angle = 90, code = 3)))
  axis(side = 1, at = c(1,2), labels = c("Long (home)", "Short (away)"), cex.axis = 1.2)
  abline(h=c(0,1), lty = 2)}

mean(long$time)
sd(long$time)
var(long$time) * h2.long$h2.corrected

mean(short$time)
sd(short$time)
var(short$time) * h2.short$h2.corrected

t.test(long$time, short$time)

##histograms
hist(long$time, col = rgb(0.1,0.1,0.1,0.5), breaks = 50, ylim = c(0,50))
hist(short$time, col = rgb(0.8,0.8,0.8,0.5), breaks = 50, add = TRUE)


### 1-way ANOVA ----
modelShort1 <- lm(time ~ Sire, data = short)
summary(modelShort1)
Anova(modelShort1, type = 3)
anova(modelShort1)
plot(modelShort1)

modelLong1 <- lm(time ~ Sire, data = long)
summary(modelLong1)
Anova(modelLong1, type = 3)
anova(modelLong1)
plot(modelLong1)
confint(modelLong1)

### Random effects Models ----
#Long day treatment
days_meanL <- long %>% summarise(mean = mean(time)) %>% as.numeric()
mixedLong1 <- lmer(time ~ (1 | Sire), data = long)
summary(mixedLong1)
# confint(mixedLong1, method = "boot", oldNames = F)
confint(mixedLong1, method = "profile", oldNames = F)
rand(mixedLong1)
days_res_L <- var_from_model_phs(mixedLong1, days_meanL)
days_res_forplot_L <- days_res_L$VarCompsLong %>%
  mutate(treatment = "long")

#Short day treatment
days_meanS <- short %>% summarise(mean = mean(time)) %>% as.numeric()
mixedShort1 <- lmer(time ~ (1 | Sire), data = short)
summary(mixedShort1)
# confint(mixedShort1, method = "boot")
confint(mixedShort1, method = "profile")
rand(mixedShort1)
days_res_S <- var_from_model_phs(mixedShort1, days_meanS)
days_res_forplot_S <- days_res_S$VarCompsLong %>%
  mutate(treatment = "short")

###Days Bootstrap confidence intervals ----
#long
boot_days_L <- bootstrap_VaVp(data = long, n_groups = n_groups, n_iter = n_iter, mean = days_meanL) %>% 
  bind_rows() %>% filter(is.na(var_sire) == F) #combine list into one dataframe, remove NAs
boot_sum_days_L <- boot_summary(boot_days_L) %>%
  mutate(treatment = "long")

#short
boot_days_S <- bootstrap_VaVp(data = short, n_groups = n_groups, n_iter = n_iter, mean = days_meanS) %>% 
  bind_rows() %>% filter(is.na(var_sire) == F) #combine list into one dataframe, remove NAs
boot_sum_days_S <- boot_summary(boot_days_S) %>%
  mutate(treatment = "short")

###Plots of Results ----
days_res_comb <- rbind(days_res_forplot_L, days_res_forplot_S) # means from REML

days_boot_sum_both <- rbind(boot_sum_days_L,boot_sum_days_S) # confidence intervals from bootstrap
days_boot_sum_both$var_comp <- factor(days_boot_sum_both$var_comp, levels = c('VP', 'VA', 'h2', 'I')) #change order of VA and VP
levels(days_boot_sum_both$var_comp) <- c('V[Phenotypic]', 'V[Additive]', 'h^2', 'I[A]') #change labels for VA and VP

Vars <- ggplot(data = days_boot_sum_both %>% filter(var_comp != 'h^2' & var_comp != 'I[A]'), aes(x = treatment, y = mean)) +
  geom_linerange(aes(ymax = mean + sd, ymin = mean - sd),
                 size = 1, lty = 'solid') +
  geom_point(data = days_res_comb %>% filter(var_comp != 'h^2' & var_comp != 'I[A]'), aes(x = treatment, y = est, shape = treatment), size = 3, fill = 'white') +
  facet_wrap( ~ var_comp, labeller = label_parsed, scales = "fixed") +
  labs(y = "Variance of\ndays until diapause", x = "Environment") +
  scale_x_discrete(labels = c('long' = "North\n(Home)", 'short' = "South\n(Away)"))+
  scale_shape_manual(labels = c('long' = "Home", 'short' = "Away"), values = c('long' = 21, 'short' = 22)) +
  theme_bw() +
  theme(text = element_text(size = 18), #axis.ticks.x = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none')
Vars

Herit <- ggplot(data = days_boot_sum_both %>% filter(var_comp == 'h^2'), aes(x = treatment, y = mean)) +
  geom_linerange(aes(ymax = mean + sd, ymin = mean - sd),
                 size = 1, lty = 'solid') +
  geom_point(data = days_res_comb %>% filter(var_comp == 'h^2'), aes(x = treatment, y = est, shape = treatment), size = 3, fill = 'white') +
  facet_wrap( ~ var_comp, labeller = label_parsed, scales = "fixed") +
  labs(y = "Heritability of\ndays until diapause", x = "Environment") +
  scale_x_discrete(labels = c('long' = "North\n(Home)", 'short' = "South\n(Away)")) +
  scale_shape_manual(values = c('long' = 21, 'short' = 22)) +
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(-0.05, 1.14)) +
  theme_bw() +
  theme(text = element_text(size = 18), #axis.ticks.x = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none')
Herit

Evolv <- ggplot(data = days_boot_sum_both %>% filter(var_comp == 'I[A]'), aes(x = treatment, y = mean)) +
  geom_linerange(aes(ymax = mean + sd, ymin = mean - sd),
                 size = 1, lty = 'solid') +
  geom_point(data = days_res_comb %>% filter(var_comp == 'I[A]'), aes(x = treatment, y = est, shape = treatment), size = 3, fill = 'white') +
  facet_wrap( ~ var_comp, labeller = label_parsed, scales = "fixed") +
  labs(y = "Evolvability of\ndays until diapause", x = "Environment") +
  scale_x_discrete(labels = c('long' = "North\n(Home)", 'short' = "South\n(Away)")) +
  scale_shape_manual(values = c('long' = 21, 'short' = 22)) +
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(-0.05, 1.14)) +
  theme_bw() +
  theme(text = element_text(size = 18), #axis.ticks.x = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none')
Evolv

ggarrange(Vars, Herit, Evolv, ncol = 3, widths = c(.8, 0.5, 0.5), labels = 'AUTO')
##Export: 900 x 450 px (9.375 x 4.69 in)

### Results for plot - OLD
# long.res <- as.data.frame(VarCorr(mixedLong1))
# short.res <- as.data.frame(VarCorr(mixedShort1))
# 
# REML.results <- data.frame(var_sire = NA, var_e = NA, V_A = NA, V_P = NA, h2 = NA) #blank data.frame to store results
# REML.results <- REML.results %>% add_row(var_sire = c(long.res[1,'vcov']), #store results and calculate VA & VP & h2
#                                          var_e = c(long.res[2,'vcov']), 
#                                          V_A = c(4 * long.res[1,'vcov']), 
#                                          V_P = c(long.res[1,'vcov'] + long.res[2,'vcov']),
#                                          h2 = c(V_A / V_P)) 
# REML.results <- REML.results %>% add_row(var_sire = short.res[1,'vcov'], #store results and calculate VA & VP & h2
#                                          var_e = short.res[2,'vcov'], 
#                                          V_A = 4 * var_sire, 
#                                          V_P = var_sire + var_e,
#                                          h2 = V_A / V_P) %>%
#   filter(is.na(var_sire) == F) %>%
#   mutate(treatment = c('long', 'short')) %>%
#   dplyr::select(c(V_A:treatment)) %>%
#   pivot_longer(cols = c(1:3), names_to = "var_comp", values_to = "est")
# options(digits = 10)
# REML.results$var_comp <- factor(REML.results$var_comp, levels = c('V_P', 'V_A', 'h2')) #change order of VA and VP
# levels(REML.results$var_comp) <- c('V[Phenotypic]', 'V[Additive]', 'h^2') #change labels for VA and VP


### Plot of variation among families ----

resim_short <- REsim(mixedShort1) %>% mutate(treat = 'short') %>% #short dataset
  mutate(addIntercept = mean + fixef(mixedShort1))
resim_long <- REsim(mixedLong1) %>% mutate(treat = 'long') %>% #long dataset
  mutate(addIntercept = mean + fixef(mixedLong1)) %>% 
  arrange(addIntercept) %>%
  mutate(groupID2 = factor(groupID, levels = groupID)) 
str(resim_long)

# plotREsim(REsim(mixedLong1))
# plotREsim(REsim(mixedShort1))

varAmongFams <- ggplot(data = resim_long, aes(x = groupID2, y = addIntercept)) +
  # geom_point(aes(fill = treat), size = 3, shape = 21) +  #For means, without sd
  # geom_point(data = resim_short, aes(x = groupID, y = addIntercept, fill = treat), size = 3, shape = 21) +
  geom_pointrange(aes(ymin = addIntercept - sd, ymax = addIntercept + sd, shape = treat), #For means with sd
                  fatten = 2, size = 1, fill = 'white') + #shape = 21) +
  geom_pointrange(data = resim_short, aes(x = groupID, ymin = addIntercept - sd, ymax = addIntercept + sd, shape = treat),
                  fatten = 2, size = 1, fill = 'white') + #, shape = 21) +
  lims(y = c(0,31)) +
  labs(x = "Rank order of families", y = "Days until diapause", shape = "Environment") +
  scale_shape_manual(labels = c('long' = "North fall (home)", 'short' = "South fall (away)"), values = c(21, 22))+# c('long' = , 'short' = '2')) +
  theme_classic() +
  theme(text = element_text(size = 18), axis.ticks.x = element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'right')
varAmongFams
##Export: 700 x 400 px (7.21 x 4.17 in)

# ggarrange(Vars, Herit, Evolv, varAmongFams, ncol = 4, widths = c(1, 0.5, 0.5, 1), labels = 'AUTO')


### Trying out other models ----
mixedLong2 <- glmer(time ~ (1 | Sire), family = "poisson", data = long)
summary(mixedLong2)

##
longSummary <- long %>% group_by(Sire) %>% 
  summarise(
  n = n(),
  meanTime = mean(time)
)

N.long = length(longSummary$Sire)
T.long = sum(longSummary$n)
n0.long = (T.long - (sum(longSummary$n^2)/T.long))/(N.long-1)

Vars.long = (MSs - MSe)/n0.long

# Mixed model with both treatments
modelboth1 <- lm(time ~ Sire * Treatment, data = y.n)
Anova(modelboth1, type = 3)

modelboth1 <- lmer(time ~ Treatment * (1 | Sire) + (1 | Dam/Sire), data = y.n)
summary(modelboth1)
rand(modelboth1)
Anova(modelboth1, type = 3)
anova(modelboth1, type = 3)
intervals(modelboth1)

####
library(fullfact)
observLmer(observ = long, dam = "Dam", sire = "Sire", response = "time", ml = F)

#######################-
# Weight at eclosion ----

#######################-

### Import weight data and clean up ----

weight_data <- read_xlsx("G:\\My Drive\\GRAD SCHOOL\\RESEARCH\\D. carinulata\\Diapause analyses and manuscript\\Diapause heritability.2019\\data entry\\weight_width-8.9.19.xlsx",
                sheet = "weight_import")

weight_data <- weight_data %>%
  filter(is.na(weight) == F ) %>% # filter out empty rows
  mutate(sex = case_when(         # Make new column for sex
    treatment == 'M' ~ 'M',
    treatment == 'L' | treatment == 'S' ~ 'F'
  ))
weight_data$sire <- as.factor(weight_data$sire)
weight_data$dam <- as.factor(weight_data$dam)
weight_data$treatment <- as.factor(weight_data$treatment)
weight_data$sex <- as.factor(weight_data$sex)
str(weight_data)

weight_data %>% group_by(sire,dam) %>% summarise(
  noff = n()) %>% summarise(
    range = range(noff),
    nSire = n(),
    mean = mean(noff)
  )

heritability.fs(weight_data$weight, weight_data$sire, weight_data$dam)

#2-way nested ANOVA
weight_anova1 <- lm(weight ~ sex + sire + sire/dam, data = weight_data)
summary(weight_anova1)
Anova(weight_anova1, type = 3)
anova(weight_anova1)
plot(weight_anova1)

weight_anovaF <- lm(weight ~ sire + sire/dam, data = weight_data %>% filter(sex == "F"))
anova(weight_anovaF)
weight_data_F <- weight_data %>% filter(sex == "F")
heritability.fs(weight_data_F$weight, weight_data_F$sire, weight_data_F$dam)
# anova(lm(weight ~ sire, data = weight_data %>% filter(treatment == "L")))

weight_anovaM <- lm(weight ~ sire, data = weight_data %>% filter(sex == "M"))
anova(weight_anovaM)
weight_data_M <- weight_data %>% filter(sex == "M")
heritability.phs(weight_data_M$weight, weight_data_M$sire)


### Random effects Models ----
weight_mixed1 <- lmer(weight ~ sex + (1 | sire/dam), data = weight_data)
summary(weight_mixed1)
# confint(mixedLong1, method = "boot", oldNames = F)
confint(weight_mixed1, method = "profile", oldNames = F)
rand(weight_mixed1)

# Females only 
weight_meanF <- weight_data %>% filter(sex == "F") %>% summarise(mean = mean(weight)) %>% as.numeric()
weight_mixedF <- lmer(weight ~ (1 | sire/dam), data = weight_data %>% filter(sex == "F"))
summary(weight_mixedF)
rand(weight_mixed1)
fem_weight_res <- var_from_model_fs(weight_mixedF, mean = weight_meanF) ##Results from REML
fem_weight_forplot <- fem_weight_res$VarCompsLong %>%
  mutate(sex = "F")

# Males only
weight_meanM <- weight_data %>% filter(sex == "M") %>% summarise(mean = mean(weight)) %>% as.numeric()
weight_mixedM <- lmer(weight ~ (1 | sire), data = weight_data %>% filter(sex == "M"))
summary(weight_mixedM)
rand(weight_mixedM)
mal_weight_res <- var_from_model_phs(weight_mixedM, mean = weight_meanM) ##Results from REML
mal_weight_forplot <- mal_weight_res$VarCompsLong %>%
  mutate(sex = "M")

### Weight - Bootstrap confidence intervals ----
#NOTE: Bootstrap parameters are set on 'bootstrap_va_vp file

#Bootstrap sd - females only
weight_data %>% filter(sex == "F") %>% group_by(sire) %>% sample_n_groups(n_groups) %>% View()
# bootstrap_VaVp_weight_fs
boot_weight_f <- bootstrap_VaVp_weight_fs(data = weight_data %>% filter(sex == "F"), n_groups = n_groups, n_iter = n_iter, mean = weight_meanF) %>% 
  bind_rows() %>% filter(is.na(var_sire) == F)
boot_sum_F <- boot_summary(boot_weight_f) %>%
  mutate(sex = "F")

#Bootstrap sd - males only
boot_weight_m <- bootstrap_VaVp_weight(data = weight_data %>% filter(sex == "M"), n_groups = n_groups, n_iter = n_iter, mean = weight_meanM) %>% 
  bind_rows() %>% filter(is.na(var_sire) == F)
boot_sum_M <- boot_summary(boot_weight_m) %>%
  mutate(sex = "M")

### Plot Weight results ----

weight_res_comb <- rbind(fem_weight_forplot, mal_weight_forplot) # means from REML

boot_sum_both <- rbind(boot_sum_F,boot_sum_M) # confidence intervals from bootstrap
boot_sum_both$var_comp <- factor(boot_sum_both$var_comp, levels = c('VP', 'VA', 'h2', 'I')) #change order of VA and VP
levels(boot_sum_both$var_comp) <- c('V[Phenotypic]', 'V[Additive]', 'h^2', 'I[A]') #change labels for VA and VP

Vars_weight <- ggplot(data = boot_sum_both %>% filter(var_comp != 'h^2' & var_comp != 'I[A]'), aes(x = sex, y = mean)) +
  geom_linerange(aes(ymax = mean + sd, ymin = mean - sd),
                 size = 1, lty = 'solid') +
  geom_point(data = weight_res_comb %>% filter(var_comp != 'h^2' & var_comp != 'I[A]'), aes(x = sex, y = est, fill = sex), size = 3, shape = 21) +
  facet_wrap( ~ var_comp, labeller = label_parsed, scales = "fixed") +
  labs(y = "Variance of body mass", x = "") +
  scale_x_discrete(labels = c('F' = "Female", 'M' = "Male"))+
  scale_fill_manual(labels = c('F' = "Female", 'M' = "Male"), values = c('F' = 'black', 'M' = 'white')) +
  theme_bw() +
  theme(text = element_text(size = 18), #axis.ticks.x = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none')
Vars_weight

Herit_weight <- ggplot(data = boot_sum_both %>% filter(var_comp == 'h^2'), aes(x = sex, y = mean)) +
  geom_linerange(aes(ymax = mean + sd, ymin = mean - sd),
                 size = 1, lty = 'solid') +
  geom_point(data = weight_res_comb %>% filter(var_comp == 'h^2'), aes(x = sex, y = est, fill = sex), size = 3, shape = 21) +
  facet_wrap( ~ var_comp, labeller = label_parsed, scales = "fixed") +
  labs(y = "Heritability of body mass", x = "") +
  scale_x_discrete(labels = c('F' = "Female", 'M' = "Male")) +
  scale_fill_manual(values = c('F' = 'black', 'M' = 'white')) +
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 18), #axis.ticks.x = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none')
Herit_weight

Evolv_weight <- ggplot(data = boot_sum_both %>% filter(var_comp == 'I[A]'), aes(x = sex, y = mean)) +
  geom_linerange(aes(ymax = mean + sd, ymin = mean - sd),
                 size = 1, lty = 'solid') +
  geom_point(data = weight_res_comb %>% filter(var_comp == 'I[A]'), aes(x = sex, y = est, fill = sex), size = 3, shape = 21) +
  facet_wrap( ~ var_comp, labeller = label_parsed, scales = "fixed") +
  labs(y = "Evolvability of body mass", x = "") +
  scale_x_discrete(labels = c('F' = "Female", 'M' = "Male")) +
  scale_fill_manual(values = c('F' = 'black', 'M' = 'white')) +
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(0, 1)) +
  theme_bw() +
  theme(text = element_text(size = 18), #axis.ticks.x = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none')
Evolv_weight

ggarrange(Vars_weight, Herit_weight, Evolv_weight, ncol = 3, widths = c(1, 0.5, 0.5))

### Plot of variation among families ----

resim_weight_M <- REsim(weight_mixedM) %>% mutate(treat = 'm') %>% #male dataset
  mutate(addIntercept = mean + fixef(weight_mixedM))
resim_weight_F <- REsim(weight_mixedF) %>% filter(groupFctr == 'sire') %>% 
  mutate(treat = 'f') %>% #female dataset
  mutate(addIntercept = mean + fixef(weight_mixedF)) %>% 
  arrange(addIntercept) %>%
  mutate(groupID2 = factor(groupID, levels = groupID)) 
str(resim_weight_F)

ggplot(data = resim_weight_F, aes(x = groupID2, y = addIntercept)) +
  # geom_point(aes(fill = treat), size = 3, shape = 21) +  #For means, without sd
  # geom_point(data = resim_weight_M, aes(x = groupID, y = addIntercept, fill = treat), size = 3, shape = 21) +
  geom_pointrange(aes(ymin = addIntercept - sd, ymax = addIntercept + sd, fill = treat), #For means with sd
                  fatten = 2, size = 1, shape = 21) +
  geom_pointrange(data = resim_weight_M, aes(x = groupID, ymin = addIntercept - sd, ymax = addIntercept + sd, fill = treat),
                  fatten = 2, size = 1, shape = 21) +
  lims(y = c(9,12.8)) +
  labs(x = "Rank order of families", y = "Body mass (mg)", fill = "") +
  scale_fill_manual(labels = c('m' = "Male", 'f' = "Female"), values = c('f' = 'black', 'm' = 'white')) +
  theme_classic() +
  theme(text = element_text(size = 18), axis.ticks.x = element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'bottom')


#######################-
# Thorax Width ----

#######################-
cor.test(width_data$weight, width_data$th_width)
cor.test(width_data_F$weight, width_data_F$th_width)
cor.test(width_data_M$weight, width_data_M$th_width)
qplot(data = width_data, x = weight, y = th_width, color = sex) +
  geom_smooth(method = 'lm')


width_data <- weight_data %>%
  filter(is.na(th_width) == F )

width_data %>% group_by(sire,dam) %>% summarise(
  noff = n()) %>% summarise(
    # range = range(noff),
    nDam = n(),
    mean = mean(noff)
  ) %>% View()

heritability.fs(width_data$th_width, width_data$sire, width_data$dam)

### 2-way nested ANOVA ----
width_anova1 <- lm(th_width ~ sex + sire + sire/dam, data = width_data)
summary(width_anova1)
Anova(width_anova1, type = 3)
anova(width_anova1)
plot(width_anova1)

width_anovaF <- lm(th_width ~ sire + sire/dam, data = width_data %>% filter(sex == "F"))
anova(width_anovaF)
width_data_F <- width_data %>% filter(sex == "F")
heritability.fs(width_data_F$th_width, width_data_F$sire, width_data_F$dam)

width_anovaM <- lm(th_width ~ sire, data = width_data %>% filter(sex == "M"))
anova(width_anovaM)
width_data_M <- width_data %>% filter(sex == "M")
heritability.phs(width_data_M$th_width, width_data_M$sire)

### Random effects Models ----
width_mixed1 <- lmer(th_width ~ sex + (1 | sire/dam), data = width_data)
summary(width_mixed1)
# confint(width_mixed1, method = "boot", oldNames = F)
confint(width_mixed1, method = "profile", oldNames = F)
rand(width_mixed1)

# Females only 
width_meanF <- width_data %>% filter(sex == "F") %>% summarise(mean = mean(th_width)) %>% as.numeric()
width_mixedF <- lmer(th_width ~ (1 | sire/dam), data = width_data %>% filter(sex == "F"))
summary(width_mixedF)
rand(width_mixedF)
fem_width_res <- var_from_model_fs(width_mixedF, width_meanF) ##Results from REML
fem_width_forplot <- fem_width_res$VarCompsLong %>%
  mutate(sex = "F")

# Males only
width_meanM <- width_data %>% filter(sex == "M") %>% summarise(mean = mean(th_width)) %>% as.numeric()
width_mixedM <- lmer(th_width ~ (1 | sire), data = width_data %>% filter(sex == "M"))
summary(width_mixedM)
rand(width_mixedM)
mal_width_res <- var_from_model_phs(width_mixedM, width_meanM) ##Results from REML
mal_width_forplot <- mal_width_res$VarCompsLong %>%
  mutate(sex = "M")

### Width - Bootstrap confidence intervals ----

#Bootstrap sd - females only
boot_width_f <- bootstrap_VaVp_width_fs(data = width_data %>% filter(sex == "F"), n_groups = n_groups, n_iter = n_iter, mean = width_meanF) %>% 
  bind_rows() %>% filter(is.na(var_sire) == F)
boot_sum_width_F <- boot_summary(boot_width_f) %>%
  mutate(sex = "F")

#Bootstrap sd - males only
boot_width_m <- bootstrap_VaVp_width(data = width_data %>% filter(sex == "M"), n_groups = n_groups, n_iter = n_iter, mean = width_meanM) %>% 
  bind_rows() %>% filter(is.na(var_sire) == F)
boot_sum_width_M <- boot_summary(boot_width_m) %>%
  mutate(sex = "M")

### Plot Width results ----

width_res_comb <- rbind(fem_width_forplot, mal_width_forplot) # means from REML

boot_sum_width_both <- rbind(boot_sum_width_F,boot_sum_width_M) # confidence intervals from bootstrap
boot_sum_width_both$var_comp <- factor(boot_sum_width_both$var_comp, levels = c('VP', 'VA', 'h2', 'I')) #change order of VA and VP
levels(boot_sum_width_both$var_comp) <- c('V[Phenotypic]', 'V[Additive]', 'h^2', 'I[A]') #change labels for VA and VP

Vars_width <- ggplot(data = boot_sum_width_both %>% filter(var_comp != 'h^2' & var_comp != 'I[A]'),
                     aes(x = sex, y = mean)) +
  geom_linerange(aes(ymax = mean + sd, ymin = mean - sd),
                 size = 1, lty = 'solid') +
  geom_point(data = width_res_comb %>% filter(var_comp != 'h^2' & var_comp != 'I[A]'), aes(x = sex, y = est, fill = sex), size = 3, shape = 21) +
  facet_wrap( ~ var_comp, labeller = label_parsed, scales = "fixed") +
  labs(y = "Variance of thorax width", x = "") +
  scale_x_discrete(labels = c('F' = "Female", 'M' = "Male"))+
  scale_fill_manual(labels = c('F' = "Female", 'M' = "Male"), values = c('F' = 'black', 'M' = 'white')) +
  theme_bw() +
  theme(text = element_text(size = 18), #axis.ticks.x = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none')
Vars_width

Herit_width <- ggplot(data = boot_sum_width_both %>% filter(var_comp == 'h^2'), aes(x = sex, y = mean)) +
  geom_linerange(aes(ymax = mean + sd, ymin = mean - sd),
                 size = 1, lty = 'solid') +
  geom_point(data = width_res_comb %>% filter(var_comp == 'h^2'), aes(x = sex, y = est, fill = sex), size = 3, shape = 21) +
  facet_wrap( ~ var_comp, labeller = label_parsed, scales = "fixed") +
  labs(y = "Heritability of thorax width", x = "") +
  scale_x_discrete(labels = c('F' = "Female", 'M' = "Male")) +
  scale_fill_manual(values = c('F' = 'black', 'M' = 'white')) +
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(-0.03, 1)) +
  theme_bw() +
  theme(text = element_text(size = 18), #axis.ticks.x = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none')
Herit_width

Evolv_width <- ggplot(data = boot_sum_width_both %>% filter(var_comp == 'I[A]'), aes(x = sex, y = mean)) +
  geom_linerange(aes(ymax = mean + sd, ymin = mean - sd),
                 size = 1, lty = 'solid') +
  geom_point(data = width_res_comb %>% filter(var_comp == 'I[A]'), aes(x = sex, y = est, fill = sex), size = 3, shape = 21) +
  facet_wrap( ~ var_comp, labeller = label_parsed, scales = "fixed") +
  labs(y = "Evolvability of thorax width", x = "") +
  scale_x_discrete(labels = c('F' = "Female", 'M' = "Male")) +
  scale_fill_manual(values = c('F' = 'black', 'M' = 'white')) +
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(-0.03, 1)) +
  theme_bw() +
  theme(text = element_text(size = 18), #axis.ticks.x = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none')
Evolv_width

ggarrange(Vars_width, Herit_width, Evolv_width, ncol = 3, widths = c(1, 0.5, 0.5))

### Plot of variation among families ----

resim_width_M <- REsim(width_mixedM) %>% mutate(treat = 'm') %>% #male dataset
  mutate(addIntercept = mean + fixef(width_mixedM))
resim_width_F <- REsim(width_mixedF) %>% filter(groupFctr == 'sire') %>% 
  mutate(treat = 'f') %>% #female dataset
  mutate(addIntercept = mean + fixef(width_mixedF)) %>% 
  arrange(addIntercept) %>%
  mutate(groupID2 = factor(groupID, levels = groupID)) 
str(resim_width_F)

ggplot(data = resim_width_F, aes(x = groupID2, y = addIntercept)) +
  # geom_point(aes(fill = treat), size = 3, shape = 21) +  #For means, without sd
  # geom_point(data = resim_width_M, aes(x = groupID, y = addIntercept, fill = treat), size = 3, shape = 21) +
  geom_pointrange(aes(ymin = addIntercept - sd, ymax = addIntercept + sd, fill = treat), #For means with sd
                  fatten = 2, size = 1, shape = 21) +
  geom_pointrange(data = resim_width_M, aes(x = groupID, ymin = addIntercept - sd, ymax = addIntercept + sd, fill = treat),
                  fatten = 2, size = 1, shape = 21) +
  lims(y = c(2,2.62)) +
  labs(x = "Rank order of families", y = "Thorax width (mm)", fill = "") +
  scale_fill_manual(labels = c('m' = "Male", 'f' = "Female"), values = c('f' = 'black', 'm' = 'white')) +
  theme_classic() +
  theme(text = element_text(size = 18), axis.ticks.x = element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'bottom')

#######################-
## Functions ##----
#######################-

### Functions to extract var comps ----
var_from_model_fs <- function(model, mean) {
  model_results <- as.data.frame(VarCorr(model))
  REML.results <- data.frame(var_sire = NA, var_dam = NA, var_e = NA, V_A = NA, V_P = NA, h2 = NA, I = NA) #blank data.frame to store results
  REML.results.all <- REML.results %>% add_row(var_sire = c(model_results[2,'vcov']), #store results and calculate VA & VP & h2
                                               var_dam = c(model_results[1,'vcov']),
                                               var_e = c(model_results[3,'vcov']), 
                                               V_A = c(4 * var_sire),
                                               V_P = c(var_sire + var_dam + var_e),
                                               h2 = c(V_A / V_P),
                                               I = c(V_A / mean^2)) %>%
    filter(is.na(var_sire) == F)
    
  REML.results2 <- REML.results.all %>%
    dplyr::select(c(V_A:I)) %>%
    pivot_longer(cols = c(1:4), names_to = "var_comp", values_to = "est")
  REML.results2$var_comp <- factor(REML.results2$var_comp, levels = c('V_P', 'V_A', 'h2', 'I')) #change order of VA and VP
  levels(REML.results2$var_comp) <- c('V[Phenotypic]', 'V[Additive]', 'h^2', 'I[A]') #change labels for VA and VP
  
  output = list(AllVarComps = REML.results.all, VarCompsLong = REML.results2)
  
  return(output)
}

var_from_model_phs <- function(model, mean) {
  model_results <- as.data.frame(VarCorr(model))
  REML.results <- data.frame(var_sire = NA, var_e = NA, V_A = NA, V_P = NA, h2 = NA, I = NA) #blank data.frame to store results
  REML.results.all <- REML.results %>% add_row(var_sire = model_results[1,'vcov'], #store results and calculate VA & VP & h2
                                               var_e = model_results[2,'vcov'],
                                               V_A = 4 * var_sire,
                                               V_P = var_sire + var_e,
                                               h2 = V_A / V_P,
                                               I = V_A / mean^2) %>%
    filter(is.na(var_sire) == F)
  
  REML.results2 <- REML.results.all %>%
    dplyr::select(c(V_A:I)) %>%
    pivot_longer(cols = c(1:4), names_to = "var_comp", values_to = "est")
  REML.results2$var_comp <- factor(REML.results2$var_comp, levels = c('V_P', 'V_A', 'h2', 'I')) #change order of VA and VP
  levels(REML.results2$var_comp) <- c('V[Phenotypic]', 'V[Additive]', 'h^2', 'I[A]') #change labels for VA and VP
  
  output = list(AllVarComps = REML.results.all, VarCompsLong = REML.results2)
  
  return(output)
}

### Function to bootsrap confidence interals for weight  ----
VarComps_of_samples_weight <- function(data, n_groups, mean) {
  varComps <- data.frame(var_sire = NA, var_e = NA, V_A = NA, V_P = NA, h2 = NA, I = NA) #blank data.frame to store results
  
  sampled <- data %>% group_by(sire) %>% sample_n_groups(n_groups) #new data frame with sample of original
  
  mixedsamp <- lmer(weight ~ (1 | unique_id), data = sampled) #Random effects model
  vcompsamp <- as.data.frame(VarCorr(mixedsamp)) #extract variance components
  
  varComps <- varComps %>% add_row(var_sire = vcompsamp[1,'vcov'], #store results and calculate VA & VP & h2
                                   var_e = vcompsamp[2,'vcov'], 
                                   V_A = 4 * var_sire, 
                                   V_P = var_sire + var_e,
                                   h2 = V_A / V_P,
                                   I = V_A / mean^2) 
}
### Function to do boostrapping procedure (combines previous functions)
bootstrap_VaVp_weight <- function(data, n_groups, n_iter, mean){
  
  datalist = vector("list", length = n_iter) #list to store results
  
  #Loop to run multiple times
  for (i in 1:n_iter) {
    VarComps2 <- VarComps_of_samples_weight(data = data, n_groups = n_groups, mean = mean)
    datalist[[i]] <- VarComps2
  }
  datalist
}

### Function to bootstrap confidence intervals for weight, full-sib design  ----
VarComps_of_samples_weight_fs <- function(data, n_groups, mean) {
  varComps <- data.frame(var_sire = NA, var_dam = NA, var_e = NA, V_A = NA, V_P = NA, h2 = NA, I = NA) #blank data.frame to store results
  
  sampled <- data %>% group_by(sire) %>% sample_n_groups(n_groups) #new data frame with sample of original
  
  mixedsamp <- lmer(weight ~ (1 | unique_id/dam), data = sampled) #Random effects model
  vcompsamp <- as.data.frame(VarCorr(mixedsamp)) #extract variance components
  
  varComps <- varComps %>% add_row(var_sire = c(vcompsamp[2,'vcov']), #store results and calculate VA & VP & h2
                                   var_dam = c(vcompsamp[1,'vcov']),
                                   var_e = c(vcompsamp[3,'vcov']), 
                                   V_A = c(4 * var_sire),
                                   V_P = c(var_sire + var_dam + var_e),
                                   h2 = c(V_A / V_P),
                                   I = V_A / mean^2) 
}
### Function to do boostrapping procedure (combines previous functions) for weight, full-sib design
bootstrap_VaVp_weight_fs <- function(data, n_groups, n_iter, mean){
  
  datalist = vector("list", length = n_iter) #list to store results
  
  #Loop to run multiple times
  for (i in 1:n_iter) {
    VarComps2 <- VarComps_of_samples_weight_fs(data = data, n_groups = n_groups, mean = mean)
    datalist[[i]] <- VarComps2
  }
  datalist
}

### Function to bootsrap confidence interals for width - phs  ----
VarComps_of_samples_width <- function(data, n_groups, mean) {
  varComps <- data.frame(var_sire = NA, var_e = NA, V_A = NA, V_P = NA, h2 = NA, I = NA) #blank data.frame to store results
  
  sampled <- data %>% group_by(sire) %>% sample_n_groups(n_groups) #new data frame with sample of original
  
  mixedsamp <- lmer(th_width ~ (1 | unique_id), data = sampled) #Random effects model
  vcompsamp <- as.data.frame(VarCorr(mixedsamp)) #extract variance components
  
  varComps <- varComps %>% add_row(var_sire = vcompsamp[1,'vcov'], #store results and calculate VA & VP & h2
                                   var_e = vcompsamp[2,'vcov'], 
                                   V_A = 4 * var_sire, 
                                   V_P = var_sire + var_e,
                                   h2 = V_A / V_P,
                                   I = V_A / mean^2) 
}
### Function to do boostrapping procedure (combines previous functions), width - phs
bootstrap_VaVp_width <- function(data, n_groups, n_iter, mean){
  
  datalist = vector("list", length = n_iter) #list to store results
  
  #Loop to run multiple times
  for (i in 1:n_iter) {
    VarComps2 <- VarComps_of_samples_width(data = data, n_groups = n_groups, mean = mean)
    datalist[[i]] <- VarComps2
  }
  datalist
}

### Function to bootstrap confidence intervals for width, full-sib design  ----
VarComps_of_samples_width_fs <- function(data, n_groups, mean) {
  varComps <- data.frame(var_sire = NA, var_dam = NA, var_e = NA, V_A = NA, V_P = NA, h2 = NA, I = NA) #blank data.frame to store results
  
  sampled <- data %>% group_by(sire) %>% sample_n_groups(n_groups) #new data frame with sample of original
  
  mixedsamp <- lmer(th_width ~ (1 | unique_id/dam), data = sampled) #Random effects model
  vcompsamp <- as.data.frame(VarCorr(mixedsamp)) #extract variance components
  
  varComps <- varComps %>% add_row(var_sire = c(vcompsamp[2,'vcov']), #store results and calculate VA & VP & h2
                                   var_dam = c(vcompsamp[1,'vcov']),
                                   var_e = c(vcompsamp[3,'vcov']), 
                                   V_A = c(4 * var_sire),
                                   V_P = c(var_sire + var_dam + var_e),
                                   h2 = c(V_A / V_P),
                                   I = V_A / mean^2) 
}
### Function to do boostrapping procedure (combines previous functions) for width, full-sib design
bootstrap_VaVp_width_fs <- function(data, n_groups, n_iter, mean){
  
  datalist = vector("list", length = n_iter) #list to store results
  
  #Loop to run multiple times
  for (i in 1:n_iter) {
    VarComps2 <- VarComps_of_samples_width_fs(data = data, n_groups = n_groups, mean = mean)
    datalist[[i]] <- VarComps2
  }
  datalist
}

### Summarize outputs from bootstrapping ----
boot_summary <- function(bootout){
  summ_A <- bootout %>% summarize(  #summarize results for VA
    mean = mean(V_A),
    sd = sd(V_A),
    se = sd(V_A)/sqrt(n()),
    CI.low = quantile(V_A, 0.025),
    CI.high = quantile(V_A, 0.975)) %>%
    mutate(var_comp = c('VA'))

  summ_P <- bootout %>% summarize(  #summarize results for VP
    mean = mean(V_P),
    sd = sd(V_P),
    se = sd(V_P)/sqrt(n()),
    CI.low = quantile(V_P, 0.025),
    CI.high = quantile(V_P, 0.975)) %>%
    mutate(var_comp = c('VP'))
  
  summ_h <- bootout %>% summarize(  #summarize results for h2
    mean = mean(h2),
    sd = sd(h2),
    se = sd(h2)/sqrt(n()),
    CI.low = quantile(h2, 0.025),
    CI.high = quantile(h2, 0.975)) %>%
    mutate(var_comp = c('h2'))
  
  summ_I <- bootout %>% summarize(  #summarize results for h2
    mean = mean(I),
    sd = sd(I),
    se = sd(I)/sqrt(n()),
    CI.low = quantile(I, 0.025),
    CI.high = quantile(I, 0.975)) %>%
    mutate(var_comp = c('I'))
  
  rbind(summ_A,summ_P,summ_h, summ_I)  
}







# library(fullfact)
# fullfact_weight1 <- observLmer(observ = weight_data, dam = "dam", sire = "sire", response = "weight", ml = F)
# fullfact_weight1




