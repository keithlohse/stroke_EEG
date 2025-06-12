# Title: Stroke and Aging EEG Study - Secondary Data Analysis
# Author: Keith Lohse, PhD, PStat
# Date: 2025-06-12


# Import data ----
library(tidyverse); library(patchwork)
library(psych); library(car); library(lme4); library(lmerTest)
library(factoextra); library(corrplot); library(caTools);
library(akima); library(viridis)

getwd()
setwd("C:/Users/lohse/Box/Lohse Brain Age EEG/data/")

list.files()

# Import the Cramer stroke data ----
STROKE <- read.csv("./data_CRAMER_IMPUTED_SUBSET.csv", header=TRUE, stringsAsFactors = TRUE)%>%
  mutate(dataset=c("cramer"), group=c("stroke")) %>% select(-X.2, -X.1, -X)
head(STROKE)

STROKE <- STROKE %>% group_by(subID) %>%
  fill(EEG_Directory_Name:group, .direction="updown")
summary((STROKE%>% group_by(subID) %>% slice(1))$age)

# 1. Descriptive Statistics and Explanation of Spec Param ----------------------
head(STROKE)
describe(STROKE %>% group_by(subID) %>% slice(1))
summary((STROKE %>% group_by(subID) %>% slice(1))$gender)
summary((STROKE %>% group_by(subID) %>% slice(1))$lesion_site)
summary((STROKE %>% group_by(subID) %>% slice(1))$lesion_hemisphere)
summary((STROKE %>% group_by(subID) %>% slice(1))$dom_hemi_affected)
summary((STROKE %>% group_by(subID) %>% slice(1))$ethnicity_race)

summary(STROKE %>% group_by(subID) %>% summarise(missing = sum(is.na(r.2)==FALSE)))
summary(STROKE$r.2)
summary(STROKE$Offset)
summary(STROKE$Exponent)

summary(STROKE %>% group_by(subID) %>% slice(1) %>% select(age, gender, days_to_enrollment, lesion_hemisphere,
                                                           dom_hemi_affected, lesion_volume))

  

# Import the combined control data ----
CONTROL <- read.csv("data_IMPUTED_CONTROLS.csv", stringsAsFactors = TRUE) %>%
  select(-X) %>% mutate(group=c("control"))
head(CONTROL)

CONTROL <- CONTROL %>% group_by(subID) %>%
  fill(age:group, .direction="updown")
summary(CONTROL$dataset)

EULER <- CONTROL %>% filter(dataset=="euler")
PATHANIA <- CONTROL %>% filter(dataset=="pathania")
MPI <- CONTROL %>% filter(dataset=="mpi")

summary(EULER %>% group_by(subID) %>% summarise(complete = sum(is.na(Error)==FALSE)))
summary(EULER$r.2) # why are these values missing for euler data?
# r^2 values were named R.2 rather than r.2 and not merged : /
list.files()
summary(read.csv("euler_MERGED.csv", stringsAsFactors = TRUE) %>% select(R.2))
summary(EULER$Offset)
summary(EULER$Exponent)

summary(EULER %>% group_by(subID) %>% slice(1) %>% select(age, sex))

summary(PATHANIA %>% group_by(subID) %>% summarise(complete = sum(is.na(Error)==FALSE)))
summary(PATHANIA$r.2) 
summary(PATHANIA$Offset)
summary(PATHANIA$Exponent)

summary(PATHANIA %>% group_by(subID) %>% slice(1) %>% select(age, sex))

summary(MPI %>% group_by(subID) %>% summarise(missing = sum(is.na(Error)==FALSE)))
summary(MPI$r.2) 
summary(MPI$Offset)
summary(MPI$Exponent)

summary(MPI %>% group_by(subID) %>% slice(1) %>% select(age, sex))




# Combined imputed data for plotting ----
head(STROKE)
head(CONTROL)


COMB <- rbind(STROKE %>% select(subID, Channels, r.2, dataset,
                                Imp_Off, Imp_Exp, CF1, PW1, BW1, CF2, PW2, BW2,
                                CF3, PW3, BW3, CF4, PW4, BW4, CF5, PW5, BW5, 
                                CF6, PW6, BW6, age, sex, group), 
              CONTROL%>% select(subID, Channels, r.2, dataset,
                                Imp_Off, Imp_Exp, CF1, PW1, BW1, CF2, PW2, BW2,
                                CF3, PW3, BW3, CF4, PW4, BW4, CF5, PW5, BW5, 
                                CF6, PW6, BW6, age, sex, group)
              ) %>%
  rename(Offset = "Imp_Off", Exponent="Imp_Exp")
head(COMB)


# Scalp Topography Figure ------------------------------------------------------
summary(COMB$Channels)
summary(COMB$age)
COMB$group <- factor(COMB$group)
summary(COMB$group)
COMB$age_group <- factor(ifelse(COMB$age < 60, "YA", "OA"))
xtabs(~COMB$group+COMB$age_group)

# mean exponent plot data
df_plot <- COMB %>%
  mutate(group = fct_recode(group, Control="control", Stroke="stroke")) %>%
  group_by(group, age_group, Channels) %>%
  summarise(mean_exponent = mean(Exponent, na.rm=TRUE)) %>%
  mutate(Group = factor(paste(group, "-", age_group, sep="")), 
         Group = fct_relevel(Group, "Control-YA", "Control-OA", "Stroke-YA", "Stroke-OA"))

levels(df_plot$Group)

# EEG Layout
electrode_layout <- tibble::tribble(
  ~Channels, ~x,   ~y,
  "Fp1",	-2.0,	3.5,	
  "Fp2", 	2.0,	3.5,
  "AF7", -2.5,	3.0,
  "AF8", 	2.5, 	3.0,
  "AF3",	-1.5, 3.0,
  "AF4",	1.5,	3.0,
  "F7",	-3.0,	2.0,	
  "F8",	3.0,	2.0,
  "F5",	-2.0,	2.0,	
  "F6",	2.0,	2.0,
  "F3",	-1.0,	2.0,	
  "F4",	1.0,	2.0,
  "F1",	-0.5,	2.0,	
  "F2",	0.5,	2.0,
  "Fz",	0.0,	2.0,	
  "FC5",	-2.0,	1.0,
  "FC6",	2.0,	1.0,	
  "FC3",	-1.0,	1.0,
  "FC4",	1.0,	1.0,	
  "FC1",	-0.5,	1.0,
  "FC2",	0.5,	1.0,	
  "FCz",	0.0,	1.0,
  "T7",	-3.0,	0.0,	
  "T8",	3.0,	0.0,
  "C5",	-2.0,	0.0,	
  "C6",	2.0,	0.0,
  "C3", -1.0,	0.0,	
  "C4",	1.0,	0.0,
  "C1",	-0.5,	0.0,	
  "C2",	0.5,	0.0,
  "Cz",	0.0,	0.0,	
  "CP5",	-2.0,	-1.0,
  "CP6",	2.0,	-1.0,	
  "CP3",	-1.0,	-1.0,
  "CP4",	1.0,	-1.0,	
  "CP1",	-0.5,	-1.0,
  "CP2",	0.5,	-1.0,	
  "CPz",	0.0,	-1.0,
  "P7",	-3.0,	-2.0,	
  "P8",	3.0,	-2.0,
  "P5",	-2.0,	-2.0,	
  "P6",	2.0,	-2.0,
  "P3",	-1.0,	-2.0,
  "P4",	1.0,	-2.0,
  "P1",	-0.5,	-2.0,	
  "P2",	0.5,	-2.0,
  "Pz",	0.0,	-2.0,	
  "PO7",	-2.5,	-3.0,
  "PO8",	2.5,	-3.0,
  "PO3",	-1.5,	-3.0,
  "PO4",	1.5,	-3.0,	
  "POz",	0.0,	-3.0,
  "O1",	-1.0,	-3.5,	
  "O2",	1.0,	-3.5,
  "Oz",	0.0,	-3.5,
)

df_plot <- df_plot %>%
  left_join(electrode_layout, by = "Channels") %>%
  mutate(Channels = factor(Channels))

levels(df_plot$Channels)


# Function to interpolate for one group
interpolate_group <- function(df_group, grid_res = 100) {
  interp_res <- with(df_group, akima::interp(
    x = x,
    y = y,
    z = mean_exponent,
    duplicate = "mean",
    nx = grid_res,
    ny = grid_res
  ))
  expand.grid(x = interp_res$x, y = interp_res$y) %>%
    mutate(z = as.vector(interp_res$z)) %>%
    filter(!is.na(z)) %>%
    mutate(Group = unique(df_group$Group))
}

# Apply to each group
interp_all <- df_plot %>%
  group_split(Group) %>%
  map_df(interpolate_group) 


ggplot(data = interp_all) +
  geom_raster(aes(x = x, y = y, fill = z), interpolate = TRUE) +
  geom_point(data = df_plot, aes(x = x, y = y), color = "black") +
  geom_text(data = df_plot, aes(x = x, y = y, label = Channels), 
            size = 2, vjust = -0.8, check_overlap = TRUE) +
  scale_fill_viridis_c(name = "Exponent", option = "plasma") +
  coord_fixed() +
  facet_wrap(~Group, ncol=2) +
  theme_minimal() +
  labs(title = "EEG Exponent Scalp Map", x = NULL, y = NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )


# We ultimately use Asher's figure from MatLab to show scalp topographies
ggsave(
  filename="./scalp_map.jpeg",
  plot = last_plot(),
  path = "C:/Users/lohse/Box/Lohse Brain Age EEG/reports/",
  width = 6,
  height =6,
  units = "in",
  dpi = 300
)


# Prepping data for Analysis ---------------------------------------------------
# Select subset of electrodes to plot:
plot_subset <- c("F3", "Fz", "F4", "C3", "Cz", "C4", "P3", "Pz", "P4", "O1", "Oz", "O2")


# Select subset of electrode to analyze (droping temporal and FP):
subset <- c("F3", "Fz", "F4", "C3", "Cz", "C4", "P3", "Pz", "P4", "O1", 
            "Oz", "O2", "CP1", "CP2", "CP5",  "CP6", "F7", "F8", "FC1", "FC2",
            "FC5", "FC6", "P7", "P8")

summary(COMB$Channels)
gsub("[[:digit:]]","", COMB$Channels) # we can extract the region of the electrode using regex

mean_age<-mean(COMB %>% group_by(subID) %>% slice(1) %>% pull(age))

COMB <- COMB %>% mutate(
  region = factor(gsub("[[:digit:]]","", Channels)),
  group = factor(group),
  age.c = age/10-mean_age/10,
  Channels = fct_relevel(Channels, "F3", "Fz", "F4", 
                         "C3", "Cz", "C4", 
                         "P3", "Pz", "P4", 
                         "O1", "Oz", "O2"))
summary(COMB$age.c)
summary(COMB$Channels)
summary(COMB$region)

str_sub(COMB$region, start=0, end=1)

# Drop Fp1 and Fp2, T7 and T8
summary(COMB$Channels)
COMB <- COMB %>% 
  filter(Channels %in% subset) %>%
  mutate(Channels = factor(Channels),
         # remove second letters so everthing is F, C, P or O
         # this means we operationally calling CP, "C", and FC, "F"
         region = str_sub(region, start=0, end=1), 
         region = factor(region),
         region = fct_relevel(region, "F", "C", "P", "O"),
         group = fct_relevel(group, "control", "stroke"))

summary(COMB$Channels)
summary(COMB$region)

head(COMB)
summary(COMB$r.2[COMB$group=="control"])
summary(COMB$r.2[COMB$group=="stroke"])
summary(COMB$r.2[COMB$region=="F"])
summary(COMB$r.2[COMB$region=="C"])
summary(COMB$r.2[COMB$region=="P"])
summary(COMB$r.2[COMB$region=="O"])

# Set Colorblind friendly pallette 
cbPalette <- c("#999999", "#56B4E9","#E69F00", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
               "#999933", "#882255", "#661100", "#6699CC")


ggplot(data=COMB %>% filter(Channels %in% plot_subset),
       aes(x=Exponent)) +
  geom_density(aes(fill=group), col="black", alpha=0.5) +
  geom_point(aes(y=0, fill=group), col="black", alpha=0.5,shape=21, 
             position=position_jitter(height=0.05))+
  facet_wrap(~Channels, ncol=3) +
  scale_x_continuous(name = "Spectral Slope (Exponent)") +
  scale_y_continuous(name = "Density")+
  theme_bw()+labs(fill="Group")+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text=element_text(size=10, color="black"), 
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size=10, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=10),
        legend.position = "bottom")

ggsave(
  filename="./rainclouds_exponents.jpeg",
  plot = last_plot(),
  path = "C:/Users/lohse/Box/Lohse Brain Age EEG/reports/",
  width = 9,
  height = 9,
  units = "in",
  dpi = 300
)


ggplot(data=COMB %>% filter(Channels %in% plot_subset),
       aes(x=Offset, 
           y=Exponent)) +
  geom_point(aes(fill=group), col="black", alpha=0.5,shape=21)+
  stat_smooth(aes(col=group), method="lm", se=FALSE)+
  facet_wrap(~Channels, ncol=3) +
  scale_x_continuous(name = "Imputed Offset") +
  scale_y_continuous(name = "Imputed Exponent")+
  theme_bw()+labs(fill="Group", col="Group")+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text=element_text(size=10, color="black"), 
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size=10, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=10),
        legend.position = "bottom")



ggsave(
  filename="./exponents_to_offsets.jpeg",
  plot = last_plot(),
  path = "C:/Users/lohse/Box/Lohse Brain Age EEG/reports/",
  width = 9,
  height = 9,
  units = "in",
  dpi = 300
)







# 2. Effects of Age and Stroke Incidence ---------------------------------------
# 2.1 Spectral slope (exponent) ----


# Linear model of age in each region -------------------------------------------
options(contrasts=c("contr.poly","contr.poly"))
COMB$stroke <- as.numeric(COMB$group) - 1
COMB$control <- (as.numeric(COMB$group) - 2)*(-1)
contrasts(COMB$group)

# Comparison of Random-Effects
sub_mod <- lmer(Exponent~1+(1|subID),
                       data=COMB,
                       REML=FALSE,
                       control=lmerControl(optimizer="bobyqa",
                                           optCtrl=list(maxfun=5e5)))

channel_mod <- lmer(Exponent~1+(1|subID)+(1|Channels),
                    data=COMB,
                    REML=FALSE,
                    control=lmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=5e5)))

full_rand <- lmer(Exponent~1+(1|subID)+(1|Channels)+(1|region),
                   data=COMB,
                   REML=FALSE,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=5e5)))

fixed_region <- lmer(Exponent~1+region+(1|subID)+(1|Channels),
                  data=COMB,
                  REML=FALSE,
                  control=lmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=5e5)))


anova(sub_mod, channel_mod, full_rand, fixed_region)



# Testing Effects of Age ---
age_linear <- lmer(Exponent~
                     # fixed effects
                     1+age.c*region+
                     # random effects
                     (1|subID)+(1|Channels),
                   data=COMB,
                   REML=FALSE,
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=5e5)))


age_quad <- lmer(Exponent~
                   # fixed effects
                   1+age.c*region+I(age.c^2)*region+
                   # random effects
                   (1|subID)+(1|Channels),
                 data=COMB,
                 REML=FALSE,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=5e5)))

age_cube <- lmer(Exponent~
                   # fixed effects
                   1+age.c*region+I(age.c^2)*region+I(age.c^3)*region+
                   # random effects
                   (1|subID)+(1|Channels),
                 data=COMB,
                 REML=FALSE,
                 control=lmerControl(optimizer="bobyqa",
                                     optCtrl=list(maxfun=5e5)))

anova(age_linear, age_quad, age_cube)
# We do find some evidence for a non-linear effect of age. We dug deeper into 
# this, which included fitting one- and two-knot linear spline models. 
# However, this non-linearity appeared to be driven by the relative lack of 
# "middle aged" adults (especially in the control data). Erring on the side of 
# caution, we focus our analyses on the simpler linear effect of age that captures
# the majority of the variance.


# Fitting the Fixed Effects Model ----
exp_mod_linear <- lmer(Exponent ~ 
                         age.c*region*group + 
                         (1|subID) + (1|Channels),
                       data = COMB,
                       REML = FALSE,
                       control = lmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 5e5)))

hist(COMB$age.c)
#confint.merMod(exp_mod_cube)
vif(exp_mod_spline)
vif(exp_mod_linear)

summary(exp_mod_linear)
anova(exp_mod_linear)


sink(file="./exp_mod_linear.txt")
summary(exp_mod_linear)
anova(exp_mod_linear)
sink(file = NULL)



plot(density(rstudent(exp_mod_linear)))
qqnorm(rstudent(exp_mod_linear))
abline(0,1)
qqnorm(scale(ranef(exp_mod_linear)$subID$`(Intercept)`))
abline(0,1)
qqnorm(scale(ranef(exp_mod_linear)$Channel$`(Intercept)`))
abline(0,1)


# Post-Hoc Test
library(emmeans)
emtrends(exp_mod_linear, pairwise ~ group+region, var = "age.c")


# Figure 4 ----
# Fit predicted spline data
PRED<-data.frame(
  expand.grid(group = c(levels(COMB$group)),
              region = c(levels(COMB$region)),
              age = c(18:86))
)

PRED$age.c <- PRED$age/10 - mean(unique(COMB$age)/10)
PRED$linear_exp <- predict(exp_mod_linear, PRED, re.form = NA)


head(COMB)

# Figure 4A Exponent by Age ---------------------
# need to stitch together with Asher's plot outside of R
ggplot(data=COMB,
       aes(x=age, 
           y=Exponent)) +
  geom_point(aes(group=region, col=region, shape=region), alpha=0.5,
             position=position_dodge(width=0.5))+
  #stat_smooth(aes(lty=region), col="black", method="lm", se=FALSE, lwd=0.5)+
  geom_line(data=PRED, aes(y=linear_exp, lty=region), col="black", lwd=0.5)+
  facet_wrap(~group) +
  scale_x_continuous(name = "Age") +
  scale_y_continuous(name = "Exponent")+
  theme_bw()+labs(shape="Region", col="Region", lty="Region")+
  #guides(col="none", shape="none")+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text=element_text(size=10, color="black"), 
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size=10, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=10),
        legend.position = "bottom")


ggsave(
  filename="./Figure4A.png",
  plot = last_plot(),
  path = "C:/Users/lohse/Box/Lohse Brain Age EEG/reports/",
  width = 6,
  height = 3.5,
  units = "in",
  dpi = 300
)





# 3. Narrow band power ---------------------------------------------------------
# 3.1 Number of Peaks ----
head(COMB)
PEAK_NUM <- COMB %>% select(subID, dataset, Channels, region, group, age_group, age, age.c, sex, 
                            CF1, PW1, CF2, PW2, CF3, PW3, CF4, PW4, CF5, 
                            PW5, CF6, PW6) %>%
  pivot_longer(CF1:PW6) %>% 
  mutate_at(c('value'), ~na_if(., 0)) %>% # Euler data has 0 where should be NAs for power
  mutate(variable = substr(name, 1,2),
         number = substr(name, 3,3)) %>%
  select(-name) %>%
  pivot_wider(names_from = variable, values_from=value)


head(PEAK_NUM)
summary(PEAK_NUM$CF)
summary(PEAK_NUM$PW)
summary(PEAK_NUM$group)

# number of peaks per person
PEAK_COUNT <- PEAK_NUM %>% group_by(subID, Channels) %>%
  summarize(dataset=dataset[1],
            group = group[1],
            age = age[1],
            age.c = age.c[1],
            region = region[1],
            sex = sex[1],
            peak_count = sum(is.na(CF)==FALSE))


# Descriptive statistics for peak counts
head(PEAK_COUNT)
PEAK_COUNT %>% 
  group_by(dataset) %>%
  summarize(ave = mean(peak_count),
            med = median(peak_count),
            s = sd(peak_count))


hist(PEAK_COUNT$peak_count)
mean(PEAK_COUNT$peak_count)
var(PEAK_COUNT$peak_count)

options(contrasts=c("contr.poly","contr.poly"))
peak_count_linear <- lmer(sqrt(peak_count+1)~
                             # fixed effects
                             age.c*group*region+
                             # random effects
                             (1|subID)+(1|Channels),
                          data=PEAK_COUNT, REML=FALSE)

anova(peak_count_linear)
summary(peak_count_linear)

sink(file="./peak_count_linear.txt")
summary(peak_count_linear)
anova(peak_count_linear)
sink(file = NULL)

plot(density(rstudent(peak_count_linear)))
qqnorm(rstudent(peak_count_linear))
abline(0,1)

# Post-Hoc Test (refit model without transformation)
emmeans(peak_count_linear, pairwise ~ group)
emtrends(peak_count_linear, pairwise ~ group, var = "age.c")


# 3.2 Central Frequencies ----
head(PEAK_NUM)

PEAK_BAND <- PEAK_NUM %>% 
  filter(is.na(CF)==FALSE) %>%
  mutate(band = factor(cut(CF, 
                           breaks=c(0,4,8,12,30),
                           labels=c("delta", "theta", "alpha", "beta"))),
         band=fct_relevel(band, "delta", "theta", "alpha", "beta")) %>%
  mutate(region = factor(substr(Channels, start=0, stop=1)),
         side = factor(substr(Channels, start=2, stop=3)),
         side = factor(ifelse(side == "z", "z",
                              ifelse(side == "1"|side=="3", "left",
                                     "right"))))



summary(PEAK_NUM$CF)
summary(PEAK_BAND$CF)
summary(PEAK_BAND$age.c)
summary(PEAK_BAND$group)

CF <- ggplot(data=PEAK_BAND,
           aes(x=age, 
               y=CF)) +
  geom_point(aes(col=band, shape=band), alpha=0.3)+
  stat_smooth(aes(lty=band), col="black", lwd=0.5, method="lm", se=FALSE)+
  facet_wrap(~group) +
  scale_x_continuous(name = "Age Group") +
  scale_y_continuous(name = "Central Frequency of Peaks (Hz)", limits=c(0,30))+
  theme_bw()+labs(lty="Band", col="Band", shape="Band")+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text=element_text(size=10, color="black"), 
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size=10, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=10),
        legend.position = "bottom")



head(PEAK_BAND)
summary(PEAK_BAND$band)


options(contrasts=c("contr.poly","contr.poly"))
band_linear <- lmer(CF~
                      # fixed effects
                      band*age.c*group*region+
                      # random effects
                      (1|subID)+(1|Channels),
                    REML=FALSE, 
                    data=PEAK_BAND %>% filter(band != "delta"), 
                    control=lmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=5e5)))
vif(band_linear)
anova(band_linear)
summary(band_linear)

sink(file="./band_linear.txt")
summary(band_linear)
anova(band_linear)
sink(file = NULL)

plot(density(rstudent(band_linear)))
qqnorm(rstudent(band_linear))
abline(0,1)


# Post-Hoc Tests
emtrends(band_linear, pairwise ~ band, var = "age.c")
emmeans(band_linear, pairwise ~ region|band)


# 3.3 Power at Central Frequencies ----
head(PEAK_BAND)
hist(PEAK_BAND$CF)

PW <- ggplot(data=PEAK_BAND, 
             aes(x=age, 
                 y=PW)) +
  geom_point(aes(col=band, shape=band), alpha=0.2)+
  stat_smooth(aes(lty=band), col="black", lwd=0.5, method="lm", se=FALSE)+
  facet_wrap(~group) +
  scale_x_continuous(name = "Age (y)") +
  scale_y_continuous(name = "Power at Central Frequency")+
  theme_bw()+labs(col="Band", shape="Band", lty="Band")+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text.y=element_text(size=10, color="black"), 
        axis.text.x=element_text(size=10, color="black"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size=10, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=10),
        legend.position = "bottom")



library(patchwork)
CF/PW

ggsave(
  filename="./Figure5.jpeg",
  plot = last_plot(),
  path = "C:/Users/lohse/Box/Lohse Brain Age EEG/reports/",
  width = 5,
  height = 6,
  units = "in",
  dpi = 600
)



options(contrasts=c("contr.poly","contr.poly"))
band_power_linear <- lmer(log(PW+1)~
                             # fixed effects
                             band*age.c*group*region+
                             # random effects
                             (1|subID)+(1|Channels),
                           REML=TRUE, 
                           data=PEAK_BAND %>% filter(band!="delta"))

anova(band_power_linear)
summary(band_power_linear)

sink(file="./band_power_linear.txt")
summary(band_power_linear)
anova(band_power_linear)
sink(file = NULL)

#confint.merMod(band_power_linear)
qqnorm(rstudent(band_power_linear))
abline(0,1)




# Post-Hoc Tests (refit model without transformation for post-hoc tests)
emmeans(band_power_linear, pairwise ~ group|band)
emtrends(band_power_linear, pairwise ~ band, var = "age.c")








# 4. Associations within the Stroke Sub-Group ----------------------------------
head(STROKE)
summary(STROKE$lesion_hemisphere)
summary(STROKE$lesion_site)
summary(STROKE$Channels)

write.csv(STROKE, "./data_STROKE_SUBSET.csv")

write.csv(STROKE %>% 
  select(subID, age, Offset, Exponent) %>%
  group_by(subID) %>%
  summarize(age= age[1],
            grand_ave_offset = mean(Offset),
            grand_ave_exponent = mean(Exponent)),
  "./data_STROKE_SUBSET_AVE.csv")
  
summary(STROKE$Channels)

# STROKE <- STROKE %>% 
#   filter(Channels %in% subset) 


STROKE <- STROKE %>% 
  mutate(
  channel_region = factor(substr(Channels, 1, 1)),
  channel_side = substr(Channels, 
                        stringi::stri_length(Channels), 
                        stringi::stri_length(Channels))) %>%
  mutate(channel_side = factor(ifelse(channel_side=="z", "z",
                                      ifelse((as.numeric(channel_side) %% 2)==0, "R",
                                             "L"))),
         channel_side = fct_relevel(channel_side, "L", "z", "R"),
         channel_region = fct_relevel(channel_region, "F", "C", "P", "O", "T")) %>%
  filter(channel_region != "T") %>%
  filter(channel_side != "z") %>%
  filter(Channels != "Fp1") %>%
  filter(Channels != "Fp2") %>%
  mutate(Channels = factor(Channels),
         lesion_hemisphere = factor(lesion_hemisphere),
         channel_side = factor(channel_side),
         channel_region = factor(channel_region),
         contra = factor(ifelse(lesion_hemisphere==channel_side,
                                "ipsilesional",
                                "contralesional")))


summary(STROKE$channel_side)
summary(STROKE$channel_region)
summary(STROKE$contra)
summary(STROKE$Channels)

# 4.1 Spatial effects of lesion location on exponent ----
head(STROKE)

# Missingness
MISSING <- STROKE %>% select(subID, Channels, lesion_hemisphere, channel_region, contra, Exponent) %>%
  mutate(missing = is.na(Exponent) == TRUE) %>% 
  group_by(subID, lesion_hemisphere, channel_region, contra) %>%
  summarize(count_missing = sum(missing))

# Exponent by Lesion Hemisphere
ggplot(data=MISSING,
       aes(x=contra, 
           y=count_missing)) +
  geom_point(aes(col=contra), alpha=0.3,
             position=position_jitterdodge(dodge.width=0.75,
                                           jitter.width=0.1, jitter.height = 0.1))+
  geom_boxplot(aes(fill=contra), alpha=0.3,
               position=position_dodge(width=0.75), outlier.shape=NA)+
  facet_wrap(~channel_region) +
  scale_x_discrete(name = "Location") +
  scale_y_continuous(name = "Amount of Missing Data")+
  theme_bw()+
  #labs(fill="Montage Region", col="Montage Region")+
  #guides(shape="none")+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text=element_text(size=10, color="black"), 
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size=10, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=10, face="bold"),
        legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=10),
        legend.position = "none")

ggsave(
  filename="./missingness_by_region.jpeg",
  plot = last_plot(),
  path = "C:/Users/lohse/Box/Lohse Brain Age EEG/reports/",
  width = 6,
  height = 5.5,
  units = "in",
  dpi = 300
)


# Exponent by Lesion Hemisphere
xtabs(~Channels+channel_side, data=STROKE)
head(STROKE)

STROKE <- STROKE %>% 
  group_by(subID, channel_region) %>%
  mutate(channel_number = substr(Channels, 
                      stringi::stri_length(Channels), 
                      stringi::stri_length(Channels)),
         elect_pair = floor((as.numeric(channel_number)-1)/2),
         elect_pair = paste(subID, channel_region, elect_pair, sep="")) %>%
  arrange(subID, Channels)


ggplot(data=STROKE,
       aes(x=contra, 
           y=Imp_Exp)) +
  geom_line(aes(group=elect_pair), lwd=0.5, col="grey", alpha=0.4,
            position=position_dodge(width=0.1))+
  geom_point(aes(group=elect_pair, col=contra), alpha=0.3,
             position=position_dodge(width=0.1))+
  geom_boxplot(aes(fill=contra), alpha=1.0, outlier.shape=NA, width=0.1,
               position= position_nudge(x=-.15))+
  facet_wrap(~channel_region) +
  scale_x_discrete(name = "Channel Side") +
  scale_y_continuous(name = "Exponent")+
  theme_bw()+
  #labs(fill="Montage Side", col="Montage Side")+
  #guides(shape="none")+
  scale_fill_manual(values=cbPalette[c(2,4)])+
  scale_colour_manual(values=cbPalette[c(2,4)])+
  theme(axis.text=element_text(size=12, color="black"), 
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "none")

ggsave(
  filename="./Figure6.jpeg",
  plot = last_plot(),
  path = "C:/Users/lohse/Box/Lohse Brain Age EEG/reports/",
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)

colnames(STROKE)
options(contrasts=c("contr.poly","contr.poly"))
exp_lesion_mod <- lmer(Imp_Exp~
                         # fixed effects
                         sex+age+days_to_enrollment+lesion_volume+
                         contra*channel_region+
                         # random effects
                         (1+contra|subID)+(1|Channels),
                       REML=TRUE, 
                       data=STROKE %>% 
                         filter(channel_side != "z"))
anova(exp_lesion_mod)
summary(exp_lesion_mod)


sink(file="./STROKE_exp_lesion_mod.txt")
summary(exp_lesion_mod)
anova(exp_lesion_mod)
sink(file = NULL)


# Level 1 residuals
qqnorm(rstudent(exp_lesion_mod))
abline(0,1)

# Level 2 random intercepts and slopes
qqnorm(scale(ranef(exp_lesion_mod)$subID$`(Intercept)`))
qqnorm(scale(ranef(exp_lesion_mod)$subID$contra))
qqnorm(scale(ranef(exp_lesion_mod)$Channels))
abline(0,1)


emmeans::emmeans(exp_lesion_mod, pairwise~contra|channel_region)






# 4.2. Lesion location on narrowband power -----
colnames(STROKE)
PEAK_NUM <- STROKE %>% select(subID, elect_pair, Channels, lesion_site, 
                              lesion_hemisphere, age, days_to_enrollment,
                              sex, lesion_volume, cst_p_injured,
                              channel_region, channel_number, channel_side, contra,
                              CF1, PW1, CF2, PW2, CF3, PW3, CF4, PW4, CF5, 
                            PW5, CF6, PW6) %>%
  pivot_longer(CF1:PW6) %>% 
  mutate_at(c('value'), ~na_if(., 0)) %>% # Euler data has 0 where should be NAs for power
  mutate(variable = substr(name, 1,2),
         number = substr(name, 3,3)) %>%
  select(-name) %>%
  pivot_wider(names_from = variable, values_from=value) %>%
  filter(channel_side != "z") %>%
  filter(channel_region != "T")

summary(PEAK_NUM$channel_side)
summary(PEAK_NUM$channel_region)
summary(PEAK_NUM$elect_pair)

head(PEAK_NUM)
summary(PEAK_NUM$CF)
summary(PEAK_NUM$PW)




# 4.2.1 Central Frequencies ----
head(PEAK_NUM)

PEAK_BAND <- PEAK_NUM %>% 
  filter(is.na(CF)==FALSE) %>%
  mutate(band = factor(cut(CF, 
                           breaks=c(0,4,8,12,30),
                           labels=c("delta", "theta", "alpha", "beta"))),
         band=fct_relevel(band, "delta", "theta", "alpha", "beta"),
         channel_region = factor(channel_region)) 


head(PEAK_BAND)

stroke_CF <- ggplot(data=PEAK_BAND,
       aes(x=contra, 
           y=CF)) +
  geom_point(aes(fill=band),
             col="black", alpha=0.5, shape=21,
             position=position_jitterdodge(dodge.width = 0.75))+
  geom_boxplot(aes(fill=band),
               col="black", alpha=0.5, outlier.shape=NA,
               position=position_dodge(width = 0.75))+
  facet_wrap(~channel_region, ncol=2) +
  scale_x_discrete(name = "Channel Side") +
  scale_y_continuous(name = "Central Frequency of Peaks")+
  theme_bw()+labs(fill="Channel Side", col="Channel Side")+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text=element_text(size=12, color="black"), 
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "none")

# ggsave(
#   filename="./central_freq.jpeg",
#   plot = last_plot(),
#   path = "C:/Users/lohse/Box/Lohse Brain Age EEG/reports/",
#   width = 6,
#   height = 3,
#   units = "in",
#   dpi = 300
# )


head(PEAK_BAND)
# Number of peaks
PEAK_COUNT <- PEAK_BAND %>% group_by(subID, Channels, band) %>%
  summarize(elect_pair = elect_pair [1],
            age = age[1],
            sex = sex[1],
            days_to_enrollment = days_to_enrollment[1],
            lesion_volume = lesion_volume[1],
            cst_p_injured = cst_p_injured[1],
            channel_region = channel_region[1],
            channel_number = channel_number[1],
            contra = contra[1],
            band = band[1],
            peak_count = sum(is.na(CF)==FALSE))


ggplot(data=PEAK_COUNT,
       aes(x=contra, 
           y=peak_count)) +
  geom_point(aes(fill=band),
             col="black", alpha=0.5, shape=21,
             position=position_jitterdodge(dodge.width = 0.75))+
  geom_boxplot(aes(fill=band),
               col="black", alpha=0.5, outlier.shape=NA,
               position=position_dodge(width = 0.75))+
  facet_wrap(~channel_region, ncol=2) +
  scale_x_discrete(name = "Channel Side") +
  scale_y_continuous(name = "# of Peaks")+
  theme_bw()+labs(fill="Channel Side", col="Channel Side")+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text=element_text(size=12, color="black"), 
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "none")


band_count <- lmer(sqrt(peak_count+1)~
                      # fixed effects
                      sex+age+days_to_enrollment+lesion_volume+
                      band*contra*channel_region+
                      # random effects
                      (1|subID), # variance for Channel RE goes to 0
                    REML=FALSE, 
                    data=PEAK_COUNT %>% filter(band != "delta"), 
                    control=lmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=5e5)))
anova(band_count)
summary(band_count)

sink(file="./STROKE_band_count.txt")
summary(band_count)
anova(band_count)
sink(file = NULL)

# Central frequency of peaks
band_linear <- lmer(CF~
                      # fixed effects
                      sex+age+days_to_enrollment+lesion_volume+
                      band*contra*channel_region+
                      # random effects
                      (1+contra|subID), # variance for Channel RE goes to 0
                    REML=FALSE, 
                    data=PEAK_BAND %>% filter(band != "delta"), 
                    control=lmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=5e5)))
anova(band_linear)
summary(band_linear)


sink(file="./STROKE_band_linear.txt")
summary(band_linear)
anova(band_linear)
sink(file = NULL)



# 4.2.2 Power at CF ----
head(PEAK_BAND)

stroke_PW <- ggplot(data=PEAK_BAND,
       aes(x=contra, 
           y=PW)) +
  geom_point(aes(fill=band),
             col="black", alpha=0.5, shape=21,
             position=position_jitterdodge(dodge.width = 0.75))+
  geom_boxplot(aes(fill=band),
               col="black", alpha=0.5, outlier.shape=NA,
               position=position_dodge(width = 0.75))+
  facet_wrap(~channel_region, ncol=2) +
  scale_x_discrete(name = "Channel Side") +
  scale_y_continuous(name = "Power at CF (uV^2)")+
  theme_bw()+labs(fill="Band", col="Band")+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text=element_text(size=12, color="black"), 
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom")

stroke_PW

ggsave(
  filename="./Figure8.jpeg",
  plot = last_plot(),
  path = "C:/Users/lohse/Box/Lohse Brain Age EEG/reports/",
  width = 6,
  height = 6,
  units = "in",
  dpi = 600
)

stroke_CF/stroke_PW


ggsave(
  filename="./stroke_CF_and_PW.jpeg",
  plot = last_plot(),
  path = "C:/Users/lohse/Box/Lohse Brain Age EEG/reports/",
  width = 6,
  height = 10,
  units = "in",
  dpi = 300
)



options(contrasts=c("contr.poly","contr.poly"))
power_band <- lmer(log(PW+1)~
                      # fixed effects
                     sex+age+days_to_enrollment+lesion_volume+
                     band*contra*channel_region+
                      # random effects
                      (1+contra|subID)+(1|Channels),
                    REML=FALSE, 
                   data=PEAK_BAND %>% filter(band != "delta"), 
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=5e5)))
anova(power_band)
summary(power_band)
#confint.merMod(power_band)

sink(file="./STROKE_power_band.txt")
summary(power_band)
anova(power_band)
sink(file = NULL)


# Level 1 residuals
qqnorm(rstudent(power_band))
abline(0,1)

# Level 2 random intercepts and slopes
qqnorm(scale(ranef(power_band)$subID$`(Intercept)`))
qqnorm(scale(ranef(power_band)$Channels))
abline(0,1)


emmeans::emmeans(power_band, pairwise~channel_region|contra)
emmeans::emmeans(power_band, pairwise~band|channel_region)




# 5.0. Relationship of exponent to box and block test ----
colnames(STROKE)
summary(STROKE$age)
summary(STROKE$bbt_affected)

mean(STROKE$age)
mean(STROKE$bbt_affected, na.rm=TRUE)

STROKE <- STROKE %>% 
  ungroup() %>%
  mutate(bbt_affected.c = bbt_affected - mean(bbt_affected, na.rm=TRUE),
         age.c = age - mean(age, na.rm=TRUE),
         age.20 = age - 20, # for testing BBT effects at different ages
         age.40 = age - 40,
         age.60 = age - 60,
         age.80 = age - 80,
         age_band = cut(age, breaks=c(20, 60, 100)),
         lesion_volume.c = lesion_volume - mean(lesion_volume, na.rm=TRUE),
         days_to_enrollment.c = days_to_enrollment - mean(days_to_enrollment, na.rm=TRUE))
  
summary(STROKE$age)
summary(STROKE$age_band)

summary(STROKE$bbt_affected.c)
summary(STROKE$age.c)
summary(STROKE$sex)

hist(STROKE$bbt_affected)
summary(STROKE$bbt_affected)
plot(x=STROKE$bbt_affected, y=STROKE$Exponent)

options(contrasts=c("contr.poly","contr.poly"))
exp_age_linear <- lmer(Imp_Exp~
                     # fixed effects
                     sex+days_to_enrollment.c+lesion_volume.c+
                     bbt_affected.c*age.c*channel_region*contra+
                     # random effects
                     (1+contra|subID)+(1|Channels),
                   REML=FALSE, 
                   data=STROKE, 
                   control=lmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=5e5)))
anova(exp_age_linear)
summary(exp_age_linear)
vif(exp_age_linear)

sink(file="./STROKE_exp_age_linear.txt")
summary(exp_age_linear)
anova(exp_age_linear)
sink(file = NULL)


# Level1 residuals
qqnorm(rstudent(exp_age_linear))
abline(0,1)

# Level 2 random intercepts and slopes
qqnorm(scale(ranef(exp_age_linear)$subID$`(Intercept)`))
qqnorm(scale(ranef(exp_age_linear)$Channels))
abline(0,1)


emmeans::emmeans(exp_age_linear, pairwise~channel_region|contra)
emtrends(exp_age_linear, pairwise ~ channel_region, var = "bbt_affected.c")


# Figure showing model predictions at different ages ---------------------------
exp_age_linear_raw <- lmer(Imp_Exp~
                         # fixed effects
                         days_to_enrollment.c+lesion_volume.c+
                         bbt_affected*age*channel_region+
                         # random effects
                         (1+contra|subID)+(1|Channels),
                       REML=FALSE, 
                       data=STROKE, 
                       control=lmerControl(optimizer="bobyqa",
                                           optCtrl=list(maxfun=5e5)))
anova(exp_age_linear_raw)
summary(exp_age_linear_raw)
vif(exp_age_linear_raw)

hist(STROKE$age)
# making simulated data for the effects of BBT, Age, and Channel Region
PREDS <- data.frame(expand.grid(
  channel_region = c(levels(STROKE$channel_region)), # use all the same channel regions as model
  age = c(20, 40, 60, 80), # select 4 ages we want to highlight
  bbt_affected = c(seq(from=0, to=50, by=5))
  ))
  
  
PREDS <- PREDS %>% 
  mutate(
    days_to_enrollment.c=0,
    lesion_volume.c=0)

PREDS$Exponent <- predict(exp_age_linear_raw, newdata=PREDS, re.form = NA)

ggplot(data=PREDS,
       aes(y=Exponent, 
           x=bbt_affected)) +
  geom_line(aes(lty=channel_region))+
  facet_wrap(~age, ncol=2) +
  scale_y_continuous(name = "Exponent") +
  scale_x_continuous(name = "BBT Affected")+
  theme_bw()+labs(lty="Region")+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text=element_text(size=12, color="black"), 
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.position = "bottom")


ggsave(
  filename="./Figure9.jpeg",
  plot = last_plot(),
  path = "C:/Users/lohse/Box/Lohse Brain Age EEG/reports/",
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)



