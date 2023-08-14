# cort analysis 
library(viridis)
library(glue)
library(lme4)
library(MASS)
library(car)
library(tidyverse)

#source
df <- read_csv("cort/FullCort.csv") 
head(df)

df<- df %>% dplyr::select(-Well, -final_conc_ng_ul) %>%unique(.)


#changing times
df$time <- ifelse(df$time == "1 hr", "70 min", df$time)
df$time <- ifelse(df$time == "25hr", "25 hr", df$time)
table(df$time)

##Changing period
head(df)
df$period <- ifelse(df$period == "Pre", 0, df$period)
df$period <- ifelse(df$period == "Post" & df$time == "70 min", 1, df$period)
df$period <- ifelse(df$period == "Post" & df$time == "25 hr", 2, df$period)

df <- df %>% 
  dplyr::select( -time)

##adding domgroup 
df$domgroup <- ifelse(df$Prerank == 3|df$Prerank==4, "Subordinate", "Dominant")  
## post domgroup 
df$domgroupPost <- ifelse(df$Postrank == 3|df$Postrank==4, "Subordinate", "Dominant")  

## Getting rid of duplicates 
df <- unique(df)


## making sure everything is a factor ( idk if you need to do this, but I think so)
df$period <- as.factor(df$period)
df$Prerank <- factor(df$Prerank, levels = c("1","3", "4"))
df$Postrank <- factor(df$Postrank, levels = c("1", "3","4"))
df$batch <- as.factor(df$batch)
df$plate <- as.factor(df$plate)
df$pre_idbatchcage <- as.factor(df$pre_idbatchcage)
df$condition<- factor(df$condition, levels = c("control", "same", "ascenders", "descenders"))
df$post_idbatch <- as.factor(df$post_idbatch)
df$group <- as.factor(df$group)


head(df)
colnames(df)
#period 0 = pre, 1 = 70min, 2 = 25hr

dfx <- df %>% select(plate, mean_con_ng_ul, batch, period, precage, post_idbatch, group, Prerank, Postrank, condition, domgroup, domgroupPost)
head(dfx)


source('functions/geom_boxjitter.R')

##### no differences between controls and reorganized group pre reorganization #####

g1 <- df %>%filter(period == "0")
g1$group <- ifelse(g1$group == "reorganized", "Reorganized","Control")
g1$group <- factor(g1$group, levels = c("Control", "Reorganized"))
head(g1)

g1$mean_con_ng_ul

#stats
hist(g1$mean_con_ng_ul) #not normal 

pre.glm <-glmer(mean_con_ng_ul~group +(1|batch)+(1|plate), data =g1,family = Gamma(link = "log"))

summary(pre.glm)
# Fixed effects:
#   Estimate Std. Error t value Pr(>|z|)    
# (Intercept)        4.4323     0.1225  36.183   <2e-16 ***
#   groupReorganized   0.1523     0.1166   1.307    0.191 
AIC(pre.glm)


#checks 
acf(resid(pre.glm))
qqPlot(resid(pre.glm))
hist(resid(pre.glm))
plot(pre.glm)

durbinWatsonTest(resid(pre.glm))
shapiro.test(resid(pre.glm)) # not normal :( 

#figure 
p1 <- ggplot(g1, aes(group, mean_con_ng_ul,color =group, fill = group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "",
       x = "",
       y = paste0("Corticosterone (ng/ml)", "\n", "Before Reorganzation"))+
  scale_color_manual(values = c("#22A884FF", "#414487FF")) +
  scale_fill_manual(values = c("#22A884FF", "#414487FF")) +theme_bw() +
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),    
        axis.title.x = element_text(hjust = 1),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        text = element_text(size = 20)
        
  )

p1

ggsave("manuscript/cort/results_figures/CORTB4_suppl.png", p1,height=4.5, width = 4.5, dpi = 600)

#### no difference between ranks pre reorganization ####
head(g1)
#stats
hist(g1$mean_con_ng_ul) #not normal 

g1$prerank2 <- factor(g1$Prerank, levels = c(3,4,1))

# compared to rank 1
prerank.glm <-glmer(mean_con_ng_ul~Prerank +(1|batch)+(1|plate), data =g1,family = Gamma(link = "log"))
summary(prerank.glm)
AIC(prerank.glm)

#checks 
acf(resid(prerank.glm))
qqPlot(resid(prerank.glm))
hist(resid(prerank.glm))
plot(prerank.glm)

durbinWatsonTest(resid(prerank.glm))
shapiro.test(resid(prerank.glm))

# compared to rank 3
prerank.glm2 <-glmer(mean_con_ng_ul~prerank2 +(1|batch)+(1|plate), data =g1,family = Gamma(link = "log"))
summary(prerank.glm2)
AIC(prerank.glm2)

#checks 
acf(resid(prerank.glm2))
qqPlot(resid(prerank.glm2))
hist(resid(prerank.glm2))
plot(prerank.glm2)

durbinWatsonTest(resid(prerank.glm2))
shapiro.test(resid(prerank.glm2))


#figure 
p2<- ggplot(g1, aes(domgroup, mean_con_ng_ul,color =domgroup, fill = domgroup))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "",
       x = "",
       y = paste0("Corticosterone (ng/ml)", "\n", "Before Reorganzation"))+
  scale_color_manual(values = viridis::viridis(2)) +
  scale_fill_manual(values = viridis::viridis(2)) +theme_bw() +
  theme(axis.text.x = element_text(vjust = 1,size = 20),
        # axis.ticks = element_blank(),
        axis.text.y = element_text(hjust = 0.5,size = 20),
        axis.text = element_text(color="#3C3C3C",size = 20),
        axis.title = element_text(size = 20),    
        axis.title.x = element_text(hjust = 1),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15),
        text = element_text(size = 20)
        
  )

p2

ggsave("manuscript/cort/results_figures/CORTB4_domgroup_suppl.png", p2,height=4.5, width = 4.5, dpi = 600)


pp <- gridExtra::grid.arrange(p1,p2, nrow=1)

ggsave("manuscript/cort/results_figures/groupDOMtogether_suppl.png", pp,height=4, width = 10, dpi = 600)

#### differences in post - pre cort for controls and reorganized groups ####
df1 <- read_csv("cort/FullCort.csv") 
head(df1)
df1<- df1 %>% dplyr::select(-Well, -final_conc_ng_ul) %>%unique(.)
df1$time <- ifelse(df1$time == "1 hr", "Reorganized 70 min", df1$time)
df1$time <- ifelse(df1$time == "25hr", "Reorganized 25 hr", df1$time)
## Time for the control data 
df1$time <- ifelse(df1$time == "Reorganized 70 min" & df1$group == "control","Control 70 min", df1$time)
df1$time <- ifelse(df1$time == "Reorganized 25 hr" & df1$group == "control","Control 25 hr", df1$time)

#Making things factors
# df1$period <- factor(df1$period, levels = c("Pre", "Post"))
# df1$time <- factor(df1$time,levels =c ("Reorganized 70 min", 'Control 70 min', "Reorganized 25 hr", "Control 25 hr"))

##adding domgroup 
df1$domgroup <- ifelse(df1$Prerank == 3|df1$Prerank==4, "Subordinate", "Dominant")  

##adding final condition groups
df1$condition1<- ifelse(df1$condition == "ascenders",  "Ascenders", df1$condition)
df1$condition1 <-ifelse(df1$condition == "descenders",  "Descenders", df1$condition1)
df1$condition1 <- ifelse(df1$condition == "same" & df1$domgroup== "Dominant", "Dominant", df1$condition1)
df1$condition1 <- ifelse(df1$condition == "same" & df1$domgroup== "Subordinate", "Subordinate", df1$condition1)
df1$condition1 <- ifelse(df1$condition == "control" & df1$time == "Control 70 min", "Control 70 min", df1$condition1) 
df1$condition1 <- ifelse(df1$condition == "control" & df1$time == "Control 25 hr", "Control 25 hr", df1$condition1)
df1$condition1 <- factor(df1$condition1, levels = c("Descenders","Dominant","Subordinate","Ascenders", "Control 70 min", "Control 25 hr"))
  
head(df1)
colnames(df1)

dp <- df1 %>% select(mean_con_ng_ul, plate, batch, period, time,pre_idbatchcage,post_idbatch, Prerank, Postrank, condition1)
head(dp)

dp <-dp %>% 
  group_by(pre_idbatchcage) %>% 
  pivot_wider(values_from = mean_con_ng_ul, names_from = period) %>% 
  mutate(diff = Post - Pre) %>% ungroup()
head(dp)

dp$time2 <- ifelse(grepl("70 min", dp$time), "70 min", "25 hr")
dp$time2 <- factor(dp$time2, levels = c("70 min", "25 hr"))

dp$condition1 <- gsub("25 hr", "", dp$condition1)
dp$condition1 <- gsub("70 min", "", dp$condition1)

dp$condition1 <- factor(dp$condition1, levels =c("Control ", "Ascenders", "Subordinate", "Descenders", "Dominant"))

ggplot(dp, aes(condition1, diff,color =condition1, fill = condition1))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.4,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  scale_color_manual(values = viridis::viridis(6)) +
  scale_fill_manual(values = viridis::viridis(6))+
  labs(y = "Difference in Corticosterone (ng/ml)", 
       x = "")+
  facet_wrap(~time2) +
  ylim(-400, 600)+
  theme_bw()+ theme(axis.text.x = element_text(vjust = 1, size = 15),
                    legend.position = 'none',
                    axis.text.y = element_text(hjust = 0.5, size = 15),
                    panel.grid.major = element_blank(),
                    # panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    plot.background = element_blank(),
                    axis.title = element_text(size = 20),
                    axis.text= element_text(color="#3C3C3C", size=20),
                    strip.text.x = element_text(size = 15),
                    strip.background = element_rect(color = "black", fill = NA, size = 1),
                    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

  
#stats 70 min 

head(dp)
dp70 <- dp %>% filter(time2 == "70 min")
hist(dp70$diff)
dp70$condition1 <- gsub(" ", "", dp70$condition1)

dp70$condition1 <- factor(dp70$condition1, levels = c('Control',"Ascenders", "Subordinate", "Descenders", "Dominant"))
dp70$condition2 <- factor(dp70$condition1, levels = c("Ascenders", "Subordinate", "Descenders", "Dominant", "Control"))
dp70$condition3<- factor(dp70$condition1, levels = c("Subordinate", "Descenders", "Dominant", "Control", "Ascenders"))
dp70$condition4<- factor(dp70$condition1, levels = c("Descenders", "Dominant", "Control", "Ascenders", "Subordinate"))


# Z-Score Standardization: (X – μ) / σ
dm <- mean(dp70$diff)
ds <- sd(dp70$diff)
dp70$zdiff <- (dp70$diff -dm)/ds
hist(dp70$zdiff)

#Min-Max Normalization: (X – min(X)) / (max(X) – min(X))
dp70$ndiff <- (dp70$diff - min(dp70$diff))/ (max(dp70$diff) - min(dp70$diff)) 
hist(dp70$ndiff)

head(dp70)
#stats that compare groups to controls
c70.lm <-glm(ndiff~condition1 +(1|batch)+(1|plate) , data =dp70)
summary(c70.lm)
AIC(c70.lm)  

# Coefficients: (2 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              69.65      26.40   2.639 0.010514 *  
#   condition1Ascenders      22.48      39.91   0.563 0.575274    
# condition1Subordinate    69.14      39.91   1.733 0.088143 .  
# condition1Descenders    136.73      39.91   3.426 0.001092 ** 
#   condition1Dominant      175.76      49.89   3.523 0.000807 ***
#   1 | batchTRUE               NA         NA      NA       NA    
# 1 | plateTRUE               NA         NA      NA       NA    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 112 on 62 degrees of freedom
# Multiple R-squared:  0.2502,	Adjusted R-squared:  0.2018 
# F-statistic: 5.171 on 4 and 62 DF,  p-value: 0.001165


#stats that compare groups to ascenders 
c70.lm2 <-lm(diff~condition2 +(1|batch)+(1|plate), data =dp70)
summary(c70.lm2)
AIC(c70.lm2)  

# Coefficients: (2 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)              92.13      29.93   3.078  0.00310 **
#   condition2Subordinate    46.66      42.33   1.102  0.27454   
# condition2Descenders    114.25      42.33   2.699  0.00895 **
#   condition2Dominant      153.28      51.84   2.957  0.00439 **
#   condition2Control       -22.48      39.91  -0.563  0.57527   
# 1 | batchTRUE               NA         NA      NA       NA   
# 1 | plateTRUE               NA         NA      NA       NA   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 112 on 62 degrees of freedom
# Multiple R-squared:  0.2502,	Adjusted R-squared:  0.2018 
# F-statistic: 5.171 on 4 and 62 DF,  p-value: 0.001165


#stats that compare groups to subs
c70.lm3 <-lm(diff~condition3 +(1|batch)+(1|plate), data =dp70)
summary(c70.lm3)
AIC(c70.lm3)  

# Coefficients: (2 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            138.80      29.93   4.637 1.87e-05 ***
#   condition3Descenders    67.59      42.33   1.597   0.1154    
# condition3Dominant     106.62      51.84   2.057   0.0439 *  
#   condition3Control      -69.14      39.91  -1.733   0.0881 .  
# condition3Ascenders    -46.66      42.33  -1.102   0.2745    
# 1 | batchTRUE              NA         NA      NA       NA    
# 1 | plateTRUE              NA         NA      NA       NA    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 112 on 62 degrees of freedom
# Multiple R-squared:  0.2502,	Adjusted R-squared:  0.2018 
# F-statistic: 5.171 on 4 and 62 DF,  p-value: 0.001165


#stats that compare groups to descenders
c70.lm4 <-lm(diff~condition4 +(1|batch)+(1|plate), data =dp70)
summary(c70.lm4)
AIC(c70.lm4)  
# 
# Coefficients: (2 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             206.38      29.93   6.895 3.26e-09 ***
#   condition4Dominant       39.03      51.84   0.753  0.45438    
# condition4Control      -136.73      39.91  -3.426  0.00109 ** 
#   condition4Ascenders    -114.25      42.33  -2.699  0.00895 ** 
#   condition4Subordinate   -67.59      42.33  -1.597  0.11542    
# 1 | batchTRUE               NA         NA      NA       NA    
# 1 | plateTRUE               NA         NA      NA       NA    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



### stats 25 hr

head(dp)
dp25 <- dp %>% filter(time2 == "25 hr")
hist(dp70$diff)
dp25$condition1 <- gsub(" ", "", dp25$condition1)

dp25$condition1 <- factor(dp25$condition1, levels = c('Control',"Ascenders", "Subordinate", "Descenders", "Dominant"))
dp25$condition2 <- factor(dp25$condition1, levels = c("Ascenders", "Subordinate", "Descenders", "Dominant", "Control"))
dp25$condition3<- factor(dp25$condition1, levels = c("Subordinate", "Descenders", "Dominant", "Control", "Ascenders"))
dp25$condition4<- factor(dp25$condition1, levels = c("Descenders", "Dominant", "Control", "Ascenders", "Subordinate"))

# compared to controls
c25.lm <-lm(diff~condition1 +(1|batch)+(1|plate) , data =dp25)
summary(c25.lm)
AIC(c25.lm)  

# no differences 


# compared to ascenders
c25.lm2 <-lm(diff~condition2 +(1|batch)+(1|plate) , data =dp25)
summary(c25.lm2)
AIC(c25.lm2)  

# no differences 


# compared to ascenders
c25.lm3 <-lm(diff~condition3 +(1|batch)+(1|plate) , data =dp25)
summary(c25.lm3)
AIC(c25.lm3)  

# no differences 

# compared to descenders
c25.lm4 <-lm(diff~condition4 +(1|batch)+(1|plate) , data =dp25)
summary(c25.lm4)
AIC(c25.lm3)  

# no differences 
