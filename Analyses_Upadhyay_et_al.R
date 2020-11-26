library(readxl)
library(dplyr)
library(emmeans)
library(ggplot2)
library(nlme)
library(multcomp)
library(caret)
library(EnvStats)
library(cowplot)
library(growth)

############## 
# UCLA data
##############

# Read in pass B2
b2bb8 <- read_xlsx("Pass B2 - All Machines, Main Stats.xlsx", sheet = "BB8")
b2r2nu <- read_xlsx("Pass B2 - All Machines, Main Stats.xlsx", sheet = "R2Nu")
b2r2sar <- read_xlsx("Pass B2 - All Machines, Main Stats.xlsx", sheet = "R2Sar")
b2253 <- read_xlsx("Pass B2 - All Machines, Main Stats.xlsx", sheet = "253")

b2bb8$Pass <- "B2"
b2bb8$Instrument <- "Config3B"

b2r2nu$Pass <- "B2"
b2r2nu$Instrument <- "Config3A"

b2r2sar$Pass <- "B2"
b2r2sar$Instrument <- "Config2"

b2253$Pass <- "B2"
b2253$Instrument <- "Config1"

# Read in pass B3
b3bb8 <- read_xlsx("Pass B3 - All Machines, Main Stats.xlsx", sheet = "BB8")
b3r2nu <- read_xlsx("Pass B3 - All Machines, Main Stats.xlsx", sheet = "R2Nu")
b3r2sar <- read_xlsx("Pass B3 - All Machines, Main Stats.xlsx", sheet = "R2Sar")
b3253 <- read_xlsx("Pass B3 - All Machines, Main Stats.xlsx", sheet = "253")

b3bb8$Pass <- "B3"
b3bb8$Instrument <- "Config3B"

b3r2nu$Pass <- "B3"
b3r2nu$Instrument <- "Config3A"

b3r2sar$Pass <- "B3"
b3r2sar$Instrument <- "Config2"

b3253$Pass <- "B3"
b3253$Instrument <- "Config1"

# Concatenate data
alldata <- Reduce(function(x, y) merge(x, y, all=TRUE), list(b2bb8, b2r2nu, b2r2sar, b2253, b3bb8, b3r2nu, b3r2sar, b3253))
colnames(alldata)[1] <- "SampleName"
colnames(alldata)[3] <- "d13CVPDBFinal"
colnames(alldata)[5] <- "d13CVPDBse"
colnames(alldata)[6] <- "d18OVPDBFinal"
colnames(alldata)[8] <- "d18OVPDBse"
colnames(alldata)[9] <- "D47CDESFinal"
colnames(alldata)[11] <- "D47CDESse"

alldata <- alldata[,-2]

# Produce a data summary
alldatsummary <- alldata %>%
  group_by(Pass, Instrument, `SampleName`) %>%
  summarize(n = n(), 
            `d13C VPDB (Final)` = mean(d13CVPDBFinal), 
            `d18O VPDB (Final)` = mean(`d18OVPDBFinal`),
            `D47 CDES (Final)` = mean(`D47CDESFinal`)
            )

####################################################

# Quick visualization
ggplot(data=alldata[alldata$d13CVPDBFinal > -40,], aes(x=SampleName, y=d13CVPDBFinal, fill=Instrument)) + geom_boxplot() + facet_wrap(~Pass)

alldata$SampleName <- as.character(alldata$SampleName)
options(contrasts = c("contr.sum","contr.poly"))
summary(model1 <- aov(d13CVPDBFinal ~ Pass*Instrument + Pass*SampleName + Instrument*SampleName + Error(SampleName), data=alldata))

# Setup for nlme and multcomp
alldata$SampleName <- as.factor(alldata$SampleName)
alldata$Instrument <- as.factor(alldata$Instrument)
alldata$Pass <- as.factor(alldata$Pass)

# Optional d13c and d18o analyses. Not included in paper but included here for informational purposes.
lme_d13c = lme(d13CVPDBFinal ~ SampleName + Pass*Instrument,
               data=alldata, random = ~1|d13CVPDBse, method = "REML")
summary(lme_d13c)

summary(glht(lme_d13c, linfct=mcp(SampleName = "Tukey", interaction_average = TRUE, covariate_average = TRUE)), 
        test = adjusted(type = "bonferroni"))
summary(glht(lme_d13c, linfct=mcp(Instrument = "Tukey", interaction_average = TRUE, covariate_average = TRUE)), 
        test = adjusted(type = "bonferroni"))
summary(glht(lme_d13c, linfct=mcp(Pass = "Tukey", interaction_average = TRUE, covariate_average = TRUE)), 
        test = adjusted(type = "bonferroni"))

# d180
lme_d18o = lme(d18OVPDBFinal ~ SampleName + Pass*Instrument, 
               data=alldata, random = ~1|d18OVPDBse, method = "REML")

summary(glht(lme_d18o, linfct=mcp(SampleName = "Tukey", interaction_average = TRUE, covariate_average = TRUE)), 
        test = adjusted(type = "bonferroni"))
summary(glht(lme_d18o, linfct=mcp(Instrument = "Tukey", interaction_average = TRUE, covariate_average = TRUE)), 
        test = adjusted(type = "bonferroni"))
summary(glht(lme_d18o, linfct=mcp(Pass = "Tukey", interaction_average = TRUE, covariate_average = TRUE)), 
        test = adjusted(type = "bonferroni"))

# These analyses included in manuscript
# Cap47
lme_D47 <- lme(D47CDESFinal ~ SampleName + Pass*Instrument, 
               data=alldata, random = ~1|D47CDESse, method = "REML")
summary(lme_D47)

modpairwise <- emmeans(lme_D47, pairwise ~ SampleName + Pass + Instrument)
modpairwise2 <- emmeans(lme_D47, pairwise ~ Pass + Instrument)
emmeans(lme_D47, pairwise ~ Pass)
emmeans(lme_D47, pairwise ~ Pass)

write.csv(modpairwise$contrasts, "contrasts.csv")
write.csv(modpairwise$emmeans, "emmeans.csv")

write.csv(modpairwise2$contrasts, "contrasts2.csv")
write.csv(modpairwise2$emmeans, "emmeans2.csv")


# Setup for precision analysis

precisiondf <- alldata[,c(8, 10:12)]
names(precisiondf)[1] <- "Standard" 

# Column of accepted standard values for D47 taken from Easotope
precisiondf <- precisiondf %>%
  mutate(D47CDESAccepted = case_when(
    grepl("ETH-1", Standard) ~ 0.263,
    grepl("ETH-2", Standard) ~ 0.26, 
    grepl("ETH-3", Standard) ~ 0.69, 
    grepl("ETH-4", Standard) ~ 0.507,
    grepl("102-GC-AZ01", Standard) ~ 0.713,
    grepl("Carmel Chalk", Standard) ~ 0.664,
    grepl("Carrara Marble", Standard) ~ 0.37,
    grepl("CMTile", Standard) ~ 0.376,
    grepl("IAEA-C1", Standard) ~ 0.362,
    grepl("IAEA-C2", Standard) ~ 0.723,
    grepl("MallinckrodtCal", Standard) ~ 0.533,
    grepl("MERCK", Standard) ~ 0.596,
    grepl("NBS 19", Standard) ~ 0.383,
    grepl("Spel 2-8-E", Standard) ~ 0.659, # Easotope notes this is a preliminary value still under testing
    grepl("SRM 88B", Standard) ~ 0.582,
    grepl("TV01", Standard) ~ 0.689, # Easotope notes this is a preliminary value still under testing
    grepl("TV03", Standard) ~ 0.7,
    grepl("Veinstrom", Standard) ~ 0.713
  )
  )

precisiondf %>%
  group_by(Standard) %>%
  summarize(n = n())

precisiondf <- precisiondf[precisiondf$Standard != "TV01",]  

# Set up models for precision tests
set.seed(280)
D47_index <- createDataPartition(precisiondf$D47CDESFinal, p = .75, list = FALSE)
D47_tr <- precisiondf[ D47_index, ]
D47_te <- precisiondf[ D47_index, ]

set.seed(7279)
lm_fit1 <- train(D47CDESFinal ~ Standard ,
                data = D47_tr, 
                method = "lm")
D47_pred1 <- predict(lm_fit1, D47_te)

lm_fit1

lm_fit2 <- train(D47CDESFinal ~ Standard + Instrument,
                 data = D47_tr, 
                 method = "lm")
D47_pred2 <- predict(lm_fit2, D47_te)

lm_fit2

lm_fit3 <- train(D47CDESFinal ~ Standard + Pass,
                 data = D47_tr, 
                 method = "lm")
D47_pred3 <- predict(lm_fit3, D47_te)

lm_fit3

lm_fit4 <- train(D47CDESFinal ~ Standard + Instrument + Pass,
                 data = D47_tr, 
                 method = "lm")
D47_pred4 <- predict(lm_fit4, D47_te)

lm_fit4

postResample(pred = D47_pred1, obs = D47_te$D47CDESFinal)
postResample(pred = D47_pred2, obs = D47_te$D47CDESFinal)
postResample(pred = D47_pred3, obs = D47_te$D47CDESFinal)
postResample(pred = D47_pred4, obs = D47_te$D47CDESFinal)

# Accuracy

# Calculate the percent error
accuracy <- precisiondf %>%
  group_by(Standard, Instrument, Pass) %>%
  mutate(error = ((D47CDESFinal - D47CDESAccepted)/D47CDESAccepted)*100
              )

# Get error summary stats
range(accuracy$error)
histogram(accuracy$error)
summary(accuracy$error)

# Make a table of the errors
accuracy %>%
  group_by(Standard, Instrument, Pass)  %>%
  summarise(meanerr = mean(error))

# Test for outliers
test <- rosnerTest(accuracy$error, k = 10)
test

# Drop nine true outilers
accuracy_outliersremoved <- accuracy[-c(5469, 1458, 1467, 1456, 3740, 1514, 4800, 3801, 4798),]

# Get summary stats again
range(accuracy_outliersremoved$error)
summary(accuracy_outliersremoved$error)

errtable <- accuracy_outliersremoved %>%
  group_by(Standard, Instrument, Pass)  %>%
  summarise(meanerr = mean(error),
            mederr = median(error),
            maxerr = max(error),
            minerr = min(error))

write.csv(errtable, "errortable.csv", row.names = FALSE)

ac1 <- ggplot(accuracy, aes(x=error)) + geom_histogram(bins = 50, color = "gray60", fill = "gray80") + theme_classic() + ylab("Count") + xlab("Error (%)") +
  xlim(-40,40) + ylim(0,1045) + geom_vline(xintercept=c(-3.24324,3.36606), linetype="dotted") + geom_vline(xintercept = 0.03292, linetype = "solid")

ac2 <- ggplot(accuracy_outliersremoved, aes(x=error)) + geom_histogram(bins = 50, color = "gray60", fill = "gray80") + theme_classic() + ylab("Count") + xlab("Error (%)") +
  xlim(-40,40)+ ylim(0,1045) + geom_vline(xintercept=c(-3.24324,3.36606), linetype="dotted") + geom_vline(xintercept = 0.05779, linetype = "solid")

cowplot::plot_grid(ac1, ac2, labels=c('A', 'B'))
ggsave("SIFig1.png", scale = 1.25, dpi = 600, width = 5.5, units = "in")

# Optional figures

# Pass B3
ggplot(alldatats[alldatats$Pass == "B3",], aes(x=DateTime, y=D47CDESFinal, color=Standard)) + geom_point() + geom_path() + facet_wrap(~Instrument) +
  scale_color_viridis_d() + theme_classic() + ylab(expression(paste(~Delta[47], " (CDES) \U2030"))) + xlab("Date")
  
ggsave("standardtsb3.tiff", dpi=600, compression = "lzw", width = 8, scale = 1.3, units = "in")

# Pass B2
ggplot(alldatats[alldatats$Pass == "B2",], aes(x=DateTime, y=D47CDESFinal, color=Standard)) + geom_point() + geom_path() + facet_wrap(~Instrument) +
  scale_color_viridis_d() + theme_classic() + ylab(expression(paste(~Delta[47], " (CDES) \U2030"))) + xlab("Date")

ggsave("standardtsb2.tiff", dpi=600, compression = "lzw", width = 8, scale = 1.3, units = "in")

