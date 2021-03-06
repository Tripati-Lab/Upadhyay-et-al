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

# Read in pass AnalysisA
AnalysisAbb8 <- read_xlsx("Pass AnalysisA - All Machines, Main Stats.xlsx", sheet = "BB8")
AnalysisAr2nu <- read_xlsx("Pass AnalysisA - All Machines, Main Stats.xlsx", sheet = "R2Nu")
AnalysisAr2sar <- read_xlsx("Pass AnalysisA - All Machines, Main Stats.xlsx", sheet = "R2Sar")
AnalysisA253 <- read_xlsx("Pass AnalysisA - All Machines, Main Stats.xlsx", sheet = "253")
additional <- read_xlsx("Additional standards Bulk and Clumped.xlsx")

AnalysisAbb8$Pass <- "AnalysisA"
AnalysisAbb8$Instrument <- "Config3B"

AnalysisAr2nu$Pass <- "AnalysisA"
AnalysisAr2nu$Instrument <- "Config3A"

AnalysisAr2sar$Pass <- "AnalysisA"
AnalysisAr2sar$Instrument <- "Config2"

AnalysisA253$Pass <- "AnalysisA"
AnalysisA253$Instrument <- "Config1"

additional$Pass <- "AnalysisA"

# Read in pass AnalysisB
AnalysisBbb8 <- read_xlsx("Pass AnalysisB - All Machines, Main Stats.xlsx", sheet = "BB8")
AnalysisBr2nu <- read_xlsx("Pass AnalysisB - All Machines, Main Stats.xlsx", sheet = "R2Nu")
AnalysisBr2sar <- read_xlsx("Pass AnalysisB - All Machines, Main Stats.xlsx", sheet = "R2Sar")
AnalysisAnalysisA53 <- read_xlsx("Pass AnalysisB - All Machines, Main Stats.xlsx", sheet = "253")

AnalysisBbb8$Pass <- "AnalysisB"
AnalysisBbb8$Instrument <- "Config3B"

AnalysisBr2nu$Pass <- "AnalysisB"
AnalysisBr2nu$Instrument <- "Config3A"

AnalysisBr2sar$Pass <- "AnalysisB"
AnalysisBr2sar$Instrument <- "Config2"

AnalysisAnalysisA53$Pass <- "AnalysisB"
AnalysisAnalysisA53$Instrument <- "Config1"

# Concatenate data
alldata <- Reduce(function(x, y) merge(x, y, all=TRUE), list(AnalysisAbb8, AnalysisAr2nu, AnalysisAr2sar, AnalysisA253, AnalysisBbb8, AnalysisBr2nu, AnalysisBr2sar, AnalysisAnalysisA53, additional))

# Standardize names
names(alldata)<-gsub("\\s","",names(alldata))
names(alldata)<-gsub("_","",names(alldata))
names(alldata) <- make.names(names(alldata))

alldata <- alldata[,-6] # Don't need the Standard column

# Standardize name
alldata$SampleName <- recode(alldata$SampleName, "Spel 2-8-E" = "SPEL-2-8-E")

# Calculate SE
alldata <- alldata %>%
  group_by(Pass, Instrument, SampleName) %>%
  mutate(D47CDES.Final.SE = sd(D47CDES.Final.)/sqrt(length(D47CDES.Final.)),
         d13CVPDB.Final.SE = sd(d13CVPDB.Final.)/sqrt(length(d13CVPDB.Final.)),
         d18OVPDB.Final.SE = sd(d18OVPDB.Final.)/sqrt(length(d18OVPDB.Final.))
  )

# Drop singular observation (Carrara Marble, Pass AnalysisA, Config 3A)
alldata <- alldata[!is.na(alldata$D47CDES.Final.SE),]

# Produce a data summary
alldatsummary <- alldata %>%
  group_by(Pass, Instrument, `SampleName`) %>%
  summarize(n = n(), 
            `d13C VPDB (Final)` = mean(`d13CVPDB.Final.`), 
            `d18O VPDB (Final)` = mean(d18OVPDB.Final.),
            `D47 CDES (Final)` = mean(D47CDES.Final.)
            )

####################################################

# Quick visualization
ggplot(data=alldata[alldata$`d13CVPDB.Final.` > -40,], aes(x=SampleName, y=`d13CVPDB.Final.`, fill=Instrument)) + geom_boxplot() + facet_wrap(~Pass)

alldata$SampleName <- as.character(alldata$SampleName)
options(contrasts = c("contr.sum","contr.poly"))
summary(model1 <- aov(`d13CVPDB.Final.` ~ Pass*Instrument + Pass*SampleName + Instrument*SampleName + Error(SampleName), data=alldata))


# Remove CIT standards ahead of analysis
alldata <- alldata[!grepl("CIT", alldata$SampleName), ]

# Setup for nlme and multcomp
alldata$SampleName <- as.factor(alldata$SampleName)
alldata$Instrument <- as.factor(alldata$Instrument)
alldata$Pass <- as.factor(alldata$Pass)

# Optional d13c and d18o analyses. Not included in paper but included here for informational purposes.
lme_d13c = lme(`d13CVPDB.Final.` ~ SampleName + Pass*Instrument,
               data=alldata, random = ~1|d13CVPDBse, method = "REML")
summary(lme_d13c)

summary(glht(lme_d13c, linfct=mcp(SampleName = "Tukey", interaction_average = TRUE, covariate_average = TRUE)), 
        test = adjusted(type = "bonferroni"))
summary(glht(lme_d13c, linfct=mcp(Instrument = "Tukey", interaction_average = TRUE, covariate_average = TRUE)), 
        test = adjusted(type = "bonferroni"))
summary(glht(lme_d13c, linfct=mcp(Pass = "Tukey", interaction_average = TRUE, covariate_average = TRUE)), 
        test = adjusted(type = "bonferroni"))

# d180
lme_d18o = lme(d18OVPDB.Final. ~ SampleName + Pass*Instrument, 
               data=alldata, random = ~1|d18OVPDBse, method = "REML")

summary(glht(lme_d18o, linfct=mcp(SampleName = "Tukey", interaction_average = TRUE, covariate_average = TRUE)), 
        test = adjusted(type = "bonferroni"))
summary(glht(lme_d18o, linfct=mcp(Instrument = "Tukey", interaction_average = TRUE, covariate_average = TRUE)), 
        test = adjusted(type = "bonferroni"))
summary(glht(lme_d18o, linfct=mcp(Pass = "Tukey", interaction_average = TRUE, covariate_average = TRUE)), 
        test = adjusted(type = "bonferroni"))

# These analyses included in manuscript
# Cap47
lme_D47 <- lme(D47CDES.Final. ~ SampleName + Pass*Instrument, 
               data=alldata, random = ~1|D47CDES.Final.SE, method = "REML")
summary(lme_D47)

modpairwise <- emmeans(lme_D47, pairwise ~ SampleName + Pass + Instrument)
modpairwise2 <- emmeans(lme_D47, pairwise ~ Pass + Instrument)
emmeans(lme_D47, pairwise ~ Pass)
emmeans(lme_D47, pairwise ~ Instrument)

write.csv(modpairwise$contrasts, "contrasts_June_revision.csv")
write.csv(modpairwise$emmeans, "emmeans_June_revision.csv")

write.csv(modpairwise2$contrasts, "contrasts2.csv")
write.csv(modpairwise2$emmeans, "emmeans2.csv")

yaxislabs <- c("Analysis 1 Config 1", "Analysis 2 Config 1", "Analysis 1 Config 2", "Analysis 2 Config 2",
               "Analysis 1 Config 3A", "Analysis 2 Config 3A", "Analysis 1 Config 3B", "Analysis 2 Config 3B")


plot(emmeans(lme_D47, pairwise ~ Pass + Instrument), comparison = TRUE, colors = "darkblue") + 
  scale_y_discrete(labels = yaxislabs) +
  xlab("Estimated marginal mean") + ylab("Analysis by instrument") +
  theme_bw()

ggsave("emmplot_Junerevision.png", dpi = 600, height = 7, width = 5, units = "in")
ggsave("SI Fig 5.tiff", dpi = 1200, height = 7, width = 5, units = "in", compression = "lzw")

# Setup for precision analysis

precisiondf <- alldata[,c(1, 3:5)]
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
    grepl("SPEL-2-8-E", Standard) ~ 0.659, # Easotope notes this is a preliminary value still under testing
    grepl("SRM 88B", Standard) ~ 0.582,
    grepl("TV01", Standard) ~ 0.689, # Easotope notes this is a preliminary value still under testing
    grepl("TV03", Standard) ~ 0.7,
    grepl("Veinstrom", Standard) ~ 0.713
  )
  )

precisiondf %>%
  group_by(Standard) %>%
  summarize(n = n())

# Set up models for precision tests
set.seed(280)
D47_index <- createDataPartition(precisiondf$D47CDES.Final., p = .75, list = FALSE)
D47_tr <- precisiondf[ D47_index, ]
D47_te <- precisiondf[ D47_index, ]

set.seed(7279)
lm_fit1 <- train(D47CDES.Final. ~ Standard ,
                data = D47_tr, 
                method = "lm")
D47_pred1 <- predict(lm_fit1, D47_te)

lm_fit1

lm_fit2 <- train(D47CDES.Final. ~ Standard + Instrument,
                 data = D47_tr, 
                 method = "lm")
D47_pred2 <- predict(lm_fit2, D47_te)

lm_fit2

lm_fit3 <- train(D47CDES.Final. ~ Standard + Pass,
                 data = D47_tr, 
                 method = "lm")
D47_pred3 <- predict(lm_fit3, D47_te)

lm_fit3

lm_fit4 <- train(D47CDES.Final. ~ Standard + Instrument + Pass,
                 data = D47_tr, 
                 method = "lm")
D47_pred4 <- predict(lm_fit4, D47_te)

lm_fit4

postResample(pred = D47_pred1, obs = D47_te$D47CDES.Final.)
postResample(pred = D47_pred2, obs = D47_te$D47CDES.Final.)
postResample(pred = D47_pred3, obs = D47_te$D47CDES.Final.)
postResample(pred = D47_pred4, obs = D47_te$D47CDES.Final.)

# Accuracy

# Calculate the percent error
accuracy <- precisiondf %>%
  group_by(Standard, Instrument, Pass) %>%
  mutate(error = ((D47CDES.Final. - D47CDESAccepted)/D47CDESAccepted)*100,
         absoluteerror = D47CDES.Final. - D47CDESAccepted
              )

# Get error summary stats
range(accuracy$error)
histogram(accuracy$error)
summary(accuracy$error)

range(accuracy$absoluteerror)
histogram(accuracy$absoluteerror)
summary(accuracy$absoluteerror)

# Make a table of the errors
accuracy %>%
  group_by(Standard, Instrument, Pass)  %>%
  summarise(meanerr = mean(error))

# Test for outliers using percent error
test <- rosnerTest(accuracy$error, k = 10)
test

# Test for outliers using absolute error
test2 <- rosnerTest(accuracy$absoluteerror, k = 10)
test2

# Drop true outilers
accuracy_outliersremoved_pcterr <- accuracy[-dput(as.numeric(test$all.stats$Obs.Num[test$all.stats$Outlier == TRUE])),]
accuracy_outliersremoved_absterr <- accuracy[-dput(as.numeric(test2$all.stats$Obs.Num[test2$all.stats$Outlier == TRUE])),]

# Get summary stats again
range(accuracy_outliersremoved_pcterr$error)
summary(accuracy_outliersremoved_pcterr$error)
summary(accuracy_outliersremoved_pcterr$absoluteerror)
summary(accuracy_outliersremoved_absterr$error)
summary(accuracy_outliersremoved_absterr$absoluteerror)

hist(accuracy_outliersremoved$absoluteerror)

# Absolute error
errtable_absterr <- accuracy_outliersremoved_absterr %>%
  group_by(Standard, Instrument, Pass)  %>%
  summarise(meanerr = mean(error),
            mederr = median(error),
            maxerr = max(error),
            minerr = min(error),
            absoluteerrormeanerr = mean(absoluteerror),
            absoluteerrormederr = median(absoluteerror),
            absoluteerrormaxerr = max(absoluteerror),
            absoluteerrorminerr = min(absoluteerror))

write.csv(errtable_absterr, "errortable_absterr.csv", row.names = FALSE)

# Percent error
errtable_pcterr <- accuracy_outliersremoved_pcterr %>%
  group_by(Standard, Instrument, Pass)  %>%
  summarise(meanerr = mean(error),
            mederr = median(error),
            maxerr = max(error),
            minerr = min(error),
            absoluteerrormeanerr = mean(absoluteerror),
            absoluteerrormederr = median(absoluteerror),
            absoluteerrormaxerr = max(absoluteerror),
            absoluteerrorminerr = min(absoluteerror))

write.csv(errtable_pcterr, "errortable_pcterr.csv", row.names = FALSE)

# Histograms of errors
ac1 <- ggplot(accuracy, aes(x=error)) + geom_histogram(bins = 50, color = "gray60", fill = "gray80") + theme_classic() + ylab("Count") + xlab("Error (%)") +
  xlim(-40,40) + ylim(0,1045) + geom_vline(xintercept=c(-3.24324,3.36606), linetype="dotted") + geom_vline(xintercept = 0.03292, linetype = "solid")

ac2 <- ggplot(accuracy, aes(x=absoluteerror)) + geom_histogram(bins = 50, color = "gray60", fill = "gray80") + theme_classic() + ylab("Count") + xlab("Error (‰)") +
  xlim(-0.17,0.17) + ylim(0,1045) + geom_vline(xintercept=c(-0.015,0.016), linetype="dotted") + geom_vline(xintercept = 0.0002312, linetype = "solid")

ac3 <- ggplot(accuracy_outliersremoved_pcterr, aes(x=error)) + geom_histogram(bins = 50, color = "gray60", fill = "gray80") + theme_classic() + ylab("Count") + xlab("Error (%)") +
  xlim(-40,40)+ ylim(0,1045) + geom_vline(xintercept=c(-3.24324,3.36606), linetype="dotted") + geom_vline(xintercept = 0.05779, linetype = "solid")

ac4 <- ggplot(accuracy_outliersremoved_pcterr, aes(x=absoluteerror)) + geom_histogram(bins = 50, color = "gray60", fill = "gray80") + theme_classic() + ylab("Count") + xlab("Error (‰)") +
  xlim(-0.17,0.17)+ ylim(0,1045) + geom_vline(xintercept=c(-0.015,0.016), linetype="dotted") + geom_vline(xintercept = 0.0003431, linetype = "solid")

ac5 <- ggplot(accuracy_outliersremoved_absterr, aes(x=error)) + geom_histogram(bins = 50, color = "gray60", fill = "gray80") + theme_classic() + ylab("Count") + xlab("Error (%)") +
  xlim(-40,40)+ ylim(0,1045) + geom_vline(xintercept=c(-3.24324,3.36606), linetype="dotted") + geom_vline(xintercept = 0.07502, linetype = "solid")

ac6 <- ggplot(accuracy_outliersremoved_absterr, aes(x=absoluteerror)) + geom_histogram(bins = 50, color = "gray60", fill = "gray80") + theme_classic() + ylab("Count") + xlab("Error (‰)") +
  xlim(-0.17,0.17)+ ylim(0,1045) + geom_vline(xintercept=c(-0.015,0.016), linetype="dotted") + geom_vline(xintercept = 0.0004328, linetype = "solid")

cowplot::plot_grid(ac1, ac2, ac3, ac4, ac5, ac6, 
                   ncol = 2, nrow = 3,
                   labels=c('A', 'B', 'C', 'D', 'E', 'F'))

ggsave("SIFig1.png", scale = 1.25, dpi = 600, width = 5.5, units = "in")
ggsave("SIFig1.tiff", scale = 1.25, dpi = 1200, width = 5.5, units = "in", compression = "lzw")

# Optional figures

# Pass AnalysisB
ggplot(alldatats[alldatats$Pass == "AnalysisB",], aes(x=DateTime, y=D47CDES.Final., color=Standard)) + geom_point() + geom_path() + facet_wrap(~Instrument) +
  scale_color_viridis_d() + theme_classic() + ylab(expression(paste(~Delta[47], " (CDES) \U2030"))) + xlab("Date")
  
ggsave("standardtsAnalysisB.tiff", dpi=600, compression = "lzw", width = 8, scale = 1.3, units = "in")

# Pass AnalysisA
ggplot(alldatats[alldatats$Pass == "AnalysisA",], aes(x=DateTime, y=D47CDES.Final., color=Standard)) + geom_point() + geom_path() + facet_wrap(~Instrument) +
  scale_color_viridis_d() + theme_classic() + ylab(expression(paste(~Delta[47], " (CDES) \U2030"))) + xlab("Date")

ggsave("standardtsAnalysisA.tiff", dpi=600, compression = "lzw", width = 8, scale = 1.3, units = "in")

