install.packages("dplyr")
install.packages("readr")
install.packages("readxl") 
install.packages("epitools")
install.packages("tidyr")
install.packages("vcd")
install.packages("MASS")
install.packages("ggplot2")

library(dplyr)
library(readr)
library(epitools)
library(tidyr)
library(vcd)
library(MASS)
library(ggplot2)

# NSDUH_2019 Dataset Import

setwd("C:/File/Path.txt)

NSDUH_2019 <- read_delim("NSDUH_2019_Tab.txt", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)



attach(NSDUH_2019)

## Transforming values
### Changing 2 to 1 for binomial regression of AMDEYR, AUUNMTXYR, AUOPTYR,and CATAG2


NSDUH_sub <- data.frame(CATAG2, AMDEYR, 
                             AUUNMTYR, AUOPTYR, ACTD2001, ADTMTHLP)


NSDUH_sub["AMDEYR"][NSDUH_sub["AMDEYR"] == 2] <- 0
NSDUH_sub["AUUNMTYR"][NSDUH_sub["AUUNMTYR"] == 2] <- 0
NSDUH_sub["AUOPTYR"][NSDUH_sub["AUOPTYR"] == 2] <- 0
NSDUH_sub["ACTD2001"][NSDUH_sub["ACTD2001"] == 2] <- 0
NSDUH_sub <- NSDUH_sub[NSDUH_sub$CATAG2 != 1,]
NSDUH_sub["CATAG2"][NSDUH_sub["CATAG2"] == 2] <- 0
NSDUH_sub["CATAG2"][NSDUH_sub["CATAG2"] == 3] <- 1



### Filtering variables ACTD2001, AUUNMTYR, and AUOPTYR



NSDUH_sub <- NSDUH_sub[NSDUH_sub$ACTD2001 != 85,]
NSDUH_sub <- NSDUH_sub[NSDUH_sub$ACTD2001 != 97,]
NSDUH_sub <- NSDUH_sub[NSDUH_sub$ACTD2001 != 98,]
NSDUH_sub <- NSDUH_sub[NSDUH_sub$ACTD2001 != 99,]

NSDUH_sub <- NSDUH_sub[NSDUH_sub$AUUNMTYR != 85,]
NSDUH_sub <- NSDUH_sub[NSDUH_sub$AUUNMTYR != 94,]
NSDUH_sub <- NSDUH_sub[NSDUH_sub$AUUNMTYR != 97,]

NSDUH_sub <- NSDUH_sub[NSDUH_sub$AUOPTYR != 3,]
NSDUH_sub <- NSDUH_sub[NSDUH_sub$AUOPTYR != 85,]
NSDUH_sub <- NSDUH_sub[NSDUH_sub$AUOPTYR != 94,]
NSDUH_sub <- NSDUH_sub[NSDUH_sub$AUOPTYR != 97,]


# Testing hypotheses

## Generating frequency table for hypothesis 1 to confirm cell counts are >= 5

attach(NSDUH_sub)

table(AMDEYR, ACTD2001)

## Binomial logistic regression for hypothesis 1

GLM1 <- glm(AMDEYR ~ factor(ACTD2001), data = NSDUH_sub, family = "binomial")
summary(GLM1)

exp(coefficients(GLM1))

GLM1_CI <- confint((GLM1))
exp(GLM1_CI)


## Generating frequency table for hypothesis 2 to confirm cell counts are >= 5

table(AUOPTYR, ACTD2001)

## Binomial logistic regression for hypothesis 2

GLM2 <- glm(AUOPTYR ~ ACTD2001, data = NSDUH_sub, family = "binomial")
summary(GLM2)

exp(coefficients(GLM2))

GLM2_CI <- confint(GLM2)

exp(GLM2_CI)


## Generating frequency table for hypothesis 3 to confirm cell counts are >= 5

table(AUUNMTYR, ACTD2001)

## Binomial logistic regression for hypothesis 3

GLM3 <- glm(AUUNMTYR ~ ACTD2001, data = NSDUH_sub, family = "binomial")
summary(GLM3)

exp(coefficients(GLM3))

GLM3_CI <- confint(GLM3)

exp(GLM3_CI)

# Creating subset for ACTD2001 = Yes

NSDUH_2019_ACTD2001 <- filter(NSDUH_sub, ACTD2001 == 1)


## Generating frequency table for hypothesis 4 to confirm cell counts are >= 5

attach(NSDUH_2019_ACTD2001)

table(AUOPTYR, CATAG2)

## Binomial logistic for hypothesis 4

GLM4 <- glm(AUOPTYR ~ CATAG2, data = NSDUH_2019_ACTD2001, family = "binomial")
summary(GLM4)

exp(coefficients(GLM4))

GLM4_CI <- confint(GLM4)

exp(GLM4_CI)


## Creating graphic for logistic regression results

### Create labels for plot
boxLabels = c("Major Depressive Episode", "Outpatient Treatment (>= 26 years old)", "Outpatient Treatment (>= 18 years old)", "Perceived Unmet Need")

# Enter OR and CI data. 


ORPlot <- data.frame(yAxis = length(boxLabels):1, 
                 boxOdds = c(2.04, 2.75, 2.01, 3.03), 
                 boxCILow = c(1.43, 1.32, 1.47, 2.01), 
                 boxCIHigh = c(2.93, 6.70, 2.76, 4.64))


ORPlot <- ggplot(ORPlot, aes(x = boxOdds, y = boxLabels)) + 
    geom_vline(aes(xintercept = 1.0), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = 
                     .2, color = "gray50") +
    geom_point(size = 3.5, color = "orange") +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("Odds ratio") +
    annotate(geom = "text", y =1.3, x = 2.0, label = "p = 8.52e-05", size = 3, hjust = 0) +
    annotate(geom = "text", y =2.3, x = 2, label = "p = 1.5e-05", size = 3, hjust = 0) +
    annotate(geom = "text", y =3.3, x = 2.75, label = "p = 0.013", size = 3, hjust = 0) +
    annotate(geom = "text", y =4.3, x = 3.05, label = "p = 2.02e-07", size = 3, hjust = 0) +
      ggtitle("Active Duty In or After 2001")

ORPlot

## Determining How Much Outpatient Mental Health Treatment Helped Depression by Percent for ACTD2001 = YES (1) and NO (0)


attach(NSDUH_sub)

PCT_HLP_COM <- data.frame(ADTMTHLP, ACTD2001)
PCT_HLP_COM <- filter(PCT_HLP_COM, ADTMTHLP != "94")
PCT_HLP_COM <- filter(PCT_HLP_COM, ADTMTHLP != "97")
PCT_HLP_COM <- filter(PCT_HLP_COM, ADTMTHLP != "98")
PCT_HLP_COM <- filter(PCT_HLP_COM, ADTMTHLP != "99")
PCT_HLP_COM$ADTMTHLP <- as.factor(PCT_HLP_COM$ADTMTHLP)
PCT_HLP_COM$ACTD2001 <- as.factor(PCT_HLP_COM$ACTD2001)
PCT_HLP_COM <- PCT_HLP_COM %>%
  group_by(ACTD2001, ADTMTHLP) %>%
  summarise(percent = 100 * n () / nrow(PCT_HLP_COM))

PCT_HLP_COM


## Creating Grouped Bar Chart for ADTMTHLP by ACTD2001 Status


ggplot(PCT_HLP_COM, aes(fill = ACTD2001,
                        y=percent, 
                        x=ADTMTHLP)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label=round(percent, digits = 2)), position=position_dodge(width=0.9), vjust=-0.25) +
  scale_fill_discrete(labels = c("No", "Yes")) +
  labs(title = "Perceived Level of Treatment Effectiveness by Active-Duty Status", x = "", y = "Percent (%)")
