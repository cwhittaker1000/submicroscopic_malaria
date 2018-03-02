library(ggplot2)
library(LMest)

data_frame <- as.data.frame(Playing_Around_15th_November_With_Walldorf_Data)
data_frame <- as.data.frame(Rough_Preliminary_Data_Exploration_For_R)

data_frame$Old_Or_New <- as.factor(data_frame$Old_Or_New)
data_frame$Global_Region <- as.factor(data_frame$Global_Region)
data_frame <- data_frame[order(data_frame$PCR_Percentage),]

subset_low <- data_frame[data_frame$PCR_Prev < 0.2,]
subset_mid <- data_frame[data_frame$PCR_Prev < 0.5 & data_frame$PCR_Prev >= 0.2,]
subset_high <- data_frame[data_frame$PCR_Prev > 0.5,]

ggplot(data_frame[1:120,], aes(x= PCR_Percentage, y=LM_Percentage)) +
  geom_point(aes(colour = Old_Or_New), size =2)+
  scale_color_manual(values = c("#00AFB5", "#FF7700"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75), plot.title = element_text(hjust = 0),
        axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
        legend.position = "bottom", legend.title = element_blank()) +
  scale_y_continuous(breaks = c(0,0.05, 0.1, 0.15,0.2)) +
  scale_x_continuous(breaks = c(0,0.2, 0.4, 0.6, 0.8, 1.0)) +
  labs(x = "% PCR Positive", y = "% Microscopy Positive")


ggplot(data_frame[1:120,], aes(x= PCR_Percentage, y=LM_Percentage, colour = Global_Region)) +
  geom_point(size =2 )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75), plot.title = element_text(hjust = 0),
        axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"), legend.position = "bottom") +
  #scale_y_continuous(breaks = c(0,0.2, 0.4, 0.6, 0.8, 1.0,1.2,1.4)) +
  #scale_x_continuous(breaks = c(0,0.2, 0.4, 0.6, 0.8, 1.0,1.2,1.4)) +
  labs(x = "% PCR Positive", y = "Proportion Detected by Microscopy")

#Prep for Fitting Various Models to First 120 Data Points
data_frame_subset <- data_frame[1:120,]

South_America <- data_frame_subset[data_frame_subset$Global_Region == "South America",]
West_Africa <- data_frame_subset[data_frame_subset$Global_Region == "West Africa",]
East_Africa <- data_frame_subset[data_frame_subset$Global_Region == "East Africa",]
Asia <- data_frame_subset[data_frame_subset$Global_Region == "Asia",]

#Fitting Linear Models To Data For First 120 Data Points
summary(lm(West_Africa$PCR_Percentage ~ 0 + West_Africa$LM_Percentage))
summary(lm(East_Africa$PCR_Percentage ~ 0 + East_Africa$LM_Percentage))
summary(lm(Asia$PCR_Percentage ~ 0 + Asia$LM_Percentage))
summary(lm(South_America$PCR_Percentage ~ 0 + South_America$LM_Percentage))

LM <- seq(0,0.2,0.001)

plot(South_America$PCR_Percentage, South_America$LM_Percentage, xlim = c(0, 0.2), ylim = c(0,0.2), pch = 20, col = "blue")
lines(4.9554*LM, LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Percentage, West_Africa$LM_Percentage, xlim = c(0, 0.2), ylim = c(0,0.2), pch = 20, col = "red")
lines(1.3616*LM, LM, col = "red", lwd = 2)
points(East_Africa$PCR_Percentage, East_Africa$LM_Percentage, xlim = c(0,0.2), ylim = c(0,0.2), pch = 20, col = "green")
lines(1.7707*LM, LM, col = "green", lwd = 2)
points(Asia$PCR_Percentage, Asia$LM_Percentage, xlim = c(0,0.2), ylim = c(0,0.2), pch = 20, col = "black")
lines(1.6157*LM, LM, col = "black", lwd = 2)

#Fitting Logistic Model To Data For First 120 Data Points
East_Africa$PCR_Percentage <- East_Africa$PCR_Percentage+0.0001
East_Africa$LM_Percentage <- East_Africa$LM_Percentage+0.0001
Asia$PCR_Percentage <- Asia$PCR_Percentage+0.0001
Asia$LM_Percentage <- Asia$LM_Percentage+0.0001
South_America$PCR_Percentage <- South_America$PCR_Percentage+0.0001
South_America$LM_Percentage <- South_America$LM_Percentage+0.0001

summary(lm(logit(West_Africa$PCR_Percentage) ~ 0 + logit(West_Africa$LM_Percentage)))
summary(lm(logit(East_Africa$PCR_Percentage) ~ 0 + logit(East_Africa$LM_Percentage)))
summary(lm(logit(Asia$PCR_Percentage) ~ 0 + logit(Asia$LM_Percentage)))
summary(lm(logit(South_America$PCR_Percentage) ~ 0 + logit(South_America$LM_Percentage)))

LM <- seq(0,0.2,0.001)
logit_LM <- logit(LM)

#Plotting on the logit scale
plot(logit(South_America$PCR_Percentage), logit(South_America$LM_Percentage), xlim = c(-10, 0), ylim = c(-10,0), pch = 20, col = "blue")
lines(0.59699*logit_LM, logit_LM, col = "blue", lwd = 2)
points(logit(West_Africa$PCR_Percentage), logit(West_Africa$LM_Percentage), xlim = c(-10, 0), ylim = c(-10,0), pch = 20, col = "red")
lines(0.80607*logit_LM, logit_LM, col = "red", lwd = 2)
points(logit(East_Africa$PCR_Percentage), logit(East_Africa$LM_Percentage), xlim = c(-10,0), ylim = c(-10,0), pch = 20, col = "green")
lines(0.51459*logit_LM, logit_LM, col = "green", lwd = 2)
points(logit(Asia$PCR_Percentage), logit(Asia$LM_Percentage), xlim = c(-10,0), ylim = c(-10,0), pch = 20, col = "black")
lines(0.59204*logit_LM, logit_LM, col = "black", lwd = 2)

#Plotting on the real scale- note the issue with the green line amongst others, as in, it's quite clearly not passing through
#in such a way that accurately reflects the underlying points.
plot(South_America$PCR_Percentage, South_America$LM_Percentage, xlim = c(0, 0.2), ylim = c(0,0.2), pch = 20, col = "blue")
lines(expit(0.59699*logit_LM), expit(logit_LM), col = "blue", lwd = 2)
points(West_Africa$PCR_Percentage, West_Africa$LM_Percentage, xlim = c(0, 0.2), ylim = c(0,0.2), pch = 20, col = "red")
lines(expit(0.80607*logit_LM), expit(logit_LM), col = "red", lwd = 2)
points(East_Africa$PCR_Percentage, East_Africa$LM_Percentage, xlim = c(0,0.2), ylim = c(0,0.2), pch = 20, col = "green")
lines(expit(0.51459*logit_LM), expit(logit_LM), col = "green", lwd = 2)
points(Asia$PCR_Percentage, Asia$LM_Percentage, xlim = c(0,0.2), ylim = c(0,0.2), pch = 20, col = "black")
lines(expit(0.59204*logit_LM), expit(logit_LM), col = "black", lwd = 2)

#Fitting model with square root transformed variables
summary(lm(sqrt(West_Africa$PCR_Percentage) ~ 0 + sqrt(West_Africa$LM_Percentage)))
summary(lm(sqrt(East_Africa$PCR_Percentage) ~ 0 + sqrt(East_Africa$LM_Percentage)))
summary(lm(sqrt(Asia$PCR_Percentage) ~ 0 + sqrt(Asia$LM_Percentage)))
summary(lm(sqrt(South_America$PCR_Percentage) ~ 0 + sqrt(South_America$LM_Percentage)))

LM <- seq(0,0.2,0.001)
sqrt_LM <- sqrt(LM)

#Plotting on the square rooted scale
plot(sqrt(South_America$PCR_Percentage), sqrt(South_America$LM_Percentage), xlim = c(0, 0.5), ylim = c(0, 0.5), pch = 20, col = "blue")
lines(2.2137*sqrt_LM, sqrt_LM, col = "blue", lwd = 2)
points(sqrt(West_Africa$PCR_Percentage), sqrt(West_Africa$LM_Percentage), xlim = c(0, 0.5), ylim = c(0, 0.5), pch = 20, col = "red")
lines(1.21180*sqrt_LM, sqrt_LM, col = "red", lwd = 2)
points(sqrt(East_Africa$PCR_Percentage), sqrt(East_Africa$LM_Percentage), xlim = c(0, 0.5), ylim = c(0, 0.5), pch = 20, col = "green")
lines(1.4943*sqrt_LM, sqrt_LM, col = "green", lwd = 2)
points(sqrt(Asia$PCR_Percentage), sqrt(Asia$LM_Percentage), xlim = c(0, 0.5), ylim = c(0, 0.5), pch = 20, col = "black")
lines(1.34205*sqrt_LM, sqrt_LM, col = "black", lwd = 2)

#Plotting on the real scale- note the issue with the green line amongst others, as in, it's quite clearly not passing through
#in such a way that accurately reflects the underlying points.
plot(South_America$PCR_Percentage, South_America$LM_Percentage, xlim = c(0, 0.2), ylim = c(0,0.2), pch = 20, col = "blue")
lines((2.2137*sqrt_LM)^2, sqrt_LM^2, col = "blue", lwd = 2)
points(West_Africa$PCR_Percentage, West_Africa$LM_Percentage, xlim = c(0, 0.2), ylim = c(0,0.2), pch = 20, col = "red")
lines((1.21180*sqrt_LM)^2, sqrt_LM^2, col = "red", lwd = 2)
points(East_Africa$PCR_Percentage, East_Africa$LM_Percentage, xlim = c(0,0.2), ylim = c(0,0.2), pch = 20, col = "green")
lines((1.4943*sqrt_LM)^2, sqrt_LM^2, col = "green", lwd = 2)
points(Asia$PCR_Percentage, Asia$LM_Percentage, xlim = c(0,0.2), ylim = c(0,0.2), pch = 20, col = "black")
lines((1.34205*sqrt_LM)^2, sqrt_LM^2, col = "black", lwd = 2)

#Fitting model with square root transformed variables BUT with LM_Percentage as the independent variable
summary(lm(sqrt(West_Africa$LM_Percentage) ~ 0 + sqrt(West_Africa$PCR_Percentage)))
summary(lm(sqrt(East_Africa$LM_Percentage) ~ 0 + sqrt(East_Africa$PCR_Percentage)))
summary(lm(sqrt(Asia$LM_Percentage) ~ 0 + sqrt(Asia$PCR_Percentage)))
summary(lm(sqrt(South_America$LM_Percentage) ~ 0 + sqrt(South_America$PCR_Percentage)))

PCR <- seq(0,0.2,0.001)
sqrt_PCR <- sqrt(PCR)

#Plotting on the square rooted scale
plot(sqrt(South_America$PCR_Percentage), sqrt(South_America$LM_Percentage), xlim = c(0, 0.5), ylim = c(0, 0.5), pch = 20, col = "blue")
lines(sqrt_PCR, 0.34972*sqrt_PCR,  col = "blue", lwd = 2)
points(sqrt(West_Africa$PCR_Percentage), sqrt(West_Africa$LM_Percentage), xlim = c(0, 0.5), ylim = c(0, 0.5), pch = 20, col = "red")
lines(sqrt_PCR, 0.8073*sqrt_PCR, col = "red", lwd = 2)
points(sqrt(East_Africa$PCR_Percentage), sqrt(East_Africa$LM_Percentage), xlim = c(0, 0.5), ylim = c(0, 0.5), pch = 20, col = "green")
lines(sqrt_PCR, 0.52634*sqrt_PCR, col = "green", lwd = 2)
points(sqrt(Asia$PCR_Percentage), sqrt(Asia$LM_Percentage), xlim = c(0, 0.5), ylim = c(0, 0.5), pch = 20, col = "black")
lines(sqrt_PCR, 0.67285*sqrt_PCR, col = "black", lwd = 2)

#Plotting on the real scale- note the issue with the green line amongst others, as in, it's quite clearly not passing through
#in such a way that accurately reflects the underlying points.
plot(South_America$PCR_Percentage, South_America$LM_Percentage, xlim = c(0, 0.2), ylim = c(0,0.2), pch = 20, col = "blue")
lines(sqrt_LM^2, (0.34972*sqrt_LM)^2,  col = "blue", lwd = 2)
points(West_Africa$PCR_Percentage, West_Africa$LM_Percentage, xlim = c(0, 0.2), ylim = c(0,0.2), pch = 20, col = "red")
lines(sqrt_LM^2, (0.8073*sqrt_LM)^2, col = "red", lwd = 2)
points(East_Africa$PCR_Percentage, East_Africa$LM_Percentage, xlim = c(0,0.2), ylim = c(0,0.2), pch = 20, col = "green")
lines(sqrt_LM^2, (0.52634*sqrt_LM)^2,  col = "green", lwd = 2)
points(Asia$PCR_Percentage, Asia$LM_Percentage, xlim = c(0,0.2), ylim = c(0,0.2), pch = 20, col = "black")
lines(sqrt_LM^2, (0.67285*sqrt_LM)^2, col = "black", lwd = 2)
x <- seq(0,0.2,0.001)
y <- seq(0,0.2,0.001)
lines(x,y)
