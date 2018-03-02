# Access relevant packages
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)

# Assign whichever dataset being used to data_frame
data_frame <- as.data.frame(Data.4th.February)

# Manipulate data_frame in prep for analysis
data_frame$Age_Group_2 <- as.factor(data_frame$Age_Group_2)

# Create individual datasets for each global region
young_children <- data_frame[data_frame$Age_Group_2 == "0-5", ]
old_children <- data_frame[data_frame$Age_Group_2 == "5-15", ]
adults <- data_frame[data_frame$Age_Group == "15+",]

# Without the zeroes and the associated transformations of them
young_children <- data_frame[data_frame$Age_Group_2 == "0-5" & data_frame$Was_Initially_Zero. == "N", ]
old_children <- data_frame[data_frame$Age_Group_2 == "5-15" & data_frame$Was_Initially_Zero. == "N", ]
adults <- data_frame[data_frame$Age_Group == "15+" & data_frame$Was_Initially_Zero. == "N", ]

# Running the linear regression, with PCR regressed onto Microscopy. 









