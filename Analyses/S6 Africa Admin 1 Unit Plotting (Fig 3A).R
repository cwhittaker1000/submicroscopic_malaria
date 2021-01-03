# Access and load various required packages for the analyses
library(rjags); library(ssa); library(binom); library(raster); library(sf); library(maptools); library(countrycode);
library(dplyr)
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/Sub-Patent Malarial Infections/Sub_Patent_Malaria_Analysis/")

# Loading In the Dataset
data_frame <- read.csv("Data/SI_Systematic_Review_Results_R_Import.csv")
data_frame <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2 & !is.na(data_frame$Transmission_Setting_15), ]
plotting_df <- data_frame[, c("Country", "Admin_1", "Transmission_Setting_15")]
countries <- unique(as.character(data_frame$Country))
admin_ones <- as.character(data_frame$Admin_1)
unique_admin_ones <- unique(as.character(data_frame$Admin_1))

indices <- seq(1:21)
countries[indices]

indices <- indices[-c(4, 5)] # removing Kenya and Malawi as I'll have to process those manually
countries[indices]

par(mfrow = c(1, 1))
polygons_list <- list()
transmission_archetype <- c()
country <- c()
counter <- 1

for (i in indices) {
  country_specific_admin_ones <- getData('GADM', country = countries[i], level = 1)
  country_specific_dataset <- data_frame[data_frame$Country == countries[i], ]
  survey_present_admin_ones <- unique(as.character(country_specific_dataset$Admin_1))
  for (j in 1:length(survey_present_admin_ones)) {
    index <- grep(paste0("^", survey_present_admin_ones[j], "$"), country_specific_admin_ones$NAME_1)
    if (i == 1 & j == 1) {
      admin_one_spatial_polygons_dataframe <- country_specific_admin_ones[index, ]
      admin_unit_specific_dataset <- country_specific_dataset[country_specific_dataset$Admin_1 == survey_present_admin_ones[j], ]
      transmission_archetype[counter] <- as.character(admin_unit_specific_dataset$Transmission_Setting_15[1])
      country[counter] <- countries[i]
      print(counter)
      counter <- counter + 1
    }
    else {
      admin_one_spatial_polygons_dataframe <- rbind(admin_one_spatial_polygons_dataframe, country_specific_admin_ones[index, ])
      admin_unit_specific_dataset <- country_specific_dataset[country_specific_dataset$Admin_1 == survey_present_admin_ones[j], ]
      transmission_archetype[counter] <- as.character(admin_unit_specific_dataset$Transmission_Setting_15[1])
      country[counter] <- countries[i]
      counter <- counter + 1
    }
  }
  print(c(i, countries[i]))
}
plot(admin_one_spatial_polygons_dataframe[28, ])
transmission_archetype
country

BF <- which(country == "Burkina Faso")
transmission_archetype[BF]

which(country == "Burkina Faso")
which(admin_one_spatial_polygons_dataframe$GID_0 == "BFA")

unique(country)
unique(admin_one_spatial_polygons_dataframe$GID_0)

# sorting Kenya
Kenya_admin_units <- c("Mombasa", "Kwale", "Kilifi", "Tana River", "Lamu", "Taita Taveta", "Siaya", "Kisumu", 
                       "Migori", "Homa Bay", "Kisii", "Nyamira", "Kakamega", "Vihiga", "Bungoma", "Busia", 
                       "Turkana", "West Pokot", "Samburu", "Trans Nzoia", "Uasin Gishu", "Elgeyo-Marakwet", 
                       "Nandi", "Baringo", "Laikipia", "Nakuru", "Narok", "Kajiado", "Kericho", "Bomet")
Kenya_admin_shapefiles <- getData('GADM', country = "Kenya", level = 1)
for (i in 1:length(Kenya_admin_units)){
  
  country_specific_dataset <- data_frame[data_frame$Country == "Kenya", ]
  index <- grep(paste0("^", Kenya_admin_units[i], "$"), Kenya_admin_shapefiles$NAME_1)
  admin_one_spatial_polygons_dataframe <- rbind(admin_one_spatial_polygons_dataframe, Kenya_admin_shapefiles[index, ])
  
  if (i <= 6) {
    admin_unit_specific_dataset <- country_specific_dataset[country_specific_dataset$Admin_1 == "Coast", ]
    transmission_archetype[counter] <- as.character(admin_unit_specific_dataset$Transmission_Setting_15[1])
    country[counter] <- "Kenya"
    counter <- counter + 1
  }
  else if (i > 6 & i <= 12) {
    admin_unit_specific_dataset <- country_specific_dataset[country_specific_dataset$Admin_1 == "Nyanza", ]
    transmission_archetype[counter] <- as.character(admin_unit_specific_dataset$Transmission_Setting_15[1])
    country[counter] <- "Kenya"
    counter <- counter + 1
  }
  else if (i > 12 & i <= 16) {
    admin_unit_specific_dataset <- country_specific_dataset[country_specific_dataset$Admin_1 == "Western", ]
    transmission_archetype[counter] <- as.character(admin_unit_specific_dataset$Transmission_Setting_15[1])
    country[counter] <- "Kenya"
    counter <- counter + 1
  } 
  else if (i > 16) {
    admin_unit_specific_dataset <- country_specific_dataset[country_specific_dataset$Admin_1 == "Rift Valley", ]
    transmission_archetype[counter] <- as.character(admin_unit_specific_dataset$Transmission_Setting_15[1])
    country[counter] <- "Kenya"
    counter <- counter + 1
  }
}
admin_one_spatial_polygons_dataframe

# sorting Malawi
Malawi_admin_units <- c("Dedza", "Dowa", "Kasungu", "Lilongwe", "Mchinji", "Nkhotakota", "Ntcheu", "Ntchisi", "Salima", 
                       "Balaka", "Blantyre", "Chikwawa", "Chiradzulu", "Machinga", "Mangochi", "Mulanje", "Mwanza", "Nsanje", "Thyolo", "Phalombe", "Zomba", "Neno")
Malawi_admin_shapefiles <- getData('GADM', country = "Malawi", level = 1)
for (i in 1:length(Malawi_admin_units)){

  country_specific_dataset <- data_frame[data_frame$Country == "Malawi", ]
  index <- grep(paste0("^", Malawi_admin_units[i], "$"), Malawi_admin_shapefiles$NAME_1)
  admin_one_spatial_polygons_dataframe <- rbind(admin_one_spatial_polygons_dataframe, Malawi_admin_shapefiles[index, ])
  
  if (i <= 9) {
    admin_unit_specific_dataset <- country_specific_dataset[country_specific_dataset$Admin_1 == "Central Region", ]
    transmission_archetype[counter] <- as.character(admin_unit_specific_dataset$Transmission_Setting_15[1])
    country[counter] <- "Malawi"
    counter <- counter + 1
  }
  else if (i > 10) {
    admin_unit_specific_dataset <- country_specific_dataset[country_specific_dataset$Admin_1 == "Southern Region", ]
    transmission_archetype[counter] <- as.character(admin_unit_specific_dataset$Transmission_Setting_15[1])
    country[counter] <- "Malawi"
    counter <- counter + 1
  }
}
admin_one_spatial_polygons_dataframe$NAME_1
transmission_archetype
country

BF <- which(country == "Burkina Faso")
transmission_archetype[BF]

list1 <- countrycode::codelist_panel %>% 
  filter(continent == "Africa")
country_codes <- unique(list1[, "wb"])
country_codes <- c(country_codes[-44], "EH") 
for (i in 1:length(country_codes)) {
  if (i == 1) {
    country_shapefile <- getData('GADM', country = country_codes[i], level = 0)
  }
  else {
    new_country_shapefile <- getData('GADM', country = country_codes[i], level = 0)
    country_shapefile <- rbind(country_shapefile, new_country_shapefile)  
  }
  print(i)
}

pdf("Figures/Figure 3 - Transmission Archetype/Figure_3A_Admin_1_Mapping.pdf", width = 12.5, height = 12.5)
palette(c("#00A600FF", "#ECB176FF", "darkgrey"))
plot(country_shapefile)
trans_arch <- as.factor(transmission_archetype)
for (i in 1:length(admin_one_spatial_polygons_dataframe)) {
  plot(admin_one_spatial_polygons_dataframe[i, ], lwd = 2, add = TRUE, col = trans_arch[i])
}
dev.off()

