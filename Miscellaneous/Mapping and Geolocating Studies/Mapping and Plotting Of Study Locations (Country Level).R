# Plotting the locations of the studies
library(rworldmap)

# Plotting by country and number of studies in each country


e2020_Countries <- E2020Countries
e2020_Countries <- e2020_Countries[, c("Country_ID", "Category")]

for (i in 1:length(e2020_Countries$Category)) {
  
  e2020_Countries$Category[i] <- rnorm(1, 20, 5)
  
}


sPDF <- joinCountryData2Map(e2020_Countries, joinCode = "ISO3", nameJoinColumn = "Country_ID" )

mapCountryData(sPDF, nameColumnToPlot="Category", catMethod = c(10, 15, 20, 25, 30), colourPalette = "terrain")

?terrain.colors

?mapCountryData

# Plotting by country and number of datapoints in each country

sPDF <- joinCountryData2Map( Number_of_studies_by_country,
                             , joinCode = "ISO3"
                             , nameJoinColumn = "Country_ID" )

mapCountryData( sPDF, nameColumnToPlot="Number_Of_Datapoints", numCats = 4,  
                catMethod = c(0, 5, 10, 15, 60), 
                colourPalette = c("#F9E500", "#FFB200", "#FF4C00" 
                                  , "#F40400"))
