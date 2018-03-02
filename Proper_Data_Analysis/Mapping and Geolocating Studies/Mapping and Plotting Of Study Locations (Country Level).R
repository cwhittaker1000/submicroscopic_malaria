# Plotting the locations of the studies
library(rworldmap)

# Plotting by country and number of studies in each country

sPDF <- joinCountryData2Map( Number_of_studies_by_country,
                             , joinCode = "ISO3"
                             , nameJoinColumn = "Country_ID" )

mapCountryData( sPDF, nameColumnToPlot="Number_Of_Studies",  numCats = 9, 
                catMethod = c(1,2,3,4,5,6,7,8,9))


# Plotting by country and number of datapoints in each country

sPDF <- joinCountryData2Map( Number_of_studies_by_country,
                             , joinCode = "ISO3"
                             , nameJoinColumn = "Country_ID" )

mapCountryData( sPDF, nameColumnToPlot="Number_Of_Datapoints", numCats = 4,  
                catMethod = c(0, 5, 10, 15, 60), 
                colourPalette = c("#F9E500", "#FFB200", "#FF4C00" 
                                  , "#F40400"))




