library(ggplot2)
library(reshape2)
library(dplyr)
library(drc)

## Contains functinos for converting strings into time.
# source("TimeConversion.R")

## Contains functions for converting the Combinations into placements, and vice versa.
## convertToPlacement(x) where x is the name of the combination
## convertToPot(x) ...
## convertToCamera(x) ...
source("CameraConversion.R")

## Supporting functions that make everything easier.
## mergeVectors(x, y) Merge 2 string vectors by a seperator "_" and returning z
## divide_by_harverst(df, combination_table) Takes a data.frame (df) and uses a combination_table which contains the harvest information
## pickone(x) Takes the first element from vector x, assuming all values in this vector are identical.
## predict_points(df) Takes the TAf table and predicts points linearly between points if there is more than 1 day between them.
source("supportfunctions.R")

##############################################
## Compiling Harvest information
##############################################

harvested <- read.csv("nchain_harvested.tab.csv", sep=";")

full_table <- read.table("clover.tab", header=TRUE, sep="\t")
full_table <- full_table[!is.na(full_table$ID),]
full_table <- full_table[!is.na(full_table$Table),]
rownames(full_table) <- full_table$ID
#full_table$ID <- NULL
full_table$Cutting <- NULL
full_table$Harvest <- NULL
full_table$Potting <- NULL
full_table$Status <- NULL
full_table %>% arrange(ID) -> full_table
full_table$Combination <- mergeVectors(full_table$Clover, full_table$Rhizobium)

HarvestData <- read.csv("harvest1.csv", sep=";")
names(HarvestData)[1] <- "barcode"
HarvestData$Combination <- full_table[HarvestData$barcode,"Combination"]
rownames(HarvestData) <- HarvestData$barcode
HarvestData$harvest_date <- as.Date(HarvestData$harvest_date, format="%d/%m/%y")
HarvestData$cutting_date <- as.Date(HarvestData$cutting_date, format="%d/%m/%y")
HarvestData$inoculation_date <- as.Date(HarvestData$inoculation_date, format="%d/%m/%y")
HarvestData$potting_date <- as.Date(HarvestData$potting_date, format="%d/%m/%y")
HarvestData <- HarvestData[as.character(sort(HarvestData$barcode)),]
rownames(HarvestData) <- HarvestData$Combination

#full_table <- full_table[sort(harvested$Pot), ]
full_table <- full_table[sort(HarvestData$barcode), ]

#row.names(harvested) <- harvested$Pot
#harvested$X.2 <- NULL
#harvested <- harvested[as.character(sort(as.numeric(rownames(harvested)))),]

#full_table$Harvest <- as.Date(harvested$date, format="%d/%m/%Y")
full_table$Harvest <- as.Date(HarvestData$harvest_date, format="%d/%m/%y")

combination_table <- full_table
row.names(combination_table) <- combination_table$Combination

##############################################
## Loading the growth data
##############################################

TA <- read.csv("CloverImageDataTA.csv")

all_names <- colnames(TA)[2:ncol(TA)]

## Filter by harvested
TA <- TA[,c("Time", full_table$Combination)]

## Convert into Date.
TA$Time <- as.Date(TA$Time, format="%Y/%m/%d - %H:%M")


## melt from matrix into continuous data.frame
TA <- melt(TA, id.vars = "Time", variable.name="Pot", value.name = "TotalArea")

## Convert pixels into cm^2
max_pixel_count <- 400*400
crop_in_plate_size <- 40*40

TA$TotalArea <- TA$TotalArea/max_pixel_count
TA$TotalArea <- TA$TotalArea*crop_in_plate_size

## Filter early days.
TA %>% filter(Time!="2017-05-23") %>% filter(Time!="2017-05-22") -> TA

## Add Camera ID and Pot ID to the data.frame
TA %>% mutate(Camera = unlist(convertToCamera(Pot)),
              Placement = unlist(convertToPot(Pot))) -> TA

## Remove all rows with missing data.
TA <- TA[!is.na(TA$TotalArea),]

## Create per day summary data
TA %>% group_by(Pot, Time) %>%
    summarise(Std=sd(TotalArea), TotalArea=median(TotalArea), Camera=pickone(Camera),
              Placement=pickone(Placement)) -> TAf

## Convert into data.frame
#TAf <- as.data.frame(TAf)

## Remove missing data if any.
TAf <- TAf[!is.na(TAf$TotalArea),]
TAf <- TAf[!is.na(TAf$Std),]

## Divided by harvest
TAf$Period <- divide_by_harvest(TAf, combination_table)

show <- function(Period){
    length(unique(Period))>1
}

TAf %>% group_by(Pot) %>% filter(show(Period)) -> TAf

##############################################
## Filtering days (Done manually)
##############################################

## Contains the manual filtering list and black list.
## manual_filter(pot, time) Uses the manual_list to filter if that day is in the list.
## blacklist_filter(pot) If pot is in black_list then filter it.
source('manual_list.R')

TAf %>% group_by(Pot, Time) %>% filter(blacklist_filter(Pot)) -> TAblacklist

in_filter <- Vectorize(function(Pot, Manual_list) {
    exists(as.character(Pot), where=Manual_list)
}, vectorize.args = "Pot")

TAf %>% group_by(Pot, Time) %>% filter(!blacklist_filter(Pot)) %>%
    filter(manual_filter(Pot, Time)) -> TAf

## Add predicted results
TAfpred <- predict_points(TAf)

TAfpred %>% group_by(Pot, Time) -> TAf

############################
# Getting the growth rates #
#####################################
## Get growth rates using the medians
#####################################

TAf %>% group_by(Pot, Period) %>%
    summarise(b=drm(TotalArea ~ Time, fct=L.4())$coefficients[1],
              c=drm(TotalArea ~ Time, fct=L.4())$coefficients[2],
              d=drm(TotalArea ~ Time, fct=L.4())$coefficients[3],
              e=drm(TotalArea ~ Time, fct=L.4())$coefficients[4]) -> GrowthRates

GrowthRates %>% arrange(b) -> GrowthRates

GrowthRates <- as.data.frame(GrowthRates)
#print(GrowthRates)
