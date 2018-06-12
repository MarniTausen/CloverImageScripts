library(zoo)
library(dplyr)

#GrowthRates <- read.csv("GrowthRates.csv")
HarvestData <- read.csv("harvest1.csv", sep=";")

names(HarvestData)[1] <- "barcode"
HarvestData$clover <- as.character(HarvestData$clover)
HarvestData$rhizobium <- as.character(HarvestData$rhizobium)

HarvestData$harvest_date <- as.Date(HarvestData$harvest_date, format="%d/%m/%y")

GrowthRates$X <- NULL

names(GrowthRates) <- c("Pot", "Period", "slope", "start", "end", "turn_numeric")

GrowthRates$turn_date <- as.Date(GrowthRates$turn_numeric)
GrowthRates$TotalGrowth <- GrowthRates$end-GrowthRates$start

GrowthRates$Ratio <- GrowthRates$TotalGrowth*((GrowthRates$slope)^2)/2

GrowthRates %>% arrange(Pot) -> GrowthRates

identify_replicates <- Vectorize(function(x, pattern) regexpr(paste(pattern,"[.]",sep=""),
                                                              paste(x,".",sep=""))[1]>0)

replicate_consistency <- function(df){
    replicates <- list()
    for(name in unique(df$Pot)){
        df %>% filter(identify_replicates(df$Pot, name)) -> replicates[[name]]
    }

    summary_df <- data.frame(Pot=vector(), Slope_difference=vector(),
                             Growth_difference=vector())


    for(replicate in replicates){
        name <- sort(unique(replicate$Pot))[1]
        if(nrow(replicate)>2){
            replicate %>% filter(Period=="growth1") -> G1
            slope_diff = sqrt((G1$slope[1]-G1$slope[2])^2)
            growth_diff = sqrt((G1$TotalGrowth[1]-G1$TotalGrowth[2])^2)
            tmp <- data.frame(Pot=name, Slope_difference=slope_diff,
                              Growth_difference=growth_diff)
            summary_df <- rbind(summary_df, tmp)
        } else {
            tmp <- data.frame(Pot=name, Slope_difference=NA, Growth_difference=NA)
            summary_df <- rbind(summary_df, tmp)
        }

    }
    summary_df
}

Replicate_statistics <- replicate_consistency(GrowthRates)

missing_replicates <- Replicate_statistics[is.na(Replicate_statistics$Slope_difference),]

cat("Number of pots with no replicates:", nrow(missing_replicates), "\n")

Replicate_statistics <- Replicate_statistics[!is.na(Replicate_statistics$Slope_difference),]

cat("Number of pots with replicates:", nrow(Replicate_statistics), "\n")

Replicate_statistics %>% arrange(Slope_difference) -> Replicate_statistics

print(Replicate_statistics)

## Comparing the harvest data, and size at harvest.

GrowthRates %>% filter(Period=="growth1") -> BeforeHarvest

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

HarvestData$Pot <- full_table[HarvestData$barcode, "Combination"]
HarvestData$Combination <- HarvestData$Pot
rownames(HarvestData) <- HarvestData$Combination

BeforeHarvest$Pot <- as.character(BeforeHarvest$Pot)

test <- inner_join(HarvestData, BeforeHarvest, by=c("Pot"))

cat("Number of Pots where there is growth data:", nrow(BeforeHarvest), "\n")

cat("Number of Pots in the HarvestData", nrow(HarvestData), "\n")

cat("Number of Pots which overlap between them", nrow(test), "\n")

match <- inner_join(HarvestData, full_table, by=c("Combination"))

cat("Number of Pots in the harvest table:", nrow(match), "\n")

get_harvest_days <- function(Time, Pot){
    if(is.na(HarvestData[as.character(Pot), "harvest_date"])) return(FALSE)
    if(Time==(HarvestData[as.character(Pot), "harvest_date"]-1)) return(TRUE)
    else return(FALSE)
}

TAf %>% filter(get_harvest_days(Time, Pot)) %>%
    mutate(weight=HarvestData[Pot, "weight"],
           growth_day=HarvestData[Pot, "growth_day"]) -> TAfweight

multi_effective_strains <- as.character(read.table("multi_effective_strains.csv", header=FALSE)$V1)
multi_effective_strains[1] <- unlist(strsplit(multi_effective_strains[1], split="\277"))[2]

for(strain in multi_effective_strains){
    TAf %>% filter(get_harvest_days(Time, Pot)) %>%
        filter(identify_replicates(Pot, strain)) %>%
        mutate(weight=HarvestData[Pot, "weight"],
               growth_day=HarvestData[Pot, "growth_day"]) -> TAfweight
    figure <- qplot(TotalArea, weight, data=TAfweight, geom="point") + theme_classic() +
        labs(title=strain)
    print(figure)
}

#qplot(TotalArea, weight, data=TAfweight, geom="point") + theme_classic()
#qplot(TotalArea, growth_day, data=TAfweight, geom="point") + theme_classic()
#qplot(weight, growth_day, data=TAfweight, geom="point") + theme_classic()

#In manual list
#for(strain in multi_effective_strains){
#    cat("Strain:", strain, "\n")
#    print(length(names(Manual_list)[identify_replicates(names(Manual_list), strain)]))
#}

#In black list
#for(strain in multi_effective_strains){
#    cat("Strain:", strain, "\n")
#    print(Black_list[identify_replicates(Black_list, strain)])
#}

#In harvest data
#for(strain in multi_effective_strains){
#    cat("Strain:", strain, "\n")
#    print(c(length(HarvestData$Pot[identify_replicates(HarvestData$Pot, strain)]),
#            length(names(Manual_list)[identify_replicates(names(Manual_list), strain)]),
#            length(names(Manual_list)[identify_replicates(names(Manual_list), strain)])+
#            length(Black_list[identify_replicates(Black_list, strain)]),
#            length(all_names[identify_replicates(all_names, strain)])))
#}


full_harvest_table <- read.table("20171006_clover_harvest.tab", sep="\t")
colnames(full_harvest_table) <- c("barcode", "harvest_id", "clover", "rhizobium",
                                  "table", "placement", "harvest", "cutting_date",
                                  "potting_date", "inoculation_date", "harvest_date",
                                  "weight", "growth_day")


full_harvest_table$Pot <- full_table[full_harvest_table$barcode, "Combination"]
full_harvest_table <- full_harvest_table[!is.na(full_harvest_table$weight),]

full_harvest_table$clover <- as.character(full_harvest_table$clover)
full_harvest_table$rhizobium <- as.character(full_harvest_table$rhizobium)
full_harvest_table$cutting_date <- as.Date(full_harvest_table$cutting_date,
                                           format="%d/%m/%Y")
full_harvest_table$inoculation_date <- as.Date(full_harvest_table$inoculation_date,
                                               format="%d/%m/%Y")
full_harvest_table$potting_date <- as.Date(full_harvest_table$potting_date,
                                           format="%d/%m/%Y")
full_harvest_table$harvest_date <- as.Date(full_harvest_table$harvest_date,
                                           format="%d/%m/%Y")

full_harvest_table %>% arrange(barcode) %>% group_by(barcode)

full_harvest_table <- full_harvest_table[!is.na(full_harvest_table$Pot),]

#In harvest data
#for(strain in multi_effective_strains){
#    cat("Strain:", strain, "\n")
#    print(unique(full_harvest_table$Pot)[identify_replicates(
#                    unique(full_harvest_table$Pot), strain)])
#}

#for(strain in multi_effective_strains){
#    cat("Strain:", strain, "\n")
#    print(unique(HarvestData$Pot)[identify_replicates(
#                    unique(HarvestData$Pot), strain)])
#}

findequals <- function(vx, vy){
    matches <- 0
    for(x in vx){
        if(any(x==vy)){
            matches <- matches+1
            #cat(x, "in", deparse(substitute(vy)), "\n")
        } else {
            cat(x, "not in", deparse(substitute(vy)), "\n")
        }
    }
    print(matches)
    print(length(vx))
}

#findequals(HarvestData$Pot, full_harvest_table$Pot)

#findequals(HarvestData$barcode, full_harvest_table$barcode)

