library(lme4)

TAarea$weight <- HarvestData[TAarea$Pot, "weight"]
TAarea$growthday <- HarvestData[TAarea$Pot, "growth_day"]

TAarea$Camera <- unlist(convertToCamera(TAarea$Pot))
TAarea$Placement <- unlist(convertToPot(TAarea$Pot))

combination_split <- Vectorize(function(x, i) unlist(strsplit(x, split=".", fixed=TRUE))[i],
                               vectorize.args = "x")

is_replicate <- function(x) c(2, 1)[is.na(combination_split(x, 3))+1]

join <- Vectorize(function(x, y) paste(x, y, sep="."))

TAarea$Clover <- combination_split(TAarea$Pot, 1)
TAarea$Rhizobium <- combination_split(TAarea$Pot, 2)
TAarea$Combination <- join(TAarea$Clover, TAarea$Rhizobium)
TAarea$Replicate <- is_replicate(TAarea$Pot)

TAarea$Inoculationdate <- HarvestData[TAarea$Pot, "inoculation_date"]
TAarea$Harvestdate <- HarvestData[TAarea$Pot, "harvest_date"]


TAarea$Period <- NULL

Get_originals <- function(Name){
    sequence <- unlist(strsplit(as.character(Name), split=""))
    n <- length(sequence)
    if(all(sequence[(n-1):n]==c(".", "1"))) return(FALSE)
    else return(TRUE)
}

Get_originals <- Vectorize(Get_originals)

Get_replicate <- function(Name){
    paste(Name, ".1", sep="")
}

any_match <- function(x, y){
    v <- vector(length=length(x))
    j <-
    for(i in 1:length(x)){
        if(any(x[i]==y)) v[i] <- TRUE
        else v[i] <- FALSE
    }
    v
}

## Score of 1 perfect consitency, 1-0 (varying degrees of constistency)
## <0 inverse consistency, meaning the replicates are worse than the average.
replicate_consistency <- function(dataset, measure){
    x <- dataset[[measure]]
    x <- x/max(x)
    originals <- TAarea$Pot[Get_originals(as.character(dataset$Pot))]
    within_variance <- 0
    for(original in originals){
        replicate <- Get_replicate(as.character(original))
        pair <- c(original, replicate)
        new_variance <- var(x[any_match(dataset$Pot, pair)])
        if(!is.na(new_variance)) within_variance <- within_variance+new_variance
    }
    within_variance <- within_variance/length(originals)
    total_variance <- var(x)
    cat("Replicate variance:", within_variance, "\n")
    cat("Total variance:", total_variance, "\n")
    (total_variance-within_variance)/total_variance
}

replicate_consistency(TAarea, "weight")
replicate_consistency(TAarea, "AreaRatio")
replicate_consistency(TAarea, "growthday")
replicate_consistency(TAarea, "Area")

randint <- function(min, max) round(runif(1, min, max))


find_minimum_time_frame <- function(Harvestdates, Inoculationdates){
    as.numeric(min(Harvestdates-Inoculationdates))
}


normalize_weight <- function(growthrates, timeframe){
    growthrates*timeframe
}

# Just transforms it into the growthday measurement, scaled by time.
t <- find_minimum_time_frame(TAarea$Harvestdate, TAarea$Inoculationdate)
TAarea$normweight <- normalize_weight(TAarea$growthday, t)
replicate_consistency(TAarea, "normweight")

normalized_replicate_time <- (TAf %>% filter(is_or_has_replicate(Pot, all_names), Period=="growth1") %>% group_by(Pot) %>% summarise(Days=sum(select_time_frame_filter(Pot, TotalArea, Time))))

normalized_replicate_time <- as.data.frame(normalized_replicate_time)

rownames(normalized_replicate_time) <- as.character(normalized_replicate_time$Pot)

as.data.frame(normalized_replicate_time %>% arrange(as.character(Pot)))

t <- as.numeric(TAarea$Harvestdate-TAarea$Inoculationdate)
TAarea$normweight <- TAarea$weight*(normalized_replicate_time[TAarea$Pot, "Days"]/t)

replicate_consistency(TAarea, "weight")
replicate_consistency(TAarea, "AreaRatio")
replicate_consistency(TAarea, "growthday")
replicate_consistency(TAarea, "Area")
replicate_consistency(TAarea, "normweight")

multi_effective_strains <- as.character(read.delim("multi_effective_strains.csv", header=FALSE)$V1)
multi_effective_strains[1] <- unlist(strsplit(multi_effective_strains[1], split="\277"))[2]

TAarea %>% filter(any_match(Rhizobium, multi_effective_strains)) %>%
    filter(is_or_has_replicate(Pot, all_names)) -> MEarea

replicate_consistency(MEarea, "weight")
replicate_consistency(MEarea, "AreaRatio")
replicate_consistency(MEarea, "growthday")
replicate_consistency(MEarea, "Area")
replicate_consistency(MEarea, "normweight")

fit <- glm(Area ~ Rhizobium, data=MEarea)
fit2 <- glm(growthday ~ Rhizobium, data=MEarea)

qplot(summary(fit2)$coefficients[,1], summary(fit)$coefficients[,1], xlim=c(-0.1, 0.3))

HarvestData$rhizobium <- as.character(HarvestData$rhizobium)

HarvestData %>% filter(any_match(rhizobium, multi_effective_strains)) -> MEharvest


summary(glm(growth_day ~ rhizobium, data=MEharvest))


fit <- lm(Area ~ Rhizobium, data=MEarea)
fit2 <- lm(growth_day ~ rhizobium, data=MEharvest)

summary(fit)
summary(fit2)

qplot(summary(fit2)$coefficients[c(-1, -21),1], summary(fit)$coefficients[-1,1], xlim=c(-0.1, 0.3))


fit <- lm(Area ~ Rhizobium + Clover, data=MEarea)
fit2 <- lm(growth_day ~ rhizobium + clover, data=MEharvest)

summary(fit)
summary(fit2)

qplot(summary(fit2)$coefficients[c(2:20, 22:23),1], summary(fit)$coefficients[2:22,1])

AreaTable <- as.data.frame(summary(fit)$coefficients[2:22,])
AreaTable$Names <- rownames(summary(fit)$coefficients[2:22,])
AreaTable <- AreaTable[, c('Names', 'Estimate')]

GrowthTable <- as.data.frame(summary(fit2)$coefficients[2:22,])
GrowthTable$Names <- rownames(summary(fit2)$coefficients[2:22,])
GrowthTable <- GrowthTable[, c('Names', 'Estimate')]

MEharvest %>% filter(any_match(Combination, MEarea$Pot)) -> MEharvestsample

MEarea$Rhizobium <- as.factor(MEarea$Rhizobium)
MEarea$Clover <- as.factor(MEarea$Clover)

MEharvest$rhizobium <- as.factor(MEharvest$rhizobium)
MEharvest$clover <- as.factor(MEharvest$clover)

fit <- glm(Area ~ Rhizobium + Clover, data=MEarea)
fit2 <- glm(growth_day ~ rhizobium + clover, data=MEharvest)

summary(fit)
summary(fit2)

summary(glht(fit2, linfct=mcp(rhizobium="Tukey")))

qplot(summary(fit2)$coefficients[c(2:20, 22:23),1], summary(fit)$coefficients[2:22,1])

AreaTable <- as.data.frame(summary(fit)$coefficients[2:22,])
AreaTable$Names <- rownames(summary(fit)$coefficients[2:22,])
AreaTable <- AreaTable[, c('Names', 'Estimate')]

GrowthTable <- as.data.frame(summary(fit2)$coefficients[2:22,])
GrowthTable$Names <- rownames(summary(fit2)$coefficients[2:22,])
GrowthTable <- GrowthTable[, c('Names', 'Estimate')]


qplot(MEarea$growthday, MEarea$Area)

cor(MEarea$growthday, MEarea$Area, method="spearman")

cor(summary(fit2)$coefficients[c(2:20, 22:23),1], summary(fit)$coefficients[2:22,1], method="spearman")


fit <- lm(Area ~ Rhizobium, data=MEarea)

summary(fit)

summary(glht(fit, linfct=mcp(Rhizobium="Tukey")))


MEarea$Clover <- factor(MEarea$Clover)

fit <- lm(growthday ~ Rhizobium + Clover, data=MEarea)

summary(fit)

filtered <- (MEarea %>% filter(!any_match(Clover, c("Aalon_0718", "Aoost_0215",
                                                       "Banna_0733", "Clfin_0213",
                                                       "Banca_0947", "Clfin_0102",
                                                       "Aalon_0617", "Aalon_0512")),
                                  !any_match(Rhizobium, c("SM88"))))

fit <- lm(Area ~ Rhizobium + Clover,
          data=(MEarea %>% filter(!any_match(Clover, c("Aalon_0718", "Aoost_0215",
                                                       "Banna_0733", "Clfin_0213",
                                                       "Banca_0947", "Clfin_0102",
                                                       "Aalon_0617", "Aalon_0512")),
                                  !any_match(Rhizobium, c("SM88")))))


summary(fit)

significance_column <- Vectorize(function(pvalue){
    if(pvalue>0.1) return(' ')
    if(pvalue<=0.1 && pvalue>0.5) return('.')
    if(pvalue<=0.5 && pvalue>0.01) return('*')
    if(pvalue<=0.01 && pvalue>0.001) return('**')
    if(pvalue<=0.001) return('***')
})

rhivrhi <- summary(glht(fit, linfct=mcp(Rhizobium="Tukey")))

rvrtable <- data.frame(estimates=rhivrhi$test$coefficients)
rvrtable$std.error <- rhivrhi$test$sigma
rvrtable$tvalue <- rhivrhi$test$tstat
rvrtable$pvalue <- rhivrhi$test$pvalues
rvrtable[[' ']] <- significance_column(rvrtable$pvalue)
rvrtable <- cbind(rownames(rvrtable), rvrtable)
colnames(rvrtable)[1] <- "Test"

rvrtable %>% filter(pvalue<0.05)



arearc <- lm(Area ~ Rhizobium + Clover, data=MEarea)
areatable <- as.data.frame(summary(arearc)$coefficients)
colnames(areatable)[4] <- "p.value"
areatable <- cbind(rownames(areatable), areatable)
colnames(areatable)[1] <- "Genotypes"


areatable %>% filter(grepl("Rhizobium", Genotypes)) %>% filter(p.value<0.05) %>%
    arrange(-Estimate) %>% mutate(' '=significance_column(p.value))

gpdrc <- lm(growth_day ~ rhizobium + clover, data=MEharvest)
gpdtable <- as.data.frame(summary(gpdrc)$coefficients)
colnames(gpdtable)[4] <- "p.value"
gpdtable <- cbind(rownames(gpdtable), gpdtable)
colnames(gpdtable)[1] <- "Genotypes"

gpdtable %>% filter(grepl("rhizobium", Genotypes)) %>% filter(p.value<0.05) %>%
    arrange(-Estimate) %>% mutate(' '=significance_column(p.value))



## Try Bayz for hertitability

source("http://www.bayz.biz/Rbayz.R")

BAYZHOME = "/usr/local/bin/"

fit <- bayz.mm(data=MEarea, resp="Area", fixmod="fac.Replicate", ranmod="fac.Rhizobium+fac.Clover", chain=c(99900, 1000, 50))

bayz.summ(fit)

fit$post.mean[181]/(fit$post.mean[181]+fit$post.mean[117])

## Replicate versus replicates

(MEarea %>% filter(Replicate==1) %>% arrange(Pot))$Area
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$Area

qplot((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$Area,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$Area) + theme_classic() +
    labs(x="Area (Replicate 1)", y="Area (Replicate 2)")

qplot((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$growthday,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$growthday) + theme_classic() +
    labs(x="Growth per day (Replicate 1)", y="Growth per day (Replicate 2)")

qplot((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$weight,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$weight) + theme_classic() +
    labs(x="Weight (Replicate 1)", y="Weight (Replicate 2)")

qplot((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$AreaRatio,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$AreaRatio) + theme_classic() +
    labs(x="Area ratio (Replicate 1)", y="Area ratio (Replicate 2)")

cat("Area replicate correlation:",
cor((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$Area,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$Area), "\n")

cat("Growth per day replicate correlation:",
cor((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$growthday,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$growthday), "\n")

cat("Weight replicate correlation:",
cor((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$weight,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$weight), "\n")

cat("Weight replicate correlation:",
cor((MEarea %>% filter(Replicate==1) %>% arrange(Pot))$AreaRatio,
(MEarea %>% filter(Replicate==2) %>% arrange(Pot))$AreaRatio), "\n")

## Calculate heritiabilities


calculate_heritability <- function(lmerMod, variable){
    RE <- as.data.frame(VarCorr(lmerMod))
    if(!class(variable)=="character") stop("Variable must be a string of Characters!")
    if(length(variable)==1) {
        return(RE[RE$grp==variable, "vcov"]/sum(RE$vcov))
    }

print(RE)

    any_match <- function(x, y){
        v <- vector(length=length(x))
        j <-
            for(i in 1:length(x)){
                if(any(x[i]==y)) v[i] <- TRUE
                else v[i] <- FALSE
            }
        v
    }

    if(length(variable)>1){
        return(sum(RE[any_match(RE$grp,variable), "vcov"])/sum(RE$vcov))
    }
}

fit <- lmer(growthday ~ factor(Replicate) + (1|Rhizobium) + (1|Clover),
            data=MEarea)

calculate_heritability(fit, c("Rhizobium", "Clover"))

fit <- lmer(Area ~ factor(Replicate) + (1|Rhizobium) + (1|Clover),
            data=MEarea)

calculate_heritability(fit, c("Rhizobium", "Clover"))

fit <- lmer(weight ~ factor(Replicate) + (1|Rhizobium) + (1|Clover),
            data=MEarea)

calculate_heritability(fit, c("Rhizobium", "Clover"))

fit <- lmer(AreaRatio ~ factor(Replicate) + (1|Rhizobium) + (1|Clover),
            data=MEarea)

calculate_heritability(fit, c("Rhizobium", "Clover"))
