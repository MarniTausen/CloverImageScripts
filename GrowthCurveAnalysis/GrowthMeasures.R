calculate_area <- function(x, y, a, b, baseline=0){
    square <- (abs(min(y, b)-baseline))*(abs(a-x))
    triangle <- ((abs(b-y))*(abs(a-x)))/2
    square+triangle
}

area_under_the_curve <- function(TotalArea, Time){
    base <- min(TotalArea)
    n <- length(TotalArea)
    total_area <- 0
    for(i in 1:(n-1)){
        total_area <- total_area + calculate_area(1, TotalArea[i], 2, TotalArea[i+1], base)
    }
    total_area
}

area_ratio <- function(TotalArea, Time){
    total_area <- area_under_the_curve(TotalArea, Time)
    max_area <- (max(TotalArea)-min(TotalArea))*length(TotalArea)
    total_area/max_area
}

select_time_frame <- function(Pot, TotalArea, Time, func){

    #print(as.character(Pot[1]))

### Figure out whether this is the replicate .1, or has one.
    sequence <- unlist(strsplit(as.character(Pot[1]), split=""))
    n <- length(sequence)
    if(all(sequence[(n-1):n]==c(".", "1"))) replicate <- TRUE
    else replicate <- FALSE

### Get name of replicate
    if(replicate) other <- paste(head(sequence, n=n-2), collapse = "")
    else other <- paste(c(sequence, ".", "1"), collapse = "")

    H1 <- max(Time)
    H2 <- max((TAf %>% filter(Pot==other, Period=="growth1"))$Time)
    H1_t <- HarvestData[as.character(Pot)[1], "inoculation_date"]
    H2_t <- HarvestData[other, "inoculation_date"]

    if(H1<=H2) { H <- H1; Ht <- H1_t }
    else { H <- H2; Ht <- H2_t }

    in1 <- HarvestData[as.character(Pot)[1], "inoculation_date"]
    in2 <- HarvestData[other, "inoculation_date"]

    start_date <- as.Date("2017-05-24", format="%Y-%m-%d")

    ### Find common regions, which are comparable.
    ## Find biggest common region, with min(HarvestDate)-max(Inoculationdate)

    if(in1<in2){
        y <- in2 - in1
        p_start_date <- start_date

        #cat(as.character(Pot[1]), Ht-H, y, "\n")

        if((Ht-H)>y){
            #cat("This is true \n")
            x <- H - p_start_date
        } else {
            x <- H - p_start_date - y
        }

        p_end_date <- p_start_date + x
    } else {
        y <- in1 - in2
        p_start_date <- start_date+y

        if((Ht-H)>y){
            #cat("This is also true \n")
            x <- H - start_date
            #x <- x + y
            #print(x)
        } else {
            #cat("Why is this one true? \n")
            x <- H - p_start_date
        }

        p_end_date <- p_start_date + x
    }

    data.frame(Time=Time, TotalArea) %>%
        filter(Time>=p_start_date, Time<=p_end_date) -> df

    Time <- df$Time
    TotalArea <- df$TotalArea

    func(TotalArea, Time)
}

is_or_has_replicate <- function(Name, names){
    sequence <- unlist(strsplit(as.character(Name), split=""))
    n <- length(sequence)
    if(all(sequence[(n-1):n]==c(".", "1"))) is_replicate <- TRUE
    else is_replicate <- FALSE
    if(is_replicate){
        non_replicate <- paste(sequence[1:(n-2)], collapse="")
        if(any(names==non_replicate)) return(TRUE)
        else return(FALSE)
    } else {
        non_replicate <- paste(c(sequence, ".1"), collapse="")
        if(any(names==non_replicate)) return(TRUE)
        else return(FALSE)
    }
}

all_names <- as.character(unique(TAf$Pot))
TAf %>% filter(is_or_has_replicate(Pot, all_names)) %>% group_by(Pot, Period) %>%
    summarise(AreaRatio=area_ratio(TotalArea, Time),
              Area=area_under_the_curve(TotalArea, Time),
              incoluation_date=HarvestData[Pot, "inoculation_date"][1],
              harvest_date=HarvestData[Pot, "harvest_date"][1]) %>%
    arrange(Pot) -> TAarea

TAf %>% filter(is_or_has_replicate(Pot, all_names), Period=="growth1") %>%
    group_by(Pot, Period) %>%
    summarise(AreaRatio=select_time_frame(Pot, TotalArea, Time, area_ratio),
              Area=select_time_frame(Pot, TotalArea, Time, area_under_the_curve)) -> TAarea

TAarea$Pot <- as.character(TAarea$Pot)

print(as.data.frame(TAarea %>% filter(Period=="growth1") %>% arrange(Pot)))

#write.table(as.data.frame(TAarea %>% arrange(Pot)),
#            file="AreaMeasurements2.csv", sep=",", row.names=FALSE)





select_time_frame_filter <- function(Pot, TotalArea, Time){

    #print(as.character(Pot[1]))

### Figure out whether this is the replicate .1, or has one.
    sequence <- unlist(strsplit(as.character(Pot[1]), split=""))
    n <- length(sequence)
    if(all(sequence[(n-1):n]==c(".", "1"))) replicate <- TRUE
    else replicate <- FALSE

### Get name of replicate
    if(replicate) other <- paste(head(sequence, n=n-2), collapse = "")
    else other <- paste(c(sequence, ".", "1"), collapse = "")

    H1 <- max(Time)
    H2 <- max((TAf %>% filter(Pot==other, Period=="growth1"))$Time)
    H1_t <- HarvestData[as.character(Pot)[1], "harvest_date"]
    H2_t <- HarvestData[other, "harvest_date"]

    if(H1<=H2) { H <- H1; Ht <- H1_t }
    else { H <- H2; Ht <- H2_t }

    in1 <- HarvestData[as.character(Pot)[1], "inoculation_date"]
    in2 <- HarvestData[other, "inoculation_date"]

    start_date <- as.Date("2017-05-24", format="%Y-%m-%d")

    ### Find common regions, which are comparable.
    ## Find biggest common region, with min(HarvestDate)-max(Inoculationdate)

    if(in1<in2){
        y <- in2 - in1
        p_start_date <- start_date

        #cat(as.character(Pot[1]), Ht-H, y, "\n")

        if((Ht-H)>y){
            cat("This is true \n")
            x <- H - p_start_date
        } else {
            x <- H - p_start_date - y
        }

        p_end_date <- p_start_date + x
    } else {
        y <- in1 - in2
        p_start_date <- start_date+y

        if((Ht-H)>y){
            cat("This is also true \n")
            x <- H - start_date
            #x <- x + y
            #print(x)
        } else {
            cat("Why is this one true? \n")
            x <- H - p_start_date
        }

        p_end_date <- p_start_date + x
    }

    v <- ((Time>=p_start_date)+(Time<=p_end_date))>1
    v
}

testGrowthCurveMini <- function(Plant, filter=TRUE){
    Plants <- c(Plant, paste(Plant, ".1", sep=""))
    print(Plants)
    if(filter){
        TAf %>% filter(any(as.character(Pot)==Plants), Period=="growth1") %>%
            group_by(Pot) %>% filter(select_time_frame_filter(Pot, TotalArea, Time)) -> TAtemp
    } else {
        TAf %>% filter(any(as.character(Pot)==Plants), Period=="growth1") -> TAtemp
    }

    cat("Incoluation date of:", Plants[1], "\n")
    print(HarvestData[Plants[1], "inoculation_date"])
    cat("Incoluation date of:", Plants[2], "\n")
    print(HarvestData[Plants[2], "inoculation_date"])
    cat("Harvest date of:", Plants[1], "\n")
    print(HarvestData[Plants[1], "harvest_date"])
    cat("Harvest date of:", Plants[2], "\n")
    print(HarvestData[Plants[2], "harvest_date"])

    print(TAarea[TAarea$Pot==Plants[1],])
    print(TAarea[TAarea$Pot==Plants[2],])

    TAplot <- TAtemp

    figure <- ggplot(TAplot, aes(x=Time, y=TotalArea, color=Pot)) + geom_line(size=1.05) +
        geom_errorbar(data=TAplot, aes(ymin=TotalArea-Std, ymax=TotalArea+Std),
                      alpha=0.6, width=0.1) + geom_point(aes(color=Predicted)) +
        labs(x="", title=paste(pickone(TAplot$Pot), pickone(TAplot$Camera),
                               pickone(TAplot$Placement), collapse=" ")) +
        scale_x_date(date_labels = "%b %d", date_breaks = "2 days") + theme_bw() +
        scale_y_continuous(limits=c(0, 1500)) +
        scale_color_manual(values=c("red", "blue", "gray", "black"))


    print(figure)

}

showtime <- function(pot){
    (TAf %>% filter(Pot==pot))$Time
}
