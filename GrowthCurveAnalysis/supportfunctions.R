contains <- function(x, m){
    if(m==x) return(TRUE)
    FALSE
}
contains <- Vectorize(contains)

mergeVectors <- function(x, y){
    if(length(x)!=length(y)) stop("Vectors are of unequal length")
    z <- vector(length=length(x))
    for(i in 1:length(x)){
        tmp <- paste(x[i], y[i], sep=".")
        j <- 0
        while(any(z==tmp)){
            j <- j+1
            tmp <- paste(x[i], y[i], as.character(j), sep=".")
        }
        z[i] <- tmp
    }
    z
}

divide_by_harvest <- function(df, combination_table){
    v <- vector(length=nrow(df))
    for(i in 1:nrow(df)){
        if(df$Time[i]<combination_table[df$Pot[i], 'Harvest']){
            v[i] <- "growth1"
        } else {
            v[i] <- "growth2"
        }
    }
    v
}

pickone <- function(x) x[1]


predict_points <- function(df){

    new_df = data.frame()
    last_day = NULL
    new_pot = TRUE
    pot_name = NULL
    period = NULL
    phantom = vector()

    for(row in 1:nrow(df)){

        if(new_pot==TRUE){
            last_day = df$Time[row]
            pot_name = as.character(df$Pot[row])
            period = df$Period[row]
            new_pot = FALSE
            if(nrow(new_df)==0) new_df = df[row,]
            else new_df = rbind(new_df, df[row,])
            phantom <- c(phantom, "true data")
        } else {

            if(as.numeric(df$Time[row]-last_day) > 1){
                Days <- last_day+1:(as.numeric(df$Time[row]-last_day)-1)

                ## Find the results inbetween the Days.
                temp <- data.frame(Time=df$Time[(row-1):row],
                                   TotalArea=df$TotalArea[(row-1):row])
                ### Do linear regression inbetween points
                fit <- lm(TotalArea ~ Time, data=temp)
                ### Predict points between.
                new_points <- predict(fit, newdata=data.frame(Time=Days))
                new_days <- data.frame(Pot=rep(df$Pot[row], times=length(Days)),
                                       Time=Days,
                                       Std=rep(0, times=length(Days)),
                                       TotalArea=new_points,
                                       Camera=rep(df$Camera[row], times=length(Days)),
                                       Placement=rep(df$Placement[row], times=length(Days)),
                                       Period=rep(df$Period[row], times=length(Days)))

                new_days %>% group_by(Pot, Time) -> new_days
                new_df = rbind(new_df, new_days)
                phantom <- c(phantom, rep("predicted data", times=length(Days)))

            }

            last_day = df$Time[row]
            new_df = rbind(new_df, df[row,])
            phantom <- c(phantom, "true data")

        }

        if((row+1)<nrow(df) && pot_name!=as.character(df$Pot[row+1])){
            new_pot = TRUE
        }
        if((row+1)<nrow(df) && period!=df$Period[row+1]){
            new_pot = TRUE
        }
    }

    new_df$Predicted = phantom
    new_df
}
