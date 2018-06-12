#TA %>% group_by(Pot, Time) %>%
#    summarise(Std=sd(TotalArea), TotalArea=median(TotalArea)) -> TAf

source("manual_list.R")

Cams <- Black_list

for(Cam in Cams){
    #TAf %>% filter(Camera==Cam) -> TAplot
    TAblacklist %>% filter(Pot==Cam) -> TAplot

    print(Cam)

    if(nrow(TAplot)==0) next

    TAplot %>% filter(Period=="growth1") -> TAplot1
    TAplot %>% filter(Period=="growth2") -> TAplot2

    figure <- ggplot(TAplot, aes(x=Time, y=TotalArea, color=Period)) + geom_line(size=1.05) +
        geom_errorbar(data=TAplot, aes(ymin=TotalArea-Std, ymax=TotalArea+Std),
                      alpha=0.6, width=0.1) + geom_point(color="black") +
        labs(x="", title=paste(pickone(TAplot$Pot), pickone(TAplot$Camera),
                               pickone(TAplot$Placement), collapse=" "),
             y="Total Coverage (cm^2)") +
        scale_x_date(date_labels = "%b %d", date_breaks = "2 days") + theme_bw() +
        scale_y_continuous(limits=c(0, 1600))

    if(nrow(TAplot1)!=0){
        logfit1 <- drm(TotalArea ~ Time, data=TAplot1, fct=L.4())
        figure <- figure + geom_line(data=TAplot1, aes(x=Time, y=predict(logfit1)))
    }
    if(nrow(TAplot2)!=0){
        logfit2 <- drm(TotalArea ~ Time, data=TAplot2, fct=L.4())
        figure <- figure + geom_line(data=TAplot2, aes(x=Time, y=predict(logfit2)))
    }

    print(figure)

}
