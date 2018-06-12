PotCam <- read.csv("PotCam.csv", row.names = 1)
#TestPotCam <- read.csv("test.csv", row.names = 1)

Pots <- colnames(PotCam)
Cams <- rownames(PotCam)

camera <- list()

for(Cam in Cams){
    for(Pot in Pots){
        Rname <- paste(unlist(strsplit(as.character(PotCam[Cam, Pot]), split="-")),
                       collapse=".")
        if(!exists(Rname, where=camera)){
            camera[[Rname]] <- paste(Pot, Cam, sep="_")
        } else {
            i <- 1
            while(exists(paste(Rname, as.character(i), sep="."), where=camera)){
                i <- i+1
            }
            camera[[paste(Rname, as.character(i), sep=".")]] <-
                paste(Pot, Cam, sep="_")
        }
    }
}

convertToPlacement <- Vectorize(function(x) unlist(camera[[as.character(x)]]))
convertToPot <- Vectorize(function(x) unlist(strsplit(
                                          as.character(
                                              unlist(convertToPlacement(x))),"_"))[1])
convertToCamera <- Vectorize(function(x) unlist(strsplit(
                                             as.character(
                                                 unlist(convertToPlacement(x))),"_"))[2])


revCamera <- list()

for(placement in names(camera)){
    revCamera[[camera[[placement]]]] <- placement
}

TotalPot <- function(x){
    TableNumber <- floor(x/10)
    PotNumber <- as.integer(((x/10)-CamNumber)*10+1)
    revCamera[[paste(c("Pot", as.character(PotNumber), "_Cam", as.character(CamNumber)), collapse="")]]
}
