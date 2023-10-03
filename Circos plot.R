###Circos plot


library(tidyverse)
library(data.table)


data <- read.csv("circos1.csv", header = T)
datax <- data.frame(data$Pathway.1,data$G.annotion) 
datax$y <- 1

colnames(datax) <- c("ID", "Loc", "y")
datax1 <- dcast(data=datax, ID ~ Loc)


for (i in 1:ncol(datax1)){
  datax1[,i][is.na(datax1[,i])] <- 0
}



biga <- read.csv("biga3.csv", header = T)

datax3 <- merge(biga, datax1, by.x = "Pathway.1", by.y = "ID", all = F)
datax4 <- datax3[order(datax3$pathwayID),]
row.names(datax4) <- datax4$pathwayID


library(circlize)
library(RColorBrewer)

grid.col <- c(
  G1_Cellular.Processes = "#FFCC33",
  G2_Cellular.Processes = "#FFCC33",
  G3_Cellular.Processes = "#FFCC33",
  G5_Cellular.Processes = "#FFCC33",
  G4_Cellular.Processes = "#FFCC33",
  G1_Environmental.Information.Processing = "#99FFFF",
  G2_Environmental.Information.Processing = "#99FFFF",
  G4_Environmental.Information.Processing = "#99FFFF",
  G1_Genetic.Information.Processing = "#FF99FF",
  G2_Genetic.Information.Processing = "#FF99FF",
  G5_Genetic.Information.Processing = "#FF99FF",
  G4_Genetic.Information.Processing = "#FF99FF",
  G6_Genetic.Information.Processing = "#FF99FF",
  G1_Human.Diseases = "#33CCFF",
  G2_Human.Diseases = "#33CCFF",
  G5_Human.Diseases = "#33CCFF",
  G4_Human.Diseases = "#33CCFF",
  G6_Human.Diseases = "#33CCFF",
  G1_Metabolism = "#FF3333",
  G2_Metabolism = "#FF3333",
  G3_Metabolism = "#FF3333",
  G5_Metabolism = "#FF3333",
  G4_Metabolism = "#FF3333",
  G6_Metabolism = "#FF3333",
  G1_Organismal.Systems = "#00FF66",
  G2_Organismal.Systems = "#00FF66",
  G3_Organismal.Systems = "#00FF66",
  G5_Organismal.Systems = "#00FF66",
  G4_Organismal.Systems = "#00FF66",
  P01 = "#999999",
  P02 = "#999999",
  P03 = "#999999",
  P04 = "#999999",
  P05 = "#999999",
  P06 = "#999999",
  P07 = "#999999",
  P08 = "#999999",
  P09 = "#999999",
  P10 = "#999999",
  P11 = "#999999",
  P12 = "#999999",
  P13 = "#999999",
  P14 = "#999999",
  P15 = "#999999",
  P16 = "#999999",
  P17 = "#999999",
  P18 = "#999999",
  P19 = "#999999",
  P20 = "#999999",
  P21 = "#999999",
  P22 = "#999999",
  P23 = "#999999",
  P24 = "#999999",
  P25 = "#999999",
  P26 = "#999999",
  P27 = "#999999",
  P28 = "#999999",
  P29 = "#999999",
  P30 = "#999999",
  P31 = "#999999",
  P32 = "#999999",
  P33 = "#999999",
  P34 = "#999999",
  P35 = "#999999",
  P36 = "#999999",
  P37 = "#999999",
  P38 = "#999999",
  P39 = "#999999",
  P40 = "#999999",
  P41 = "#999999",
  P42 = "#999999",
  P43 = "#999999",
  P44 = "#999999",
  P45 = "#999999",
  P46 = "#999999",
  P47 = "#999999",
  P48 = "#999999",
  P49 = "#999999",
  P50 = "#999999",
  P51 = "#999999",
  P52 = "#999999",
  P53 = "#999999",
  P54 = "#999999",
  P55 = "#999999",
  P56 = "#999999",
  P57 = "#999999",
  P58 = "#999999",
  P59 = "#999999",
  P60 = "#999999",
  P61 = "#999999",
  P62 = "#999999",
  P63 = "#999999",
  P64 = "#999999",
  P65 = "#999999",
  P66 = "#999999",
  P67 = "#999999",
  P68 = "#999999",
  P69 = "#999999",
  P70 = "#999999",
  P71 = "#999999",
  P72 = "#999999",
  P73 = "#999999",
  P74 = "#999999",
  P75 = "#999999",
  P76 = "#999999",
  P77 = "#999999",
  P78 = "#999999",
  P79 = "#999999",
  P80 = "#999999",
  P81 = "#999999",
  P82 = "#999999",
  P83 = "#999999",
  P84 = "#999999",
  P85 = "#999999",
  P86 = "#999999")



set.seed(123)



chordDiagram(t(datax4[,4:32]), 
             grid.col = grid.col,
             annotationTrack = c("grid","axis"), 
             #directional = 1, 
             transparency = 0.5,
             preAllocateTracks = list( track.height = uh(3,"mm"), 
             track.margin = c(uh(6, "mm"), 0) ),
             annotationTrackHeight = mm_h(c(2, 4)),
             scale = FALSE,
             big.gap = 5,
             link.border = NA)




circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  #circos.text(mean(xlim), mean(ylim), sector.index, col = "white", cex = 0.6, niceFacing = TRUE)
})




Addiciton.vs.Control <- c('G1_Metabolism',
                          'G1_Genetic.Information.Processing',
                          'G1_Environmental.Information.Processing',
                          'G1_Cellular.Processes', 
                          'G1_Organismal.Systems',
                          'G1_Human.Diseases'
                          )

Withdrawal.vs.Control <- c('G2_Metabolism',
                           'G2_Genetic.Information.Processing',
                           'G2_Environmental.Information.Processing',
                           'G2_Cellular.Processes', 
                           'G2_Organismal.Systems',
                           'G2_Human.Diseases'
                           )

Methadone.vs.Control <- c('G3_Metabolism',
                          'G3_Cellular.Processes',
                          'G3_Organismal.Systems'
                          )

Addiciton.vs.Withdrawal <- c('G4_Metabolism',
                             'G4_Genetic.Information.Processing',
                             'G4_Environmental.Information.Processing',
                             'G4_Cellular.Processes', 
                             'G4_Organismal.Systems',
                             'G4_Human.Diseases'
                             )

Addiciton.vs.Methadone <- c('G5_Metabolism',
                            'G5_Genetic.Information.Processing',
                            'G5_Cellular.Processes', 
                            'G5_Organismal.Systems',
                            'G5_Human.Diseases'
                            )

Withdrawal.vs.Methadone <- c('G6_Metabolism',
                             'G6_Genetic.Information.Processing',
                             'G6_Human.Diseases'
                             )


highlight.sector(Addiciton.vs.Control, track.index = 1, col = "#0099FF", text = "Addiciton.vs.Control", 
                 cex = 0.8, text.col = "white", niceFacing = TRUE)


highlight.sector(Withdrawal.vs.Control, track.index = 1, col = "#CC1200", text = "Withdrawal.vs.Control", 
                 cex = 0.8, text.col = "white", niceFacing = TRUE)

highlight.sector(Methadone.vs.Control, track.index = 1, col = "#003399", text = "Methadone.vs.Control", 
                 cex = 0.8, text.col = "white", niceFacing = TRUE)

highlight.sector(Addiciton.vs.Withdrawal, track.index = 1, col = "#FF6600", text = "Addiciton.vs.Withdrawal", 
                 cex = 0.8, text.col = "white", niceFacing = TRUE)

highlight.sector(Addiciton.vs.Methadone, track.index = 1, col = "#CC00FF", text = "Addiciton.vs.Methadone", 
                 cex = 0.8, text.col = "white", niceFacing = TRUE)

highlight.sector(Withdrawal.vs.Methadone, track.index = 1, col = "#339999", text = "Withdrawal.vs.Methadone", 
                 cex = 0.8, text.col = "white", niceFacing = TRUE)

Human.Diseases <- biga$pathwayID[74:86]
Organismal.Systems <- biga$pathwayID[65:73]
Cellular.Processes <- biga$pathwayID[56:64]
Environmental.Information.Processing <- biga$pathwayID[50:55]
Genetic.Information.Processing <- biga$pathwayID[42:49]
Metabolism <- biga$pathwayID[1:41]


highlight.sector(Cellular.Processes, track.index = 1, col = "#FFCC33", text = "Cellular.Processes", 
                 cex = 0.8, text.col = "white", niceFacing = TRUE)

highlight.sector(Environmental.Information.Processing, track.index = 1, col = "#99FFFF", text = "Environmental.Information.Processing", 
                 cex = 0.8, text.col = "white", niceFacing = TRUE)

highlight.sector(Genetic.Information.Processing, track.index = 1, col = "#FF99FF", text = "Genetic.Information.Processing", 
                 cex = 0.8, text.col = "white", niceFacing = TRUE)

highlight.sector(Human.Diseases, track.index = 1, col = "#33CCFF", text = "Human.Diseases", 
                 cex = 0.8, text.col = "white", niceFacing = TRUE)

highlight.sector(Metabolism, track.index = 1, col = "#FF3333", text = "Metabolism", 
                 cex = 0.8, text.col = "white", niceFacing = TRUE)

highlight.sector(Organismal.Systems, track.index = 1, col = "#00FF66", text = "Organismal.Systems", 
                 cex = 0.8, text.col = "white", niceFacing = TRUE)



