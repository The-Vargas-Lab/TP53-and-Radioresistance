#Load necessary libraries
library("reshape2")
library("openxlsx")
library("readxl")
library("tidyverse")
library("BIGL")


#Analyzing Drug + radiation combination experiments (Figure 4, S4)


##Data prep for synergy analyses
###Make list of sheetnames from Excel file and loop through to create single dataframe for each sheet 
sheetnumber <- seq(1:9)
sheet_names <- c()
for (i in 1:length(sheetnumber)){
  sheet_names <- append(sheet_names, (paste0("Sheet",i)))
}
AUC_1 <- list()
AUC_2 <- list()
filename_1 <- "JHUEM2, Hec108, Hec1B -Nutlin + Rad, AP47, 5-14-21"
filename_2 <- "JHUEM2, Hec108, Hec1B -AMG232 + Rad, AP46, 5-7-21"
Drug_1 <- "Nutlin-3"
Drug_2 <- "AMG-232"
cell_lines <- c("JHUEM2", "Hec108", "Hec1B")

###Read in Nutlin + Radiation data file
for (i in 1:length(sheet_names)){
  AUC_1[[i]] <-read_xlsx(paste0(filename_1,".xlsx"), sheet_names[i])
}

###Read in AMG-232 + Radiation data file
for (i in 1:length(sheet_names)){
  AUC_2[[i]] <-read_xlsx(paste0(filename_2,".xlsx"), sheet_names[i])
}

###Extract list of radiation doses and reverse order for use later
Rad <- as.character(rev(unique(AUC_1[[1]]$Radiation)))

###Loop through each dataframe, take the mean of 6 replicates, and restructure, making cell names into rownames 
for (i in 1:length(AUC_1)){
  AUC_1[[i]] <- melt(AUC_1[[i]], id = c("Dose", "Radiation"))
  AUC_1[[i]] <- dcast(AUC_1[[i]], Dose ~ Radiation, mean)
  row.names(AUC_1[[i]]) <- AUC_1[[i]]$Dose
  AUC_1[[i]] <- AUC_1[[i]][,-1]
  AUC_2[[i]] <- melt(AUC_2[[i]], id = c("Dose", "Radiation"))
  AUC_2[[i]] <- dcast(AUC_2[[i]], Dose ~ Radiation, mean)
  row.names(AUC_2[[i]]) <- AUC_2[[i]]$Dose
  AUC_2[[i]] <- AUC_2[[i]][,-1]
  
  
  ###Loop through data tables and normalize each response value to Untreated sample
  
  AUC_1[[i]] <- AUC_1[[i]]/AUC_1[[i]]["0", "0"]
  AUC_2[[i]] <- AUC_2[[i]]/AUC_2[[i]]["0", "0"]
  
  ###For Radiation doses of 2, 4, 8 Gy, divide fraction surviving by 2 because plated at twice the density
  
  for (j in Rad){
    if (j == "2" | j == "4" | j == "8"){
      AUC_1[[i]][j] <- (AUC_1[[i]][j])/2
      AUC_2[[i]][j] <- (AUC_2[[i]][j])/2
    }
    else {
      AUC_1[[i]][j] <- AUC_1[[i]][j]
      AUC_2[[i]][j] <- AUC_2[[i]][j]
    }
  }
  
  
  ###Transpose each table and extract to separate dataframe
  AUC_1[[i]] <- as.data.frame(t(AUC_1[[i]]))
  AUC_2[[i]] <- as.data.frame(t(AUC_2[[i]]))
  
  ###Create new Radiation column for Excel labels and pivot to long form
  AUC_1[[i]]$Radiation <- rownames(AUC_1[[i]])
  AUC_1[[i]] <- melt(AUC_1[[i]], id= c("Radiation"))
  AUC_2[[i]]$Radiation <- rownames(AUC_2[[i]])
  AUC_2[[i]] <- melt(AUC_2[[i]], id= c("Radiation"))
}


###For each cell line, combine 3 replicate tables, take the mean, normalize all to the maximum viability, and rename columns
k <- 1
while(k <= length(AUC_1)){
  AUC_1[[k]] <- rbind(AUC_1[[k]], AUC_1[[k+1]], AUC_1[[k+2]])
  AUC_1[[k]] <- dcast(AUC_1[[k]], Radiation ~ variable, mean)
  AUC_1[[k]] <- melt(AUC_1[[k]], id= c("Radiation"))
  AUC_1[[k]]$value <- (AUC_1[[k]]$value)/(max(AUC_1[[k]]$value))
  names(AUC_1[[k]])[1:3] <- c("d2", "d1", "effect")
  AUC_2[[k]] <- rbind(AUC_2[[k]], AUC_2[[k+1]], AUC_2[[k+2]])
  AUC_2[[k]] <- dcast(AUC_2[[k]], Radiation ~ variable, mean)
  AUC_2[[k]] <- melt(AUC_2[[k]], id= c("Radiation"))
  AUC_2[[k]]$value <- (AUC_2[[k]]$value)/(max(AUC_2[[k]]$value))
  names(AUC_2[[k]])[1:3] <- c("d2", "d1", "effect")
  k <- k + 3
}

###Add column with cell line name
AUC_1[[1]]$Cell <- cell_lines[1]
AUC_1[[4]]$Cell <- cell_lines[2]
AUC_1[[7]]$Cell <- cell_lines[3]
AUC_2[[1]]$Cell <- cell_lines[1]
AUC_2[[4]]$Cell <- cell_lines[2]
AUC_2[[7]]$Cell <- cell_lines[3]

###Convert to dataframes
Mean_1 <- as.data.frame(AUC_1[[1]])
Mean_2 <- as.data.frame(AUC_1[[4]])
Mean_3 <- as.data.frame(AUC_1[[7]])
Mean_4 <- as.data.frame(AUC_2[[1]])
Mean_5 <- as.data.frame(AUC_2[[4]])
Mean_6 <- as.data.frame(AUC_2[[7]])

###Convert Drug (d1) and Radiation (d2) doses to numeric 
####d1 is factor, so first needs changed to character
Mean_1$d1 <- as.character(Mean_1$d1)
Mean_2$d1 <- as.character(Mean_2$d1)
Mean_3$d1 <- as.character(Mean_3$d1)
Mean_4$d1 <- as.character(Mean_4$d1)
Mean_5$d1 <- as.character(Mean_5$d1)
Mean_6$d1 <- as.character(Mean_6$d1)
Mean_1$d1 <- as.numeric(Mean_1$d1)
Mean_2$d1 <- as.numeric(Mean_2$d1)
Mean_3$d1 <- as.numeric(Mean_3$d1)
Mean_4$d1 <- as.numeric(Mean_4$d1)
Mean_5$d1 <- as.numeric(Mean_5$d1)
Mean_6$d1 <- as.numeric(Mean_6$d1)
Mean_1$d2 <- as.numeric(Mean_1$d2)
Mean_2$d2 <- as.numeric(Mean_2$d2)
Mean_3$d2 <- as.numeric(Mean_3$d2)
Mean_4$d2 <- as.numeric(Mean_4$d2)
Mean_5$d2 <- as.numeric(Mean_5$d2)
Mean_6$d2 <- as.numeric(Mean_6$d2)




##Synergy analysis with BIGL
###Fit marginal dose response curves to data using fitMarginal function in BIGL
Marg1Fit <- fitMarginals(Mean_1, method= "nlslm", names= c(Drug_1, "Radiation"))
Marg2Fit <- fitMarginals(Mean_2, method= "nlslm", names= c(Drug_1, "Radiation"))
Marg3Fit <- fitMarginals(Mean_3, method= "nlslm", names= c(Drug_1, "Radiation"))
Marg4Fit <- fitMarginals(Mean_4, method= "nlslm", names= c(Drug_2, "Radiation"))
Marg5Fit <- fitMarginals(Mean_5, method= "nlslm", names= c(Drug_2, "Radiation"))
Marg6Fit <- fitMarginals(Mean_6, method= "nlslm", names= c(Drug_2, "Radiation"))

###Use fitSurface function in BIGL to fit predicted surface using Highest Single Agent (HSA)(JHUEM2, Hec108) and Loewe models (Hec1B)
Cell1Surf <- fitSurface(Mean_1, Marg1Fit, null_model = "hsa", B.CP = 50, statistic = "both")
Cell2Surf <- fitSurface(Mean_2, Marg2Fit, null_model = "hsa", B.CP = 50, statistic = "both")
Cell3Surf <- fitSurface(Mean_3, Marg3Fit, null_model = "loewe", B.CP = 50, statistic = "both")
Cell4Surf <- fitSurface(Mean_4, Marg4Fit, null_model = "hsa", B.CP = 50, statistic = "both")
Cell5Surf <- fitSurface(Mean_5, Marg5Fit, null_model = "hsa", B.CP = 50, statistic = "both")
Cell6Surf <- fitSurface(Mean_6, Marg6Fit, null_model = "loewe", B.CP = 50, statistic = "both")

###Generate isobolograms of predicted surfaces using isobologram function in BIGL and save as pdfs
pdf(paste0(cell_lines[1], " + ", Drug_1, ", HSA.pdf"), width= 8.00, height= 5.00)
isobologram(Cell1Surf)
dev.off()

pdf(paste0(cell_lines[2], " + ", Drug_1, ", HSA.pdf"), width= 8.00, height= 5.00)
isobologram(Cell2Surf)
dev.off()

pdf(paste0(cell_lines[3], " + ", Drug_1, ", HSA.pdf"), width= 8.00, height= 5.00)
isobologram(Cell3Surf)
dev.off()

pdf(paste0(cell_lines[1], " + ", Drug_2, ", HSA.pdf"), width= 8.00, height= 5.00)
isobologram(Cell4Surf)
dev.off()

pdf(paste0(cell_lines[2], " + ", Drug_2, ", HSA.pdf"), width= 8.00, height= 5.00)
isobologram(Cell5Surf)
dev.off()

pdf(paste0(cell_lines[3], " + ", Drug_2, ", HSA.pdf"), width= 8.00, height= 5.00)
isobologram(Cell6Surf)
dev.off()


###Generate contour plots of synergy scores using contour function in BIGL and save as pdfs
pdf(paste0(cell_lines[1], " + ", Drug_1, " contour, HSA.pdf"), width= 8.00, height= 5.00)
contour(Cell1Surf, main= paste0("Contour Plot for ", cell_lines[1], " MaxR"))
dev.off()

pdf(paste0(cell_lines[2], " + ", Drug_1, " contour, HSA.pdf"), width= 8.00, height= 5.00)
contour(Cell2Surf, main= paste0("Contour Plot for ", cell_lines[2], " MaxR"))
dev.off()

pdf(paste0(cell_lines[3], " + ", Drug_1, " contour, HSA.pdf"), width= 8.00, height= 5.00)
contour(Cell3Surf, main= paste0("Contour Plot for ", cell_lines[3], " MaxR"))
dev.off()

pdf(paste0(cell_lines[1], " + ", Drug_2, " contour, HSA.pdf"), width= 8.00, height= 5.00)
contour(Cell4Surf, main= paste0("Contour Plot for ", cell_lines[1], " MaxR"))
dev.off()

pdf(paste0(cell_lines[2], " + ", Drug_2, " contour, HSA.pdf"), width= 8.00, height= 5.00)
contour(Cell5Surf, main= paste0("Contour Plot for ", cell_lines[2], " MaxR"))
dev.off()

pdf(paste0(cell_lines[3], " + ", Drug_2, " contour, HSA.pdf"), width= 8.00, height= 5.00)
contour(Cell6Surf, main= paste0("Contour Plot for ", cell_lines[3], " MaxR"))
dev.off()






