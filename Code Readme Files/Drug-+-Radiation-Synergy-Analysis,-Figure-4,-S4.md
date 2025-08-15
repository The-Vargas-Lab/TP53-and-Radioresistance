Drug + Radiation Synergy Analysis, Figure 4, S4
================
Aaron Petty
2025-08-12

### *Data prep for synergy analyses*

#### Load packages

``` r
library("reshape2")
library("readxl")
library("tidyverse")
library("BIGL")
```

#### Make list of sheetnames from Excel file and loop through to create single dataframe for each sheet

``` r
sheetnumber <- seq(1:9)
sheet_names <- c()
for (i in 1:length(sheetnumber)){
  sheet_names <- append(sheet_names, (paste0("Sheet",i)))
}
Nutlin_1 <- list()
AMG_2 <- list()
filename_1 <- "JHUEM2, Hec108, Hec1B -Nutlin + Rad, AP47, 5-14-21"
filename_2 <- "JHUEM2, Hec108, Hec1B -AMG232 + Rad, AP46, 5-7-21"
Drug_1 <- "Nutlin-3"
Drug_2 <- "AMG-232"
cell_lines <- c("JHUEM2", "Hec108", "Hec1B")
```

#### Read in Nutlin + Radiation data file

``` r
for (i in 1:length(sheet_names)){
  Nutlin_1[[i]] <-read_xlsx(paste0(filename_1,".xlsx"), sheet_names[i])
}
```

#### Read in AMG-232 + Radiation data file

``` r
for (i in 1:length(sheet_names)){
  AMG_2[[i]] <-read_xlsx(paste0(filename_2,".xlsx"), sheet_names[i])
}
```

#### Extract list of radiation doses and reverse order for use later

``` r
Rad <- as.character(rev(unique(Nutlin_1[[1]]$Radiation)))
```

#### Loop through each dataframe to normalize and restructure data

``` r
for (i in 1:length(Nutlin_1)){
  # Take the mean of 6 replicates, and restructure, making drug doses into     rownames 
  Nutlin_1[[i]] <- melt(Nutlin_1[[i]], id = c("Dose", "Radiation"))
  Nutlin_1[[i]] <- dcast(Nutlin_1[[i]], Dose ~ Radiation, mean)
  row.names(Nutlin_1[[i]]) <- Nutlin_1[[i]]$Dose
  Nutlin_1[[i]] <- Nutlin_1[[i]][,-1]
  AMG_2[[i]] <- melt(AMG_2[[i]], id = c("Dose", "Radiation"))
  AMG_2[[i]] <- dcast(AMG_2[[i]], Dose ~ Radiation, mean)
  row.names(AMG_2[[i]]) <- AMG_2[[i]]$Dose
  AMG_2[[i]] <- AMG_2[[i]][,-1]

  # Normalize each response value to Untreated sample
  
  Nutlin_1[[i]] <- Nutlin_1[[i]]/Nutlin_1[[i]]["0", "0"]
  AMG_2[[i]] <- AMG_2[[i]]/AMG_2[[i]]["0", "0"]
  
  # In JHUEM2 and Hec108 at Radiation doses of 2, 4, 8 Gy, divide fraction surviving by 2 because plated at twice the density
  
    for (j in Rad){
    
      if(i==7 | i==8 | i==9){
      Nutlin_1[[i]][j] <- Nutlin_1[[i]][j]
      AMG_2[[i]][j] <- AMG_2[[i]][j]
    }
      
    else if (j=="0" | j=="0.5" | j=="1"){
      Nutlin_1[[i]][j] <- Nutlin_1[[i]][j]
     AMG_2[[i]][j] <- AMG_2[[i]][j]
    }
    
    else {
      Nutlin_1[[i]][j] <- (Nutlin_1[[i]][j])/2
      AMG_2[[i]][j] <- (AMG_2[[i]][j])/2
    }
  }
  
  # Transpose each table and extract to separate dataframe
  Nutlin_1[[i]] <- as.data.frame(t(Nutlin_1[[i]]))
  AMG_2[[i]] <- as.data.frame(t(AMG_2[[i]]))
  
  # Create new Radiation column for Excel labels and pivot to long form
  Nutlin_1[[i]]$Radiation <- rownames(Nutlin_1[[i]])
  Nutlin_1[[i]] <- melt(Nutlin_1[[i]], id= c("Radiation"))
  AMG_2[[i]]$Radiation <- rownames(AMG_2[[i]])
  AMG_2[[i]] <- melt(AMG_2[[i]], id= c("Radiation"))
}
```

#### For each cell line, combine 3 replicate tables, take the mean, normalize all to the maximum viability, and rename columns

``` r
k <- 1
while(k <= length(Nutlin_1)){
  Nutlin_1[[k]] <- rbind(Nutlin_1[[k]], Nutlin_1[[k+1]], Nutlin_1[[k+2]])
  Nutlin_1[[k]] <- dcast(Nutlin_1[[k]], Radiation ~ variable, mean)
  Nutlin_1[[k]] <- melt(Nutlin_1[[k]], id= c("Radiation"))
  Nutlin_1[[k]]$value <- (Nutlin_1[[k]]$value)/(max(Nutlin_1[[k]]$value))
  names(Nutlin_1[[k]])[1:3] <- c("d2", "d1", "effect")
  AMG_2[[k]] <- rbind(AMG_2[[k]], AMG_2[[k+1]], AMG_2[[k+2]])
  AMG_2[[k]] <- dcast(AMG_2[[k]], Radiation ~ variable, mean)
  AMG_2[[k]] <- melt(AMG_2[[k]], id= c("Radiation"))
  AMG_2[[k]]$value <- (AMG_2[[k]]$value)/(max(AMG_2[[k]]$value))
  names(AMG_2[[k]])[1:3] <- c("d2", "d1", "effect")
  k <- k + 3
}
```

#### Add column with cell line name

``` r
Nutlin_1[[1]]$Cell <- cell_lines[1]
Nutlin_1[[4]]$Cell <- cell_lines[2]
Nutlin_1[[7]]$Cell <- cell_lines[3]
AMG_2[[1]]$Cell <- cell_lines[1]
AMG_2[[4]]$Cell <- cell_lines[2]
AMG_2[[7]]$Cell <- cell_lines[3]
```

#### Convert to dataframes

``` r
Nutlin_JH <- as.data.frame(Nutlin_1[[1]])
Nutlin_108 <- as.data.frame(Nutlin_1[[4]])
Nutlin_1B <- as.data.frame(Nutlin_1[[7]])
AMG_JH <- as.data.frame(AMG_2[[1]])
AMG_108 <- as.data.frame(AMG_2[[4]])
AMG_1B <- as.data.frame(AMG_2[[7]])
```

#### Convert Drug (d1) and Radiation (d2) doses to numeric

``` r
# d1 is factor, so first needs changed to character
Nutlin_JH$d1 <- as.character(Nutlin_JH$d1)
Nutlin_108$d1 <- as.character(Nutlin_108$d1)
Nutlin_1B$d1 <- as.character(Nutlin_1B$d1)
AMG_JH$d1 <- as.character(AMG_JH$d1)
AMG_108$d1 <- as.character(AMG_108$d1)
AMG_1B$d1 <- as.character(AMG_1B$d1)
Nutlin_JH$d1 <- as.numeric(Nutlin_JH$d1)
Nutlin_108$d1 <- as.numeric(Nutlin_108$d1)
Nutlin_1B$d1 <- as.numeric(Nutlin_1B$d1)
AMG_JH$d1 <- as.numeric(AMG_JH$d1)
AMG_108$d1 <- as.numeric(AMG_108$d1)
AMG_1B$d1 <- as.numeric(AMG_1B$d1)
Nutlin_JH$d2 <- as.numeric(Nutlin_JH$d2)
Nutlin_108$d2 <- as.numeric(Nutlin_108$d2)
Nutlin_1B$d2 <- as.numeric(Nutlin_1B$d2)
AMG_JH$d2 <- as.numeric(AMG_JH$d2)
AMG_108$d2 <- as.numeric(AMG_108$d2)
AMG_1B$d2 <- as.numeric(AMG_1B$d2)
```

### *Synergy analysis with BIGL*

#### Fit individual dose response curves to data using fitMarginal function in BIGL

``` r
Nutlin_JH_Marg <- fitMarginals(Nutlin_JH, method= "nlslm", names= c(Drug_1, "Radiation"))
Nutlin_108_Marg <- fitMarginals(Nutlin_108, method= "nlslm", names= c(Drug_1, "Radiation"))
Nutlin_1B_Marg <- fitMarginals(Nutlin_1B, method= "nlslm", names= c(Drug_1, "Radiation"))
AMG_JH_Marg <- fitMarginals(AMG_JH, method= "nlslm", names= c(Drug_2, "Radiation"))
AMG_108_Marg <- fitMarginals(AMG_108, method= "nlslm", names= c(Drug_2, "Radiation"))
AMG_1B_Marg <- fitMarginals(AMG_1B, method= "nlslm", names= c(Drug_2, "Radiation"))
```

#### Use fitSurface function in BIGL to fit predicted surface using Highest Single Agent (HSA)(JHUEM2, Hec108) and Loewe models (Hec1B)

``` r
Nutlin_JH_Surf <- fitSurface(Nutlin_JH, Nutlin_JH_Marg, null_model = "hsa", B.CP = 50, statistic = "both")

Nutlin_108_Surf <- fitSurface(Nutlin_108, Nutlin_108_Marg, null_model = "hsa", B.CP = 50, statistic = "both")

Nutlin_1B_Surf <- fitSurface(Nutlin_1B, Nutlin_1B_Marg, null_model = "loewe", B.CP = 50, statistic = "both")

AMG_JH_Surf <- fitSurface(AMG_JH, AMG_JH_Marg, null_model = "hsa", B.CP = 50, statistic = "both")

AMG_108_Surf <- fitSurface(AMG_108, AMG_108_Marg, null_model = "hsa", B.CP = 50, statistic = "both")

AMG_1B_Surf <- fitSurface(AMG_1B, AMG_1B_Marg, null_model = "loewe", B.CP = 50, statistic = "both")
```

#### Generate isobolograms of predicted surfaces using isobologram function in BIGL and save as pdfs

``` r
isobologram(Nutlin_JH_Surf)
```

![](Drug-+-Radiation-Synergy-Analysis,-Figure-4,-S4_files/figure-gfm/isobol-1.pdf)<!-- -->

``` r
isobologram(Nutlin_108_Surf)
```

![](Drug-+-Radiation-Synergy-Analysis,-Figure-4,-S4_files/figure-gfm/isobol-2.pdf)<!-- -->

``` r
isobologram(Nutlin_1B_Surf)
```

![](Drug-+-Radiation-Synergy-Analysis,-Figure-4,-S4_files/figure-gfm/isobol-3.pdf)<!-- -->

``` r
isobologram(AMG_JH_Surf)
```

![](Drug-+-Radiation-Synergy-Analysis,-Figure-4,-S4_files/figure-gfm/isobol-4.pdf)<!-- -->

``` r
isobologram(AMG_108_Surf)
```

![](Drug-+-Radiation-Synergy-Analysis,-Figure-4,-S4_files/figure-gfm/isobol-5.pdf)<!-- -->

``` r
isobologram(AMG_1B_Surf)
```

![](Drug-+-Radiation-Synergy-Analysis,-Figure-4,-S4_files/figure-gfm/isobol-6.pdf)<!-- -->

#### Generate contour plots of synergy scores using contour function in BIGL and save as pdfs

``` r
contour(Nutlin_JH_Surf, main= paste0("Contour Plot for ", cell_lines[1], " MaxR"))
```

![](Drug-+-Radiation-Synergy-Analysis,-Figure-4,-S4_files/figure-gfm/contour-1.pdf)<!-- -->

``` r
contour(Nutlin_108_Surf, main= paste0("Contour Plot for ", cell_lines[2], " MaxR"))
```

![](Drug-+-Radiation-Synergy-Analysis,-Figure-4,-S4_files/figure-gfm/contour-2.pdf)<!-- -->

``` r
contour(Nutlin_1B_Surf, main= paste0("Contour Plot for ", cell_lines[3], " MaxR"))
```

![](Drug-+-Radiation-Synergy-Analysis,-Figure-4,-S4_files/figure-gfm/contour-3.pdf)<!-- -->

``` r
contour(AMG_JH_Surf, main= paste0("Contour Plot for ", cell_lines[1], " MaxR"))
```

![](Drug-+-Radiation-Synergy-Analysis,-Figure-4,-S4_files/figure-gfm/contour-4.pdf)<!-- -->

``` r
contour(AMG_108_Surf, main= paste0("Contour Plot for ", cell_lines[2], " MaxR"))
```

![](Drug-+-Radiation-Synergy-Analysis,-Figure-4,-S4_files/figure-gfm/contour-5.pdf)<!-- -->

``` r
contour(AMG_1B_Surf, main= paste0("Contour Plot for ", cell_lines[3], " MaxR"))
```

![](Drug-+-Radiation-Synergy-Analysis,-Figure-4,-S4_files/figure-gfm/contour-6.pdf)<!-- -->

------------------------------------------------------------------------

#### **Built with Version 4.4.3**
