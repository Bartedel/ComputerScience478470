install.packages("stringr")
install.packages("tidyverse")
install.packages("stringi")
install.packages("stringdist")
install.packages("cluster")
install.packages("ade4")
install.packages("NMOF")
install.packages("rjson")
library(tidyverse)
library(stringr)
library(stringi)
library(stringdist)
library(cluster)
library(ade4)
library(NMOF)
library(rjson)
#importing data
Alldata <- unlist(fromJSON(file = "TVs-all-merged.json"), recursive = FALSE)
#setting parameters
keyTreshold <- 0.768
cosTreshold <- 0.634
titleTreshold <- 0
titleWeightPar <- 0.63
clusterTreshold <- 0.57
#cleaning data
test <- Alldata
for (i in 1: length(Alldata)){
test[[i]][["title"]] <- str_replace(test[[i]][["title"]], " inches", "inch")
test[[i]][["title"]] <- str_replace(test[[i]][["title"]], " inches", "inch")                   
test[[i]][["title"]] <- str_replace(test[[i]][["title"]], " hertz", "hz")
test[[i]][["title"]] <- str_replace(test[[i]][["title"]], " hz", "hz")
}

#making the model Word set
MWTitle = ""


for (i in 1:length(Alldata)){
     MWTitle <- c(MWTitle, strsplit(test[[i]][["title"]], " "))
}

MWTitle <- MWTitle[-1]

MWTitle2 <- ""

for (i in 1:length(MWTitle)){
  MWTitle2 <- c(MWTitle2, MWTitle[[i]])
}

MWTitle2 <- tolower(MWTitle2)
MWTitle2 <- str_replace(MWTitle2, "\"", "inch")
MWTitle2 <- str_replace(MWTitle2, "inches", "inch")
MWTitle2 <- str_replace(MWTitle2, "-inches", "inch")
MWTitle2 <- str_replace(MWTitle2, "-inch", "inch")
MWTitle2 <- str_replace(MWTitle2, "-in", "inch")
MWTitle2 <- str_replace(MWTitle2, "\\sin", "inch")
MWTitle2 <- str_replace(MWTitle2, "hertz", "hz")
MWTitle2 <- str_replace(MWTitle2, "-hertz", "hz")
MWTitle2 <- str_replace(MWTitle2, "-hz", "hz")
MWTitle2 <- str_replace(MWTitle2, " hz$", "hz")
MWTitle2 <- str_replace(MWTitle2, " in$", "inch")
MWTitle2 <- unique(MWTitle2)


for (i in 1:length(MWTitle2)){
  if (!grepl("[a-zA-Z0-9]*(([0-9]+[^0-9, ]+)|([^0-9, ]+[0-9]+))[a-zA-Z0-9]*", MWTitle2[i])){ 
    MWTitle2[i] <- ""
  }
}
MWTitle2 <- stri_remove_empty(MWTitle2)
#MWTitle2 <- c(brandList, MWTitle2)

binMatrix <- matrix(0, length(MWTitle2), length(MWTitle))
n <- 1


#filling in the Binary Matrix
for (i in 1: length(Alldata)){
    
    
    wordsTData <- strsplit(test[[i]][["title"]], " ")
    wordsT <- ""
    for (m in 1:length(wordsTData)){
      wordsT <- c(wordsT, wordsTData[[m]])
    }
    wordsT <- tolower(wordsT)
    wordsT <- str_replace(wordsT, "\"", "inch")
    wordsT <- str_replace(wordsT, "inches", "inch")
    wordsT <- str_replace(wordsT, "-inches", "inch")
    wordsT <- str_replace(wordsT, "-inch", "inch")
    wordsT <- str_replace(wordsT, "-in", "inch")
    wordsT <- str_replace(wordsT, "\\sin", "inch")
    wordsT <- str_replace(wordsT, "hertz", "hz")
    wordsT <- str_replace(wordsT, "-hertz", "hz")
    wordsT <- str_replace(wordsT, "-hz", "hz")
    wordsT <- str_replace(wordsT, " hz$", "hz")
    wordsT <- str_replace(wordsT, " in$", "inch")
    wordsT <- unique(wordsT)
    wordsT <- str_remove_all(wordsT, "\\W$|^\\W")
    wordsT <- stri_remove_empty(wordsT)
    for (k in 1:length(MWTitle2)){

      if (is.element(MWTitle2[k], wordsT)){
        binMatrix[k,n] <- 1
      }
    }
    n <- n + 1
}

#Computing the signature matrix using minhashes
rowNr <- 1:length(MWTitle2)
minHash <- matrix(NA, length(MWTitle2), 684)
for (i in 1:684){
  a <- sample(1:1348, 1, replace=TRUE)
  b <- sample(1:1348, 1, replace=TRUE)
  minHash[,i] <- (a*rowNr+b)%%length(MWTitle2)
}

sigMatrix <- matrix((length(MWTitle2)+1), 684, 1624)



for (i in 1:length(rowNr)){
  for (j in 1:length(MWTitle)){
    if (binMatrix[i,j] == 1){
      for (k in 1:684){
        sigMatrix[k,j] <- min(sigMatrix[k,j],minHash[i,k])
      }
    }
  }
}


#LSH
r <- 18
b <- 684/r 
n <- 1
candidatePairs <- matrix(0, 1624, 1624)

for (i in 1:b){
  print(i)
  hashValue <- matrix(NA, 1624, 1)
  for (j in 1:length(MWTitle)){
    hashValue[j] <- paste(sigMatrix[n:(n+r-1),j], collapse = "")
  }
  simValues <- hashValue[duplicated(hashValue)]
  for (k in 1: length(simValues)){
    simValues2 <- which(hashValue == simValues[k])
    for (l in 1:length(simValues2)){
      for (m in 1:length(simValues2)){
        candidatePairs[simValues2[l], simValues2[m]] <- 1
      }
    }
  }
  n <- n + r
}
counter1 <- 0
counter2 <- 0

for(i in 1: 1623){
  for (j in i:1624){
    if(candidatePairs[i,j] == 1 && Alldata[[i]][["modelID"]] == Alldata[[j]][["modelID"]] && i != j){
      counter1 <- counter1 + 1
    }
    if (candidatePairs[i,j] == 1 && i != j){
      counter2 <- counter2 + 1
    }
  }
  
}

#Evaluation LSH
Ratio <- counter2/1317876
print(Ratio)
PQlsh <- counter1/counter2
PClsh <- counter1/362
Fstar <- (PQlsh*PClsh*2)/(PQlsh+PClsh)
print(PQlsh)
print(PClsh)
print(Fstar)


#check if there are different brands
difBrands <- matrix(1, 1624, 1624)

for (i in 1:1624){
  for (j in 1:1624){
    if (identical(is.element(brandList, tolower(c(MWTitle[[i]]))), is.element(brandList, tolower(c(MWTitle[[j]]))))){
      difBrands[i,j] <- 0
    }
  }
}
# check if a shop is the same
sameShop <-matrix(0, 1624, 1624)
webShop <- list()

for (i in 1:length(Alldata)){
    webShop <- c(webShop, Alldata[[i]][["shop"]])
}

for (i in 1:1624){
  for (j in 1:1624){
    if (identical(webShop[i], webShop[j])){
      sameShop[i,j] <- 1
    }
  }
}


#MSM
distanceFinal <- matrix(NA, 1624, 1624)

for (i in 1:(length(Alldata))){
  print(i)
  for (j in 1:length(Alldata)){
    if (candidatePairs[i,j] == 0 || sameShop[i,j] == 1 || difBrands == 1){
      distanceFinal[i,j] <- 10000
    }else{
      sim = 0
      avgSim = 0
      z = 0
      w = 0
      nmki <- tolower(Alldata[[i]][["featuresMap"]])
      nmkj <- tolower(Alldata[[j]][["featuresMap"]])
      nmki <- str_remove_all(nmki, "\\W")
      nmkj <- str_remove_all(nmkj, "\\W")
      
      for (k in 1: length(nmki)){
        for (l in 1: length(nmkj)){
          simKey <- stringsim(str_pad(str_remove_all(tolower(names(Alldata[[i]][["featuresMap"]]))[k], "\\W"),3), str_pad(str_remove_all(tolower(names(Alldata[[j]][["featuresMap"]][l])), "\\W"),3), method = "qgram", q=3)
          if (simKey > keyTreshold){
            valueSim <- stringsim(str_pad(nmki[k],3), str_pad(nmkj[l],3), method = "qgram", q=3)
            weight <- simKey
            sim <- sim + weight * valueSim
            z <- z + 1
            w = w + weight
            nmki2 <- nmki[-k]
            nmkj2 <- nmkj[-l]
            
          }
        }
      }
      if (w > 0){
        avgSim = sim/w
      }
      mwPerc <- length(intersect(nmki2, nmkj2))/(length(nmki2)+length(nmkj2))
      titleSimCos <- stringsim(str_remove_all(tolower(Alldata[[i]][["title"]]), "\\W"),str_remove_all(tolower(Alldata[[j]][["title"]]), "\\W"), method = "cosine")
      if(titleSimCos > cosTreshold){
        titleSim <- 1
      }else{
        wordsTiData <- strsplit(Alldata[[i]][["title"]], " ")
        wordsTi <- ""
        
        for (m in 1:length(wordsTiData)){
          wordsTi <- c(wordsTi, wordsTiData[[m]])
        }
        wordsTi <- c(tolower(wordsTi))
        for (m in 1:length(wordsTi)){
          if (!grepl("^[0-9]+[a-z]+$", wordsTi[m])){ 
            wordsTi[m] <- ""
          }
        }
        
        wordsTi <- stri_remove_empty(wordsTi)
        
        wordsTjData <- strsplit(Alldata[[j]][["title"]], " ")
        wordsTj <- ""
        
        for (m in 1:length(wordsTjData)){
          wordsTj <- c(wordsTj, wordsTjData[[m]])
        }
        wordsTj <- c(tolower(wordsTj))
        for (m in 1:length(wordsTj)){
          if (!grepl("^[0-9]+[a-z]+$", wordsTj[m])){ 
            wordsTj[m] <- ""
          }
        }
        
        wordsTj <- stri_remove_empty(wordsTj)
        levSim <- 0
        counter <- 0
        for (m in 1:length(wordsTi)){
          for (n in 1:length(wordsTj)){
            if(!is_empty(wordsTj) && !is_empty(wordsTi) && stringdist(str_extract(wordsTi[m], "[a-z]+"), str_extract(wordsTj[n], "[a-z]+"), method = "lv") < 0.5 && str_extract(wordsTi[m], "[0-9]+") == str_extract(wordsTj[n], "[0-9]+")){
              counter <- 1
            }
            levSim <- levSim + stringsim(wordsTi[m], wordsTj[n], method = "lv")
          }
        }
        aveLevSim <- levSim/(length(wordsTi) * length(wordsTj))
        preTitleSim <- (aveLevSim + titleSimCos)/ 2
        if (counter == 1 && preTitleSim > titleTreshold){
          titleSim <- preTitleSim
        }else{
          titleSim <- -1
        }

      }
      if (titleSim == -1){
        theta1 <- z/min(length(nmki), length(nmkj))
        theta2 <- 1- theta1
        distanceFinal[i,j] <- 1 - (theta1 * avgSim + theta2 * mwPerc)
      }else{
        theta1 <- (1-titleWeightPar)* (z/min(length(nmki), length(nmkj)))
        theta2 <- 1 - titleWeightPar - theta1
        distanceFinal[i,j] <- 1 - (theta1 * avgSim + theta2 * mwPerc + titleWeightPar * titleSim)
      }
    }
  }
}

#Hierarchal Clustering  
clustDist <- as.dist(distanceFinal)
cluster <- hclust(clustDist, method = "complete")
clusterFinal <- cutree(cluster, h=clusterTreshold)
duplCluster <- clusterFinal[duplicated(clusterFinal)]
clusteredTV <- vector(mode = "list", length = length(duplCluster))
for (i in 1:length(duplCluster)){
  clusteredTV[[i]] <- which(clusterFinal %in% duplCluster[i])
}
clusteredTV <- unique(clusteredTV)

FP <- 0
counterCluster <- matrix(0, length(clusteredTV), 1)
for(i in 1: length(clusteredTV)){
  for ( j in 1: length(clusteredTV[[i]])){
    for ( l in 1:length(clusteredTV[[i]])){
      if(Alldata[[clusteredTV[[i]][j]]][["modelID"]] == Alldata[[clusteredTV[[i]][l]]][["modelID"]] && j != l && length(clusteredTV[[i]] == 2)){
        counterCluster[i] <- counterCluster[i] + 1
      }
      if(j != l && length(clusteredTV[[i]] == 2)){
        FP <- FP + 1
      }
    }
  }
}

#Evaluation MSM
FP <- FP/2
TP <- sum(counterCluster)/2
Precision <- TP/(TP+FP)
Recall <- TP/362
F <- (2*Precision*Recall)/(Precision+Recall)
print(F)

