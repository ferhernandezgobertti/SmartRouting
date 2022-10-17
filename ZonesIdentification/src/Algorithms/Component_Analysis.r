
# MAIN FUNCTION: processAnalysis
# Generates Principal Component Analysis results file with the identifier's variational
# spatial information based on processing eigenvectors.
# Inputs: inputFilename (existent CSV File), outputFilename (non-existent CSV File) - complete addresses
# Output: CSV file on specified outputfilename

new.getZonesData <- function(filename) {
    zonesData <- read.csv(file=filename, sep=",", header=T)
    return(zonesData)
}

new.splitZones <- function(table) {
    #splitTable %>% separate(table$zone, c("Alpha", "Beta", "Gamma", "Delta"))
    #splitTable %>% extract(table$zone, c("Alpha", "Beta", "Gamma", "Delta"), "([^-]+), ([^.]+), ([^*]+), ([^[0-9A-Z]]+)")
    splitTable <- str_split_fixed(table$zone, c("Alpha", "rest"), "-", 2)
    splitTable <- str_split_fixed(splitTable$rest, c("Beta", "rest"),".", 2)
    splitTable <- str_split_fixed(splitTable$rest, c("Gamma", "Delta"),"*", 2)
    return(splitTable)
}

new.getVariationalIds <- function(table, colNum, refChar) {
    varIdsTable <- table %>% mutate(colnames(table)[colNum] = as.numeric(as.character(colnames(table)[colNum])) - as.numeric(refChar))
    return(varIdsTable)
}

new.setVariationalColumns <- function(table) {
    table <- getVariationalIds(table, 3, 'A') #alpha_variation_pos
    table <- getVariationalIds(table, 3, 'Z') #alpha_variation_neg
    table <- getVariationalIds(table, 4, '1') #beta_variation_pos
    table <- getVariationalIds(table, 4, '9') #beta_variation_neg
    table <- getVariationalIds(table, 5, '1') #gamma_variation_pos
    table <- getVariationalIds(table, 5, '9') #gamma_variation_neg
    table <- getVariationalIds(table, 6, 'A') #delta_variation_pos
    table <- getVariationalIds(table, 6, 'Z') #delta_variation_neg
    return(table)
}

new.getVarData <- function(varTable) {
    alpha_variation_pos <- c(varTable[1,1]*varTable[1,3], varTable[2,1]*varTable[3,1])
    alpha_variation_neg <- c(varTable[1,1]*varTable[1,4], varTable[2,1]*varTable[4,1])
    beta_variation_pos  <- c(varTable[1,1]*varTable[1,5], varTable[2,1]*varTable[5,1])
    beta_variation_neg  <- c(varTable[1,1]*varTable[1,6], varTable[2,1]*varTable[6,1])
    gamma_variation_pos <- c(varTable[1,1]*varTable[1,7], varTable[2,1]*varTable[7,1])
    gamma_variation_neg <- c(varTable[1,1]*varTable[1,8], varTable[2,1]*varTable[8,1])
    delta_variation_pos <- c(varTable[1,1]*varTable[1,9], varTable[2,1]*varTable[9,1])
    delta_variation_neg <- c(varTable[1,1]*varTable[1,10], varTable[2,1]*varTable[10,1])
    varData <- data.frame(alpha_variation_pos, alpha_variation_neg, beta_variation_pos, beta_variation_neg, gamma_variation_pos, gamma_variation_neg, delta_variation_pos, delta_variation_neg)
    return(varData)
}

new.getPrincipalAnalysis <- function(table) {
    #varResults <- prcomp(table[,c(1,2,7:14)], center = TRUE, scale. = FALSE)
    varResults <- prcomp(table[,c(1:10)], center = TRUE, scale. = FALSE)
    summary(varResults)
    str(varResults)
    varData <- getVarData(varResults$rotation)
    return(varData)
}

new.writeAnalysisResult <- function(table, filename) {
    write.csv(table,filename, row.names = FALSE)
}

new.processAnalysis <- function(inputFilename, outputFilename) {
    library(dplyr)
    library(tidyr)
    library(stringr)

    zonesData <- getZonesData(inputFilename)
    zonesData <- splitZones(zonesData)
    zonesData <- setVariationalColumns(zonesData)

    variationalData <- getPrincipalAnalysis(zonesData)
    writeAnalysisResult(variationalData, outputFilename)
}