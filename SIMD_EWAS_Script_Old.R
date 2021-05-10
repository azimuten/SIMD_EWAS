setwd("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/")

############## Raw Model ####################################



##Creating phenotype file for SIMD
View(genscot_merged)
SIMD_pheno <- genscot_merged[,c("ID", "rank")]
colnames(SIMD_pheno)[1] <- "id"
View(SIMD_pheno)
write.csv(SIMD_pheno, file = "SIMD_pheno.csv")

wave1_modelE_cov <- readRDS("./wave1_modelE_covariates_5087_20181127.rds")
wave3_modelE_cov <- readRDS("./wave3_modelE_covariates_4450_20191120.rds")



## Creating a custom covariate file for wave 3 EWAS
setwd("C:/Users/ajesp/Documents/PhD/SIMD_EWAS")
genscot_merged <- read.csv("genscot_merged.csv")
cellcountinfo <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/cellcountinfo.csv", header=TRUE, sep=",")
cellcountinfo_w3 <- read.table("C:/Users/ajesp/Documents/PhD/MDD_inflam_EWAS//samplesheet.final.csv", header=TRUE, sep=",")
IDinfo <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/IDinfo.csv", header=TRUE, sep=",")
Batch_w3 <- read.csv("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/Batch_w3.csv")
colnames(Batch_w3)[1] <- "ID"



#Match IDs of cellcountinfo and IDinfo
IDinfo$Sentrix_ID <- sub("^", "X", IDinfo$Sentrix_ID )
IDinfo$Sentrix_Position <- sub("^", "_", IDinfo$Sentrix_Position )
IDinfo$SampleID <- paste0(IDinfo$Sentrix_ID, IDinfo$Sentrix_Position)

#Merge IDinfo and cellountinfo
cellcountinfo <- merge(cellcountinfo, IDinfo, by = "SampleID", all = TRUE)

#Rename SampleID
colnames(cellcountinfo)[8] <- "ID"
colnames(cellcountinfo_w3)[1] <- "ID"

#concatenate cellcount info from wave 1 and wave 3
cellcount_all <- rbind(cellcountinfo[,c(2:8)], cellcountinfo_w3[,c(1,8:13)])

#Rename cellcount_all ID column
colnames(cellcount_all)[5] <- "id"

#Merge cellcount, age, smoking, cell counts and batch number
SIMD_cov <- merge(genscot_merged[,c(2,3,4,16,17)], cellcount_all[,c(1:4,6,7)], by = "ID")
SIMD_cov <- merge(SIMD_cov, Batch_w3, by = "ID")


#rename ID column to id
colnames(SIMD_cov)[1] <- "id"

#save as RDS
saveRDS(SIMD_cov, file = "SIMD_cov.rds")



setwd("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/")

##opening Ewas output
SIMD_W1_TopTables <- read.delim("SIMD_rank_W1.toptable.txt")
SIMD_W3_TopTables <- read.delim("SIMD_rank_w3.toptable.txt")

View(SIMD_W1_TopTables)
View(SIMD_W3_TopTables)



## Creating a presentation table
SIMD_pres_table <- head(SIMD_W1_TopTables[order(SIMD_W1_TopTables$P.Value), ], 10)
View(SIMD_pres_table)

write.csv(SIMD_pres_table, file = "SIMD_W1_Top10.csv")


SIMD_pres_table_w3 <- head(SIMD_W3_TopTables[order(SIMD_W3_TopTables$P.Value), ], 10)
View(SIMD_pres_table_w3)

write.csv(SIMD_pres_table_w3, file = "SIMD_W3_Top10.csv")


## Manhattan plot
SIMD_man_W1 <- SIMD_W1_TopTables[,c(1,2,3,4,10)]
colnames(SIMD_man_W1)[1] <- "SNP"
colnames(SIMD_man_W1)[3] <- "CHR"
colnames(SIMD_man_W1)[4] <- "BP"
colnames(SIMD_man_W1)[5] <- "P"

SIMD_man_W1$CHR <- gsub("chr", "", SIMD_man_W1$CHR)
SIMD_man_W1 <- transform(SIMD_man_W1, CHR = as.numeric(CHR))



manhattan(SIMD_man_W1, col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-9)), main = "SIMD W1 (Smoking)", cex.main = 1.5) #export with width 800 and keep aspect ratio



SIMD_man_W3 <- SIMD_W3_TopTables[,c(1,2,3,4,10)]
colnames(SIMD_man_W3)[1] <- "SNP"
colnames(SIMD_man_W3)[3] <- "CHR"
colnames(SIMD_man_W3)[4] <- "BP"
colnames(SIMD_man_W3)[5] <- "P"

SIMD_man_W3$CHR <- gsub("chr", "", SIMD_man_W3$CHR)
SIMD_man_W3 <- transform(SIMD_man_W3, CHR = as.numeric(CHR))



manhattan(SIMD_man_W3, col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-09)), main = "SIMD W3 (Smoking)", cex.main = 1.5) #export with width 800 and keep aspect ratio






##Meta analysis

#creating wave 1 dataframes for METAL meta analysis
Metal_SIMD_w1 <- SIMD_W1_TopTables[, c("ID", "P.Value", "N", "beta")]
View(Metal_SIMD_w1)

#adding reference allele columns to the metal data frame
Metal_SIMD_w1$A1 <- rep("A",nrow(Metal_SIMD_w1))
Metal_SIMD_w1$A2 <- rep("C",nrow(Metal_SIMD_w1))

#adding a direction of effect column
Metal_SIMD_w1$Effect <- ifelse(Metal_SIMD_w1$beta >= 0, "+", "-")

#save dataframes 
write.csv(Metal_SIMD_w1, file = "Metal_SIMD_w1.csv", row.names = FALSE, quote = FALSE)



#creating wave 3 dataframes for METAL meta analysis
Metal_SIMD_w3 <- SIMD_W3_TopTables[, c("ID", "P.Value", "N", "beta")]
View(Metal_SIMD_w3)

#adding reference allele columns to the metal data frame
Metal_SIMD_w3$A1 <- rep("A",nrow(Metal_SIMD_w3))
Metal_SIMD_w3$A2 <- rep("C",nrow(Metal_SIMD_w3))

#adding a direction of effect column
Metal_SIMD_w3$Effect <- ifelse(Metal_SIMD_w3$beta >= 0, "+", "-")

#save dataframes 
write.csv(Metal_SIMD_w3, file = "Metal_SIMD_w3.csv", row.names = FALSE, quote = FALSE)


#opening metal output
SIMD_meta_Output <- read.table("SIMD_meta1.TBL", header = TRUE, fill = TRUE)
View(SIMD_meta_Output)

#merging in gene annotation, chromosome and map info from both w1 and w3 dataframes
SIMD_Meta_TopTable <- merge(SIMD_meta_Output[,c("MarkerName", "Weight", "Zscore", "P.value", "Direction")], 
                                      SIMD_W1_TopTables[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_TopTable <- merge(SIMD_Meta_TopTable, SIMD_W3_TopTables[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_TopTable$geneSymbol.x <- as.character(SIMD_Meta_TopTable$geneSymbol.x)
SIMD_Meta_TopTable$geneSymbol.y <- as.character(SIMD_Meta_TopTable$geneSymbol.y)

SIMD_Meta_TopTable$geneSymbol.x[ is.na(SIMD_Meta_TopTable$geneSymbol.x) ] <- SIMD_Meta_TopTable$geneSymbol.y[ is.na(SIMD_Meta_TopTable$geneSymbol.x) ]
SIMD_Meta_TopTable$geneSymbol.y =NULL

SIMD_Meta_TopTable$CHR.x[ is.na(SIMD_Meta_TopTable$CHR.x) ] <- SIMD_Meta_TopTable$CHR.y[ is.na(SIMD_Meta_TopTable$CHR.x) ]
SIMD_Meta_TopTable$CHR.y =NULL

SIMD_Meta_TopTable$MAPINFO.x[ is.na(SIMD_Meta_TopTable$MAPINFO.x) ] <- SIMD_Meta_TopTable$MAPINFO.y[ is.na(SIMD_Meta_TopTable$MAPINFO.x) ]
SIMD_Meta_TopTable$MAPINFO.y =NULL

SIMD_Meta_TopTable$FEATURE.x <- as.character(SIMD_Meta_TopTable$FEATURE.x)
SIMD_Meta_TopTable$FEATURE.y <- as.character(SIMD_Meta_TopTable$FEATURE.y)

SIMD_Meta_TopTable$FEATURE.x[ is.na(SIMD_Meta_TopTable$FEATURE.x) ] <- SIMD_Meta_TopTable$FEATURE.y[ is.na(SIMD_Meta_TopTable$FEATURE.x) ]
SIMD_Meta_TopTable$FEATURE.y =NULL

SIMD_Meta_TopTable$CpGISLAND.x[ is.na(SIMD_Meta_TopTable$CpGISLAND.x) ] <- SIMD_Meta_TopTable$CpGISLAND.y[ is.na(SIMD_Meta_TopTable$CpGISLAND.x) ]
SIMD_Meta_TopTable$CpGISLAND.y =NULL

colnames(SIMD_Meta_TopTable)[6] <- "geneSymbol"
colnames(SIMD_Meta_TopTable)[7] <- "CHR"
colnames(SIMD_Meta_TopTable)[8] <- "MAPINFO"
colnames(SIMD_Meta_TopTable)[9] <- "FEATURE"
colnames(SIMD_Meta_TopTable)[10] <- "CpGISLAND"

View(SIMD_Meta_TopTable)

write.csv(SIMD_Meta_TopTable, file = "SIMD_Meta_TopTable.csv")
SIMD_Meta_TopTable <- read.csv("SIMD_Meta_TopTable.csv")


SIMD_meta_pres_table <- head(SIMD_Meta_TopTable[order(SIMD_Meta_TopTable$P.value), ], 20)
View(SIMD_meta_pres_table)

write.csv(SIMD_meta_pres_table, file = "SIMD_Meta_Top20.csv")

#creating a manhattan plot of the meta analysis
SIMD_Meta_Man <- SIMD_Meta_TopTable[,c(1,4,6,7,8)]
colnames(SIMD_Meta_Man)[1] <- "SNP"
colnames(SIMD_Meta_Man)[2] <- "P"
colnames(SIMD_Meta_Man)[4] <- "CHR"
colnames(SIMD_Meta_Man)[5] <- "BP"


SIMD_Meta_Man$CHR <- gsub("chr", "", SIMD_Meta_Man$CHR)
SIMD_Meta_Man <- transform(SIMD_Meta_Man, CHR = as.numeric(CHR))
SIMD_Meta_Man <- transform(SIMD_Meta_Man, BP = as.numeric(BP))
SIMD_Meta_Man <- na.omit(SIMD_Meta_Man)

View(SIMD_Meta_Man)

manhattan((SIMD_Meta_Man), col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-14)), main = "SIMD Meta (Smoking)", cex.main = 1.5) #export with width 800 and keep aspect ratio









################### Never Smoked RAW Model ##############################################################################

setwd("C:/Users/ajesp/Documents/PhD/SIMD_EWAS")
genscot_merged <- read.csv("genscot_merged.csv")


##Creating a custom covariate file for NeverSmoked wave 3 EWAS
cellcountinfo <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/cellcountinfo.csv", header=TRUE, sep=",")
cellcountinfo_w3 <- read.table("C:/Users/ajesp/Documents/PhD/MDD_inflam_EWAS//samplesheet.final.csv", header=TRUE, sep=",")
IDinfo <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/IDinfo.csv", header=TRUE, sep=",")
Batch_w3 <- read.csv("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/Batch_w3.csv")
colnames(Batch_w3)[1] <- "ID"

#Match IDs of cellcountinfo and IDinfo
IDinfo$Sentrix_ID <- sub("^", "X", IDinfo$Sentrix_ID )
IDinfo$Sentrix_Position <- sub("^", "_", IDinfo$Sentrix_Position )
IDinfo$SampleID <- paste0(IDinfo$Sentrix_ID, IDinfo$Sentrix_Position)

#Merge IDinfo and cellountinfo
cellcountinfo <- merge(cellcountinfo, IDinfo, by = "SampleID", all = TRUE)

#Rename SampleID
colnames(cellcountinfo)[8] <- "ID"
colnames(cellcountinfo_w3)[1] <- "ID"

#concatenate cellcount info from wave 1 and wave 3
cellcount_all <- rbind(cellcountinfo[,c(2:8)], cellcountinfo_w3[,c(1,8:13)])

#Merge cellcount, age, cell counts and batch number
SIMD_cov_NS_w3 <- merge(genscot_merged[,c(2,16,17)], cellcount_all[,c(1:4,6,7)], by = "ID")
SIMD_cov_NS_w3 <- merge(SIMD_cov_NS_w3, Batch_w3, by = "ID")


#rename ID column to id
colnames(SIMD_cov_NS_w3)[1] <- "id"

#save as RDS
saveRDS(SIMD_cov_NS_w3, file = "SIMD_cov_NS_w3.rds")

##creating a custom Phenotype file for NS EWAS
SIMD_NeverSmoke <- genscot_merged[,c("ID", "rank", "ever_smoke")]
SIMD_NeverSmoke$never_smoke <- ifelse(SIMD_NeverSmoke$ever_smoke == 4,1,0)
SIMD_NeverSmoke <- na.omit(smoketest[!(SIMD_NeverSmoke$never_smoke == 0),])
                           
names(SIMD_NeverSmoke)[1] <- "id"
write.csv(SIMD_NeverSmoke, file = "SIMD_NeverSmoke.csv")


##Custom covariate file for never smoked w1
SIMD_NS_Cov_w1 <- genscot_merged[, c("ID", "age", "sex")]
names(SIMD_NS_Cov_w1)[1] <- "id"
saveRDS(SIMD_NS_Cov_w1, file = "SIMD_NS_Cov.rds")


#Opening EWAS toptables of w1 and w3
SIMD_NS_W1_TopTables <- read.delim("SIMD_NS_W1.toptable.txt")
SIMD_NS_W3_TopTables <- read.delim("SIMD_NS_w3.toptable.txt")
View(SIMD_NS_W1_TopTables)
View(SIMD_NS_W3_TopTables)

##Meta analysis

#creating wave 1 dataframes for METAL meta analysis
METAL_SIMD_NS_W1 <- SIMD_NS_W1_TopTables[, c("ID", "P.Value", "N", "beta")]
View(METAL_SIMD_NS_W1)

#adding reference allele columns to the metal data frame
METAL_SIMD_NS_W1$A1 <- rep("A",nrow(METAL_SIMD_NS_W1))
METAL_SIMD_NS_W1$A2 <- rep("C",nrow(METAL_SIMD_NS_W1))

#adding a direction of effect column
METAL_SIMD_NS_W1$Effect <- ifelse(METAL_SIMD_NS_W1$beta >= 0, "+", "-")

#creating wave 3 dataframes for METAL meta analysis
METAL_SIMD_NS_W3 <- SIMD_NS_W3_TopTables[, c("ID", "P.Value", "N", "beta")]
View(METAL_SIMD_NS_W3)

#adding reference allele columns to the metal data frame
METAL_SIMD_NS_W3$A1 <- rep("A",nrow(METAL_SIMD_NS_W3))
METAL_SIMD_NS_W3$A2 <- rep("C",nrow(METAL_SIMD_NS_W3))

#adding a direction of effect column
METAL_SIMD_NS_W3$Effect <- ifelse(METAL_SIMD_NS_W3$beta >= 0, "+", "-")

#save dataframes 
write.csv(METAL_SIMD_NS_W1, file = "METAL_SIMD_NS_W1.csv", row.names = FALSE, quote = FALSE)
write.csv(METAL_SIMD_NS_W3, file = "METAL_SIMD_NS_W3.csv", row.names = FALSE, quote = FALSE)

#opening metal output
SIMD_NS_METAL_Output <- read.table("SIMD_rank_NS_Meta.TBL", header = TRUE)
View(SIMD_NS_METAL_Output)

#merging in gene annotation, chromosome and map info from both w1 and w3 dataframes
SIMD_Meta_NS_TopTable <- merge(SIMD_NS_METAL_Output[,c("MarkerName", "Weight", "Zscore", "P.value", "Direction")], 
                                SIMD_NS_W1_TopTables[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_NS_TopTable <- merge(SIMD_Meta_NS_TopTable, SIMD_NS_W3_TopTables[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_NS_TopTable$geneSymbol.x <- as.character(SIMD_Meta_NS_TopTable$geneSymbol.x)
SIMD_Meta_NS_TopTable$geneSymbol.y <- as.character(SIMD_Meta_NS_TopTable$geneSymbol.y)

SIMD_Meta_NS_TopTable$geneSymbol.x[ is.na(SIMD_Meta_NS_TopTable$geneSymbol.x) ] <- SIMD_Meta_NS_TopTable$geneSymbol.y[ is.na(SIMD_Meta_NS_TopTable$geneSymbol.x) ]
SIMD_Meta_NS_TopTable$geneSymbol.y =NULL

SIMD_Meta_NS_TopTable$CHR.x[ is.na(SIMD_Meta_NS_TopTable$CHR.x) ] <- SIMD_Meta_NS_TopTable$CHR.y[ is.na(SIMD_Meta_NS_TopTable$CHR.x) ]
SIMD_Meta_NS_TopTable$CHR.y =NULL

SIMD_Meta_NS_TopTable$MAPINFO.x[ is.na(SIMD_Meta_NS_TopTable$MAPINFO.x) ] <- SIMD_Meta_NS_TopTable$MAPINFO.y[ is.na(SIMD_Meta_NS_TopTable$MAPINFO.x) ]
SIMD_Meta_NS_TopTable$MAPINFO.y =NULL

SIMD_Meta_NS_TopTable$FEATURE.x <- as.character(SIMD_Meta_NS_TopTable$FEATURE.x)
SIMD_Meta_NS_TopTable$FEATURE.y <- as.character(SIMD_Meta_NS_TopTable$FEATURE.y)

SIMD_Meta_NS_TopTable$FEATURE.x[ is.na(SIMD_Meta_NS_TopTable$FEATURE.x) ] <- SIMD_Meta_NS_TopTable$FEATURE.y[ is.na(SIMD_Meta_NS_TopTable$FEATURE.x) ]
SIMD_Meta_NS_TopTable$FEATURE.y =NULL

SIMD_Meta_NS_TopTable$CpGISLAND.x[ is.na(SIMD_Meta_NS_TopTable$CpGISLAND.x) ] <- SIMD_Meta_NS_TopTable$CpGISLAND.y[ is.na(SIMD_Meta_NS_TopTable$CpGISLAND.x) ]
SIMD_Meta_NS_TopTable$CpGISLAND.y =NULL

colnames(SIMD_Meta_NS_TopTable)[6] <- "geneSymbol"
colnames(SIMD_Meta_NS_TopTable)[7] <- "CHR"
colnames(SIMD_Meta_NS_TopTable)[8] <- "MAPINFO"
colnames(SIMD_Meta_NS_TopTable)[9] <- "FEATURE"
colnames(SIMD_Meta_NS_TopTable)[10] <- "CpGISLAND"

View(SIMD_Meta_NS_TopTable)

SIMD_NS_meta_pres_table <- head(SIMD_Meta_NS_TopTable[order(SIMD_Meta_NS_TopTable$P.value), ], 10)
View(SIMD_pres_table_w3)

write.csv(SIMD_NS_meta_pres_table, file = "SIMD_NS_meta_Top10.csv")

#creating a manhattan plot of the meta analysis
SIMD_NS_Meta_Man <- SIMD_Meta_NS_TopTable[,c(1,4,6,7,8)]
colnames(SIMD_NS_Meta_Man)[1] <- "SNP"
colnames(SIMD_NS_Meta_Man)[2] <- "P"
colnames(SIMD_NS_Meta_Man)[4] <- "CHR"
colnames(SIMD_NS_Meta_Man)[5] <- "BP"


SIMD_NS_Meta_Man$CHR <- gsub("chr", "", SIMD_NS_Meta_Man$CHR)
SIMD_NS_Meta_Man <- transform(SIMD_NS_Meta_Man, CHR = as.numeric(CHR))



manhattan(na.omit(SIMD_NS_Meta_Man), col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-10)), main = "SIMD rank in never smoked (Meta)", cex.main = 1.5) #export with width 800 and keep aspect ratio








################### Never Smoked BMI_Alc Model ##############################################################################

setwd("C:/Users/ajesp/Documents/PhD/SIMD_EWAS")
genscot_merged <- read.csv("genscot_merged.csv")
cellcountinfo <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/cellcountinfo.csv", header=TRUE, sep=",")
cellcountinfo_w3 <- read.table("C:/Users/ajesp/Documents/PhD/MDD_inflam_EWAS//samplesheet.final.csv", header=TRUE, sep=",")
IDinfo <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/IDinfo.csv", header=TRUE, sep=",")
Alcohol <- read.csv("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/alcohol.csv")
colnames(Alcohol)[1] <- "ID"
Batch_w3 <- read.csv("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/Batch_w3.csv")
colnames(Batch_w3)[1] <- "ID"


##For Wave 1
#Merge age, sex, BMI, and alcohol
SIMD_NS_Cov_BMI_Alc_w1 <- merge(genscot_merged[,c(2,16,17,28)], Alcohol[,c(1,4)], by = "ID")

#rename ID column to id
colnames(SIMD_NS_Cov_BMI_Alc_w1)[1] <- "id"

#save as RDS
saveRDS(SIMD_NS_Cov_BMI_Alc_w1, file = "SIMD_NS_Cov_BMI_Alc_w1.rds")


##For Wave 3 
#Match IDs of cellcountinfo and IDinfo
IDinfo$Sentrix_ID <- sub("^", "X", IDinfo$Sentrix_ID )
IDinfo$Sentrix_Position <- sub("^", "_", IDinfo$Sentrix_Position )
IDinfo$SampleID <- paste0(IDinfo$Sentrix_ID, IDinfo$Sentrix_Position)

#Merge IDinfo and cellountinfo
cellcountinfo <- merge(cellcountinfo, IDinfo, by = "SampleID", all = TRUE)

#Rename SampleID
colnames(cellcountinfo)[8] <- "ID"
colnames(cellcountinfo_w3)[1] <- "ID"

#concatenate cellcount info from wave 1 and wave 3
cellcount_all <- rbind(cellcountinfo[,c(2:8)], cellcountinfo_w3[,c(1,8:13)])

#Merge smoking, age, sex, BMI, cellcount, alcohol and batch number
SIMD_NS_Cov_BMI_Alc_w3 <- merge(genscot_merged[,c(2,16,17,28)], cellcount_all[,c(1:4,6,7)], by = "ID")
SIMD_NS_Cov_BMI_Alc_w3 <- merge(SIMD_NS_Cov_BMI_Alc_w3, Alcohol[,c(1,4)], by = "ID")
SIMD_NS_Cov_BMI_Alc_w3 <- merge(SIMD_NS_Cov_BMI_Alc_w3, Batch_w3, by = "ID")

#rename ID column to id
colnames(SIMD_NS_Cov_BMI_Alc_w3)[1] <- "id"

View(SIMD_NS_Cov_BMI_Alc_w3)

#save as RDS
saveRDS(SIMD_NS_Cov_BMI_Alc_w3, file = "SIMD_NS_Cov_BMI_Alc_w3.rds")



SIMD_NS_BMI_Alc_W1_TopTables <- read.delim("SIMD_NS_BMI_Alc_W1.toptable.txt")
SIMD_NS_BMI_Alc_W3_TopTables <- read.delim("SIMD_NS_BMI_Alc_w3.toptable.txt")
View(SIMD_NS_BMI_Alc_W1_TopTables)
View(SIMD_NS_BMI_Alc_W3_TopTables)

##Meta analysis

#creating wave 1 dataframes for METAL meta analysis
METAL_SIMD_NS_BMI_Alc_W1 <- SIMD_NS_BMI_Alc_W1_TopTables[, c("ID", "P.Value", "N", "beta")]

#adding reference allele columns to the metal data frame
METAL_SIMD_NS_BMI_Alc_W1$A1 <- rep("A",nrow(METAL_SIMD_NS_BMI_Alc_W1))
METAL_SIMD_NS_BMI_Alc_W1$A2 <- rep("C",nrow(METAL_SIMD_NS_BMI_Alc_W1))

#adding a direction of effect column
METAL_SIMD_NS_BMI_Alc_W1$Effect <- ifelse(METAL_SIMD_NS_BMI_Alc_W1$beta >= 0, "+", "-")

View(METAL_SIMD_NS_BMI_Alc_W1)


#creating wave 3 dataframes for METAL meta analysis
METAL_SIMD_NS_BMI_Alc_W3 <- SIMD_NS_BMI_Alc_W3_TopTables[, c("ID", "P.Value", "N", "beta")]

#adding reference allele columns to the metal data frame
METAL_SIMD_NS_BMI_Alc_W3$A1 <- rep("A",nrow(METAL_SIMD_NS_BMI_Alc_W3))
METAL_SIMD_NS_BMI_Alc_W3$A2 <- rep("C",nrow(METAL_SIMD_NS_BMI_Alc_W3))

#adding a direction of effect column
METAL_SIMD_NS_BMI_Alc_W3$Effect <- ifelse(METAL_SIMD_NS_BMI_Alc_W3$beta >= 0, "+", "-")

View(METAL_SIMD_NS_BMI_Alc_W3)


#save dataframes 
write.csv(METAL_SIMD_NS_BMI_Alc_W1, file = "METAL_SIMD_NS_BMI_Alc_W1.csv", row.names = FALSE, quote = FALSE)
write.csv(METAL_SIMD_NS_BMI_Alc_W3, file = "METAL_SIMD_NS_BMI_Alc_W3.csv", row.names = FALSE, quote = FALSE)

#opening metal output
SIMD_NS_BMI_Alc_METAL_Output <- read.table("Meta_SIMD_NS_BMI_Alc1.tbl", header = TRUE)
View(SIMD_NS_BMI_Alc_METAL_Output)

#merging in gene annotation, chromosome and map info from both w1 and w3 dataframes
SIMD_Meta_NS_BMI_Alc_TopTable <- merge(SIMD_NS_BMI_Alc_METAL_Output[,c("MarkerName", "Weight", "Zscore", "P.value", "Direction")],
                                       SIMD_NS_BMI_Alc_W1_TopTables[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_NS_BMI_Alc_TopTable <- merge(SIMD_Meta_NS_BMI_Alc_TopTable, SIMD_NS_BMI_Alc_W3_TopTables[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_NS_BMI_Alc_TopTable$geneSymbol.x <- as.character(SIMD_Meta_NS_BMI_Alc_TopTable$geneSymbol.x)
SIMD_Meta_NS_BMI_Alc_TopTable$geneSymbol.y <- as.character(SIMD_Meta_NS_BMI_Alc_TopTable$geneSymbol.y)

SIMD_Meta_NS_BMI_Alc_TopTable$geneSymbol.x[ is.na(SIMD_Meta_NS_BMI_Alc_TopTable$geneSymbol.x) ] <- SIMD_Meta_NS_BMI_Alc_TopTable$geneSymbol.y[ is.na(SIMD_Meta_NS_BMI_Alc_TopTable$geneSymbol.x) ]
SIMD_Meta_NS_BMI_Alc_TopTable$geneSymbol.y =NULL

SIMD_Meta_NS_BMI_Alc_TopTable$CHR.x[ is.na(SIMD_Meta_NS_BMI_Alc_TopTable$CHR.x) ] <- SIMD_Meta_NS_BMI_Alc_TopTable$CHR.y[ is.na(SIMD_Meta_NS_BMI_Alc_TopTable$CHR.x) ]
SIMD_Meta_NS_BMI_Alc_TopTable$CHR.y =NULL

SIMD_Meta_NS_BMI_Alc_TopTable$MAPINFO.x[ is.na(SIMD_Meta_NS_BMI_Alc_TopTable$MAPINFO.x) ] <- SIMD_Meta_NS_BMI_Alc_TopTable$MAPINFO.y[ is.na(SIMD_Meta_NS_BMI_Alc_TopTable$MAPINFO.x) ]
SIMD_Meta_NS_BMI_Alc_TopTable$MAPINFO.y =NULL

SIMD_Meta_NS_BMI_Alc_TopTable$FEATURE.x <- as.character(SIMD_Meta_NS_BMI_Alc_TopTable$FEATURE.x)
SIMD_Meta_NS_BMI_Alc_TopTable$FEATURE.y <- as.character(SIMD_Meta_NS_BMI_Alc_TopTable$FEATURE.y)

SIMD_Meta_NS_BMI_Alc_TopTable$FEATURE.x[ is.na(SIMD_Meta_NS_BMI_Alc_TopTable$FEATURE.x) ] <- SIMD_Meta_NS_BMI_Alc_TopTable$FEATURE.y[ is.na(SIMD_Meta_NS_BMI_Alc_TopTable$FEATURE.x) ]
SIMD_Meta_NS_BMI_Alc_TopTable$FEATURE.y =NULL

SIMD_Meta_NS_BMI_Alc_TopTable$CpGISLAND.x[ is.na(SIMD_Meta_NS_BMI_Alc_TopTable$CpGISLAND.x) ] <- SIMD_Meta_NS_BMI_Alc_TopTable$CpGISLAND.y[ is.na(SIMD_Meta_NS_BMI_Alc_TopTable$CpGISLAND.x) ]
SIMD_Meta_NS_BMI_Alc_TopTable$CpGISLAND.y =NULL

colnames(SIMD_Meta_NS_BMI_Alc_TopTable)[6] <- "geneSymbol"
colnames(SIMD_Meta_NS_BMI_Alc_TopTable)[7] <- "CHR"
colnames(SIMD_Meta_NS_BMI_Alc_TopTable)[8] <- "MAPINFO"
colnames(SIMD_Meta_NS_BMI_Alc_TopTable)[9] <- "FEATURE"
colnames(SIMD_Meta_NS_BMI_Alc_TopTable)[10] <- "CpGISLAND"

View(SIMD_Meta_NS_BMI_Alc_TopTable)

SIMD_Meta_NS_BMI_Alc_pres_table <- head(SIMD_Meta_NS_BMI_Alc_TopTable[order(SIMD_Meta_NS_BMI_Alc_TopTable$P.value), ], 10)
View(SIMD_Meta_NS_BMI_Alc_pres_table)

write.csv(SIMD_Meta_NS_BMI_Alc_pres_table, file = "SIMD_Meta_NS_BMI_Alc_Top10.csv")

#creating a manhattan plot of the meta analysis
SIMD_Meta_NS_BMI_Alc_Man <- SIMD_Meta_NS_BMI_Alc_TopTable[,c(1,4,6,7,8)]
colnames(SIMD_Meta_NS_BMI_Alc_Man)[1] <- "SNP"
colnames(SIMD_Meta_NS_BMI_Alc_Man)[2] <- "P"
colnames(SIMD_Meta_NS_BMI_Alc_Man)[4] <- "CHR"
colnames(SIMD_Meta_NS_BMI_Alc_Man)[5] <- "BP"


SIMD_Meta_NS_BMI_Alc_Man$CHR <- gsub("chr", "", SIMD_Meta_NS_BMI_Alc_Man$CHR)
SIMD_Meta_NS_BMI_Alc_Man <- transform(SIMD_Meta_NS_BMI_Alc_Man, CHR = as.numeric(CHR))
SIMD_Meta_NS_BMI_Alc_Man <- transform(SIMD_Meta_NS_BMI_Alc_Man, BP = as.numeric(BP))
SIMD_Meta_NS_BMI_Alc_Man <- na.omit(SIMD_Meta_NS_BMI_Alc_Man)

View(SIMD_Meta_NS_BMI_Alc_Man)

manhattan((SIMD_Meta_NS_BMI_Alc_Man), col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-09)), main = "Never Smoke SIMD Meta (BMI, Alcohol)", cex.main = 1.5) #export with width 800 and keep aspect ratio







##################################### W1 and W3 differences #####################################################################################################

setwd("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/")

cellcountinfo <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/cellcountinfo.csv", header=TRUE, sep=",")
IDinfo<- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/IDinfo.csv", header=TRUE, sep=",")
genscot_merged <- read.csv("./genscot_merged.csv")
alcohol <- read.csv("alcohol.csv")
colnames(IDinfo)[1] <-"ID"
SIMD_pheno_W1 <- merge(genscot_merged [,c("ID", "ever_smoke", "pack_years", "area", "rank", "age", "bmi", "dep_status")], IDinfo, by = "ID") 
colnames(pedigree)[2] <- "ID"
SIMD_pheno_W1 <- merge(SIMD_pheno_W1, pedigree[,c("famid", "ID")], by = "ID")

SIMD_pheno_W1[,c(9,10)] = NULL

SIMD_pheno_W1$Wave <- rep("W1",nrow(SIMD_pheno_W1))
SIMD_pheno_W3$Wave <- rep("W3",nrow(SIMD_pheno_W3))
SIMD_pheno <- rbind(SIMD_pheno_W1, SIMD_pheno_W3)

famid_W1 <- as.data.frame(table(SIMD_pheno_W1$famid))
colnames(famid_W1)[1] <- "famid"


install.packages("readxl")
library("readxl")
W3_sample_list <- read_excel("Wave3_sample_list_from_AC.xlsx")
colnames(W3_sample_list)[7] <-"ID"

SIMD_pheno_W3 <- merge(genscot_merged [,c("ID", "ever_smoke", "pack_years", "area", "rank", "age", "bmi", "dep_status")], W3_sample_list[,"ID"], by = "ID") 
SIMD_pheno_W3 <- merge(SIMD_pheno_W3, pedigree[,c("famid", "ID")], by = "ID")


famid_W3 <- as.data.frame(table(SIMD_pheno_W3$famid))
colnames(famid_W3)[1] <- "famid"


SIMD_famid_freq_W1 <- as.data.frame(table(famid_W1$Freq))
SIMD_famid_freq_W3 <- as.data.frame(table(famid_W3$Freq))

SIMD_famid_freq_W1$Wave <- rep("W1",nrow(SIMD_famid_freq_W1))
SIMD_famid_freq_W3$Wave <- rep("W3",nrow(SIMD_famid_freq_W3))
SIMD_famid_freq <- rbind(SIMD_famid_freq_W1, SIMD_famid_freq_W3)


famid_W1$Freq1 <- ifelse(famid_W1$Freq == 1, "", famid_W1$Freq)
famid_W1$Wave <- rep("W1",nrow(famid_W1))
famid_W3$Wave <- rep("W3",nrow(famid_W3))
famid_freq <- rbind(famid_W1, famid_W3)
famid_freq <- famid_freq[!famid_freq$Freq==1,]

ggplot(famid_freq, aes(x = Freq, fill = Wave)) + geom_density(alpha=.3)

ggplot(famid_freq, aes(x=Freq, fill=Wave)) +
  geom_histogram(binwidth=1, alpha=.5, position = "identity") +
  labs(x = "Number of family ID repeats", y = "Frequency", title = "Relatedness")

SIMD_famid_freq <- merge(SIMD_famid_freq_W1, SIMD_famid_freq_W3, by = "Var1", all = TRUE)

qplot(famid_W1$Freq, geom = "histogram")

#Plotting distriution of SIMD ranks
range(genscot_merged$rank)
breaks_rank = seq(28, 6600, by = 100)
rank_cut = cut(genscot_merged$rank, breaks_rank, right=FALSE)
barplot(table(rank_cut), main = "SIMD rank", xlab = "SIMD rank", ylab = "frequency")

#plotting distribution of SIMD rank in W1
range(SIMD_pheno_W1$rank, na.rm =TRUE)
breaks_rank = seq(1, 6700, by = 100)
rank_cut = cut(SIMD_pheno_W1$rank, breaks_rank, right=FALSE)
barplot(table(rank_cut), main = "SIMD rank W1", xlab = "SIMD rank", ylab = "frequency")

#plotting distribution of SIMD rank in W3
range(SIMD_pheno_W3$rank, na.rm =TRUE)
breaks_rank = seq(1, 6700, by = 100)
rank_cut = cut(SIMD_pheno_W3$rank, breaks_rank, right=FALSE)
barplot(table(rank_cut), main = "SIMD rank W3", xlab = "SIMD rank", ylab = "frequency")


library(ggplot2)

#Pack years
tail(sort(SIMD_pheno_W1$pack_years), 10)
table(SIMD_pheno_W1$pack_years)
qplot(SIMD_pheno_W1$pack_years, geom = "histogram")

ggplot(SIMD_pheno, aes(x = pack_years, fill = Wave)) + geom_density(alpha=.3, position = "identity")+
  labs(x = "Pack Years", y = "Frequency", title = "Distribution of Smoking")


tail(sort(SIMD_pheno_W3$pack_years), 10)
table(SIMD_pheno_W3$pack_years)
qplot(SIMD_pheno_W3$pack_years, geom = "histogram")

#BMI
qplot(SIMD_pheno_W1$bmi,
      geom="histogram",
      binwidth = 0.5,
      main = "Histogram for BMI (W1)")

qplot(SIMD_pheno_W3$bmi,
      geom="histogram",
      binwidth = 0.5,
      main = "Histogram for BMI (W1)")

ggplot(SIMD_pheno, aes(x = bmi, fill = Wave)) + geom_density(alpha=.3, position = "identity")+
  labs(x = "BMI", y = "Frequency", title = "Distribution of BMI")

#Area
area_dist_w1 <- as.data.frame(SIMD_pheno_W1$area)
area_dist_w1$Wave <- rep("W1",nrow(area_dist_w1))
area_dist_w3 <- as.data.frame(SIMD_pheno_W3$area)
area_dist_w3$Wave <- rep("W3",nrow(area_dist_w3))
colnames(area_dist_w1)[1] <- "area"
colnames(area_dist_w3)[1] <- "area"
area_dist <- rbind(area_dist_w1,area_dist_w3)

ggplot(area_dist, aes(x=area, fill=Wave)) +
  geom_histogram(binwidth=1, alpha=.5, position = "identity", stat = "count") +
  labs(x = "Area", y = "Frequency", title = "Participation by area")


#SIMD rank
qplot(SIMD_pheno_W1$rank,
      geom="histogram",
      binwidth = 100,
      main = "Histogram for SIMD rank (W1)")

qplot(SIMD_pheno_W3$rank,
      geom="histogram",
      binwidth = 100,
      main = "Histogram for SIMD rank (W3)")

ggplot(SIMD_pheno, aes(x=rank, fill=Wave)) +
  geom_histogram(binwidth=1, alpha=.5, position = "identity") +
  labs(x = "SIMDrank", y = "Frequency", title = "Distribution of SIMD rank")

ggplot(SIMD_pheno, aes(x = rank, fill = Wave)) + geom_density(alpha=.3)+
  labs(x = "SIMD rank", y = "Frequency", title = "Distribution of SIMD rank")+
  scale_fill_manual(values = c("red2", "mediumblue"))


##calculating realtedness (not actually relatedness tho since these could be non-blood relatives)

pedigree <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/pedigree.csv", header=TRUE, sep=",")
View(pedigree)




############################# Model including BMI and Alcohol #########################################################

setwd("C:/Users/ajesp/Documents/PhD/SIMD_EWAS")
genscot_merged <- read.csv("genscot_merged.csv")
cellcountinfo <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/cellcountinfo.csv", header=TRUE, sep=",")
cellcountinfo_w3 <- read.table("C:/Users/ajesp/Documents/PhD/MDD_inflam_EWAS//samplesheet.final.csv", header=TRUE, sep=",")
IDinfo <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/IDinfo.csv", header=TRUE, sep=",")
Alcohol <- read.csv("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/alcohol.csv")
colnames(Alcohol)[1] <- "ID"
Batch_w3 <- read.csv("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/Batch_w3.csv")
colnames(Batch_w3)[1] <- "ID"


##For Wave 1
#Merge smoking, age, sex, BMI, and alcohol
SIMD_cov_bmi_alcohol_w1 <- merge(genscot_merged[,c("ID","bmi")], Alcohol[,c("ID","units")], by = "ID")
colnames(SIMD_cov_bmi_alcohol_w1)[1] <- "id"
SIMD_cov_bmi_alcohol_w1 <- merge(wave1_modelE_cov, SIMD_cov_bmi_alcohol_w1, by = "id")


#save as RDS
saveRDS(SIMD_cov_bmi_alcohol_w1, file = "SIMD_cov_bmi_alcohol_w1.rds")

##For Wave 3 
#Match IDs of cellcountinfo and IDinfo
IDinfo$Sentrix_ID <- sub("^", "X", IDinfo$Sentrix_ID )
IDinfo$Sentrix_Position <- sub("^", "_", IDinfo$Sentrix_Position )
IDinfo$SampleID <- paste0(IDinfo$Sentrix_ID, IDinfo$Sentrix_Position)

#Merge IDinfo and cellountinfo
cellcountinfo <- merge(cellcountinfo, IDinfo, by = "SampleID", all = TRUE)

#Rename SampleID
colnames(cellcountinfo)[8] <- "ID"
colnames(cellcountinfo_w3)[1] <- "ID"

#concatenate cellcount info from wave 1 and wave 3
cellcount_all <- rbind(cellcountinfo[,c(2:8)], cellcountinfo_w3[,c(1,8:13)])

#Merge smoking, age, sex, BMI, cellcount, alcohol and batch number
#SIMD_cov_bmi_alcohol_w3 <- merge(genscot_merged[,c("ID", "bmi")], cellcount_all[,c(1:4,6,7)], by = "ID")
SIMD_cov_bmi_alcohol_w3 <- merge(genscot_merged[,c("ID", "bmi")], Alcohol[,c("ID", "units")], by = "ID")
colnames(SIMD_cov_bmi_alcohol_w3)[1] <- "id"
SIMD_cov_bmi_alcohol_w3 <- merge(SIMD_cov_bmi_alcohol_w3, wave3_modelE_cov, by = "id")
View(SIMD_cov_bmi_alcohol_w3)


#save as RDS
saveRDS(SIMD_cov_bmi_alcohol_w3, file = "SIMD_cov_bmi_alcohol_w3.rds")

### EWAS including BMI and Alcohol
SIMD_BMI_Alc_W1_TopTable <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_BMI_Alc_w1.toptable.txt")
SIMD_BMI_Alc_W3_TopTable <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_BMI_Alc_w3.toptable.txt")
View(SIMD_BMI_Alc_W1_TopTable)
View(SIMD_BMI_Alc_W3_TopTable)


## Manhattan plot of W1_BMI_Alc
SIMD_man_BMI_Alc_W1 <- SIMD_BMI_Alc_W1_TopTable[,c(1,2,3,4,10)]
colnames(SIMD_man_BMI_Alc_W1)[1] <- "SNP"
colnames(SIMD_man_BMI_Alc_W1)[3] <- "CHR"
colnames(SIMD_man_BMI_Alc_W1)[4] <- "BP"
colnames(SIMD_man_BMI_Alc_W1)[5] <- "P"

SIMD_man_BMI_Alc_W1$CHR <- gsub("chr", "", SIMD_man_BMI_Alc_W1$CHR)
SIMD_man_BMI_Alc_W1 <- transform(SIMD_man_BMI_Alc_W1, CHR = as.numeric(CHR))



manhattan(SIMD_man_BMI_Alc_W1, col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-9)), main = "SIMD W1 (BMI, Alcohol, Smoking)", cex.main = 1.5) #export with width 800 and keep aspect ratio


## Creating a presentation table for BMI_Alc_W1
SIMD_BMI_Alc_W1_Top10 <- head(SIMD_BMI_Alc_W1_TopTable[order(SIMD_BMI_Alc_W1_TopTable$P.Value), ], 10)
View(SIMD_BMI_Alc_W1_Top10)

write.csv(SIMD_BMI_Alc_W1_Top10, file = "SIMD_BMI_Alc_W1_Top10.csv")


##Manhattan plot for BMI_Alc w3
SIMD_man_BMI_Alc_W3 <- SIMD_BMI_Alc_W3_TopTable[,c(1,2,3,4,10)]
colnames(SIMD_man_BMI_Alc_W3)[1] <- "SNP"
colnames(SIMD_man_BMI_Alc_W3)[3] <- "CHR"
colnames(SIMD_man_BMI_Alc_W3)[4] <- "BP"
colnames(SIMD_man_BMI_Alc_W3)[5] <- "P"

SIMD_man_BMI_Alc_W3$CHR <- gsub("chr", "", SIMD_man_BMI_Alc_W3$CHR)
SIMD_man_BMI_Alc_W3 <- transform(SIMD_man_BMI_Alc_W3, CHR = as.numeric(CHR))
SIMD_man_BMI_Alc_W3 <- transform(SIMD_man_BMI_Alc_W3, BP = as.numeric(BP))
SIMD_man_BMI_Alc_W3 <- na.omit(SIMD_man_BMI_Alc_W3)



manhattan(SIMD_man_BMI_Alc_W3, col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-9)), main = "SIMD W3 (BMI, Alcohol, Smoking)", cex.main = 1.5) #export with width 800 and keep aspect ratio


## Creating a presentation table for BMI_Alc_W1
SIMD_BMI_Alc_W3_Top10 <- head(SIMD_BMI_Alc_W3_TopTable[order(SIMD_BMI_Alc_W3_TopTable$P.Value), ], 10)
View(SIMD_BMI_Alc_W3_Top10)

write.csv(SIMD_BMI_Alc_W3_Top10, file = "SIMD_BMI_Alc_W3_Top10.csv")

##Meta analysis

#creating wave 1 dataframes for METAL meta analysis
Metal_SIMD_BMI_Alc_W1 <- SIMD_BMI_Alc_W1_TopTable[, c("ID", "P.Value", "N", "beta")]
View(Metal_SIMD_BMI_Alc_W1)

#adding reference allele columns to the metal data frame
Metal_SIMD_BMI_Alc_W1$A1 <- rep("A",nrow(Metal_SIMD_BMI_Alc_W1))
Metal_SIMD_BMI_Alc_W1$A2 <- rep("C",nrow(Metal_SIMD_BMI_Alc_W1))

#adding a direction of effect column
Metal_SIMD_BMI_Alc_W1$Effect <- ifelse(Metal_SIMD_BMI_Alc_W1$beta >= 0, "+", "-")

#creating wave 3 dataframes for METAL meta analysis
Metal_SIMD_BMI_Alc_W3 <- SIMD_BMI_Alc_W3_TopTable[, c("ID", "P.Value", "N", "beta")]
View(Metal_SIMD_BMI_Alc_W3)

#adding reference allele columns to the metal data frame
Metal_SIMD_BMI_Alc_W3$A1 <- rep("A",nrow(Metal_SIMD_BMI_Alc_W3))
Metal_SIMD_BMI_Alc_W3$A2 <- rep("C",nrow(Metal_SIMD_BMI_Alc_W3))

#adding a direction of effect column
Metal_SIMD_BMI_Alc_W3$Effect <- ifelse(Metal_SIMD_BMI_Alc_W3$beta >= 0, "+", "-")

#save dataframes 
write.csv(Metal_SIMD_BMI_Alc_W1, file = "Metal_SIMD_BMI_Alc_W1.csv", row.names = FALSE, quote = FALSE)
write.csv(Metal_SIMD_BMI_Alc_W3, file = "Metal_SIMD_BMI_Alc_W3.csv", row.names = FALSE, quote = FALSE)

#Opening meta analysis
SIMD_BMI_Alc_meta_Output <- read.table("SIMD1BMI_Alc_2", header = TRUE, fill = TRUE)
View(SIMD_BMI_Alc_meta_Output)

#merging in gene annotation, chromosome and map info from both w1 and w3 dataframes
SIMD_Meta_BMI_Alc_TopTable <- merge(SIMD_BMI_Alc_meta_Output[,c("MarkerName", "Weight", "Zscore", "P.value", "Direction")], 
                                    SIMD_BMI_Alc_W1_TopTable[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_BMI_Alc_TopTable <- merge(SIMD_Meta_BMI_Alc_TopTable, SIMD_BMI_Alc_W3_TopTable[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_BMI_Alc_TopTable$geneSymbol.x <- as.character(SIMD_Meta_BMI_Alc_TopTable$geneSymbol.x)
SIMD_Meta_BMI_Alc_TopTable$geneSymbol.y <- as.character(SIMD_Meta_BMI_Alc_TopTable$geneSymbol.y)

SIMD_Meta_BMI_Alc_TopTable$geneSymbol.x[ is.na(SIMD_Meta_BMI_Alc_TopTable$geneSymbol.x) ] <- SIMD_Meta_BMI_Alc_TopTable$geneSymbol.y[ is.na(SIMD_Meta_BMI_Alc_TopTable$geneSymbol.x) ]
SIMD_Meta_BMI_Alc_TopTable$geneSymbol.y =NULL

SIMD_Meta_BMI_Alc_TopTable$CHR.x[ is.na(SIMD_Meta_BMI_Alc_TopTable$CHR.x) ] <- SIMD_Meta_BMI_Alc_TopTable$CHR.y[ is.na(SIMD_Meta_BMI_Alc_TopTable$CHR.x) ]
SIMD_Meta_BMI_Alc_TopTable$CHR.y =NULL

SIMD_Meta_BMI_Alc_TopTable$MAPINFO.x[ is.na(SIMD_Meta_BMI_Alc_TopTable$MAPINFO.x) ] <- SIMD_Meta_BMI_Alc_TopTable$MAPINFO.y[ is.na(SIMD_Meta_BMI_Alc_TopTable$MAPINFO.x) ]
SIMD_Meta_BMI_Alc_TopTable$MAPINFO.y =NULL

SIMD_Meta_BMI_Alc_TopTable$FEATURE.x <- as.character(SIMD_Meta_BMI_Alc_TopTable$FEATURE.x)
SIMD_Meta_BMI_Alc_TopTable$FEATURE.y <- as.character(SIMD_Meta_BMI_Alc_TopTable$FEATURE.y)

SIMD_Meta_BMI_Alc_TopTable$FEATURE.x[ is.na(SIMD_Meta_BMI_Alc_TopTable$FEATURE.x) ] <- SIMD_Meta_BMI_Alc_TopTable$FEATURE.y[ is.na(SIMD_Meta_BMI_Alc_TopTable$FEATURE.x) ]
SIMD_Meta_BMI_Alc_TopTable$FEATURE.y =NULL

SIMD_Meta_BMI_Alc_TopTable$CpGISLAND.x[ is.na(SIMD_Meta_BMI_Alc_TopTable$CpGISLAND.x) ] <- SIMD_Meta_BMI_Alc_TopTable$CpGISLAND.y[ is.na(SIMD_Meta_BMI_Alc_TopTable$CpGISLAND.x) ]
SIMD_Meta_BMI_Alc_TopTable$CpGISLAND.y =NULL

colnames(SIMD_Meta_BMI_Alc_TopTable)[6] <- "geneSymbol"
colnames(SIMD_Meta_BMI_Alc_TopTable)[7] <- "CHR"
colnames(SIMD_Meta_BMI_Alc_TopTable)[8] <- "MAPINFO"
colnames(SIMD_Meta_BMI_Alc_TopTable)[9] <- "FEATURE"
colnames(SIMD_Meta_BMI_Alc_TopTable)[10] <- "CpGISLAND"

View(SIMD_Meta_BMI_Alc_TopTable)

write.csv(SIMD_Meta_BMI_Alc_TopTable, file = "SIMD_Meta_BMI_Alc_TopTable.csv")

SIMD_Meta_BMI_Alc_TopTable <- read.csv("SIMD_Meta_BMI_Alc_TopTable.csv")

SIMD_Meta_BMI_Alc_pres_table <- head(SIMD_Meta_BMI_Alc_TopTable[order(SIMD_Meta_BMI_Alc_TopTable$P.value), ], 10)
View(SIMD_Meta_BMI_Alc_pres_table)

write.csv(SIMD_Meta_BMI_Alc_pres_table, file = "SIMD_BMI_Alc_Meta_Top10.csv")


#creating a manhattan plot of the meta analysis
SIMD_Meta_BMI_Alc_Man <- SIMD_Meta_BMI_Alc_TopTable[,c(1,4,6,7,8)]
colnames(SIMD_Meta_BMI_Alc_Man)[1] <- "SNP"
colnames(SIMD_Meta_BMI_Alc_Man)[2] <- "P"
colnames(SIMD_Meta_BMI_Alc_Man)[4] <- "CHR"
colnames(SIMD_Meta_BMI_Alc_Man)[5] <- "BP"


SIMD_Meta_BMI_Alc_Man$CHR <- gsub("chr", "", SIMD_Meta_BMI_Alc_Man$CHR)
SIMD_Meta_BMI_Alc_Man <- transform(SIMD_Meta_BMI_Alc_Man, CHR = as.numeric(CHR))
SIMD_Meta_BMI_Alc_Man <- transform(SIMD_Meta_BMI_Alc_Man, BP = as.numeric(BP))
SIMD_Meta_BMI_Alc_Man <- na.omit(SIMD_Meta_BMI_Alc_Man)

View(SIMD_Meta_BMI_Alc_Man)

manhattan((SIMD_Meta_BMI_Alc_Man), col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-09)), main = "SIMD Meta (BMI, Alcohol, Smoking)", cex.main = 1.5) #export with width 800 and keep aspect ratio





##################### Correlating effect sizes ###############################################################################

SIMD_correlation_model12 <- merge(SIMD_corrected_meta_Output[,c(1,5)], SIMD_BMI_Alc_meta_Output[,c(1,5)], by = "MarkerName")
cor.test(SIMD_correlation_model12$Zscore.x, SIMD_correlation_model12$Zscore.y, method = "pearson")

SIMD_correlation_model12 <- merge(SIMD_corrected_meta_Output[,c(1,6)], SIMD_BMI_Alc_meta_Output[,c(1,6)], by = "MarkerName")
cor.test(SIMD_correlation_model12$P.value.x, SIMD_correlation_model12$P.value.y, method = "pearson")

#betas

SIMD_correlation_w1 <- merge(SIMD_W1_TopTables[,c(1,14)], SIMD_BMI_Alc_W1_TopTable[,c(1,14)], by = "ID")
cor.test(SIMD_correlation_w1$beta.x, SIMD_correlation_w1$beta.y, method = "pearson")

SIMD_correlation_w3 <- merge(SIMD_corrected_w3_TopTables[,c(1,14)], SIMD_BMI_Alc_W3_TopTable[,c(1,14)], by = "ID")
cor.test(SIMD_correlation_w3$beta.x, SIMD_correlation_w3$beta.y, method = "pearson")

SIMD_correlation_model1 <- merge(SIMD_W1_TopTables[,c(1,14)], SIMD_corrected_w3_TopTables[,c(1,14)], by = "ID")
cor.test(SIMD_correlation_model1$beta.x, SIMD_correlation_model1$beta.y, method = "pearson")

SIMD_correlation_model2 <- merge(SIMD_BMI_Alc_W1_TopTable[,c(1,14)], SIMD_BMI_Alc_W3_TopTable[,c(1,14)], by = "ID")
cor.test(SIMD_correlation_model2$beta.x, SIMD_correlation_model2$beta.y, method = "pearson")


#B
SIMD_correlation_w1 <- merge(SIMD_W1_TopTables[,c(1,12)], SIMD_BMI_Alc_W1_TopTable[,c(1,12)], by = "ID")
cor.test(SIMD_correlation_w1$B.x, SIMD_correlation_w1$B.y, method = "pearson")

SIMD_correlation_w3 <- merge(SIMD_corrected_w3_TopTables[,c(1,12)], SIMD_BMI_Alc_W3_TopTable[,c(1,12)], by = "ID")
cor.test(SIMD_correlation_w3$B.x, SIMD_correlation_w3$B.y, method = "pearson")

SIMD_correlation_model1 <- merge(SIMD_W1_TopTables[,c(1,12)], SIMD_corrected_w3_TopTables[,c(1,12)], by = "ID")
cor.test(SIMD_correlation_model1$B.x, SIMD_correlation_model1$B.y, method = "pearson")

SIMD_correlation_model2 <- merge(SIMD_BMI_Alc_W1_TopTable[,c(1,12)], SIMD_BMI_Alc_W3_TopTable[,c(1,12)], by = "ID")
cor.test(SIMD_correlation_model2$B.x, SIMD_correlation_model2$B.y, method = "pearson")


#AveExpr
SIMD_correlation_w1 <- merge(SIMD_W1_TopTables[,c(1,8)], SIMD_BMI_Alc_W1_TopTable[,c(1,8)], by = "ID")
cor.test(SIMD_correlation_w1$AveExpr.x, SIMD_correlation_w1$AveExpr.y, method = "pearson")

SIMD_correlation_w3 <- merge(SIMD_corrected_w3_TopTables[,c(1,8)], SIMD_BMI_Alc_W3_TopTable[,c(1,8)], by = "ID")
cor.test(SIMD_correlation_w3$AveExpr.x, SIMD_correlation_w3$AveExpr.y, method = "pearson")

SIMD_correlation_model1 <- merge(SIMD_W1_TopTables[,c(1,8)], SIMD_corrected_w3_TopTables[,c(1,8)], by = "ID")
cor.test(SIMD_correlation_model1$AveExpr.x, SIMD_correlation_model1$AveExpr.y, method = "pearson")

SIMD_correlation_model2 <- merge(SIMD_BMI_Alc_W1_TopTable[,c(1,8)], SIMD_BMI_Alc_W3_TopTable[,c(1,8)], by = "ID")
cor.test(SIMD_correlation_model2$AveExpr.x, SIMD_correlation_model2$AveExpr.y, method = "pearson")


#P.value
SIMD_correlation_model2 <- merge(SIMD_BMI_Alc_W1_TopTable[,c(1,10)], SIMD_BMI_Alc_W3_TopTable[,c(1,10)], by = "ID")
cor.test(SIMD_correlation_model2$P.Value.x, SIMD_correlation_model2$P.Value.y, method = "pearson")



##checking covariate files for BMI and Alc model
SIMD_cov_bmi_alcohol_w1 <- readRDS("SIMD_cov_bmi_alcohol_w1.rds")
View(SIMD_cov_bmi_alcohol_w1)

SIMD_cov_bmi_alcohol_w3 <- readRDS("SIMD_cov_bmi_alcohol_w3.rds")
View(SIMD_cov_bmi_alcohol_w3)

SIMD_BMI_Alc_W1_TopTable <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_BMI_Alc_w1.toptable.txt")
SIMD_BMI_Alc_W3_TopTable <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_BMI_Alc_w3.toptable.txt")
View(SIMD_BMI_Alc_W1_TopTable)
View(SIMD_BMI_Alc_W3_TopTable)

## BMI test
setwd("c:/Users/ajesp/Documents/PhD/SIMD_EWAS/BMI test/")
BMI_w1 <- read.delim("./BMI_W1.toptable.txt")
BMI_w3 <- read.delim("./BMI_W3.toptable.txt")
View(BMI_w1)
View(BMI_w3)
BMI_correlation <- merge(BMI_w1[,c(1,14)], BMI_w3[,c(1,14)], by = "ID")
cor.test(BMI_correlation$beta.x, BMI_correlation$beta.y, method = "spearman")

BMI_correlation <- merge(BMI_w1[,c(1,10)], BMI_w3[,c(1,10)], by = "ID")
cor.test(BMI_correlation$P.Value.x, BMI_correlation$P.Value.y, method = "spearman")

##alcohol test
alcohol_w1 <- read.delim("./units_ts.toptable.txt", header = TRUE)
alcohol_w3 <- read.delim("./units_w3_newcov.toptable.txt", header = TRUE)
View(alcohol_w1)
View(alcohol_w3)

alcohol_correlation <- merge(alcohol_w1[,c(1,14)], alcohol_w3[,c(1,14)], by = "ID")
cor.test(alcohol_correlation$beta.x, alcohol_correlation$beta.y, method = "spearman")

alcohol_correlation <- merge(alcohol_w1[,c(1,10)], alcohol_w3[,c(1,10)], by = "ID")
cor.test(alcohol_correlation$P.Value.x, alcohol_correlation$P.Value.y, method = "spearman")









################### Pathway enrichment analysis ###############################################################################################

library(org.Hs.eg.db) 
library(KEGGREST)
library(Category)
library(GO.db)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(stringr)
library(dplyr)
library(AnnotationDbi) 
library(qusage)

# Get annotation table for 850k
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# Extract all gene symbols from the 850k Annotation table
allgenes <- anno$UCSC_RefGene_Name
allgenes <- as.character(allgenes)
allgenes_split <- list()
for(i in 1:length(allgenes)) { 
  allgenes_split[i] <- str_split(allgenes[i], ";")
}
allgenes_split2 <- list()
for(i in 1:length(allgenes_split)) { 
  allgenes_split2[i] <- unique(allgenes_split[[i]][[1]])
}
allgenes <- unique(unlist(allgenes_split2))
allgenes <- allgenes[-which(allgenes == "")]


# Get KEGG Pathways
pathways <- keggList("pathway", "hsa")

# make them into KEGG-style human pathway identifiers
human.pathways <- sub("path:", "", names(pathways))
pathway.ids <- human.pathways

terms.list <- list()
for(i in seq(1, length(pathway.ids), 10)) {
  for(j in i:(i+9)){
    ind <- which(i:(i+9) == j)
  }
  tmp <- pathway.ids[i:(i+9)]
  tmp.pathways <- setNames(keggGet(tmp), tmp)
  for(j in i:(i+9)){
    ind <- which(i:(i+9) == j)
    terms.list[[tmp.pathways[[ind]]$ENTRY]] <-  tmp.pathways[[ind]]$GENE[c(TRUE, FALSE)]
  }
}



# Get EntrezIDs for probes             
entrez <- AnnotationDbi::select(org.Hs.eg.db, allgenes, "ENTREZID", "SYMBOL")
# Gene symbols with multiple ENTREZ IDs
doubles <- subset(entrez, SYMBOL %in% names(which(table(entrez$SYMBOL)>1)))$ENTREZID
which(as.numeric(doubles) %in% as.numeric(unlist(terms.list)))
#[1] 3   # Keep 3, remove 4. Can remove 1 or 2 and 5 or 6 
doubles2 <- subset(entrez, SYMBOL %in% names(which(table(entrez$SYMBOL)>1)))

double_rm <- doubles[c(1,4,5)]
entrez <- entrez[-which(entrez$ENTREZID %in% double_rm)]
all.geneIDs <- na.omit(entrez[,2])

# Limit the KEGG pathway genes to those in the EPIC array universe
for(i in 1:length(terms.list)) {
  terms.list[[i]] <- terms.list[[i]][which(terms.list[[i]] %in% all.geneIDs)]
}

# Run enrichment analysis on SIMD W1 EWAS using suggestive p = 1 x 10-5
significant.genes <- SIMD_W1_TopTables[SIMD_W1_TopTables$P.Value < 3.6e-08, "geneSymbol"] # Character vector of your significant gene IDs from annotation table/toptable
significant.genes_split <- list()
for(i in 1:length(significant.genes)) { 
  significant.genes_split[i] <- str_split(significant.genes[i], ";")
}
significant.genes_split2 <- list()
for(i in 1:length(significant.genes_split)) { 
  significant.genes_split2[i] <- unique(significant.genes_split[[i]][[1]])
}
significant.genes <- unique(unlist(significant.genes_split2))
# significant.genes <- significant.genes[which(significant.genes == "")]
significant.genes <- entrez[which(entrez$SYMBOL %in% significant.genes), "ENTREZID"] 
significant.genes <- significant.genes[which(!is.na(significant.genes))]  # Converts gene symbols to Entrez IDs where possible  



## Enrichment analysis including pathways with at least two genes of interest ##

# Hypergeometric test for enrichment (KEGG terms)
# phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
# q = size of overlap between significant genes and your pathway genes -1;
# m = number of genes in pathway;
# n = (total number of genes on platform - m);
# k = number of significant genes. 

hyper <- list()
for(i in 1:length(terms.list)){
  hyper[[names(terms.list)[i]]] <- 
    ifelse(length(intersect(significant.genes, terms.list[[i]])) > 1,
           (1 - phyper((length(intersect(significant.genes, terms.list[[i]]))-1), 
                       length(terms.list[[i]]),
                       (length(all.geneIDs)-length(significant.genes)),
                       length(significant.genes))), 
           NA)
}

hyperdf <- data.frame(PathID = names(hyper), P = unlist(hyper), Pathway =NA)
for(i in 1:nrow(hyperdf)){
  hyperdf$Pathway[i] <- setNames(keggGet(hyperdf$PathID[i]), hyperdf$PathID[i])[[1]]$PATHWAY
}
hyperdf <- hyperdf[order(hyperdf$P),]
hyperdf$Pbon <- p.adjust(hyperdf$P, "bonferroni")
kegg_hyper <- na.omit(hyperdf)
write.table(kegg_hyper, file="KEGG_pathway_enrichment.xls", sep='\t', quote=F, col.names=NA)




# GO analysis
bp <- read.gmt("C:/Users/ajesp/Documents/PhD/MDD_inflam_EWAS/c5.bp.v6.2.entrez.gmt")
for(i in 1:length(bp)){
  bp[[i]] <- bp[[i]][which(bp[[i]] %in% all.geneIDs)]
}
mf <- read.gmt("C:/Users/ajesp/Documents/PhD/MDD_inflam_EWAS/c5.mf.v6.2.entrez.gmt")
for(i in 1:length(mf)){
  mf[[i]] <- mf[[i]][which(mf[[i]] %in% all.geneIDs)]
}
cc <- read.gmt("C:/Users/ajesp/Documents/PhD/MDD_inflam_EWAS/c5.cc.v6.2.entrez.gmt")
for(i in 1:length(cc)){
  cc[[i]] <- cc[[i]][which(cc[[i]] %in% all.geneIDs)]
}

hyper <- list()
for(i in c("bp","mf","cc")){
  terms.list <- get(i)
  hyper[[i]] <- list()
  for(j in names(terms.list)){
    hyper[[i]][[j]] <- ifelse(length(intersect(significant.genes, terms.list[[j]])) > 1,
                              (1 - phyper((length(intersect(significant.genes, terms.list[[j]]))-1), 
                                          length(terms.list[[j]]),
                                          (length(all.geneIDs)-length(significant.genes)),
                                          length(significant.genes))), 
                              NA)
  }
}



hyperdf_mf <- data.frame(Term = names(hyper[["mf"]]), P = unlist(hyper[["mf"]]), Description =names(hyper[["mf"]]), Ontology = "MF")
hyperdf_cc <- data.frame(Term = names(hyper[["cc"]]), P = unlist(hyper[["cc"]]), Description =names(hyper[["cc"]]), Ontology = "CC")
hyperdf_bp <- data.frame(Term = names(hyper[["bp"]]), P = unlist(hyper[["bp"]]), Description =names(hyper[["bp"]]), Ontology = "BP")
go_hyper <- rbind(hyperdf_mf, hyperdf_cc, hyperdf_bp)
go_hyper <- na.omit(go_hyper)
go_hyper <- go_hyper[order(go_hyper$P),]
go_hyper$Pbon <- p.adjust(go_hyper$P, "bonferroni")
write.table(go_hyper, file="GO_term_enrichment", sep='\t', quote=F, col.names=NA)






# GWAS catalog extract function
# This takes a list of toptable data frames and queries the GWAS catalog 
# for genes associated with traits at a given pval threshold
# NB: depends on dplyr and stringr



gwascat <- function(x, meth_Pcol="adj.P.Val", meth_Pthresh=0.05, gwas_Pthresh=5e-8, study="ewas", tt_gene_col="UCSC_RefGene_Name") { 
  x$GeneSymbol <- x[,tt_gene_col]
  x$Meth_Pcol <- x[,meth_Pcol]
  
  gwas <- read.delim("C:/Users/ajesp/Documents/PhD/MDD_inflam_EWAS/gwas_catalog_v1.0-associations_e95_r2019-03-01.tsv", sep='\t', header=T)
  gwas = gwas[which(gwas$P.VALUE < gwas_Pthresh), ]
  
  # query new gwas catalog
  rows <- list()
  for(row in 1:nrow(gwas)){
    rows[[row]] <- str_split(as.character(gwas[row, "REPORTED.GENE.S."]), ", ") %>% unlist
  }
  
  mygenes <- as.character(x[which(x[,"Meth_Pcol"]<meth_Pthresh), "GeneSymbol"]) %>% unique
  if(length(which(is.na(mygenes))) > 0){
    mygenes <- mygenes[-which(is.na(mygenes))]} else{
      mygenes <- mygenes[-which(mygenes=="")]
    }
  if(length(grep(";", mygenes)) > 0) {
    mygenes <- str_split(mygenes, ";") %>% unlist %>% unique
  }
  
  #   
  gwas_mat <- matrix(ncol=5)
  colnames(gwas_mat) <- c("Gene", "Trait", "SNP",  "P-Value", "Study Link")
  for(gene in mygenes) { 
    tmpgene <- paste0("\\b", gene, "\\b") # This should prevent matching to substrings in GWAS table
    # e.g. DISC1 matching to DISC1FP1
    gwas_pvals <- gwas[grep(tmpgene, rows),"P.VALUE"]
    gwas_traits <- as.character(gwas[grep(tmpgene, rows),"DISEASE.TRAIT"])
    gwas_link <- as.character(gwas[grep(tmpgene, rows),"LINK"])
    gwas_snps <- as.character(gwas[grep(tmpgene, rows),"SNPS"])
    gwas_mat <- rbind(gwas_mat, cbind(rep(gsub("\\\\b", "", tmpgene), length(gwas_pvals)), gwas_traits, gwas_snps, gwas_pvals, gwas_link))
  }
  gwas_mat <- gwas_mat[- 1,]
  
  write.table(gwas_mat, file=paste0("gwas_catalog_hits_", study, ".xls"), sep='\t', quote=F, row.names=F)
}


gwascat(SIMD_W1_TopTables, meth_Pcol="P.Value", meth_Pthresh=1e-5, gwas_Pthresh=5e-8, study="SIMD_EWAS", tt_gene_col="geneSymbol")




library(org.Hs.eg.db) 
library(KEGGREST)
library(Category)
library(GO.db)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(stringr)
library(dplyr)
library(AnnotationDbi) 
library(qusage)

# Get annotation table for 850k
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# Extract all gene symbols from the 850k Annotation table
allgenes <- anno$UCSC_RefGene_Name
allgenes <- as.character(allgenes)
allgenes_split <- list()
for(i in 1:length(allgenes)) { 
  allgenes_split[[i]] <- str_split(allgenes[i], ";")
}
allgenes_split2 <- list()
for(i in 1:length(allgenes_split)) { 
  allgenes_split2[[i]] <- unique(allgenes_split[[i]])
}
allgenes <- unique(unlist(allgenes_split2))
allgenes <- allgenes[-which(allgenes == "")]


# Get KEGG Pathways
pathways <- keggList("pathway", "hsa")

# make them into KEGG-style human pathway identifiers

human.pathways <- sub("path:", "", names(pathways))
pathway.ids <- human.pathways

terms.list <- list()
for(i in seq(1, length(pathway.ids), 10)) {
  for(j in i:(i+9)){
    ind <- which(i:(i+9) == j)
  }
  tmp <- pathway.ids[i:(i+9)]
  tmp.pathways <- setNames(keggGet(tmp), tmp)
  for(j in i:(i+9)){
    ind <- which(i:(i+9) == j)
    terms.list[[tmp.pathways[[ind]]$ENTRY]] <-  tmp.pathways[[ind]]$GENE[c(TRUE, FALSE)]
  }
}



# Get EntrezIDs for probes             
entrez <- AnnotationDbi::select(org.Hs.eg.db, allgenes, "ENTREZID", "SYMBOL")
# Gene symbols with multiple ENTREZ IDs
doubles <- subset(entrez, SYMBOL %in% names(which(table(entrez$SYMBOL)>1)))$ENTREZID
which(as.numeric(doubles) %in% as.numeric(unlist(terms.list)))
#[1] 3   # Keep 3, remove 4. Can remove 1 or 2 and 5 or 6 
doubles2 <- subset(entrez, SYMBOL %in% names(which(table(entrez$SYMBOL)>1)))

double_rm <- doubles[c(1,4,5)]
entrez <- entrez[-which(entrez$ENTREZID %in% double_rm)]
all.geneIDs <- na.omit(entrez[,2])

# Limit the KEGG pathway genes to those in the EPIC array universe
for(i in 1:length(terms.list)) {
  terms.list[[i]] <- terms.list[[i]][which(terms.list[[i]] %in% all.geneIDs)]
}

significant.genes <- TopTables[[2]][TopTables[[2]]$P.Value < 1e-05, "geneSymbol"] # Character vector of your significant gene IDs from annotatuib table/toptable
significant.genes_split <- list()
for(i in 1:length(significant.genes)) { 
  significant.genes_split[i] <- str_split(significant.genes[i], ";")
}
significant.genes_split2 <- list()
for(i in 1:length(significant.genes_split)) { 
  significant.genes_split2[[i]] <- unique(significant.genes_split[[i]])
}
significant.genes <- unique(unlist(significant.genes_split2))
significant.genes <- significant.genes[-which(significant.genes == "")]
significant.genes <- entrez[which(entrez$SYMBOL %in% significant.genes), "ENTREZID"] 
significant.genes <- significant.genes[which(!is.na(significant.genes))]  # Converts gene symbols to Entrez IDs where possible  


# Hypergeometric test for enrichment (KEGG terms)
# phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
# q = size of overlap between significant genes and your pathway genes -1;
# m = number of genes in pathway;
# n = (total number of genes on platform - m);
# k = number of significant genes. 

hyper <- list()
for(i in 1:length(terms.list)){
  hyper[[names(terms.list)[i]]] <- 1 - phyper((length(intersect(significant.genes, terms.list[[i]]))-1), 
                                              length(terms.list[[i]]),
                                              (length(all.geneIDs)-length(significant.genes)),
                                              length(significant.genes))
}

hyperdf <- data.frame(PathID = names(hyper), P = unlist(hyper), Pathway =NA)
for(i in 1:nrow(hyperdf)){
  hyperdf$Pathway[i] <- setNames(keggGet(hyperdf$PathID[i]), hyperdf$PathID[i])[[1]]$PATHWAY
}
kegg_hyper <- hyperdf <- hyperdf[order(hyperdf$P),]




### GO analysis ### 

# Using c5 gene sets. F
bp <- read.gmt("c5.bp.v6.2.entrez.gmt")
for(i in 1:length(bp)){
  bp[[i]] <- bp[[i]][which(bp[[i]] %in% all.geneIDs)]
}
mf <- read.gmt("c5.mf.v6.2.entrez.gmt")
for(i in 1:length(mf)){
  mf[[i]] <- mf[[i]][which(mf[[i]] %in% all.geneIDs)]
}
cc <- read.gmt("c5.cc.v6.2.entrez.gmt")
for(i in 1:length(cc)){
  cc[[i]] <- cc[[i]][which(cc[[i]] %in% all.geneIDs)]
}

hyper <- list()
for(i in c("bp","mf","cc")){
  terms.list <- get(i)
  for(j in names(terms.list)){
    hyper[[i]][[j]] <- 1 - phyper((length(intersect(significant.genes, terms.list[[j]]))-1), 
                                  length(terms.list[[j]]),
                                  (length(all.geneIDs)-length(significant.genes)),
                                  length(significant.genes))
  }
}



hyperdf_mf <- data.frame(Term = names(hyper[["mf"]]), P = unlist(hyper[["mf"]]), Description =names(hyper[["mf"]]), Ontology = "MF")
hyperdf_cc <- data.frame(Term = names(hyper[["cc"]]), P = unlist(hyper[["cc"]]), Description =names(hyper[["cc"]]), Ontology = "CC")
hyperdf_bp <- data.frame(Term = names(hyper[["bp"]]), P = unlist(hyper[["bp"]]), Description =names(hyper[["bp"]]), Ontology = "BP")
go_hyper <- rbind(hyperdf_mf, hyperdf_cc, hyperdf_bp)
go_hyper <- go_hyper[order(go_hyper$P),]
go_hyper$fdr <- p.adjust(go_hyper$P, method="fdr")



###############testing Marks Cov file#############################
SIMD_rank_w3_test <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_rank_w3_test.toptable.txt")
SIMD_rank_w3 <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_rank_w3.toptable.txt")
View(SIMD_rank_w3_test)
View(SIMD_rank_w3)

MarkTest_correlation <- merge(SIMD_rank_w3_test[,c(1,14)], SIMD_rank_w3[,c(1,14)], by = "ID")
cor.test(MarkTest_correlation$beta.x, MarkTest_correlation$beta.y, method = "pearson")





########### Model including Passive Smoking #######################################################
tobacco_v2 <- read.csv("./tobaccov2.csv")
tobacco_v5 <- read.csv("./tobaccov5.csv")
View(tobacco_v2)
View(tobacco_v5)

tobacco_all <- rbind(tobacco_v2[,c(1:2, 11:15)], tobacco_v5[,c(1:2, 10:14)])
View(tobacco_all)

#For wave 1
SIMD_cov_PS_BMI_Alc_w1 <- readRDS("./SIMD_cov_bmi_alcohol_w1.rds")
View(SIMD_cov_PS_BMI_Alc_w1)
SIMD_cov_PS_BMI_Alc_w1 <- merge(SIMD_cov_PS_BMI_Alc_w1, tobacco_all[,c(1,7)], by.x = "id", by.y = "ID")

#save as RDS
saveRDS(SIMD_cov_PS_BMI_Alc_w1, file = "SIMD_cov_PS_BMI_Alc_w1.rds")


#For wave 3
SIMD_cov_PS_BMI_Alc_w3 <- readRDS("./SIMD_cov_bmi_alcohol_w3.rds")
View(SIMD_cov_PS_BMI_Alc_w3)
SIMD_cov_PS_BMI_Alc_w3 <- merge(SIMD_cov_PS_BMI_Alc_w3, tobacco_all[,c(1,7)], by.x = "id", by.y = "ID")
saveRDS(SIMD_cov_PS_BMI_Alc_w3, file = "SIMD_cov_PS_BMI_Alc_w3.rds")

##DO I need this??
setwd("C:/Users/ajesp/Documents/PhD/SIMD_EWAS")
genscot_merged <- read.csv("genscot_merged.csv")
cellcountinfo <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/cellcountinfo.csv", header=TRUE, sep=",")
cellcountinfo_w3 <- read.table("C:/Users/ajesp/Documents/PhD/MDD_inflam_EWAS//samplesheet.final.csv", header=TRUE, sep=",")
IDinfo <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/IDinfo.csv", header=TRUE, sep=",")
Alcohol <- read.csv("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/alcohol.csv")
colnames(Alcohol)[1] <- "ID"
Batch_w3 <- read.csv("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/Batch_w3.csv")
colnames(Batch_w3)[1] <- "ID"


## Opening EWAS output
SIMD_PS_BMI_Alc_w1_TopTable <- read.delim("./SIMD_PS_BMI_Alc_w1.toptable.txt")
SIMD_PS_BMI_Alc_w3_TopTable <- read.delim("./SIMD_PS_BMI_Alc_w3.toptable.txt")
View(SIMD_PS_BMI_Alc_w1_TopTable)
View(SIMD_PS_BMI_Alc_w3_TopTable)

##Meta analysis

#creating wave 1 dataframes for METAL meta analysis
Metal_SIMD_PS_BMI_Alc_W1 <- SIMD_PS_BMI_Alc_w1_TopTable[, c("ID", "P.Value", "N", "beta")]
View(Metal_SIMD_PS_BMI_Alc_W1)

#adding reference allele columns to the metal data frame
Metal_SIMD_PS_BMI_Alc_W1$A1 <- rep("A",nrow(Metal_SIMD_PS_BMI_Alc_W1))
Metal_SIMD_PS_BMI_Alc_W1$A2 <- rep("C",nrow(Metal_SIMD_PS_BMI_Alc_W1))

#adding a direction of effect column
Metal_SIMD_PS_BMI_Alc_W1$Effect <- ifelse(Metal_SIMD_PS_BMI_Alc_W1$beta >= 0, "+", "-")

#creating wave 3 dataframes for METAL meta analysis
Metal_SIMD_PS_BMI_Alc_W3 <- SIMD_PS_BMI_Alc_w3_TopTable[, c("ID", "P.Value", "N", "beta")]
View(Metal_SIMD_PS_BMI_Alc_W3)

#adding reference allele columns to the metal data frame
Metal_SIMD_PS_BMI_Alc_W3$A1 <- rep("A",nrow(Metal_SIMD_PS_BMI_Alc_W3))
Metal_SIMD_PS_BMI_Alc_W3$A2 <- rep("C",nrow(Metal_SIMD_PS_BMI_Alc_W3))

#adding a direction of effect column
Metal_SIMD_PS_BMI_Alc_W3$Effect <- ifelse(Metal_SIMD_PS_BMI_Alc_W3$beta >= 0, "+", "-")

#save dataframes 
write.csv(Metal_SIMD_PS_BMI_Alc_W1, file = "Metal_SIMD_PS_BMI_Alc_W1.csv", row.names = FALSE, quote = FALSE)
write.csv(Metal_SIMD_PS_BMI_Alc_W3, file = "Metal_SIMD_PS_BMI_Alc_W3.csv", row.names = FALSE, quote = FALSE)

Meta_SIMD_PS_BMI_Alc1.tbl

#Opening meta analysis
SIMD_PS_BMI_Alc_meta_Output <- read.table("./Meta_SIMD_PS_BMI_Alc1.tbl", header = TRUE, fill = TRUE)
View(SIMD_PS_BMI_Alc_meta_Output)

#merging in gene annotation, chromosome and map info from both w1 and w3 dataframes
SIMD_Meta_PS_BMI_Alc_TopTable <- merge(SIMD_PS_BMI_Alc_meta_Output[,c("MarkerName", "Weight", "Zscore", "P.value", "Direction")], 
                                    SIMD_PS_BMI_Alc_w1_TopTable[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_PS_BMI_Alc_TopTable <- merge(SIMD_Meta_PS_BMI_Alc_TopTable, SIMD_PS_BMI_Alc_w3_TopTable[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_PS_BMI_Alc_TopTable$geneSymbol.x <- as.character(SIMD_Meta_PS_BMI_Alc_TopTable$geneSymbol.x)
SIMD_Meta_PS_BMI_Alc_TopTable$geneSymbol.y <- as.character(SIMD_Meta_PS_BMI_Alc_TopTable$geneSymbol.y)

SIMD_Meta_PS_BMI_Alc_TopTable$geneSymbol.x[ is.na(SIMD_Meta_PS_BMI_Alc_TopTable$geneSymbol.x) ] <- SIMD_Meta_PS_BMI_Alc_TopTable$geneSymbol.y[ is.na(SIMD_Meta_PS_BMI_Alc_TopTable$geneSymbol.x) ]
SIMD_Meta_PS_BMI_Alc_TopTable$geneSymbol.y =NULL

SIMD_Meta_PS_BMI_Alc_TopTable$CHR.x[ is.na(SIMD_Meta_PS_BMI_Alc_TopTable$CHR.x) ] <- SIMD_Meta_PS_BMI_Alc_TopTable$CHR.y[ is.na(SIMD_Meta_PS_BMI_Alc_TopTable$CHR.x) ]
SIMD_Meta_PS_BMI_Alc_TopTable$CHR.y =NULL

SIMD_Meta_PS_BMI_Alc_TopTable$MAPINFO.x[ is.na(SIMD_Meta_PS_BMI_Alc_TopTable$MAPINFO.x) ] <- SIMD_Meta_PS_BMI_Alc_TopTable$MAPINFO.y[ is.na(SIMD_Meta_PS_BMI_Alc_TopTable$MAPINFO.x) ]
SIMD_Meta_PS_BMI_Alc_TopTable$MAPINFO.y =NULL

SIMD_Meta_PS_BMI_Alc_TopTable$FEATURE.x <- as.character(SIMD_Meta_PS_BMI_Alc_TopTable$FEATURE.x)
SIMD_Meta_PS_BMI_Alc_TopTable$FEATURE.y <- as.character(SIMD_Meta_PS_BMI_Alc_TopTable$FEATURE.y)

SIMD_Meta_PS_BMI_Alc_TopTable$FEATURE.x[ is.na(SIMD_Meta_PS_BMI_Alc_TopTable$FEATURE.x) ] <- SIMD_Meta_PS_BMI_Alc_TopTable$FEATURE.y[ is.na(SIMD_Meta_PS_BMI_Alc_TopTable$FEATURE.x) ]
SIMD_Meta_PS_BMI_Alc_TopTable$FEATURE.y =NULL

SIMD_Meta_PS_BMI_Alc_TopTable$CpGISLAND.x[ is.na(SIMD_Meta_PS_BMI_Alc_TopTable$CpGISLAND.x) ] <- SIMD_Meta_PS_BMI_Alc_TopTable$CpGISLAND.y[ is.na(SIMD_Meta_PS_BMI_Alc_TopTable$CpGISLAND.x) ]
SIMD_Meta_PS_BMI_Alc_TopTable$CpGISLAND.y =NULL

colnames(SIMD_Meta_PS_BMI_Alc_TopTable)[6] <- "geneSymbol"
colnames(SIMD_Meta_PS_BMI_Alc_TopTable)[7] <- "CHR"
colnames(SIMD_Meta_PS_BMI_Alc_TopTable)[8] <- "MAPINFO"
colnames(SIMD_Meta_PS_BMI_Alc_TopTable)[9] <- "FEATURE"
colnames(SIMD_Meta_PS_BMI_Alc_TopTable)[10] <- "CpGISLAND"

View(SIMD_Meta_PS_BMI_Alc_TopTable)

write.csv(SIMD_Meta_PS_BMI_Alc_TopTable, file = "SIMD_Meta_PS_BMI_Alc_TopTable.csv")

SIMD_Meta_PS_BMI_Alc_TopTable <- read.csv("SIMD_Meta_PS_BMI_Alc_TopTable.csv")

SIMD_Meta_PS_BMI_Alc_pres_table <- head(SIMD_Meta_PS_BMI_Alc_TopTable[order(SIMD_Meta_PS_BMI_Alc_TopTable$P.value), ], 10)
View(SIMD_Meta_PS_BMI_Alc_pres_table)

write.csv(SIMD_Meta_PS_BMI_Alc_pres_table, file = "SIMD_PS_BMI_Alc_Meta_Top10.csv")


#creating a manhattan plot of the meta analysis
SIMD_Meta_PS_BMI_Alc_Man <- SIMD_Meta_PS_BMI_Alc_TopTable[,c(1,4,6,7,8)]
colnames(SIMD_Meta_PS_BMI_Alc_Man)[1] <- "SNP"
colnames(SIMD_Meta_PS_BMI_Alc_Man)[2] <- "P"
colnames(SIMD_Meta_PS_BMI_Alc_Man)[4] <- "CHR"
colnames(SIMD_Meta_PS_BMI_Alc_Man)[5] <- "BP"


SIMD_Meta_PS_BMI_Alc_Man$CHR <- gsub("chr", "", SIMD_Meta_PS_BMI_Alc_Man$CHR)
SIMD_Meta_PS_BMI_Alc_Man <- transform(SIMD_Meta_PS_BMI_Alc_Man, CHR = as.numeric(CHR))
SIMD_Meta_PS_BMI_Alc_Man <- transform(SIMD_Meta_PS_BMI_Alc_Man, BP = as.numeric(BP))
SIMD_Meta_PS_BMI_Alc_Man <- na.omit(SIMD_Meta_PS_BMI_Alc_Man)

View(SIMD_Meta_PS_BMI_Alc_Man)

manhattan((SIMD_Meta_PS_BMI_Alc_Man), col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-09)), main = "SIMD Meta (BMI, Alcohol, Smoking, Passive Smoking)", cex.main = 1.5) #export with width 800 and keep aspect ratio




############ Model including urban vs rural living (UvR) #################################################

#Opening the UvR data file
library(readxl)
SIMD_urban_rural <- read_excel("~/PhD/SIMD_EWAS/SIMD urban rural.xlsx")

#For wave 1
SIMD_cov_UvR_BMI_Alc_w1 <- readRDS("./SIMD_cov_bmi_alcohol_w1.rds")
View(SIMD_cov_UvR_BMI_Alc_w1)
SIMD_cov_UvR_BMI_Alc_w1 <- merge(SIMD_cov_UvR_BMI_Alc_w1, SIMD_urban_rural[,c(1,4)], by = "id")

#save as RDS
saveRDS(SIMD_cov_UvR_BMI_Alc_w1, file = "SIMD_cov_UvR_BMI_Alc_w1.rds")


#For wave 3
SIMD_cov_UvR_BMI_Alc_w3 <- readRDS("./SIMD_cov_bmi_alcohol_w3.rds")
View(SIMD_cov_UvR_BMI_Alc_w3)
SIMD_cov_UvR_BMI_Alc_w3 <- merge(SIMD_cov_UvR_BMI_Alc_w3, SIMD_urban_rural[,c(1,4)], by = "id")
saveRDS(SIMD_cov_UvR_BMI_Alc_w3, file = "SIMD_cov_UvR_BMI_Alc_w3.rds")


## Opening EWAS output
SIMD_UvR_BMI_Alc_w1_TopTable <- read.delim("./SIMD_UvR_BMI_Alc_w1.toptable.txt")
SIMD_UvR_BMI_Alc_w3_TopTable <- read.delim("./SIMD_UvR_BMI_Alc_w3.toptable.txt")
View(SIMD_UvR_BMI_Alc_w1_TopTable)
View(SIMD_UvR_BMI_Alc_w3_TopTable)

##Meta analysis

#creating wave 1 dataframes for METAL meta analysis
Metal_SIMD_UvR_BMI_Alc_W1 <- SIMD_UvR_BMI_Alc_w1_TopTable[, c("ID", "P.Value", "N", "beta")]
View(Metal_SIMD_UvR_BMI_Alc_W1)

#adding reference allele columns to the metal data frame
Metal_SIMD_UvR_BMI_Alc_W1$A1 <- rep("A",nrow(Metal_SIMD_UvR_BMI_Alc_W1))
Metal_SIMD_UvR_BMI_Alc_W1$A2 <- rep("C",nrow(Metal_SIMD_UvR_BMI_Alc_W1))

#adding a direction of effect column
Metal_SIMD_UvR_BMI_Alc_W1$Effect <- ifelse(Metal_SIMD_UvR_BMI_Alc_W1$beta >= 0, "+", "-")

#creating wave 3 dataframes for METAL meta analysis
Metal_SIMD_UvR_BMI_Alc_W3 <- SIMD_UvR_BMI_Alc_w3_TopTable[, c("ID", "P.Value", "N", "beta")]
View(Metal_SIMD_UvR_BMI_Alc_W3)

#adding reference allele columns to the metal data frame
Metal_SIMD_UvR_BMI_Alc_W3$A1 <- rep("A",nrow(Metal_SIMD_UvR_BMI_Alc_W3))
Metal_SIMD_UvR_BMI_Alc_W3$A2 <- rep("C",nrow(Metal_SIMD_UvR_BMI_Alc_W3))

#adding a direction of effect column
Metal_SIMD_UvR_BMI_Alc_W3$Effect <- ifelse(Metal_SIMD_UvR_BMI_Alc_W3$beta >= 0, "+", "-")

#save dataframes 
write.csv(Metal_SIMD_UvR_BMI_Alc_W1, file = "Metal_SIMD_UvR_BMI_Alc_W1.csv", row.names = FALSE, quote = FALSE)
write.csv(Metal_SIMD_UvR_BMI_Alc_W3, file = "Metal_SIMD_UvR_BMI_Alc_W3.csv", row.names = FALSE, quote = FALSE)


#Open meta analysis
SIMD_UvR_BMI_Alc_meta_Output <- read.table("./Meta_SIMD_UvR_BMI_Alc1.tbl", header = TRUE, fill = TRUE)
View(SIMD_UvR_BMI_Alc_meta_Output)

#merging in gene annotation, chromosome and map info from both w1 and w3 dataframes
SIMD_Meta_UvR_BMI_Alc_TopTable <- merge(SIMD_UvR_BMI_Alc_meta_Output[,c("MarkerName", "Weight", "Zscore", "P.value", "Direction")], 
                                       SIMD_UvR_BMI_Alc_w1_TopTable[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_UvR_BMI_Alc_TopTable <- merge(SIMD_Meta_UvR_BMI_Alc_TopTable, SIMD_UvR_BMI_Alc_w3_TopTable[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_UvR_BMI_Alc_TopTable$geneSymbol.x <- as.character(SIMD_Meta_UvR_BMI_Alc_TopTable$geneSymbol.x)
SIMD_Meta_UvR_BMI_Alc_TopTable$geneSymbol.y <- as.character(SIMD_Meta_UvR_BMI_Alc_TopTable$geneSymbol.y)

SIMD_Meta_UvR_BMI_Alc_TopTable$geneSymbol.x[ is.na(SIMD_Meta_UvR_BMI_Alc_TopTable$geneSymbol.x) ] <- SIMD_Meta_UvR_BMI_Alc_TopTable$geneSymbol.y[ is.na(SIMD_Meta_UvR_BMI_Alc_TopTable$geneSymbol.x) ]
SIMD_Meta_UvR_BMI_Alc_TopTable$geneSymbol.y =NULL

SIMD_Meta_UvR_BMI_Alc_TopTable$CHR.x[ is.na(SIMD_Meta_UvR_BMI_Alc_TopTable$CHR.x) ] <- SIMD_Meta_UvR_BMI_Alc_TopTable$CHR.y[ is.na(SIMD_Meta_UvR_BMI_Alc_TopTable$CHR.x) ]
SIMD_Meta_UvR_BMI_Alc_TopTable$CHR.y =NULL

SIMD_Meta_UvR_BMI_Alc_TopTable$MAPINFO.x[ is.na(SIMD_Meta_UvR_BMI_Alc_TopTable$MAPINFO.x) ] <- SIMD_Meta_UvR_BMI_Alc_TopTable$MAPINFO.y[ is.na(SIMD_Meta_UvR_BMI_Alc_TopTable$MAPINFO.x) ]
SIMD_Meta_UvR_BMI_Alc_TopTable$MAPINFO.y =NULL

SIMD_Meta_UvR_BMI_Alc_TopTable$FEATURE.x <- as.character(SIMD_Meta_UvR_BMI_Alc_TopTable$FEATURE.x)
SIMD_Meta_UvR_BMI_Alc_TopTable$FEATURE.y <- as.character(SIMD_Meta_UvR_BMI_Alc_TopTable$FEATURE.y)

SIMD_Meta_UvR_BMI_Alc_TopTable$FEATURE.x[ is.na(SIMD_Meta_UvR_BMI_Alc_TopTable$FEATURE.x) ] <- SIMD_Meta_UvR_BMI_Alc_TopTable$FEATURE.y[ is.na(SIMD_Meta_UvR_BMI_Alc_TopTable$FEATURE.x) ]
SIMD_Meta_UvR_BMI_Alc_TopTable$FEATURE.y =NULL

SIMD_Meta_UvR_BMI_Alc_TopTable$CpGISLAND.x[ is.na(SIMD_Meta_UvR_BMI_Alc_TopTable$CpGISLAND.x) ] <- SIMD_Meta_UvR_BMI_Alc_TopTable$CpGISLAND.y[ is.na(SIMD_Meta_UvR_BMI_Alc_TopTable$CpGISLAND.x) ]
SIMD_Meta_UvR_BMI_Alc_TopTable$CpGISLAND.y =NULL

colnames(SIMD_Meta_UvR_BMI_Alc_TopTable)[6] <- "geneSymbol"
colnames(SIMD_Meta_UvR_BMI_Alc_TopTable)[7] <- "CHR"
colnames(SIMD_Meta_UvR_BMI_Alc_TopTable)[8] <- "MAPINFO"
colnames(SIMD_Meta_UvR_BMI_Alc_TopTable)[9] <- "FEATURE"
colnames(SIMD_Meta_UvR_BMI_Alc_TopTable)[10] <- "CpGISLAND"

View(SIMD_Meta_UvR_BMI_Alc_TopTable)

write.csv(SIMD_Meta_UvR_BMI_Alc_TopTable, file = "SIMD_Meta_UvR_BMI_Alc_TopTable.csv")

SIMD_Meta_UvR_BMI_Alc_TopTable <- read.csv("SIMD_Meta_UvR_BMI_Alc_TopTable.csv")

SIMD_Meta_UvR_BMI_Alc_pres_table <- head(SIMD_Meta_UvR_BMI_Alc_TopTable[order(SIMD_Meta_UvR_BMI_Alc_TopTable$P.value), ], 10)
View(SIMD_Meta_UvR_BMI_Alc_pres_table)

write.csv(SIMD_Meta_UvR_BMI_Alc_pres_table, file = "SIMD_UvR_BMI_Alc_Meta_Top10.csv")


#creating a manhattan plot of the meta analysis
SIMD_Meta_UvR_BMI_Alc_Man <- SIMD_Meta_UvR_BMI_Alc_TopTable[,c(1,4,6,7,8)]
colnames(SIMD_Meta_UvR_BMI_Alc_Man)[1] <- "SNP"
colnames(SIMD_Meta_UvR_BMI_Alc_Man)[2] <- "P"
colnames(SIMD_Meta_UvR_BMI_Alc_Man)[4] <- "CHR"
colnames(SIMD_Meta_UvR_BMI_Alc_Man)[5] <- "BP"


SIMD_Meta_UvR_BMI_Alc_Man$CHR <- gsub("chr", "", SIMD_Meta_UvR_BMI_Alc_Man$CHR)
SIMD_Meta_UvR_BMI_Alc_Man <- transform(SIMD_Meta_UvR_BMI_Alc_Man, CHR = as.numeric(CHR))
SIMD_Meta_UvR_BMI_Alc_Man <- transform(SIMD_Meta_UvR_BMI_Alc_Man, BP = as.numeric(BP))
SIMD_Meta_UvR_BMI_Alc_Man <- na.omit(SIMD_Meta_UvR_BMI_Alc_Man)

View(SIMD_Meta_UvR_BMI_Alc_Man)

manhattan((SIMD_Meta_UvR_BMI_Alc_Man), col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-10)), main = "SIMD Meta (BMI, Alcohol, Smoking, Urban vs Rural)", cex.main = 1.5) #export with width 800 and keep aspect ratio



######################### Model including AHRR probe as covar. ##########################################
SIMD_Top10_CpGS_W1 <- readRDS("./SIMD_Top10_CpGs_W1.rds")
SIMD_Top10_CpGS_W3 <- readRDS("./SIMD_Top10_CpGs_W3.rds")
wave3_sentrix <- readRDS("./wave3_sentrix.rds")
wave1_sentrix <- readRDS("./sentrix_5101.rds")
SIMD_cov_bmi_alcohol_w1 <- readRDS("./SIMD_cov_bmi_alcohol_w1.rds")
SIMD_cov_bmi_alcohol_w3 <- readRDS("./SIMD_cov_bmi_alcohol_w3.rds")


View(SIMD_Top10_CpGS_W1)

#transpose data so CpFs become columns
SIMD_Top10_CpGS_W1 <- t(SIMD_Top10_CpGS_W1)

#make rownames into a column
SIMD_Top10_CpGS_W1 <- cbind(rownames(SIMD_Top10_CpGS_W1), data.frame(SIMD_Top10_CpGS_W1, row.names=NULL))

#rename ID column
colnames(SIMD_Top10_CpGS_W1)[1] <- "Sample_Sentrix_ID"

#merge with sentrix IDs
SIMD_Top10_CpGS_W1 <- merge(SIMD_Top10_CpGS_W1, wave1_sentrix[,c(1,4)], by = "Sample_Sentrix_ID")


#delete Sample_Sentrix_ID column
SIMD_Top10_CpGS_W1[,"Sample_Sentrix_ID"] <- NULL

#merge into the existing BMI, Alc covar file
SIMD_cov_BMI_Alc_AHRR_w1 <- merge(SIMD_cov_bmi_alcohol_w1, SIMD_Top10_CpGS_W1[,c(1, 11)], by = "id")
saveRDS(SIMD_cov_BMI_Alc_AHRR_w1, file = "SIMD_cov_BMI_Alc_AHRR_w1.rds")

#transpose data so CpFs become columns
SIMD_Top10_CpGS_W3 <- t(SIMD_Top10_CpGS_W3)

#make rownames into a column
SIMD_Top10_CpGS_W3 <- cbind(rownames(SIMD_Top10_CpGS_W3), data.frame(SIMD_Top10_CpGS_W3, row.names=NULL))

#rename ID column
colnames(SIMD_Top10_CpGS_W3)[1] <- "Sample_Sentrix_ID"

#merge with sentrix IDs
SIMD_Top10_CpGS_W3 <- merge(SIMD_Top10_CpGS_W3, wave3_sentrix[,c(1,4)], by = "Sample_Sentrix_ID")


#delete Sample_Sentrix_ID column
SIMD_Top10_CpGS_W3[,"Sample_Sentrix_ID"] <- NULL

#save dataframe
write.csv(SIMD_Top10_CpGS_W1, file = "SIMD_Top10_CpGS_W1.csv")
write.csv(SIMD_Top10_CpGS_W3, file = "SIMD_Top10_CpGS_W3.csv")

#merge into the existing BMI, Alc covar file
SIMD_cov_BMI_Alc_AHRR_w3 <- merge(SIMD_cov_bmi_alcohol_w3, SIMD_Top10_CpGS_W3[,c(1, 11)], by = "id")
saveRDS(SIMD_cov_BMI_Alc_AHRR_w3, file = "SIMD_cov_BMI_Alc_AHRR_w3.rds")


#Opening w1 EWAS output
SIMD_BMI_Alc_AHRR_W1_TopTable <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_BMI_Alc_AHRR_w1.toptable.txt")
View(SIMD_BMI_Alc_AHRR_W1_TopTable)



## Manhattan plot of W1_BMI_Alc
SIMD_man_BMI_Alc_AHRR_W1 <- SIMD_BMI_Alc_AHRR_W1_TopTable[,c(1,2,3,4,10)]
colnames(SIMD_man_BMI_Alc_AHRR_W1)[1] <- "SNP"
colnames(SIMD_man_BMI_Alc_AHRR_W1)[3] <- "CHR"
colnames(SIMD_man_BMI_Alc_AHRR_W1)[4] <- "BP"
colnames(SIMD_man_BMI_Alc_AHRR_W1)[5] <- "P"

SIMD_man_BMI_Alc_AHRR_W1$CHR <- gsub("chr", "", SIMD_man_BMI_Alc_AHRR_W1$CHR)
SIMD_man_BMI_Alc_AHRR_W1 <- transform(SIMD_man_BMI_Alc_AHRR_W1, CHR = as.numeric(CHR))



manhattan(SIMD_man_BMI_Alc_AHRR_W1, col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-9)), main = "SIMD W1 (AHRR, BMI, Alcohol, Smoking)", cex.main = 1.5) #export with width 800 and keep aspect ratio


## Creating a presentation table for BMI_Alc_W1
SIMD_BMI_Alc_AHRR_W1_Top10 <- head(SIMD_BMI_Alc_AHRR_W1_TopTable[order(SIMD_BMI_Alc_AHRR_W1_TopTable$P.Value), ], 10)
View(SIMD_BMI_Alc_AHRR_W1_Top10)

write.csv(SIMD_BMI_Alc_AHRR_W1_Top10, file = "SIMD_BMI_Alc_AHRR_W1_Top10.csv")




#Opening w3 EWAS output
SIMD_BMI_Alc_AHRR_W3_TopTable <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_BMI_Alc_AHRR_w3.toptable.txt")
View(SIMD_BMI_Alc_AHRR_W3_TopTable)



#creating wave 1 dataframes for METAL meta analysis
Metal_SIMD_BMI_Alc_AHRR_W1 <- SIMD_BMI_Alc_AHRR_W1_TopTable[, c("ID", "P.Value", "N", "beta")]
View(Metal_SIMD_BMI_Alc_AHRR_W1)

#adding reference allele columns to the metal data frame
Metal_SIMD_BMI_Alc_AHRR_W1$A1 <- rep("A",nrow(Metal_SIMD_BMI_Alc_AHRR_W1))
Metal_SIMD_BMI_Alc_AHRR_W1$A2 <- rep("C",nrow(Metal_SIMD_BMI_Alc_AHRR_W1))

#adding a direction of effect column
Metal_SIMD_BMI_Alc_AHRR_W1$Effect <- ifelse(Metal_SIMD_BMI_Alc_AHRR_W1$beta >= 0, "+", "-")

#creating wave 3 dataframes for METAL meta analysis
Metal_SIMD_BMI_Alc_AHRR_W3 <- SIMD_BMI_Alc_AHRR_W3_TopTable[, c("ID", "P.Value", "N", "beta")]
View(Metal_SIMD_BMI_Alc_AHRR_W3)

#adding reference allele columns to the metal data frame
Metal_SIMD_BMI_Alc_AHRR_W3$A1 <- rep("A",nrow(Metal_SIMD_BMI_Alc_AHRR_W3))
Metal_SIMD_BMI_Alc_AHRR_W3$A2 <- rep("C",nrow(Metal_SIMD_BMI_Alc_AHRR_W3))

#adding a direction of effect column
Metal_SIMD_BMI_Alc_AHRR_W3$Effect <- ifelse(Metal_SIMD_BMI_Alc_AHRR_W3$beta >= 0, "+", "-")

#save dataframes 
write.csv(Metal_SIMD_BMI_Alc_AHRR_W1, file = "Metal_SIMD_BMI_Alc_AHRR_W1.csv", row.names = FALSE, quote = FALSE)
write.csv(Metal_SIMD_BMI_Alc_AHRR_W3, file = "Metal_SIMD_BMI_Alc_AHRR_W3.csv", row.names = FALSE, quote = FALSE)

#Opening meta analysis
SIMD_BMI_Alc_AHRR_meta_Output <- read.table("Meta_SIMD_BMI_Alc_AHRR1.tbl", header = TRUE, fill = TRUE)
View(SIMD_BMI_Alc_AHRR_meta_Output)

#merging in gene annotation, chromosome and map info from both w1 and w3 dataframes
SIMD_Meta_BMI_Alc_AHRR_TopTable <- merge(SIMD_BMI_Alc_AHRR_meta_Output[,c("MarkerName", "Weight", "Zscore", "P.value", "Direction")], 
                                    SIMD_BMI_Alc_AHRR_W1_TopTable[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_BMI_Alc_AHRR_TopTable <- merge(SIMD_Meta_BMI_Alc_AHRR_TopTable, SIMD_BMI_Alc_AHRR_W3_TopTable[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_BMI_Alc_AHRR_TopTable$geneSymbol.x <- as.character(SIMD_Meta_BMI_Alc_AHRR_TopTable$geneSymbol.x)
SIMD_Meta_BMI_Alc_AHRR_TopTable$geneSymbol.y <- as.character(SIMD_Meta_BMI_Alc_AHRR_TopTable$geneSymbol.y)

SIMD_Meta_BMI_Alc_AHRR_TopTable$geneSymbol.x[ is.na(SIMD_Meta_BMI_Alc_AHRR_TopTable$geneSymbol.x) ] <- SIMD_Meta_BMI_Alc_AHRR_TopTable$geneSymbol.y[ is.na(SIMD_Meta_BMI_Alc_AHRR_TopTable$geneSymbol.x) ]
SIMD_Meta_BMI_Alc_AHRR_TopTable$geneSymbol.y =NULL

SIMD_Meta_BMI_Alc_AHRR_TopTable$CHR.x[ is.na(SIMD_Meta_BMI_Alc_AHRR_TopTable$CHR.x) ] <- SIMD_Meta_BMI_Alc_AHRR_TopTable$CHR.y[ is.na(SIMD_Meta_BMI_Alc_AHRR_TopTable$CHR.x) ]
SIMD_Meta_BMI_Alc_AHRR_TopTable$CHR.y =NULL

SIMD_Meta_BMI_Alc_AHRR_TopTable$MAPINFO.x[ is.na(SIMD_Meta_BMI_Alc_AHRR_TopTable$MAPINFO.x) ] <- SIMD_Meta_BMI_Alc_AHRR_TopTable$MAPINFO.y[ is.na(SIMD_Meta_BMI_Alc_AHRR_TopTable$MAPINFO.x) ]
SIMD_Meta_BMI_Alc_AHRR_TopTable$MAPINFO.y =NULL

SIMD_Meta_BMI_Alc_AHRR_TopTable$FEATURE.x <- as.character(SIMD_Meta_BMI_Alc_AHRR_TopTable$FEATURE.x)
SIMD_Meta_BMI_Alc_AHRR_TopTable$FEATURE.y <- as.character(SIMD_Meta_BMI_Alc_AHRR_TopTable$FEATURE.y)

SIMD_Meta_BMI_Alc_AHRR_TopTable$FEATURE.x[ is.na(SIMD_Meta_BMI_Alc_AHRR_TopTable$FEATURE.x) ] <- SIMD_Meta_BMI_Alc_AHRR_TopTable$FEATURE.y[ is.na(SIMD_Meta_BMI_Alc_AHRR_TopTable$FEATURE.x) ]
SIMD_Meta_BMI_Alc_AHRR_TopTable$FEATURE.y =NULL

SIMD_Meta_BMI_Alc_AHRR_TopTable$CpGISLAND.x[ is.na(SIMD_Meta_BMI_Alc_AHRR_TopTable$CpGISLAND.x) ] <- SIMD_Meta_BMI_Alc_AHRR_TopTable$CpGISLAND.y[ is.na(SIMD_Meta_BMI_Alc_AHRR_TopTable$CpGISLAND.x) ]
SIMD_Meta_BMI_Alc_AHRR_TopTable$CpGISLAND.y =NULL

colnames(SIMD_Meta_BMI_Alc_AHRR_TopTable)[6] <- "geneSymbol"
colnames(SIMD_Meta_BMI_Alc_AHRR_TopTable)[7] <- "CHR"
colnames(SIMD_Meta_BMI_Alc_AHRR_TopTable)[8] <- "MAPINFO"
colnames(SIMD_Meta_BMI_Alc_AHRR_TopTable)[9] <- "FEATURE"
colnames(SIMD_Meta_BMI_Alc_AHRR_TopTable)[10] <- "CpGISLAND"

View(SIMD_Meta_BMI_Alc_AHRR_TopTable)

write.csv(SIMD_Meta_BMI_Alc_AHRR_TopTable, file = "SIMD_Meta_BMI_Alc_AHRR_TopTable.csv")

SIMD_Meta_BMI_Alc_AHRR_TopTable <- read.csv("SIMD_Meta_BMI_Alc_AHRR_TopTable.csv")

SIMD_Meta_BMI_Alc_AHRR_pres_table <- head(SIMD_Meta_BMI_Alc_AHRR_TopTable[order(SIMD_Meta_BMI_Alc_AHRR_TopTable$P.value), ], 10)
View(SIMD_Meta_BMI_Alc_AHRR_pres_table)

write.csv(SIMD_Meta_BMI_Alc_AHRR_pres_table, file = "SIMD_BMI_Alc_AHRR_Meta_Top10.csv")


#creating a manhattan plot of the meta analysis
SIMD_Meta_BMI_Alc_AHRR_Man <- SIMD_Meta_BMI_Alc_AHRR_TopTable[,c(1,4,6,7,8)]
colnames(SIMD_Meta_BMI_Alc_AHRR_Man)[1] <- "SNP"
colnames(SIMD_Meta_BMI_Alc_AHRR_Man)[2] <- "P"
colnames(SIMD_Meta_BMI_Alc_AHRR_Man)[4] <- "CHR"
colnames(SIMD_Meta_BMI_Alc_AHRR_Man)[5] <- "BP"


SIMD_Meta_BMI_Alc_AHRR_Man$CHR <- gsub("chr", "", SIMD_Meta_BMI_Alc_AHRR_Man$CHR)
SIMD_Meta_BMI_Alc_AHRR_Man <- transform(SIMD_Meta_BMI_Alc_AHRR_Man, CHR = as.numeric(CHR))
SIMD_Meta_BMI_Alc_AHRR_Man <- transform(SIMD_Meta_BMI_Alc_AHRR_Man, BP = as.numeric(BP))
SIMD_Meta_BMI_Alc_AHRR_Man <- na.omit(SIMD_Meta_BMI_Alc_AHRR_Man)

View(SIMD_Meta_BMI_Alc_AHRR_Man)

manhattan((SIMD_Meta_BMI_Alc_AHRR_Man), col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-09)), main = "SIMD Meta (BMI, Alcohol, Smoking, AHRR)", cex.main = 1.5) #export with width 800 and keep aspect ratio







##heatmap of top 10 cpgs in BMI, Alc model w1 ####################

correlation_CpGs <- SIMD_Top10_CpGS_W1[,c(1:10)]

CpG_cor <- round(cor(correlation_CpGs), 2)

library(gplots)
heatmap.2(CpG_cor, Rowv=TRUE, Colv=TRUE,
          col=colorRampPalette(c("blue","white","red"))(n = 299), #colours-range and the number of individual colours in the palette (n=x)
          dendrogram="none", tracecol=F, key=F, 
          cellnote=CpG_cor, #cell labels taken from x
          notecol="black", #colour of cell labels
          cexRow=1, cexCol=1, #size of rows and coloumns
          margins = c(10,10), #size of margins
          lhei = c(1,12), lwid = c(0.5,6)) #change the ratio of the window that the graph covers




############## AHRR mediation analysis ##############################################################
install.packages("standardize")
library("standardize")
install.packages("lavaan")
library("lavaan")

SIMD_Top10_CpGS_W1 <- readRDS("./SIMDmodel_Top10_CpGs.rds")
SIMD_Top10_CpGS_W3 <- readRDS("./SIMD_Top10_CpGs_W3.rds")

genscot_merged <- read.csv("./genscot_merged.csv")

SIMD_cov_bmi_alcohol_w1 <- readRDS("./SIMD_cov_bmi_alcohol_w1.rds")
SIMD_cov_bmi_alcohol_w3 <- readRDS("./SIMD_cov_bmi_alcohol_w3.rds")

W1_PCs <- readRDS("./Factominer_PCs_5087_resid_PC100.rds")
W3_PCs <- readRDS("./wave3_pcs.rds")

#Getting W1 PCs ready
  #make rownames into a column
  W1_PCs <- cbind(rownames(W1_PCs), data.frame(W1_PCs, row.names=NULL))
  #rename ID column
  colnames(W1_PCs)[1] <- "ID"
  colnames(wave1_sentrix) [4] <- "ID"
  #merge with sentrix IDs
  W1_PCs <- merge(W1_PCs, wave1_sentrix[,c(1,4)], by = "ID")
  #delete Sample_Sentrix_ID column
  W1_PCs[,c(1,22:102)] <- NULL
  W1_PCs$PC1 <-scale(W1_PCs$PC1)
  W1_PCs$PC2 <-scale(W1_PCs$PC2)
  W1_PCs$PC3<-scale(W1_PCs$PC3)
  W1_PCs$PC4 <-scale(W1_PCs$PC4)
  W1_PCs$PC5 <-scale(W1_PCs$PC5)
  W1_PCs$PC6 <-scale(W1_PCs$PC6)
  W1_PCs$PC7 <-scale(W1_PCs$PC7)
  W1_PCs$PC8 <-scale(W1_PCs$PC8)
  W1_PCs$PC9 <-scale(W1_PCs$PC9)
  W1_PCs$PC10 <-scale(W1_PCs$PC10)
  W1_PCs$PC11 <-scale(W1_PCs$PC11)
  W1_PCs$PC12 <-scale(W1_PCs$PC12)
  W1_PCs$PC13 <-scale(W1_PCs$PC13)
  W1_PCs$PC14 <-scale(W1_PCs$PC14)
  W1_PCs$PC15 <-scale(W1_PCs$PC15)
  W1_PCs$PC16 <-scale(W1_PCs$PC16)
  W1_PCs$PC17 <-scale(W1_PCs$PC17)
  W1_PCs$PC18 <-scale(W1_PCs$PC18)
  W1_PCs$PC19 <-scale(W1_PCs$PC19)
  W1_PCs$PC20 <-scale(W1_PCs$PC20)
  
#Getting W3 PCs ready
  #merge with sentrix IDs
  W3_PCs <- merge(W3_PCs, wave3_sentrix[,c(1,4)], by = "Sample_Sentrix_ID")
  #delete Sample_Sentrix_ID column
  W3_PCs[,"Sample_Sentrix_ID"] <- NULL
  W3_PCs[,c(21:100)] <- NULL
  W3_PCs$PC1 <-scale(W3_PCs$PC1)
  W3_PCs$PC2 <-scale(W3_PCs$PC2)
  W3_PCs$PC3<-scale(W3_PCs$PC3)
  W3_PCs$PC4 <-scale(W3_PCs$PC4)
  W3_PCs$PC5 <-scale(W3_PCs$PC5)
  W3_PCs$PC6 <-scale(W3_PCs$PC6)
  W3_PCs$PC7 <-scale(W3_PCs$PC7)
  W3_PCs$PC8 <-scale(W3_PCs$PC8)
  W3_PCs$PC9 <-scale(W3_PCs$PC9)
  W3_PCs$PC10 <-scale(W3_PCs$PC10)
  W3_PCs$PC11 <-scale(W3_PCs$PC11)
  W3_PCs$PC12 <-scale(W3_PCs$PC12)
  W3_PCs$PC13 <-scale(W3_PCs$PC13)
  W3_PCs$PC14 <-scale(W3_PCs$PC14)
  W3_PCs$PC15 <-scale(W3_PCs$PC15)
  W3_PCs$PC16 <-scale(W3_PCs$PC16)
  W3_PCs$PC17 <-scale(W3_PCs$PC17)
  W3_PCs$PC18 <-scale(W3_PCs$PC18)
  W3_PCs$PC19 <-scale(W3_PCs$PC19)
  W3_PCs$PC20 <-scale(W3_PCs$PC20)
  

mediation_W1 <- merge(SIMD_Top10_CpGS_W1[,c(1,3,11)],  genscot_merged[,c("ID","rank")], by.x = "id", by.y = "ID")
mediation_W3 <- merge(SIMD_Top10_CpGS_W3[,c(1,3,11)],  genscot_merged[,c("ID","rank")], by.x = "id", by.y = "ID")

mediation_W1$cg05575921 <- scale(mediation_W1$cg05575921)
mediation_W1$cg15446156 <- scale(mediation_W1$cg15446156)

mediation_W3$cg05575921 <- scale(mediation_W3$cg05575921)
mediation_W3$cg15446156 <- scale(mediation_W3$cg15446156)

mediation_W1 <- merge(mediation_W1, SIMD_cov_bmi_alcohol_w1, by = "id")
mediation_W3 <- merge(mediation_W3, SIMD_cov_bmi_alcohol_w3, by = "id")

mediation_W1 <- merge(mediation_W1, W1_PCs, by = "id")
mediation_W3 <- merge(mediation_W3, W3_PCs, by = "id")


## Wave 1
model <- 'cg15446156 ~ c * rank + ever_smoke + pack_years + sex + age + bmi + usual + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19

cg05575921 ~ a * rank + ever_smoke + pack_years + sex + age + bmi + usual + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 +
PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19

cg15446156 ~ b * cg05575921
ab := a * b
total := c + (a * b)'
fit <- sem(model, data = mediation_W1)
mediation_result_W1 <- as.data.frame(unclass(parameterEstimates(fit)))
View(mediation_result_W1)



## Wave 3
model <- 'cg15446156 ~ c * rank + ever_smoke + pack_years + sex + age + bmi + CD8T + CD4T + NK + Bcell + Gran + usual + Batch +
PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +
PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20

cg05575921 ~ a * rank + ever_smoke + pack_years + sex + age + bmi + CD8T + CD4T + NK + Bcell + Gran + usual + Batch +
PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +
PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20

cg15446156 ~ b * cg05575921
ab := a * b
total := c + (a * b)'
fit <- sem(model, data = mediation_W3)
mediation_result_W3 <- as.data.frame(unclass(parameterEstimates(fit)))
View(mediation_result_W3)


#Meta analysis




#model <- ' # direct effect
#Y ~ c*X
#Y is the outcome - DNAm in this case
#X is the exposure - SIMD rank in this case

# mediator
#M ~ a*X
#M is the mediator - AHRR probe
#Y ~ b*M

# indirect effect (a*b)
#ab := a*b

# total effect
#total := c + (a*b)
#'
#fit <- sem(model, data = Data)
#summary(fit)


########### Correlation of smoking variabels ################################

tobacco_v2 <- read.csv("./tobaccov2.csv")
tobacco_v5 <- read.csv("./tobaccov5.csv")

tobacco_all <- rbind(tobacco_v2[,c(1:2, 11:15)], tobacco_v5[,c(1:2, 10:14)])
View(tobacco_all)

correlation_smoking_pheno <- merge(genscot_merged[,c(2,4)], tobacco_all[,c(1,2,7)], by = "ID")
View(correlation_smoking_pheno)

read.csv("./SIMD_Top10_CpGS_W1.csv")
read.csv("./SIMD_Top10_CpGS_W3.csv")

SIMD_Top10_CpGs <- rbind(SIMD_Top10_CpGS_W1, SIMD_Top10_CpGS_W3)
View(SIMD_Top10_CpGs)
colnames(SIMD_Top10_CpGs)[11] <- "ID"

correlation_smoking_pheno <- merge(correlation_smoking_pheno, SIMD_Top10_CpGs[,c(1,11)], by = "ID")

smoking_cor <- na.omit(correlation_smoking_pheno[,c(2:5)])

smoking_heat <- round(cor(smoking_cor), 2)
View(smoking_heat)




install.packages("ggcorrplot")
library(ggcorrplot)
ggcorrplot(smoking_heat, #DF with correlations
           hc.order = TRUE, #orders the variables according to size of cor.
           method = "square", #the shape of the tiles, can also be "circle"
           type = "lower", #only shows half the plot
           lab = TRUE) #adds the correlation coefficient to each tile






###################Dealing with ever_smoke as a factor etc##################################################################




SIMD_cov_BMI_Alc_AHRR_w3 <- readRDS("./SIMD_cov_BMI_Alc_AHRR_w3.rds")
View(SIMD_cov_BMI_Alc_AHRR_w3)
table(SIMD_cov_BMI_Alc_AHRR_w3$ever_smoke)

SIMD_cov_BMI_Alc_AHRR_w1 <- readRDS("./SIMD_cov_BMI_Alc_AHRR_w1.rds")
View(SIMD_cov_BMI_Alc_AHRR_w1)
table(SIMD_cov_BMI_Alc_AHRR_w1$ever_smoke)






wave1_modelE_cov <- readRDS("./wave1_modelE_covariates_5087_20181127.rds")
wave3_modelE_cov <- readRDS("./wave3_modelE_covariates_4450_20191120.rds")
View(wave1_modelE_cov)
View(wave3_modelE_cov)





#################### Making master cov file for normalised EWAS ###################################################

SIMD_cov_BMI_Alc_AHRR_w1 <- readRDS("./SIMD_cov_bmi_Alc_AHRR_w1.rds")
SIMD_cov_BMI_Alc_AHRR_w3 <- readRDS("./SIMD_cov_BMI_Alc_AHRR_w3.rds")
SIMD_cov_bmi_alcohol_w3 <- readRDS("./SIMD_cov_bmi_alcohol_w3.rds")
SIMD_cov_bmi_alcohol_w1 <- readRDS("./SIMD_cov_bmi_alcohol_w1.rds")

View(SIMD_cov_BMI_Alc_AHRR_w1)
View(SIMD_cov_BMI_Alc_AHRR_w3)

SIMD_cov_BMI_Alc_AHRR <- rbind(SIMD_cov_BMI_Alc_AHRR_w1[,c("id","ever_smoke","pack_years", "bmi", "units", "cg05575921")], SIMD_cov_BMI_Alc_AHRR_w3[,c("id","ever_smoke","pack_years", "bmi", "usual", "cg05575921")])
View(SIMD_cov_BMI_Alc_AHRR)
saveRDS(SIMD_cov_BMI_Alc_AHRR, file = "SIMD_cov_BMI_Alc_AHRR.rds")


SIMD_cov_BMI_Alc <- rbind(SIMD_cov_bmi_alcohol_w1[,c("id","ever_smoke","pack_years", "bmi", "units")], SIMD_cov_bmi_alcohol_w3[,c("id","ever_smoke","pack_years", "bmi", "units")])
View(SIMD_cov_BMI_Alc)
saveRDS(SIMD_cov_BMI_Alc, file = "SIMD_cov_BMI_Alc.rds")


#Opening EWAS output from normalised data
SIMD_BMI_Alc_norm_TopTable <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/New SIMD EWAS/SIMD_BMI_Alc_noprune.toptable.txt")
View(SIMD_BMI_Alc_norm_TopTable)

#Manhattan plot of SIMD EWAS in normalised data
SIMD_man_norm <- SIMD_BMI_Alc_norm_TopTable[,c(1,2,3,4,10)]
colnames(SIMD_man_norm)[1] <- "SNP"
colnames(SIMD_man_norm)[3] <- "CHR"
colnames(SIMD_man_norm)[4] <- "BP"
colnames(SIMD_man_norm)[5] <- "P"

SIMD_man_norm$CHR <- gsub("chr", "", SIMD_man_norm$CHR)
SIMD_man_norm <- transform(SIMD_man_norm, CHR = as.numeric(CHR))



manhattan(SIMD_man_norm, col = c(rgb(0,50,95,maxColorValue = 255), rgb(193,0,67, maxColorValue = 255)), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-11)), main = "SIMD (normalised data)", cex.main = 1.5) #export with width 800 and keep aspect ratio





#making a presentation table of norm. data EWAS
SIMD_BMI_Alc_norm_pres_table <- head(SIMD_BMI_Alc_norm_TopTable[order(SIMD_BMI_Alc_norm_TopTable$P.Value), ], 36)
View(SIMD_BMI_Alc_norm_pres_table)

write.csv(SIMD_BMI_Alc_norm_pres_table, file = "./New SIMD EWAS/SIMD_BMI_Alc_norm_Top36.csv")

#Correlating findings from the two analyses
SIMD_new_norm_correlation <- merge(SIMD_BMI_Alc_norm_TopTable[,c("ID", "beta")], SIMD_Meta_BMI_Alc_new_TopTable[,c("MarkerName","Zscore")], by.x = "ID", by.y = "MarkerName")
cor.test(SIMD_new_norm_correlation$beta, SIMD_new_norm_correlation$Zscore, method = "pearson")
View(SIMD_new_norm_correlation)

#Correlating top 10 CpGs
SIMD_new_norm_Top10_correlation <- merge(head(SIMD_BMI_Alc_norm_TopTable[order(SIMD_BMI_Alc_norm_TopTable$P.Value), c("ID", "beta")], 36), head(SIMD_Meta_BMI_Alc_new_TopTable[order(SIMD_Meta_BMI_Alc_new_TopTable$P.value), c(1, 3)], 36), by.x = "ID", by.y = "MarkerName")
cor.test(SIMD_new_norm_Top10_correlation$beta, SIMD_new_norm_Top10_correlation$Zscore, method = "pearson")
View(SIMD_new_norm_Top10_correlation)

#Correlating top 100 CpGs
SIMD_new_norm_Top100_correlation <- merge(head(SIMD_BMI_Alc_norm_TopTable[order(SIMD_BMI_Alc_norm_TopTable$P.Value), c("ID", "beta")], 100), head(SIMD_Meta_BMI_Alc_new_TopTable[order(SIMD_Meta_BMI_Alc_new_TopTable$P.value), c(1, 3)], 100), by.x = "ID", by.y = "MarkerName")
cor.test(SIMD_new_norm_Top100_correlation$beta, SIMD_new_norm_Top100_correlation$Zscore, method = "pearson")
View(SIMD_new_norm_Top100_correlation)

#Correlating top 1,000 CpGs
SIMD_new_norm_Top1000_correlation <- merge(head(SIMD_BMI_Alc_norm_TopTable[order(SIMD_BMI_Alc_norm_TopTable$P.Value), c("ID", "beta")], 1000), head(SIMD_Meta_BMI_Alc_new_TopTable[order(SIMD_Meta_BMI_Alc_new_TopTable$P.value), c(1, 3)], 1000), by.x = "ID", by.y = "MarkerName")
View(SIMD_new_norm_Top1000_correlation)
cor.test(SIMD_new_norm_Top1000_correlation$beta, SIMD_new_norm_Top1000_correlation$Zscore, method = "pearson")


#Correlating top 10,000 CpGs
SIMD_new_norm_Top10000_correlation <- merge(head(SIMD_BMI_Alc_norm_TopTable[order(SIMD_BMI_Alc_norm_TopTable$P.Value), c("ID", "beta")], 10000), head(SIMD_Meta_BMI_Alc_new_TopTable[order(SIMD_Meta_BMI_Alc_new_TopTable$P.value), c(1, 3)], 10000), by.x = "ID", by.y = "MarkerName", all = TRUE)
View(SIMD_new_norm_Top10000_correlation)
cor.test(SIMD_new_norm_Top10000_correlation$beta, SIMD_new_norm_Top10000_correlation$Zscore, method = "pearson")

#Correlating top 100,000 CpGs
SIMD_new_norm_Top100000_correlation <- merge(head(SIMD_BMI_Alc_norm_TopTable[order(SIMD_BMI_Alc_norm_TopTable$P.Value), c("ID", "beta")], 100000), head(SIMD_Meta_BMI_Alc_new_TopTable[order(SIMD_Meta_BMI_Alc_new_TopTable$P.value), c(1, 3)], 100000), by.x = "ID", by.y = "MarkerName", all = TRUE)
View(SIMD_new_norm_Top100000_correlation)
cor.test(SIMD_new_norm_Top100000_correlation$beta, SIMD_new_norm_Top100000_correlation$Zscore, method = "pearson")


#Finding top 10 Meta ouput in all of norm output
SIMD_new_Top10_in_norm <- merge(SIMD_BMI_Alc_norm_TopTable[,c("ID", "beta")], head(SIMD_Meta_BMI_Alc_new_TopTable[order(SIMD_Meta_BMI_Alc_new_TopTable$P.value), c(1, 3)], 10), by.x = "ID", by.y = "MarkerName")
cor.test(SIMD_new_Top10_in_norm$beta, SIMD_new_Top10_in_norm$Zscore, method = "pearson")
View(SIMD_new_Top10_in_norm)

#Finding top 36 norm ouput in all of meta output
SIMD_norm_Top36_in_new <- merge(head(SIMD_BMI_Alc_norm_TopTable[order(SIMD_BMI_Alc_norm_TopTable$P.Value), c("ID", "beta")], 36), SIMD_Meta_BMI_Alc_new_TopTable[, c(1, 3)], by.x = "ID", by.y = "MarkerName")
cor.test(SIMD_norm_Top36_in_new$beta, SIMD_norm_Top36_in_new$Zscore, method = "pearson")
View(SIMD_norm_Top36_in_new)




################# New BMI-Alc-Smoking EWAS with factor in cov ##############################
### EWAS including BMI and Alcohol
SIMD_BMI_Alc_W1_TopTable_new <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/New SIMD EWAS/SIMD_BMI_Alc_w1.toptable.txt")
SIMD_BMI_Alc_W3_TopTable_new <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/New SIMD EWAS/SIMD_BMI_Alc_w3.toptable.txt")
View(SIMD_BMI_Alc_W1_TopTable_new)
View(SIMD_BMI_Alc_W3_TopTable_new)






##Meta analysis

#creating wave 1 dataframes for METAL meta analysis
Metal_SIMD_BMI_Alc_W1_new <- SIMD_BMI_Alc_W1_TopTable_new[, c("ID", "P.Value", "N", "beta")]
View(Metal_SIMD_BMI_Alc_W1_new)

#adding reference allele columns to the metal data frame
Metal_SIMD_BMI_Alc_W1_new$A1 <- rep("A",nrow(Metal_SIMD_BMI_Alc_W1_new))
Metal_SIMD_BMI_Alc_W1_new$A2 <- rep("C",nrow(Metal_SIMD_BMI_Alc_W1_new))

#adding a direction of effect column
Metal_SIMD_BMI_Alc_W1_new$Effect <- ifelse(Metal_SIMD_BMI_Alc_W1_new$beta >= 0, "+", "-")

#creating wave 3 dataframes for METAL meta analysis
Metal_SIMD_BMI_Alc_W3_new <- SIMD_BMI_Alc_W3_TopTable_new[, c("ID", "P.Value", "N", "beta")]
View(Metal_SIMD_BMI_Alc_W3_new)

#adding reference allele columns to the metal data frame
Metal_SIMD_BMI_Alc_W3_new$A1 <- rep("A",nrow(Metal_SIMD_BMI_Alc_W3_new))
Metal_SIMD_BMI_Alc_W3_new$A2 <- rep("C",nrow(Metal_SIMD_BMI_Alc_W3_new))

#adding a direction of effect column
Metal_SIMD_BMI_Alc_W3_new$Effect <- ifelse(Metal_SIMD_BMI_Alc_W3_new$beta >= 0, "+", "-")

#save dataframes 
write.csv(Metal_SIMD_BMI_Alc_W1_new, file = "./New SIMD EWAS/Metal_SIMD_BMI_Alc_W1.csv", row.names = FALSE, quote = FALSE)
write.csv(Metal_SIMD_BMI_Alc_W3_new, file = "./New SIMD EWAS/Metal_SIMD_BMI_Alc_W3.csv", row.names = FALSE, quote = FALSE)

#Opening meta analysis
SIMD_BMI_Alc_meta_new_Output <- read.table("./New SIMD EWAS/Meta_SIMD_BMI_Alc_new1.tbl", header = TRUE, fill = TRUE)
View(SIMD_BMI_Alc_meta_new_Output)

#merging in gene annotation, chromosome and map info from both w1 and w3 dataframes
SIMD_Meta_BMI_Alc_new_TopTable <- merge(SIMD_BMI_Alc_meta_new_Output[,c("MarkerName", "Weight", "Zscore", "P.value", "Direction")], 
                                        SIMD_BMI_Alc_W1_TopTable_new[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_BMI_Alc_new_TopTable <- merge(SIMD_Meta_BMI_Alc_new_TopTable, SIMD_BMI_Alc_W3_TopTable_new[,c("ID","geneSymbol",  "CHR", "MAPINFO", "FEATURE", "CpGISLAND")], by.x = "MarkerName", by.y = "ID", all = TRUE)

SIMD_Meta_BMI_Alc_new_TopTable$geneSymbol.x <- as.character(SIMD_Meta_BMI_Alc_new_TopTable$geneSymbol.x)
SIMD_Meta_BMI_Alc_new_TopTable$geneSymbol.y <- as.character(SIMD_Meta_BMI_Alc_new_TopTable$geneSymbol.y)

SIMD_Meta_BMI_Alc_new_TopTable$geneSymbol.x[ is.na(SIMD_Meta_BMI_Alc_new_TopTable$geneSymbol.x) ] <- SIMD_Meta_BMI_Alc_new_TopTable$geneSymbol.y[ is.na(SIMD_Meta_BMI_Alc_new_TopTable$geneSymbol.x) ]
SIMD_Meta_BMI_Alc_new_TopTable$geneSymbol.y =NULL

SIMD_Meta_BMI_Alc_new_TopTable$CHR.x[ is.na(SIMD_Meta_BMI_Alc_new_TopTable$CHR.x) ] <- SIMD_Meta_BMI_Alc_new_TopTable$CHR.y[ is.na(SIMD_Meta_BMI_Alc_new_TopTable$CHR.x) ]
SIMD_Meta_BMI_Alc_new_TopTable$CHR.y =NULL

SIMD_Meta_BMI_Alc_new_TopTable$MAPINFO.x[ is.na(SIMD_Meta_BMI_Alc_new_TopTable$MAPINFO.x) ] <- SIMD_Meta_BMI_Alc_new_TopTable$MAPINFO.y[ is.na(SIMD_Meta_BMI_Alc_new_TopTable$MAPINFO.x) ]
SIMD_Meta_BMI_Alc_new_TopTable$MAPINFO.y =NULL

SIMD_Meta_BMI_Alc_new_TopTable$FEATURE.x <- as.character(SIMD_Meta_BMI_Alc_new_TopTable$FEATURE.x)
SIMD_Meta_BMI_Alc_new_TopTable$FEATURE.y <- as.character(SIMD_Meta_BMI_Alc_new_TopTable$FEATURE.y)

SIMD_Meta_BMI_Alc_new_TopTable$FEATURE.x[ is.na(SIMD_Meta_BMI_Alc_new_TopTable$FEATURE.x) ] <- SIMD_Meta_BMI_Alc_new_TopTable$FEATURE.y[ is.na(SIMD_Meta_BMI_Alc_new_TopTable$FEATURE.x) ]
SIMD_Meta_BMI_Alc_new_TopTable$FEATURE.y =NULL

SIMD_Meta_BMI_Alc_new_TopTable$CpGISLAND.x[ is.na(SIMD_Meta_BMI_Alc_new_TopTable$CpGISLAND.x) ] <- SIMD_Meta_BMI_Alc_new_TopTable$CpGISLAND.y[ is.na(SIMD_Meta_BMI_Alc_new_TopTable$CpGISLAND.x) ]
SIMD_Meta_BMI_Alc_new_TopTable$CpGISLAND.y =NULL

colnames(SIMD_Meta_BMI_Alc_new_TopTable)[6] <- "geneSymbol"
colnames(SIMD_Meta_BMI_Alc_new_TopTable)[7] <- "CHR"
colnames(SIMD_Meta_BMI_Alc_new_TopTable)[8] <- "MAPINFO"
colnames(SIMD_Meta_BMI_Alc_new_TopTable)[9] <- "FEATURE"
colnames(SIMD_Meta_BMI_Alc_new_TopTable)[10] <- "CpGISLAND"

View(SIMD_Meta_BMI_Alc_new_TopTable)

write.csv(SIMD_Meta_BMI_Alc_new_TopTable, file = "SIMD_Meta_BMI_Alc_new_TopTable.csv")

SIMD_Meta_BMI_Alc_new_TopTable <- read.csv("SIMD_Meta_BMI_Alc_new_TopTable.csv")

SIMD_Meta_BMI_Alc_new_pres_table <- head(SIMD_Meta_BMI_Alc_new_TopTable[order(SIMD_Meta_BMI_Alc_new_TopTable$P.value), ], 10)
View(SIMD_Meta_BMI_Alc_new_pres_table)

write.csv(SIMD_Meta_BMI_Alc_new_pres_table, file = "./New SIMD EWAS/SIMD_BMI_Alc_new_Meta_Top10.csv")


#creating a manhattan plot of the meta analysis
SIMD_Meta_BMI_Alc_Man <- SIMD_Meta_BMI_Alc_TopTable[,c(1,4,6,7,8)]
colnames(SIMD_Meta_BMI_Alc_Man)[1] <- "SNP"
colnames(SIMD_Meta_BMI_Alc_Man)[2] <- "P"
colnames(SIMD_Meta_BMI_Alc_Man)[4] <- "CHR"
colnames(SIMD_Meta_BMI_Alc_Man)[5] <- "BP"


SIMD_Meta_BMI_Alc_Man$CHR <- gsub("chr", "", SIMD_Meta_BMI_Alc_Man$CHR)
SIMD_Meta_BMI_Alc_Man <- transform(SIMD_Meta_BMI_Alc_Man, CHR = as.numeric(CHR))
SIMD_Meta_BMI_Alc_Man <- transform(SIMD_Meta_BMI_Alc_Man, BP = as.numeric(BP))
SIMD_Meta_BMI_Alc_Man <- na.omit(SIMD_Meta_BMI_Alc_Man)

View(SIMD_Meta_BMI_Alc_Man)

manhattan((SIMD_Meta_BMI_Alc_Man), col = c("black", "grey"), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-09)), main = "SIMD Meta (BMI, Alcohol, Smoking)", cex.main = 1.5) #export with width 800 and keep aspect ratio




##################### Running EWAS in norm data with wave as covariate '###################
SIMD_cov_BMI_Alc_AHRR_w1 <- readRDS("./SIMD_cov_bmi_Alc_AHRR_w1.rds")
SIMD_cov_BMI_Alc_AHRR_w3 <- readRDS("./SIMD_cov_BMI_Alc_AHRR_w3.rds")
SIMD_cov_bmi_alcohol_w3 <- readRDS("./SIMD_cov_bmi_alcohol_w3.rds")
SIMD_cov_bmi_alcohol_w1 <- readRDS("./SIMD_cov_bmi_alcohol_w1.rds")

View(SIMD_cov_BMI_Alc_AHRR_w1)
View(SIMD_cov_BMI_Alc_AHRR_w3)

#add a column with wave number
SIMD_cov_BMI_Alc_AHRR_w1_with_wave <- SIMD_cov_BMI_Alc_AHRR_w1
SIMD_cov_BMI_Alc_AHRR_w1_with_wave$wave <- 1

#add a column with wave number
SIMD_cov_BMI_Alc_AHRR_w3_with_wave <- SIMD_cov_BMI_Alc_AHRR_w3
SIMD_cov_BMI_Alc_AHRR_w3_with_wave$wave <- 2

#make a wave number and ID df
SIMD_wave_ID <- rbind(SIMD_cov_BMI_Alc_AHRR_w1_with_wave[,c("id","wave")], SIMD_cov_BMI_Alc_AHRR_w3_with_wave[,c("id","wave")])
View(SIMD_wave_ID)

#merge into SIMD normalised data, bmi and alcohol covariate file
SIMD_cov_BMI_Alc_wave <- merge(SIMD_cov_BMI_Alc, SIMD_wave_ID, by ="id")
View(SIMD_cov_BMI_Alc_wave)

#make wave column a factor
SIMD_cov_BMI_Alc_wave$wave <- as.factor(SIMD_cov_BMI_Alc_wave$wave)

#saving covariate file as rds file
saveRDS(SIMD_cov_BMI_Alc_wave, file = "SIMD_cov_BMI_Alc_wave.rds")




#Running the ewas in EDDIE

ssh s0951790@@eddie3.ecdf.ed.ac.uk

PATH=$PATH:/exports/igmm/eddie/GenScotDepression/local/bin
source ~/.bash_profile
module load R
stradl_ewas --pdata phenotypes.csv --pheno trait1 --out stradl_trait1_ewas
/exports/igmm/eddie/GenScotDepression/local/EWAS/nocpg_carmen_ewas --pdata SIMD.csv --pheno rank --cov /exports/eddie/scratch/s0951790/SIMD_cov_BMI_Alc_wave.rds --no-prune --out SIMD_BMI_Alc_Wave


#opening ewas pipeline results
SIMD_BMI_Alc_wave_TopTable <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/New SIMD EWAS/SIMD_BMI_Alc_Wave.toptable.txt")
View(SIMD_BMI_Alc_wave_TopTable)

#saving top results as .csv file
SIMD_BMI_Alc_wave_pres_table <- head(SIMD_BMI_Alc_wave_TopTable[order(SIMD_BMI_Alc_wave_TopTable$P.Value), ], 15)
View(SIMD_BMI_Alc_wave_pres_table)

write.csv(SIMD_BMI_Alc_wave_pres_table, file = "./New SIMD EWAS/SIMD_BMI_Alc_wave_Top10.csv")

#Manhattan plot of SIMD EWAS w. wave as cov
SIMD_man_wave <- SIMD_BMI_Alc_wave_TopTable[,c(1,2,3,4,10)]
colnames(SIMD_man_wave)[1] <- "SNP"
colnames(SIMD_man_wave)[3] <- "CHR"
colnames(SIMD_man_wave)[4] <- "BP"
colnames(SIMD_man_wave)[5] <- "P"

SIMD_man_wave$CHR <- gsub("chr", "", SIMD_man_wave$CHR)
SIMD_man_wave <- transform(SIMD_man_wave, CHR = as.numeric(CHR))



manhattan(SIMD_man_wave, col = c(rgb(0,50,95,maxColorValue = 255), rgb(193,0,67, maxColorValue = 255)), suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-11)), main = "SIMD (normalised data w. wave as covariate)", cex.main = 1.5) #export with width 800 and keep aspect ratio

library(aplot)
library(calibrate)
install.packages("ggrepel")
library(ggrepel)
install.packages("qqman")
library(qqman)

manhattan(SIMD_man_wave, col = c(rgb(0,50,95,maxColorValue = 255), rgb(193,0,67, maxColorValue = 255)), 
          suggestiveline = -log10(1e-5), genomewideline=-log10(3.6e-8), 
          ylim = c(0, -log10(1e-11)), main = "Scottish Index of Multiple Deprivation", cex.main = 1.5, 
          annotatePval =  3.6e-8, annotateTop = FALSE)+
   #export with width 800 and keep aspect ratio


geom_label_repel(data=topHits[topHits$P], aes(label=topHits$SNP, alpha=0.7), size=5, force=1.3))

textxy(pos, -log10(P), offset = 0.625, labs = topHits$SNP, cex = 1), ...)


sig = 3.6e-8
sugg = 1e-5

gg.manhattan(SIMD_man_wave, hlight = NA, threshold = sig,
             col = c("blue", "green"),
             ylims = c(0,12), title = "Scottish Index of Multiple Deprivation")


#correlate with findings wihtout wave as cov
SIMD_wave_correlation <- merge(SIMD_BMI_Alc_norm_TopTable[,c("ID", "beta")], SIMD_BMI_Alc_wave_TopTable[,c("ID", "beta")], by = "ID")
View(SIMD_wave_correlation)
cor.test(SIMD_wave_correlation$beta.x, SIMD_wave_correlation$beta.y, method = "pearson")

#paired t-test to see if the betas are significantly different
t.test(SIMD_wave_correlation$beta.x, SIMD_wave_correlation$beta.y, paired = TRUE, alternative = "two.sided")


SIMD_BMI_Alc_wave_pres_table

#correlate with findings wihtout wave as cov in top 15 
SIMD_wave_15_correlation <- merge(SIMD_BMI_Alc_wave_pres_table[,c("ID", "beta")], SIMD_BMI_Alc_norm_TopTable[,c("ID", "beta")], by = "ID")
View(SIMD_wave_15_correlation)
cor.test(SIMD_wave_15_correlation$beta.x, SIMD_wave_15_correlation$beta.y, method = "pearson")

#paired t-test to see if the betas are significantly different
t.test(SIMD_wave_15_correlation$beta.x, SIMD_wave_15_correlation$beta.y, paired = TRUE, alternative = "two.sided")

######## SIMD and depression regression
genscot_merged$dep_status <- as.factor(genscot_merged$dep_status)

SIMD_MDD <- genscot_merged[,c("ID", "ever_smoke", "pack_years", "dep_status", "rank", "quintile", "bmi", "age", "sex")]
SIMD_MDD <- merge(SIMD_MDD, Alcohol[,c("ID", "units")], by = "ID")
colnames(SIMD_MDD) [10] <- "Units"
base.model <- glm(dep_status ~ rank, data = SIMD_MDD, family=binomial())
summary(base.model)
full.model <- glm(dep_status ~ rank + Units + pack_years + bmi + age + sex, data = SIMD_MDD, family=binomial())
summary(full.model)

full.model.q <- glm(dep_status ~ quintile + Units + pack_years + bmi + age + sex, data = SIMD_MDD, family=binomial())
summary(full.model.q)


########## dmeographics table ############
SIMD_Demo <- na.omit(SIMD_cov_BMI_Alc_wave)
View(SIMD_Demo)
install.packages(psych)
library(psych)

###########  DMR Analysis   ##########################################################################################################################

library(devtools)
install_github("perishky/dmrff")
library(dmrff)

##creating a EWAS summary stats df w. estimate (logFC), standard error, p-value, chromosome and position
SIMD_BMI_Alc_wave_TopTable <- read.delim("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/New SIMD EWAS/SIMD_BMI_Alc_Wave.toptable.txt") #opening EWAS output
SIMD_EWAS_sumstat <- SIMD_BMI_Alc_wave_TopTable[,c("ID", "geneSymbol", "CHR", "MAPINFO", "logFC", "P.Value", "se")]
colnames(SIMD_EWAS_sumstat)[3] <- "chr"
colnames(SIMD_EWAS_sumstat)[4] <- "pos"
colnames(SIMD_EWAS_sumstat)[5] <- "estimate"
colnames(SIMD_EWAS_sumstat)[6] <- "p.value"
SIMD_EWAS_sumstat$chr <- sub("chr", "", SIMD_EWAS_sumstat$chr)
#Saving summary statistic df
write.csv(SIMD_EWAS_sumstat, file = "C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_EWAS_sumstat.csv")

##creating df of sentrix IDs for participants included in EWAS
DMR_Cov <- na.omit(SIMD_cov_BMI_Alc_wave)
colnames(DMR_Cov)[1] <- "ID"
colnames(sentrix_ID_w1_w3)[2] <- "ID"
colnames(cellcountinfo)[1] <- "Sample_Sentrix_ID"
sentrix_ID_w1_w3 <- rbind(cellcountinfo[,c(1,8,9,10)], cellcountinfo_w3[,c(1:3,7)])
DMR_Cov <- merge(DMR_Cov[,c(1,6)], sentrix_ID_w1_w3[,c(1,2)], by ="ID")
write.csv(DMR_Cov, file = "C:/Users/ajesp/Documents/PhD/SIMD_EWAS/DMR_IDs.csv")

#Copying files in Eddie
scp s0951790@eddie.ecdf.ed.ac.uk:/exports/igmm/eddie/GenScotDepression/miruna/EWAS_Jan/carmen_data/mvalues_txt/all_mvalues.txt s0951790@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/s0951790/

scp C:/Users/ajesp/Documents/PhD/SIMD_EWAS/DMR_IDs.csv s0951790@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/s0951790/
scp C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_EWAS_sumstat.csv s0951790@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/s0951790/
  
  
#requesting more memory
qlogin -l h_vmem=128G
#Opening R
Module load R
R


scp s0951790@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/s0951790/all_mvalues.txt C:/Users/ajesp/Documents/PhD/SIMD_EWAS/ 
  

#opening required packages
library(tidyverse)
library(tidyr)
library(dplyr)
library(devtools)
install_github("perishky/dmrff")
library(dmrff)

#Opening m-value file
m_vals <- read.delim("./all_mvalues.txt")
#saving as RDS file
#saveRDS(mvals, file = "./mvals.rds")
#Opening ID file
DMR_IDs <- read.csv("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/DMR_IDs.csv")
#Removing IDs from mvals file that are NOT part of the EWAS
mvals_DMR <- mvals[,which(colnames(mvals) %in% DMR_Cov$Sample_Sentrix_ID)]
#convert to matrix
mvals_DMR <- as.matrix(mvals_DMR)
#Opening summary stats file
SIMD_EWAS_sumstat <- read.csv("C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_EWAS_sumstat.csv")
#running dmrff function
dmrs <- dmrff(estimate=SIMD_EWAS_sumstat$estimate,
              se=SIMD_EWAS_sumstat$se,
              p.value=SIMD_EWAS_sumstat$p.value,
              methylation=mvals_DMR,
              chr=SIMD_EWAS_sumstat$chr,
              pos=SIMD_EWAS_sumstat$pos,
              maxgap=500,
              verbose=T)


#opening DMR file
DMRs_SIMD_adj <- read.csv("./dmrs_simd_adj.csv")
SIMD_EWAS_sumstat <- read.csv("./SIMD_EWAS_sumstat.csv")
View(DMRs_SIMD_adj)

#annotating DMRs
install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

data(list="IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
data(Locations)
data(Other)
annotation <- cbind(as.data.frame(Locations), as.data.frame(Other))
View(annotation)

dmrs_anno_test <- dmrff.annotate(DMRs_SIMD_adj, annotation)

library(ggplot2)
dmrff.manhattan.plot(DMRs_SIMD_adj, chr = DMRs_SIMD_adj$chr, pos = DMRs_SIMD_adj$, title = "Manhattan plot")



sites <- dmrff.sites(DMRs_SIMD_adj, SIMD_EWAS_sumstat$chr, SIMD_EWAS_sumstat$pos)
sites <- cbind(sites, SIMD_EWAS_sumstat[sites$site, c("ID", "geneSymbol", "estimate", "p.value")])
sites <- cbind(sites, dmr=DMRs_SIMD_adj[sites$region,c("start","end","z","p.adjust")])

sites$geneSymbol <- gsub(";.*", "", sites$geneSymbol)

sites <- aggregate(estimate ~ region*chr*geneSymbol*dmr.start*dmr.end*dmr.z*dmr.p.adjust, mean, data = sites)
View(sites)
write.csv(sites, file="./SIMD_DMR_genes.csv")

##DMR FUMA output
DMR_GeneSets <- read.delim("./GS.txt")

##Overlap of DMP CpGs in DMRs
SIMD_BMI_Alc_wave_Top15 <- read.csv("./New SIMD EWAS/SIMD_BMI_Alc_wave_Top10.csv")
CpG_overlap <- merge(sites, SIMD_BMI_Alc_wave_Top15, by = "ID")
View(CpG_overlap)




###### LASSO regression Generation Scotland W1 & W3 ######################################################################

SIMD_lasso_pheno_w1 <- readRDS("./SIMD_cov_bmi_alcohol_w1.rds")
PCs_w1 <- readRDS("./Factominer_PCs_5087_resid_PC100.rds")

IDinfo <- read.table("C:/Users/ajesp/Documents/PhD/1st year/Block 2/Mini project/Anders_PhD_Mini_Project/IDinfo.csv", header=TRUE, sep=",")
IDinfo$Sentrix_Position <- sub("^", "_", IDinfo$Sentrix_Position )
IDinfo$SampleID <- paste0(IDinfo$Sentrix_ID, IDinfo$Sentrix_Position)
PCs_w1 <- merge(PCs_w1, IDinfo[,c("Sample_Name", "SampleID")], by.x = "Sample_Sentrix_ID", by.y = "SampleID")
colnames(PCs_w1)[102] <- "ID"


colnames(SIMD_lasso_pheno_w1)[1] <- "ID"
SIMD_lasso_pheno_w1 <- merge(SIMD_lasso_pheno_w1,  genscot_merged[,c("ID", "rank")], by = "ID")
SIMD_lasso_pheno_w1 <- merge(SIMD_lasso_pheno_w1, PCs_w1[,c(2:11, 102)], by = "ID")
SIMD_lasso_pheno_w1_na <-na.omit(SIMD_lasso_pheno_w1)

SIMD_lasso_w1_lm <- lm(rank ~ ever_smoke + pack_years + bmi + usual + sex + age + PC1:PC10, data = SIMD_lasso_pheno_w1_na)
SIMD_lasso_w1_res <- data.frame(resid(SIMD_lasso_w1_lm))
SIMD_lasso_pheno_w1_na$resid <- resid(SIMD_lasso_w1_lm)
SIMD_lasso_pheno_w1 <-SIMD_lasso_pheno_w1_na

colnames(IDinfo)[1] <- "ID"
SIMD_lasso_pheno_w1 <- merge(SIMD_lasso_pheno_w1, IDinfo[,c(1,4)], by = "ID")
write.csv(SIMD_lasso_pheno_w1, file = "SIMD_lasso_pheno_aggr_w1.csv")

SIMD_lasso_pheno_w1 <- SIMD_lasso_pheno_w1[,c(1,19,20)]
colnames(SIMD_lasso_pheno_w1)[1] <- "gwas"
colnames(SIMD_lasso_pheno_w1)[2] <- "pheno"
colnames(SIMD_lasso_pheno_w1)[3] <- "ID"
write.csv(SIMD_lasso_pheno_w1, file = "SIMD_lasso_pheno_w1.csv")


#save penotype df as a SIMD_lasso_pheno_w1.txt file
write.table(SIMD_lasso_pheno_w1, file = "SIMD_lasso_pheno_w1.txt")

##In Eddie
  #Copy files into Eddie scratch space
  scp C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_lasso_pheno_w1.txt s0951790@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/s0951790/

  #requesting more memory
  qlogin -l h_vmem=128G

  #Opening R
  module load igmm/apps/R/3.6.3

#####Lasso script from Miruna (edited)

#-	FIRST, you need to open an R session and install.packages("glmnet"); 
  #this is because you need to interact with the installation (choose your own library 
  #from a drop down list supplied by R); then close the R session 
  #(but save the workspace when prompted by R); this will allow you to "library(glmnet)" 
  #when the actual Rscript runs

install.packages("glmnet")

## RSCRIPT FOR LASSO
## LASSO prediction in GS 5K 1st wave
# Pre-normalized beta values have been cleaned for unsuitable probes using Prep_LASSO_TKC.R
install.packages("methods")
library(glmnet)
library(methods)
args=commandArgs(TRUE)

clean_beta  <- read.table("/exports/igmm/eddie/GenScotDepression/miruna/EWAS_Jan/carmen_data/mvalues_txt/mvalues_450/all_450mvalues.txt")

names(clean_beta) <- substring(names(clean_beta),2,20) 

phenotype_file  <- read.table("/exports/eddie/scratch/s0951790/SIMD_lasso_pheno_w1.txt",header=T)

# input files
# clean_beta <- readRDS("/exports/eddie/scratch/tclarke2/pre_norm_beta_vals_clean450.Rdata")
# phenotype_file <- read.table("/exports/eddie/scratch/tclarke2/residual_phenos.txt", header=T)

# subset methylation data according to whether they have the predictor of interest
# phenotypes should be in the SENTRIXPLATE_SENTRIXPOSITION format and not GWAS id
pheno <- phenotype_file[!is.na(phenotype_file$pheno),]
a = which(colnames(clean_beta) %in% pheno$ID)
clean_beta1 = clean_beta[,a]
rm(clean_beta)
tclean_beta = t(clean_beta1)
rm(clean_beta1)

ids = rownames(tclean_beta)
dep1 = pheno[match(ids, pheno$ID),]
y = dep1$pheno
# remove columns with NAs as glmnet can't use missing data
x1 = tclean_beta[,apply(tclean_beta, 2, function(tar) !any(is.na(tar)))]
rm(tclean_beta) 

set.seed(1.234)

# Cross-validation: find the best shrinkage value - the alpha=1 means it's a LASSO model
lasso.cv <- cv.glmnet(x1, y, alpha=1, nfolds=10)
save(lasso.cv, file="/exports/eddie/scratch/s0951790/Lasso_output_SIMD_w1.RData") 

lambda.min = lasso.cv$lambda.min
m = glmnet(x1,y,lambda=lambda.min)

# create dataframe with weights for each probe
test = data.frame(coef.name = dimnames(coef(m))[[1]], coef.value = matrix(coef(m)))

# only take predictors where the co-efficient value is not 0
coef = test[test$coef.value!=0,]
save(coef, file="/exports/eddie/scratch/s0951790/Lasso_coef_SIMD_w1.RData") 

#End R and save workspace
q()
y


#copy files to own computer
scp s0951790@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/s0951790/Lasso_output_SIMD_w1.RData C:/Users/ajesp/Documents/PhD/SIMD_EWAS/


#make phenotype (rds) file with phenotype, ID and mehtylation ID
SIMD_predict_pheno_w3 <- readRDS("./SIMD_cov_bmi_alcohol_w3.rds")
colnames(SIMD_predict_pheno_w3)[1] <- "ID"

wave3_sentrix <- readRDS("./wave3_sentrix.rds")
colnames(wave3_sentrix)[1] <- "ID"

SIMD_predict_pheno_w3 <- merge(SIMD_predict_pheno_w3[,1], wave3_sentrix[,c(1,4)], by = "ID")
SIMD_predict_pheno_w3 <- SIMD_predict_pheno_w3[,c(1,14)]

SIMD_predict_pheno_w3 <- merge(genscot_merged[,c("ID", "rank")], SIMD_predict_pheno_w3, by = "ID")
colnames(SIMD_predict_pheno_w3)[1] <- "gwas"
colnames(SIMD_predict_pheno_w3)[2] <- "pheno"
colnames(SIMD_predict_pheno_w3)[3] <- "ID"

saveRDS(SIMD_predict_pheno_w3, file = "SIMD_predict_pheno_w3.rds")

##In Eddie
#Copy files into Eddie scratch space
scp C:/Users/ajesp/Documents/PhD/SIMD_EWAS/SIMD_predict_pheno_w3.rds s0951790@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/s0951790/
  
  #requesting more memory
  qlogin -l h_vmem=128G

#Opening R
module load igmm/apps/R/3.6.3
R

data  <- read.table("/exports/igmm/eddie/GenScotDepression/miruna/EWAS_Jan/carmen_data/mvalues_txt/mvalues_450/all_450mvalues.txt")
names(data) <- substring(names(data),2,20)

dep_test  <- readRDS("/exports/eddie/scratch/s0951790/SIMD_predict_pheno_w3.rds")

a = which(colnames(data) %in% dep_test$ID)
meth = data[,a]
rm(data)
dat = meth
rm(meth)

meth = t(dat)
rm(dat)
meth1 = as.data.frame(meth)
meth1$id = as.character(rownames(meth1)) 

## Read in LASSO coefficients, which contain coef.name (CpG site); coef.value (CpG weight) --> this was derived in Generation Scotland wave 1

load("/exports/eddie/scratch/s0951790/Lasso_coef_SIMD_w1.RData") 
a = which(names(meth1) %in% coef$coef.name)
meth2 = meth1[,a]
meth3 = t(meth2)
probes <- intersect(coef$coef.name, rownames(meth3))
rownames(coef) = coef$coef.name

b = meth3[probes,]
p = coef[probes,]

for (i in probes) {
  b[i,]= b[i,]*p[i,"coef.value"]
}

predicted_dep=colSums(b) + coef[1,2] 
pred_dep = as.data.frame(predicted_dep)
pred_dep$ID = rownames(pred_dep)

dep = merge(dep_test, pred_dep, by="ID")

## Save this dataset with antidep_wave3_scores.Rdata (or whichever name you want); this should have: participant ID; depression status; DNAm risk score as columns

save(dep, file="/exports/eddie/scratch/s0951790/SIMD_w3_DNAm_scores.RData") 
scp s0951790@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/s0951790/SIMD_w3_DNAm_scores.RData C:/Users/ajesp/Documents/PhD/SIMD_EWAS/

load("./SIMD_w3_DNAm_scores.RData")
View(dep)


##ROC tutorial

install.packages("pROC")
library(pROC)
install.packages("randomForest")
library(randomForest)

#use the roc function with known value first and estimate second.
roc(dep$pheno, dep$predicted_dep, plot=TRUE)
#didnt work because this is not a binary variable.

#Trying to use multiclass.roc instead
multiclass.roc(dep$pheno, dep$predicted_dep, plot = F)
#


################ Rural living and depression ################################
install.packages("readxl")
library(readxl)
SIMD_urban_rural <- read_excel("~/PhD/SIMD_EWAS/SIMD urban rural.xlsx")
colnames(SIMD_urban_rural)[1] <- "ID"



urban_depression <- genscot_merged[,c("ID", "dep_status", "inflam_disease", "rank", "bmi", "age", "sex")]

urban_depression <- merge(urban_depression, SIMD_urban_rural[,c(1,4)], by = "ID")
urban_depression$dep_status <- as.factor(urban_depression$dep_status)
urban_depression$urban <- as.factor(urban_depression$urban)



fit <- glm(inflam_disease ~ urban + rank + age + sex, data = urban_depression)
summary(fit)

fit <- glm(inflam_disease ~ urban + age + sex, data = urban_depression, family = binomial())
summary(fit)

fit <- glm(dep_status ~ urban + rank + age + sex, data = urban_depression, family = binomial())
summary(fit)

fit <- glm(dep_status ~ inflam_disease + rank + age + sex, data = urban_depression, family = binomial())
summary(fit)


cor.test(urban_depression$rank, urban_depression$urban, method = "pearson")
cor.test(urban_depression$urban, urban_depression$, method = "pearson")


cor.test(genscot_merged$rank, genscot_merged$Bcell, method = "pearson")



######## SIMD and depression regression
genscot_merged$dep_status <- as.factor(genscot_merged$dep_status)

SIMD_MDD <- genscot_merged[,c("ID", "ever_smoke", "pack_years", "dep_status", "rank", "quintile", "bmi", "age", "sex")]
SIMD_MDD <- merge(SIMD_MDD, Alcohol[,c("ID", "units")], by = "ID")
colnames(SIMD_MDD) [10] <- "Units"
base.model <- glm(dep_status ~ rank, data = SIMD_MDD, family=binomial())
summary(base.model)
full.model <- glm(dep_status ~ rank + Units + pack_years + bmi + age + sex, data = SIMD_MDD, family=binomial())
summary(full.model)

full.model.q <- glm(dep_status ~ quintile + Units + pack_years + bmi + age + sex, data = SIMD_MDD, family=binomial())
summary(full.model.q)

