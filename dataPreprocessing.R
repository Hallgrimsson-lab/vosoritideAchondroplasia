# load in, symmetrize, and register achondroplasia drug participants, FB2 achondroplasia, and 3Dfn data

# Achondroplasia drug data####
#covariates
# drugCovariates <- read.csv("~/postdoc/achondroplasiaDrug/achondroplasia_drug_covariates.csv")
drugCovariates <- read.csv("~/postdoc/achondroplasiaDrug/ACH_covariates.csv")

# orderedDrugCovariates <- drugCovariates[match(tools::file_path_sans_ext(basename(drugTrial)), drugCovariates$Image.ID),]
orderedDrugCovariates <- drugCovariates[match(tools::file_path_sans_ext(basename(drugTrial)), drugCovariates$Image_ID),]

drugTrial <- list.files("/Volumes/Storage4/FaceBase_3/Data/Images/Misc/Achondroplasia_Drug/2_Registered/Dense_5k/Scans/Original/", pattern = "*.ply", full.names = T)
drugTrialReflected <- list.files("/Volumes/Storage4/FaceBase_3/Data/Images/Misc/Achondroplasia_Drug/2_Registered/Dense_5k/Scans/Reflected/", pattern = "*.ply", full.names = T)
drugScans <- array(NA, dim = c(ncol(atlas$vb), 3, length(drugTrial)))
sample1k <- sample(1:ncol(atlas$vb), 200)

for(i in 1:dim(drugScans)[3]){

  tmpScan <- file2mesh(drugTrial[i])
  tmpScan <- rotmesh.onto(tmpScan, t(tmpScan$vb[-4, sample1k]), t(atlas$vb[-4, sample1k]), scale =
                            T)$mesh

  tmpScanReflected <- file2mesh(drugTrialReflected[i])
  tmpScanReflected <- rotmesh.onto(tmpScanReflected, t(tmpScanReflected$vb[-4, sample1k]), t(atlas$vb[-4, sample1k]), scale =
                                     T)$mesh

  #symmetrize
  tmpScan$vb <- (tmpScanReflected$vb + tmpScan$vb)/2

  drugScans[,,i] <- t(rotmesh.onto(tmpScan, t(tmpScan$vb[-4, sample1k]), sampleMeanShape[sample1k,], scale =
                                     T)$mesh$vb[-4,])

}

dimnames(drugScans)[[3]] <- tools::file_path_sans_ext(basename(drugTrial))

drugScans <- drugScans[,,is.na(orderedDrugCovariates$Image_ID) == F]
finalDrugCovariates <- orderedDrugCovariates[is.na(orderedDrugCovariates$Image_ID) == F,]


# 3d facial norms####
filteredFbmd <- fbMd[fbMd$Syndrome %in% c("Achondroplasia", "Unaffected Unrelated"), ]
filteredFbmd <- filteredFbmd[match(filteredHdrda$Family, filteredFbmd$Family_ID),]


pittScans <- list.files("/Volumes/Storage4/FaceBase_3/Data/Images/Pitt/2_Registered/Dense_5k/Meshes/Original/", pattern = "*.ply", full.names = T)
pittScansReflected <- list.files("/Volumes/Storage4/FaceBase_3/Data/Images/Pitt/2_Registered/Dense_5k/Meshes/Reflected/", pattern = "*.ply", full.names = T)
pittScanNames <- tools::file_path_sans_ext(basename(pittScans))

fb2ShapeSymmetrized <- array(NA, dim = c(nverts(atlas), 3, length(pittScans)))

for(i in 1:nrow(filteredHdrda)){
  tmpScan <- file2mesh(pittScans[i])
  tmpScan <- rotmesh.onto(tmpScan, t(tmpScan$vb[-4, sample1k]), t(atlas$vb[-4, sample1k]), scale =
                            T)$mesh

  tmpScanReflected <- file2mesh(pittScansReflected[i])
  tmpScanReflected <- rotmesh.onto(tmpScanReflected, t(tmpScanReflected$vb[-4, sample1k]), t(atlas$vb[-4, sample1k]), scale =
                                     T)$mesh

  #symmetrize
  tmpScan$vb <- (tmpScanReflected$vb + tmpScan$vb)/2

  fb2ShapeSymmetrized[,,i] <- t(rotmesh.onto(tmpScan, t(tmpScan$vb[-4, sample1k]), sampleMeanShape[sample1k,], scale =
                                               T)$mesh$vb[-4,])

  if(i %% 500 == F) print(i)
}

pittShapeSymmetrized <- fb2ShapeSymmetrized
# filter for pitt that we have metadata for
matchedPittMd <- fbMd[match(pittScanNames, fbMd$Image_Name),]
View(cbind(pittScanNames, matchedPittMd))

# fb2 achondroplasia####
achondScans <- list.files("/Volumes/Storage4/FaceBase_3/Data/Images/AchSubset/2_Registered/Dense_5k/Meshes/Original/", pattern = "*.ply", full.names = T)
achondScanNames <- tools::file_path_sans_ext(basename(achondScans))

achondMd <- fbMd[fbMd$Syndrome == "Achondroplasia",]

matchedAchondMd <- achondMd[na.omit(match(achondScanNames, achondMd$Image_Name)),]
View(matchedAchondMd)

achondScans <- achondScans[match(matchedAchondMd$Image_Name, achondScanNames)]
achondScansReflected <- list.files("/Volumes/Storage4/FaceBase_3/Data/Images/AchSubset/2_Registered/Dense_5k/Meshes/Reflected/", pattern = "*.ply", full.names = T)
#matching assumes there are no failed reflections--beware
achondScansReflected <- achondScansReflected[match(matchedAchondMd$Image_Name, achondScanNames)]


achondShapeSymmetrized <- array(NA, dim = c(nverts(atlas), 3, length(achondScans)))

for(i in 1:length(achondScans)){
  tmpScan <- file2mesh(achondScans[i])
  tmpScan <- rotmesh.onto(tmpScan, t(tmpScan$vb[-4, sample1k]), t(atlas$vb[-4, sample1k]), scale =
                            T)$mesh

  tmpScanReflected <- file2mesh(achondScansReflected[i])
  tmpScanReflected <- rotmesh.onto(tmpScanReflected, t(tmpScanReflected$vb[-4, sample1k]), t(atlas$vb[-4, sample1k]), scale =
                                     T)$mesh

  #symmetrize
  tmpScan$vb <- (tmpScanReflected$vb + tmpScan$vb)/2

  achondShapeSymmetrized[,,i] <- t(rotmesh.onto(tmpScan, t(tmpScan$vb[-4, sample1k]), sampleMeanShape[sample1k,], scale =
                                                  T)$mesh$vb[-4,])
}

#save out all metadata and landmarks
save(finalDrugCovariates, drugScans, matchedAchondMd, achondShapeSymmetrized, matchedPittMd, pittShapeSymmetrized, file = "~/postdoc/achondroplasiaDrug/achondDrugStartingShapeMd.Rdata")

