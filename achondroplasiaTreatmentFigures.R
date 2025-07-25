#### Achondroplasia drug severity figure
# this script assumes you've already symmetrized all the data and prepped the covariate info
# we also used parallelPCA on 3Dfn to generate the PC axes for this analysis

# load in 3dfn PCA
load("~/postdoc/achondroplasiaDrug/symmetricPittPCA.Rdata")
load("~/postdoc/achondroplasiaDrug/achondDrugStartingShapeMd.Rdata")

nPCs <- symmetricPCA$n


nsCol <-  "#3F51B5"
syndCol <- "#F44336"

filteredHdrda <- rbind.data.frame(matchedAchondMd, matchedPittMd)
filteredHdrda$Syndrome <- as.character(filteredHdrda$Syndrome)
finalDrugControls <- finalDrugCovariates[finalDrugCovariates$Group == "Control",]
finalDrugDemographics <- finalDrugCovariates[finalDrugCovariates$Visit == 1,]

achondDataFull <- data.frame(Sex = c((filteredHdrda$Sex == "F") - 1, (finalDrugDemographics$Sex == "F") - 1), Syndrome = c(filteredHdrda$Syndrome, finalDrugDemographics$Group), Age = c(filteredHdrda$Age, finalDrugDemographics$Age), Group = c(filteredHdrda$Syndrome, finalDrugDemographics$Phase))
achondDataFull <- achondDataFull[is.na(achondDataFull$Age) == F,]
achondDataFull$Group[is.na(achondDataFull$Group)] <- "Achondroplasia"
achondDataFull$Group[achondDataFull$Group == "Achondroplasia"] <- "Achondroplasia Control"
achondDataFull$Group[achondDataFull$Group == "3" & achondDataFull$Syndrome == "Drug"] <- "Achondroplasia Phase 3 (Drug)"
achondDataFull$Group[achondDataFull$Group == "3" & achondDataFull$Syndrome == "Placebo"] <- "Achondroplasia Phase 3 (Placebo)"
achondDataFull$Group[achondDataFull$Group == "2"] <- "Achondroplasia Phase 2"
achondDataFull$Group <- factor(achondDataFull$Group, levels = c("Unaffected Unrelated", "Achondroplasia Control", "Achondroplasia Phase 2", "Achondroplasia Phase 3 (Drug)", "Achondroplasia Phase 3 (Placebo)"))

# figure 1 - demographics####
# a) count by sex

demographicDf <- data.frame(table(achondDataFull$Group, achondDataFull$Sex))
demographicDf$Freq[demographicDf$Freq > 55] <- 55
demographicDf$Var1 <- factor(demographicDf$Var1, levels = c("Unaffected Unrelated", "Achondroplasia Control", "Achondroplasia Phase 2", "Achondroplasia Phase 3 (Drug)", "Achondroplasia Phase 3 (Placebo)"))


p <- ggplot(demographicDf, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(labels = c("Male", "Female"), values = c("#0E084C","#3D5FA6")) +
  ggbreak::scale_y_break(c(7,90)) +
  labs(x = "", y = "Sample Size", fill = "Sex") +
  theme_minimal(base_line_size = 0) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1.02, size = 8.5),
        axis.text.y = element_text(size = 6),
        legend.position = "bottom",  # Position legend at top right (relative coordinates)
        legend.direction = "horizontal",  # Set legend to horizontal
        legend.title = element_text(size = 8),  # Adjust legend title size
        legend.text = element_text(size = 8))

pdf(file = "~/postdoc/achondroplasiaDrug/sampleSizes.pdf", height = 4, width = 4)
p
dev.off()


# b) age density by group
p <- ggplot(achondDataFull, aes(x = Age, y = Group)) +
  ggridges::geom_density_ridges(scale = .8, rel_min_height = 0.01, fill = "black") +
  theme_minimal(base_line_size = 0) +
  theme(axis.text.x = element_text(hjust = 1, size = 8.5, angle = 60)) +
  labs(title = "", x = "Age (years)", y = "") +
  coord_flip()

pdf(file = "~/postdoc/achondroplasiaDrug/ageRidgePlot.pdf", width = 4, height = 4)
p
dev.off()

# project into pc space
drugScores <- getPCscores(drugScans, as.matrix(symmetricPCA$original$loadings[,1:nPCs]), as.matrix(sampleMeanShape))
achondScores <- getPCscores(achondShapeSymmetrized, as.matrix(symmetricPCA$original$loadings[,1:nPCs]), as.matrix(sampleMeanShape))
pittScores <- getPCscores(pittShapeSymmetrized, as.matrix(symmetricPCA$original$loadings[,1:nPCs]), as.matrix(sampleMeanShape))

plot(symmetricPCA$original$rotated[,1:2])
points(drugScores[,1:2], col = 2)
points(achondScores[,1:2], col = 4)

#model
achondData <- data.frame(Sex = c((filteredHdrda$Sex == "F") - 1, (finalDrugControls$Sex == "F") - 1), Syndrome = c(filteredHdrda$Syndrome, (rep("Achondroplasia", nrow(finalDrugControls)))), Age = c(filteredHdrda$Age, finalDrugControls$Age), Scores = rbind(achondScores, pittScores, drugScores[finalDrugCovariates$Group == "Control",]))
achondData <- achondData[is.na(achondData$Age) == F,]
achondData$Syndrome <- factor(achondData$Syndrome, levels = c("Unaffected Unrelated", "Achondroplasia"))

View(achondData)

cubicAge <- poly(achondData$Age, degree = 3)
achondModel <- lm(as.matrix(achondData[, -c(1:3)]) ~ cubicAge[,1] + cubicAge[,2] + cubicAge[,3] + achondData$Sex +  achondData$Syndrome + achondData$Syndrome:cubicAge[,1])
severityScores <- Morpho::RegScore(achondModel)

#project achondroplasia drug participants into score space
# severity score distributions####
severityDrugs <- RegScore(achondModel, drugScores)[,5]

syndromeOfInterest <- "Achondroplasia"

syndromeScoresF <- density(severityScores[achondData$Syndrome == syndromeOfInterest & achondData$Sex == 0, 5]/sd(severityScores[achondData$Syndrome == syndromeOfInterest, 5]))
syndromeScoresM <- density(severityScores[achondData$Syndrome == syndromeOfInterest & achondData$Sex == -1, 5]/sd(severityScores[achondData$Syndrome == syndromeOfInterest, 5]))
syndromeScoresM$y <- -syndromeScoresM$y

syndromeScoresUUF <- density(severityScores[achondData$Syndrome == "Unaffected Unrelated" & achondData$Sex == 0, 5]/sd(severityScores[achondData$Syndrome == syndromeOfInterest, 5]))
syndromeScoresUUM <- density(severityScores[achondData$Syndrome == "Unaffected Unrelated" & achondData$Sex == -1, 5]/sd(severityScores[achondData$Syndrome == syndromeOfInterest, 5]))
syndromeScoresUUM$y <- -syndromeScoresUUM$y

# Standardize scores and get density for achondroplasia drug individuals

syndromeScoresDrugsTreatmentF <- density(severityDrugs[finalDrugCovariates$Sex == "F" & finalDrugCovariates$Group == "Drug"]/sd(severityScores[achondData$Syndrome == syndromeOfInterest, 5]))
syndromeScoresDrugsTreatmentM <- density(severityDrugs[finalDrugCovariates$Sex == "M" & finalDrugCovariates$Group == "Drug"]/sd(severityScores[achondData$Syndrome == syndromeOfInterest, 5]))
syndromeScoresDrugsTreatmentM$y <- -syndromeScoresDrugsTreatmentM$y

#make the pdf bar for the severity plot
pdf(file = glue("~/postdoc/achondroplasiaDrug/SeveritySegments.pdf"), width = 7, height = 4)

unaffectedCol <- nsCol
plot(1:10, type = "n", xlim = c(-3.5,7), ylim = c(-1,2.1), xlab = glue("{syndromeOfInterest} score (sd)"), ylab = "", axes = F, cex = .4)

# draw x axis line in the middle
segments(-3, 0, 7, 0)

# draw actual x axis with ticks 10% lower than the minimum density function value
axis(1, pos = 1.5 * min(min(syndromeScoresM$y), min(syndromeScoresUUM$y)), lwd = .75)  # X-axis
text(cbind(c(-2.5,-2.5), c(.075, -.075)), labels = c("Female", "Male"), cex = .6)

polygon(syndromeScoresUUF, col = adjustcolor(unaffectedCol, .3))
polygon(syndromeScoresUUM, col = adjustcolor(unaffectedCol, .3))

polygon(syndromeScoresF, col = adjustcolor(syndCol, .3))
polygon(syndromeScoresM, col = adjustcolor(syndCol, .3))

#for every ind id that was treated, plot the beginning and end severity scores####
treatedIndex <- which(finalDrugCovariates$Group == "Drug")
scoresTreated <- severityDrugs[treatedIndex]
covariateTreated <- finalDrugCovariates[treatedIndex,]

treatedIndex2 <- which(finalDrugCovariates$Group == "Placebo")
scoresTreated2 <- severityDrugs[treatedIndex2]
covariateTreated2 <- finalDrugCovariates[treatedIndex2,]

yIterator <- seq(170, 70, length.out = length(unique(covariateTreated2$Individual)) + length(unique(covariateTreated$Individual)))/100

for(i in 1:length(unique(covariateTreated$Individual))){
  indIndex <- which(covariateTreated$Individual == unique(covariateTreated$Individual)[i])

    firstVisit <- which.min(covariateTreated$Visit[indIndex])
    lastVisit <- which.max(covariateTreated$Visit[indIndex])
    #plot both points
    tmpPoints <- matrix(c(scoresTreated[indIndex][c(firstVisit, lastVisit)]/sd(severityScores[achondData$Syndrome == syndromeOfInterest, 5]), c(yIterator[i+4], yIterator[i+4])), ncol = 2)
    points(tmpPoints, pch = 19, col = syndCol, cex = .6)
    arrows(tmpPoints[1,1], tmpPoints[1,2], tmpPoints[2,1], tmpPoints[2,2], lty = 1, col = 1, length = .05, angle = 20)
    print(glue("{unique(covariateTreated$Individual)[i]} diff: {diff(c(scoresTreated[indIndex][c(firstVisit, lastVisit)]/sd(severityScores[achondData$Syndrome == syndromeOfInterest, 5])))}"))
    }


for(i in 1:length(unique(covariateTreated2$Individual))){
    indIndex <- which(covariateTreated2$Individual == unique(covariateTreated2$Individual)[i])
    firstVisit <- which(covariateTreated2$Visit[indIndex] == 2)
    lastVisit <- which.max(covariateTreated2$Visit[indIndex])
    #plot both points
    tmpPoints <- matrix(c(scoresTreated2[indIndex][c(firstVisit, lastVisit)]/sd(severityScores[achondData$Syndrome == syndromeOfInterest, 5]), c(yIterator[i], yIterator[i])), ncol = 2)
    points(tmpPoints, pch = 19, col = syndCol, cex = .6)
    arrows(tmpPoints[1,1], tmpPoints[1,2], tmpPoints[2,1], tmpPoints[2,2], lty = 1, col = 1, length = .05, angle = 20)
    print(glue("{unique(covariateTreated2$Individual)[i]} diff: {diff(c(scoresTreated2[indIndex][c(firstVisit, lastVisit)]/sd(severityScores[achondData$Syndrome == syndromeOfInterest, 5])))}"))

  }


dev.off()

#BK, TH, AR are the biggest changes in score

# severity spectrum morphs####
# convert numbers into z-scores
tmpAge <- predict(cubicAge, 9)
tmpSex <- -.5

tmpSyndrome <- factor("Achondroplasia", levels = c("Unaffected Unrelated", "Achondroplasia"))

# build syndrome model dataframe
datamod <- ~ tmpAge[1] + tmpAge[2] + tmpAge[3] + tmpSex + tmpSyndrome + tmpSyndrome:tmpAge[1]

# predict PC scores from model
mat <- model.matrix(datamod)
# shapePred <- predPcScores(achondModel, mat, syndrome = syndrome, severity = severity)

syndromeVector <- achondModel$model$`achondData$Syndrome`
# severitySD <- apply(Morpho::RegScore(syndromeModel), MARGIN = 2, FUN = sd)
severityScores <- Morpho::RegScore(achondModel)

achondScores <- severityScores[,5]
nsMean <- mean(achondScores[achondData$Syndrome == "Unaffected Unrelated"])
achondMean <- mean(achondScores[achondData$Syndrome == "Achondroplasia"])

normalizedAchondScores <- (achondScores - nsMean) / (achondMean - nsMean) #* (achondMean - nsMean) + nsMean

#sanity check on the normalization
mean(normalizedAchondScores[achondData$Syndrome == "Achondroplasia"])
boxplot(normalizedAchondScores ~ achondData$Syndrome)

syndromeColumn <- 5
mOrF <- -.5

severity <- 6
severityMag <- (severity * (sd(normalizedAchondScores)))
severitySeq <- seq(min(normalizedAchondScores), max(normalizedAchondScores), length.out = 5)

for(i in 1:length(severitySeq)){
  severityMag <- severitySeq[i]
  mat[6] <- 0 + severityMag # main effect
  mat[7] <- mat[2] * (0 + severityMag) #age interaction

  shapePred <- as.numeric(mat) %*% achondModel$coef

  # rescale predicted pc scores
  shapePredRescaled <- rep(NA, length(shapePred))
  for(i in 1:length(shapePredRescaled)) shapePredRescaled[i] <- (shapePred[i]) #* sd(pcScores[,i])) + mean(pcScores[,i])

  finalEstimate <- Morpho::restoreShapes(shapePredRescaled, symmetricPCA$original$loadings[, 1:numPCs], sampleMeanShape)
  # finalEstimate2 <- finalEstimate + severityComponent

  finalEstimateMesh <- atlas
  finalEstimateMesh$vb[-4,] <- t(finalEstimate)
  finalEstimateMesh <- Rvcg::vcgUpdateNormals(finalEstimateMesh)

  # nsMesh <- finalEstimateMesh # i=2
  # achondMesh <- finalEstimateMesh # i=4

  plotSyndromeModel(finalEstimateMesh)
  par3d(userMatrix = diag2)
  rgl.snapshot(glue("~/postdoc/achondroplasiaDrug/AchondScore{round(severityMag, 2)}.png"))

}

#is change in score related to change in age?
# needs to be rewritten to accomodate for variable numbers of visits
# ageChange <- rep(NA, length(unique(covariateTreated)) - 1)
# scoreChange <- rep(NA, length(unique(covariateTreated)) - 1)
#
# for(i in 1:length(unique(covariateTreated$Individual.ID))){
#   indIndex <- which(covariateTreated$Individual.ID == unique(covariateTreated$Individual.ID)[i])
#   if(length(indIndex) > 1){
#
#     scoreChange[i] <- diff(scoresTreated[indIndex]/sd(severityScores[achondData$Syndrome == syndromeOfInterest, 5]))
#     ageChange[i] <- diff(covariateTreated$Age[indIndex])
#
#   }
# }

plot(scoreChange ~ ageChange)

# per vertex achondroplasia scores####
#use drugVertexScores to compare the change in per vertex achondroplasia scores
# find the first visit for drug, generate heatmap

#find the last visit for drug, subtract from first
treatedIndex <- which(finalDrugCovariates$Group == "Drug")
scoresTreated <- severityDrugs[treatedIndex]
covariateTreated <- finalDrugCovariates[treatedIndex,]
finalDiffVector <- rep(NA, length(unique(covariateTreated$Individual)))

for(i in 1:length(unique(covariateTreated$Individual))){

  syndromeEffect <- achondModel$coefficients[5 ,] + achondModel$coefficients[1,]

  interceptShape <- matrix(restoreShapes(achondModel$coefficients[1,], symmetricPCA$original$loadings[,1:nPCs], sampleMeanShape), ncol = 3, byrow = F)

  syndVector <- matrix(restoreShapes(syndromeEffect, symmetricPCA$original$loadings[,1:nPCs], sampleMeanShape), ncol = 3, byrow = F) - interceptShape

  indIndex <- which(covariateTreated$Individual == unique(covariateTreated$Individual)[i])

  firstVisit <- which.min(covariateTreated$Visit[indIndex])
  lastVisit <- which.max(covariateTreated$Visit[indIndex])
  firstVisitVerts <- restoreShapes(drugScores[indIndex,][firstVisit,], symmetricPCA$original$loadings[,1:nPCs], sampleMeanShape)
  lastVisitVerts <- restoreShapes(drugScores[indIndex,][lastVisit,], symmetricPCA$original$loadings[,1:nPCs], sampleMeanShape)

  personAligned <- rotmesh.onto(makeMesh(firstVisitVerts), firstVisitVerts, interceptShape, scale = T)

  personVector <- personAligned$yrot - interceptShape
  syndPersonVector <- lastVisitVerts - interceptShape
  # open3d()
  par3d(windowRect = c(0, 0, 800, 800), zoom = .75, userMatrix = diag1)
  # orient face, then save front.face <- par3d()$userMat
  firstProjection <- perVertexHeatmap(atlas, syndVector, personVector, type = "projection")
  rgl.snapshot(glue("~/postdoc/achondroplasiaDrug/individualScores/{unique(covariateTreated$Individual)[i]}_firstVisitScore_diag1.png"))

  par3d(userMatrix = diag2)
  rgl.snapshot(glue("~/postdoc/achondroplasiaDrug/individualScores/{unique(covariateTreated$Individual)[i]}_firstVisitScore_diag2.png"))

    #then plot for syndromic
  lastProjection <- perVertexHeatmap(atlas, syndVector, syndPersonVector, type = "projection")
  par3d(userMatrix = diag1)

  finalDiff <- lastProjection - firstProjection
  finalDiffVector[i] <- sum(sqrt((finalDiff)^2))
  finalDiff[finalDiff < -1] <- -1
  finalDiff[finalDiff > 1] <- 1

  meshDist(atlas, distvec = finalDiff, rampcolors = c(colMinus, colPlus), from = -1, to = 1)

  rgl.snapshot(glue("~/postdoc/achondroplasiaDrug/individualScores/{unique(covariateTreated$Individual)[i]}_lastVisitScore_diag1.png"))

  par3d(userMatrix = diag2)
  rgl.snapshot(glue("~/postdoc/achondroplasiaDrug/individualScores/{unique(covariateTreated$Individual)[i]}_lastVisitScore_diag2.png"))


}

treatedIndex2 <- which(finalDrugCovariates$Group == "Placebo")
scoresTreated2 <- severityDrugs[treatedIndex2]
covariateTreated2 <- finalDrugCovariates[treatedIndex2,]
finalDiffVector2 <- rep(NA, length(unique(covariateTreated2$Individual)))

for(i in 1:length(unique(covariateTreated2$Individual))){

  syndromeEffect <- achondModel$coefficients[5 ,] + achondModel$coefficients[1,]

  interceptShape <- matrix(restoreShapes(achondModel$coefficients[1,], symmetricPCA$original$loadings[,1:nPCs], sampleMeanShape), ncol = 3, byrow = F)

  syndVector <- matrix(restoreShapes(syndromeEffect, symmetricPCA$original$loadings[,1:nPCs], sampleMeanShape), ncol = 3, byrow = F) - interceptShape

  indIndex <- which(covariateTreated2$Individual == unique(covariateTreated2$Individual)[i])

  firstVisit <- which.min(covariateTreated2$Visit[indIndex])
  lastVisit <- which.max(covariateTreated2$Visit[indIndex])
  firstVisitVerts <- restoreShapes(drugScores[indIndex,][firstVisit,], symmetricPCA$original$loadings[,1:nPCs], sampleMeanShape)
  lastVisitVerts <- restoreShapes(drugScores[indIndex,][lastVisit,], symmetricPCA$original$loadings[,1:nPCs], sampleMeanShape)

  personAligned <- rotmesh.onto(makeMesh(firstVisitVerts), firstVisitVerts, interceptShape, scale = T)

  personVector <- personAligned$yrot - interceptShape
  syndPersonVector <- lastVisitVerts - interceptShape
  # open3d()
  par3d(windowRect = c(0, 0, 800, 800), zoom = .75, userMatrix = diag1)
  # orient face, then save front.face <- par3d()$userMat
  firstProjection <- perVertexHeatmap(atlas, syndVector, personVector, type = "projection")
  rgl.snapshot(glue("~/postdoc/achondroplasiaDrug/individualScores/{unique(covariateTreated2$Individual)[i]}_firstVisitScore_diag1.png"))

  par3d(userMatrix = diag2)
  rgl.snapshot(glue("~/postdoc/achondroplasiaDrug/individualScores/{unique(covariateTreated2$Individual)[i]}_firstVisitScore_diag2.png"))

  #then plot for syndromic
  lastProjection <- perVertexHeatmap(atlas, syndVector, syndPersonVector, type = "projection")
  par3d(userMatrix = diag1)

  finalDiff <- lastProjection - firstProjection
  finalDiffVector2[i] <- sum(sqrt((finalDiff)^2))

  finalDiff[finalDiff < -1] <- -1
  finalDiff[finalDiff > 1] <- 1

  meshDist(atlas, distvec = finalDiff, rampcolors = c(colMinus, colPlus), from = -1, to = 1)

  rgl.snapshot(glue("~/postdoc/achondroplasiaDrug/individualScores/{unique(covariateTreated2$Individual)[i]}_lastVisitScore_diag1.png"))

  par3d(userMatrix = diag2)
  rgl.snapshot(glue("~/postdoc/achondroplasiaDrug/individualScores/{unique(covariateTreated2$Individual)[i]}_lastVisitScore_diag2.png"))

}

unique(covariateTreated$Individual)[sort(finalDiffVector, index.return = T, decreasing = T)$ix]
unique(covariateTreated2$Individual)[sort(finalDiffVector2, index.return = T, decreasing = T)$ix]


# statistical tests####
# multivariate paired t-test

# Differences
pcDiffs <- matrix(NA, ncol = ncol(drugScores), nrow = length(unique(treatedCovariates$Individual)))
pcFirst <- matrix(NA, ncol = ncol(drugScores), nrow = length(unique(treatedCovariates$Individual)))
pcLast <- matrix(NA, ncol = ncol(drugScores), nrow = length(unique(treatedCovariates$Individual)))
treatedIndex <- which(finalDrugCovariates$Group == "Placebo" | finalDrugCovariates$Group == "Drug")
treatedCovariates <- finalDrugCovariates[treatedIndex,]

for(i in 1:length(unique(treatedCovariates$Individual))){

  indIndex <- which(treatedCovariates$Individual == unique(treatedCovariates$Individual)[i])

  firstVisit <- which.min(treatedCovariates$Visit[indIndex])
  lastVisit <- which.max(treatedCovariates$Visit[indIndex])
  firstScan <- pcFirst[i,] <- drugScores[indIndex,][firstVisit,]
  lastScan <- pcLast[i,] <- drugScores[indIndex,][lastVisit,]

  pcDiffs[i,] <- lastScan - firstScan

}

# Hotellings
pcPaired <- MVTests::Mpaired(pcFirst[,1:13], pcLast[,1:13])

#no diff for full face
pcPaired$p.value

# same answer
nDims <- 13
ICSNP::HotellingsT2(pcDiffs[,1:nDims], mu = rep(0, nDims), test = "f")

# mixed model on score ####
View(drugCovariates)

# for placebo, find the oldest non-treated day or find the first visit with treatment and pick the visit before

drugScores
drugCovariates

drugCovariates$Visit[drugCovariates$Treated == 1]

# for each placebo ind, find min visit where treated and subtract 1

covariateTreated2

library(dplyr)

df <- data.frame(
  ID = finalDrugCovariates$Individual,
  visit = finalDrugCovariates$Visit,
  date = as.Date(finalDrugCovariates$DOV, tryFormats = "%d-%m-%Y"),
  treated = finalDrugCovariates$Treated
)

# mutate a count of days from the first treatment
daysAfterTreatment <- df %>%
  mutate(original_order = row_number()) %>%
  group_by(ID) %>%
  arrange(date) %>%
  mutate(
    first_treatment_date = if_else(any(treated == 1), min(date[treated == 1], na.rm = TRUE), NA),
    date_before_first_treatment = if_else(!is.na(first_treatment_date),
                                          lag(date, default = NA)[date == first_treatment_date],
                                          NA),
    days_from_before_first_treatment = if_else(!is.na(date_before_first_treatment),
                                               as.numeric(date - date_before_first_treatment),
                                               NA)
  ) %>%
  ungroup() %>%
  arrange(original_order)

View(daysAfterTreatment)

# merge to existing covariate data and model with scores
library(lme4)

finalDrugCovariates$treatmentLag <- daysAfterTreatment$days_from_before_first_treatment/365
treatedIndex <- which(finalDrugCovariates$Group == "Placebo" | finalDrugCovariates$Group == "Drug")
treatedCovariates <- finalDrugCovariates[treatedIndex,]
treatedCovariates$score <- severityDrugs[treatedIndex]

# colsToKeep <- c(3,6,16,18,29,30)
# treatedCovariates <- treatedCovariates[,colsToKeep]

mmAchondroplasia <- lmer(data = treatedCovariates, formula = score ~ Group + treatmentLag + Group:treatmentLag + (1 | Individual))
nullAchondroplasia <- lmer(data = treatedCovariates, formula = score ~ treatmentLag + (1 | Individual))
summary(mmAchondroplasia)
summary(nullAchondroplasia)

anova(mmAchondroplasia, nullAchondroplasia)


car::Anova(mmAchondroplasia)

# Figure 3 - Severity over time####
treatedCovariates$PhaseFactor <- treatedCovariates$Phase
treatedCovariates$PhaseFactor[treatedCovariates$Phase == 2] <- "Phase 2"
treatedCovariates$PhaseFactor[treatedCovariates$Phase == 3 & treatedCovariates$Group == "Placebo"] <- "Phase 3 - Placebo"
treatedCovariates$PhaseFactor[treatedCovariates$Phase == 3 & treatedCovariates$Group == "Drug"] <- "Phase 3 - Drug"


pdf("~/postdoc/achondroplasiaDrug/scoreOverTime.pdf", height = 5, width = 10)
ggplot(treatedCovariates, aes(x = treatmentLag/365, y = score, group = Phase)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "lightgrey") +
  geom_line(aes(group = Individual)) +
  geom_point() +
  xlim(c(-1.5, 8)) +
  # scale_color_manual(values = phaseColors) +
  labs(x = "Years after first treatment", y = "Achondroplasia score") +
  theme_minimal(base_line_size = 0) +
  facet_wrap(~PhaseFactor) #+
  # theme(legend.position = "none")
dev.off()
