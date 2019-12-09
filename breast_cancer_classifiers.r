# R-3.5.2
#####################################################################
# Classification Van't Veer et al., Nature 2002 70-gene 'MammaPrint'
#####################################################################
`MammaPrint` <- function(DATA, ID, CM=TRUE){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# CM, adjust data by median center genes, default = TRUE
#
	#### Matching DATA w/ model
	intsectSymbol <- setdiff(intersect(ID, .mammaprint_model$Symbol), which(is.na(.mammaprint_model$Symbol)))
	if (length(intsectSymbol)==0){ stop("\nno match\n") }
	dat <- DATA[match(intsectSymbol, ID),]
	ref <- .mammaprint_model[match(intsectSymbol, .mammaprint_model$Symbol),]
	if (CM) { dat <- dat - apply(dat,1,median,na.rm=T) }
	#### Classification
	Classif_MammaPrint <- list(Score70g=apply(dat,2, function(x){cor.test(x, ref$Profil.Good.70g,na.rm=T)$estimate}))
	Classif_MammaPrint$Grp70g <- cut(Classif_MammaPrint$Score70g, c(-Inf, 0.4, Inf), c("good","poor"))   # cut-off : Fig2 b/c Van't Veer et al., Nature 2002
	#### output
	return(data.frame(Classif_MammaPrint))
}
##########################################################################
# Classification Agendia/BluePrint, Krijgsman et al., BCRT 2012 'MSP-80g'
##########################################################################
`BluePrint` <- function(DATA, ID, CM=TRUE){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# CM, adjust data by median center genes, default = TRUE
#
	#### Matching DATA w/ model
	intsectSymbol <- setdiff(intersect(ID, .blueprint_model$ref.sel$Symbol), which(is.na(.blueprint_model$ref.sel$Symbol)))
	if (length(intsectSymbol)==0){ stop("\nno match\n") }
	dat <- DATA[match(intsectSymbol, ID),]
	if (CM) { dat <- dat - apply(dat,1,median,na.rm=T) }
	ref <- .blueprint_model$ref.sel[match(intsectSymbol, .blueprint_model$ref.sel$Symbol),]
	# Scores, centroid correlations
	Score_mat <- data.frame(list(
		Basal = apply(dat,2,function(x){cor.test(x,ref$Basal,na.rm=T)$estimate[[1]]}),
		HER2 = apply(dat,2,function(x){cor.test(x,ref$HER2,na.rm=T)$estimate[[1]]}),
		Luminal = apply(dat,2,function(x){cor.test(x,ref$Luminal,na.rm=T)$estimate[[1]]}) ))
	######## def Grp 'part 1'
	Score_mat$Grp_80GS 	<- unlist(apply(Score_mat,1,function(x){colnames(Score_mat)[which.max(x)[1]]}))
	######## def Grp 'part 2, Luminal discretization' w/ Mammaprint/70-gene risk classifier
	# Mammaprint/70g
	Prog70gMammaPrint <- MammaPrint(DATA, ID, CM=TRUE)
    	Score_mat$Mammaprint70g <- as.character(Prog70gMammaPrint$Grp70g)
	Score_mat$grp <- as.character(Score_mat$Grp_80GS)
	Score_mat$grp[which(Score_mat$grp %in% "Luminal" & Score_mat$Mammaprint70g %in% "good")] <- "LuminalA"
	Score_mat$grp[which(Score_mat$grp %in% "Luminal" & Score_mat$Mammaprint70g %in% "poor")] <- "LuminalB"
	# output
	return(Score_mat)
}
#################################################
# Classification Parker et al., JCO 2009 'PAM50'
#################################################
`PAM50` <- function(DATA, ID, CM=TRUE){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# CM, adjust data by median center genes, default = TRUE
#
	#### Matching DATA w/ model
	intsectSymbol <- setdiff(intersect(ID, .pam50_model$Symbol), which(is.na(.pam50_model$Symbol)))
	if (length(intsectSymbol)==0){ stop("\nno match\n") }
	dat <- DATA[match(intsectSymbol, ID),]
	ref <- .pam50_model[match(intsectSymbol, .pam50_model$Symbol),]
	if (CM) { dat <- dat - apply(dat,1,median,na.rm=T) }
	#### Classification, nearest centroid correlations
	selref <- which(colnames(ref) %in% c("Basal","Her2","LumA","LumB","Normal"))
	Score_mat <- data.frame(apply(ref[,selref],2, function(x){
			unlist(apply(dat,2,function(y){ cor.test(x,y,na.rm=T)$estimate[[1]] }))
		}))
	Score_mat$PAM50 <- unlist(apply(Score_mat,1,function(x){colnames(Score_mat)[which.max(x)[1]]}))
	# output
	return(Score_mat)
}
#################################################
# Classification Hess et al., JCO 2006 'DLDA-30'
#################################################
`DLDA30` <- function(DATA, ID, CM=TRUE){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# CM, adjust data by median center genes, default = TRUE
#
	require(sfsmisc)
	#### Matching DATA w/ model
	intsectSymbol <- setdiff(intersect(ID, .dlda30_model$ref$Symbol), which(is.na(.dlda30_model$ref$Symbol)))
	if (length(intsectSymbol)==0){ stop("\nno match\n") }
	datsel <- DATA[match(intsectSymbol, ID),]
	modelsel <- .dlda30_model$learn$x[match(intsectSymbol, .dlda30_model$ref$Symbol),]
	if (CM) {
		datsel <- datsel - apply(datsel,1,median,na.rm=T)
		modelsel <- modelsel - apply(modelsel,1,median,na.rm=T)
	 }
	refsel <- .dlda30_model$ref[match(intsectSymbol, .dlda30_model$ref$Symbol),]
	# ctrl/remove gene w/ NA values
	selNA <- apply(datsel, 1, function(x){sum(is.na(x))})
	datsel <- datsel[which(selNA==0),]
	modelsel <- modelsel[which(selNA==0),]
	# DLDA classification
	# Learning model
	GRP_Model <- ifelse(.dlda30_model$learn$y=="RD",0,1)
	DLDA_Model <- dDA(t(modelsel), GRP_Model, pool=FALSE)
	# Classification
	PredTest <- list(Grp=diagDA(t(modelsel), GRP_Model, t(datsel), pool=FALSE))
	PredTest$GrpTxt <- ifelse(PredTest$Grp==1, "pCR_like", ifelse(PredTest$Grp==0, "RD_like",NA))
	# output
	return(data.frame(PredTest))
}
###################################################################
# Classification Guerrero-Zotana et al., CCR 2018, 24 genes 'E2F4'
###################################################################
`E2F4` <- function(DATA, ID, CM=TRUE){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# CM, adjust data by median center genes, default = TRUE
#
	#### Matching DATA w/ model
	intsectSymbol <- setdiff(intersect(ID, .E2F4_model$ref$Symbol), which(is.na(.E2F4_model$ref$Symbol)))
	if (length(intsectSymbol)==0){ stop("\nno match\n") }
	dat <- DATA[match(intsectSymbol, ID),]
	if (CM) { dat <- dat - apply(dat,1,median,na.rm=T) }
	ref <- .E2F4_model$ref[match(intsectSymbol, .E2F4_model$ref$Symbol),]
	# Metagene
	Score_mat <- list(score_raw=apply(dat,2,median,na.rm=T))
	Score_mat$score_scaled <- (Score_mat$score_raw - mean(Score_mat$score_raw,na.rm=T))/sd(Score_mat$score_raw,na.rm=T)
	Score_mat$grp <- cut(Score_mat$score_scaled, c(-Inf, 0, Inf), c("low","high"))
	# output
	return(data.frame(Score_mat))
}
############################################################################
# Classification Bindea et al., Immunity 2013, 28 Immunes types, 'Immunome'
############################################################################
`Immunome28` <- function(DATA, ID, CM=TRUE){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# CM, adjust data by median center genes, default = TRUE
#
	#### Matching DATA w/ model
	intsectSymbol <- setdiff(intersect(ID, .immunome28_model$ref$Symbol), which(is.na(.immunome28_model$ref$Symbol)))
	if (length(intsectSymbol)==0){ stop("\nno match\n") }
	dat <- DATA[match(intsectSymbol, ID),]
	if (CM) { dat <- dat - apply(dat,1,median,na.rm=T) }
	ref <- .immunome28_model$ref[match(intsectSymbol, .immunome28_model$ref$Symbol),]
	# Metagene
    	Metagene_sellist <- tapply(seq_along(ref$Symbol), ref$CellType, function(x){x})
	Metagene_list <- lapply(Metagene_sellist, function(x){apply(data.frame(dat)[x,],2,mean,na.rm=T)})
	# ctrl Metagene
	MetageneID <- unique(as.character(.immunome28_model$ref$CellType))
	if (length(Metagene_list)<length(MetageneID)){
     	selMiss <- which(!MetageneID %in% names(Metagene_list))
		N <- length(Metagene_list)
		for (i in seq_along(selMiss)){
          	MetageneID[[N+i]] <- rep(NA, length(MetageneID[[N]]))
		}
		names(MetageneID)[(N+1):(N+length(selMiss))] <- MetageneID[selMiss]
	}
  	#output
  	return(data.frame(Metagene_list))
}
######################################################################
# Classification Bertucci et al., BJC 2018, 20 genes, 'ICR'
######################################################################
`ICR` <- function(DATA, ID, CM=TRUE){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# CM, adjust data by median center genes, default = TRUE
#
	#### Matching DATA w/ model
	intsectSymbol <- setdiff(intersect(ID, .icr_model$ref$Symbol), which(is.na(.icr_model$ref$Symbol)))
	if (length(intsectSymbol)==0){ stop("\nno match\n") }
	dat <- DATA[match(intsectSymbol, ID),]
	if (CM) { dat <- dat - apply(dat,1,median,na.rm=T) }
	ref <- .icr_model$ref[match(intsectSymbol, .icr_model$ref$Symbol),]
	# NMF classification
	require(ConsensusClusterPlus)
     ConsClust <- ConsensusClusterPlus(dat, distance="euclidean", innerLinkage="ward.D2", finalLinkage="complete",    # params ref papier Hendrickx et al.
                                    maxK=4, reps=500, seed=123)[[4]]
  	# ICR groups attribution
	Metagene_ctrl_score_groups <- tapply(apply(dat,2,mean,na.rm=T), ConsClust$consensusClass, median, na.rm=T)
	order_Ctrl <- order(unlist(Metagene_ctrl_score_groups))
	ICR <- rep(NA, length(ConsClust$consensusClass))
	for (i in seq_along(order_Ctrl)){
		ICR[which(ConsClust$consensusClass %in% order_Ctrl[i])] <- paste("ICR", i, sep="_")
	}
	ICR <- data.frame(ICR)
	ICR$ICR_2K <- ifelse(ICR$ICR %in% "ICR_4", "ICR_4", ifelse(ICR$ICR %in% c("",NA), NA, "ICR_1-3"))
	# output
	return(data.frame(ICR))
}
#####################################################################
# Classification Creighton et al., PNAS 2009, 493 pbs, 'CD44+/CD24-'
#####################################################################
`CD44CD24` <- function(DATA, ID, CM=TRUE){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# CM, adjust data by median center genes, default = TRUE
#
	#### Matching DATA w/ model
	intsectSymbol <- setdiff(intersect(ID, .cd44cd24_model$ref$Symbol), which(is.na(.cd44cd24_model$ref$Symbol)))
	if (length(intsectSymbol)==0){ stop("\nno match\n") }
	dat <- DATA[match(intsectSymbol, ID),]
	if (CM) { dat <- dat - apply(dat,1,median,na.rm=T) }
	ref <- .cd44cd24_model$ref[match(intsectSymbol, .cd44cd24_model$ref$Symbol),]
	# Classification, Pearson Correlation
	Score_mat <- list(score=apply(dat,2, function(x){cor.test(x, ref$Pattern, na.rm=T, method="pearson")$estimate[[1]]}))
	Score_mat$grp <- cut(Score_mat$score, c(-Inf,0,Inf), c("noCD44+/CD24-_MS-like", "CD44+/CD24-_MS-like"))		# natural r0 cut-off
	# output
	return(data.frame(Score_mat))
}
###################################################################################################
# Classification Lim et al., Nature Medicine 2009, 3+1 of 'Breast Epithelial Hierachy' populations
###################################################################################################
`breastHierarch` <- function(DATA, ID, CM=TRUE){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# CM, adjust data by median center genes, default = TRUE
#
	#### Matching DATA w/ model
	intsectSymbol_list <- lapply(.breastHierarch_model, function(x){setdiff(intersect(ID, x$Symbol), which(is.na(x$Symbol))) })
	dat_list <- ref_list <- list()
	for (i in seq_along(intsectSymbol_list)){
     	if (length(intsectSymbol_list[[i]]) < 1) { stop("\nno match\n") }
		dat_list[[i]] <- DATA[match(intsectSymbol_list[[i]], ID),]
		ref_list[[i]] <- .breastHierarch_model[[i]][match(intsectSymbol_list[[i]],.breastHierarch_model[[i]]$Symbol),]
	}
	names(dat_list) <- names(ref_list) <- names(.breastHierarch_model)
	if (CM) {for(i in seq_along(dat_list)){ dat_list[[i]] <- dat_list[[i]]-apply(dat_list[[i]],1,median,na.rm=T) }}
	# Classification score
	score_mat <- list()
	for(i in seq_along(dat_list)){
     	score_mat[[i]] <- apply(dat_list[[i]]*ref_list[[i]]$Average,2,mean,na.rm=T)
	}
	names(score_mat) <- names(dat_list)
	# output
	return(data.frame(score_mat))
}
######################################################################
# Classification Filipits et al., CCR 2011, 8+3 genes, 'EndoPredict'
######################################################################
`EndoPredict` <- function(DATA, ID, ER=NULL, HER2=NULL, Tc=NULL, Nc=NULL){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# ER, vector of samples estrogen receptor statut (bin, 1|0)
# HER2, vector of samples HER2 statut (bin, 1|0)
# Tc : clinical tumor size. code= 1: <=1 cm, 2: >1 to <=2 cm, 3: >2 to <=5 cm, and 4: >5 cm
# Nc : clinical nodal status. code= 1: negative, 2: 1–3 positive nodes, 3: 4–10 positive nodes, and 4: >10 positive nodes
# CM, adjust data by median center genes, default = TRUE
#
	#### Matching DATA w/ model
	if (is.null(ER)|is.null(HER2)){stop("\nER & HER2 status required\n")}
	ER <- as.numeric(as.character(ER))
	HER2 <- as.numeric(as.character(HER2))
	intsectSymbol <- setdiff(intersect(ID, .endopred_model$ref$Symbol), which(is.na(.endopred_model$ref$Symbol)))
	if (length(intsectSymbol)==0){ stop("\nno match\n") }
	dat <- DATA[match(intsectSymbol, ID),]
	ref <- .endopred_model$ref[match(intsectSymbol, .endopred_model$ref$Symbol),]
	selERposHER2neg <- which(ER == 1 & HER2 == 0)
	if (length(selERposHER2neg) < 2){stop("\nER+/HER2- samples missing\n") }
	dat <- (((dat - apply(dat[,selERposHER2neg],1,mean,na.rm=T))/apply(dat[,selERposHER2neg],1,sd,na.rm=T))*ref$SD_COI)+ref$Mean_COI # data normalized / ref data on ER+/HER2- samples
	# Classification score
	sel_ctrl <- which(ref$type %in% "norm")
	dat_score <- (dat[-sel_ctrl,] - ref$b[-sel_ctrl]) / ref$a[-sel_ctrl]
	ref_score <- ref[-sel_ctrl,]
	Su <- (apply(dat_score,2,function(x){ sum((x*ref_score$coef_EndoPredict),na.rm=T)-2.63 })*1.5)+18.95
	Su_rescaled <- ifelse(Su<0,0,ifelse(Su>15,15,Su))
	EP <- list(Su=Su_rescaled)
	EP$grp <- cut(Su_rescaled, c(-Inf, 4.9999, Inf), c("low-risk", "high-risk"))
	if (!is.null(Tc) & !is.null(Nc)){
     	Tc <- as.numeric(as.character(Tc))
        	Nc <- as.numeric(as.character(Nc))
        	EP$clin_score <- as.numeric(EP$Su*0.28 + 0.35*Tc + 0.64*Nc)
        	EP$clin_grp <- cut(EP$clin_score, c(-Inf, 3.3, Inf), c("low-risk", "high-risk"))
    	}else{
		EP$clin_score <- EP$clin_grp <- rep(NA, length(Su))
	}
	# output
	return(data.frame(EP))
}
#####################################################################
# Classification Paik et al., NEJM 2004, 21 genes, 'Recurrence Score'
#####################################################################
`RecurrenceScore` <- function(DATA, ID, CM=TRUE){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# CM, adjust data by median center genes, default = TRUE
#
	#### Matching DATA w/ model
	intsectSymbol <- setdiff(intersect(ID, .recurrenceScore_model$ref$Symbol), which(is.na(.recurrenceScore_model$ref$Symbol)))
	if (length(intsectSymbol)==0){ stop("\nno match\n") }
	dat <- DATA[match(intsectSymbol, ID),]
	if (CM) { dat <- dat - apply(dat,1,median,na.rm=T) }
	ref <- .recurrenceScore_model$ref[match(intsectSymbol, .recurrenceScore_model$ref$Symbol),]
	# Classification score
	# ref center
	dat.test <- t(dat [ref$type != "Ref",])
	dat.test <- t((dat.test) - apply( t(dat[ref$type == "Ref",]),1,mean, na.rm=T))
	ref.test <- ref[ref$type != "Ref",]
	# norm 0-15
	dat.test.015 <- dat.test
	for (i in 1:ncol(dat.test)){
	    temp0 <- (as.numeric(dat.test[,i])-min(as.numeric(dat.test[,i])))
	    coef.max <- 15/max(temp0)
	    dat.test.015[,i] <- coef.max * temp0
	}
	# weight coef Paik
	dat.test.015.coef <- dat.test.015 * ref.test$coef.Paik
	## markers
	# GRB7
	if (sum(ref.test$type == "HER2") >1){
		GRB7 <- as.matrix(apply(dat.test.015.coef[ref.test$type=="HER2",],2,sum))
	}else{
		GRB7 <- (as.matrix(dat.test.015.coef[ref.test$type=="HER2",]))
	}
	GRB7[GRB7<8]<-8
	# ER
	ER <- as.matrix(apply(dat.test.015.coef[ref.test$type=="Estrogen",],2,sum))
	# Prolif
	Prolif <- as.matrix(apply(dat.test.015.coef[ref.test$type=="Prolif",],2,sum))
	Prolif[Prolif<6.5]<-6.5
	# Invasion
	if (sum(ref.test$type == "Invasion") >1){
		Invasion <- as.matrix(apply(dat.test.015.coef[ref.test$type=="Invasion",],2,sum))
	}else{
		Invasion <- (as.matrix(dat.test.015.coef[ref.test$type=="Invasion",]))
	}
	# Neutre & intermediate score
	selNeutre <- which(ref.test$type=="neutre")
	if (length(selNeutre) == 0){
		RS.int <- cbind(GRB7*0.47, ER*(-0.34), Prolif*1.04, Invasion*0.1)
	}else{
		if (length(selNeutre) == 1){
	     	RS.int<- cbind(GRB7*0.47, ER*(-0.34), Prolif*1.04, Invasion*0.1, (dat.test.015.coef[ref.test$type=="neutre",]))
		}else{
	      	RS.int<- cbind(GRB7*0.47, ER*(-0.34), Prolif*1.04, Invasion*0.1, t(dat.test.015.coef[ref.test$type=="neutre",]))
	  	}
	}
	# RSu score
	RSu.raw <- apply(RS.int,1, sum, na.rm=T)
	# RS score
	RS <- RSu.raw
	RS[RS>100] <- 100
	RS[RS<=0] <- 0
	RS[RSu.raw <= 100 & RSu.raw >0] <- 20*(RS[RSu.raw <= 100 & RSu.raw >0]-6.7)
	# RS scaled score
	mat_score <- list(score= (RS - min(RS,na.rm=T)) * (100 / (max(RS,na.rm=T) - min(RS,na.rm=T))))
	mat_score$grp <- cut(mat_score$score, c(-Inf,18,30.9999,Inf), c("low", "intermediate", "high"))
	# output
	return(data.frame(mat_score))
}
############################################################################################################
# Classification Turner et al., JCO 2019, E2Fregulon metagene (from 'E2F targets' Hallmark MySig/Base v.6.2)
############################################################################################################
`E2Fregulon` <- function(DATA, ID, CM=TRUE){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# CM, adjust data by median center genes, default = TRUE
#
	#### Matching DATA w/ model
	intsectSymbol <- setdiff(intersect(ID, .E2Fregulon$ref$Symbol), which(is.na(.E2Fregulon$ref$Symbol)))
	if (length(intsectSymbol)==0){ stop("\nno match\n") }
	dat <- DATA[match(intsectSymbol, ID),]
	if (CM) { dat <- (dat - apply(dat,1,median,na.rm=T)) }
	ref <- .E2Fregulon$ref[match(intsectSymbol, .E2Fregulon$ref$Symbol),]
	# Metagene
	Score_mat <- list(Zscore=apply(dat,2,mean,na.rm=T)/apply(dat,2,sd,na.rm=T))
	# output
	return(data.frame(Score_mat))
}
############################################################################################################
# Classification Malorni et al., OncoTarget 2016, E2F Z-score metagene (87 genes)
############################################################################################################
`RBsig` <- function(DATA, ID, CM=TRUE){
# DATA, data matrix (log2-transformed)
# ID, vector of Symbol
# CM, adjust data by median center genes, default = TRUE
#
	#### Matching DATA w/ model
	intsectSymbol <- setdiff(intersect(ID, .RBsig$ref$Symbol), which(is.na(.RBsig$ref$Symbol)))
	if (length(intsectSymbol)==0){ stop("\nno match\n") }
	dat <- DATA[match(intsectSymbol, ID),]
	if (CM) { dat <- dat - apply(dat,1,median,na.rm=T) }
	ref <- .RBsig$ref[match(intsectSymbol, .RBsig$ref$Symbol),]
	# Metagene
	Score_mat <- list(score_raw=apply(dat,2,median,na.rm=T))
	Score_mat$Zscore <- (Score_mat$score_raw - mean(Score_mat$score_raw,na.rm=T))/sd(Score_mat$score_raw,na.rm=T)
	# output
	return(data.frame(Score_mat))
}
