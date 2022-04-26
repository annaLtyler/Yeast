#This function gets the mean phenotype effects for
#two markers and their combined effect using three
#different patterns of genotype coding
#To use a binning other than the ones offered by 
#additive, dominant, and recessive, use geno.bins
#to define custom binning


get.pheno.effects <- function(data.obj, geno.obj, marker1.name, 
	marker2.name, pheno.name, covar = NULL, 
	scan.what = c("normalized.traits", "raw.traits", "eigentraits"), 
	collapsed.net = FALSE, geno.coding = c("Additive", "Dominant", "Recessive"), 
	geno.bins = NULL, present.thresh = 0.3){

	pheno <- get_pheno(data.obj, scan_what = scan.what, covar = covar)
	geno <- get_geno(data.obj, geno.obj)

	covar.info <- get_covar(data.obj)
		if(!is.null(covar.info$covar_table)){
			not.na.locale <- which(!is.na(rowSums(covar.info$covar_table)))
			}else{
			not.na.locale <- 1:nrow(pheno)	
			}

	#=========================================================
	#internal functions
	#=========================================================
	get.marker.name <- function(marker.name){
		if(length(grep("Chr", marker.name)) > 0){
			marker.locale <- which(names(data.obj$linkage_blocks_full) == marker.name)
			marker.label <- data.obj$linkage_blocks_full[[marker.locale]]
			return(marker.label)
			}else{
			if(is.null(allele)){
				full.name <- marker.name
				}else{
				full.name <- paste0(marker.name, "_", allele)
				}
			return(marker.name)
			}
		}

	collapsed_to_full <- function(collapsed_marker_name){
		block.locale <- which(names(data.obj$linkage_blocks_collapsed) == collapsed_marker_name)
		marker.names <- data.obj$linkage_blocks_collapsed[[block.locale]]
		marker.locale <- match(marker.names, unlist(data.obj$linkage_blocks_full))
		collapsed_marker_names <- names(data.obj$linkage_blocks_full)[marker.locale]
		return(collapsed_marker_names)
	}
	#=========================================================


	if(collapsed.net){
		#if we are interested in the collapsed network, 
		#find the strongest interaction among the marker 
		#pairs making up the linkage blocks.
		
		m1.blocks <- collapsed_to_full(marker1.name)
		m2.blocks <- collapsed_to_full(marker2.name)
		all.int <- data.obj$full_net[m1.blocks, m2.blocks, drop=FALSE]
		max.int <- which(abs(all.int) == max(abs(all.int)), arr.ind = TRUE)
		m1.name <- rownames(all.int)[max.int[1,1]]
		m2.name <- colnames(all.int)[max.int[1,2]]

		max.marker1 <- get.marker.name(m1.name)
		max.marker2 <- get.marker.name(m2.name)

		split.m1 <- strsplit(max.marker1, "_")
		m1.name <- sapply(split.m1, function(x) x[1])
		m1.allele <- sapply(split.m1, function(x) x[2])

		split.m2 <- strsplit(max.marker2, "_")
		m2.name <- sapply(split.m2, function(x) x[1])
		m2.allele <- sapply(split.m2, function(x) x[2])
		
		m1.marker.ind <- match(m1.name, dimnames(geno)[[3]])
		m1.allele.ind <- match(m1.allele, dimnames(geno)[[2]])
		
		m2.marker.ind <- match(m2.name, dimnames(geno)[[3]])
		m2.allele.ind <- match(m2.allele, dimnames(geno)[[2]])

		m1.geno <- geno[not.na.locale,m1.allele.ind[max.int[1,1]], m1.marker.ind[max.int[1,1]]]
		m2.geno <- geno[not.na.locale,m2.allele.ind[max.int[1,1]], m2.marker.ind[max.int[1,2]]]

		}else{
			
		m1.marker <- get.marker.name(marker1.name)
		split.m1 <- strsplit(m1.marker, "_")
		m1.name <- sapply(split.m1, function(x) x[1])
		m1.allele <- sapply(split.m1, function(x) x[2])

		m2.marker <- get.marker.name(marker2.name)
		split.m2 <- strsplit(m2.marker, "_")
		m2.name <- sapply(split.m2, function(x) x[1])
		m2.allele <- sapply(split.m2, function(x) x[2])

		m1.geno <- geno[not.na.locale,which(dimnames(geno)[[2]] == m1.allele), which(dimnames(geno)[[3]] == m1.name)]
		m2.geno <- geno[not.na.locale,which(dimnames(geno)[[2]] == m2.allele), which(dimnames(geno)[[3]] == m2.name)]
		}

	
	if(is.null(geno.bins)){
		if(geno.coding == "Dominant"){
			geno.bins <- c(0, 1)
			binned.m1  <- m1.geno
			binned.m1[which(binned.m1 >= present.thresh)] <- 1
			binned.m1[which(binned.m1 < present.thresh)] <- 0
			
			binned.m2 <- m2.geno
			binned.m2[which(binned.m2 >= present.thresh)] <- 1
			binned.m2[which(binned.m2 < present.thresh)] <- 0			
			}
			
		if(geno.coding == "Recessive"){
			geno.bins <- c(0, 1)
			binned.m1  <- m1.geno
			binned.m1[which(binned.m1 >= (present.thresh*2))] <- 1
			binned.m1[which(binned.m1 < (present.thresh*2))] <- 0
			
			binned.m2 <- m2.geno
			binned.m2[which(binned.m2 >= (present.thresh*2))] <- 1
			binned.m2[which(binned.m2 < (present.thresh*2))] <- 0					
			}
	
		if(geno.coding == "Additive"){
			geno.bins <- c(0, 0.5, 1)
		}
	}

	pheno.locale <- which(colnames(pheno) == pheno.name)
	pheno.vals <- pheno[,pheno.locale]

	m1.min.locale <- which(binned.m1 == min(geno.bins))
	m1.max.locale <- which(binned.m1 == max(geno.bins))
	m2.min.locale <- which(binned.m2 == min(geno.bins))
	m2.max.locale <- which(binned.m2 == max(geno.bins))

	no.m1.no.m2 <- intersect(m1.min.locale, m2.min.locale)
	baseline.val <- mean(pheno.vals[no.m1.no.m2], na.rm = TRUE)

	just.m1 <- intersect(m1.max.locale, m2.min.locale)
	m1.coef <- mean(pheno.vals[just.m1], na.rm = TRUE) - baseline.val

	just.m2 <- intersect(m1.min.locale, m2.max.locale)
	m2.coef <- mean(pheno.vals[just.m2], na.rm = TRUE) - baseline.val
	
	m1.m2 <- intersect(m1.max.locale, m2.max.locale)
	double.coef <- mean(pheno.vals[m1.m2], na.rm = TRUE) - baseline.val
		
	result <- c("m1.name" = m1.name, "m2.name" = m2.name, "m1.effect" = m1.coef, "m2.effect" = m2.coef, "double.effect" = double.coef) 
	return(result)
	}