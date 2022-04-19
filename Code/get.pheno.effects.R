#This function gets the mean phenotype effects for
#two markers and their combined effect using three
#different patterns of genotype coding
#To use a binning other than the ones offered by 
#additive, dominant, and recessive, use geno.bins
#to define custom binning


get.pheno.effects <- function(data.obj, geno.obj, marker1.name, marker2.name, pheno.name, covar = NULL, 
scan.what = c("normalized.traits", "raw.traits", "eigentraits"), 
collapsed.net = FALSE, geno.coding = c("Additive", "Dominant", "Recessive"), geno.bins = NULL){

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
	#=========================================================


	if(collapsed.net){
		m1.markers <- data.obj$linkage_blocks_collapsed[[which(names(data.obj$linkage_blocks_collapsed) == marker1.name)]]
		split.m1 <- strsplit(m1.markers, "_")
		m1.name <- sapply(split.m1, function(x) x[1])
		m1.allele <- sapply(split.m1, function(x) x[2])

		m2.markers <- data.obj$linkage_blocks_collapsed[[which(names(data.obj$linkage_blocks_collapsed) == marker2.name)]]
		split.m2 <- strsplit(m1.markers, "_")
		m2.name <- sapply(split.m2, function(x) x[1])
		m2.allele <- sapply(split.m2, function(x) x[2])
		
		m1.marker.ind <- match(m1.name, dimnames(geno)[[3]])
		m1.allele.ind <- match(m1.allele, dimnames(geno)[[2]])
		
		m2.marker.ind <- match(m2.name, dimnames(geno)[[3]])
		m2.allele.ind <- match(m2.allele, dimnames(geno)[[2]])

		all.int <- data.obj$full_net[m1.marker.ind, m2.marker.ind, drop=FALSE]
		max.int <- which(abs(all.int) == max(abs(all.int)), arr.ind = TRUE)

		m1.geno <- geno[not.na.locale,m1.allele.ind[max.int[1,1]], m1.marker.ind[max.int[1,1]]]
		m2.geno <- geno[not.na.locale,m2.allele.ind[max.int[1,1]], m2.marker.ind[max.int[1,2]]]

		m1.name <- 	m1.markers[max.int[1,1]]
		m2.name <- 	m2.markers[max.int[1,2]]
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
			}
			
		if(geno.coding == "Recessive"){
			geno.bins <- c(0, 1)
			}
	
		if(geno.coding == "Additive"){
			geno.bins <- c(0, 0.5, 1)
		}
	}

	binned.m1 <- bin.vector(m1.geno, geno.bins)
	binned.m2 <- bin.vector(m2.geno, geno.bins)

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