plot_clinical_effects_query <- function(data_obj, geno_obj, 
    pheno_type = c("pheno", "norm_pheno", "ET"), p_or_q = 0.05, 
    scale.coord = 1, verbose = FALSE){

    query_genotype <- data_obj$query_genotype
    trait_cols <- categorical_pal(8)
    pheno_type <- pheno_type[1]

    if(pheno_type == "pheno"){
        pheno <- data_obj$pheno
    }
    if(pheno_type == "norm_pheno"){
        pheno <- apply(data_obj$pheno, 2, rankZ)
    }
    if(pheno_type == "ET"){
        pheno <- data_obj$ET
    }

    common.ind <- Reduce("intersect", list(rownames(pheno), rownames(geno_obj), rownames(query_genotype)))
    common.pheno.locale <- match(common.ind, rownames(pheno))
    common.geno.locale <- match(common.ind, rownames(geno_obj))
    common.query.locale <- match(common.ind, rownames(query_genotype))
    matched.pheno <- pheno[common.pheno.locale,]
    matched.geno <- geno_obj[common.geno.locale,,]
    matched.query <- query_genotype[common.query.locale,]

    var_int <- write_variant_influences(data_obj, p_or_q = p_or_q, 
        include_main_effects = FALSE, mark_covar = FALSE, 
        write_file = FALSE)

    #each interaction in query_cape has the query genotype as
    #the source or the target.
    #plot these types of interactions separately.
    source.chr <- as.numeric(var_int[,2])
    source.pos <- as.numeric(var_int[,3])
    target.chr <- as.numeric(var_int[,5])
    target.pos <- as.numeric(var_int[,6])

    u_chr <- sort(unique(c(source.chr, target.chr)))
    u_chr <- u_chr[-1] #take off the 0 chromosome. this is the query marker
        
    get_effect <- function(phenoV, geno1, geno2){
        geno_pairs <- pair.matrix(c(0,1), self.pairs = TRUE, ordered = TRUE)
        pheno.groups <- apply(geno_pairs, 1, function(x) phenoV[intersect(which(geno1 == x[1]), which(geno2 == x[2]))])
        #boxplot(pheno.groups);abline(h = 0)
        group.means <- sapply(pheno.groups, function(x) mean(x, na.rm = TRUE))
        #group.se <- sapply(pheno.groups, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
        ref.effect <- group.means[1]
        centered.means <- group.means - ref.effect
        source.effect <- centered.means[4]
        target.effect <- centered.means[2]
        actual.effect <- centered.means[3]

        add.pred <- source.effect + target.effect

        if(is.finite(add.pred)){
			if(add.pred < 0){
				int.expect <- "expect.negative"
				}else{
				int.expect <- "expect.positive"
				}
			
			if(abs(actual.effect) < abs(add.pred)){
				int.effect <- "alleviating"
				}else{
				int.effect <- "aggravating"
				}
			}else{
				int.expect <- "none"
				int.effect <- "none"
				}
		
		pheno.result <- c("main1" = source.effect, "main2" = target.effect, 
        "additive" = add.pred, "actual" = actual.effect, "additive.direction" = int.expect, 
        "effect" = int.effect)
		return(pheno.result)
    }

    source.effects <- target.effects <- vector(mode = "list", length = length(u_chr))
    names(source.effects) <- names(target.effects) <- u_chr

    for(ch in u_chr){
        if(verbose){report.progress(ch, length(u_chr))}
        target.chr.locale <- which(target.chr == ch)
        source.chr.locale <- which(source.chr == ch)

        target.markers <- var_int[target.chr.locale,"Target"]
        source.markers <- var_int[source.chr.locale,"Source"]
        
        split.target <- strsplit(target.markers, "_")
        target.marker.names <- sapply(split.target, function(x) x[1])
        target.allele.names <- sapply(split.target, function(x) x[2])

        split.source <- strsplit(source.markers, "_")
        source.marker.names <- sapply(split.source, function(x) x[1])
        source.allele.names <- sapply(split.source, function(x) x[2])

        source.allele[[ch]] <- source.allele.names
        target.allele[[ch]] <- target.allele.names

        target.geno <- sapply(1:length(target.marker.names), function(x) matched.geno[,target.allele.names[x], target.marker.names[x]])
        colnames(target.geno) <- target.markers
        source.geno <- sapply(1:length(source.marker.names), function(x) matched.geno[,source.allele.names[x], source.marker.names[x]])
        colnames(source.geno) <- source.markers
        
        target.int <- lapply(1:ncol(matched.pheno), function(y) t(apply(target.geno, 2, function(x) get_effect(matched.pheno[,y], matched.query, x))))
        source.int <- lapply(1:ncol(matched.pheno), function(y) t(apply(source.geno, 2, function(x) get_effect(matched.pheno[,y], x, matched.query))))

        target.effects[[ch]] <- target.int
        source.effects[[ch]] <- source.int

    }

    all.target.effects <- Reduce("rbind", target.effects)
    all.source.effects <- lapply(1:ncol(pheno), function(y) Reduce("rbind", sapply(source.effects, function(x) x[[y]])))

    classify_motif <- function(marker.effects){
        if(sign(marker.effects[1]) == sign(marker.effects[2])){
            main <- "coherent"
        }else{
            main <- "incoherent"
        }
    }

    for(p in 1:length(all.source.effects)){
        all.num.effects <- apply(all.source.effects[[p]][,1:4], 2, as.numeric)
        #dim(all.num.effects)
        plot.new()
        plot.window(xlim = c(1,4), ylim = c(min(all.num.effects), max(all.num.effects)))
        apply(all.num.effects, 1, function(x) points(x, type = "b", col = "gray"))    
        axis(2)
        mtext(colnames(pheno)[p], side = 3, line = -1.5)
    }

}
