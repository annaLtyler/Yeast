plot_variant_effects_query <- function(data_obj, geno_obj, 
    pheno_type = c("pheno", "ET"), p_or_q = 0.05, verbose = FALSE){

    query_genotype <- data_obj$query_genotype
    trait_cols <- categorical_pal(8)
    pheno_type <- pheno_type[1]

    if(pheno_type == "pheno"){
        pheno <- data_obj$pheno
    }else{
        pheno <- data_obj$ET
    }

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
        
    get_int <- function(phenoV, geno1, geno2){
        geno_pairs <- pair.matrix(c(0,1), self.pairs = TRUE, ordered = TRUE)
        pheno.groups <- apply(geno_pairs, 1, function(x) phenoV[intersect(which(geno1 == x[1]), which(geno2 == x[2]))])
        #boxplot(pheno.groups)
        group.means <- sapply(pheno.groups, function(x) mean(x, na.rm = TRUE))
        #group.se <- sapply(pheno.groups, function(x) sd(x, na.rm = TRUE)/sqrt(length(x)))
        centered.means <- group.means - group.means[1]
        add.pred <- centered.means[2] + centered.means[4]

        result <- c("main1" = centered.means[4], "main2" = centered.means[2], 
            "add" = add.pred, "int" = centered.means[3])

    }

    source.deviation <- target.deviation <- vector(mode = "list", length = length(u_chr))
    names(source.deviation) <- names(target.deviation) <- u_chr

    source.ch.pos <- target.ch.pos <- vector(mode = "list", length = length(u_chr))
    names(source.ch.pos) <- names(target.ch.pos) <- u_chr
    
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

        target.geno <- sapply(1:length(target.marker.names), function(x) geno_obj[,target.allele.names[x], target.marker.names[x]])
        colnames(target.geno) <- target.markers
        source.geno <- sapply(1:length(source.marker.names), function(x) geno_obj[,source.allele.names[x], source.marker.names[x]])
        colnames(source.geno) <- source.markers
        
        target.int <- lapply(1:ncol(pheno), function(y) t(apply(target.geno, 2, function(x) get_int(pheno[,y], query_genotype, x))))
        source.int <- lapply(1:ncol(pheno), function(y) t(apply(source.geno, 2, function(x) get_int(pheno[,y], query_genotype, x))))

        source.ch.pos[[ch]] <- source.pos[source.chr.locale]
        target.ch.pos[[ch]] <- target.pos[target.chr.locale]

        target.effects <- sapply(target.int, function(x) x[,"int"] - x[,"add"])
        colnames(target.effects) <- colnames(pheno)

        source.effects <- sapply(source.int, function(x) x[,"int"] - x[,"add"])
        colnames(source.effects) <- colnames(pheno)

        target.deviation[[ch]] <- target.effects
        source.deviation[[ch]] <- source.effects
    }


    ylim = c(min(c(min(unlist(target.deviation), na.rm = TRUE),min(unlist(source.deviation), na.rm = TRUE))),
           max(c(max(unlist(target.deviation), na.rm = TRUE),max(unlist(source.deviation), na.rm = TRUE))))


    for(ch in 1:length(source.deviation)){
        #quartz(width = 10, height = 6)
        xlim <- c(1, max(c(source.ch.pos[[ch]], target.ch.pos[[ch]])))

        par(mfrow = c(2,1), mar = c(2,4,2,2))

        plot.new()
        plot.window(xlim = xlim, ylim = ylim)
        for(p in 1:ncol(pheno)){
            points(source.ch.pos[[ch]], source.deviation[[ch]][,p], type = "h",
            col = trait_cols[p])
        } 
        axis(1); axis(2)
        abline(h = 0)
        legend("topright", col = trait_cols[1:ncol(pheno)], lty = 1, lwd = 2,
            legend = colnames(pheno))
        mtext("Deviation of Trait from Additive", side = 2, line = 2.5)
        mtext("Query as Target", side = 3)

        plot.new()
        plot.window(xlim = xlim, ylim = ylim)
        par(mar = c(4,4,0,2))
        for(p in 1:ncol(pheno)){
            points(target.ch.pos[[ch]], target.deviation[[ch]][,p], type = "h",
            col = trait_cols[p])
        } 
        axis(1); axis(2)
        abline(h = 0)
        legend("topright", col = trait_cols[1:ncol(pheno)], lty = 1, lwd = 2,
            legend = colnames(pheno))
        mtext("Deviation of Trait from Additive", side = 2, line = 2.5)
        mtext("Position", side = 1, line = 2.5)
        mtext("Query as Source", side = 3, line = -2.5)

        mtext(paste("Chr", ch), side = 3, outer = TRUE, line = -1)
    }

    results <- list("deviation_with_query_as_source" = target.deviation,
        "deviation_with_query_as_target" = source.deviation)
    invisible(results)

}
