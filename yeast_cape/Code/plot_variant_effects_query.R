plot_variant_effects_query <- function(data_obj, geno_obj, query_genotype, 
    pheno_type = c("pheno", "ET"), p_or_q = 0.05){

    trait_cols <- categorical_pal(8)

    gene <- get_query_geno(data_obj, geno_obj, query_genotype)

    if(pheno_type == "pheno"){
        pheno <- data_obj$pheno[match(rownames(gene), rownames(data_obj$pheno)),]
    }else{
        pheno <- data_obj$ET[match(rownames(gene), rownames(data_obj$ET)),]
    }
    query <- query_genotype[match(rownames(gene), rownames(query_genotype)),1,drop=FALSE]

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
        group.means <- sapply(pheno.groups, mean)
        #group.se <- sapply(pheno.groups, function(x) sd(x)/sqrt(length(x)))
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
        target.chr.locale <- which(target.chr == ch)
        source.chr.locale <- which(source.chr == ch)

        target.markers <- var_int[target.chr.locale,"Target"]
        source.markers <- var_int[source.chr.locale,"Source"]
        
        target.main.idx <- match(target.markers, colnames(gene))
        source.main.idx <- match(source.markers, colnames(gene))

        target.int <- lapply(1:ncol(pheno), function(y) t(sapply(target.main.idx, function(x) get_int(pheno[,y], query, gene[,x]))))
        source.int <- lapply(1:ncol(pheno), function(y) t(sapply(source.main.idx, function(x) get_int(pheno[,y], query, gene[,x]))))

        source.ch.pos[[ch]] <- source.pos[source.chr.locale]
        target.ch.pos[[ch]] <- target.pos[target.chr.locale]

        target.deviation[[ch]] <- sapply(target.int, function(x) x[,"int"] - x[,"add"])
        source.deviation[[ch]] <- sapply(source.int, function(x) x[,"int"] - x[,"add"])
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
}