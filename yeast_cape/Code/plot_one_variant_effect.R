plot_one_variant_effect <- function(data_obj, geno_obj, marker_name,
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

    test.source.idx <- grep(marker_name, var_int[,"Source"]) 
    test.target.idx <- grep(marker_name, var_int[,"Target"])

    if(length(test.source.idx) == 0 && length(test.target.idx) == 0){
        return("No interactions found for this marker")
    }
        
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

    source.deviation <- vector(mode = "list", length = length(test.source.idx))
    names(source.deviation) <- var_int[test.source.idx,"Source"]
    
    if(length(test.source.idx) > 0){
        for(s in 1:length(test.source.idx)){
            split.source <- strsplit(var_int[test.source.idx[s],"Source"], "_")
            source.marker.names <- sapply(split.source, function(x) x[1])
            source.allele.names <- sapply(split.source, function(x) x[2])

            source.geno <- sapply(1:length(source.marker.names), function(x) geno_obj[,source.allele.names[x], source.marker.names[x]])
            colnames(source.geno) <- source.marker.names
            
            source.int <- lapply(1:ncol(pheno), function(y) t(apply(source.geno, 2, function(x) get_int(pheno[,y], x, query_genotype))))

            for(p in 1:length(source.int)){
                barplot(source.int[[p]], 
                main = paste(colnames(pheno)[p], names(source.deviation)[s], sep = "\n"), 
                names = c(names(source.deviation)[s], "query", "Additive", "Actual"))
                abline(h = 0)
            }

            source.effects <- sapply(source.int, function(x) x[,"int"] - x[,"add"])
            if(!is.null(nrow(source.effects))){
                colnames(source.effects) <- colnames(pheno)
            }else{
                names(source.effects) <- colnames(pheno)
            }
            source.deviation[[s]] <- source.effects
        }
    }

    target.deviation <- vector(mode = "list", length = length(test.target.idx))
    names(target.deviation) <- var_int[test.target.idx,"Target"]

    if(length(test.target.idx) > 0){
        for(s in 1:length(test.target.idx)){
            split.target <- strsplit(var_int[test.target.idx[s],"Target"], "_")
            target.marker.names <- sapply(split.target, function(x) x[1])
            target.allele.names <- sapply(split.target, function(x) x[2])

            target.geno <- sapply(1:length(target.marker.names), function(x) geno_obj[,target.allele.names[x], target.marker.names[x]])
            colnames(target.geno) <- target.marker.names
            
            target.int <- lapply(1:ncol(pheno), function(y) t(apply(target.geno, 2, function(x) get_int(pheno[,y], query_genotype, x))))

            for(p in 1:length(target.int)){
                barplot(target.int[[p]], 
                main = paste(colnames(pheno)[p], names(target.deviation)[s], sep = "\n"), 
                names = c(names(target.deviation)[s], "query", "Additive", "Actual"))
                abline(h = 0)
            }

            target.effects <- sapply(target.int, function(x) x[,"int"] - x[,"add"])
            if(!is.null(nrow(target.effects))){
                colnames(target.effects) <- colnames(pheno)
            }else{
                names(target.effects) <- colnames(pheno)
            }
            target.deviation[[s]] <- target.effects
        }

    }
    results <- list("deviation_with_query_as_source" = target.deviation,
        "deviation_with_query_as_target" = source.deviation)
    invisible(results)

}
