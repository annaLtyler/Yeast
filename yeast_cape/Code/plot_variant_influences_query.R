plot_variant_influences_query <- function(data_obj, geno_obj, p_or_q = 0.05){

    query_genotype <- data_obj$query_genotype
    trait_cols <- categorical_pal(8)

    var_inf <- write_variant_influences(data_obj, p_or_q = p_or_q, 
        include_main_effects = TRUE, mark_covar = FALSE, 
        write_file = FALSE)
    just_main <- var_inf[which(is.na(var_inf[,6])),]
    
    main_chr <- as.numeric(just_main[,2])
    main_pos <- as.numeric(just_main[,3])
    main_trait <- just_main[,4]
    u_trait <- unique(main_trait)
    main_col <- rep(NA, nrow(just_main))
    for(p in 1:length(u_trait)){
        main_col[which(main_trait == u_trait[p])] <- trait_cols[p]
    }
    main_effect <- as.numeric(just_main[,"Effect"])
    main.lim <- c(min(main_effect), max(main_effect))

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

    all.effect <- as.numeric(var_int[,"Effect"])
    int.lim <- c(min(all.effect), max(all.effect))

    u_chr <- sort(unique(c(source.chr, target.chr)))
    u_chr <- u_chr[-1] #take off the 0 chromosome. this is the query marker
    
    #quartz(width = 12, height = 6)
    
    #query as source
    for(ch in u_chr){
        #quartz(width = 8, height = 6)
        par(mfrow = c(2,1), mar = c(0,4,2,4))
        target.chr.locale <- which(target.chr == ch)
        source.chr.locale <- which(source.chr == ch)

        target.chr.pos <- target.pos[target.chr.locale]
        source.chr.pos <- source.pos[source.chr.locale]
        target.chr.effect <- all.effect[target.chr.locale]
        source.chr.effect <- all.effect[source.chr.locale]

        #plot main effects
        main.chr.locale <- which(main_chr == ch)
        main.chr.pos <- main_pos[main.chr.locale]
        main.chr.effect <- main_effect[main.chr.locale]
        main.chr.col <- main_col[main.chr.locale]
        
        chr.boundaries <- c(0, max(c(target.chr.pos, source.chr.pos, main.chr.pos)))
        
        plot.new()
        plot.window(xlim = chr.boundaries, ylim = main.lim)
        points(main.chr.pos, main.chr.effect, col = main.chr.col, type = "h") 
        legend("topright", col = trait_cols[1:length(u_trait)], lty = 1, lwd = 3,
            legend = u_trait, cex = 0.7)
        abline(h = 0)       
        axis(2)
        mtext(side = 2, "Test Marker Main Effect", line = 2.5)
        mtext(side = 3, paste("Chr", u_chr[ch]))

        par(mar = c(4,4,0,4))
        plot.new()
        plot.window(xlim = chr.boundaries, ylim = int.lim)
        points(target.chr.pos, target.chr.effect, type = "h", col = "#a6611a")
        points(source.chr.pos, source.chr.effect, type = "h", col = "#018571")
        abline(h = 0)
        axis(2);axis(1)
        mtext(side = 2, "Interaction Effect Size", line = 2.5)
        mtext(side = 1, "Position", line = 2.5)
        legend("topright", legend = c("Query as Source", "Query as Target"), 
            col = c("#a6611a", "#018571"), lty = 1, lwd = 3, cex = 0.7)
    }

}