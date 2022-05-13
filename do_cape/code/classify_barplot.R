#classify barplot results
classify_barplot <- function(barplot_results){
    pheno.max <- barplot_results[[1]] + barplot_results[[2]]
    pheno.min <- barplot_results[[1]] - barplot_results[[2]]

    compare_main <- function(min1, max1, min2, max2){
        do.overlap <- segments.overlap(min1, max1, min2, max2)
        if(!do.overlap){
            if(min2 > max1){
                return("increases trait")
            }
            if(min2 < min1){
                return("decreases trait")
            }
        }else{
            return("none")
        }
    }

    compare_int <- function(min1, max1, min2, max2){
        do.overlap <- segments.overlap(min1, max1, min2, max2)
        if(!do.overlap){
            if(abs(max2) > abs(max1)){
                if(max2 < 0){
                    return("more negative")
                }else{
                    return("more positive")
                }
            }
            if(abs(max2) < abs(max1)){
                return("less than additive")
            }
        }else{
            return("none")
        }
    }


    add_min <- sum(barplot_results[[1]][2:3]) - sum(barplot_results[[2]][2:3])
    add_max <- sum(barplot_results[[1]][2:3]) + sum(barplot_results[[2]][2:3])

    #get categories for two main effects
    main1 <- compare_main(pheno.min[1], pheno.max[1], pheno.min[2], pheno.max[2])
    main2 <- compare_main(pheno.min[1], pheno.max[1], pheno.min[3], pheno.max[3])

    #get category for interaction
    int <- compare_int(add_min, add_max, pheno.min[4], pheno.max[4])

    result <- c("main1" = main1, "main2" = main2, "interaction" = int)
    return(result)

}
