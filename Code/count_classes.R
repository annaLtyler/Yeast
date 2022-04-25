count_classes <- function(class_table, plot.result = TRUE){
    
    u_classes <- unique(class_table)
    class.tally <- apply(apply(u_classes, 1, function(x) apply(class_table, 1, function(y) identical(x,y))), 2, function(z) which(z))
    class.count <- sapply(class.tally, length)
    names(class.count) <- names(class.tally) <- apply(u_classes, 1, function(x) paste(x, collapse = "_"))

    if(plot.result){
        par(mar = c(4,20,4,2))
        barplot(sort(class.count), las = 2, horiz = TRUE, cex.names = 0.7)
    }

    result <- list("tally" = class.tally, "count" = class.count)
    return(result)

}