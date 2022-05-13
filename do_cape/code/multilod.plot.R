multilod.plot <- function(scan1.results, map, lod.thresh = 0, 
color.scheme = c("blue", "green", "purple", "red", "orange", "brown", "yellow", "gray"),
border.col = "darkgray", border.lwd = 3, row.names = colnames(scan1.results), 
row.name.shift = 0, chr.name.y = 0.2, color.bar.cex = 1, color.bar.axis.lin = 0, 
color.fun = c("linear", "exponential"), steepness = 1){

    val.mat <- as.matrix(scan1.results)
    low.lod <- which(val.mat < lod.thresh)
    if(length(low.lod) > 0){
        val.mat[low.lod] <- NA
        has.vals <- which(apply(val.mat, 1, function(x) !all(is.na(x))))
        val.mat <- val.mat[has.vals,]
    }
    col.scale <- color.scheme[1]

    layout(matrix(c(1,1,2,0), ncol = 2), widths = c(1, 0.1), heights = c(0.4, 1))
    par(mar = c(2, 6, 2, 0))
    imageWithText(t(val.mat), show.text = FALSE, col.names = NULL,
    col.scale = col.scale, color.fun = "exponential", exp.steepness = steepness,
    row.text.shift = row.name.shift, row.names = row.names)
    
    #demarcate the chromosomes on the plot
    markers <- rownames(val.mat)
    all.markers <- lapply(map, names)
    
    for(i in 1:length(all.markers)){
        chr.locale <- which(markers %in% all.markers[[i]])
        chr.min <- min(chr.locale)
        chr.max <- max(chr.locale)
        chr.mid <- mean(c(chr.min, chr.max))
        draw.rectangle(chr.min, chr.max, 0.5, ncol(val.mat)+0.5, 
        border.col = border.col, lwd = border.lwd)
        text(x = chr.mid, y = chr.name.y, names(all.markers)[i])
    }

    par(mar = c(2,2,2,2))
    imageWithTextColorbar(val.mat, col.scale = col.scale, cex = color.bar.cex,
    axis.line = color.bar.axis.lin, color.fun = "exponential", exp.steepness = steepness)
}