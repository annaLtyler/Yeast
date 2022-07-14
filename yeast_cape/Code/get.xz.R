#This function aligns two matrices based on their rownames
#We use this to get X and Z matrices for CCA.

get.xz  <- function(x.mat, z.mat){
    common.ind <- intersect(rownames(x.mat), rownames(z.mat))
    common.x.locale <- match(common.ind, rownames(x.mat))
    common.z.locale <- match(common.ind, rownames(z.mat))

    ordered.x <- x.mat[common.x.locale,]
    ordered.z  <- z.mat[common.z.locale,]

    x.na.locale <- which(is.na(ordered.x), arr.ind = TRUE)
    z.na.locale <- which(is.na(ordered.z), arr.ind = TRUE)

    na.locale <- rbind(x.na.locale, z.na.locale)
    if(nrow(na.locale) > 0){
        ordered.x <- ordered.x[-na.locale[,1],]
    }
    if(nrow(na.locale) > 0){
        ordered.z <- ordered.z[-na.locale[,1],]
    }
    return(list("X" = ordered.x, "Z" = ordered.z))
}
