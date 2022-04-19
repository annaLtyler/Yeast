#This function looks at the main effects in different motifs

motif.effects <- function(data.obj, geno.obj, interaction = c("enhancing", "suppressing"), 
	main = c("coherent", "incoherent"), source.sign = c("pos", "neg", "any"), color.scheme = "blue",
	geno.bins = c(0, 0.5, 1)){
	
	# data.obj <- get.network(data.obj, p.or.q = p.or.q, collapse.linked.markers = TRUE)
	motif.obj <- find.motifs(data.obj)
	
	#collect the unique motifs for all phenotypes
	motifs <- NULL
	for(i in 1:length(motif.obj[[2]])){
		if(interaction == "enhancing"){
			motif.locale.int <- which(motif.obj[[2]][[i]][,1] == 1)
			}else{
			motif.locale.int <- which(motif.obj[[2]][[i]][,1] == -1)	
			}
		if(main == "coherent"){
			motif.locale.main <- which(motif.obj[[2]][[i]][,2] == motif.obj[[2]][[i]][,3])
			}else{
			motif.locale.main <- which(motif.obj[[2]][[i]][,2] != motif.obj[[2]][[i]][,3])	
			}
		if(source.sign == "pos"){
			sign.locale <- which(motif.obj[[2]][[i]][,2] == 1)
			}
		if(source.sign == "neg"){
			sign.locale <- which(motif.obj[[2]][[i]][,2] == -1)
			}
		if(source.sign == "any"){
			sign.locale <- 1:dim(motif.obj[[2]][[i]])[1]
			}
		motif.locale <- intersect(intersect(motif.locale.int, motif.locale.main), sign.locale)
		motifs <- rbind(motifs, motif.obj[[1]][[i]][motif.locale,])
		}
	

	u_pheno <- unique(motifs[,3])
	motif.pheno.effects <- vector(mode = "list", length = length(u_pheno))
	names(motif.pheno.effects) <- u_pheno
	for(ph in 1:length(u_pheno)){
		pheno.motifs <- motifs[which(motifs[,3] == u_pheno[ph]),1:2]
		u_motifs <- pheno.motifs
		pheno.motif.effects <- t(apply(pheno.motifs, 1, function(x) get.pheno.effects(data.obj, geno.obj, x[1], x[2], u_pheno[ph], data.obj$p_covar, geno.bins = geno.bins)))

		#get the phenotype effects for all motifs in the class
	
		collapsed.net <- data.obj$collapsed_net
		num.markers <- min(dim(collapsed.net))
		just.pheno <- collapsed.net[,(num.markers+1):dim(collapsed.net)[2]]
		
		source.locale <- match(u_motifs[,1], rownames(collapsed.net))
		source.mat <- just.pheno[source.locale,]
		source.pheno <- source.mat
		source.pheno[which(source.pheno == 0)] <- NA
		#boxplot(source.pheno)
		
		target.locale <- match(u_motifs[,2], rownames(collapsed.net))
		target.mat <- just.pheno[target.locale,]
		target.pheno <- target.mat
		target.pheno[which(target.pheno == 0)] <- NA	
		#boxplot(target.pheno)
	
	#cluster by the target phenotype effects, since this is a large part of how directionality is determined
	tm <- heatmap(target.mat)
	
	plot.mat <- function(mat){
		ymin <- min(mat);ymax <- max(mat)
		xmin <- 1;xmax <- ncol(mat)
		plot.height <- ymax - ymin
		par(mar = c(2, 2, 2, 2))
		plot.new()
		plot.window(xlim = c(xmin, xmax), ylim = c(ymin, ymax))
		for(i in 1:ncol(mat)){
			points(x = jitter(rep(i, nrow(mat))), mat[,i])
			}
		axis(2)
		par(xpd = TRUE)
		text(x = 1:ncol(mat), y = rep(ymin-plot.height*0.1, ncol(mat)), labels = colnames(mat))
		par(xpd = FALSE)
		}
	
	# quartz()
	layout.mat <- matrix(c(1:4), ncol = 2)
	layout(layout.mat, heights = c(1,0.5))
	par(mar = c(4,2,3,2))
	
	imageWithText(source.mat[tm$rowInd,], show.text = FALSE, 
	split.at.vals = TRUE, col.names = names(motif.obj[[1]]), col.text.adj = 0, 
	col.text.shift = 0, main = paste("Source\n", interaction, main, "\nSource =", source.sign), 
	col.scale = color.scheme, col.text.rotation = 0)
	par(mar = c(3,2,0,2))
	boxplot(source.mat, las = 2)
	#plot.mat(source.mat)
	par(mar = c(4,2,3,2))
	imageWithText(target.mat[tm$rowInd,], show.text = FALSE, split.at.vals = TRUE, col.names = names(motif.obj[[1]]), col.text.adj = 1, col.text.shift = 0, main = paste("Target\n", interaction, main, "\nSource =", source.sign), col.scale = color.scheme)
	par(mar = c(3,2,0,2))
	boxplot(target.mat, las = 2)
	#plot.mat(target.mat)

	ymin <- min(c(source.mat, target.mat))
	ymax <- max(c(source.mat, target.mat))

		
	# layout.mat <- get.layout.mat(ncol(just.pheno))
	# layout(layout.mat)
	# for(ph in 1:dim(source.mat)[2]){
		# plot.new()
		# plot.window(xlim = c(0,1), ylim = c(ymin, ymax))
		# segments(x0 = rep(0.2, dim(source.mat)[1]), y0 = source.mat[,ph], x1 = rep(0.8, dim(source.mat)[1]), y1 = target.mat[,ph])
		# segments(0,0,1,0, col = "red")
		# mtext(names(motif.obj[[2]])[ph], cex = 2)
		# text(x = 0.2, y = ymin-0.15, "Source")
		# text(x = 0.8, y = ymin-0.15, "Target")
		# mtext(paste(interaction, main, "\nSource =", source.sign), outer = TRUE, line = -2)
		# }

	# quartz()
	# imageWithText(cbind(source.sum, target.sum, (source.sum-target.sum)), show.text = FALSE, split.at.vals = TRUE)	
	
	}
}