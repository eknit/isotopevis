#' Standard Species
#' 
#' Takes a dataframe and identifies the species column (either "species", "Species", 
#' "taxon" or "Taxon") and matches the listed species against a dictionary list of 
#' standardized English names in the "species_lookup_table" The standardized names 
#' are then looked up in a standardized plotting table (plot_lookup_table) to 
#' find the plotting parameters. 
#' 
#' @param df The data.frame in which the isotope ratios 
#'  for all species are stored. 
#' @return A data.frame with the same number of rows as the original data frame,
#'    with the standardized English species name and other standardized grouping 
#'    information for each entry
#' @examples 
#' x <- standard_species(test_df)
#'

standard_species <- function(df){
  names <- colnames(df)
  species_names <-  df[,grepl("species", names) |
                         grepl("Species", names) |
                         grepl("taxon", names) |
                         grepl("Taxon", names)]
  species_lookup <- species
  plot_params <- plot_params
  matched_species <- species_lookup[match(
      species_names, species_lookup[, 2]), 1]
  matched_plot_params <- plot_params[match(
      matched_species, plot_params[, 1]), 1:ncol(plot_params)]
  standard <- cbind(df, matched_plot_params)
  data.frame(standard[, 2:ncol(standard)], row.names=1:nrow(standard))
}

#' Add.alpha
#' 
#' Converts colours to transparent.
#' 
#' @param col A color. 
#' @return the transparent version of the colour
#' @examples 
#' add.alpha("blue")
#'
add.alpha <- function(col, alpha=1){ # Converts colours to transparent
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#' d13C
#' 
#' Convenient text for plotting.
d13C <- expression(paste(delta^{13}, "C (‰)"))

#' d15N
#' 
#' Convenient text for plotting.
d15N <- expression(paste(delta^{15}, "N (‰)"))

#' plotmyxy
#' 
#'  Makes a standard isotope plot using the standard plotting charactes from the standard.species()
#' function. Can take all the normal graphical parameters for further customization. It can also
#' use big D13C or little d13C. The default (bigD=FALSE) is for little d13C. If bigD is selected, the 
#' D13C column must be included in the original data.frame. 
plotmyxy <- function(df, bigD = FALSE, dair=NULL, ...){
  if (bigD==TRUE){
    x <- df$D13C
    xlab <- expression(paste(Delta^{13}, "C (‰)"))
  }
  else {
    x <- df$normd13C
    xlab <- expression(paste(delta^{13}, "C (‰)"))
  }
  
  ylab <- expression(paste(delta^{15}, "N (‰)"))
  y <- df$normd15N
  
  #Plotting parameters
  par(mar=c(5,5,2,10))
  plot(x, y, col=paste(df$col.list),
       pch=as.numeric(as.character(df$pch.list)),
       bg=paste(df$bg.list),
       xlab=xlab, ylab=ylab,
       ...)
  
  #Legend
  legtab <-  df[!duplicated(df$English),]
  legtab <- legtab[order(as.numeric(as.character(legtab$Order))),]
  par(new=T)
  par(mar=c(1,1,1,1))
  plot(1:10, 1:10, axes=F, xlab="", ylab="", type="n")
  par(xpd=T)
  legend("right", paste(legtab$English), pch=as.numeric(as.character(legtab$pch.list)), 
         col=paste(legtab$col.list), pt.bg=paste(legtab$bg.list), cex=0.9, pt.cex=1.3, bty="n") 
}

#' Get model names
#' 
#' Extracts the category names from the list of those provided in the data frame df

get_model_names <- function(df){
  types <- as.character(df$Type[!duplicated(df$model_groups)])
  group_names <- df$model_groups[!duplicated(df$model_groups)]
  model_names <- cbind(types, group_names)
  model_names <- model_names[complete.cases(model_names),]
  model_names
}

#'  Make summary table

make_summary_table <- function(df){
  s <- ddply(df, .(model_groups), summarize, md13C=mean(normd13C), 
             d13Csd=sd(normd13C), 
             md15N=mean(normd15N),
             d15Nsd=sd(normd15N),
            Type=unique(Type)[1],
            Group=unique(Group)[1],
            model_groups=unique(model_groups)[1])
  s 
}

#' Make endpoints
#' 
#' Takes isotope values of the selected groups and divides them by type (Cereal, Pulse or Animal)
#' and addes the appropriate discrimination factor from the params table, and pools the two standard
#' deviations. It also assigns the appropriate digestable C and N proportions to each type of mixing 
#' group based on the params table.


make_endpoints <- function(df, params){
  
  s <- make_summary_table(df)
  model_names <- get_model_names(df)
  #Cereals
  model_names[match(
    species_names, species_lookup[, 2]), 1]
  ecer <- subset(s, model_groups %in% c("Barley", "Cereal", "Millet"))
  ecer$md13C <- ecer$md13C + params[1,2]
  ecer$md15N <- ecer$md15N + params[1,4]
  ecer$d13Csd <- sqrt(ecer$d13Csd^2 + params[1,3]^2)
  ecer$d15Nsd <- sqrt(ecer$d15Nsd^2 + params[1,5]^2)
  other_params <- data.frame(matrix(rep(params[1, 6:9], 3),  ncol=4, byrow=T))
  names(other_params)<- names(params[6:9])
  ecer <- cbind(ecer, other_params)
  
  #Pulses
  #Cereals
  epul <- subset(s, model_groups %in% "Pulse")
  epul$md13C <- epul$md13C + params[2,2]
  epul$md15N <- epul$md15N + params[2,4]
  epul$d13Csd <- sqrt(epul$d13Csd^2 + params[2,3]^2)
  epul$d15Nsd <- sqrt(epul$d15Nsd^2 + params[2,5]^2)
  other_params <- data.frame(matrix(rep(params[2, 6:9], 1), ncol=4, byrow=T))
  names(other_params)<- names(params[6:9])
  epul <- cbind(epul, other_params)
  
  #Pulses
  #Cereals
  ean <- subset(s, model_groups %in% c("Cattle", "Fish", "Sheep/Goat/Pig", "Wild Canopy"))
  ean$md13C <- ean$md13C + params[3,2]
  ean$md15N <- ean$md15N + params[3,4]
  ean$d13Csd <- sqrt(ean$d13Csd^2 + params[3,3]^2)
  ean$d15Nsd <- sqrt(ean$d15Nsd^2 + params[3,5]^2)
  other_params <- data.frame(matrix(rep(params[3, 6:9], 1), ncol=4, byrow=T))
  names(other_params)<- names(params[6:9])
  ean <- cbind(ean, other_params)
  
  e <- rbind(ecer, epul, ean)
  e$order <- c(1,2,6,3,4,8,5,7)
  e <- e[order(e$order),]
  
  e$pcDC <- unlist(e$DigestC)/(unlist(e$DigestC)+unlist(e$DigestN))
  e$pcDN <- unlist(e$DigestN)/(unlist(e$DigestC)+unlist(e$DigestN))
  e
}

#' Calculates the hypothetical d13C and d15N ratios for the selected set of proportions (model)
#' The proportions from the model are modified by the C/N ratio of the digestible C and N.
#' From the endpoints calculated by the endpoints() function 10,000 normally distributed values
#' are selected, for each of the endpoint groups. Each of these 10,000 ranomly selected values
#' are averaged, weighted by the digestible C/N ratio-corrected proportions from the selected model
#' The distribution of these 10,000 weighted mixtures is reported for both x and y.

make_xyvals <- function(e, model){
  props <- model/100
  xconc <- (props*e$pcDC)/sum(props*e$pcDC)
  yconc <- (props*e$pcDN)/sum(props*e$pcDN)
  B <- 10000
  sample.matrix <- matrix(nrow=B, ncol=length(xconc))
  for (i in 1:length(xconc)){
    sample.matrix[,i] <- rnorm(B, e$md13C[i], e$d13Csd[i])
  }
  prop.matrix <- sweep(sample.matrix, MARGIN=2, as.numeric(xconc), "*")
  xvals <- rowSums(prop.matrix)
  for (i in 1:length(xconc)){
    sample.matrix[,i] <- rnorm(B, e$md15N[i], e$d15Nsd[i])
  }
  prop.matrix <- sweep(sample.matrix, MARGIN=2, as.numeric(xconc), "*")
  yvals <- rowSums(prop.matrix)
  x <- data.frame(xvals=xvals[1:B-1], yvals=yvals[1:B-1])
  x
}

make_xyvals2 <- function(e, model){
  props <- model/100
  xconc <- (props*e$pcDC)/sum(props*e$pcDC)
  yconc <- (props*e$pcDN)/sum(props*e$pcDN)
  B <- 10000
  sample.matrix <- matrix(nrow=B, ncol=length(xconc))
  for (i in 1:length(xconc)){
    sample.matrix[,i] <- rnorm(B, e$md13C[i], e$d13Csd[i])
  }
  prop.matrix <- sweep(sample.matrix, MARGIN=2, as.numeric(xconc), "*")
  xvals <- rowSums(prop.matrix)
  for (i in 1:length(xconc)){
    sample.matrix[,i] <- rnorm(B, e$md15N[i], e$d15Nsd[i])
  }
  prop.matrix <- sweep(sample.matrix, MARGIN=2, as.numeric(xconc), "*")
  yvals <- rowSums(prop.matrix)
  x <- data.frame(xvals=xvals[1:B-1], yvals=yvals[1:B-1])
  x
}

#' This function plots the resulting x and y values as a 2d histogram using the MASS package, and then
#' overlays a barplot of the relative proportions for the selected model, as well as the human 
#' d13C and d15N values for reference. The cex values are adjusted so that it appears best as a 3x4
#' matrix with other plots of the same type (implemented below)

#' @param colval A number between 1 and 360
model_plot <- function(colval, x, this.model=this.model){
  if (this.model$Works=="Yes"){
    r <- rev(sequential_hcl(32, h=colval, c=c(80, 0), l=c(30, 100)))
  }
  else {
    r <- rev(sequential_hcl(32, h=colval, c=c(0, 0), l=c(30, 100)))
  }
  xlims=c(-24, -10)
  ylims=c(4, 18)
  par(mar=c(4,5,2,0))
  par(cex=1)
  k <- kde2d(x$xvals, x$yvals, n=50, lims=c(xlims+c(-0.2,0.2), ylims+c(-0.2, 0.2)))
  image(k, col=r, axes=F)
  par(new=T)
  plot(humans$normd13C, humans$normd15N, col="black", pch=8, xlim=xlims, ylim=ylims, xlab=d13C,
       ylab=d15N, main=this.model$ModelName, cex.main=0.92)
  par(new=T)
  barplot(as.numeric(this.model[,2:9]), col=r[32], ylim=c(-300,100),
          xlim=c(-12, 9), axes=F, cex.axis=0.8)
  par(xpd=T)
  text(((1:8)*1.2)-0.3, -5, c("Barley", "Cereal", "Pulses", "Cattle", "SGP", "Millet", "Wild", "Fish"), cex=0.8, srt=60, adj=1)
  axis(2, at=seq(0,100, by=20), pos=0, cex.axis=0.8)
}

#' This function combines all of the above plots - note that you can select which parameter table
#' you'd like to use here, for easy comparison
full_plot <- function(modelno, params, s){
  this.model <- models[modelno,]
  e <- make_endpoints(s, params)
  x <- make_xyvals(e, this.model[2:9])
  if (this.model$Works=="Yes"){
    colval <- 360/(nrow(models[models$Works=="Yes",])) + (360/(nrow(models[models$Works=="Yes",]))*modelno)
  }
  else {colval <- 260}
  model_plot(colval, x, this.model)
}

