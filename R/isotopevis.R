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