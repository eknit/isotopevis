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
  species_names <-  df[,grepl("species", names, ignore.case=T) |
                         grepl("Species", names, ignore.case=T)]
  if (ncol(species_names)>1) warning("Multiple species match columns found, using only first")
  species_names <- species_names[,1]
  species_lookup <- species
  plot_params <- plot_params
  matched_species <- species_lookup[match(
    species_names, species_lookup[, 2]), 1]
  matched_plot_params <- plot_params[match(
    matched_species, plot_params[, 1]), 1:ncol(plot_params)]
  standard <- cbind(df, matched_plot_params)
  data.frame(standard[, 1:ncol(standard)], row.names=1:nrow(standard))
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

#' make_tint
#' Automatically makes tinted shade for plotting. Returns RGB colour.
#' 
make_tint <- function(col, tint_factor){
  apply(x <- sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1]+(1-x[1])*tint_factor,
              x[2]+(1-x[2])*tint_factor,
              x[3]+(1-x[3])*tint_factor))
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

#' generate_mvnorm
#' 
#' Generates a bivariate normal distribution of the dietary group supplied 
#' (batch) using the mvnorm function
#' 
generate_mvnorm <- function(batch, N){
  Cvar <- var(batch$normd13C)
  Nvar <- var(batch$normd15N)
  Cmu <- mean(batch$normd13C)
  Nmu <- mean(batch$normd15N)
  batch_cov <- cov(batch$normd13C, batch$normd15N)
  Sigma <- matrix(c(Cvar, batch_cov, batch_cov, Nvar), 2,2)
  result <- mvrnorm(N, c(Cmu, Nmu), Sigma)
  result
}

#' bi_endpoints_sim
#' 
#' Simulates bivariate endpoints by adding the distribution of the C and N discrimination
#' factors to the bivariate normal distribution of the food source generated by the
#' make_mvnorm function.
#' 
#'
 bi_endpoints_sim <- function(batch, N){
   if (batch$Type[1]=="Plant" & s$Group[1] != "Pulse"){
     parCmu <- params[1,2]
     parCvar <- (params[1,3])^2
     parNmu <- params[1,4]
     parNvar <- (params[1,5])^2
   }
   else if (batch$Group[1] == "Pulse"){
     parCmu <- params[2,2]
     parCvar <- (params[2,3])^2
     parNmu <- params[2,4]
     parNvar <- (params[2,5])^2
   }
   else if (batch$Type[1] == "Animal"){
     parCmu <- params[3,2]
     parCvar <- (params[3,3])^2
     parNmu <- params[3,4]
     parNvar <- (params[3,5])^2
   }
   Sigma2 <- matrix(c(parCvar, 0, 0, parNvar), 2,2)
   result2 <- mvrnorm(N, c(parCmu, parNmu), Sigma2)
   result2
}

#' model_plot function - version 2 using kde instead of histogram, based on multivariate distribution
#' 
model_plot <- function(s, model.no, ...){
  this.model <- models[model.no, ]
  if (this.model$Works=="No"){
    colval <- "grey20"
  }
  else {
    colval <- model.no
  }
  print(colval)
  xlims=c(min(s$md13C)+1, max(s$md13C)+(abs(min(s$md13C)-max(s$md13C))))
  ylims=c(min(s$md15N)+1, max(s$md15N)+(abs(min(s$md15N)-max(s$md15N))*1.8))
  par(cex=1)
  endpoints <- make_xyvals2(s, model=this.model[2:(ncol(this.model)-2)])
  H <- Hpi(x=endpoints)
  fhat <- kde(x=endpoints, H=H)
  plot(fhat, display="filled.contour2", cont=c(95, 67), col=c('white', 		make_tint(colval, 0.5), colval), xlim=xlims, ylim=ylims, axes=F, xlab="", ylab="")
  points(humans$normd13C, humans$normd15N, pch=8)
  par(new=T)
  barplot(as.numeric(this.model[,2:(ncol(this.model)-2)]), col=colval, ylim=c(-300,100),
          xlim=c((ncol(this.model)-2)*-1.3, (ncol(this.model)-2)), axes=F, cex.axis=0.6)
  par(xpd=T)
  model_names <- get_model_names(df)[,2]
  text(((1:length(model_names))*1.2)-0.3, -5, model_names, cex=0.6, srt=60, adj=1)
  axis(2, at=seq(0,100, by=20), pos=0, cex.axis=0.6)
  par(new=T)
  plot(1:10, 1:10, xlim=xlims, ylim=ylims, axes=F, xlab="", ylab="")
}

#' full_plot makes grid of all plots
full_plot <- function(models, s, ncol=3, nrow=4){
  palette(c("#FFFF54","#C8E64C", "#8CD446", "#4DC742", "#45D2B0", "#46ACD3", "#438CCB", "#4262C7", "#5240C3", "#8C3FC0", "#D145C1", "#E64C8D", "#FF5454", "#FF8054", "#FFA054", "#FFB554"))
  N <- nrow(models)
  par(mfrow=c(nrow,ncol))
  side <- seq(1, N, by=ncol)
  bottom <- (N-ncol+1):N
  par(oma=c(5,5,1,1))
  par(mar=c(0,0,0,0))
  for (i in 1:N){
    model_plot2(s, i)
    if (i %in% side){
      axis(2)
      mtext(d15N, 2, 2, cex=0.9)
    }
    if (!i %in% side){
      axis(2, tck=0.02, labels=NA, cex=0.7)
    }
    if (i %in% bottom){
      axis(1)
      mtext(d13C, 1, 2, cex=0.9)
      
    }
    if (!i %in% bottom){
      axis(1, tck=0.02, labels=NA, cex=0.8)
    }
    box()
  }
}