# helper function to filter genes out of cummeRbund fpkm matrix based on status flag and/or 
# mean fpkm across samples
filterCuff <- function(cuff, dat, filter, OK.only, fpkm) {
  # build SQLite query
  q <- 'select distinct gene_id from geneExpDiffData WHERE' # base query
  if(OK.only) q <- paste(q, "status NOT LIKE 'OK'")
  else {
    filt <- vector("character", length(filter)) # initialize vector to hold formatted filter terms
    for(i in 1:length(filter)) { filt[i] <- paste0("status LIKE '", filter[i], "'")} # format filter terms
    q <- paste(q, filt[1]) # add first filter term to base query; if more terms loop through and add them
    if(length(filt) > 1) {
      for(item in filt[2:length(filt)]) q <- paste(q, item, sep=' OR ')
    }
  }

  # execute query on the cuffDb object
  genesToRemove <- dbGetQuery( DB(cuff), q)

  # remove rows containing genes returned by query
  dat <- dat[!(rownames(dat) %in% genesToRemove$gene_id), ]
  
  # if fpkm >= 0 then filter rows with mean expression across all samples < fpkm
  if(fpkm >= 0){
    dat <- dat[rowMeans(dat) > fpkm, ]
  }

  # print dimensions of matrix to give indication of how many genes remain and return the modified matrix
  print( paste('Number of Genes:', dim(dat)[1]) )
  print( paste('Number of Samples:', dim(dat)[2]) )

  # return results
  dat

}


# Helper function for ggplot2 plotting of MDS plots
# Input is a dataframe with minimum columns M1, M2, and Condition which will be mapped to x, y, and color, respectively.
# Code also checks for an optional "Differentiation" column and maps that to shape when present (useful for cell culture expts).
# Color scale is currently optimized for a specific 5 sample experiment. Future version should take color scale as an
# argument to allow more flexible reusability.
# Returns value is the ggplot object.
publishPlot <- function(df) {
  p <- ggplot(df)
  if('Differentiation' %in% colnames(df)) {
    p <- p + geom_point(aes(x = M1, y = M2, color = Condition, shape=Differentiation), size = 5)
  } else {
    p <- p + geom_point(aes(x = M1, y = M2, color = Condition), size = 5) 
  }
  p <- p + scale_color_manual(values=c('#231F20', '#4E4E4F', '#868686', '#F1142B', '#C0142B')) +
          theme_bw()
  return(p)
}

# modification of cummeRbund MDSplot function for class CuffData to use only genes passing user-defined
# filter criteria based on status flag and mean fpkm across all samples
# also has argument altPlot which produces a plot without text labels when TRUE
MDSplot.filter <-  function (object, filter=c('HIDATA', 'FAIL'), OK.only=FALSE, fpkm = -1,
  replicates = FALSE, logMode = TRUE, altPlot = FALSE, sampleMap=NULL, batchMap=NULL, ...) 
  {
    if(altPlot & is.null(sampleMap)) stop('altPlot =T but no valid sampleMap found')
    .local <- function (object, replicates = FALSE, logMode = TRUE, 
                        pseudocount = 1) 
    { 
      if (replicates) {
        dat <- repFpkmMatrix(object)
      }
      else {
        dat <- fpkmMatrix(object)
      }

      # call filterCuff helper function to perform the actual filtering
      dat <- filterCuff(object, dat, filter, OK.only, fpkm)

      if (logMode) {
        dat <- log10(dat + pseudocount)
      }
      d <- JSdist(makeprobs(dat))
      fit <- cmdscale(d, eig = TRUE, k = 2)
      res <- data.frame(names = rownames(fit$points), 
                        M1 = fit$points[, 1], 
                        M2 = fit$points[, 2])
      p <- ggplot(res)

      if(altPlot) {
        if(is.null(names(sampleMap))) {
          warning('Cannot verify that sampleMap matches data because names(sampleMap) = NULL')
        } else {
            stopifnot(colnames(dat) == names(sampleMap))
        }
        
        res$Condition <- sampleMap
        if(!is.null(batchMap)) {res$Differentiation <- as.factor(batchMap)}
        p <- publishPlot(res)
      } else {
          p <- p + geom_point(aes(x = M1, y = M2, color = names), size=3) + 
            geom_text(aes(x = M1, y = M2, label = names, color = names)) + 
            theme_bw()
      }
      p
    }
    .local(object, replicates, logMode, ...)
  }


# modification of cummeRbund csDendro code for CuffData class that uses only genes passing user-defined filter criteria
# based upon status flag and mean fpkm across all samples
csDendro.filter <- function (object, logMode = TRUE, pseudocount = 1, replicates = FALSE, filter=c('HIDATA', 'FAIL'),
	OK.only=FALSE, fpkm = -1, ret.hclust=F, ...) 
{
  if (replicates) {
    fpkmMat <- repFpkmMatrix(object)
  }
  else {
    fpkmMat <- fpkmMatrix(object)
  }

  # call filterCuff helper function to perform the actual filtering
  fpkmMat <- filterCuff(object, fpkmMat, filter, OK.only, fpkm)
  
  if (logMode) {
    fpkmMat <- log10(fpkmMat + pseudocount)
  }
  res <- JSdist(makeprobs(fpkmMat))
  res2 <- as.dendrogram(hclust(res))
  plot(res2, main = paste("All", deparse(substitute(object)), 
                         sep = " "), ...)
  if(ret.hclust) {
    return (res2)
  }
  res
}



# csDensity cummeRbund method altered to change plotting output -- no fill and no grid to make
# visual comparison of overlap easier
csDensity_custom <- function (object, logMode = TRUE, pseudocount = 1, labels, features = FALSE,
                              removeZeroes=F, ...) 
{
  .local <- function (object, logMode = TRUE, labels, features = FALSE, replicates = FALSE, ...) 
  {
    if (is(object, "CuffData")) {
      if (replicates) {
        dat <- repFpkm(object, features = features)
        colnames(dat)[colnames(dat) == "rep_name"] <- "condition"
      }
      else {
        dat <- fpkm(object, features = features)
        colnames(dat)[colnames(dat) == "sample_name"] <- "condition"
      }
    }
    else {
      stop("Un-supported class of object.")
    }
    if (logMode) 
      dat$fpkm <- dat$fpkm + pseudocount
    p <- ggplot(dat)
    if (logMode) {
      p <- p + geom_density(aes(x = log10(fpkm), group = condition, 
                                color=condition), alpha = I(1/6), linetype='dashed')
    }
    else {
      p <- p + geom_density(aes(x = fpkm, group = condition, 
                                color=condition), alpha = I(1/6), linetype='dashed')
    }
    p <- p + labs(title = object@tables$mainTable)
    p <- p + scale_fill_hue(l = 50, h.start = 200) + scale_color_hue(l = 50, 
                                                                     h.start = 200)
    if(removeZeroes) {
      if(pseudocount == 0) {
        warning('Pseudocount = 0 so argument removeZeroes ignored')
      } else {
        p <- p + xlim(0.8*(log10(pseudocount)) + pseudocount, max(log10(dat$fpkm)))
        }
    }
    p + theme_classic()
  }
  .local(object, logMode, labels, features, ...)
}


# Modification of cummeRbund PCA plotting function that enables filtering of genes based on Cuffdiff status flag
# and/or fpkm expression level. Also modifies plot appearance to include percent variation explained in the axis
# labels.
PCAplot <-  function (object, replicates = FALSE, logMode = TRUE, filter=c('HIDATA', 'FAIL'), 
  OK.only=FALSE, fpkm = -1, altPlot = FALSE, ...) 
{
  .local <- function (object, replicates = FALSE, logMode = TRUE, 
                      pseudocount = 1) 
  { 
    if (replicates) {
      dat <- repFpkmMatrix(object)
    }
    else {
      dat <- fpkmMatrix(object)
    }

    # call filterCuff helper function to perform the actual filtering
    dat <- filterCuff(object, dat, filter, OK.only, fpkm)
    
    if (logMode) {
      dat <- log10(dat + pseudocount)
    }
    fit <- prcomp(t(dat))
    pve <- fit$sdev^2 / sum(fit$sdev^2)
    res <- data.frame(names = rownames(fit$x),
                      PC1 = fit$x[, 1], PC2 = fit$x[, 2])
    p <- ggplot(res)
    if(altPlot) {
      p <- p + geom_point(aes(x = PC1, y = PC2, color = names), size = 5) +
        xlab(paste('PC1 (', round(pve[1] * 100, 1), '%)', sep='')) +
        ylab(paste('PC2 (', round(pve[2] * 100, 1), '%)', sep='')) +
        theme_bw()
    } else {
      p <- p + geom_point(aes(x = PC1, y = PC2, color = names), size=3) + 
        geom_text(aes(x = PC1, y = PC2, label = names, color = names)) +
        xlab(paste('PC1 (', round(pve[1] * 100, 1), '%)', sep='')) +
        ylab(paste('PC2 (', round(pve[2] * 100, 1), '%)', sep='')) +
        theme_bw()
    }
    p
  }
  .local(object, replicates, logMode, ...)
}


# function to generate a .rnk file based on log2-fold-change for use with GSEA
# usage: specify the gene_exp.diff dataframe containing data and the contrast of interest
# contrast must be a valid contrast based on the cuffdiff run
# returns a formatted .rnk list which needs to be written to a file for use with GSEA:
# e.g., write.table(rnkFile, 'testRnk.rnk', quote=F, sep='\t', row.names=F)
makeRnkFile <- function(data.diff, contrast, pseudocount = 0.0001) {
  dat <- subset(data.diff, sample_1 == contrast[1] & sample_2 == contrast[2])
  # if invalid contrast an empty dataframe will be returned; check for this and exit
  # with a message if necessary
  if(nrow(dat) == 0) {
    print('No data for specified contrast. Check that contrast is valid for this experiment.')
    return
  }
  # filter out genes with status FAIL, LOWDATA, HIDATA, and NOTEST
  dat <- subset(dat, status == 'OK')
  dat$fc <- log2(dat$value_2 + pseudocount) - log2(dat$value_1 + pseudocount)
  rnk <- data.frame(gene=dat$gene, log2fc=dat$fc)
  colnames(rnk)[1] <- "#gene"
  rnk <- rnk[order(rnk$log2fc, decreasing=T), ]
  return(rnk)
}

