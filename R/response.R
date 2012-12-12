# $Id: response.R 283 2011-12-21 12:15:27Z jbao $

library(sn)
library(fields)
library(splancs)
library(plyr)
library(doMC)
#library(Rsge)
library(ggplot2)

#fontsize <- 20
#theme_set(theme_bw(fontsize))

#if (FALSE) {
registerDoMC()
#print(Sys.getenv('R_LIBS_SITE'))

# Read command-line arguments
args=(commandArgs(TRUE))
for (i in 1:length(args)) {
    eval(parse(text=args[[i]]))
}
n <- as.numeric(Sys.getenv("SGE_TASK_ID"))
#filename <- '~/data/DREAM/gnw/scalefree2/gnw/Size1000/norm_perturbation/mds/pca_scalefree2-1_perturbation-1_1000_normed_euc'
#prefix <- 'test'
#outdir <- './'

# Load and fit the MDS data
#files <- c('~/data/hgf/pca_hgf_17726_normed_euc', '~/data/data_aceton_w_mds.txt', 
#    '~/data/cfu/mds_cfu_epo_21464_euc', '~/data/yeast/pca_yeast_alpha_6178_normed_euc')
mds.fit <- function(filename) {
    coord <- read.delim(filename, header=F, sep=' ')
    coord <- coord[,1:2]
    #radius <- sqrt(rowSums((coord-kronecker(matrix(1,dim(coord)[1],1),
    #    t(mean(coord))))^2))
    
    # Fit a skew Gaussian distribution to the Data
    skewt_fit <- msn.mle(y=coord)
    # Extract the parameters of the distribution
    dp <- skewt_fit$dp

    # Create a regular grid with p values
    pp <- 501
    offset <- 0.1

    meshgrid <- function(a,b) {
        list(
            x=outer(b*0,a,FUN="+"),
            y=outer(b,a*0,FUN="+")
        )
    } 

    x <- seq(min(coord[,1])-offset, max(coord[,1])+offset, length=pp)
    y <- seq(min(coord[,2])-offset, max(coord[,2])+offset, length=pp)
    # the values of the probability density is now in the variable out
    #out <- dsn2.plot(x,y,dp=dp,nlevels=10,col="red")
    grid <- meshgrid(x, y)
    x_len <- as.vector(grid$x)
    y_len <- as.vector(grid$y)
    z_len <- dmsn(cbind(x_len,y_len), dp=dp)
    # Spatial grid distance
    dx <- x[2]-x[1]
    dy <- y[2] - y[1]

    # Interpolate the points from the regular grid to our grid points
    # this is the probability density for each data point
    obj<- list( x=x, y=y, z= matrix(z_len,nrow=pp,byrow=TRUE))
    pdensity <- interp.surface( obj, coord)
    fit <- list(coord=coord, pp=pp, x=x, y=y, x_len=x_len, y_len=y_len, z_len=z_len, dx=dx,
        dy=dy, pdensity=pdensity)
}
#all.fit <- NULL
#for (i in 1:length(files)) 
#    all.fit <- c(all.fit, list(mds.fit(files[i])))
fit <- mds.fit(filename)
#print(n)

#stop()

# for each data point, calculate the isocline of equal probability density
# and then integrate over all probability density grid points that
# lie inside the isocline
# 1- the integrated density is then the probability of each data point
# this is the most time consuming process
#pvals <- NULL
#pdensity[pdensity<1e-5] = 0
#print('Calculating contour lines...')
#for (i in 1:nrow(coord)) {
#pd <- fit$pdensity
getPval <- function(pd, fit) {
    #print(paste('density = ',pd,sep=''))
    mycontour <- contourLines(fit$x, fit$y, matrix(log10(fit$z_len),nrow=fit$pp,
        byrow=TRUE),levels=log10(pd+.Machine$double.eps))
    #browser()
	poly <- cbind(mycontour[[1]]$x,mycontour[[1]]$y)
	# if the contour is NOT closed
    if (!any(poly[1,]==poly[nrow(poly),])) {
        center = colMeans(poly)
        if (center[1] > 0) {
            if (center[2] > 0) # upper right
                poly <- rbind(poly, c(max(fit$x),min(fit$y)), 
                    c(min(fit$x),min(fit$y)), c(min(fit$x),max(fit$y)))
            else # lower right
                poly <- rbind(poly, c(min(fit$x),min(fit$y)), 
                    c(min(fit$x),max(fit$y)), c(max(fit$x),max(fit$y)))
        }
        else {
            if (center[2] > 0) # upper left
                poly <- rbind(poly, c(max(fit$x),max(fit$y)), 
                    c(max(fit$x),min(fit$y)), c(min(fit$x),min(fit$y)))
            else # lower left
                poly <- rbind(poly, c(min(fit$x),max(fit$y)), 
                    c(max(fit$x),max(fit$y)), c(max(fit$x),min(fit$y)))
        }
    }
    inside <- inout(cbind(fit$x_len,fit$y_len),poly, bound=TRUE)
	pvals <- sum(fit$z_len[!inside])*fit$dx*fit$dy # avoid float point subtraction
    #print(paste('p = ',pvals,sep=''))
}
#print('Calculating p-values...')
pvals <- laply(fit$pdensity, getPval, fit, .parallel=T)
#print(fit$pdensity[n])
#pvals <- getPval(fit$pdensity[n], fit)
#write.table(cbind(fit$coord,fit$pdensity,pvals), 
#    file=paste(outdir,prefix,'.dat',sep=''), quote=F, 
#    col.names=F, row.names=F)
write.table(pvals, 
    file=paste(outdir,prefix,'.dat',sep=''), quote=F, 
    col.names=F, row.names=F)
#}
#stop()
#for (i in 1:length(files))
    #pvals <- c(pvals, list(unlist(sge.parLapply(all.fit[[i]]$pdensity, getPval,
    #    all.fit[[i]], njobs=length(all.fit[[i]]$pdensity)))))

#Plot everything and calculalte the FDR corrected p values for every data point
#qvals <- p.adjust(pvals,method="BH")

#pdf('~/data/hgf/mds.pdf')
#df <- as.data.frame(coord)
#df$pvals <- pvals
#m <- ggplot(df, aes(x=V1,y=V2))
#print(m + geom_point() + geom_contour(aes(x=x_len,y=y_len,z=z_len,colour=..level..), 
#    bins=5) +
#    scale_colour_gradient("Density", low = "blue", high = "red") +
#    xlab("PC 1") + ylab("PC 2") +
#    coord_equal() +
#    opts(legend.position=c(0.2,0.85),
#    axis.title.x=theme_text(size=fontsize,vjust=0,hjust=0.5),title="MDS",
#    plot.title=theme_text(size=fontsize,vjust=1)))
#dev.off()

#pdf('~/data/yeast/response_dist.pdf')
#df <- as.data.frame(-log10(pvals))
#names(df) <- 'pvals'
#m <- ggplot(df, aes(x = pvals, y=..ndensity..))
#print(m + #geom_density() +
#    geom_histogram(colour = "darkgreen", fill = "white")+
#    #breaks=seq(0,max(df$pvals),length=31)) +
#    #binwidth=diff(range(df$pvals))/10) + 
#    scale_y_sqrt() + scale_x_sqrt() +
#    xlab("Response strength") + ylab("Density") + 
#    opts(axis.title.x = theme_text(size=fontsize, vjust = 0, hjust=0.5),
#    title=expression(paste("Yeast + ",alpha,"-factor",sep="")), 
#    plot.title=theme_text(size=fontsize,vjust=1)))
#dev.off()

