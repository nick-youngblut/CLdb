### start-up
rm(list=ls())

### packages
suppressPackageStartupMessages(library(Hmisc))
pkgs <- Cs(
	optparse,
	ape,
	genoPlotR
	)
for(i in 1:length(pkgs)){
	tmp <- pkgs[i]
	suppressPackageStartupMessages(library(pkgs[i], character.only=TRUE))
	}

### I/O
# initialize options
option_list <- list(
	make_option(c("-d", "--dna_segs"), type="character", help="dna_segs table"),
	make_option(c("-x", "--xlims"), type="character", help="xlims table"),
	make_option(c("-c", "--comparisons"), type="character", help="comparisons table [optional]"),
	make_option(c("-t", "--tree"), type="character", help="tree file plotted next to the loci [optional]"),	
	make_option(c("-o", "--outname"), type="character", help="Output file name. [default: modified dna_segs file name]"),
	make_option(c("-v", "--verbose"), action="store_false", default=TRUE, help="Print extra output")
	)
# get command line options, if help option encountered print help and exit, # otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))



#--- functions ---#
ext.edit <- function(file, ext){
	file <- gsub("\\.[^\\.]+$|$", ext, file, perl=TRUE)
	return(file)
	}

nwk2phylog <- function(file){
# loading a newick tree and converting to phylog format #
	nwk.str <- scan(file, what="complex", nmax=1, sep="\r")
	nwk.str <- gsub(" ", "_", nwk.str)
	tree <- newick2phylog(nwk.str)
	return(tree)
	}

edit.phylog.names <- function(phylog){
# phylog conversion edits names; this function edits back #
	## editing phylog names ##
	phylog$tre <- gsub("\\(X([0-9])", "(\\1", phylog$tre)
	phylog$tre <- gsub(",X([0-9])", ",\\1", phylog$tre)
	
	### leaves 
	names(phylog$leaves) <- gsub("^X([0-9])", "\\1", names(phylog$leaves))
	
	### nodes 
	names(phylog$nodes) <- gsub("^X([0-9])", "\\1", names(phylog$nodes))
	
	### parts
	for(i in names(phylog$parts)){
		phylog$parts[[i]] <- gsub("^X([0-9])", "\\1", phylog$parts[[i]])
		}
	names(phylog$parts) <- gsub("^X([0-9])", "\\1", names(phylog$parts))
	
	### paths
	for(i in names(phylog$paths)){
		phylog$paths[[i]] <- gsub("^X([0-9])", "\\1", phylog$paths[[i]])
		}
	names(phylog$paths) <- gsub("^X([0-9])", "\\1", names(phylog$paths))
	
	### phylog$droot
	names(phylog$droot) <- gsub("^X([0-9])", "\\1", names(phylog$droot))
	
	return(phylog)
	}

format.dna.segs.xlims <- function(x){
# formating the xlims table #
	x$dna_segs_id <- gsub(" ", "_", x$dna_segs_id)
	x$dna_segs_id <- gsub("\\.", "_", x$dna_segs_id)
	x$dna_segs_id <- gsub("-", "_", x$dna_segs_id)
	x$start <- as.numeric(x$start)
	x$end <- as.numeric(x$end)
	return(x)
	}

df2list.dna_segs <- function(x){
# converting dna_segs to a list by dna_segs_id #
	## convert to list ##
	x.l <- list()

	u.dna_segs_id <- unique(x$dna_segs_id)
	for(i in 1:length(u.dna_segs_id) ){
		# dna_segs #
		ii <- as.character(u.dna_segs_id[i])
		x.l[[ii]] <- dna_seg(x[x$dna_segs_id==ii, ])
		}
	
	return(x.l)
	}
	
df2list.xlims <- function(x){
# converting dna_segs to a list by dna_segs_id #
	## convert to list ##
	x.l <- list()

	u.dna_segs_id <- unique(x$dna_segs_id)
	for(i in 1:length(u.dna_segs_id) ){
		ii <- as.character(u.dna_segs_id[i])
		se <- x[ x$dna_segs_id==ii, c(1,2)]
		x.l[[ii]] <- as.vector(as.matrix(se[1,]))
		}
	
	return(x.l)
	}
	
edit.compare.names <- function(x){
# editing names in comparison table #
	x$dna_segs_id1 <- gsub(" ", "_", x$dna_segs_id1)
	x$dna_segs_id1 <- gsub("\\.", "_", x$dna_segs_id1)
	x$dna_segs_id1 <- gsub("-", "_", x$dna_segs_id1)
	return(x)
	}

df2list.compare <- function(compare, dna_segs.tbl, color_scheme="grey"){
# converting dna_segs to a list by dna_segs_id #
		# dna_segs.tbl <- DNAsegs.tbl
		# compare <- compare.tbl
	u.dna_segs_id <- unique(dna_segs.tbl$dna_segs_id)
	
	compare.l <- list();
	for(i in 1:(length(u.dna_segs_id)-1) ){
		# dna_segs.tbl #
		ii <- as.character(u.dna_segs_id[i])
		compare.l[[i]] <- compare[compare$dna_segs_id1==ii, ]
		
		if(nrow(compare.l[[i]]) == 0){
			compare.l[[i]][1,] <- c(rep(0, 4), 1, rep(NA, 5)) 
			}
	
		compare.l[[i]]$col <- apply_color_scheme(compare.l[[i]]$col, color_scheme=color_scheme)
		compare.l[[i]] <- comparison(compare.l[[i]])
		}
	return(compare.l)
	}
	

set.plot.size <- function(xlims, dna_segs){
# setting total plot size based on number of loci to plot #
	xlim.diff <- vector();
	for (i in names(xlims)){
		for (ii in seq(1,length(xlims[[i]]), 2)){
			xlim.diff <- append(xlim.diff, abs(xlims[[i]][ii] - xlims[[i]][ii+1]))
			}
		}
		
	height <- length(dna_segs)
	width <- log(max(xlim.diff), 1.8)
	
	return(list(width=width, height=height))
	}


#--- main ---#
# I/O error #
if(is.null(opt$dna_segs)){ stop("ERROR: provide a dna_segs file!") }
if(is.null(opt$xlims)){ stop("ERROR: provide an xlims file!") }
outfile <- ext.edit(opt$dna_segs, '.pdf')

# loading the pruned tree output by CLdb_dna_segs_orderByTree.pl #
## will be plotted with loci
## not needed
if(! is.null(opt$tree)){
	tree <- nwk2phylog(opt$tree)			# newick tree 
	tree <- edit.phylog.names(tree)				
	} else {
	tree = NULL
	}

# loading the dna_segs & xlims tables #
## both should be ordered by the tree if providing a tree
DNAsegs.tbl <- read.delim( opt$dna_segs )
xlims.tbl <- read.delim( opt$xlims )

# formatting the dna_segs & xlims tables for genoPlotR #
## editing names
DNAsegs.tbl <- format.dna.segs.xlims(DNAsegs.tbl)
xlims.tbl <- format.dna.segs.xlims(xlims.tbl)
## converting to lists
dna_segs <- df2list.dna_segs(DNAsegs.tbl)
xlims <- df2list.xlims(xlims.tbl)

# loading and formatting the comparisons table #
## not needed if no comparisions necessary 
if(! is.null(opt$comparisons)){
	compare.tbl <- read.delim(opt$comparisons)
	compare.tbl <- edit.compare.names(compare.tbl)
	compare <- df2list.compare(compare.tbl, DNAsegs.tbl)
	} else {
	compare = NULL
	}

# setting plot dimensions based on the number & length of loci #
p.dim <- set.plot.size(xlims, dna_segs)


# plotting #
## comment out any portions that were not provided (e.g. 'tree' or 'compare')
#quartz(width=p.dim[["width"]], height=p.dim[["height"]])		# plot size
pdf(file=outfile, width=p.dim[["width"]], height=p.dim[["height"]])
plot_gene_map(dna_segs=dna_segs,
		xlims=xlims,
		tree=tree,
		comparisons=compare,
		main=NULL,
		dna_seg_scale=TRUE,
		scale=FALSE,
		tree_branch_labels=0,
		seg_plot_height=20)
dev.off()
		
		
# status
message("loci plot written")
