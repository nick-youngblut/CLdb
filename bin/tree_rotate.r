# Script for rotating internal tree nodes to order tips in a particular manner

#--- start-up ---#
rm(list=ls())


#--- packages ---#
suppressPackageStartupMessages(require(Hmisc))
pkgs <- Cs(
  optparse,
  ape
)
for(i in 1:length(pkgs)){
  tmp <- pkgs[i]
  suppressPackageStartupMessages(require(pkgs[i], character.only=TRUE))
}

#--- Option parsing ---#
option_list <- list(
  make_option(c("-t", "--tree"), type="character", help="Tree file"),
  make_option(c("-f", "--format"), type="character", default="newick", help="Tree file format (newick or nexus)"),
  make_option(c("-l", "--leaves"), type="character", help="File listing desired order of tree leaves; 1 per line ('-' if STDIN)"),
  make_option(c("-o", "--outname"), type="character", help="Output file name. [default: modified input file name]"),
  make_option(c("-v", "--verbose"), action="store_false", default=TRUE, help="Print extra output"),
  make_option(c("-d", "--description"), action="store_true", default=FALSE, help="Script description")
)
# get command line options, if help option encountered print help and exit
## otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))


#--- I/O error check ---#
# Description
if(opt$description == TRUE){
  message("
Rotating internal tree nodes to match the user-defined leaf order.

Rotating the tree to completely match the user-defined leaf order may not be possible 
depending on the defined leaf order.

All possible combinations of internal nodes (only those that affect leaf order)
will be rotated and the tree with closest match to the user-defined order will be written. 

Spearman's rank correlation is used to compare the leaf ordering of each rotated tree 
vs the user-defined ordering. The tree with the best correlation is written.
")
  quit(save='no')
}

# File extension
ext.edit <- function(file, ext){
  file <- gsub("\\.[^\\.]+$|$", ext, file, perl=TRUE)
  return(file)
}

if(grepl("^new", opt$format, ignore.case=TRUE, perl=TRUE) == TRUE ) { ext <- ".nwk"} else
  if(grepl("^nex", opt$format, ignore.case=TRUE, perl=TRUE) == TRUE ) { ext <- ".tre"} else{
    stop(" ERROR: tree format must be nexus or newick")
  }

# Required opts/defaults
if(is.null(opt$tree)){ stop(" ERROR: provide a tree file (-t)")}
if(is.null(opt$leaves)){ stop(" ERROR: provide a leaves file (-l)")}
if(is.null(opt$outname)){ opt$outname <- ext.edit(opt$tree, paste(c("_rot", ext), collapse="")) }



#--- Functions ---#
# FUNC: number of possible rotation combinations for provided tip ordering
tree.rotate.ncomb = function(tree, x){
  # tree = phylo object
  # x = vector of new leaf ordering (based on leaf labels)
  
  # user-defined leaf ordering (by leaf IDs, not leaf labels)
  tip.index = as.data.frame(list(index = 1:length(tree$tip)), row.names=tree$tip )
  user.tip.order = tip.index[x, ]
  
  # determining which taxa have changed order
  new.order = which( tree$tip.label != x )
  if(length(new.order) == 0){
    message("No rotation needed!")
    return(tree)
  }
  
  # determine most topologically distant of taxa
  ## matrix of all leaf-leaf distances
  leaf.dists = cophenetic(compute.brlen(tree))
  ## all combintations of focal leaves
  cmb = combn(new.order, 2)
  ## pairwise distances of just those leaves
  p.dists = apply(cmb, 2, function(x, dists=leaf.dists){
    return( leaf.dists[x[1], x[2]] )
  })
  ## pair with max distance
  leaves.max.dist = cmb[, which(p.dists == max(p.dists))]
  
  # determine number of nodes in tree of lca for max-dist pair
  ## getting lca
  lca = getMRCA(tree, leaves.max.dist)  
  ## getting all internal node children of lca
  int.node.ids = getDesc(tree, lca)
  
  # getting all possible combinations
  #int.node.ids = unique(subtree$edge[,1])
  if( length(int.node.ids) > 1 ){
    combs = sapply(1:length(int.node.ids), function(x){
      combn(int.node.ids, x)
    })
    # making sure it it a list (n = 2 = list)
    if(is.matrix(combs)==TRUE){  # if just matrix, making list
      combs = list(combs)
    }
  } else {   # n = 1
    combs = list()
  }
  ## adding just single node rotation
  combs[[length(combs) + 1]] = matrix(int.node.ids, nrow=1)
  
  comb.count = length(unlist(combs))
  
  return(comb.count)
}


# FUNC: getting children nodes from a phylo object
getDesc = function(tree, nodes, desc=c(), leaves=FALSE){
  # tree = phylo object
  # nodes = basal node
  # desc = leave blank
  # leaves = include leaves
  # appending to list of descendents (includes 1st provided node)
  desc = append(desc, nodes)
  # getting all child nodes  
  child.nodes = tree$edge[tree$edge[,1]==nodes, 2]
  # determining which child nodes are internal
  if(leaves==FALSE){
    int.nodes = tree$edge[tree$edge[,1]==child.nodes, 1]
    child.nodes = int.nodes
  }  
  if( length(child.nodes) > 0 ){
    getDesc(tree, nodes=child.nodes, desc=desc)
  } else {
    return(desc)
  }
}


# FUNC: rotation of nodes to get defined ordering
rotate.tree = function(tree, x){
  # tree = phylo object
  # x = vector of new leaf ordering (based on leaf labels)
  
  # user-defined leaf ordering (by leaf IDs, not leaf labels)
  tip.index = as.data.frame(list(index = 1:length(tree$tip)), row.names=tree$tip )
  user.tip.order = tip.index[x, ]
  
  # determining which taxa have changed order
  new.order = which( tree$tip.label != x )
  
  # determine most topologically distant of taxa
  ## matrix of all leaf-leaf distances
  leaf.dists = cophenetic(compute.brlen(tree))
  ## all combintations of focal leaves
  cmb = combn(new.order, 2)
  ## pairwise distances of just those leaves
  p.dists = apply(cmb, 2, function(x, dists=leaf.dists){
    return( leaf.dists[x[1], x[2]] )
  })
  ## pair with max distance
  leaves.max.dist = cmb[, which(p.dists == max(p.dists))]
  
  # determine number of nodes in tree of lca for max-dist pair
  ## getting lca
  lca = getMRCA(tree, leaves.max.dist)  
  ## getting all internal node children of lca
  int.node.ids = getDesc(tree, lca)
  
  # getting all possible combinations
  #int.node.ids = unique(subtree$edge[,1])
  if( length(int.node.ids) > 1 ){
    combs = sapply(1:length(int.node.ids), function(x){
      combn(int.node.ids, x)
    })
    # making sure it it a list (n = 2 = list)
    if(is.matrix(combs)==TRUE){  # if just matrix, making list
      combs = list(combs)
    }
  } else {   # n = 1
    combs = list()
  }
  ## adding just single node rotation
  combs[[length(combs) + 1]] = matrix(int.node.ids, nrow=1)
  
  
  #--- running through combinations to determine closest to correct order ---#
  ## function for getting tree with max 'rot.sim' value
  get.max.rot.sim.tree = function(x){
    # x = list of trees 
    # max rot.sim value among trees
    max.rot.sim = max(
      unlist(
        lapply(x, function(x){ x$rot.sim })
      )
    )
    
    # getting tree with max value
    tree.max.rot.sim = function(x, max.rot.sim){
      for( i in 1:length(x) ){
        if( x[[i]]['rot.sim'] == max.rot.sim ){
          return(x[[i]])
        }
      }
    }
    max.rot.sim.tree = tree.max.rot.sim(x, max.rot.sim)
    return(max.rot.sim.tree)
  }
  
  # getting best trees for each n-combs (best = closest to correct ordering)
  max.rot.sim.trees = lapply(combs, tree=tree, user.tip.order=user.tip.order, FUN = function(x, tree, ...){    
    ret2 = apply(x, 2, tree=tree, FUN = function(y, tree, ...){
      # applying all rotations for node combination
      tree.rot = tree    
      for (node in y){
        tree.rot = rotate(tree.rot, node)
      }
      
      # compare to non-rotated tree (range: 0-1, 0 = no matching; 1 = all matching)
      ## get tip ordering 
      tree.rot.order = setdiff(tree.rot$edge[,2], tree.rot$edge[,1])
      
      
      ## sanity check
      if(length(user.tip.order) != length(tree.rot.order)){
        stop("user.tip.order != tree.rot.order")
      }
      
      # ranking distance between user-defined and rotated tree tips
      #rot.sim = length( which( user.tip.order == tree.rot.order) ) /
      #  length(tree.rot.order)   # fraction that match in their order
      rot.sim = spearman(user.tip.order, tree.rot.order)
      tree.rot$rot.sim = rot.sim
      return(tree.rot) 
    })
    # ret2 = list    
    
    # getting max rotation similarity
    max.rot.sim.tree = get.max.rot.sim.tree(ret2)
    
    #print(max.rot.sim.tree)        
    #stop('#--- done ---#')
    
    return(max.rot.sim.tree)
  })
  
  # finding max tree
  best.rot.tree = get.max.rot.sim.tree(max.rot.sim.trees)
  
  # status
  message("Spearman's rho for best tree's order vs user-defined order: ", round(best.rot.tree$rot.sim, 3))
  
  return(best.rot.tree)
}


#--- MAIN ---#
# I/O
## tree
if(ext == ".nwk"){ tree <- read.tree(opt$tree) } else
if(ext == ".tre"){ tree <- read.nexus(opt$tree) }

## leaves
if(opt$leaves == "-"){
  user.leaves = read.delim(pipe('cat /dev/stdin'), header=F)
} else {
  user.leaves = read.delim(opt$leaves, header=F)
}
user.leaves = as.character(user.leaves[,1])  # vector of leaves

# sanity checks
## check that user.leaves is proper length
if( length(user.leaves) != length(tree$tip.label)){
  stop("'--leaves' must be same length as tree leaves!")
}
not.in = user.leaves[! user.leaves %in% tree$tip.label]
if( length(not.in) > 0 ){
  for(i in not.in){
    message("Error: '", i, "'' not in provided tree")
  }
  stop("'--leaves' must match the provided tree!")
}

# status
ptm <- proc.time()[1]


# n-combinatinos
n.comb = tree.rotate.ncomb(tree, user.leaves)
message("Number of combinations to try: ", n.comb)

# best rotated tree
best.rot.tree = rotate.tree(tree, user.leaves)

# writing best tree
if(ext == ".nwk"){ write.tree(best.rot.tree, file=opt$outname) } else
if(ext == ".tre"){ write.nexus(best.rot.tree, file=opt$outname) }

# status
message("Best tree file written: '", opt$outname, "'")
message("Time to complete (system time): ", round(proc.time()[1] - ptm, 3), " secs\n")

