# 05/07/2015
# Dom Bennett
# Plotting trees

# LIBS
library (ape)

# FUNCTIONS
addSupportIndicators <- function (support.thresh=95, pch=19, ...) {
  # Add node dots to indicate support value
  nodes <- which(as.numeric (tree$node.label) > support.thresh) + length (tree$tip.label)
  nodelabels(node=nodes, pch=19, ...) 
}

removeUnsupportedNodes <- function (tree, support.thresh=50) {
  # Return tree with unsupported nodes removed
  # get all nodes with less than 50% support
  unresolved <- which (as.numeric (tree$node.label) < support.thresh)
  for (node in unresolved) {
    # get supporting edge of node
    sedge <- which (tree$edge[ ,2] == node + length (tree$tip.label))
    # get desecending edges
    dedge <- which (tree$edge[ ,1] == node + length (tree$tip.label))
    # add length of sedge to dedge
    tree$edge.length[dedge] <- tree$edge.length[dedge] + tree$edge.length[sedge]
    # set sedge length to 0
    tree$edge.length[sedge] <- 0
  }
  # drop unsupported nodes
  tree <- di2multi (tree)
  return(tree)
}


# PROCESS
wd <- 'final_trees'
treefiles <- list.files (wd, pattern='\\.tre')
for (f in treefiles) {
  tree <- read.tree(file.path (wd, f))
  print(f)
  print(length(tree$tip.label))
  tree <- removeUnsupportedNodes (tree)
  # re-root and drop outgroup
  if ('outgroup' %in% tree$tip.label) {
    tree <- root (tree, 'outgroup')
    tree <- drop.tip (tree, 'outgroup')
  } else {
    # node 76 is the node that matches with Haukisalmi et al.
    tree <- root (tree, node=76)
  }
  # rename and colour sample
  tip.labels <- sub("__", " (", tree$tip.label)
  tip.labels <- sub("_", " ", tip.labels)
  tip.labels <- sub("$", ")", tip.labels)
  tip.labels <- gsub("_", ", ", tip.labels)
  tree$tip.label <- tip.labels
  i <- which (grepl('sample', tree$tip.label))
  tree$tip.label[i] <- '"2009 captive dormouse sample"'
  tip.sizes <- rep (1.5, length (tree$tip.label))
  tip.sizes[i] <- 1.7
  tip.cols <- rep ('gray15', length (tree$tip.label))
  tip.cols[i] <- "black"
  # plot
  filename <- paste0 (sub ('\\.tre', '', f), '.pdf')
  pdf (file.path (wd, filename), h=30, w=14)
  plot (tree, cex=tip.sizes, lwd=2, tip.color=tip.cols)
  addSupportIndicators (support.thresh=75, cex=1.5)
  addSupportIndicators (support.thresh=75, col="white")
  addSupportIndicators (support.thresh=90, cex=1.5)
  addSupportIndicators (support.thresh=90, col='grey')
  addSupportIndicators (cex=1.5)
  axisPhylo()
  dev.off()
}