# 05/07/2015
# Dom Bennett
# Plotting trees

# LIBS
library (ape)

# FUNCTIONS
addSupportIndicators <- function (support.thresh=95, pch=19, col='black') {
  # Add node dots to indicate support value
  nodes <- which(as.numeric (tree$node.label) > support.thresh) + length (tree$tip.label)
  nodelabels(node=nodes, pch=pch, col=col)
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
  i <- which (tree$tip.label == 'sample')
  tree$tip.label[i] <- 'Dormouse_sample'
  tip.colours <- rep ('black', length (tree$tip.label))
  tip.colours[i] <- 'firebrick3'
  # plot
  filename <- paste0 (sub ('\\.tre', '', f), '.pdf')
  pdf (file.path (wd, filename), h=21, w=14)
  plot (tree, tip.color=tip.colours)
  addSupportIndicators (support.thresh=75, pch=19, col='grey65')
  addSupportIndicators (support.thresh=90, pch=19, col='grey30')
  addSupportIndicators ()
  axisPhylo()
  dev.off()
}