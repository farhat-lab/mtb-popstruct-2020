#!/usr/bin/env Rscript

# parameters: <tree_file> <dir_vcf> <threshold_fst> <min_num_isolates_group> <tag_out_files>

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 5){
  stop("please provide the following arguments: <tree_file> <path_to_vcf> <threshold_fst> <min_num_isolates_group> <tag_out_files>")
}

library("ape")
library("phytools")
library("data.tree")
library("phangorn")
library("PopGenome")

cat("* I parse the arguments\n")
tre = args[1]
cat("  - tree: ",tre,"\n",sep="")
path_vcf = args[2]
cat("  - vcf: ",path_vcf,"\n",sep="")
thr_fst = as.numeric(args[3])
cat("  - threshold_fst: ",thr_fst,"\n",sep="")
min_n_isol_group = as.numeric(args[4])
cat("  - min_n_isolates_group: ",min_n_isol_group,"\n",sep="")
tag_out = args[5]
cat("  - tag_out_files: ",tag_out,"\n",sep="")


# I read the tree
cat("* I load the tree\n")
tree_all = read.tree(tre)
## I do the midpoint rooting
m_tree=midpoint.root(tree_all)
write.tree(m_tree,file=paste(tag_out,"_rooted.tree",sep=""))
## I assegna label to each node
m_tree$node.label=seq(1:length(m_tree$node.label))

#I check if the tree is binary
if(! is.binary(m_tree)){
  stop(cat("* the tree you provided is not binary\n"))
}

# I convert the tree to a data.tree object
cat("* I convert the tree to a data.tree object\n")
dt=as.Node(m_tree)

cat("* I load the genetic data (VCF)\n")
# I load the vcf
gpop_data <- readData(path_vcf, format="VCF", progress_bar_switch=FALSE)


cat("* I detect the nodes where I will calculate the F_st values\n")
list_comparisons=list()
get_groups = function(node, min_num_isolates_group){
  cat("* analyzing node", node$name,"\n")
  vect_leaves1=c()
  vect_leaves2=c()
  count1=node$children[[1]]$leafCount
  count2=node$children[[2]]$leafCount
  if((count1 > min_num_isolates_group) && (count2 > min_num_isolates_group)){
    vect_leaves1 = unname(sapply(node$children[[1]]$leaves, function(x){x$name}))
    vect_leaves2 = unname(sapply(node$children[[2]]$leaves, function(x){x$name}))
    list_comparisons[[node$name]] <<- list(node$name,vect_leaves1,vect_leaves2)
  }
}

dt$Do(get_groups, min_num_isolates_group = min_n_isol_group, filterFun = isNotLeaf)

calc_fst=function(list_nodes_children, genpop_data){
  list_node_fst=list()
  for(i in 1:length(list_nodes_children)){
    cat("* calculating F_st for node", list_nodes_children[[i]][[1]],"\n")
    genpop_data  <- set.populations(genpop_data, list(list_nodes_children[[i]][[2]],list_nodes_children[[i]][[3]]))
    cat("\n")
    fst_data <- PopGenome::F_ST.stats(genpop_data,FAST = TRUE, mode = "nucleotide")
    cat("\n")
    list_node_fst[[list_nodes_children[[i]][[1]]]]=fst_data@nucleotide.F_ST

  }
  return(list_node_fst)
}

cat("* I calculate the F_st values\n")
# I calculate the Fst at the nodes of interest
data_fst=calc_fst(list_comparisons, gpop_data)

#I build the vectors with the information to display in the final tree.
vect_fst_values=unlist(data_fst)
pos_to_modify=as.numeric(names(vect_fst_values))

cat("* I write down the final trees\n")
## tree #1 -- FST values
out_tree1=m_tree
out_tree1$node.label=rep("-1",length(out_tree1$node.label))
out_tree1$node.label[pos_to_modify]=unname(vect_fst_values)
write.tree(out_tree1,file=paste(tag_out,"_rooted_fst_raw.tree",sep=""))

## tree #2 -- binarized values
out_tree2=m_tree
out_tree2$node.label=rep("0",length(out_tree2$node.label))
set_to_one=as.numeric(names(vect_fst_values[vect_fst_values > thr_fst]))
out_tree2$node.label[set_to_one]=1
write.tree(out_tree2,file=paste(tag_out,"_rooted_fst_binary.tree",sep=""))
