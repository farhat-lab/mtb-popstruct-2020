#!/usr/bin/env Rscript

# parameters: <tree_file> <dir_vcf> <threshold_fst> <min_num_isolates_group> <tag_out_files>

args = commandArgs(trailingOnly=TRUE)
# this is useful to test the code interactively
#args = c("results/trees/lineage1_sensitive/tree_lin1_sensitive.tree","./vcf_test/","0.33","10","THX")

if(length(args) != 5){
  stop("please provide the following arguments: <tree_file> <path_to_vcf> <threshold_fst> <min_num_isolates_group> <tag_out_files>")
}

library("ape")
library("phytools")
library("data.tree")
library("phangorn")
library("PopGenome")
library("RColorBrewer")

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
## I get the original node_labeles
orig_node_lab = m_tree$node.label
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
list_data_all_nodes=list()

get_groups = function(node, min_num_isolates_group){
  cat("* analyzing node", node$name,"\n")
  vect_leaves1=c()
  vect_leaves2=c()
  count1=node$children[[1]]$leafCount
  count2=node$children[[2]]$leafCount
  vect_leaves1 = unname(sapply(node$children[[1]]$leaves, function(x){x$name}))
  vect_leaves2 = unname(sapply(node$children[[2]]$leaves, function(x){x$name}))
  list_data_all_nodes[[node$name]] <<- list(node$name,vect_leaves1,vect_leaves2)
  if((count1 > min_num_isolates_group) && (count2 > min_num_isolates_group)){
    list_comparisons[[node$name]] <<- list(node$name,vect_leaves1,vect_leaves2)
  }
}

dt$Do(get_groups, min_num_isolates_group = min_n_isol_group, filterFun = isNotLeaf)
# I write down the list with all nides and their descendents
save(list_data_all_nodes,file=paste(tag_out,"_list_data_all_nodes.Rdata",sep=""))

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
vect_fst_values=unlist(data_fst)
udata=data.frame(node=names(vect_fst_values), fst=vect_fst_values)
write.csv(udata, file = paste(tag_out,"_table_fst_values.tsv",sep=""), quote=F, row.names=F)

cat("* I generate the FigTree color annotation file\n")
# I get the list of the nodes that have significant Fst
nodes_fst_above_thr = names(vect_fst_values[vect_fst_values > thr_fst])
# Then I count their children and I rank them by the number of children
num_children = c()
for(node in nodes_fst_above_thr){
  # I count the number of children
  num_children = c(num_children, length(list_comparisons[[node]][[2]]) +  length(list_comparisons[[node]][[3]]))
}
names(num_children) = c(nodes_fst_above_thr)
ranked_nodes = sort(num_children, decreasing = TRUE)
# Now I apply the colors
## I get the leaves of the tree
vect_leaves_tree = m_tree$tip.label
num_leaves = length(vect_leaves_tree)
df_colors = data.frame("isolate" = vect_leaves_tree, "!color" = rep("#000000",num_leaves), check.names = FALSE, stringsAsFactors = FALSE)
## I generate an appropriate palette of colors
### I use RColorBrewer
coul = brewer.pal(8, "Dark2")
### I extend the palette -- thanks to: https://www.r-graph-gallery.com/40-rcolorbrewer-get-a-longer-palette/
#### How many colors do I need?
num_colors = 2 * length(ranked_nodes)
colors = colorRampPalette(coul)(num_colors)

# I assign the colors to the leaves and I generate at the same time a data frame that cointains the lists of the children of each node having F_st > thr_fst
num_nodes_fst_greater_thr = length(names(ranked_nodes))
df_groups = data.frame("node" = names(ranked_nodes), "children_1" = rep("-", num_nodes_fst_greater_thr), "children_2" = rep("-", num_nodes_fst_greater_thr), stringsAsFactors= FALSE)

counter = 1
names_pie = rep("-", length(colors))
for(node in names(ranked_nodes)){
color_1 = colors[counter]
names_pie[counter] = paste(node,"_c1")
color_2 = colors[(length(colors)/2)+counter]
names_pie[(length(colors)/2)+counter] = paste(node,"_c2")

isolates_1 = list_comparisons[[node]][[2]]
isolates_2 = list_comparisons[[node]][[3]]
# I select the lines that contain isolates that are children of "node" and I assign the colors
df_colors[df_colors$isolate %in% isolates_1,]$"!color" = rep(color_1, length(df_colors[df_colors$isolate %in% isolates_1,"!color"]))
df_colors[df_colors$isolate %in% isolates_2,]$"!color" = rep(color_2, length(df_colors[df_colors$isolate %in% isolates_2,"!color"]))
counter = counter + 1

# I work on the data frame nodes -> children (the tree is a binary tree)
df_groups[df_groups$node==node, "children_1"] = paste(list_comparisons[[node]][[2]], collapse = ",")
df_groups[df_groups$node==node, "children_2"] = paste(list_comparisons[[node]][[3]], collapse = ",")
}
# I save the annotation files (isolates -> colors)
write.table(df_colors, file = paste(tag_out,"_color_annotation_FigTree.tsv",sep=""), sep ="\t", quote = FALSE, row.names = FALSE)
write.table(df_groups, file = paste(tag_out,"_nodes_fst_above_thr_children.tsv",sep=""), sep ="\t", quote = FALSE, row.names = FALSE)

#### I plot the palette
pdf(file=paste(tag_out,"_palette.pdf",sep=""))
pie(rep(1,length(colors)), col = colors , main="", labels = names_pie)
dev.off()

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
df_nodes = data.frame("node" = 1:length(out_tree2$node.label), "fst_greater_thr" = out_tree2$node.label)
write.table(df_nodes, file = paste(tag_out,"_node_annotation_FigTree.tsv",sep=""), sep ="\t", quote = FALSE, row.names = FALSE)

## tree #3 -- tree with the original node numbers
write.tree(m_tree,file=paste(tag_out,"_rooted_node_numbers.tree",sep=""))

## tree #4 -- nexus tree using metacomments
out_tree4=m_tree
vect_new_labels = c()
for(i in 1:length(orig_node_lab)){
  if(i %in% set_to_one){
    current_label = paste(out_tree4$node.label[i],"[&fst_greater_thr=1,support=", orig_node_lab[i],"]", sep="" )
  }
  else{
    current_label = paste(out_tree4$node.label[i],"[&fst_greater_thr=0,support=", orig_node_lab[i],"]", sep="" )
  }
vect_new_labels = c(vect_new_labels, current_label)
}
out_tree4$node.label = vect_new_labels
print(out_tree4$node.label)
writeNexus(out_tree4,file=paste(tag_out,"_rooted.nex",sep=""))
