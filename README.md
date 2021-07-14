## Background

This repository includes two scripts (`define_groups.R` and `define_groups_colorize.R`) that were used to explore the phylogeny and guide the definitions of the sub-lineages. 

Such scripts take a binary tree in Newick format and a combined  `.vcf`  file as an input, determine the FST (Fixation Index) values for each of the nodes of the tree (the user can also decide a threshold `min_num_isolates_group` to stop this process) and return a tree which is rooted (midpoint rooting) and annotated with the FST values, that can be visualized with [FigTree](http://tree.bio.ed.ac.uk/software/figtree/).

Note:  a combined  `.vcf`  file is a `.vcf` file containing variant calls for all the samples present in the tree

What differs between the two scripts are the outputs:  

* `define_groups.R` returns (1) a Newick tree, rooted and annotated with all the calculated FST values; (2) a rooted Newick tree where the FST values are binarized, according to a threshold decided by the user; 
* `define_groups_colorize.R` returns (1) a Newick tree, rooted and annotated with all the calculated FST values; (2) a rooted Newick tree where the FST values are binarized, according to a threshold decided by the user; (3) a rooted Newick tree with the original node numbers; (4) a Nexus tree with meta-comments that can be visualized with FigTree.



## Usage

The general syntax to execute the scripts is the following:

```
define_groups_colorize.R [orig_newick_tree] [dest_folder] [fst_threshold] [min_num_isolates_group] 
```

Note: `define_groups.R` has the same syntax, so it is sufficient to change the name of the script.
