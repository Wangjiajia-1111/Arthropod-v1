library(ggplot2)
library(ggtree)
library(treeio)
args <- commandArgs(trailingOnly = TRUE)
tree = args[1]
group = args[2]
output_prefix = args[3]

tree <- read.tree(tree)  # 读取树文件
group_info <- read.delim(group, sep='\t', row.names = 1, header = F) # 读取物种分类文件
# 将物种分类信息与tree文件合并
groupInfo <- split(row.names(group_info), group_info$V2)
tree <- groupOTU(tree, groupInfo)

# rectangular_tree_plot
pdf(paste(output_prefix, ".rectangular_tree.pdf", sep = ""))
ggtree(tree,layout = "rectangular",color="black",size=0.2) +
    geom_tiplab(aes(color=group), size=0.53, hjust=-0.01) +
    theme(legend.position = "none") + 
    xlim(NA, 100) +
    scale_color_manual(values=c("0"="black","1"="red"))
dev.off()

# circular_tree_plot
pdf(paste(output_prefix, ".circular_tree.pdf", sep = ""))
ggtree(tree,layout = "circular",color="black",size = 0.2) +
    geom_tiplab2(aes(color=group), size=0.55, hjust=-0.01) + 
    xlim(NA, 100) +
    theme(legend.position = "none") +
    scale_color_manual(values=c("0"="black","1"="red"))
dev.off()

