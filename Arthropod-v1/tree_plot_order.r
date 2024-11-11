library(ggplot2)
library(ggtree)
library(treeio)
args <- commandArgs(trailingOnly = TRUE)
tree = args[1]
group = args[2]
output_prefix = args[3]

tree <- read.tree(tree)  # 读取树文件
group_info <- read.delim(group, sep='\t', row.names = 1) # 读取物种分类文件
#colnames(group_info)=c("Group")
# 将物种分类信息与tree文件合并
groupInfo <- split(row.names(group_info), group_info$order)
tree1 <- groupOTU(tree, groupInfo)

# 创建颜色向量
color_vector <- setNames(group_info$order_color, group_info$order)

# rectangular_tree_plot
pdf(paste(output_prefix, ".rectangular_tree.pdf", sep = ""))
ggtree(tree1,layout = "rectangular", aes(color=group), size=0.1) +
    geom_tiplab(aes(color=group), size=0.6, hjust=-0.01) +
    theme(legend.position = "none") + 
    xlim(NA, 80) + 
    scale_color_manual(values=color_vector)
dev.off()

# circular_tree_plot
pdf(paste(output_prefix, ".circular_tree.pdf", sep = ""))
ggtree(tree1,layout = "circular", aes(color=group), size = 0.1) +
    geom_tiplab2(aes(color=group), size=0.6, hjust=-0.01) + 
    xlim(NA, 80) +
    theme(legend.position = "none") + 
    scale_color_manual(values=color_vector)
dev.off()

