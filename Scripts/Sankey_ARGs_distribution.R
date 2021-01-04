######Distribution of resistance genes at different taxonomic levels in bacteria####

table <- read.delim('Bacteria_level_for_Sankey.xls', stringsAsFactors = FALSE, check.names = FALSE)
#Integrate the nested relation of classification and ARGs number to construct the link list of Sankey graph
genus_species <- table[c('Genus','Species','Freq')]
names(genus_species) <- c('source', 'target', 'Freq')
family_genus <- aggregate(table$Freq, by = list(table$Family, table$Genus), FUN = sum)
names(family_genus) <- c('source', 'target', 'Freq')
order_family <- aggregate(table$Freq, by = list(table$Order, table$Family), FUN = sum)
names(order_family) <- c('source', 'target', 'Freq')
class_order <- aggregate(table$Freq, by = list(table$Class, table$Order), FUN = sum)
names(class_order) <- c('source', 'target', 'Freq')
phylum_class <- aggregate(table$Freq, by = list(table$Phylum, table$Class), FUN = sum)
names(phylum_class) <- c('source', 'target', 'Freq')
kingdom_phylum <- aggregate(table$Freq, by = list(table$Kingdom, table$Phylum), FUN = sum)
names(kingdom_phylum) <- c('source', 'target', 'Freq')

link_list <- rbind(kingdom_phylum,phylum_class, class_order, order_family, family_genus,genus_species)

#Build the node list and assign ID  to the category names in the Link list
node_list <- reshape2::melt(table, id = 'Freq')
node_list <- node_list[!duplicated(node_list$value), ]
head(node_list)

link_list$IDsource <- match(link_list$source, node_list$value) - 1 
link_list$IDtarget <- match(link_list$target, node_list$value) - 1
head(link_list)

#plot sankey graph by networkD3 package
library(networkD3)
p <- sankeyNetwork(Links = link_list, Nodes = node_list,
                   Source = 'IDsource', Target = 'IDtarget', Value = 'Freq', 
                   NodeID = 'value', NodeGroup = 'variable', 
                   fontSize = 12, sinksRight = FALSE,width = 1200,height = 700)

p
#save as .html
saveNetwork(p,"Bacteria_level_Sankey.html")

library(webshot)
webshot("Bacteria_level_Sankey.html","Bacteria_level_Sankey.pdf",vwidth = 1200,vheight=700)
