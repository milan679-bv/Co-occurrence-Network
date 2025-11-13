
library(phyloseq)
library(plyr)
library(igraph)
library(ggplot2)
library(Hmisc)

otufile="HE_network/bac_rare_p.biom"
otufile="HE_network/fungi_rare_rel_ab.biom"

mapfile="HE_network/mapping_network.txt"

fungi_biom = import_biom(otufile) ## import of OTU table ##
colnames(tax_table(fungi_biom)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") ## setting up taxonomic levels ##
fungi_map<-import_qiime_sample_data(mapfile) ## import of metadata ##
fungi_raw<-merge_phyloseq(fungi_biom, fungi_map) ## creats main working file - OTu table with metadata ##


physeq_fungi_filter= prune_taxa(taxa_sums(fungi_raw)>50,fungi_raw)#1552 OTUs, filter taxa whose abundance is lower than 5 reads per sample
fungi_filter = filter_taxa(fungi_raw, function(x) sum(x > 3) > (0.05*length(x)), TRUE) # Remove taxa not seen more than 3 times in at least 20% of the samples.
fungi_filter=genefilter_sample(fungi_raw, filterfun_sample(function(x) x >= 3),A = 0.05*nsamples(fungi_raw))

fungi_Genus=tax_glom(fungi_raw, "Genus")
cryoOM_genus=subset_samples(fungi_Genus, horizon=="cryoOM")
top100=names(sort(taxa_sums(cryoOM_genus), TRUE)[1:200])
top100_prune=prune_taxa(top100, cryoOM_genus)

cryoOM_genus_filter <- prune_taxa(taxa_sums(cryoOM_genus)>55,cryoOM_genus) #38 (topsoil) or 14 (cryoOM) comes from the sum(taxa_sums(cryoOM_genus) = 38000 or 14000 which is 100%), so to filter at 0.01% we need to 14*100/14000=0.01% 


otu_table_fungi <- data.frame(otu_table(top100_prune))  #get otu table
taxa_table_fungi <- data.frame(tax_table(top100_prune),stringsAsFactors = FALSE)  #get taxa table
sample_data_fungi <- data.frame(sample_data(top100_prune))  #sample data/mapping file

taxa_table_fungi$abundance <- rowSums(data.frame(otu_table(top100_prune)))/sum(taxa_sums(cryoOM_genus))*100  #calculate relative abundance, 75000 is number of otus presents in the sample.
sum(taxa_table_fungi$abundance) #92.19 % for cryoOM and 92.31 for topsoil

otu_table_fungi <- data.frame(t(otu_table_fungi))

#cryoOM
otu_table_bac_cryoOM <- otu_table_fungi
sort(colSums(otu_table_bac_cryoOM))

bac_cryoOM<-as.matrix(otu_table_bac_cryoOM)
srcm1<-rcorr(bac_cryoOM, type="spearman")   #pairwise correlations
rval<-abs(srcm1$r)
rcrit<-0.6

Pcrit<-0.01

Pval<-srcm1$P
Padj<-p.adjust(Pval, "fdr")
correls_cryoOM <- as.data.frame.table(t(ifelse(test=rval>rcrit&rval<1&Padj<Pcrit,yes=srcm1$r,no="0"))) #to filter the correlations according to the cutoffs
correls_cryoOM1 <-correls_cryoOM[which(correls_cryoOM$Freq != 0),] 
length(which(correls_cryoOM$Freq !=0))

hist(srcm1$r)  #to look at the distribution of coeffients

network1 <- graph_from_data_frame(d=correls_cryoOM1,directed = F)  #creat igraph object from data frame

#load network from files

load("/media/sf_Shared_Folder/Data/Herschel/seq/fungi/fungi_communities_composition/Network/Network_topsoil_5reads_absolute_abun_r06_igraph.RData")
load("/media/sf_Shared_Folder/Data/Herschel/seq/fungi/fungi_communities_composition/Network/Network_cryoOM_5reads_absolute_abun_r06_igraph.RData")

load("/media/sf_Shared_Folder/Data/Herschel/seq/fungi/fungi_communities_composition/Network/Network_site1_5reads_absolute_abun_r06_igraph.RData")
load("/media/sf_Shared_Folder/Data/Herschel/seq/fungi/fungi_communities_composition/Network/Network_site2_5reads_absolute_abun_r06_igraph.RData")
load("/media/sf_Shared_Folder/Data/Herschel/seq/fungi/fungi_communities_composition/Network/Network_site3_5reads_absolute_abun_r06_igraph.RData")


network1 <- simplify(network1, remove.multiple = T, remove.loops = T)  #to get rid of doubled edges
plot(network1,vertex.label=NA,vertex.size=5,layout=layout_nicely) #simply visualize the network


edges1 <- ends(network1, E(network1)) #get edge table
edges1 <- data.frame(edges1)

edges1$cor <- paste(edges1$X1,edges1$X2)
correls_cryoOM1$cor <- paste(correls_cryoOM1$Var1,correls_cryoOM1$Var2)
correls_cryoOM1.new <-correls_cryoOM1[match(edges1$cor, correls_cryoOM1$cor),]

correls_cryoOM1.new$direction<-"black"  #assign postive or negative to edges
correls_cryoOM1.new$direction[grep("-",correls_cryoOM1.new$Freq)]<-"red"

vertex.data1 <- as_data_frame(network1,what = "vertices")  #get vertex table
vertex.data1 <- cbind(vertex.data1,taxa_table_fungi[match(rownames(vertex.data1),rownames(taxa_table_fungi)),])


#Add color for each varibles seperatly for fungi

vertex.data1$phylum_name=factor(vertex.data1$Phylum, levels = c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Entomophthoromycota", "Entorrhizomycota", "Glomeromycota", "GS19", "Mortierellomycota", "Mucoromycota", "Olpidiomycota", "Rozellomycota", "unidentified"))
phylum_col_name=c("lightsalmon","deepskyblue1","red","gray","lightpink4","lightslateblue","darkgreen","chartreuse","yellow","chocolate1","cornflowerblue","darkmagenta")[vertex.data1$phylum_name]
vertex.data1$phylum_col=paste(phylum_col_name)

vertex.data1$guild_name=factor(vertex.data1$Order, levels = c("Symbiotroph", "Saprotroph", "Pathotroph", "Saprotroph-Symbiotroph", "Pathotroph-Symbiotroph", "Pathotroph-Saprotroph", "Pathotroph-Saprotroph-Symbiotroph", "Pathogen-Saprotroph-Symbiotroph", "Unassigned"))#convert guild into the factor
guild_col_name=c("lightsalmon","deepskyblue1","red","darkmagenta","lightpink4","lightslateblue","darkgreen","chartreuse","yellow")[vertex.data1$guild_name] #assign factor to different levels
vertex.data1$guild_col=paste(guild_col_name) #create new color with defined color

vertex.data1$module_name=factor(vertex.data1$Module, levels = c("Module 1","Module 22", "Module 21","Module 5","Module 28"))
module_col_name=c("red","gold","deepskyblue1","magenta","darkgreen")[vertex.data1$module_name]
vertex.data1$module_col=paste(module_col_name)

#Add color for each varibles seperatly for bacteria

vertex.data1$phylum_name=factor(vertex.data1$Phylum, levels = c("Acidobacteria","Actinobacteria","Alphaproteobacteria","Armatimonadetes","Bacteroidetes","Chloroflexi","Deltaproteobacteria","Firmicutes","Gammaproteobacteria","Gemmatimonadetes","Kiritimatiellaeota","Nitrospirae","Patescibacteria","Planctomycetes","Verrucomicrobia"))
phylum_col_name=c("violet","chartreuse3","chocolate1","beige","lightpink","aquamarine3","cyan1","brown1","blue1","darkgoldenrod1","tomato1","bisque3","coral1","black","darkgreen","deepskyblue")[vertex.data1$phylum_name]
vertex.data1$phylum_col=paste(phylum_col_name)


#creat a new network with the new edge and vertex tables
network1 <- graph_from_data_frame(d=correls_cryoOM1.new,vertices = vertex.data1,directed = F)

V(network1)$color <- vertex.data1$phylum_col#assign colors
V(network1)$size <- vertex.data1$abundance  #assign sizes

E(network1)$weight<- 0.3 #increase the width of the edge

E(network1)$color <- "grey8"  #change the edge line type
E(network1)$lty <- E(network1)$direction
E(network1)$lty <-gsub("pos","solid",as.character(E(network1)$lty))
E(network1)$lty <-gsub("neg","dashed",as.character(E(network1)$lty))

#ploting
l=layout_with_fr(network1) #Saving the layout in  l allows us to get the exact same result multiple times
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1) #By default, the coordinates of the plots are rescaled to the [-1,1] interval for both x and y. You can change that with the parameter rescale=FALSE and rescale your plot manually by multiplying the coordinates by a scalar. You can use norm_coords to normalize the plot with the boundaries you want.
plot(network1, vertex.size=log(vertex.data1$abundance*100)*2, edge.color=correls_cryoOM1.new$direction,vertex.label=NA, rescale=F, layout=l*1.2, edge.curved=curves, edge.width = E(network1)$weight) #curves function is located in script folder

#legend in the graphs
phylum_col=c("lightsalmon","deepskyblue1","red","gray","lightpink4","lightslateblue","darkgreen","chartreuse","yellow","chocolate1","cornflowerblue","darkmagenta")
names(phylum_col)=c("Ascomycota", "Basidiomycota", "Chytridiomycota", "Entomophthoromycota", "Entorrhizomycota", "Glomeromycota", "GS19", "Mortierellomycota", "Mucoromycota", "Olpidiomycota", "Rozellomycota", "unidentified")

guild_col=c("lightsalmon","deepskyblue1","red","darkmagenta","lightpink4","lightslateblue","darkgreen","chartreuse","yellow")# create new object with same color as above
names(guild_col)=c("Symbiotroph", "Saprotroph", "Pathotroph", "Saprotroph-Symbiotroph", "Pathotroph-Symbiotroph", "Pathotroph-Saprotroph", "Pathotroph-Saprotroph-Symbiotroph", "Pathogen-Saprotroph-Symbiotroph", "Undefined")

module_col=c("red","gold","deepskyblue1","magenta","darkgreen")
names(module_col)=c("Module 3","Module 16","Module 9","Module 5","Module 6")

phylum_col=c("violet","chartreuse3","chocolate1","beige","lightpink","aquamarine3","cyan1","brown1","blue1","darkgoldenrod1","tomato1","bisque3","coral1","black","darkgreen","deepskyblue")
names(phylum_col)=c("Acidobacteria","Actinobacteria","Alphaproteobacteria","Armatimonadetes","Bacteroidetes","Chloroflexi","Deltaproteobacteria","Firmicutes","Gammaproteobacteria","Gemmatimonadetes","Kiritimatiellaeota","Nitrospirae","Patescibacteria","Planctomycetes","Verrucomicrobia")


legend("bottomright", legend=levels(as.factor(vertex.data1$phylum_name)), fill = phylum_col, pt.cex = 1, cex = 0.5 , horiz = FALSE)
legend("bottomright", legend=levels(as.factor(vertex.data1$guild_name))  , fill = guild_col  , pt.cex = 3, cex = 0.6, horiz = FALSE)
legend("bottomright", legend=names(module_col)  , fill = module_col  , pt.cex = 3, cex = 0.7, horiz = FALSE)

legend(x=-1.5, y=-1.1, c("Newspaper","Television", "Online News"), pch=21, col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

write.graph(graph = network1, file = "CryoOM_filter_0.01OTUs_at_genus.gml", format = "gml")

#caculate topological features
edge_number=ecount(network1)#number of edged
node_number=vcount(network1)#number of node
average_path= round(average.path.length(network1),digit=2)# average path lengh

average_degree=round(mean(degree(network1)),digit=2)
sd_degree=round(sd(degree(network1)),digit=2)

Average_degree=paste(average_degree, sd_degree, sep = " ± ")

nodes_closeness = as.data.frame(closeness(network1, mode="all", weights=NA,normalized = F))
average_closness = round(mean(log10(nodes_closeness[,1])),digit=2)
sd_closness=round(sd(log10(nodes_closeness[,1])),digit=2)

Average_closness=paste(average_closness, sd_closness, sep = " ± ")

nodes_betweenness = as.data.frame(betweenness(network1, directed=T,weights=NA,normalized = F))
average_betweenness = round(mean(nodes_betweenness[,1]),digit=2)
sd_betweenness=round(sd(nodes_betweenness[,1]),digit=2)
Average_betweenness=paste(average_betweenness, sd_betweenness, sep = " ± ")

e_density = round(edge_density(network1,loops = F),digit=4)
diameter = diameter(network1)

nodes_transitivity = as.data.frame(transitivity(network1,type="local",vids=V(network1), weights = NA,isolates = "zero")) #clustering coefficient
Clustering_coefficient = round(mean(nodes_transitivity[,1]),digit=3)
Abundance=sum(vertex.data1$abundance)

module1 <- cluster_edge_betweenness(network1)
modularity_number=round(modularity(module1),digit=2)
module_number="38 (1,22,21,5,28)"
horizon="cryoOM"

module_count=sort(summary(as.factor(module1$membership)), decreasing=T)[1:5]
assort <- assortativity_degree(network1, directed = TRUE)

df=data.frame(horizon,edge_number,node_number,average_path,Average_degree, Average_closness,Average_betweenness, e_density,diameter, Clustering_coefficient,modularity_number, module_number, Abundance)
write.csv(t(df), file="cryoOM_AA_network.csv")
write.csv(module_count, file="Module_count_topsoil_top_100_OTUs_at_genus.csv")

#Modularity
module1 <- cluster_edge_betweenness(network1) #Girvan-Newman algorithm
vertex.data1$module_name="Module"
vertex.data1$Module=paste(vertex.data1$module_name, module1$membership)
sort(summary(as.factor(module1$membership)), decreasing=T)[1:5]#find top modules
write.csv(vertex.data1, file="Vertex_cryoOM_r07.csv")

#Module bar chart
module_subset=subset(vertex.data1, Module=="Module 1"|Module=="Module 7"|Module=="Module 2"|Module=="Module 16"|Module=="Module 6")
module_subset$module_relative_abundance=module_subset$abundance*100/sum(module_subset$abundance)#convert in relative abundance

ggplot(module_subset, aes(x=Module,y=module_relative_abundance, fill=Order))+geom_bar(stat = "identity", position = "fill")+scale_fill_manual(values = guild_col, name="Functional guilds")+Axis_manupulation_cleanup+scale_x_discrete(limits=c("Module 3","Module 13","Module 17","Module 10","Module 11"))+scale_y_continuous(labels = scales::percent)+xlab("")

#Add topology to the vertex data frame
nodes.degree <- as.data.frame(degree(network1,mode = "all"))
nodes.closeness <- as.data.frame(closeness(network1, mode="all", weights=NA,normalized = F))
nodes.betweenness <- as.data.frame(betweenness(network1, directed=T,weights=NA,normalized = F))
nodes.transitivity <- as.data.frame(transitivity(network1,type="local",vids=V(network1), weights = NA,isolates = "zero"))

vertex_df=data.frame(vertex.data1,nodes.degree,nodes.closeness, nodes.betweenness, nodes.transitivity)

module_col=c("red","gold","deepskyblue1","magenta","darkgreen")


no.3=c("3") #create new color varibles for each modules # topsoil
no.4=c("15")
no.16=c("16")
no.22=c("8")
no.7=c("1")
no.6=c("6")
no.9=c("9")
no.14=c("14")

no.3=c("2") #create new color varibles for each modules #cryoOM
no.4=c("8")
no.16=c("3")
no.22=c("1")
no.7=c("6")
no.6=c("5")
no.9=c("10")
no.14=c("15")


number_df=c(20,20,20,20,20)
module_df=c("Module 1","Module 2","Module 3","Module 4","Module 5")
color_df=c("red","gold","deepskyblue1","magenta","darkgreen")
pie(number_df, labels = module_df, col=color_df)

V(network1)$color <- "gray"
V(network1)$color[vertex.data1$membership %in% no.3] <- "red" #assinghed each modules to new color
V(network1)$color[vertex.data1$membership %in% no.4] <- "gold"
V(network1)$color[vertex.data1$membership %in% no.16] <- "deepskyblue1"
V(network1)$color[vertex.data1$membership %in% no.22] <- "magenta"
V(network1)$color[vertex.data1$membership %in% no.7] <- "darkgreen"
V(network1)$color[vertex.data1$membership %in% no.6] <- "chartreuse"
V(network1)$color[vertex.data1$membership %in% no.9] <- "seagreen1"
V(network1)$color[vertex.data1$membership %in% no.14] <- "yellow"



# Letâ€™s take a look at all available layouts in igraph
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 

# Remove layouts that do not apply to our graph.

layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]



par(mfrow=c(3,3), mar=c(1,1,1,1))

for (layout in layouts) {
    print(layout)
    l <- do.call(layout, list(network1)) 
    plot(network1, edge.arrow.mode=0, layout=l, main=layout) }

#color different phylum with randon selected color

library(RColorBrewer) 
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]  
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))  #creat a color library which contains 74 colors

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
col_vector=getPalette(16)

count(taxa_table_fungi$Phylum) #22
x = data.frame(unique(taxa_table_fungi$Phylum),sample(col_vector, 16),stringsAsFactors =F)  #randomly select 22 colors for each phylum
pie(rep(1,8),labels = x$unique.taxa_table_fungi.Phylum.,col = x$sample.col_vector..8.,cex=1.5)

vertex.data1$tax.color <- vertex.data1$Phylum

#assign color to each vertex based on Phylum taxanomy
for (i in 1:16) {
  vertex.data1$tax.color <- gsub(x[i,1],x[i,2],vertex.data1$tax.color)
}


GP = filter_taxa(GlobalPatterns, function(x) sum(x > 3) > (0.2*length(x)), TRUE) # Remove taxa not seen more than 3 times in at least 20% of the samples.

GPr  = transform_sample_counts(fungi_raw, function(x) x / sum(x) )
GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-4, TRUE)



# Create table, number of features for each phyla
table(tax_table(fungi_filter)[, "Phylum"], exclude = NULL)

# How many genera would be present after filtering?
length(get_taxa_unique(fungi_filter, taxonomic.rank = "Genus"))


# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(fungi_raw),
               MARGIN = ifelse(taxa_are_rows(fungi_raw), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(fungi_raw),
                    tax_table(fungi_raw))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})


#key stone species are those who have highest betweenness and closeness centrality