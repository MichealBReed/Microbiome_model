library('igraph')
library('psych')
library('MCL')


##################################################
MATRIX<-function(table){
  columnName<-colnames(table)
  RowName<-row.names(table)
  row.names(table)<-NULL
  colnames(table)<-NULL
  mat = matrix(data = as.numeric(as.matrix(table)),nrow = length(RowName), ncol = length(columnName))
}


SpecNet<-function(table, threshold){
  
  correlation_calculation<-cor(t(table), method = 'spearman')
  
  correlation_matrix <- as.matrix(correlation_calculation)
  
  p_matrix <- corr.test(correlation_matrix, method = 'spearman', use = 'pairwise', adjust = 'none')$p
  
  binary_matrix <- ifelse(p_matrix < .001, 1, 0)
  
  sig_matrix <- correlation_matrix * binary_matrix
  
  SM1<-sig_matrix
  
  SM2<-sig_matrix
  
  SM1[SM1 < threshold] = 0
  
  SM2[SM2 > -threshold] = 0
  
  NEWMAT<-SM1 + SM2
  
  set.seed(123)
  
  network <- igraph::simplify(graph_from_adjacency_matrix(NEWMAT, mode = c('undirected'), weighted = TRUE), remove.loops = TRUE)
  
  network <- delete.vertices(network, degree(network) == 0)
}



radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

####################################################

manuscript<-read.csv('C:\\Users\\nemsr\\OneDrive\\Manuscript\\coral_microbiome.csv', sep = ',', header = TRUE)
name<-manuscript[,1]
name<-name[4:length(name)]
manuscript<-manuscript[,2:ncol(manuscript)]

Columns<-colnames(manuscript)

manuscript<-manuscript[4:nrow(manuscript),]

row.names(manuscript)<-name

manuscript<-as.matrix(manuscript)

OuterTable<-manuscript[,grep('Outer',colnames(manuscript))]
InnerTable<-manuscript[,grep('Inner', colnames(manuscript))]



manuscript<-matrix(data = as.numeric(manuscript), nrow = length(name), ncol = length(Columns))


row.names(manuscript)<-NULL
colnames(manuscript)<-NULL





OuterNumeric<-MATRIX(OuterTable)
InnerNumeric<-MATRIX(InnerTable)


OUTER<-OuterTable[-which(rowSums(OuterNumeric) == 0),]
INNER<-InnerTable[-which(rowSums(InnerNumeric) == 0),]




OUTER2<-MATRIX(OUTER)
colnames(OUTER2)<-colnames(OUTER)
rownames(OUTER2)<-rownames(OuterTable)[-which(rowSums(OuterNumeric) == 0)]






INNER2<-MATRIX(INNER)
colnames(INNER2)<-colnames(INNER)
rownames(INNER2)<-rownames(InnerTable)[-which(rowSums(InnerNumeric) == 0)]









#Combined matrix
CombinedNumeric<-(OuterNumeric + InnerNumeric) / 2
REMOVE<-which(rowSums(CombinedNumeric) == 0)
COMBINED<-CombinedNumeric[-REMOVE,]
rownames(COMBINED)<-rownames(InnerTable)[-REMOVE]


COMBINEDCOR<-cor(t(COMBINED), method = 'spearman')






#NETWORK

###############################################
##Network Analysis
#Inner Network
INNERNET<- SpecNet(INNER2,0.7)

V(INNERNET)$vertex_degree <- degree(INNERNET)

E(INNERNET)$color[E(INNERNET)$weight > 0.7] <- 'green'
E(INNERNET)$color[E(INNERNET)$weight < -0.7] <- 'red'



INNEREIGEN<-sort(eigen_centrality(INNERNET)$vector, decreasing = TRUE)




#Scale Network according to degree
inner_lay<- layout_in_circle(INNERNET, order = V(INNERNET))
set.seed(123)

V(INNERNET)$color <- ifelse(eigen_centrality(INNERNET)$vector > 0.75, 'blue','orange')

lab.locs <- radian.rescale(x=1:length(V(INNERNET)), direction=-1, start=0)

inner_label <- 1:length(V(INNERNET))

plot(INNERNET, main = 'Inner Coral', edge.curved = .75,
     layout = inner_lay, vertex.label.degree = lab.locs, vertex.label.dist = 2, vertex.label = inner_label)





InnerB2 <- sort(betweenness(INNERNET, directed = FALSE, weights = NA), decreasing = TRUE)


#Outer Network

OUTERNET <- SpecNet(OUTER2, 0.7)

V(OUTERNET)$vertex_degree <- degree(OUTERNET)

E(OUTERNET)$color[E(OUTERNET)$weight > 0.7] <- 'green'
E(OUTERNET)$color[E(OUTERNET)$weight < -0.7] <- 'red'



#Scale network according to degree
outer_lay <- layout_in_circle(OUTERNET, order = V(OUTERNET))

V(OUTERNET)$color <- ifelse(eigen_centrality(OUTERNET)$vector > .75, 'blue', 'orange')

outer.lab.locs <- radian.rescale(x=1:length(V(OUTERNET)), direction=-1, start=0)

outer_label <- 1:length(V(OUTERNET))


outer.lab.locs <- radian.rescale(x=1:length(V(OUTERNET)), direction=-1, start=0)

set.seed(123)
plot(OUTERNET, main = 'Outer Coral', edge.curved = .75,
     layout = outer_lay, vertex.label = outer_label, vertex.label.degree = outer.lab.locs,
     vertex.label.dist = 2)

OuterB2 <- sort(betweenness(OUTERNET, directed = FALSE, weights = NA), decreasing = TRUE)



OUTEREIGEN <- sort(eigen_centrality(OUTERNET)$vector, decreasing = TRUE)





