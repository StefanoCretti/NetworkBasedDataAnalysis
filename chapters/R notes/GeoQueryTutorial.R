install.packages("BiocManager") #install the package
BiocManager::install("GEOquery") #install GEOquery
library("GEOquery")
# load data
gse<- getGEO("GSE28489") # Get a specific Geo entry, it accesses the GEO database and unpacks the data
View(gse) #the first element of the list is the data that I want
# gse<- getGEO("GSE28489", destdir = ".", getGPL = FALSE)
gse <- gse[[1]]
View(gse)
gse #expression set, a special data type designed to contain several informations
# inspect matrix of expression values

gse$title



ex <- exprs(gse) #matrix of values, specialized function

View(ex)
dim(ex)
colnames(ex)

#in ex is present the data matrix of the file, 
# un'altra funzione è phenoData, che estrae le phenotypic informations.


pheno <- phenoData(gse)
pheno 

ex2 <- log2(ex) #transform the data

# analyze value distributions
boxplot(ex)
boxplot(ex2) #log transformed version

# The original data was really asymmetric, while instead the second
# is quite symmetric. For simmetricity, it is intended the 
# length of the box up and down, which in simmetrical conditions
# are comparable. The median is quite on the same line,
# so the normalization is not needed

dim(ex2)
ex2 <- na.omit(as.matrix(ex2))
# removes all the rows with NA values.
dim(ex2) #no NA  values were found



# PCA
pca <- prcomp(t(ex2)) #transpose of the matrix,
# as the function requires one sample per row and one
# per column, datastructure with different pieces of informations

View(pca)
plot(pca)
summary(pca)
screeplot(pca)
# draw PCA plot
grpcol <- c(rep("blue",5), rep("red",5), rep("green",5), rep("yellow",5) , rep("pink",3), rep("purple",5), rep("orange",5) )
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=10, col=grpcol) #first and second component
text(pca$x[,1], pca$x[,2], rownames(pca$data), cex=0.75) 