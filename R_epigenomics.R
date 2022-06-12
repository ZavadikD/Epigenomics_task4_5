datasigmoid = read.csv("sigmoid_colon_regulatoryElements.genes.distances.tsv", header = FALSE, sep = "\t")
#median(datasigmoid[,2])
#mean(datasigmoid[,2])
print(paste("the median distance between the gene and the regulatory element in sigmoid colon is:", median(datasigmoid[,2]) ,"bp;", "the mean distance is:", mean(datasigmoid[,2]), "bp"))

datastomach = read.csv("stomach_regulatoryElements.genes.distances.tsv", header = FALSE, sep = "\t")
#median(datastomach[,2])
#mean(datastomach[,2])
print(paste("the median distance between the gene and the regulatory element in stomach is:", median(datastomach[,2]) ,"bp;", "the mean distance is:", mean(datastomach[,2]), "bp"))
