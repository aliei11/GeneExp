# GeneExp
Gene Expression data derived from different Infectous Diseases

#The aim will be to re-analyse the data and compare the results with those obtaind in the paper and/or complement the analysis with novel insights.

###Project description###

#At the beginning, we read the tables required for our project. 
#The first table is consisted of gene expression values of different #samples. 
#The second one is the metadata for samples we want to study. 
#Both tables require headers and each, have their different separators.
Conditions=read.table(file = 'conditions.csv', header = TRUE, row.names = 1, sep = ',' )
RawData=read.table(file = 'GSE147507_RawReadCounts_Human.tsv', header = TRUE, row.names = 1, sep = '\t' )

•1~
#Filtering data (extracting our desired samples):
#After reading tables, since we only want genes that are available in our metadata table, we use a logical test and filter our genes.
#So, we use dim() function to see the dimension of our filtered data to see if the number of samples corresponds to the number of samples in meta data table.
SelectedCol=colnames(RawData)%in%rownames(Conditions)
dim(RawData[,SelectedCol])

#It is now time for extracting our desired samples from the raw data by sub selecting.
#Then, we have our desired data; 
#We use another logical test to see if all the samples are in the correct order based on our metadata table.
Extracted_Data=RawData[,SelectedCol]                  
colnames(Extracted_Data)==rownames(Conditions)

•2~
#In order to find how many genes are actually expressed in every sample we use the function below; we apply summation to our logical vector to sum the expressed genes for every sample:
#Genes that have at least 10 reads in at least one sample:
ExpressedGene=Extracted_Data>=10
HowMany=apply(ExpressedGene,1,sum)

#Since we also want to see how many genes are expressed at least one time, we use another logical test 
#and then use dim() function to see its dimension;
#How many genes are expressed at least 1 time: 
FilteredGenes= Extracted_Data[HowMany>0,]
dim(FilteredGenes)

#Let us see how many genes are expressed in different values; How many genes are never expressed? # 5006 How many genes are always expressed that is expressed in all the sample: #693
table(HowMany)

•2.1~ 
#We use a barplot to see how many genes are expressed in Each Sample; However, first we do our summation on each sample and then use a vector of our samples names just for simplicity of our bar plot. Then we plot number of genes expressed in every sample; How many genes expressed in each different sample:
NumExpinSamples=apply(ExpressedGene,2,sum)
barplotnames=as.character(colnames(Extracted_Data))
barplot(NumExpinSamples, ylab = "Rate of expression", col = rainbow(44),
        names.arg = barplotnames, las=2, space = 1, main='Expression Level')

•2.2~ ????? 
#Now we want to use a scatterplot to compare the total number of reads with only the genes that are considered expressed in our samples as we will see there is a correlation between the total number of reads and the expressed genes meaning that the more number of read the higher the chance of a gene to be expressed:
Totalreads=apply(Extracted_Data,2,sum)
log_Totalreads= log2(Totalreads)
plot(log_Totalreads, NumExpinSamples, pch=20 ,cex=1.5, col=rainbow(44), 
     main='Expressed Samples & Total Reads Correlation',
     xlab='Total Reads', ylab = 'Number of Expression')    
legend('topleft', legend = row.names(Conditions),col = rainbow(44), pch =20, cex =0.5 )

•3~ 
#Experimental conditions: In this part we add a new column to our metadata table called “Experimental Conditions”. To do this, first we use paste function to create names consisted of conditions requested by the question. Then we add them as a new column to our metadata table. We also use table function to see how many replicates are assigned to each condition:
Cond=paste(Conditions$virus,Conditions$cell_type,Conditions$series,Conditions$Drug.Treatment,sep="_")
Merged_design=cbind(Conditions,Experimental_condition=Cond)     #newly data frame that specifies different condition
table(Merged_design$Experimental_condition)

•3.1~ 
#Now we put our new condition table and filtered expressed genes in a DGElist by using edgeR package. First we run edgeR.
library(edgeR)
y=DGEList(counts = FilteredGenes, group = Merged_design$Experimental_condition)

#By using the first command below, we can understand that our data is not normalized. So we normalize our data and estimate the negative binomial(NB) models by using edgeR library exclusive functions:
attributes(y)
#Preforming Normalization:
y=calcNormFactors(y)
#Then NB:
y=estimateDisp(y)

•3.2~ 
#To understand the variability of the data, we need our data in a multi-dimensional scaling analysis which is known as MDS plot. Since we have different number of factors for each columns (Experimental condition,Cell_series,Cell_type,Kit) so we assign a vector of colors based on the number of factors for each column and then we plot them so we get similar colors for every same factor. MDS~PCA:
Rainbow_Colors=rainbow(16)                               ###????Change colors!!!???#####
ColorsExp=rep(Rainbow_Colors, c(3,2,3,2,3,3,3,2,3,3,3,3,3,3,2,3))
plotMDS(y,col=ColorsExp,labels = Merged_design$Experimental_condition,main='Experimental Condition')
Rainbow_Colors_series=rainbow(7)
Colors_series=rep(Rainbow_Colors_series,c(6,4,9,4,6,6,9))
plotMDS(y,col=Colors_series,labels = Merged_design$series,main='Cell Series')
RainbowColor_cell=rainbow(4)
Colors_Cell=rep(RainbowColor_cell,c(28,6,4,6))
plotMDS(y,col=Colors_Cell,labels = Merged_design$cell_type,main='Cell Type')
Rainbow_color_kit=rainbow(3)
Colors_kit=rep(Rainbow_color_kit,c(9,8,27))
plotMDS(y,col=Colors_kit,labels = Merged_design$Kit,main='Used Kit')


•TRACK1:BATCH EFFECTS•
###Description###
#Performing Differential Expression (DE) analysis of different cell line models:
#we concatenate group of y$samples
#S4 ~> IAV infected:
DE_S4=exactTest(y,pair = c('IAV_infected_A549_S4_none',
                                    'Mock_A549_S4_none'))     
#S5 ~> SARS-CoV-2 infected:
DE_S5=exactTest(y,pair = c('SARS-CoV-2 infected_A549_S5_none', 
                                     'Mock_A549_S5_none'))  
#S8 ~> RSV infected:
DE_S8=exactTest(y,pair = c('RSV infected_A549_S8_none', 
                                  'Mock_A549_S8_none')) 

#FDR is the value that shows us the chances of having False Positive in our experiment. According to question, we consider FDR=0.01 and we do this by using toptags() function.
FDR_S4=topTags(DE_S4,n=nrow(DE_S4))
FDR_S5=topTags(DE_S5,n=nrow(DE_S5))
FDR_S8=topTags(DE_S8,n=nrow(DE_S8))

•T1.1~ 
#For every comparison we assign each gene to each four possible categories then we use boxplot: DE_UP: gene is UPregulated and FDR is < 0.01; DE_DOWN: gene is DOWNregulated and FDR is < 0.01; not_DEup: gene is UPregulated and FDR is > 0.01; notDE_down: gene is DOWNregulated and FDR is >0.01

#FOR: FDR_S4
DE_UP_S4=FDR_S4$table$logFC>0 & FDR_S4$table$FDR<=0.01
UP_S4=FDR_S4$table[DE_UP_S4,]

DE_DOWN_S4=FDR_S4$table$logFC<0 & FDR_S4$table$FDR<=0.01
DOWN_S4=FDR_S4$table[DE_DOWN_S4,]

NO_DEup_S4=FDR_S4$table$logFC>0 & FDR_S4$table$FDR>=0.01
not_DEup_S4=FDR_S4$table[NO_DEup_S4,]

NO_DEdown_S4=FDR_S4$table$logFC<0 & FDR_S4$table$FDR>=0.01
not_DEdown_S4=FDR_S4$table[NO_DEdown_S4,]


boxplot(UP_S4$logFC,DOWN_S4$logFC,not_DEup_S4$logFC,not_DEdown_S4$logFC,
        names=c('DE UP','DE DOWN','NO DE UP','NO DE DOWN'), las=2, space = 1,
        main='IAV infected', outline = FALSE,col = c('Darkblue','green','darkred','purple')) 

#FOR: FDR_S5
DE_UP_S5=FDR_S5$table$logFC>0 & FDR_S5$table$FDR<=0.01
UP_S5=FDR_S5$table[DE_UP_S5,]

DE_DOWN_S5=FDR_S5$table$logFC<0 & FDR_S5$table$FDR<=0.01
DOWN_S5=FDR_S5$table[DE_DOWN_S5,]

NO_DEup_S5=FDR_S5$table$logFC>0 & FDR_S5$table$FDR>=0.01
not_DEup_S5=FDR_S5$table[NO_DEup_S5,]

NO_DEdown_S5=FDR_S5$table$logFC<0 & FDR_S5$table$FDR>=0.01
not_DEdown_S5=FDR_S5$table[NO_DEdown_S5,]


boxplot(UP_S5$logFC,DOWN_S5$logFC,not_DEup_S5$logFC,not_DEdown_S5$logFC,
        names=c('DE UP','DE DOWN','NO DE UP','NO DE DOWN'), las=2, space = 1,
        main='SARS-CoV-2 infected', outline = FALSE,col = c('pink','cyan','cadetblue','bisque'))

#FOR: FDR_S8
DE_UP_S8=FDR_S8$table$logFC>0 & FDR_S8$table$FDR<=0.01
UP_S8=FDR_S8$table[DE_UP_S8,]

DE_DOWN_S8=FDR_S8$table$logFC<0 & FDR_S8$table$FDR<=0.01
DOWN_S8=FDR_S8$table[DE_DOWN_S8,]

NO_DEup_S8=FDR_S8$table$logFC>0 & FDR_S8$table$FDR>=0.01
not_DEup_S8=FDR_S8$table[NO_DEup_S8,]

NO_DEdown_S8=FDR_S8$table$logFC<0 & FDR_S8$table$FDR>=0.01
not_DEdown_S8=FDR_S8$table[NO_DEdown_S8,]


boxplot(UP_S8$logFC,DOWN_S8$logFC,not_DEup_S8$logFC,not_DEdown_S8$logFC,
        names=c('DE UP','DE DOWN','NO DE UP','NO DE DOWN'), las=2, space = 1,
        main='RSV infected', outline= FALSE,col = c('orange','darkolivegreen','deeppink4','deepskyblue3'))

•T1.2~ 
#Now we use a venn diagram to compare and intersection of differential expressed genes (DEGs) by using down regulated genes for each sample and also using up regulated genes for each sample.
library(VennDiagram)
#Venn Diagram for UP-regulated genes:
UP_List=list(rownames(UP_S4),rownames(UP_S5),rownames(UP_S8))
venn.diagram(UP_List, category.names = c('IAV','Covid19','RSV'),
             filename = 'UPreg.png',
             imagetype= 'png', main = 'UP Regulation between infected samples', 
             col= c('blue', 'red', 'green'))

#Venn Diagram for DOWN-regulated genes:
DOWN_List=list(rownames(DOWN_S4),rownames(DOWN_S5),rownames(DOWN_S8))
venn.diagram(DOWN_List, category.names = c('IAV','Covid19','RSV'),
             filename = 'DOWNreg.png',
             imagetype= 'png', main = 'DOWN Regulation between  infected samples', 
             col= c('blue', 'red', 'green'))
             
•T1.3~ 
#Now that we have samples that we applied our criteria on, it is time to compare our samples with the original samples from Blanco_melo et al paper; So, to do this we need their table.
Paper_S1=read.table(file = 'Table_S1_Blanco-Melo_et.al.csv', header = TRUE, row.names = 1, sep = ','  )

#Now we apply the same criteria we did to our samples on the original paper samples and compare the results with each other. We just consider deregulated genes with |logfc|>1 & p_adjusted_value<0.05.
blanco_S1=Paper_S1[Paper_S1$padj_SARS.CoV.2.A549. < 0.05 & (Paper_S1$SARS.CoV.2.A549._L2FC > 1 | Paper_S1$SARS.CoV.2.A549._L2FC < -1),]

•T1.3.1~ 
#Now we plot our data to see the differences.
DEGs=cbind(sum(nrow(blanco_S1)),sum(nrow(UP_S4)),sum(nrow(UP_S5)),sum(nrow(UP_S8)))
barplot(DEGs,col =c('palevioletred3'),
        names=c('Blanco','UPREG IAV','UPREG Covid-19','UPREG RSV'),
        main="Blanco melo vs MY DATA")
apply(DEGs,1,sum)

•1.3.2~ 
#Install and load the ggVennDiagram package
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")

#UP REG comparison
UP_List=list(rownames(blanco_S1),rownames(UP_S4),rownames(UP_S5),rownames(UP_S8))
ggVennDiagram(UP_List, label_alpha = 0,
              category.names = c('Blanco','IAV','Covid19','RSV')
)
             
#DOWN REG comparison
DOWN_List=list(rownames(blanco_S1),rownames(DOWN_S4),rownames(DOWN_S5),rownames(DOWN_S8))
ggVennDiagram(DOWN_List, label_alpha = 0,
              category.names = c('Blanco','IAV','Covid19','RSV')
)












