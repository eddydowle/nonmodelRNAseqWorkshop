#Eddy Dowle 2020

#using R 3.6.2

#for today we are going to be using the edgeR library to do our differential expression analysis
#you can view the library here:
#https://bioconductor.org/packages/release/bioc/html/edgeR.html

#this is meant to be an introduction to differential expression. Dont get caught up in the R code. If you dont understand the code just run it put a comment above it to note it and try to focus on the point of what we are doing rather than the code itself. We wont have time this afternoon to go through the in's and out's of the code itself, so please focus on the broader picture. You can come back later on in your own time and work through any code bits that trouble you (google is your best friend for this). 

#disclaimer: this is not a comprehensive differential expression analysis that you can push your data through. This is just an introduction. You are expected to read the manuals and go through the examples there yourself. I have cut some big corners in this code, simply to save time for this practical season - it is therefore not suitable to be used in a black box manner (and you should not do your statistical analysis in a black box manner ~ regardless of what you are doing). 

#to install bioconductor packages:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")

library(edgeR)
library(tidyverse)

#we are going to be looking at a time series of data across 5 months

#setwd to whereever you have put the data
setwd("~/Documents/AdminStuff/Otago/Teaching/RNAseq/data/")

##########################
#getting your data into R#
##########################
#draw on board final table aim

#we are going to start of looking at gene level counts so lets have a look at these counts. These are to the gene using RSEM which is a older program (you probably will likely use something newer like featurecounts or salmon) 

#explain data

#each samples file looks like:
#gene_id	transcript_id(s)	length	effective_length	expected_count	TPM	FPKM
#gene0	rna0	2361.00	2262.63	62.00	2.47	2.29
#gene1	rna1	102.00	3.63	0.00	0.00	0.00
#gene10	rna16	1832.00	1733.63	1.00	0.05	0.05
#gene100	rna192,rna193,rna194,rna195	1444.22	1345.85	2161.00	144.62	134.13
#....

#because we are doing differential expression analysis remember we are after the counts not TPM or FPKM

#we want to build a table that looks like gene_id, counts for sampleA, counts for smapleB ..... 

#there is a few ways to get your data into R
#im going to show you three ways as getting data in often seems to be an issue

#you can also just do this yourself by looping through the files:
file_list<-list.files(pattern='*genes.results')
file_list
#have a look at one of those files
tmp <- read.table(file_list[1], header =TRUE)
head(tmp)

#create empty table
table<-c()
library(stringr)
for (file in file_list){
  sampleID<-str_match(file,"(^[A-Za-z0-9]*_[A-Za-z0-9]*_[A-Za-z0-9]*_[0-9]*)")[,1] #grab sample name
  print(sampleID)
  df<-read.table(file,header=T,sep="\t",row.names=NULL,stringsAsFactors = F) #bring in file at a time
  df2<-df[,c(1,5)] #get columns we want
  colnames(df2)[colnames(df2)=="expected_count"] <- paste(sampleID,sep="") #paste sample name into column name
  if (length(table) ==0){table<-df2 
  } else {
    table<-merge(table,df2,by="gene_id") } #merge samples together
}  

head(table)

#another easy way tximport
#I mentioned this package earlier as a way to move from transcript differential expression to gene level transcript differential expression
#we are just going to use it as a basic importer here
#BiocManager::install("tximport")

library(tximport)
file_list<-list.files(pattern='*genes.results')
file_list
#lets have a look at one of those files
tmp <- read.table(file_list[1], header =TRUE)
head(tmp)
#getting out sample names out of the file names and assigning them to our dataset
sampleID<-str_match(file_list,"(^[A-Za-z0-9]*_[A-Za-z0-9]*_[A-Za-z0-9]*_[0-9]*)")[,1]
sampleID
names(file_list) <- sampleID
#importing using tximport here we specify the program we used to generate the counts (in this case rsem) but tximport can deal with salmon/kallisto and others
?tximport
txi<-tximport(file_list, type='rsem')
summary(txi) #creates a list
head(txi$counts)

#another option is to use trinity's abundance_estimates_to_matrix perl script
#its a perl script so we are just going to run it and it will output a couple of files to our working directory
x <- 'perl abundance_estimates_to_matrix.pl'
system(x)
x <- paste('perl abundance_estimates_to_matrix.pl --est_method rsem --gene_trans_map none',paste(file_list,collapse=" "))
system(x)
#this trinity script is going to output a matrix here its going to call it isoform ~ its not its gene level but its because we didnt give it a gene_trans_map file
#if you did a trinity assembly and wanted to collapse to 'trinity genes' just change --gene_trans_map none to --gene_trans_map Trinity.fasta.gene_trans_map and it will output one to the transcript level and one matrix to the 'trinity gene' level 
counts<-read.table('rsem.isoform.counts.matrix',header=T,row.names=NULL,sep='\t')
x<-str_match(colnames(counts),"(^[A-Za-z0-9]*_[A-Za-z0-9]*_[A-Za-z0-9]*_[0-9]*)")[,1]
x[1]<-'gene_id'
colnames(counts)<-x

#if we have a look at these we can see that they produce the same thing
head(counts[order(counts$gene_id),])
head(table)
head(txi$counts)

#of course if you do your counts with featurecounts from the Rsubread package your data will already be in R and you wont have to import it. 
#featurecounts counts to the exon then sums to the gene
#library(Rsubread)
#?featureCounts

#if you have a pairwise comparison (e.g. treatment control) there is a package called bigPint https://lindsayrutter.github.io/bigPint which has some pretty cool functions for checking your data and graphing later on. Its less useful for time series stuff but still worth looking at. 

######################################
#annotations direct to genome example#
######################################

#our annotation file is going to depend on whether we have a genome or transcriptome and how we annotated.
#today we are going to look at some fly data becaue it is a fly our annotations are going to be back to drosophila our closest model species

#lets have a look at it
pom_annotation<-read.table('fly_annotation.txt',header=T,row.names=NULL,sep='\t',stringsAsFactors = F)
#this file is generated a protein blast between our flies genes and drosophila
head(pom_annotation)
#if we go look at flybase.org what we actually want is flybase gene ids (note that I have already converted out species protein codes to gene codes here to save time)
#remember we always blast to the protein database so you always end up with protein rather gene codes
#there is a handy package we can use to convert between FBpp (proteins) to FBgn (gene names)
#you can use the library biomart to get the gene names
#but this is quite slow so we wont actually run the fetch part
library(biomaRt)
#set up for drosophila
ENSEMBL_Dm=useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
listAttributes(ENSEMBL_Dm)
#lots of different things to select from
#this is what we would run to transfer and you can add other attributes on as you like
#as this is slow we are not going to run this
#annotations_with_Dm <- getBM(attributes=c('flybase_transcript_id','external_gene_name', 'ensembl_gene_id', 'ensembl_peptide_id'),
#                             filters='flybase_translation_id',
#                             values = pom_annotation$Flybase_blast,
#                             mart = ENSEMBL_Dm)
#but that code produced this file:
annotations_with_Dm <- read.table('annotations_with_Dm',row.names=NULL,header=T,sep='\t')
head(annotations_with_Dm)

#the other way you can do this is download the file yourself 
#in some ways this is safer because the genome version for model species update and this will keep your files inline with the version of the genome you annotated too
#for instance to flybase you can download the conversion file here:
#ftp://ftp.flybase.net/releases/
#choose the genome you blasted too
#under the precomputedfiles/gene folder there is a file labeled fbgn_fbtr_fbpp_fb_xxx.tsv.gz 
#we are going to use this today as the genome version this was annotated to was from a few years ago
#were these files are will vary depending on the database you are using ~ so you may have to hunt around
#for flybase these files also report other drosophila species codes ~ hence why the row numbers are very high

flybase_gn<-read.table('fbgn_fbtr_fbpp_fb_2016_05.tsv',header=T,row.names=NULL,sep='\t',stringsAsFactors = F)
head(flybase_gn)
flybase_gn$FlyBase_FBtr<-NULL

#the other one we are going to bring in is the common name annotation file, which I have striped down slightly but it gives the common name
flybase_cn<-read.table('fbgn_annotation_ID.tsv',header=T,row.names=NULL,sep='\t',quote='',stringsAsFactors = F)
head(flybase_cn)

flybase<-left_join(flybase_gn,flybase_cn,by=c('FlyBase_FBgn' = 'primary_FBgn'))
head(flybase)
pom_annotation_full<-left_join(pom_annotation,flybase,by=c('Flybase_blast'='FlyBase_FBpp'))
head(pom_annotation_full)


###############################
#annotations trinotate example#
###############################

#trinotates format is a special form of torture
#this is just the first 1000 lines of a trinotate file for us to play with

#bring in trinotate file
trinotate_example<-read.table('trinotate_example.txt', sep="\t", header=TRUE, row.names=NULL, quote="", na.strings=".",stringsAsFactors = F) 
head(trinotate_example)
colnames(trinotate_example)
#in this set SPU was our custom database the rest of it are the standard trinotate columns

#I find that its easiest to simplfy down the annotations to make them human readable
#to grab the bits of annotation I am interested in I use regular expressions. You dont need to use regular expressions for most of these as you can just split the string on a character. However I have found that the flexibility of regular expressions make them very helpful for dealing with trinotate files ~but you should use whatever works best for you. 

#I recommend you use a regular expression tester if you choose to go down this route e.g. regex101.com

annotationsSimplified = trinotate_example %>% mutate(sprot_BLASTX_gene=gsub('_.*',"",sprot_Top_BLASTX_hit))
annotationsSimplified = annotationsSimplified %>% mutate(sprot_BLASTP_gene=gsub('_.*',"",sprot_Top_BLASTP_hit))

#this gsub will definitly change depending on the custom database used
annotationsSimplified = annotationsSimplified %>% mutate(SPU_Top_BLASTX=gsub('\\^.*',"",SPU_BLASTX))
annotationsSimplified = annotationsSimplified %>% mutate(SPU_Top_BLASTP=gsub('\\^.*',"",SPU_BLASTP))

#once you have pulled your custom gene out you can then do like we did before with flybase to turn this into a gene name / common name (here you can imagine the SPU_xxxx code is similar to the FBppXXXXX code in the fly example ~ you will need to find the map file)
#for instance this species (SPU) has a file that is downloadable that looks like:

#SPU_000622	Sp-Sap130	Sin3A associated protein
#SPU_001985	Sp-Hypp_1361	hypothetical protein-1361
#SPU_007826	Sp-Hypp_1748	hypothetical protein-1748
#SPU_016202	Sp-Hypp_2252	hypothetical protein-2252
#SPU_022728	Sp-TpaseL_14	transposase-like-14
#SPU_003767	Sp-Hypp_397	hypothetical protein-397
#SPU_018712	Sp-Ipo7-2	importin 7-2
#SPU_000628	Sp-Ankrd53	ankyrin repeat domain 53

#which we could map in to give a semi-sane gene name to the SPU name

#you can then also pull the GO terms out - having them in a column by themself will allow you to plug them straight into a GO enrichment analysis
annotationsSimplified = annotationsSimplified %>% mutate(GO_code=gsub("\\^.*","",gene_ontology_blast))
annotationsSimplified = annotationsSimplified %>%mutate(GO_function=gsub(".*\\^","",gene_ontology_blast))
#you can do the same thing for the kegg terms etc

#there is one more quirk with trinotate files that I have found. Usually I have found that there is a small subset of trinity genes that have more than one protein prediction and therefore blast hit. This is really annoying as for edgeR you need to have one line per gene for the annotation. If you want to keep this information you will need to collapse them into a string. However by collapsing into a string you will lose the ablity to blindly copy and paste into enrichment analysis (although there is a solid chance none will be significant so it may never be an issue)
#you can see this issue here where of the same transcript you have two protein predictions
annotationsSimplified %>% dplyr::filter(transcript_id=='TRINITY_DN21438_c0_g1_i8')

#a work around for this issue
colnames(annotationsSimplified)
annotationsSimplified_cutdown<-annotationsSimplified[c(1,2,19:24)]
head(annotationsSimplified_cutdown)

annotationsCleanSingle <- annotationsSimplified_cutdown %>% group_by(gene_id,transcript_id) %>%
  summarize(sprot_BLASTX_gene = paste(unique(sprot_BLASTX_gene[!is.na(sprot_BLASTX_gene)]), collapse = ","),
            sprot_BLASTP_gene = paste(unique(sprot_BLASTP_gene[!is.na(sprot_BLASTP_gene)]), collapse = ","),
            SPU_Top_BLASTX = paste(unique(SPU_Top_BLASTX[!is.na(SPU_Top_BLASTX)]), collapse = ","),
            SPU_Top_BLASTP = paste(unique(SPU_Top_BLASTP[!is.na(SPU_Top_BLASTP)]), collapse = ","),
            GO_code = paste(unique(GO_code[!is.na(GO_code)]), collapse = ","),
            GO_function = paste(unique(GO_function[!is.na(GO_function)]), collapse = ","))
#and you also need to do the same thing to collapse to a trinity gene prediction rather than a transcript preduction (hopefully all the transcripts blasted to the same thing but some probably didnt)

##################################
#differential expression in edgeR#
##################################
#edgeR has a very detailed user guide found here:
#https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#it has many different examples of different analyses etc, if you use edgeR read the manual and go through the examples to find one similar to your project
#DESEQ2 also has a very through user manual too.

#we have our counts and our annotations now we can progress to a differential expression analysis

table_annotation<-right_join(pom_annotation_full,table, by="gene_id",all=TRUE)
head(table_annotation)
colnames(table_annotation)

#in this dataset we have five samples per time point and five time points
#we want to create our groupings for the samples - here I am just going to label them by month
grouping<-c(rep('M2',5),rep('M3',5),rep('M4',5),rep('M5',5),rep('M6',5))
grouping

#now we need to create an edgeR object that includes our data and annotations
colnames(table_annotation)
?DGEList
#we are going to use edgeR to normalise so the options we are interested in are counts, genes and group
#this command creates a dataset called a DGEList which is what edgeR uses for all its analyses - we will save things like the dispersion estimates etc back into this object as we go through the script
y<-DGEList(counts = table_annotation[,5:29], genes = table_annotation[,1:4], group = grouping)
summary(y)
#it contains our raw counts:
head(y$counts)
#our sample infomation with groupings (sanity check your groupings are correct), because we havent normalised yet the norm.factor is 1
y$samples
#and any annotations we added in
head(y$genes)

#filtering
#we want to remove genes that dont have sufficiently large counts for a statistical analysis
#not only are these low genes not biologicaly relavent to us (you must have some gene expression to be of interest), but they also interfer with the statistics later on
#there quite a few different ways to do this
#one way is to use the inbuilt function in edgeR. This function uses the minimum group size in your samples (in this case 5) and then finds which rows have counts in at least that minumum number regardless of group. 
?filterByExpr
keep <- filterByExpr(y)
#pre-filtering we have
nrow(y$counts)
y <- y[keep, , keep.lib.sizes=FALSE]
#post-filtering we have
nrow(y$counts) #13414

#another way to filter that doesnt involve the group size factor is to just use cpm 
y2<-DGEList(counts = table_annotation[,5:29], genes = table_annotation[,1:4], group = grouping)
keep<-rowSums(cpm(y2)>1) >=2
y2 <- y2[keep, , keep.lib.sizes=FALSE]
nrow(y2$counts) #13175

#normalising
#to normalise by library size we use edgeR's calcNormFactor - here we are just going to stick with the default TMM method
#remember we are normalising to remove sample-specific effects such as library size
?calcNormFactors
y$samples
y<-calcNormFactors(y)
y$samples
#now when we look at the samples we can see that normalisation factor is no longer just 1 it varies based on library size and composition

#before we go into our model design and analysis we can look at our raw data for trends using a MDS (deseq2 does a PCA) plot
#in this plot the distances between samples correspond to the leading log-fold-changes between sample pairs (unsupervised clustering) 
?plotMDS
plotMDS(y)
#make it pretty
#remember this is data from samples that have no morphological change through time so we are not expecting to see amazing changes like I showed you in the presentation with the sex-changing fish

#you can make that looks slightly prettier in plotMDS
cols=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00')
shps = c(15, 16, 17,18,19)
plotMDS(y, col = cols[as.factor(grouping)], pch = shps[as.factor(grouping)])
legend('topleft', col=cols, legend=levels(as.factor(grouping)), pch = shps, cex = 0.7)

#or you can just save our the co-ordinates and plot in your preferred plotting software
mds_plot<-plotMDS(y)
plot_mds<-as.data.frame(mds_plot$cmdscale.out)
plot_mds
ggplot(plot_mds, aes(x=V1,y=V2, colour=grouping,shape=grouping)) +
  geom_point( size=2)+
  theme_bw() +
  xlab('Leading logFC dim1') + 
  ylab('Leading logFC dim2')
  
#now building our model and calculating the all important dispersion 
#firstly lets build our model
grouping
#we could use our earliest time point, 2M, as our intercept 
design<-model.matrix(~grouping)
design
#or we could have set that to zero by going:
design<-model.matrix(~0+grouping)
design #we will use this one for today but setting the intercept is really handy if you have a control/contrast and doesnt limit your contrasts in any way
#just like for any other model design if you have a second factor you can add that in, lets just pretend we had a males and females as well
sex<-c(rep('F',15),rep('M',10))
model.matrix(~0+grouping*sex)

#just renaming to make our model looks sensible
colnames(design) <- levels(y$samples$group)
rownames(design) <- colnames(y)
#always sanity check that your design is correct (i.e. the 1's and 0's make sense)
design

#now we can estimate dispersion
#I mentioned earlier edgeR has three calculations for this, the function estimateDisp will do all three in one go
?estimateDisp
y <- estimateDisp(y, design)
y$common.dispersion
#the squareroot of the common dispersion gives the coefficient of biological variation
sqrt(0.066) #0.26
#which we can see if we plot the dispersion out using the plotBCV function in edgeR
plotBCV(y)
#you can see we have higher dispersion at our lower counts and then it dwindles down at the higher counts ~this is pretty normal

################
#glm's in edgeR#
################

#there is a couple of glms in edgeR
#the traditional models is:
?glmFit
#the newer model is:
?glmQLFit
#the Ql model uses a quasi-liklihood test and is more conservative (better type I error control, or less false positives)
#however from my experience the QL model is sensitive to replicate numbers and can be quite ruthless on datasets
#by default they also use different dispersion estimates glmFit uses genewise and glmQL uses trended

#to fit the bionomical generalized linear models edgeR needs the DGElist an the model design
fit <- glmQLFit(y, design)
#or
fit <- glmFit(y, design) #im just going to stick with this one for this example ~ as I want as many DE genes to come through so we can look at them in various ways, keep in mind the QL model is the recomended model. 
#this creates a large DGEGLM object
summary(fit)
#which we then need to use to test for our various comparisons
?glmLRT
?glmQLFTest #glmQLFit uses this for the comparisons

#to test for comprisons we need to tell edgeR what contrasts we want to make. There is three ways to do this. 
#test the null hypothesis that time point x - time point y = 0
#make contrasts method one
#use makecontrasts function in edgeR
?makeContrasts

#contrasts need to be based on models column names
#draw time two time series on board
my.contrasts <- makeContrasts(
  M3vsM2 = M3-M2,
  M4vsM2 = M4-M2,
  M5vsM2 = M5-M2,
  M6vsM2 = M6-M2,
  M4vsM3 = M4-M3,
  M5vsM4 = M5-M4,
  M6vsM5 = M6-M5,
  levels = design)

#we can then use these contrasts to look for differentially expressed genes in our analyses
M3vsM2_contrast<-glmLRT(fit, contrast=my.contrasts[,"M3vsM2"])
#this creates and LRT object
M3vsM2_contrast
#to get the bits we want out we can use edgeRs topTags function, the FDR method by default is Benjamini-Hochberg (BH) 
?topTags
topTags(M3vsM2_contrast)
#by default it will give you to the top 10 to get all out:
M3vsM2_contrast_toptags <- data.frame(topTags(M3vsM2_contrast,n=nrow(y$counts),sort="none"))
head(M3vsM2_contrast_toptags)
M3vsM2_contrast_toptags %>% filter(FDR < 0.05) %>% nrow()

#contrasts method two:
#now lets see what that looks like:
my.contrasts
#lets take M3vsM2 as our example, the resulting contrast is -1,1,0,0,0
#if we go back to our design
head(design)
#if we delete 1,0,0,0,0 from 0,1,0,0,0  we get -1,1,0,0,0 as I showed you earlier on the board with the more complicated model
#we can call the exact same contrast by plugging in those numbers
M3vsM2_contrast2 <- glmLRT(fit,contrast=c(-1,1,0,0,0))
M3vsM2_contrast2_toptags <- data.frame(topTags(M3vsM2_contrast2,n=nrow(y$counts),sort="none"))
head(M3vsM2_contrast2_toptags)
M3vsM2_contrast2_toptags %>% filter(FDR < 0.05) %>% nrow()

#above we are comparing the null hypothesis that M3-M2 is equal to zero 
#in the above things that are a positive fold change represent genes that are upregulated in 3M compared to 2M 
#note that this is not the same test: 1,-1,0,0,0 comparing whether M2-M3 is equal to zero

#contrasts method 3:
#the other way to call a model requires a design that is not zeros, lets have a quick look at it
design_2Mintercept<-model.matrix(~grouping)
colnames(design_2Mintercept) <- levels(y$samples$group)
rownames(design_2Mintercept) <- colnames(y)
design_2Mintercept
y2 <- estimateDisp(y, design_2Mintercept)
fit2 <- glmFit(y2, design_2Mintercept)

M3vsM2_contrast3 <- glmLRT(fit2,coef=2)
M3vsM2_contrast3_toptags <- data.frame(topTags(M3vsM2_contrast3,n=nrow(y$counts),sort="none"))
head(M3vsM2_contrast3_toptags)
M3vsM2_contrast3_toptags %>% filter(FDR < 0.05) %>% nrow()

#so we have three ways of calling the same contrast to test the null hypothesis that M3-M2 is equal to zero
#using make.contrasts M3vsM2
topTags(M3vsM2_contrast)
#using contrast=c(-1,1,0,0,0)
topTags(M3vsM2_contrast2)
#or if we have a 0 intercept model we can use coef=2
topTags(M3vsM2_contrast3)

#use which ever is easist for you
#note that with a set intercept you can still use the +/- formula to compare say 5M to 6M in this case it would be 1,0,0,0,1 - 1,0,0,1,0 = 0,0,0,-1,1
design_2Mintercept
M6vsM5_contrast2 <- glmLRT(fit2,contrast=c(0,0,0,-1,1))
M6vsM5_contrast2_toptags <- data.frame(topTags(M6vsM5_contrast2,n=nrow(y$counts),sort="none"))
head(M6vsM5_contrast2_toptags)
M6vsM5_contrast2_toptags %>% filter(FDR < 0.05) %>% nrow()

#which is the same as (note im using our orginal 'fit' here which has a 0 intercept)
M6vsM5_contrast<-glmLRT(fit, contrast=my.contrasts[,"M6vsM5"])
M6vsM5_contrast_toptags <- data.frame(topTags(M6vsM5_contrast,n=nrow(y$counts),sort="none"))
head(M6vsM5_contrast_toptags)
M6vsM5_contrast_toptags %>% filter(FDR < 0.05) %>% nrow()

#its just a glm so its endlessless flexible. Im not going to go into more complicated constrasts here but you can do basically whatever you want so if you want to compare the average of months 3-6 to 2 you can do that ~ you just end up with fractions in the contrasts rather than ones and zeros. If you want to compare anything against the control that is easy too 'coef=2:6'. There is a lot of examples in the edgeR manual-as always read the manual throughly.

##################
#building a table#
##################

#what we really need to generate at this point is a table. This table is going to have our annotation, logfc and FDR values for our comparisons of interest
#we are going to do all the contrasts in the table and stick the bits we want into a table
my.contrasts

#M3vsM2
M3vsM2<-glmLRT(fit, contrast=my.contrasts[,"M3vsM2"])
M3vsM2.toptags <- data.frame(topTags(M3vsM2,n=nrow(y$counts),sort="none"))
head(M3vsM2.toptags)
contrast_table<-M3vsM2.toptags[,c(1:5,9)]
colnames(contrast_table)[colnames(contrast_table)=="logFC"] <- "M3vsM2_logFC"
colnames(contrast_table)[colnames(contrast_table)=="FDR"] <- "M3vsM2_FDR"
head(contrast_table)
contrast_table %>% filter(M3vsM2_FDR < 0.05) %>% nrow()


#M4vsM2
M4vsM2<-glmLRT(fit, contrast=my.contrasts[,"M4vsM2"])
M4vsM2.toptags <- data.frame(topTags(M4vsM2,n=nrow(y$counts),sort="none"))
contrast_table$M4vsM2_logFC<-M4vsM2.toptags[,5]
contrast_table$M4vsM2_FDR<-M4vsM2.toptags[,9]
head(contrast_table)
contrast_table %>% filter(M4vsM2_FDR < 0.05) %>% nrow()

#M5vsM2
M5vsM2<-glmLRT(fit, contrast=my.contrasts[,"M5vsM2"])
M5vsM2.toptags <- data.frame(topTags(M5vsM2,n=nrow(y$counts),sort="none"))
contrast_table$M5vsM2_logFC<-M5vsM2.toptags[,5]
contrast_table$M5vsM2_FDR<-M5vsM2.toptags[,9]
contrast_table %>% filter(M5vsM2_FDR < 0.05) %>% nrow()

#M6vsM2
M6vsM2<-glmLRT(fit, contrast=my.contrasts[,"M6vsM2"])
M6vsM2.toptags <- data.frame(topTags(M6vsM2,n=nrow(y$counts),sort="none"))
contrast_table$M6vsM2_logFC<-M6vsM2.toptags[,5]
contrast_table$M6vsM2_FDR<-M6vsM2.toptags[,9]
contrast_table %>% filter(M6vsM2_FDR < 0.05) %>% nrow()

#M4vsM3
M4vsM3<-glmLRT(fit, contrast=my.contrasts[,"M4vsM3"])
M4vsM3.toptags <- data.frame(topTags(M4vsM3,n=nrow(y$counts),sort="none"))
contrast_table$M4vsM3_logFC<-M4vsM3.toptags[,5]
contrast_table$M4vsM3_FDR<-M4vsM3.toptags[,9]
contrast_table %>% filter(M4vsM3_FDR < 0.05) %>% nrow()

#M5vsM4
M5vsM4<-glmLRT(fit, contrast=my.contrasts[,"M5vsM4"])
M5vsM4.toptags <- data.frame(topTags(M5vsM4,n=nrow(y$counts),sort="none"))
contrast_table$M5vsM4_logFC<-M5vsM4.toptags[,5]
contrast_table$M5vsM4_FDR<-M5vsM4.toptags[,9]
contrast_table %>% filter(M5vsM4_FDR < 0.05) %>% nrow()

#M6vsM5
M6vsM5<-glmLRT(fit, contrast=my.contrasts[,"M6vsM5"])
M6vsM5.toptags <- data.frame(topTags(M6vsM5,n=nrow(y$counts),sort="none"))
contrast_table$M6vsM5_logFC<-M6vsM5.toptags[,5]
contrast_table$M6vsM5_FDR<-M6vsM5.toptags[,9]
contrast_table %>% filter(M6vsM5_FDR < 0.05) %>% nrow()

#hint R functions are good for repetative tasks

#our final table:
head(contrast_table)

#how you analyse this data is ENTIRELY dependent on your question of interest
#I have provided some examples of things you can do below. But this is NOT a comprehensive guide on all the things you can do with this data and definitly not what you should do with yours. 
#think about your question and how best you can answer your question - only you can decide what that approach is. 
#think about how you are going to display your data to your audience in an effective manner. Read quality papers on similar questions to yours

#we need to do a couple of things and then we are good to go
#one of the first things is often with time series data we are interested in any gene that is DE across the time series
contrast_table$Lowest_FDR_from2M.TimeSeries<-apply(contrast_table[,c(6,8,10,12)],1,min)
contrast_table$Lowest_FDR_serial.TimeSeries<-apply(contrast_table[,c(6,14,16,18)],1,min)

head(contrast_table)
write.table(contrast_table,'Final_table_contrasts.txt',row.names=F,quote=F,sep='\t')

#we have 4923 significant genes in our time series from 2M
contrast_table %>% filter(Lowest_FDR_from2M.TimeSeries < 0.05) %>% nrow()
#we have 2028 significant genes in our serial time series
contrast_table %>% filter(Lowest_FDR_serial.TimeSeries < 0.05) %>% nrow()

#lets have a look at those genes in general 
df <- contrast_table %>% filter(Lowest_FDR_from2M.TimeSeries < 0.05) %>% dplyr::select(gene_id,M3vsM2_logFC,M4vsM2_logFC,M5vsM2_logFC,M6vsM2_logFC) %>% gather(., variable, value, -gene_id)
head(df)
ggplot(df, aes(variable, value,group = gene_id)) +
  geom_line()

df <- contrast_table %>% filter(Lowest_FDR_serial.TimeSeries < 0.05) %>% dplyr::select(gene_id,M3vsM2_logFC,M4vsM3_logFC,M5vsM4_logFC,M6vsM5_logFC) %>% gather(., variable, value, -gene_id)
ggplot(df, aes(variable, value,group = gene_id)) +
  geom_line()

#maybe we have a time point we are particularily interested in
df <- contrast_table %>% filter(M5vsM2_FDR < 0.05) %>% dplyr::select(gene_id,M3vsM2_logFC,M4vsM2_logFC,M5vsM2_logFC,M6vsM2_logFC) %>% gather(., variable, value, -gene_id)
ggplot(df, aes(variable, value,group = gene_id)) +
  geom_line()

#in general you can see because we have so much significant gene expression going on this is not a very useful way to look at this data
#what would be cool would be to cluster the dataset into genes with similar trajectories through time so we can see in general what do the major gene trajectories look like

#######################
#heatmaps & clustering#
#######################

#heatmaps are often a good way to visualise gene expresion changes through a time series

library(RColorBrewer)
library(gplots)
mat<-data.matrix(contrast_table %>% filter(Lowest_FDR_from2M.TimeSeries < 0.05) %>% dplyr::select(M3vsM2_logFC,M4vsM2_logFC,M5vsM2_logFC,M6vsM2_logFC))
heatmap.2(mat,col=brewer.pal(11,"BrBG"),trace="none",Colv=F,labRow=F,cexCol=0.5)

#heatmaps can often look quite washed out when there is a large spread of fold changes one option is to chop them:
mat[mat < -2] <- -2
mat[mat > 2]  <-  2
heatmap.2(mat,col=brewer.pal(11,"BrBG"),trace="none",Colv=F,labRow=F,margins=c(8,8),cexCol=1)

#this easily allows us to see that half our genes are generally upregulated through the time series and half are down regulated. It also looks like there is more going on in the later months (stronger colours) than the earlier months (lighter colours). But to really see this we need to look at the average trajectories of the genes.

#we cluster genes into trajectories as we kind of assume a guilty by association correlation. That if they are doing similar things at similar times they are linked to one another in some manner

#for this we are going to have a look at the R package WGCNA. This is not really what this package is designed to do. WGCNA is a designed to be run on the individual samples rather than the final fold changes ~which we are going to do here. Because of this Im not going to run everything you would normally run in WGCNA and you should obviously be reading a the tutorials so that for your dataset you run things in the most appropriate manner. 

#there is now more appropriate tools for what we are about to do and if you want to do something similar I suggest you dont do it this way and use one of the more appropriate tools.
#I was actually going to rewrite this for DP_DG_cluster https://github.com/PrincetonUniversity/DP_GP_cluster but I ran out of time so I've stolen some old code from back in the day before we had so many cool tools. Bascially here is an idea of an analysis you could do - but dont do it this way. 
#https://lindsayrutter.github.io/bigPint/

#in saying that the WGCNA package has a great set of tutorials which if you use it you should go through. But its not really designed for what we going to use it for. Its more about finding interesting clusters or modules within rather than trajectories across
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/

#######
#WGCNA#
#######

#there is an alternative to WGCNA called https://cemitool.sysbio.tools/analysis which has an online toolset as well as a R package 

library(WGCNA)
head(contrast_table)

#WGCNA expects a transformed numeric matrix
mat<-t(contrast_table %>% filter(Lowest_FDR_from2M.TimeSeries < 0.05) %>% select(gene_id,M3vsM2_logFC,M4vsM2_logFC,M5vsM2_logFC,M6vsM2_logFC))
colnames(mat) <- mat[1,]
mat <- mat[-1, ] 
mat<-as.matrix(mat)
class(mat) <- "numeric" 

#WGCNA weighted correlation network analysis
#here we are going to cluster genes into those with similar logfc trajectories we are going to run WGCNA without the 'W', i.e., no soft thresholding (this is very non-standard and is because I'm not running it 'properly' ~ normally you give it normalised expression values) 

#first we need to build the adjacency matrix or correlation network of logfc between genes (if we were running this properly there is steps before hand to select the softpower threshold and control the noise)
?adjacency
softPower = 1;
adjacency = adjacency(mat, power = softPower);
# Turn adjacency into topological overlap matrix (measure of similarily) and calculate dissimilarity matrix (subtracts it from 1)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function to cluster our genes
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

#set the miimum number of genes to cluster into a module
minModuleSize = 30;
# Module identification using dynamic tree cut across our precious dendrogram
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

#label those clusters with a colour and visualise it beneath the dendrogram
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(mat, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(mat, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;
mergedColors

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

#seperate modules by those going up and down to plot both sides
modNums<-length(unique(mergedColors))
upGenes<-matrix(nrow=ncol(t(mat)),ncol=modNums)
downGenes<-matrix(nrow=ncol(t(mat)),ncol=modNums)

cn=1
for (j in unique(mergedColors)) {
  for (i in 1:ncol(t(mat))) {
    submat<-t(mat)[mergedColors==j,i]
    upGenes[i,cn]<-mean(submat[submat > 0])
    downGenes[i,cn]<-mean(submat[submat < 0])
  }
  cn<-cn+1
}

colnames(upGenes)<-unique(mergedColors)
colnames(downGenes)<-unique(mergedColors)
upGenes
downGenes

#plot modules out - we are going to set 2M to zero for visualisation purpose
plotTraj<-function(x,y,main) {
  range<-c(min(y,na.rm=T),max(x,na.rm=T))
  plot(2:6,c(0,x[1:4]),ylim=c(-1.1,1.2),col='black',type='l',main=main,las=1) #apple, alt version with set y limits
  lines(2:6,c(0,(y[1:4])),col='black' )
  lines(2:6,rep(0,5),lty=3)
}

upGenes
downGenes

#plot modules out
for (i in 1:ncol(upGenes)) {
    plotTraj(upGenes[,i],downGenes[,i],paste(colnames(upGenes)[i],"_n=",sum(mergedColors==colnames(upGenes)[i]),sep=""))
}

class(mergedColors)
out<-as.data.frame(t(mat))
head(out)
out$module<-mergedColors
head(out)
write.table(out,"WGCNA_modules.txt",sep="\t",row.names=F,quote=F)

BiocManager::install("gage")
library(gage)
gage

#exon analysis
