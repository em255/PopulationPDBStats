library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(Biostrings)
library(ggpubr)
AAFold <- fread("Protein_popstats.txt.gz")
gene_folding <- AAFold %>% group_by(transcript) %>% sample_n(1) %>%summarise(Gene_pLDDT=Gene_pLDDT)
AAFold <- AAFold[AAFold$freqs>0.05 & AAFold$freqs<0.95]
table(AAFold$Effect)
countdf <- cbind.data.frame(pLDDT=rep(seq(5,100,5),3),Effect=c(rep("HIGH",20),rep("MODERATE",20),rep("LOW",20)),Count=NA)

for(x in seq(5,100,5)){
  for (y in c("LOW","MODERATE","HIGH")) {
    countdf[countdf$pLDDT==x & countdf$Effect==y,"Count"] <- sum(AAFold$mean_pLDDT > x-5 & AAFold$mean_pLDDT < x & AAFold$Effect==y)
  }
}

countdf <- countdf %>% group_by(Effect) %>% mutate(Count=Count/sum(Count))
countdf$Effect <- factor(countdf$Effect,levels = c("LOW","MODERATE","HIGH") )

ggplot(countdf,aes(x=pLDDT,y=Count))+
  geom_bar(aes(fill=Effect),stat="identity",position = "dodge" )+
  ylab("Proportion of SNPs in Class")+
  ggtitle("MAF>5%")

AAFold$Effect <- factor(AAFold$Effect,levels = c("LOW","MODERATE","HIGH") )
AAFold <- AAFold[!is.na(AAFold$Effect),]

ggplot(AAFold,aes(x=Effect,y=mean_pLDDT))+
  geom_boxplot(aes(fill=Effect))


ggplot(AAFold, aes(x = AA_ProportionPosition, fill = Effect)) +                       # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.4, bins = 50)+
  xlab("% Protein Length")



ggplot(AAFold,aes(x=rASA,y=mean_pLDDT))+
  geom_point(alpha=0.005)+
  ggtitle(paste0("pLDDT vs rASA r=",round(cor(AAFold$rASA,AAFold$mean_pLDDT),2)))


ggplot(AAFold, aes(x = MAF, fill = Effect)) +                       # Draw overlaying histogram
  geom_histogram(aes(y = after_stat(density)),position = "identity", alpha = 0.2, bins = 50)+
  xlab("MAF")+
  ylab("% of Sites")


p1_pLDDT <- ggplot(AAFold,aes(x = MAF,y=mean_pLDDT,group=Effect))+
  geom_smooth(aes(col=Effect))+
  ylab("mean_pLDDT")
p2_rASA <- ggplot(AAFold,aes(x = MAF,y=rASA,group=Effect))+
  geom_smooth(aes(col=Effect))+
  ylab("rASA")
ggarrange(p1_pLDDT,p2_rASA)



AAFold_Nonsyn <- AAFold %>% filter(Effect=="MODERATE") %>% filter(nchar(alleles)==1) %>% mutate(AA_sub_ClassChange=ifelse(Ref_FunctionalClass==Alt_FunctionalClass,"Conservative","Radical")) %>% filter(!is.na(AA_sub_ClassChange))
p1_pLDDT <- ggplot(AAFold_Nonsyn,aes(x = MAF,y=mean_pLDDT,group=AA_sub_ClassChange))+
  geom_smooth(aes(col=AA_sub_ClassChange))

p2_rASA <- ggplot(AAFold_Nonsyn,aes(x = MAF,y=rASA,group=AA_sub_ClassChange))+
  geom_smooth(aes(col=AA_sub_ClassChange))
ggarrange(p1_pLDDT,p2_rASA)


data(BLOSUM62)
substitute_abbreviations <- function(abbreviations) {
  # Define a mapping of single-letter abbreviations to three-letter amino acid abbreviations
  abbreviation_map <- c("A" = "Ala", "R" = "Arg", "N" = "Asn", "D" = "Asp", "C" = "Cys",
                        "Q" = "Gln", "E" = "Glu", "G" = "Gly", "H" = "His", "I" = "Ile",
                        "L" = "Leu", "K" = "Lys", "M" = "Met", "F" = "Phe", "P" = "Pro",
                        "S" = "Ser", "T" = "Thr", "W" = "Trp", "Y" = "Tyr", "V" = "Val")

  # Substitute the abbreviations
  if (is.vector(abbreviations)) {
    # Vector input
    output <- sapply(abbreviations, function(x) abbreviation_map[[x]])
  } else {
    # Single abbreviation input
    output <- abbreviation_map[[abbreviations]]
  }

  return(output)
}

#Get down to 20 standard AA
BLOSUM62 <- BLOSUM62[1:20,1:20]
colnames(BLOSUM62) <- toupper(substitute_abbreviations(colnames(BLOSUM62)))
rownames(BLOSUM62) <- toupper(substitute_abbreviations(rownames(BLOSUM62)))



AAFold_Nonsyn <- AAFold_Nonsyn %>% mutate(BLOSOM_Score=NA)
AAFold_Nonsyn$BLOSOM_Score <- apply(AAFold_Nonsyn,1, function(x) BLOSUM62[x[7],x[20]])
AAFold_Nonsyn <- AAFold_Nonsyn[AAFold_Nonsyn$BLOSOM_Score>-4,]
AAFold_Nonsyn <- AAFold_Nonsyn[AAFold_Nonsyn$BLOSOM_Score<4,]
AAFold_Nonsyn$BLOSOM_Score <- as.factor(AAFold_Nonsyn$BLOSOM_Score)




p1_pLDDT <- ggplot(AAFold_Nonsyn,aes(x =MAF,y=mean_pLDDT,group=BLOSOM_Score))+
  geom_smooth(aes(col=BLOSOM_Score))+
  ylab("mean_pLDDT")

p2_rASA <- ggplot(AAFold_Nonsyn,aes(x =MAF,y=rASA,group=BLOSOM_Score))+
  geom_smooth(aes(col=BLOSOM_Score))+
  ylab("rASA")
ggarrange(p1_pLDDT,p2_rASA)




AAFold <- fread("Protein_pLDDT.txt.gz")
colnames(AAFold)
AAFold <- AAFold %>% select(mean_pLDDT, rASA, AA_ProportionPosition, FourFoldDegenerate, AminoAcid_1let)


pdf(file="FoldingConfidence_byposition.pdf")
ggplot(AAFold,aes(x = AA_ProportionPosition))+
  geom_smooth(aes(y=rASA),col="blue")+
  geom_smooth(aes(y=mean_pLDDT/100),col="red")+
  scale_y_continuous(
    # Features of the first axis
    name = "rASA",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.,name="pLDDT")
  ) + 
  theme_ipsum() +
  theme(
    axis.title.y = element_text(color = "blue", size=13),
    axis.title.y.right = element_text(color = "red", size=13)
  ) +
  ggtitle("pLDDT&rASA by AA position")+
  xlab("AA position")
dev.off()


#Disorder Tail search
AAFold <- head(AAFold,100000)
AAFold <- AAFold %>% mutate(fivep_disord=0) %>% mutate(threep_disord=0)
Flag_DisordTail <- function(transc){
  fivepdisord <- 1
  fivepdisord <- 1
  
  transc_df <- AAFold %>% filter(transcript==transc)
  for(x in 1:max(transc_df$AA_Pos)){
    if(mean(transc_df$mean_pLDDT[transc_df$AA_Pos==x]) < 70){transc_df$fivep_disord[transc_df$AA_Pos==x] <- 1}
    else{break}
  }
  for(y in max(transc_df$AA_Pos):1){
    if(mean(transc_df$mean_pLDDT[transc_df$AA_Pos==y]) < 70){transc_df$threep_disord[transc_df$AA_Pos==y] <- 1}
    else{break}
  }
  transc_df <- transc_df %>% select(Chr,Pos,transcript,AminoAcid,fivep_disord,threep_disord)
  return(transc_df)
}

Disorder_DF <- bind_rows(lapply(unique(AAFold$transcript),function(T) Flag_DisordTail(T)))

write.table(Disorder_DF,"DisorderedTails.txt",sep="\t",col.names = T,row.names = F,quote=F)


