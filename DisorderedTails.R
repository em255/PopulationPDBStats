library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)


AAFold <- fread("Protein_Structure.txt")
Orientation <- read.table("Transcript_Orientation.txt")
colnames(Orientation) <- c("transcript","orientation")
AAFold <- left_join(AAFold,Orientation)
AAFold <- AAFold %>% mutate(fivep_disord=0) %>% mutate(threep_disord=0) %>% filter(AminoAcid!="STOP")
AAFold$orientation[is.na(AAFold$orientation)] <- "+"

Flag_DisordTail <- function(transc){
  transc_df <- AAFold %>% filter(transcript==transc)
  for(x in 1:max(transc_df$AA_Pos)){
    if(nrow(transc_df[transc_df$AA_Pos==x,])>=1){
    if(mean(transc_df$mean_pLDDT[transc_df$AA_Pos==x],na.rm=T) < 70){transc_df$fivep_disord[transc_df$AA_Pos==x] <- 1}
    else if(mean(transc_df$rASA[transc_df$AA_Pos==x],na.rm=T) > 0.5){transc_df$fivep_disord[transc_df$AA_Pos==x] <- 1}
    else{break}
    }
  }
  for(y in max(transc_df$AA_Pos):1){
    if(nrow(transc_df[transc_df$AA_Pos==y,])>=1){
    if(mean(transc_df$mean_pLDDT[transc_df$AA_Pos==y],na.rm=T) < 70){transc_df$threep_disord[transc_df$AA_Pos==y] <- 1}
    else if(mean(transc_df$rASA[transc_df$AA_Pos==y],na.rm=T) > 0.5){transc_df$threep_disord[transc_df$AA_Pos==y] <- 1}
    else{break}
    }
  }

  transc_df <- transc_df %>% select(Chr,Pos,transcript,AminoAcid,orientation,fivep_disord,threep_disord)
  ###Now add STOP codon, and 5'/3' UTR zones of 12 bp for cushion
  if(transc_df$orientation[1]=="+"){
  	transc_df <- rbind.data.frame(
					cbind.data.frame(Chr=transc_df$Chr[1],Pos=min(transc_df$Pos,na.rm=T)-1:12,transcript=transc_df$transcript[1],AminoAcid="UTR",orientation="+",fivep_disord=1,threep_disord=0),
					transc_df,
					cbind.data.frame(Chr=transc_df$Chr[1],Pos=max(transc_df$Pos,na.rm=T)+1:3,transcript=transc_df$transcript[1],AminoAcid="STOP",orientation="+",fivep_disord=0,threep_disord=1),
					cbind.data.frame(Chr=transc_df$Chr[1],Pos=max(transc_df$Pos,na.rm=T)+4:15,transcript=transc_df$transcript[1],AminoAcid="UTR",orientation="+",fivep_disord=0,threep_disord=1)
					)
  					}
else if(transc_df$orientation[1]=="-"){
	transc_df <- rbind.data.frame(
                                        cbind.data.frame(Chr=transc_df$Chr[1],Pos=min(transc_df$Pos,na.rm=T)-4:15,transcript=transc_df$transcript[1],AminoAcid="UTR",orientation="-",fivep_disord=0,threep_disord=1),
					cbind.data.frame(Chr=transc_df$Chr[1],Pos=min(transc_df$Pos,na.rm=T)-1:3,transcript=transc_df$transcript[1],AminoAcid="STOP",orientation="-",fivep_disord=0,threep_disord=1),
                                        transc_df,
                                        cbind.data.frame(Chr=transc_df$Chr[1],Pos=max(transc_df$Pos,na.rm=T)+1:12,transcript=transc_df$transcript[1],AminoAcid="UTR",orientation="-",fivep_disord=1,threep_disord=0)
                                        )
					}

  return(transc_df)
}

Disorder_DF <- bind_rows(lapply(unique(AAFold$transcript),function(Tr) Flag_DisordTail(Tr)))

write.table(Disorder_DF,"DisorderedTails.txt",sep="\t",col.names = T,row.names = F,quote=F)

