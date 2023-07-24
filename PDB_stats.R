args<- commandArgs(trailingOnly = T)

library(dplyr)
library(data.table)
library(seqinr)
library(stringr)

print(args[1])
##Correct for the fact that PDB files are missing "-" option
PDBfile <- args[1]
#Reading in PDB file with extraneous header is difficult.  I just chop it down to the important ATOM part
system(paste0("zcat PDB/",PDBfile,".pdb.gz | grep ATOM > PDB/",PDBfile,".pdb")) 
PDB <- fread(paste0("PDB/",PDBfile,".pdb"),skip=0,fill=T)
print("PDB_files read in")
colnames(PDB) <- c("Type","AtomNum","AtomChain","AminoAcid","INFO1","AA_Pos","X","Y","Z","INFO2","pLDDT","element")

#Just in case that didn't grab everything correctly, filtere down to ATOM part of the table
PDB <- PDB[PDB$Type=="ATOM",]
#Summarize pLDDT for all atoms in each amino acid
PDB <- PDB %>% group_by(AA_Pos,AminoAcid) %>% summarize(mean_pLDDT=mean(as.numeric(pLDDT))) 
#Add a column for the total average for the gene
PDB <- cbind.data.frame(PDB,Gene_pLDDT=mean(PDB$mean_pLDDT))


###Let's run dssp to get rASA
system(paste0("mkdssp -i PDB/",PDBfile,".pdb | grep -A 100000 \"RESIDUE AA\" | cut -c 1-5,13-14,35-39 | grep -v \"\\!\" | awk '{print $1,$2,$3}' > dssp/",args[1],".dssp"))
dsspfile <- read.table(paste0("dssp/",args[1],".dssp"),skip=1) 
colnames(dsspfile) <- c("AA_Pos","AminoAcid_1let", "ASA")
#Need to calculate relative ASA https://en.wikipedia.org/wiki/Relative_accessible_surface_area
#Will use Tien et al Theoretical maximum ASA 2013
Max_ASA_table <- cbind.data.frame(AA=c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'),
		  max_ASA=c(129.0, 274.0, 195.0, 193.0, 167.0, 223.0, 225.0, 104.0, 224.0, 197.0, 201.0, 236.0, 224.0, 240.0, 159.0, 155.0, 172.0, 285.0, 263.0, 174.0))

Get_MaxASA <- function(protAA,MASA_Table){
	return(MASA_Table$max_ASA[as.character(MASA_Table$AA)==as.character(protAA)])
	}	
dsspfile <- cbind.data.frame(dsspfile, maxASA=unlist(sapply(dsspfile$AminoAcid_1let,function(x) Get_MaxASA(x,Max_ASA_table))))
dsspfile <- dsspfile %>% mutate(rASA=as.numeric(ASA)/as.numeric(maxASA))


##Read in CDS information from gff and adjust for 1-based coordinates
gff <- read.table(paste0("gff/",args[1],".gff"))
gff_adj <- gff
gff_adj[,5] <- gff_adj[,5] - min(gff[,4]) + 1
gff_adj[,4] <- gff_adj[,4] - min(gff[,4]) + 1

#initiate vector for codon and basepair positions to match up with amino acid positions
basepair_vec <- rep(1,nrow(PDB)*3+3)
CodonPos_vec <- rep(1,nrow(PDB)*3+3)
AAPos_vec <- rep(NA,nrow(PDB)*3+3)

#Reads in CDS sequence
FastaCDS <- as.vector(unlist(read.fasta(paste0("CDS/",args[1],".cds.fasta"))))




##Use gff to count out genomic positions of genes, AA, and codons
Codon_pos <- 1
AA_pos <- 1
if(gff_adj[1,7]=="+")
{
        CDS_gff <- gff[gff[,3]=="CDS",]
        bp_pos <- 1
        for(x in 1:nrow(CDS_gff))
        {
                for(gffpos in CDS_gff[x,4]:CDS_gff[x,5])
                {
		AAPos_vec[bp_pos] <- AA_pos
                basepair_vec[bp_pos] <- gffpos
                CodonPos_vec[bp_pos] <- Codon_pos %% 3
		if(Codon_pos %% 3==0){AA_pos <- AA_pos +1}
                Codon_pos <- Codon_pos + 1
                bp_pos <- bp_pos + 1
                }
        }
}


if(gff_adj[1,7]=="-")
{
        CDS_gff <- gff[gff[,3]=="CDS",]
        bp_pos <- 1
        for(x in 1:nrow(CDS_gff))
        {
                for(gffpos in CDS_gff[x,5]:CDS_gff[x,4])
                {
		AAPos_vec[bp_pos] <- AA_pos
                basepair_vec[bp_pos] <- gffpos
                CodonPos_vec[bp_pos] <- Codon_pos %% 3
		if(Codon_pos %% 3==0){AA_pos <- AA_pos +1}
                Codon_pos <- Codon_pos + 1
                bp_pos <- bp_pos + 1
                }
        }
}
        

	CodonPos_vec[CodonPos_vec==0] <- 3
#Create df with CDS information paired with AA information
	BPDF <- cbind.data.frame(gff[1,1],basepair_vec,args[1],RefAllele=toupper(FastaCDS),AAPos_vec,CodonPos_vec)
	colnames(BPDF) <- c("Chr","Pos","transcript","RefAllele","AA_Pos","CodonPos")
# Join with PDB file to correlate bp and AA to folding confidence
	BPDF <- left_join(BPDF,PDB)
#Add in stop codons
	BPDF$AminoAcid[is.na(BPDF$AminoAcid)] <- "STOP"
	BPDF$mean_pLDDT[is.na(BPDF$mean_pLDDT)] <- 100
	BPDF$mean_pLDDT <- round(BPDF$mean_pLDDT,2)
	BPDF$Chr <- as.character(BPDF$Chr)
#Join with rASA calculations
	BPDF <- inner_join(BPDF,dsspfile)
#Label sites as 4 fold degenerate
isdegenerate <- function(BPDF){
FourDsites <- rep(0,nrow(BPDF))
        for(p_site in 1:(nrow(BPDF)/3)){
                      if(BPDF$RefAllele[3*p_site-2]=="G" & BPDF$RefAllele[3*p_site-1]=="C") {FourDsites[3*p_site] <- 1}
                 else if(BPDF$RefAllele[3*p_site-2]=="C" & BPDF$RefAllele[3*p_site-1]=="G") {FourDsites[3*p_site] <- 1}
                 else if(BPDF$RefAllele[3*p_site-2]=="G" & BPDF$RefAllele[3*p_site-1]=="G") {FourDsites[3*p_site] <- 1}
                 else if(BPDF$RefAllele[3*p_site-2]=="C" & BPDF$RefAllele[3*p_site-1]=="T") {FourDsites[3*p_site] <- 1}
                 else if(BPDF$RefAllele[3*p_site-2]=="T" & BPDF$RefAllele[3*p_site-1]=="C") {FourDsites[3*p_site] <- 1}
                 else if(BPDF$RefAllele[3*p_site-2]=="A" & BPDF$RefAllele[3*p_site-1]=="C") {FourDsites[3*p_site] <- 1}
                 else if(BPDF$RefAllele[3*p_site-2]=="G" & BPDF$RefAllele[3*p_site-1]=="T") {FourDsites[3*p_site] <- 1}
            }
return(FourDsites)
}

FourD_vector <- isdegenerate(BPDF)
BPDF <- cbind.data.frame(BPDF,FourFoldDegenerate=FourD_vector)
BPDF <- BPDF %>% mutate(AA_ProportionPosition= AA_Pos/max(AA_Pos))

#We write this file out because you may interested in the invarable portion of the proteins as well
write.table(BPDF,paste0("prot_files/",args[1],".prt"),quote=F,col.names=T,row.names=F,sep="\t")
### Read in Freq File
### This Freq File is going to be modified to include snpEFF annotation as a column, and for simplicity 
### be of the form "Chr, Pos, snpEFF, alternate Alleles frequencies. . . .
freq <- fread(paste0("Freq/",args[1],".frq"), fill=T, sep="\t")
freq[freq==""] <- NA
freq[freq=="NA"] <- NA


#Turn multiallelic sites into a long format dataframe
spreadFreq <- function(freqsite){
  
  alleles <- sapply(freqsite[4:length(freqsite)], function(x) strsplit(x,":")[[1]][1])
  alleles <- alleles[!is.na(alleles)]
  freqs <- sapply(freqsite[4:length(freqsite)], function(x) strsplit(x,":")[[1]][2])
  freqs <- freqs[!is.na(freqs)]
  freqdf <- as.data.frame(matrix(NA,nrow=0,ncol=5))
  for(allelenum in 1:length(alleles)){
    freqdf <- rbind.data.frame(freqdf,c(freqsite[1:3],alleles[allelenum],freqs[allelenum]))
  }
  colnames(freqdf) <- c("Chr",  "Pos","snpEFF" ,"alleles","freqs")
  return(freqdf)
}
freq <- bind_rows(apply(freq,1,function(y) spreadFreq(y)), .id = "column_label")
freq <- freq %>% select(-column_label) %>% filter(!grepl("\\*",alleles)) %>% filter(!is.na(alleles)) %>% filter(freqs!=0 & freqs!=1) %>% mutate(MAF=ifelse(as.numeric(as.character(freqs)) < 0.5,as.numeric(as.character(freqs)),1-as.numeric(as.character(freqs))))


#These functions are used to parse the alternate allele and its snpEFF effect
ParseAllele <- function(SNPEFFstring){
  AlleleIdent <- gsub(".*ANN=","",SNPEFFstring)
  AlleleIdent <- strsplit(AlleleIdent,"\\|")[[1]][1]
  AlleleIdent <- unname(gsub("*ANN=","",AlleleIdent))
  return(AlleleIdent)
}

Get_snpEFF_effect <- function(freqsite_parsed,gene){
  LoF_Vec <- strsplit(freqsite_parsed[3],split = ",")[[1]]
  LoF_Vec_gene <- LoF_Vec[grepl(gene,LoF_Vec)]
  AllIdent <- sapply(LoF_Vec_gene, function(x) ParseAllele(x))
  AllMatch <- match(freqsite_parsed[4],AllIdent)
  LoF_del <- LoF_Vec_gene[AllMatch] 
  Effect <- strsplit(LoF_del,split="\\|")[[1]][3]
  return(Effect)
}

Get_Alt_AA <- function(freqsite_parsed,gene){
  LoF_Vec <- strsplit(freqsite_parsed[3],split = ",")[[1]]
  LoF_Vec_gene <- LoF_Vec[grepl(gene,LoF_Vec)]
  AllIdent <- sapply(LoF_Vec_gene, function(x) ParseAllele(x))
  AllMatch <- match(freqsite_parsed[4],AllIdent)
  LoF_del <- LoF_Vec_gene[AllMatch] 
  Effect <- strsplit(LoF_del,split="\\|")[[1]][3]
  if(is.na(Effect)){return(NA)}
  else if(Effect=="MODERATE"){
  AltAA <- strsplit(LoF_del,split="\\|")[[1]][11]
  return(toupper(str_sub(AltAA,-3,-1))  )
  }
  else{return(NA)}
}
freq <- cbind.data.frame(freq,Effect=apply(freq,1,function(x) Get_snpEFF_effect(x,args[1])),Alt_AA=apply(freq,1,function(x) Get_Alt_AA(x,args[1])))
freq$Pos <- as.numeric(as.character(freq$Pos))

#Join PDB/AA info with frequency/snpEFF information
BPDF <- inner_join(BPDF,freq)

# Create a matrix with uppercase amino acid abbreviations and functional classes
amino_acids <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
                 "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
amino_acids <- toupper(amino_acids)  # Convert amino acid abbreviations to uppercase
functional_classes <- c("Nonpolar", "Basic", "Polar", "Acidic", "Polar", "Polar", "Acidic",
                        "Nonpolar", "Basic", "Nonpolar", "Nonpolar", "Basic", "Nonpolar",
                        "Nonpolar", "Nonpolar", "Polar", "Polar", "Nonpolar", "Polar",
                        "Nonpolar")

amino_acid_matrix <- cbind.data.frame(AminoAcid=as.character(amino_acids), Ref_FunctionalClass=as.character(functional_classes))
BPDF$AminoAcid <- as.character(BPDF$AminoAcid)
BPDF <- inner_join(BPDF,amino_acid_matrix)
colnames(amino_acid_matrix) <- c("Alt_AA","Alt_FunctionalClass")
BPDF <- left_join(BPDF,amino_acid_matrix)
#remove snpEFF string to improve readability
BPDF <- BPDF %>% select(-snpEFF)
write.table(BPDF,paste0("pLDDT/",args[1],".AA.txt"),sep="\t",row.names=F,col.names=T,quote=F)

