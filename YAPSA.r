library(BSgenome.Hsapiens.UCSC.hg19)
library(YAPSA)
word_length <- 3

#准备 COSMIC_signature
Alex_COSMIC_signatures_path <- paste0("http://cancer.sanger.ac.uk/cancergenome/","assets/signatures_probabilities.txt")
AlexCosmicValid_sig_df <- read.csv(Alex_COSMIC_signatures_path,header=TRUE,sep="\t")
#load("D:\\Script\\R\\DataBase\\AlexCosmicValid.RData")
Alex_COSMIC_rownames <- paste(AlexCosmicValid_sig_df[,1],
                             AlexCosmicValid_sig_df[,2],sep=" ")
COSMIC_select_ind <- grep("Signature",names(AlexCosmicValid_sig_df))
AlexCosmicValid_sig_df <- AlexCosmicValid_sig_df[,COSMIC_select_ind]
number_of_Alex_COSMIC_sigs <- dim(AlexCosmicValid_sig_df)[2]
names(AlexCosmicValid_sig_df) <- gsub("Signature\\.","AC",
                                        names(AlexCosmicValid_sig_df))
rownames(AlexCosmicValid_sig_df) <- Alex_COSMIC_rownames


YASPA<-function(data,my_cutoff,key=key,dir=dir){
	data=unique(data)
	data<- translate_to_hg19(data,"CHROM")
	data1 <- 
	create_mutation_catalogue_from_df(
         data,
         this_seqnames.field = "CHROM", this_start.field = "POS",
         this_end.field = "POS", this_PID.field = "PID",
         this_subgroup.field = "SUBGROUP",
         this_refGenome = BSgenome.Hsapiens.UCSC.hg19,
         this_wordLength = word_length)
	data_df <- as.data.frame(data1$matrix)


	#（COSMIC_signature）
	COSMIC_signature_colour_vector <- c("green","pink","goldenrod",
                                    "lightblue","blue","orangered","yellow",
                                    "orange","brown","purple","red",
                                    "darkblue","magenta","maroon",
                                    "yellowgreen","violet","lightgreen",
                                    "sienna4","deeppink","darkorchid",
                                    "seagreen","grey","darkgrey",
                                    "black","yellow4","coral2","chocolate2",
                                    "navyblue","plum","springgreen")
	COSMIC_bio_process_vector <- c("spontaneous deamination","APOBEC",
                               "defect DNA DSB repair hom. recomb.",
                               "tobacco mutatgens, benzo(a)pyrene",
                               "unknown",
                               "defect DNA MMR, found in MSI tumors",
                               "UV light exposure","unknown","POL eta and SHM",
                               "altered POL E",
                               "alkylating agents, temozolomide",
                               "unknown","APOBEC","unknown",
                               "defect DNA MMR","unknown","unknown",
                               "unknown","unknown",
                               "associated w. small indels at repeats",
                               "unknown","aristocholic acid","unknown",
                               "aflatoxin","unknown","defect DNA MMR",
                               "unknown","unknown","tobacco chewing","unknown")
	AlexCosmicValid_sigInd_df <- data.frame(sig=colnames(AlexCosmicValid_sig_df))
	AlexCosmicValid_sigInd_df$index <- seq_len(dim(AlexCosmicValid_sigInd_df)[1])
	AlexCosmicValid_sigInd_df$colour <- COSMIC_signature_colour_vector
	AlexCosmicValid_sigInd_df$process <- COSMIC_bio_process_vector



	current_sig_df <- AlexCosmicValid_sig_df
	current_sigInd_df <- AlexCosmicValid_sigInd_df
	COSMICExposures_df <-LCD(data_df,current_sig_df)

	COSMIC_subgroups_df <-  make_subgroups_df(data,    COSMICExposures_df)

	 #换颜色绘图
	 pdf(paste0(dir,"/",key,".pdf"),width=10, height=7)
	#all
	p_all=exposures_barplot(in_exposures_df = COSMICExposures_df,   in_signatures_ind_df = current_sigInd_df,   in_subgroups_df = COSMIC_subgroups_df,in_title = key)


 
	#cutoff绘图
	my_cutoff <- my_cutoff
	general_cutoff_vector <- rep(my_cutoff,dim(current_sig_df)[2])
	CosmicValid_cutoffGen_LCDlist <- LCD_complex_cutoff(
		in_mutation_catalogue_df = data_df,
		in_signatures_df = current_sig_df,
		in_cutoff_vector = general_cutoff_vector,
		in_sig_ind_df = current_sigInd_df)
	
	#标准化,Unknown,6*5
	p_cut_bar_stand=exposures_barplot(
		in_exposures_df = CosmicValid_cutoffGen_LCDlist$norm_exposures,
		in_signatures_ind_df = CosmicValid_cutoffGen_LCDlist$out_sig_ind_df,
		in_subgroups_df = COSMIC_subgroups_df,in_title = key)
	
	
	dev.off()
}



data=read.csv(file="group.csv",header=T)
YASPA(data=data,my_cutoff=0.03,key="group",dir="./")

