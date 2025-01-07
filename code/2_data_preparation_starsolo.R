##setup----
# libraries
library(tidyverse) # all-purpose data analysis tools
library(fs) # filepath manipulation
library(Matrix) # for working with matrices (especially sparse ones)
library(cowplot) #ggplot add on (save_plot)
library(readxl) #for reading excel docs
library(scales) #for sci notation

# directories
wd <- getwd()
input_dir <- "/Volumes/nchevrier/projects/katej/projects/polysome-seq/2024_12_CtuKO_polyseq/SCRB/GeneFull/Filtered" # where bcbio output and meta data matrices live on project2 server
output_dir <- fs::path(wd, "input/data/2_formatted_data/star_solo")# where prepped data will live (metadata/counts csvs etc)

##read in sample info for de-multiplexing----
sample_info <- read_xlsx(fs::path(wd, "input/other/202412_CtuKO_polysome-seq_barcodes.xlsx"))
#create meta data table
meta <- sample_info %>%
        separate(Sample, 
                 into = c("Target", "Guide", "Replicate", "Pool"),
                 sep = "-", remove = F) %>%
        filter(Pool != "input_lo")
                 
#read in counts matrices for DCs and T cells (combined matrix)---
#read in tagcounts matrices
counts = readMM(fs::path(input_dir, 'matrix.mtx')) 
# add rownames and colnames
rownames(counts) = read_lines(fs::path(input_dir, 'features.tsv'))
colnames(counts) = read_lines(fs::path(input_dir, 'barcodes.tsv'))

#Shorten rownames to only gene name (remove ENSEMBL ID info)
rownames(counts) = gsub(".*\t","",gsub("\tGene.*","",rownames(counts)))

#Replace barcodes with sample names----
#convert to non-sparse format
NS_counts = as.matrix(counts)

# reorder to match metadata & pull out barcodes used in exp
NS_counts = NS_counts[, meta$Barcode]

# switch colnames to full name instead of barcode
colnames(NS_counts) = plyr::mapvalues(colnames(NS_counts),  from = meta$Barcode, to = meta$Sample)


#reorder samples so replicates are clustered together for easier viz in morpheus
meta_order = meta %>% 
             mutate(Target = factor(Target, levels = c("sgCtrl", "sgCtu1", "sgCtu2"))) %>%
             arrange(Target, Pool, Guide, Replicate)
  
NS_counts_order = NS_counts[,meta_order$Sample_Name]
  
#save separate count & meta data in case its useful later----
write.csv(as.matrix(NS_counts_order), fs::path(output_dir, "20250107_CtuKO_polyseq_counts.csv"))
write.csv(as.matrix(meta_order), fs::path(output_dir, "20250107_CtuKO_polyseq_metadata.csv"))


#----checking unfiltered matrix to see # of reads for unused barcodes
NS_counts_unfilt = as.matrix(counts) 
 
full_meta = data.frame(Barcode = colnames(NS_counts_unfilt)) %>%
            mutate(Barcode_used = ifelse(Barcode %in% meta$Barcode, T, F))

read_sum = data.frame(Reads = colSums(NS_counts)) %>%
           rownames_to_column("Sample") %>%
           arrange(desc(Reads)) %>%
           mutate(Sample_Name = factor(Sample, levels = unique(Sample)),
                  Barcode_used = T)

#raw only has 24 barcodes in it so must be already filtered matrix...           
#save barcode QC plots
p<-ggplot(read_sum, aes(x = Sample_Name, y = Reads)) +
  scale_y_continuous(labels = scales::scientific) +
  xlab(NULL) +
  geom_col(color = "black", fill= "white") +
  geom_hline(color = "red", linetype=2, aes(yintercept = 1e6)) +
  #scale_fill_manual(values = c("F" ="black", "T" ="white")) +
  theme(axis.text.x = element_text(angle = 45, hjust =1),
        legend.position = "bottom")

p
#reads might be too low to get meaninful DE data...
save_plot(fs::path(output_dir, "../../../../output/QC/read_by_sample_used.pdf"), p,
          base_aspect_ratio = 1.75)
