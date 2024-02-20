##setup----
# libraries
library(tidyverse) # all-purpose data analysis tools
library(fs) # filepath manipulation
library(Matrix) # for working with matrices (especially sparse ones)
library(cowplot) #ggplot add on (save_plot)

# directories
wd <- getwd()
input_dir <- "/Volumes/nchevrier/katej/projects/polysome-seq/" # where tagcounts matrix live on project2 server
output_dir <- fs::path(wd, "input/data/2_formatted_data")# where prepped data will live (metadata/counts csvs etc)


#2024-02-16_sgMax----
counts_dir = fs::path(input_dir, "2024_02_15_MAXko_polyseq/bcbio/final/2024-02-18_bcbio")
counts = readMM(fs::path(counts_dir, 'tagcounts.mtx'))

# add rownames and colnames
rownames(counts) = read_lines(fs::path(counts_dir, 'tagcounts.mtx.rownames'))
colnames(counts) = read_lines(fs::path(counts_dir, 'tagcounts.mtx.colnames'))

##make metadata file----
# read in table with barcodes/samples used
barcode_table = read_csv(fs::path(output_dir, "../../other/2024_02_sgMax_polysome_barcodes.csv"))

# start with a table containing the colnames of the counts table
metadata = tibble(sample = colnames(counts)) %>%
  # trim to get just the barcode
  mutate(barcode = gsub(".*:","",sample)) %>%
  # rename columns to be compatible with barcode_table
  rename_all(~str_to_title(.)) %>%
  # merge with barcode_table to get only the barcodes used (36)
  inner_join(barcode_table, by = "Barcode") %>%
  # make treatment, rep, pool, cols
  separate(Name, into = c("treatment", "pool", "replicate"), sep = "_", remove = F) %>%
  # update Name col to remove spaces to make future manipulations easier
  mutate(Name = gsub(" ", "_", Name)) %>%
  # reorder in order of sample number
  arrange(Sample)

##barcode qc----
#make qc df
all_barcodes = data.frame(Reads = counts %>% colSums(),
                          Sample = colnames(counts)) %>%
  # create logical column indicating if barcode is used
  mutate(Used = Sample %in% metadata$Sample,
         Name = plyr::mapvalues(Sample, from = metadata$Sample, to = metadata$Name )) %>%
  arrange(desc(Reads)) %>%
  mutate(Sample = factor(Sample, levels = Sample))


p2<-ggplot(subset(all_barcodes, gsub(".*:", "", Sample) %in% metadata$Barcode), 
       aes(x = Name, y = Reads)) +
  geom_bar(stat = "identity") +
  #scale_y_log10() +
  labs(x = "Barcode") +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        legend.position = "none")

#save plot
save_plot(fs::path(output_dir,
                   "../../../output/2024_02_16_sgMax_LSG/QC_plots/20240220_reads_per_barcode_used_only.pdf"),p2, dpi=90)

##subset counts and convert to gene symbols----
# reorder to match metadata
counts = counts[, metadata$Sample]
# convert to non-sparse matrix
counts = as.matrix(counts)
# switch colnames to full name instead of barcode
colnames(counts) = metadata$Name

# read symbols information
symbols = read_csv(fs::path(output_dir, '../../other/bcbio_symbol_id.csv'))
# set up conversion vector
to_symbol = symbols$r_names
names(to_symbol) = symbols$id
# convert rownames to symbols
rownames(counts) = to_symbol[rownames(counts)]

## write out datafiles----
# save metadata
write.csv(metadata, fs::path(output_dir, "20240220_LSG_MaxKO_polysome_metadata.csv"))
# save counts
write.csv(counts, fs::path(output_dir, "20240220_LSG_MaxKO_polysome_counts.csv"))




#2023-10-31_LSG----
##read in tagcounts matrices----
counts_dir = path(input_dir, "2023_10_23_LSG_polyseq/bcbio/final/2023-10-31_bcbio")
counts = readMM(path(counts_dir, 'tagcounts.mtx'))

# add rownames and colnames
rownames(counts) = read_lines(path(counts_dir, 'tagcounts.mtx.rownames'))
colnames(counts) = read_lines(path(counts_dir, 'tagcounts.mtx.colnames'))

##make metadata file----
# read in table with barcodes/samples used
barcode_table = read_csv(path(output_dir, "../../other/2023_10_polysome_barcodes.csv"))

# start with a table containing the colnames of the counts table
metadata = tibble(sample = colnames(counts)) %>%
  # trim to get just the barcode
  mutate(barcode = gsub(".*:","",sample)) %>%
  # rename columns to be compatible with barcode_table
  rename_all(~str_to_title(.)) %>%
  # merge with barcode_table to get only the barcodes used (36)
  inner_join(barcode_table, by = "Barcode") %>%
  # make treatment, rep, pool, cols
  separate(Name, into = c("treatment", "replicate", "pool"), sep = " ", remove = F) %>%
  # create polysome fraction description
  mutate(polysome_frac = ifelse(pool == "Input", "Total",
                               ifelse(pool == "Pool1", "Light",
                                      ifelse(pool == "Pool2", "Medium", "Heavy")))) %>%
  # update Name col to remove spaces to make future manipulations easier
  mutate(Name = gsub(" ", "_", Name)) %>%
  # reorder in order of sample number
  arrange(Sample)

##barcode qc----
#make qc df
all_barcodes = data.frame(Reads = counts %>% colSums(),
                          Sample = colnames(counts)) %>%
                          # create logical column indicating if barcode is used
                          mutate(Used = Sample %in% metadata$Sample) %>%
                          arrange(desc(Reads)) %>%
                          mutate(Sample = factor(Sample, levels = Sample))
               
  
p<-ggplot(all_barcodes, aes(x = Sample, y = Reads, fill = Used)) +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  labs(x = "Barcode") +
  theme(axis.text.x = element_blank())
  
#save plot
save_plot(path(output_dir, "../../../output/QC_plots/reads_per_barcode.pdf"),p, dpi=90)

##subset counts and convert to gene symbols----
# reorder to match metadata
counts = counts[, metadata$Sample]
# convert to non-sparse matrix
counts = as.matrix(counts)
# switch colnames to full name instead of barcode
colnames(counts) = metadata$Name

# read symbols information
symbols = read_csv(path(output_dir, '../../other/bcbio_symbol_id.csv'))
# set up conversion vector
to_symbol = symbols$r_names
names(to_symbol) = symbols$id
# convert rownames to symbols
rownames(counts) = to_symbol[rownames(counts)]

## write out datafiles----
# save metadata
write.csv(metadata, path(output_dir, "20231101_LSG_polysome_metadata.csv"))
# save counts
write.csv(counts, path(output_dir, "20231101_LSG_polysome_counts.csv"))
