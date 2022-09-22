#' --- 
#' Title: Circos plots to show parallel evolution along with allele frequency
#' Author: Lewis Fisher
#' Date: 19/08/2022
#' ---


# install.packages("ggplot2")
# install.packages("circlize")
# Important to use version 0.4.10 or 0.4.15
# Versions inbetween do not show labels correctly.

library(reshape)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(stringr)
library(dplyr)
library(gtools)
packageVersion("circlize")


# Import files containing coverage data
cov_dir<- "C:\\Users\\lew61\\OneDrive - The University of Liverpool\\PhD project - pseudomonas resistance\\Genomics\\increment_increase_evolution\\coverage_plots"
setwd(cov_dir)
# list.files()
cov_data<- read.csv("increment_coverage_medians.bed", header = F, sep = "\t")

# upload Breseq output
mut_dir<- "C:\\Users\\lew61\\OneDrive - The University of Liverpool\\PhD project - pseudomonas resistance\\Genomics\\increment_increase_evolution\\increment_gd_files"
setwd(mut_dir)
list.files()
mut_data<- read.csv("Increment_breseq_reannotation_PA.txt", header = T, sep = "\t")

# using regular expressions to construct metadata
mut_data$title = str_remove(mut_data$title, "_breseq")
mut_data$replicate<- str_split_fixed(mut_data$title, "_", 2)[,1]
mut_data$transfer<- str_split_fixed(mut_data$title, "_", 2)[,2]
mut_data$condition = str_extract(mut_data$replicate, "(EDTA)|(NAB)|(NE)|(TSB)")


# Not all sampples were easilt identified as intergentic. 
# Flanking genes of an intergenic mutation are delimited with a "/" 
# Mutations were labelled as "intergnic" or "gene" in metadata
intergenic_TF<- c()
for (i in 1:nrow(mut_data)){
  intergenic_i<- grepl("/", mut_data$locus_tag)
  intergenic_i<- intergenic_i[i]
  if (intergenic_i == TRUE){
    intergenic_TF<- append(intergenic_TF, "Intergenic")
  } else {intergenic_TF<-append(intergenic_TF, "Gene")}
}

# Only plotting genes/non-intergenic mutations
mut_data$Intergenic_Gene = intergenic_TF
mut_data<- subset(mut_data, Intergenic_Gene != "Intergenic")
unique(mut_data$pa_gene)

# Where possile use a gene name, if gene name not present use locus tag
combined_annotation<- c()
for (i in 1:nrow(mut_data)) {
  PA_gene<- mut_data$pa_gene[i]
  PA_number<- mut_data$pa_number[i]
  if (PA_gene == "Unknown" | PA_gene == "" | PA_gene == "unknown") {
    combined_annotation<- append(combined_annotation,PA_number)
  }else {combined_annotation<- append(combined_annotation, PA_gene)
  }
}
length(combined_annotation)
nrow(mut_data)
unique(mut_data$combined_annotation)
mut_data$combined_annotation = combined_annotation

# Save data 
write.table(mut_data, "Increment_breseq_combined_annotations_PA.csv", sep = ",", row.names = FALSE)


#########################
# Preparing circos plots both as individual timepoints and plots contianing all three timepoints

# EDTA<- subset(cov_data, V1 == "EDTA1 T1" | V1 == "EDTA1 T7" | V1 == "EDTA1 T14" )
# Subsetting time point T14 for plotting (endpoint)
EDTA<- subset(cov_data, V1 == "EDTA1 T14" )


sample_factors<- str_extract(unique(EDTA$V1), "(T[1-9]{1,2})")
# sectors<- c("T1", "T7","T14")
sectors<- c("T14")
sample_factors<- factor(sample_factors, levels = sectors)

# Obtaining the genome size of each sample and making a matrix for plots
sample_max_EDTA<- c()
for (i in 1:length(unique(EDTA$V1))) {
  sample_list<- unique(EDTA$V1)
  sample_sub<- subset(EDTA, V1 == sample_list[i])
  listtmp<- max(sample_sub$V3)
  
  sample_max_EDTA<- append(sample_max_EDTA, listtmp)
}

# Set directory to save plots
setwd(mut_dir)

png("EDTA1_para_plot.png", res = 1200, height = 8, width = 8, units = "in")

# Subset the necessary data of mutations, with a minimum frequency of 0.05
# mut_EDTA<- subset(mut_data, condition == "EDTA")
mut_EDTA<- subset(mut_data, condition == "EDTA" & transfer == "T14")
mut_EDTA<- subset(mut_EDTA, frequency >= 0.05)
mut_EDTA<-  mut_EDTA[c("transfer", "position_start","position_end","frequency", "combined_annotation", "type", "replicate")]


# Rename columns to become compatible with the circlize package
colnames(mut_EDTA)= c("chr", "start", "end", "value1", "gene", "type", "replicate")

# Ensure plots have been cleared before starting new plot
circos.clear()


# making a matrix with range in genome size for each timepoint
ref<- matrix(c(rep(0, length(sample_factors) ), sample_max_EDTA), ncol=2)

# Laying out how timepoints are laid out in the circle
circos.par("track.height"=0.8, "gap.degree"= 4, cell.padding=c(0, 0 ,0, 0), "start.degree" = 90)
circos.initialize(factors = sample_factors, xlim=ref)
col_text <- "grey40"
# Colours of each mutation type
border_colours<- c("SNP" = "darkorchid2", "INS" = "chartreuse3", "DEL" = "darkorange1", "SUB" = "blue4")

# Labelling each mutation at the locus in the genome
circos.genomicLabels(mut_EDTA, labels= mut_EDTA$gene, cex=0.3, side="outside", labels_height=0.1, padding = 0.5, niceFacing = TRUE, connection_height = 0.1, line_lwd = 0.1, 
                     col = mut_EDTA$type %>% str_replace_all(border_colours), 
                     line_col = mut_EDTA$type %>% str_replace_all(border_colours))

# Options
track_size_mut<- 0.06
pch_opt<- 16
cex_opt<- 0.3
labels_cex_opt<-0.2 
lwd_opt<- 0.21

# addition of 4000 bp to make lines visible on the plot
region_plus<- 4000

# Title in the centre of the plot
text(0, 0, "EDTA", cex = 1.5)
# Each track represents a single lineage.
# Mutations being coloured by mutation type 
# Each mutation was placed on a graph showing allele frequencies between 0 and 1 "value1" represents allele frequency
circos.genomicTrack(data=subset(mut_EDTA,  replicate == "EDTA1"), panel.fun = function(region, value, ...) {
  region$start = region$start +region_plus
  circos.genomicRect(region, value, type="l", col= "gray12", border = NA)
  circos.genomicPoints(region, value, pch = pch_opt, cex = cex_opt, col= value$type %>% str_replace_all(border_colours))
  circos.yaxis(side = "left", labels = TRUE, labels.cex = labels_cex_opt, lwd = lwd_opt)
}, track.height = track_size_mut, bg.border=F, bg.col= "gray95", ylim = c(0,1))

circos.genomicTrack(data=subset(mut_EDTA,  replicate == "EDTA2"), panel.fun = function(region, value, ...) {
  region$start = region$start +region_plus
  circos.genomicRect(region, value, type="l", col= "gray12", border = NA)
  circos.genomicPoints(region, value, pch = pch_opt, cex = cex_opt, col= value$type %>% str_replace_all(border_colours))
  circos.yaxis(side = "left", labels = TRUE, labels.cex = labels_cex_opt, lwd = lwd_opt)
}, track.height = track_size_mut, bg.border=F, bg.col= "gray95", ylim = c(0,1))

circos.genomicTrack(data=subset(mut_EDTA,  replicate == "EDTA3"), panel.fun = function(region, value, ...) {
  region$start = region$start +region_plus
  circos.genomicRect(region, value, type="l", col= "gray12", border = NA)
  circos.genomicPoints(region, value, pch = pch_opt, cex = cex_opt, col= value$type %>% str_replace_all(border_colours))
  circos.yaxis(side = "left", labels = TRUE, labels.cex = labels_cex_opt, lwd = lwd_opt)
}, track.height = track_size_mut, bg.border=F, bg.col= "gray95", ylim = c(0,1))

circos.genomicTrack(data=subset(mut_EDTA,  replicate == "EDTA4"), panel.fun = function(region, value, ...) {
  region$start = region$start +region_plus
  circos.genomicRect(region, value, type="l", col= "gray12", border = NA)
  circos.genomicPoints(region, value, pch = pch_opt, cex = cex_opt, col= value$type %>% str_replace_all(border_colours))
  circos.yaxis(side = "left", labels = TRUE, labels.cex = labels_cex_opt, lwd = lwd_opt)
}, track.height = track_size_mut, bg.border=F, bg.col= "gray95", ylim = c(0,1))


# genome scale and aesthetic
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=1, col=col_text, 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="slategray3", bg.border=F, track.height=0.06)

gen_scale <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.3)*10^6
circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
  circos.axis(h="bottom", major.at=gen_scale, labels=round(gen_scale/10^6, 1), labels.cex=0.4, 
              col=col_text, labels.col=col_text, lwd=0.7, labels.facing="outside", direction = "inside")
  circos.text(mean(CELL_META$xlim), -0.8, "Size (MB)", cex=0.6, col=col_text,facing="bending.inside", niceFacing=TRUE)
}, bg.border=F)

# Legend of mutation types by colour
lgd = Legend(at = c("SNP","INS","DEL", "SUB"), title = "Mutation Type", legend_gp = gpar(fill = c("darkorchid2","chartreuse3","darkorange1","blue4")))
lgd_list_vertical = packLegend(lgd)
draw(lgd_list_vertical, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
dev.off()




