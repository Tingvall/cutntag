#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run cutntag.nf --samples 20210121_CnT_WT_TH_BM_proB/samples_subset.txt --outdir 20210121_CnT_WT_TH_BM_proB/test --bowtie2_index /Users/johti53/bioinformatics_tools/reference_files/index/bowtie2/mm10/mm10.GRCm38.p5_gencode
  Mandatory arguments:
    --input [file]                  Comma-separated file containing information about the samples in the experiment (see docs/usage.md) (Default: './design.csv')

  """.stripIndent()
  }


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * VALIDATE INPUTS
 */

if (params.samples)     { ch_samples = Channel.fromPath(params.samples, checkIfExists: true) } else { exit 1, 'Samples not specified' }
    ch_samples
        .splitCsv(header:true, sep:'\t')
        .map { row -> [ row.sample_id,  row.group, row.replicate, [file(row.sample_path+'/*R1*'), file(row.sample_path+'/*R2*')]] }
        .into { ch_samples_split1; ch_samples_split2}

ch_samples_split1
  .map { n -> tuple(n[0], n[3]) }
  .set{ch_samples_fastq}

ch_samples_split2
  .map { n -> tuple(n[0], n[1], n[2]) }
  .into{ch_samples_info_1; ch_samples_info_2}


ch_bowtie2_index = Channel
        .fromPath("${params.bowtie2_index}*", checkIfExists: true)
        .ifEmpty { exit 1, "Bowtie2 index directory not found: " }

if (params.design)     { ch_design = Channel.fromPath(params.design, checkIfExists: true) } else { exit 1, 'Experiment design not specified' }
    ch_design
        .splitCsv(header:true, sep:'\t')
        .map { row -> [ row.sample_id, row.control_id, row.group, row.replicate] }
        .into { ch_design_seacr; ch_design_peak_stats}

if (params.chromsize)     { ch_chromsize = Channel.fromPath(params.chromsize, checkIfExists: true) } else { exit 1, 'Chromosome sizes not specified' }



println ("""
        ===========================================================================================
                                                  CUT&RUN
        ===========================================================================================
        Samples: ${params.samples}
        Design: ${params.design}
        Outdir: ${params.outdir}
        ===========================================================================================
        """)



/*
 * 1. Merge, unzip and move fastq files
 */
process PREPARE_FASTQ {

    input:
    set val(sample_id), val(fastq) from ch_samples_fastq

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq"), path("${sample_id}_R2.fastq")  into ch_fastq_untrimmed_fastqc, ch_fastq_untrimmed_trimgalore, ch_fastq_untrimmed_bowtie2


    script:
    """
    #Unzip and cat the files.
    gunzip -c ${fastq[0].join(' ')} > ${sample_id}_R1.fastq
    gunzip -c ${fastq[1].join(' ')} > ${sample_id}_R2.fastq
    """
}


/*
* 2. FastQC - Untrimmed fastq files
*/
process FASTQC {
publishDir "${params.outdir}/FastQC/Untrimmed", mode: 'copy'

when:
!params.skip_fastqc

input:
tuple val(sample_id), path(fastq_R1), path(fastq_R2) from ch_fastq_untrimmed_fastqc

output:
path '*.{zip,html}' into ch_fastqc_untrimmed

script:
"""
fastqc -q -t $task.cpus ${fastq_R1}
fastqc -q -t $task.cpus ${fastq_R2}
"""
}


/*
* 3. TrimGalore - Read trimming
*/
process TRIM {
publishDir "${params.outdir}/FastQC/Trimmed", mode: 'copy', pattern: '*.{zip,html}'

when:
!params.skip_trim

input:
tuple val(sample_id), path(fastq_R1), path(fastq_R2) from ch_fastq_untrimmed_trimgalore

output:
tuple val(sample_id), path('*R1*.gz'), path('*R2*.gz') into ch_fastq_trimmed_bowtie2
path '*.txt' into ch_trimgalore_report
path '*.{zip,html}' into ch_fastqc_trimmed

script:
def cores = 1
if (task.cpus) {
	cores = (task.cpus as int) -4
	if (cores < 1) cores = 1
	if (cores > 4) cores = 4
        }
"""
trim_galore --paired --gzip --cores $cores --fastqc --clip_R1 ${params.clip_R1} --clip_R2 ${params.clip_R2} --three_prime_clip_R1  ${params.three_prime_clip_R1} --three_prime_clip_R2 ${params.three_prime_clip_R2} ${fastq_R1} ${fastq_R2}
"""
}

/*
* 4. Bowtie2 - Alignment
*/
process ALIGN {
publishDir "${params.outdir}/Alignment/Bowtie2_reports", mode: 'copy', pattern: '*_bowtie2_report.txt'

input:
path bowtie2_index from  ch_bowtie2_index.collect()
tuple val(sample_id), path(fastq_R1), path(fastq_R2) from ch_fastq_trimmed_bowtie2

output:
tuple val(sample_id), path("${sample_id}_bowtie2.sam") into ch_sam
tuple val(sample_id), path("${sample_id}_bowtie2_report.txt") into ch_bowtie_report
//path "${sample_id}_bowtie2_report.txt" into ch_bowtie_report

script:
index_base = bowtie2_index[0].toString() - ~/.rev.\d.bt2?/ - ~/.\d.bt2?/
"""
bowtie2  --end-to-end --very-sensitive --no-discordant -I 10 -X 700  --no-mixed --phred33 -p $task.cpus -x ${index_base} -1 ${fastq_R1} -2 ${fastq_R2} -S ${sample_id}_bowtie2.sam &> ${sample_id}_bowtie2_report.txt

"""
}

/*
* 4. Remove duplicates
*/
process REMOVE_DUPLICATES {
publishDir "${params.outdir}/Sample_info/Duplicate_removal", mode: 'copy', pattern: '*_picard_rmDup_metrics.txt'

input:
tuple val(sample_id), path(sam) from ch_sam

output:
tuple val(sample_id), path("${sample_id}_bowtie2_rmDup.sam") into ch_sam_rmDup
tuple val(sample_id), path("*_picard_rmDup_metrics.txt") into ch_picard_metrics
tuple val(sample_id), path("${sample_id}_fragment_size.txt") into ch_fragment_size


script:
"""
# Sort
picard SortSam I=${sam} O=${sample_id}_bowtie2_sorted.sam SORT_ORDER=coordinate

# Remove duplicates
picard MarkDuplicates I=${sample_id}_bowtie2_sorted.sam O=${sample_id}_bowtie2_rmDup.sam REMOVE_DUPLICATES=true METRICS_FILE=${sample_id}_picard_rmDup_metrics.txt

# Fragment size distribution
samtools view -F 0x04 ${sample_id}_bowtie2_rmDup.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs(\$9)}' | sort | uniq -c | awk -v OFS="\t" '{print \$2, \$1/2}' > ${sample_id}_fragment_size.txt

"""
}


/*
* 5. Alignment plots & Fragment length
*/
process SAMPLE_INFO {
publishDir "${params.outdir}/Sample_info", mode: 'copy'

when:
!params.skip_plots

input:
val align_report2 from ch_samples_info_1.join(ch_bowtie_report).join(ch_picard_metrics).join(ch_fragment_size).collect()

output:
path "Align_RmDup_summary.txt" into ch_align_redup_summary
path "Alignment_RmDuplicates.pdf" into ch_align_plot
path "Fragment_Length.pdf" into ch_fraglen_plot


script:
"""
#!/usr/bin/env Rscript

packages <- c("dplyr", "stringr", "ggplot2", "viridis", "cowplot")
lapply(packages, library, character.only = TRUE)

# Read sample info
sample_info <-   print("${align_report2}")
sample_info <- str_replace_all(sample_info, "\\\\[|\\\\]", "")
sample_info<- read.table(text=sample_info,col.names=c('Sample','Group','Replicate', 'Align_stats', 'RmDup_stats', 'Frgament_size'), sep=",")

# Summarize alignemnt stats
align_result=c()
for(sample in sample_info\$Sample){
bowtie_stats = read.table(gsub(" ", "", toString(sample_info[sample_info\$Sample== sample,'Align_stats'])), header = FALSE, fill = TRUE)
align_rate = substr(bowtie_stats\$V1[6], 1, nchar(as.character(bowtie_stats\$V1[6]))-1)
align_result = data.frame(Sample = sample_info[sample_info\$Sample== sample,'Sample'],
                          Group = sample_info[sample_info\$Sample== sample,'Group'],
                          Replicate= as.character(sample_info[sample_info\$Sample== sample,'Replicate']),
                          Sequencing_Depth = bowtie_stats\$V1[1] %>% as.character %>% as.numeric,
                          Mapped_Fragments = bowtie_stats\$V1[4] %>% as.character %>% as.numeric + bowtie_stats\$V1[5] %>% as.character %>% as.numeric,
                          Alignment_Rate = align_rate %>% as.numeric)  %>% rbind(align_result, .)
}
align_result\$Sample = factor(align_result\$Sample)
align_result %>% mutate(Alignment_Rate = paste0(Alignment_Rate, "%"))


# Summarize duplication stats
rmdup_result2 = c()
for(sample in sample_info\$Sample){
  rmdup_result = read.table(gsub(" ", "", toString(sample_info[sample_info\$Sample==sample,'RmDup_stats'])), header = TRUE, fill = TRUE)
  rmdup_result2 = data.frame(Sample = sample_info[sample_info\$Sample== sample,'Sample'],
                            Mapped_Fragments = rmdup_result\$READ_PAIRS_EXAMINED[1] %>% as.character  %>% as.numeric,
                            Duplication_Rate = rmdup_result\$PERCENT_DUPLICATION[1] %>% as.character  %>% str_replace(",", ".")   %>% as.numeric * 100,
                            Estimated_LibrarySize = rmdup_result\$ESTIMATED_LIBRARY_SIZE[1] %>% as.character %>%  as.numeric) %>%
                            mutate(Unique_Fragments = Mapped_Fragments * (1-Duplication_Rate/100)) %>% rbind(rmdup_result2, .)
}
rmdup_result2\$Sample = factor(rmdup_result2\$Sample)


## Merge alignemtn and dpucliation stats
align_rmdup_result = left_join(align_result, rmdup_result2, by = c("Sample", "Mapped_Fragments")) %>% mutate(Duplication_Rate = paste0(Duplication_Rate, "%"))
write.table(align_rmdup_result, "Align_RmDup_summary.txt", col.names = T, row.names = F, sep = "\t", quote = FALSE)

# Plot
fig_aln1 = align_rmdup_result %>% ggplot(aes(x = Group, y = Sequencing_Depth/1000000, fill = Group)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_discrete(labels = NULL,breaks=NULL)+
  theme_bw(base_size = 10) +
  theme(legend.position="none")+
  ylab("Sequencing Depth per Million") +
  xlab("") +
  ggtitle("A. Sequencing Depth")


fig_aln2 = align_rmdup_result %>% ggplot(aes(x = Group, y = Mapped_Fragments/1000000, fill = Group)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_discrete(labels = NULL,breaks=NULL)+
  theme_bw(base_size = 10) +
  theme(legend.position="none")+
  ylab("Mapped Fragments per Million") +
  xlab("") +
  ggtitle("B. Alignable Fragment")

fig_aln3 = align_rmdup_result %>% ggplot(aes(x = Group, y = Alignment_Rate, fill = Group)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_discrete(labels =NULL,breaks=NULL)+
  theme_bw(base_size = 10) +
  theme(legend.position="none")+
  ylab("% of Mapped Fragments") +
  xlab("") +
  ggtitle("C. Alignment Rate")

fig_dup1 = align_rmdup_result %>% ggplot(aes(x = Group, y = Duplication_Rate, fill = Group)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_discrete(labels = NULL,breaks=NULL)+
  theme_bw(base_size = 10) +
  theme(legend.position="none")+
  ylab("Duplication Rate (%)") +
  xlab("")+
  ggtitle("D. Duplication rate")

fig_dup2 = align_rmdup_result %>% ggplot(aes(x = Group, y = Estimated_LibrarySize, fill = Group)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_discrete(labels = NULL,breaks=NULL)+
  theme_bw(base_size = 10) +
  theme(legend.position="none")+
  ylab("Estimated Library Size") +
  xlab("")+
  ggtitle("E. Library Size")

fig_dup3 = align_rmdup_result %>% ggplot(aes(x = Group, y = Unique_Fragments, fill = Group)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_discrete(labels = NULL, breaks=NULL)+
  theme_bw(base_size = 10) +
  theme(legend.position="none")+
  ylab("Number of Unique Fragments") +
  xlab("")+
  ggtitle("F. Unique Fragments")


legend <- get_legend( fig_aln1 + guides(color = guide_legend(nrow = 1)) +theme(legend.position = "bottom"))

plots <- plot_grid(
  fig_aln1 + theme(legend.position="none"),
  fig_aln2 + theme(legend.position="none"),
  fig_aln3 + theme(legend.position="none"),
  fig_dup1 + theme(legend.position="none"),
  fig_dup2 + theme(legend.position="none"),
  fig_dup3 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 2
)
pdf(file="Alignment_RmDuplicates.pdf", width=12, height=8)
plot_grid(plots, legend, ncol = 1, rel_heights = c(1, .1))
dev.off()


## Fragemtn length
Fragment_Length = c()
for(sample in sample_info\$Sample){
  Fragment_Length = read.table(gsub(" ", "", toString(sample_info[sample_info\$Sample == sample,'Frgament_size'])), header = FALSE) %>%
      mutate(Fragment_Length = V1 %>% as.numeric, Fragment_Count = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)),
      Sample = sample_info[sample_info\$Sample== sample,'Sample'],
      Group = sample_info[sample_info\$Sample== sample,'Group'],
      Replicate = sample_info[sample_info\$Sample== sample,'Replicate']) %>% rbind(Fragment_Length, .)
}

Fragment_Length\$Sample = factor(Fragment_Length\$Sample)
Fragment_Length\$Group = factor(Fragment_Length\$Group)
Fragment_Length\$Replicate = factor(Fragment_Length\$Replicate)

## Generate the fragment size density plot (violin plot)
fig_frlen1 = Fragment_Length %>% ggplot(aes(x = Sample, y = Fragment_Length, weight = Weight, fill = Group)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_discrete(labels = NULL, breaks=NULL)+
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")

fig_frlen2 = Fragment_Length %>% ggplot(aes(x = Fragment_Length, y = Fragment_Count, color = Group, group = Group, linetype = Replicate)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length (bp)") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

plots2 <- plot_grid(
  fig_frlen1,
  fig_frlen2,
  align = 'vh',
  hjust = -1,
  nrow = 2
)
pdf(file="Fragment_Length.pdf", width=12, height=12)
plot_grid(plots2, ncol = 1, rel_heights = c(1, .1))
dev.off()
"""
}


/*
* 6. Filter mapped reads & conversion
*/
process FILTER {
publishDir "${params.outdir}/Alignment/Bam", mode: 'copy',pattern: '*.bam'
publishDir "${params.outdir}/Alignment/Bedgraph", mode: 'copy',pattern: '*.bedgraph'

input:
tuple val(sample_id), path(sam_rmDup) from ch_sam_rmDup

output:
tuple val(sample_id), path("${sample_id}_bowtie2_rmDup_q${params.quality_score}.bam") into ch_bam_peakstat, ch_bam_bw
tuple val(sample_id), path("${sample_id}_bowtie2_rmDup_q${params.quality_score}_filt_fragemnts.bed") into ch_fragment_bed
tuple path("${sample_id}_bowtie2_rmDup_q${params.quality_score}.bedgraph"), val(sample_id) into ch_bedgraph_1, ch_bedgraph_2
tuple val(sample_id), path("${sample_id}_bowtie2_rmDup_q${params.quality_score}_fragments_count_bin${params.bin_size}.bed") into ch_correlation

script:
"""
# Filter based on quality score
samtools view -@ $task.cpus -h -q ${params.quality_score} ${sam_rmDup} -o ${sample_id}_bowtie2_rmDup_q${params.quality_score}.sam

# Keep the mapped read pairs and convert to bam
samtools view -@ $task.cpus -b -F 0x04 ${sample_id}_bowtie2_rmDup_q${params.quality_score}.sam -o ${sample_id}_bowtie2_rmDup_q${params.quality_score}.bam

# Convert into bed file format
bedtools bamtobed -i ${sample_id}_bowtie2_rmDup_q${params.quality_score}.bam -bedpe > ${sample_id}_bowtie2_rmDup_q${params.quality_score}.bed

# Keep read pairs on the same chromosome and fragment length < 1000bp.
awk '\$1==\$4 && \$6-\$2 < 1000 {print \$0}' ${sample_id}_bowtie2_rmDup_q${params.quality_score}.bed > ${sample_id}_bowtie2_rmDup_q${params.quality_score}_filt.bed

# Only extract the fragment related columns
cut -f 1,2,6 ${sample_id}_bowtie2_rmDup_q${params.quality_score}_filt.bed | sort -k1,1 -k2,2n -k3,3n  > ${sample_id}_bowtie2_rmDup_q${params.quality_score}_filt_fragemnts.bed

# Convert to bedgraph information
bedtools genomecov -bg -i ${sample_id}_bowtie2_rmDup_q${params.quality_score}_filt_fragemnts.bed -g ${params.chromsize} > ${sample_id}_bowtie2_rmDup_q${params.quality_score}.bedgraph

# For correlation plot
awk -v w=${params.bin_size} '{print \$1, int((\$2 + \$3)/(2*w))*w + w/2}' ${sample_id}_bowtie2_rmDup_q${params.quality_score}_filt_fragemnts.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print \$2, \$3, \$1}' |  sort -k1,1V -k2,2n  > ${sample_id}_bowtie2_rmDup_q${params.quality_score}_fragments_count_bin${params.bin_size}.bed
"""
}


/*
* 7. Sample correlation
*/
process CORRELATION {
publishDir "${params.outdir}/Sample_info", mode: 'copy'

when:
!params.skip_correlation

input:
val correlation from ch_samples_info_2.join(ch_correlation).collect()

output:

path "Sample_Correlation.pdf" into ch_correlation_plot


script:
"""
#!/usr/bin/env Rscript

packages <- c("stringr", "corrplot", "dplyr")
lapply(packages, library, character.only = TRUE)

# Read sample info
sample_info <-   print("${correlation}")
sample_info <- str_replace_all(sample_info, "\\\\[|\\\\]", "")
sample_info<- read.table(text=sample_info,col.names=c('Sample','Group','Replicate', 'Fragment_Count'), sep=",")

# Plot sample correlation
corr = c()
Fragment_Count = NULL
for(sample in sample_info\$Sample){

  if(is.null(Fragment_Count)){

    Fragment_Count = read.table(gsub(" ", "", toString(sample_info[sample_info\$Sample == sample,'Fragment_Count'])), header = FALSE)
    colnames(Fragment_Count) = c("chrom", "bin", sample)

  }else{

    Fragment_Count_tmp = read.table(gsub(" ", "", toString(sample_info[sample_info\$Sample == sample,'Fragment_Count'])), header = FALSE)
    colnames(Fragment_Count_tmp) = c("chrom", "bin", sample)
    Fragment_Count = full_join(Fragment_Count, Fragment_Count_tmp, by = c("chrom", "bin"))

  }
}

M = cor(Fragment_Count %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs")

pdf(file="Sample_Correlation.pdf")
corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,
         cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1,
         col = colorRampPalette(c("midnightblue","white","darkred"))(100))
dev.off()
"""
}


/*
* 8.Peak calling using SEACR
*/
process PEAK_CALLING {
publishDir "${params.outdir}/Peak_calling/Peaks_control", mode: 'copy', pattern: '*_seachr_peaks_ctrl*'
publishDir "${params.outdir}/Peak_calling/Peaks_top001", mode: 'copy', pattern: '*_seachr_peaks_top0.01*'

input:
set val(sample_id), val(control_id), val(group), val(replicate), path(sample_bedgraph), path(control_bedgraph) from ch_design_seacr.join(ch_bedgraph_2, by: 1).join(ch_bedgraph_1, by:1)

output:
tuple val(sample_id), path("${sample_id}_seachr_peaks_ctrl*") into ch_peak_ctrl
tuple val(sample_id), path("${sample_id}_seachr_peaks_top0.01*") into ch_peak_top001
script:
"""
SEACR_1.3.sh ${sample_bedgraph} ${control_bedgraph} norm stringent ${sample_id}_seachr_peaks_ctrl
SEACR_1.3.sh ${sample_bedgraph} 0.01 norm stringent ${sample_id}_seachr_peaks_top0.01
"""
}


/*
* 9.Peak stats
*/
process PEAK_STATS {
publishDir "${params.outdir}/Peak_calling/Peak_info", mode: 'copy'

when:
!params.skip_plots

input:
val(peak_info) from ch_design_peak_stats.join(ch_peak_ctrl).join(ch_peak_top001).collect()
val(bam) from ch_bam_peakstat.collect()
path align_rmdup_info from ch_align_redup_summary

output:
path "FRIP.txt" into ch_frip
path "Peak_reproducibility.txt" into ch_peak_reprod
path "Peak_plots.pdf" into ch_peak_plot



script:
"""
#!/usr/bin/env Rscript

packages <- c("dplyr", "stringr", "ggplot2", "viridis", "cowplot", "chromVAR", "GenomicRanges")
lapply(packages, library, character.only = TRUE)


# Read peak info
peak_info <-   print("${peak_info}")
peak_info <- str_replace_all(peak_info, "\\\\[|\\\\]", "")
peak_info <- gsub(', ', ',', peak_info)
peak_info<- read.table(text=peak_info,col.names=c('Sample','Control', 'Group','Replicate', 'Peak_path_ctrl', 'Peak_path_top001'), sep=",")

peak_stat=c()
peak_width=c()
for(sample in peak_info\$Sample){
      peaks_ctrl = read.table(toString(peak_info[peak_info\$Sample== sample,'Peak_path_ctrl']), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peaks_top001 = read.table(toString(peak_info[peak_info\$Sample== sample,'Peak_path_top001']), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peak_stat_ctrl = data.frame(Peak_id = paste0(sample, "_Control", sep=""),Sample = sample, Control = peak_info[peak_info\$Sample== sample,'Control'],  Group = peak_info[peak_info\$Sample== sample,'Group'],Replicate = peak_info[peak_info\$Sample== sample,'Replicate'], Peak_type = "Control", Number_of_peaks = nrow(peaks_ctrl), Peak_file = toString(peak_info[peak_info\$Sample== sample,'Peak_path_ctrl']))
      peak_stat_top001 = data.frame(Peak_id = paste0(sample, "_Top001", sep=""),Sample = sample, Control = peak_info[peak_info\$Sample== sample,'Control'],  Group = peak_info[peak_info\$Sample== sample,'Group'], Replicate = peak_info[peak_info\$Sample== sample,'Replicate'], Peak_type = "Top001", Number_of_peaks = nrow(peaks_top001), Peak_file =  toString(peak_info[peak_info\$Sample== sample,'Peak_path_top001']))
      peak_width = rbind(data.frame(Width = peaks_ctrl\$width, Peak_type = "Control", Sample = sample, Group=peak_info[peak_info\$Sample== sample,'Group'], Replicate=peak_info[peak_info\$Sample== sample,'Replicate']), data.frame(Width = peaks_top001\$width, Peak_type = "Top001", Sample = sample, Group=peak_info[peak_info\$Sample== sample,'Group'], Replicate=peak_info[peak_info\$Sample== sample,'Replicate']))  %>% rbind(peak_width, .)
      peak_stat <- rbind(peak_stat_ctrl,peak_stat_top001 ) %>% rbind(peak_stat, .)
}

# Peak reproducibility
groups = unique(peak_stat\$Group)
reps = unique(peak_stat\$Replicate)
peak_type = unique(peak_stat\$Peak_type)
peak_overlap = c()
for(type in peak_type){
  for(group in groups){
    overlap.gr = GRanges()
    for(rep in reps){
      peak_repr = read.table(peak_stat[peak_stat\$Group== group & peak_stat\$Peak_type==type & peak_stat\$Replicate==rep,'Peak_file'], header = FALSE, fill = TRUE)
      peak_repr.gr = GRanges(peak_repr\$V1, IRanges(start = peak_repr\$V2, end = peak_repr\$V3), strand = "*")
      if(length(overlap.gr) >0){
        overlap.gr = overlap.gr[findOverlaps(overlap.gr, peak_repr.gr)@from]
      }else{
        overlap.gr = peak_repr.gr

      }
    }
    peak_overlap = data.frame(Reproducible_peaks = length(overlap.gr), Group = group, Peak_type = type) %>% rbind(peak_overlap, .)
  }
}
peak_reprod = left_join(peak_stat, peak_overlap, by = c("Group", "Peak_type")) %>% mutate(Reproducible_rate = Reproducible_peaks/Number_of_peaks * 100)
write.table(peak_reprod[,c(1,2,3,4,5,6,7,9,10)], "Peak_reproducibility.txt", col.names = T, row.names = F, sep = "\t", quote = FALSE)

# FRIP
# Read peak info
bam <-   print("${bam}")
bam <- str_replace_all(bam, "\\\\[|\\\\]", "")
bam <- gsub(', ', ',', bam)
bam <- read.table(text=bam,col.names=c('Sample','Bam_file'), sep=",")


frip_info = c()
for(bam_sample in bam[bam\$Sample %in% peak_reprod\$Sample,"Sample"]){
    frip_peak = read.table(toString(peak_info[peak_reprod\$Sample== bam_sample & peak_reprod\$Peak_type=="Control",'Peak_path_ctrl']), header = FALSE, fill = TRUE)
    frip_peak.gr = GRanges(seqnames = frip_peak\$V1, IRanges(start = frip_peak\$V2, end = frip_peak\$V3), strand = "*")
    frip_fragment_counts <- getCounts(toString(bam[bam\$Sample== bam_sample,'Bam_file']), frip_peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
    frip_inpeak = counts(frip_fragment_counts)[,1] %>% sum
    frip_info = rbind(frip_info, data.frame(Counts_in_peak = frip_inpeak, Sample = bam_sample))
  }
aln_rmdup <- read.table("${align_rmdup_info}", header=T)
frip = left_join(frip_info, aln_rmdup, by = c("Sample" )) %>% mutate(FRIP = Counts_in_peak/Mapped_Fragments * 100)
write.table(frip, "FRIP.txt", col.names = T, row.names = F, sep = "\t", quote = FALSE)

# PLOTS
peak_reprod\$Replicate <- as.factor(peak_reprod\$Replicate)
fig_peak1 = peak_reprod %>% ggplot(aes(x = Group, y = Number_of_peaks, fill = Group)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  facet_grid(~Peak_type) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_discrete(labels = NULL,breaks=NULL)+
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")+
  ggtitle("A. Peak numbers")

peak_width\$Replicate <- as.factor(peak_width\$Replicate)
fig_peak2 = peak_width %>% ggplot(aes(x = Group, y = Width, fill = Group)) +
  geom_violin() +
  facet_grid(Replicate~Peak_type) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_y_continuous(trans = "log", breaks = c(400, 3000, 22000)) +
  scale_x_discrete(labels = NULL,breaks=NULL)+
  theme_bw(base_size = 18) +
  ylab("Width of Peaks") +
  xlab("")+
  ggtitle("B. Peak width")

fig_peak3 = peak_reprod %>% ggplot(aes(x = Group, y = Reproducible_rate, fill = Group, label = round(Reproducible_rate, 2))) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0.1) +
  facet_grid(Replicate~Peak_type) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_discrete(labels = NULL,breaks=NULL)+
  theme_bw(base_size = 18) +
  ylab("% of Peaks Reproduced") +
  xlab("")+
  ggtitle("C. Peak reproducibility")

frip\$Replicate <- as.factor(frip\$Replicate)
fig_peak4 = frip %>% ggplot(aes(x = Group, y = FRIP, fill = Group, label = round(FRIP, 2))) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15))+
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_x_discrete(labels = NULL,breaks=NULL)+
  theme_bw(base_size = 18) +
  ylab("% of Fragments in Peaks") +
  xlab("")+
  ggtitle("D. FRIP")

legend_peak <- get_legend( fig_peak1 + guides(color = guide_legend(nrow = 1)) +theme(legend.position = "bottom"))
plots_peaks <- plot_grid(
  fig_peak1+ theme(legend.position="none"),
  fig_peak2+ theme(legend.position="none"),
  fig_peak3+ theme(legend.position="none"),
  fig_peak4+ theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 2
)

pdf(file="Peak_plots.pdf", width = 12, height=8)
plot_grid(plots_peaks, legend_peak, ncol = 1, rel_heights = c(1, .1))
dev.off()
"""
}


/*
* 10. BogWig genreation
*/
process BIGWIG {
publishDir "${params.outdir}/BigWig", mode: 'copy'

when:
!params.skip_bigwig

input:
tuple val(sample_id), path(bam) from ch_bam_bw

output:
path "${sample_id}_bowtie2_rmDup_q${params.quality_score}.bw" into ch_bw


script:
"""
samtools sort -@ $task.cpus  -o ${sample_id}_bowtie2_rmDup_q${params.quality_score}_sorted.bam ${bam}
samtools index ${sample_id}_bowtie2_rmDup_q${params.quality_score}_sorted.bam
bamCoverage -b ${sample_id}_bowtie2_rmDup_q${params.quality_score}_sorted.bam --effectiveGenomeSize ${params.genome_size} --normalizeUsing RPGC --centerReads -o ${sample_id}_bowtie2_rmDup_q${params.quality_score}.bw
"""
}
