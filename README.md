**Script Description: 16S Amplicon Data Analysis using cutadapt and Dada2**

This script is designed for the analysis of 16S amplicon data using a combination of cutadapt, Dada2, and knitr. The script is tailored for processing raw sequencing data, denoising, and generating a comprehensive analysis report. The workflow includes:

    cutadapt:
        Trims adapter sequences and filters low-quality reads from raw sequencing data.

    Dada2:
        Performs denoising, dereplication, and generates a table of amplicon sequence variants (ASVs) to represent the microbial community.

    knitr:
        Integrates R code and text in a dynamic report, facilitating reproducibility and providing a detailed overview of the analysis process.

How to Run:

    Ensure that cutadapt, Dada2, and knitr are installed on your system.
    Customize input parameters within the script (e.g., file paths, quality thresholds).
    Execute the script in a terminal using the command: bash script_name.sh

Requirements:

    cutadapt
    Dada2 (R package)
    knitr (R package)

Output:

    Processed data files
    ASV table
    Quality assessment plots
    Comprehensive analysis report

This script streamlines the 16S amplicon data analysis pipeline, providing a robust and reproducible workflow for researchers studying microbial communities.





**Script Description:** **16S Amplicon Data Processing and Plotting usign 16S_Amplicon_Data_Processing_and_Plotting.Rmd**

This R script is tailored for the analysis of 16S amplicon data, focusing on data processing and the generation of informative plots. The script utilizes various R packages, including readODS, tidyverse, kableExtra, ggthemes, RColorBrewer, and vegan.

Key Steps:

    Load Libraries:
        Imports necessary R packages to facilitate data processing and plotting.

    Read Files:
        Reads data from an ODS file (Metadata_ASV_Taxa.ods) containing information about metadata, ASVs, and filtered ASV taxa.

    Normalization:
        Normalizes ASV counts for each peat soil sample, accounting for differences in sequencing depth.

    Transform to Long Format:
        Converts data from the "ASV_Taxa_Filter" sheet to long format, joining with metadata and summarizing.

    Group and Summarize:
        Groups data by taxonomic levels (Phylum, Class, Order, Family, Genus), summarizes relative abundance, and creates plots for each taxonomic level.

    Custom Plots:
        Produces customized plots, arranged and color-coded based on taxonomic information.

    ANOVA and Tukey HSD:
        Performs ANOVA and Tukey's Honest Significant Difference (HSD) test for a specific taxon (Genus "KS41") to explore differences among treatment groups.

How to Use:

    Ensure R is installed on your system along with required packages.
    Customize file paths or parameters if necessary.
    Run the script in an R environment or using Rscript in the terminal.

Outputs:

    Processed data files.
    Plots representing the relative abundance of taxa at various levels.
    Summary results from ANOVA and Tukey's HSD.

This script provides a comprehensive workflow for the analysis of 16S amplicon data, aiding researchers in understanding the microbial composition within different treatment groups.





**Here are the descriptions for the additional scripts:**

Script: Indicator Species Analysis (ISA)

This R script performs Indicator Species Analysis (ISA) on amplicon data. It reshapes and transposes the ASV taxa filter data, runs the multipatt() function, and saves the summary of the results in a text file named "Species_Perc_summary.txt." The ISA identifies ASVs that are indicative of specific site descriptions, aiding in the identification of microbial indicators for different environmental conditions.
Script: MOB β-diversity indices and Ordination Plot

This R script focuses on β-diversity analysis for microbial communities, specifically targeting MOB (Methane Oxidizing Bacteria). It reads data from an ODS file, performs normalization, executes Non-Metric Multidimensional Scaling (NMDS) ordination, and generates plots depicting the microbial community structure. The resulting plots are saved as "MOB_anaerobic_Ordination_Plot.svg" and "Anaerobic_MOB_Ordination_Plot.svg," offering insights into the differentiation of microbial communities across different treatment groups.
Script: β-diversity Indices and ANOSIM Analysis

This script calculates β-diversity indices and conducts statistical analyses such as ANOSIM (Analysis of Similarities) and Adonis (PERMANOVA). Utilizing the vegan and ggforce libraries, it generates NMDS plots and evaluates the dissimilarity between microbial communities based on Bray-Curtis distances. The results of ANOSIM and Adonis are printed, providing insights into the microbial community structure across different treatment groups.
Script: Microbial Community Analysis

This script analyzes the overall microbial community, reading metadata, ASV, and Taxa data from an ODS file. It performs normalization, transforms data to long format, and groups and summarizes the microbial community by Phylum and Genus. The top Phyla and Genera are visualized in bar plots, illustrating their relative abundances across different treatment groups. The resulting plots are saved as "All_Top_Phylum_Plot.svg" and "All_Top_Genera_Plot1.svg" or "All_Top_Genera_Plot22.svg," providing valuable insights into the dominant microbial taxa.
Script: All Microbial β-diversity Indices and Ordination Plot

This R script extends the analysis to all microbial taxa, calculating β-diversity indices and generating an NMDS ordination plot. It reads ASV and Taxa data, creates a phyloseq object, and performs NMDS ordination using the Bray-Curtis dissimilarity. The resulting plot, "All_anaerobic_Microbial3_Ordination_Plot.svg," visualizes the dissimilarity of microbial communities among different treatment groups. Additionally, "All_anaerobic_Microbial_Ordination2_Plot.svg" presents a detailed ordination plot with customized coloring, aiding in the interpretation of microbial community structures.
