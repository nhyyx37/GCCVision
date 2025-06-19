# GCCVision: Genome Contribution Calculator and Visualizer

**Version: 1.0**

GCCVision is an integrated toolkit for calculating and visualizing parental genome contributions in breeding populations. It consists of two main components:
1.  A Python script (`GCCVision.py`) for processing VCF files to calculate genome contributions.
2.  An interactive HTML-based tool (`GCCVision.html`) for visualizing genotype distributions along chromosomes.

## Features

-   **Informative Site Identification**: Automatically finds informative homozygous sites from a VCF file based on two parents.
-   **Genotype Classification**: Classifies sites in progeny as homozygous from parent A (`Hom_A`), homozygous from parent B (`Hom_B`), or heterozygous (`Het_AB`).
-   **Noise Reduction**: Optional sliding window filtering to clean up potential genotyping errors.
-   **Contribution Calculation**: Calculates genome-wide and per-chromosome contribution percentages for each parent.
-   **Interactive Visualization**: A feature-rich, browser-based interface to visualize genotype maps.
-   **High-Quality Export**: Export visualizations to PNG, JPEG, or SVG formats with customizable resolution.
-   **No Installation Required**: The visualizer is a single HTML file that runs in any modern web browser like Chrome, Edge or Firefox.

---

## Part 1: `GCCVision.py` - Genome Contribution Calculator

This command-line tool processes multi-sample VCF files to identify informative sites and calculate parental genome contributions for specified samples.

### Usage

The script can be run in two modes:

1.  **Standard Mode** (from a VCF file):
    ```bash
    python GCCVision.py -v input.vcf.gz -a Parent1 -b Parent2 -s Sample1,Sample2 -o results
    ```
2.  **Using a Pre-existing Informative Sites File**:
    ```bash
    python GCCVision.py -i informative_sites.tsv -v input.vcf -s Sample1,Sample2 -o results
    ```

### Parameters

| Parameter                 | Short | Description                                                                                                        | Required | Default |
| ------------------------- | ----- | ------------------------------------------------------------------------------------------------------------------ | -------- | ------- |
| `--vcf`                   | `-v`  | Input multi-sample VCF file. Can be gzipped (`.gz`).                                                               | Yes      | -       |
| `--informative-sites`     | `-i`  | Path to a pre-existing informative sites file. If provided, the site identification step is skipped.               | No       | -       |
| `--parent-a`              | `-a`  | Sample name of Parent A in the VCF file.                                                                           | Yes*     | -       |
| `--parent-b`              | `-b`  | Sample name of Parent B in the VCF file.                                                                           | Yes*     | -       |
| `--samples`               | `-s`  | Progeny samples to analyze. Can be a comma-separated list (e.g., "Sample1,Sample2") or a path to a `.txt` file with one sample per line. | Yes      | -       |
| `--output-prefix`         | `-o`  | Prefix for all output files.                                                                                       | Yes      | -       |
| `--threads`               | `-t`  | Number of threads for parallel processing.                                                                         | No       | 1       |
| `--filter-times`          | `-f`  | Number of filtering rounds using the sliding window method. Set to 0 to disable.                                   | No       | 0       |
| `--filter-window-size`    | `-w`  | The size of the sliding window for filtering.                                                                      | No       | 5       |

_*`--parent-a` and `--parent-b` are required unless `--informative-sites` is used._

### Output Files

The script generates the following files, prefixed with the string provided via `-o`:

1.  **`{prefix}_informative_sites.tsv`**: Lists all sites where parents are homozygous and different.
2.  **`{prefix}_summary.tsv`**: A summary of contribution rates for each sample, both genome-wide and per-chromosome.
3.  **`{prefix}_{sample}_details.tsv`**: For each sample, this file contains the classified genotype (`Hom_A`, `Hom_B`, `Het_AB`) for each informative site. **This file is used as the "Site Data File" in `GCCVision.html`**.
4.  **`{prefix}_{sample}_details_filtered_N.tsv`**: Filtered version of the details file, if `-f` is greater than 0.

---

## Part 2: `GCCVision.html` - Interactive Visualizer

This tool allows you to create detailed and customizable visualizations of genotype maps. To use it, simply open the `GCCVision.html` file in a modern web browser like Chrome, Edge or Firefox.

### Input Data Formats

The visualizer requires specific tab-separated formats for its input files. You can click the **"Download Sample"** button for each file type to see an example.

1.  **Chromosome Data File (`*Required`)**: A file listing chromosome names and their total lengths.
    *   **Format**: `ChromosomeName\tLength`
    *   **Example**:
        ```
        Chr01	56831624
        Chr02	48577505
        ```

2.  **Site Data File (`Optional`)**: A file listing the positions and types of sites. The `{prefix}_{sample}_details.tsv` or `{prefix}_{sample}_details_filtered_N.tsv` files generated by `GCCVision.py` can be used directly for this.
    *   **Format**: `ChromosomeName\tPosition\tSiteType`
    *   **Example**:
        ```
        Chr01	1959	Hom_A
        Chr01	2528	Het_AB
        Chr01	112623	Hom_B
        ```

3.  **Gene Data File (`Optional`)**: A file listing genes or markers to be displayed on the chromosomes.
    *   **Format**: `GeneName\tChromosomeName\tStartPosition\tEndPosition\tColor`
    *   **Example**:
        ```
        Satt300	Chr01	28572124	28572367	red
        Satt429	Chr01	47217842	47218105	blue
        ```

### How to Use the Visualizer

1.  **Load Data**:
    -   Use the **"Choose File"** buttons to select your local data files.
    -   Click **"Load Uploaded Data"** to render the visualization.
    -   Alternatively, click **"Load Sample Data"** to see a pre-loaded example.

2.  **Customize Display Settings**:
    -   **Site Filtering**: Toggle the visibility of `Hom_A`, `Hom_B`, and `Het_AB` sites. You can also toggle the genetic background analysis view.
    -   **Legend Settings**: Change the text labels for different elements in the legend.
    -   **Display Parameters**: Adjust a wide range of visual elements:
        -   **Histogram Window Size**: Change the bin size for site aggregation.
        -   **Chromosome Dimensions**: Control the width, height, and spacing of chromosomes.
        -   **Font Sizes**: Adjust the font size for gene names and background ratios.
        -   **Colors & Styles**: Set colors for chromosomes, genes, and border widths.

3.  **Export the Visualization**:
    -   In the **Export Settings** panel, choose a format (`PNG`, `JPEG`, `SVG`).
    -   Adjust the quality and scale factor for the output image.
    -   Click **"Export Image"** to download the final visualization.

## Recommended Workflow

1.  Run `GCCVision.py` with your VCF data to generate the analysis files.
    ```bash
    python GCCVision.py -v my_data.vcf.gz -a P1 -b P2 -s Progeny_1 -o Progeny_1_analysis -f 2 -w 10
    ```
2.  Prepare a **Chromosome Data File** that lists the names and lengths of the chromosomes in your reference genome.
3.  Open `GCCVision.html` in your web browser.
4.  Load your prepared **Chromosome Data File** and the **`Progeny_1_analysis_Progeny_1_details_filtered_2.tsv`** file (as the Site Data File).
5.  Customize the visualization to your liking.
6.  Export the final image for your publication or report.

