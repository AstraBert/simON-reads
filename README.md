# simON-reads

## General purpose and application
simON-reads ("Simulate Oxford Nanopore Reads") is a simple yet powerful tool to generate fastq files containing MiniON-like long reads: this python script generates artificial DNA sequencing reads with specified variations and errors for given reference sequences. It utilizes the BioPython library for handling DNA sequences and provides options to introduce single nucleotide polymorphisms (SNPs) and sequencing errors.
simON-reads, thus, represent a flexible and customizable tool for generating artificial DNA sequencing data, making it valuable for testing and validating bioinformatics softwares and pipelines. The introduced variations and errors simulate real-world scenarios, allowing for thorough testing of downstream analysis pipelines.

## Installation

1. Clone the current GitHub repository by running:
   ```bash
   git init
   git clone https://github.com/AstraBert/simON-reads
   ```
2. Move to the cloned directory and run install.sh: it will set up a new conda environment where the two necessary dependencies, [biopython](https://biopython.org/) and [matplotlib](https://matplotlib.org/) will be stored:
   ```bash
   cd ./simON-reads
   bash ./scripts/install.sh
   ```
3. Test the installation by activating the environment:
   ```bash
   conda activate ./scripts/environments/simON-reads
   simON_reads.py -h
   conda deactivate
   ```
## Options and testing
The script comes with three viable option (one is required, the other two are optional)
```
simON_reads.py -i, --infile INFILE [-snp, --single_nucleotide_polymorphism SAMPLE:POS:REF>ALT,SAMPLE:POS:REF>ALT,...] [-n, --nreads READS_NUMBER]

-i or --infile: Path to the input FASTA file containing the reference sequence(s).

-snp or --single_nucleotide_polymorphism: Insert single nucleotide variants.
Insert a single nucleotide variant; the syntax of this option should be SAMPLE:POS:REF>ALT,SAMPLE:POS:REF>ALT,...,SAMPLE:POS:REF>ALT
(it should be separated by commas without blank spaces) where SAMPLE is the header of the sequence (withouth \">\") in the original
fasta file, POS is an integer that indicates the position (0-based) of the polymorphic site, REF is the reference allele,
ALT is the alternative allele you want to be put. This will generate a diploid-like distribution of variants. (Default is "NO_SNP")

-n or --nread: Number of reads to generate for each reference sequence (default is 2000).

Input simON_reads.py -h,--help to show the help message
```

You will find a test sample of reference sequences in the test folder; to try the script, you can run:
```bash
cd ./test
simON_reads.py -i reference.fasta -snp 28S-rRNA:3:G>T,28S-rRNA:41:C>G,S7:0:A>G,S7:63:C>T -n 1000 > test.fastq
```
Always remember to redirect the stream to your desidered file, unless you want it to be printed on the standard output of your terminal.

## Output

The script generates artificial DNA sequencing reads based on the provided parameters. It outputs the average quality distribution of the generated reads as a histogram and saves it as `avg_quality_distribution.png` in the current directory.

## How does it work?

### Core Functions:

#### `seqs_to_file(genomes_dict, snp_string, nreads)`

This function takes three parameters:

- `genomes_dict`: A dictionary containing reference sequences, where keys are sequence headers, and values are corresponding sequences.
- `snp_string`: A string specifying single nucleotide polymorphisms (SNPs) in the format `SAMPLE:POS:REF>ALT,SAMPLE:POS:REF>ALT,...`.
- `nreads`: The number of artificial reads to generate for each reference sequence.

The function iterates over the reference sequences and generates artificial reads with variations:

1. **Reverse Complement and Chimeric Reads:**
   - It creates a list (`revcomp`) with 5% of the reads being reverse complemented.
   - It creates another list (`chimers`) with 0.5% of the reads being chimeric (merged from two sequences).

2. **Generate Reads:**
   - For each read:
     - If the read is marked for reverse complementation, it generates a reverse complement of the reference sequence.
     - If the read is marked as chimeric, it merges two random reference sequences.
     - Otherwise, it generates a normal read with potential SNPs and random sequencing errors.

3. **Quality Calculation:**
   - The quality of each read is calculated using the `quality_string` function, and the average quality is recorded.

4. **Output:**
   - The function prints the generated reads and their qualities.

#### `quality_string(seq)`

This function generates a random quality string for a given DNA sequence. It assigns ASCII phred score codes (from 1 to 30) to each nucleotide in the sequence.

### Main Execution:

The main execution part of the script utilizes the functions mentioned above:

1. **Argument Parsing:**
   - It uses `ArgumentParser` to parse command-line arguments, including the input FASTA file (`-i`), SNP string (`-snp`), and the number of reads to generate (`-n`).

2. **Loading Reference Sequences:**
   - The script loads reference sequences from the specified input FASTA file using the `load_data` function.

3. **Generating artificial Reads:**
   - It calls the `seqs_to_file` function to generate artificial reads based on the reference sequences, SNPs, and the specified number of reads.

4. **Plotting and Saving:**
   - The script plots the average read quality distribution as a histogram using `matplotlib` and saves it as `avg_quality_distribution.png` in the current directory.

5. **Execution Duration:**
   - The script prints the duration of execution.

### Overall Workflow:

1. **Loading Data:**
   - The reference sequences are loaded from the input FASTA file.

2. **Generating Reads:**
   - artificial reads are generated for each reference sequence.
   - SNPs and sequencing errors are introduced.
   - Reverse complementation and chimeric reads are handled.

3. **Quality Assessment:**
   - The average quality distribution of the generated reads is calculated.

4. **Output:**
   - The generated reads and their qualities are printed.

5. **Plotting:**
   - The script creates a histogram showing the distribution of average read qualities.

6. **Duration and File Saving:**
   - The duration of execution is printed.
   - The histogram is saved as an image (`avg_quality_distribution.png`).


## License and rights of usage
Please note that simON-reads is still experimental and may contain errors or may output not-100%-reliable results, so always check them and pull issues whenever you feel it to be the case, we'll be on your back as soon as possible to fix/implement/enhance whatever you suggest!

The code is protected by the GNU v.3 license. As the license provider reports: "Permissions of this strong copyleft license are conditioned on making available complete source code of licensed works and modifications, which include larger works using a licensed work, under the same license. Copyright and license notices must be preserved. Contributors provide an express grant of patent rights".

If you are using simON-reads for you project, please consider to cite the author of this code (Astra Bertelli) and this GitHub repository.
