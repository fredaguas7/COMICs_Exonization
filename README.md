# COMICs_Exonization
Pipeline For TE Exonization Analysis

# *WORK IN PROGRESS*

# Step 0 - Alignment 

Start with STAR alignments. We will use the alignments that we do for the TE counts but we can discuss if it's worth it to use other parameters.

# Step 1 - Transcriptome assembly with StringTie 

Prepare environment

```
conda create -n string_tie
conda activate string_tie
conda install bioconda::stringtie
```

Run StringTie - This will assemble the transcripts, taking a refference in account - we will then look for novel transcripts that do not match the refference

```
stringtie -o output.gtf -G my.annotation.gtf my.bam
```

Loop example

```
for bam_file in Unique_align/*.bam; do
    
    # 1. Extract the filename (e.g., SAMN11569826_Aligned.sortedByCoord.out.bam)
    filename=$(basename "$bam_file")
    
    # 2. Extract the Sample ID by removing everything after the first underscore
    # This turns "SAMN11569826_Aligned..." into "SAMN11569826"
    sample_id="${filename%%_*}"

    echo "Processing sample: $sample_id"

    # 3. Run StringTie
    stringtie -o "StringTie/${sample_id}_assembly.gtf" -p 8 \
        -G /mnt/comics-data/COMICSlab/genomes/human/GRCh38.primary_assembly/GeneAnnotations/gencode.v44.annotation.gtf \
        "$bam_file"
        
done
```

# Step 2 - Compare with the refference

Prepare environment
```
conda install bioconda::gffcompare
```

Run GFFCompare

```
gffcompare -R -r my.annotation.gtf -o output.gtf -V prefix_for_other_outputs
```

-r: Annotation file
-R: Ignore transcripts that do not overlap with any other in annotation
-o: Output name
-V : Verbose

Loop example

```
for gtf_file in *_assembly.gtf; do

    # 1. Extract the filename (e.g., SAMN11569826_assembly.gtf)
    filename=$(basename "$gtf_file")

    # 2. Extract the Sample ID by removing the suffix "_assembly.gtf"
    sample_id="${filename%_assembly.gtf}"

    echo "Comparing sample: $sample_id"

    # 3. Run GffCompare
    # -o sets the prefix for the output files (e.g., GffCompare/SAMN11569826)
    gffcompare -R -r /mnt/comics-data/COMICSlab/genomes/human/GRCh38.primary_assembly/GeneAnnotations/gencode.v44.annotation.gtf \
               -o "GffCompare/${sample_id}" -V \
               "$gtf_file"

done
```

# Step 3 - Filter transcripts overlaping with TEs

```
# 1. Set the path to your TE BED file here
TE_BED="/mnt/comics-data/f.agua/Genomes/Human_TEIndividual.bed"

# 2. Create output directory
mkdir -p TE_exon_overlap

# 3. Loop through the assemblies
for gtf_file in GffCompare/*.annotated.gtf; do

    # Extract ID
    filename=$(basename "$gtf_file")
    sample_id="${filename%.annotated.gtf}"

    echo "Intersecting TEs for: $sample_id"

    # 4. Run the intersection
    # Step A: awk '$3 == "exon"' -> Filters the GTF to keep only exon lines
    # Step B: bedtools intersect -> Finds overlaps between exons (stdin) and TEs (-b)
    # -wo -> Reports the full entry from both files PLUS the number of overlapping bases
    
    awk '$3 == "exon"' "$gtf_file" | \
    bedtools intersect -a stdin -b "$TE_BED" -wo -F 0.8 \
    > "TE_exon_overlap/${sample_id}_TE_exon_overlap.bed"

done

###### That gives a bed file for the exons. Now Filter the GTF to keep only transcripts with exons that overlap TEs

# Create a directory for the final filtered GTFs
mkdir -p TE_Filtered_GTFs

# Loop through the TE overlap files
for overlap_file in TE_exon_overlap/*_TE_exon_overlap.bed; do

    # 1. Get the Sample ID
    # Removes the directory path and the suffix
    filename=$(basename "$overlap_file")
    sample_id="${filename%_TE_exon_overlap.bed}"

    echo "Filtering GTF for sample: $sample_id"

    # Define file paths
    # The annotated GTF created by gffcompare
    input_gtf="GffCompare/${sample_id}.annotated.gtf"
    # The temporary list of IDs to keep
    id_list="TE_Filtered_GTFs/${sample_id}_ids.tmp"
    # The final output
    output_gtf="TE_Filtered_GTFs/${sample_id}_TE_containing.gtf"

    # 2. Extract Transcript IDs from the overlap file
    # We grep strictly for 'transcript_id "STRG.xxxx";' to ensure exact matching later
    grep -o 'transcript_id "[^"]*";' "$overlap_file" | sort | uniq > "$id_list"

    # 3. Filter the GffCompare GTF
    # -F: Interpret patterns as fixed strings (faster)
    # -f: Read patterns from the file (our ID list)
    grep -F -f "$id_list" "$input_gtf" > "$output_gtf"

    # 4. Clean up the temporary ID file
    rm "$id_list"

done
```

# Step 3.5 - Meta summary for all samples
Short python script

```
import pandas as pd
import glob
import os
import re

# --- CONFIGURATION ---
OVERLAP_DIR = "TE_exon_overlap"            # Bedtools output
ORIGINAL_GTF_DIR = ""         # Original GTFs (with TPMs)
FILTERED_GTF_DIR = "TE_Filtered_GTFs"  # Optional: For Class Codes
OUTPUT_FILE = "Meta_TE_Events.csv"

# Regex to find Transcript ID in the attributes string
id_pattern = re.compile(r'transcript_id "([^"]+)"')
code_pattern = re.compile(r'class_code "([^"]+)"')

def load_te_overlaps(overlap_file):
    """
    Returns a dictionary: { 'STRG.1.1': 'L1PA2 (chr1:500-600)' }
    If a transcript hits multiple TEs, we join them with '; '
    """
    tid_to_te = {}
    try:
        with open(overlap_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 16: continue

                # 1. Get Transcript ID (Col 8)
                attribs = parts[8]
                match = id_pattern.search(attribs)
                if match:
                    tid = match.group(1)
                    
                    # 2. Get TE Info (Bedtools cols 9+)
                    # BED cols: 9=Chr, 10=Start, 11=End, 12=Name
                    te_loc = f"{parts[12]} ({parts[9]}:{parts[10]}-{parts[11]})"
                    
                    if tid in tid_to_te:
                        # Avoid duplicates if multiple exons hit the SAME TE
                        if te_loc not in tid_to_te[tid]:
                            tid_to_te[tid].append(te_loc)
                    else:
                        tid_to_te[tid] = [te_loc]
    except Exception as e:
        print(f"Error reading overlap file: {e}")
        
    # Join lists into strings
    return {k: "; ".join(v) for k, v in tid_to_te.items()}

def load_class_codes(filtered_file):
    """Returns { 'STRG.1.1': 'j' }"""
    tid_to_code = {}
    if os.path.exists(filtered_file):
        with open(filtered_file, 'r') as f:
            for line in f:
                if "\ttranscript\t" in line:
                    match = id_pattern.search(line)
                    code = code_pattern.search(line)
                    if match and code:
                        tid_to_code[match.group(1)] = code.group(1)
    return tid_to_code

# --- NEW MAIN EXECUTION ---
all_rows = []
# We add "TPM_Val" as a specific column so you can sort Exons by expression later
gtf_columns = ["Sample", "Seqname", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes", "TE_Info", "Class_Code", "TPM_Val"]

print("Starting Meta-File generation...")
overlap_files = glob.glob(os.path.join(OVERLAP_DIR, "*_TE_exon_overlap.bed"))

for o_file in overlap_files:
    # 1. Get Sample ID
    fname = os.path.basename(o_file)
    sample_id = fname.replace("_TE_exon_overlap.bed", "")
    print(f"Processing: {sample_id}")
    
    # 2. Load Maps
    # te_map: { 'STRG.1.1': 'L1PA2 (chr1:500-600)' }
    te_map = load_te_overlaps(o_file)
    
    # class_map: { 'STRG.1.1': 'j' }
    class_map = load_class_codes(os.path.join(FILTERED_GTF_DIR, f"{sample_id}_TE_containing.gtf"))
    
    # 3. Read Original GTF
    original_gtf = os.path.join(ORIGINAL_GTF_DIR, f"{sample_id}_assembly.gtf")
    
    # Helper dict to remember TPMs for exons
    # { 'STRG.1.1': 5.45 }
    local_tpm_cache = {} 

    if os.path.exists(original_gtf):
        # FIRST PASS: Quick scan to grab TPMs for our candidates
        # (Because exons often appear after the transcript line, or sometimes before, 
        #  it's safer to pre-load the TPMs so every exon gets a value).
        with open(original_gtf, 'r') as f:
            for line in f:
                if "\ttranscript\t" in line:
                    match = id_pattern.search(line)
                    tpm_match = re.search(r'TPM "([^"]+)"', line)
                    if match and tpm_match:
                        local_tpm_cache[match.group(1)] = tpm_match.group(1)

        # SECOND PASS: Actually build the rows
        with open(original_gtf, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                
                # Safety: Ensure line has columns (skip headers)
                if len(parts) < 9: continue

                # Extract Transcript ID from this line (works for exons AND transcripts)
                attribs = parts[8]
                match = id_pattern.search(attribs)
                
                if match:
                    tid = match.group(1)
                    
                    # --- THE CORE CHANGE ---
                    # Check if this ID is in our TE Whitelist
                    if tid in te_map:
                        te_info = te_map[tid]
                        class_code = class_map.get(tid, ".")
                        
                        # Grab the TPM we cached earlier (so Exons get a TPM value too!)
                        tpm_val = local_tpm_cache.get(tid, "0")

                        # Construct the row
                        row = [
                            sample_id,
                            parts[0], # Seqname
                            parts[1], # Source
                            parts[2], # Feature (transcript OR exon)
                            parts[3], # Start
                            parts[4], # End
                            parts[5], # Score
                            parts[6], # Strand
                            parts[7], # Frame
                            parts[8], # Attributes
                            te_info,  # The TE Name
                            class_code,
                            tpm_val   # Explicit TPM column
                        ]
                        all_rows.append(row)

# --- SAVE ---
if all_rows:
    df = pd.DataFrame(all_rows, columns=gtf_columns)
    
    # Optional: Parse TPM out of Attributes into its own column for easier sorting?
    # Uncomment next 2 lines if you want a dedicated TPM column
    df['TPM'] = df['Attributes'].str.extract(r'TPM "([^"]+)"').astype(float)
    # df = df.sort_values(by=['Sample', 'TPM'], ascending=[True, False])
    
    df.to_csv(OUTPUT_FILE, index=False)
    print("-" * 30)
    print(f"Success! Created {OUTPUT_FILE}")
    print(f"Total TE-containing transcripts found: {len(df)}")
else:
    print("No matching transcripts found. Check paths.")
```

Then just

```
python meta_summary.py
```

# Step 4 - Look for ORFs on the transcripts with TransDecoder2

Prepare environment

```
conda create -n TD2 python=3.11
conda activate TD2
pip install TD2
```

Prepare the string tie outputs - CONSIDER ALREADY THE FILTERED FOR EXONIZATION
Note: This util files don't actually come when you pip install. Download them directly from the github and when you run you have to give the path to the .pl file.
```
util/gtf_genome_to_cdna_fasta.pl transcripts.gtf my.genome.fasta > transcripts.fasta 

util/gtf_to_alignment_gff3.pl transcripts.gtf > transcripts.gff3
```

Example loop

```
STR_DIR="/mnt/comics-data/f.agua/PDAC/PublicData/GSE130688/StringTie/TE_Filtered_GTFs"
GENOME="/mnt/comics-data/COMICSlab/genomes/human/GRCh38.primary_assembly/Genome/GRCh38.primary_assembly.genome.fa"
TD_UTIL_DIR="/mnt/comics-data/f.agua/tools/TD2/TD2/util" # <--- UPDATE THIS to your TransDecoder util path

# Create Output Folders
mkdir -p TD2/Transcripts_fasta
mkdir -p TD2/Transcripts_gff

echo "Starting TransDecoder preparation..."

# Loop through the StringTie assemblies
for gtf in "$STR_DIR"/*_TE_containing.gtf; do
    
    # Extract Sample ID (e.g., SAMN11569826)
    filename=$(basename "$gtf")
    sample_id="${filename%_TE_containing.gtf}"
    
    echo "------------------------------------------------"
    echo "Processing Sample: $sample_id"

    # 1. Generate cDNA FASTA
    # This uses the exons in the GTF to extract sequences from the genome
    "$TD_UTIL_DIR/gtf_genome_to_cdna_fasta.pl" \
        "$gtf" "$GENOME" > "TD2/Transcripts_fasta/${sample_id}.fasta"

    # 2. Generate Alignment GFF3 
    # This is needed later to map ORF coordinates back to the genome
    "$TD_UTIL_DIR/gtf_to_alignment_gff3.pl" \
        "$gtf" > "TD2/Transcripts_gff/${sample_id}.gff3"

done

echo "Done! Check the TD2/ directory for your files."
```

Run TD2

```
TD2.LongOrfs -t transcripts.fasta -O OUTPUT_DIR  -m 40 --complete-orfs-only -S
```

Example loop

```
FASTA_DIR="TD2/Transcripts_fasta"
OUT_BASE="TD2/LongOrfs_Output"

mkdir -p "$OUT_BASE"

echo "Starting TransDecoder.LongOrfs loop..."

for fasta in "$FASTA_DIR"/*.fasta; do
    
    # 1. Get Sample ID
    filename=$(basename "$fasta")
    sample_id="${filename%.fasta}"
    
    # 2. Define specific output directory for this sample
    # (Important so temp files from different samples don't clash)
    sample_out_dir="$OUT_BASE/$sample_id"
    mkdir -p "$sample_out_dir"

    echo "Processing LongOrfs for: $sample_id"

    # 3. Run TransDecoder.LongOrfs
    # -t: Input Fasta
    # -O: Output Directory
    # -m 40: Minimum protein length
    # -S: Strand-specific (good for StringTie data)
    TD2.LongOrfs \
        -t "$fasta" \
        -O "$sample_out_dir" \
        -m 40 \
        --complete-orfs-only \
        -S

done

echo "All samples processed."
```

# Step 5 - Run Blastp to compare the ORFs with annotated proteins

Prepare environment 

```
conda create -n blast
conda activate blast
conda install biocore::blast-plus
```

Get proteome fasta from UniProt and create the database 

```
makeblastdb -in uniprotkb_proteome_UP000005640_2026_01_28.fasta -parse_seqids -out uniprot_human -title "UniProt_Human" -dbtype prot
```

Keep the database path

Run Blastp

```
blastp -query peptide_file -db database_path -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads N_threads > output.file
```

-outfmt 6: Standard tabular format needed by TransDecoder
-max_target_seqs 1: We only need the best hit
-evalue 1e-5: Standard cutoff for significance

Example loop

```
LONGEST_ORFS_DIR="TD2/LongOrfs_Output"

# !!! UPDATE THIS PATH !!!
# Path to your protein database (e.g., SwissProt)
DB_PATH="/mnt/comics-data/f.agua/Genomes/Human_proteome/uniprot_human"

# Number of threads to use (BLAST is CPU intensive)
THREADS=64

echo "Starting BLASTP loop..."

for sample_dir in "$LONGEST_ORFS_DIR"/*; do
    
    # Check if it's a directory
    if [ ! -d "$sample_dir" ]; then continue; fi

    # Get Sample ID for logging
    sample_id=$(basename "$sample_dir")
    
    # Input: The peptide file generated by LongOrfs
    pep_file="${sample_dir}/longest_orfs.pep"
    
    # Output: The BLAST result file
    blast_out="${sample_dir}/blastp.outfmt6"

    echo "Processing: $sample_id"

    if [ -f "$pep_file" ]; then
        
        # Run BLASTP
        # -outfmt 6: Standard tabular format needed by TransDecoder
        # -max_target_seqs 1: We only need the best hit
        # -evalue 1e-5: Standard cutoff for significance
        blastp -query "$pep_file" \
            -db "$DB_PATH" \
            -max_target_seqs 1 \
            -outfmt 6 \
            -evalue 1e-5 \
            -num_threads "$THREADS" \
            > "$blast_out"
            
        echo "  > BLAST complete. Saved to $blast_out"
    else
        echo "  > WARNING: No pep file found in $sample_dir"
    fi

done
```

# Step 6 - Filter only the ones with blastp hits (using TD2)

```
TD2.Predict -t transcript_fasta -O my.output --retain-blastp_hits blast.output.fmt6 --complete-orfs-only
```
-t: The original transcript sequences
-O: The directory where LongOrfs ran (must be the same!)
--retain-blastp_hits: Rescue valid proteins found by BLAST
--complete-orfs-only --retain-long-orfs-mode RETAIN_LONG_ORFS_MODE
    
Loop example

```
BASE_DIR="/mnt/comics-data/f.agua/PDAC/PublicData/GSE130688"
TRANSCRIPTS_DIR="$BASE_DIR/TD2/Transcripts_fasta"  # Where the .fasta files are
LONGEST_ORFS_DIR="$BASE_DIR/TD2/LongOrfs_Output"   # Where the intermediate folders are

echo "Starting TransDecoder.Predict..."

for sample_dir in "$LONGEST_ORFS_DIR"/*; do
    
    if [ ! -d "$sample_dir" ]; then continue; fi

    # 1. Identify Sample ID and Inputs
    sample_id=$(basename "$sample_dir")
    original_fasta="${TRANSCRIPTS_DIR}/${sample_id}.fasta"
    blast_out="${sample_dir}/blastp.outfmt6"
    
    echo "------------------------------------------------"
    echo "Processing: $sample_id"

    # 2. Check for required files
    if [ ! -f "$original_fasta" ]; then
        echo "  Error: Original FASTA not found at $original_fasta"
        continue
    fi
    
    if [ ! -f "$blast_out" ]; then
        echo "  Warning: BLAST output not found. Running Predict without homology evidence."
        # You might want to skip this 'else' if you insist on BLAST results
    fi

    # 3. Run TransDecoder.Predict
    # -t: The original transcript sequences
    # -O: The directory where LongOrfs ran (must be the same!)
    # --retain-blastp_hits: Rescue valid proteins found by BLAST
    # --complete-orfs-only --retain-long-orfs-mode RETAIN_LONG_ORFS_MODE
    
    TD2.Predict \
        -t "$original_fasta" \
        -O "$sample_dir" \
        --retain-blastp_hits "$blast_out" \
        --complete-orfs-only
    
    echo "  > Prediction complete."

done

echo "All samples processed."
```

This gives a gtf with relative coordinates. We want to convert back to real genome coordinates. TD2 also has a tool for this

```
BASE_DIR="/mnt/comics-data/f.agua/PDAC/PublicData/GSE130688"
TD_UTIL_DIR="/mnt/comics-data/f.agua/tools/TD2/TD2/util"

# Inputs
PREDICT_DIR="$BASE_DIR/TD2/Predict_Output"
ALIGN_DIR="$BASE_DIR/TD2/Transcripts_gff"
FASTA_DIR="$BASE_DIR/TD2/Transcripts_fasta"

echo "Mapping TransDecoder coordinates back to the Genome..."

for sample_dir in "$PREDICT_DIR"/*; do
    [ -d "$sample_dir" ] || continue
    
    sample_id=$(basename "$sample_dir")
    
    # Define Inputs
    # 1. The relative GFF3 from TransDecoder
    td2_gff="${sample_dir}/${sample_id}.fasta.TD2.gff3"
    
    # 2. The mapping GFF3 (Transcript -> Genome)
    align_gff="${ALIGN_DIR}/${sample_id}.gff3"
    
    # 3. The Sequence file
    fasta="${FASTA_DIR}/${sample_id}.fasta"
    
    # Output file
    genome_gff="${sample_dir}/${sample_id}.genome.gff3"

    if [[ -f "$td2_gff" && -f "$align_gff" ]]; then
        
        "$TD_UTIL_DIR/cdna_alignment_orf_to_genome_orf.pl" \
            "$td2_gff" \
            "$align_gff" \
            "$fasta" \
            > "$genome_gff"
            
        echo "  > Created genomic map: $genome_gff"
    else
        echo "  > Missing files for $sample_id"
    fi
done
```
