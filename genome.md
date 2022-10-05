This is the pipeline used to assembly, remove contaminants and annotate
the *Macrocystis pyrifera* from Tasmania (Iha et al. in prep). This
pipeline may be discontinuous, but the important commands are presented
here.

## Workflow

<div id="htmlwidget-685675ac01a2bb4418cc" style="width:1344px;height:1152px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-685675ac01a2bb4418cc">{"x":{"diagram":"\n        digraph {\n\n        \n        subgraph cluster_0 {\n           rankdir = TB\n           node [shape = box]\n           A [label = \"Genomic raw reads\"]\n           B [label = \"Genomic filtered reads\"]\n           C [label = \"Assembled genome&#92;n4713 scaffolds\"]\n           \n           #edges\n           A->B [label = \"Trimmomatic\"]\n           B->C [label = \"MaSurCa\"]\n           A->C [label = \"Bowtie2&#92;nDepth of coverage\"]\n           \n           label = \"Genome assembly\"\n           labeljust = \"l\"\n        }\n          \n        subgraph cluster_1 {\n           node [shape = box]\n           TE [label = \"Taxonomic designation: Eukaryote\"]\n           TB [label = \"Taxonomic designation: Bacteria\"]\n           TO [label = \"No taxonomic designation\"]\n           FS [label = \"Filtered scaffolds\"]\n           PC [label = \"Putative contaminant\"]\n           MT [label = \"Mapped transcripts \\n with evidence of introns?\"]\n           OT [label = \"Outliers?\"]\n           ST [label = \"Stramenopiles?\"]\n           \n           node [shape=none, width=0, height=0, label=\"\"]\n           p1 -> TE\n           p1 -> TB\n           p1 -> TO\n           TE -> ST\n           ST -> PC [label = \"No\"]\n           ST -> FS [label = \"Yes\"]\n           TB->MT\n           TO->MT\n           MT->FS [label = \"Yes\"]\n           MT->PC [label = \"No\"]\n           TB->OT\n           TO->OT\n           OT->PC [label = \"Yes\"]\n           OT->FS [label = \"No\"]\n\n           edge [dir=none]\n           C->p1 [label = \"BlobTools based on BLASTX (e-value < 10-20) search \\n against Uniprot, Silva and Custom database\"]\n\n          label = \"Removing scaffolds of putative contaminants\"\n          labeljust = \"l\"\n        }\n        \n        subgraph cluster_2 {\n           rankdir = TB\n           node [shape = box]\n           T [label = \"SRA Transcriptomes&#92;nPRJNA322132 and PRJNA353611\"]\n           DA [label = \"de novo assembly\"]\n           GG [label = \"genome-guided assembly\"]\n           AT [label = \"all transcriptomes\"]\n           TC [label = \"clean transcriptomes\"]\n           \n           #edges\n           T->GG [label = \"Trinity\"]\n           T->DA [label = \"Trinity\"]\n           GG -> AT\n           DA -> AT\n           AT -> TC [label = \"SeqClean\"]\n           DA -> MT\n           \n          label = \"Transcriptome\"\n          labeljust = \"l\"\n        }\n        \n        subgraph cluster_3 {\n        rankdir = TB\n           node [shape = box]\n           RM [label = \"Repeat Masker/Repeat Modeler\"]\n           GM [label = \"GeneMark-ES\"]\n           SNAP [label = \"SNAP\"]\n           AG [label = \"AUGUSTUS\"]\n           EVM [label = \"Evidence Modeler\"]\n           \n           PS [label=\"PASA\"]\n           B2D [label=\"BLASTX\"]\n           TPSI [label=\"Transposon-PSI\"]\n           CD [label=\"CD-HITS\"]\n           GGenes [label=\"Golden Genes\"]\n           \n           FS->RM\n           FS->PS\n           TC->PS\n           PS->B2D\n           B2D->TPSI\n           TPSI->CD\n           CD->GGenes\n           RM->GM\n           RM->SNAP\n           RM->AG\n           GM->EVM\n           SNAP->EVM\n           AG->EVM\n           GGenes->SNAP\n           \n           label = \"Genome Annotation\"\n           labeljust = \"l\"\n        }\n        \n        \n        }\n        ","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script>

## Assembly genome with MaSurCa

#### Make config file

``` bash
# quick run configuration file
DATA
PE = pe 350 50 Mp_R1.fastq.gz Mp_R2.fastq.gz
END
PARAMETERS
EXTEND_JUMP_READS=0
GRAPH_KMER_SIZE = auto
USE_LINKING_MATES = 1
USE_GRID=0
GRID_ENGINE=SLURM
GRID_QUEUE=m3tb
LHE_COVERAGE=25
CA_PARAMETERS =  cgwErrorRate=0.15
CLOSE_GAPS=1
NUM_THREADS = 48
JF_SIZE = 250000000000
SOAP_ASSEMBLY=1
FLYE_ASSEMBLY=0
END
```

Run MaSurca

``` bash
module load masurca/4.0.5

# The job command(s)

masurca config.txt
./assemble.sh
```

Total scaffolds: 4713

Mapping raw reads

``` bash
bowtie/2.3.4
samtools/1.9.0
bbtools/38.37

#Index bowtie2
bowtie2-build --threads 16 scaffolds.fasta masurca

#Run mapping
bowtie2 --no-unal -p 16 -x masurca -1 Mp_R1.fastq.gz -2 Mp_R2.fastq.gz -S masurca.sam &>bowtie.log

#index scaffolds
samtools faidx scaffolds.fasta

#Convert sam to bam
samtools view -@ 16 -bt scaffolds.fasta masurca.sam | samtools sort -@ 16 -o masurca.bam

#coverage info per scaffold
pileup.sh -Xmx30g in=bmasurca.bam ref=scaffolds.fasta out=mapping_stats_masurca.txt
```

## Cleaning genome

#### Cleaning with taxon classification using Blobtools v1

The filter steps were adapted from Iha et al. 2021

**Make custom database**

I used non-redundant (NR) database from GenBank.

Download
[prot.accession2taxid.gz](https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz)

Get taxonID for the selected lineages:

-   2 -\> bacteria

-   2157 -\> archeae

-   4751 -\> fungi

-   10239 -\> virus

-   33634 -\> Stramenopiles

    -   2696291 -\> Ochrophyta
    -   2870 -\> Phaeophyceae
    -   2836 -\> Bacillariophyta

-   2763 -\> Rhodophyta

-   3041 -\> Chlorophyta

``` bash
#Using TaxonKit (https://bioinf.shenwei.me/taxonkit/)

# Create a list of taxids
taxonkit list --ids 2,2157,4751,10239,33634,2763,3041 --indent "" > lineages.taxid.txt

# Total taxIDs:
wc -l lineages.taxid.txt

# Retrieving target accessions
pigz -dc prot.accession2taxid.gz | csvtk grep -t -f taxid -P lineages.taxid.txt | csvtk cut -t -f accession.version,taxid | sed 1d > lineages.acc2taxid.txt

cut -f 1 lineages.acc2taxid.txt > lineages.acc.txt

# Target retrieved:
wc -l lineages.acc.txt
```

Split lineages.acc.txt in many files

Retrieving FASTA from NR - using array in SLURM

``` bash
#SBATCH --output=array_%A-%a.out
#SBATCH --array=0-1744

# Retrieving FASTA

blastdbcmd -db nr -entry_batch lineages.acc_${SLURM_ARRAY_TASK_ID} -out - | pigz -c > nr_${SLURM_ARRAY_TASK_ID}.db.fa.gz
```

Concatenate all sequences to create a unique database

``` bash
cat nr_*.db.fa.gz > nr.db.fa.gz

# Counting sequences:
pigz -dc nr.db.fa.gz | grep '>' | wc -l
```

**Prepare data to Blobtools v1:**

``` bash
blast+/2.12.0
diamond/2.0.4
samtools/1.9.0

#Create database with Uniprot
diamond makedb --in uniprot_ref_proteomes.fasta -d uniprot_ref_proteomes.diamond

#Running diamond to Uniprot
diamond blastx --query scaffolds.fasta --max-target-seqs 10 --sensitive --threads 16 --db uniprot_ref_proteomes.diamond.dmnd --evalue 1e-25 --outfmt 6 --out masurca.vs.uniprot.1e25.diam.out

#Running Blast to Silva
blastn -task megablast -query scaffolds.fasta -db database/Silva/silva.rDNA.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads 16 -evalue 1e-65 -out masurca.vs.silva.1e65.megablast.out

#Running Blast to Custom database
#Making db to custom database
gunzip database/Custom_db/nr.db.fa.gz
diamond makedb -p 16 --in database/Custom_db/nr.clean.db.fa -d database/Custom_db/nr.db.fa.diamond

#Running Diamond to Custom database
diamond blastx --query scaffolds.fasta --max-target-seqs 10 --sensitive --threads 16 --db database/Custom_db/nr.db.fa.diamond.dmnd --evalue 1e-20 --outfmt 6 --out masurca.vs.customdb.1e20.diamond.out

#taxify for Uniprot
blobtools taxify -f masurca.vs.uniprot.1e25.diam.out -m database/Uniprot/uniprot_ref_proteomes.taxids -s 0 -t 2

#taxify for Silva
blobtools taxify -f masurca.vs.silva.1e65.megablast.out -m database/Silva/silva_ref_rRNA.taxids -s 0 -t 2

#taxify for Custom database
blobtools taxify -f masurca.vs.customdb.1e20.diamond.out -m database/Custom_db/lineages.acc2taxid.txt -s 0 -t 1

#Concat all hits
cat masurca.vs.uniprot.1e25.diam.taxified.out masurca.vs.silva.1e65.megablast.taxified.out masurca.vs.customdb.1e20.diamond.taxified.out > hits.tsv
```

Run Blobtools

``` bash
blobtools create -i scaffolds.fasta -t hits.tsv --nodes blobtools/nodes/nodes.dmp --names blobtools/nodes/names.dmp -b masurca.bam -o masurca_blobtools

#Running Blobtools view
blobtools view -i masurca_blobtools.blobDB.json -o masurca_blobview

#Running Blobtool plot
blobtools plot -i  masurca_blobtools.blobDB.json
```

![BlobTools
plot](Masurca/masurca_blobtools14102021.blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png)

**Removing Eukaryotes scaffolds that are distant of Stramenopiles**

``` bash
cut -f6 masurca_blobview.blobdb.blobDB.bestsum.table.txt | sort -u > 00_taxa_to_remove
#manipulate on vi
```

Kept only the bold taxa:

**Actinobacteria** Arthropoda **Bacillariophyta** **Bacteria-undef**
**Bacteroidetes** **Balneolaeota** Basidiomycota **Candidatus
Tectomicrobia** **Chloroflexi** Chlorophyta **Cyanobacteria**
**Eukaryota-undef** Evosea **Firmicutes** **Gemmatimonadetes**
**no-hit** Oomycota **Planctomycetes** **Proteobacteria** Streptophyta
**undef** **Verrucomicrobia**

``` bash
grep -f 00_taxa_to_remove masurca_blobview.blobdb.blobDB.bestsum.table.txt | cut -f1 > 00_scaffolds_to_remove_byTaxa.ID
```

19 scaffolds to remove

``` bash
grep '>' scaffolds.fasta | grep -v -w -f 00_scaffolds_to_remove_byTaxa.ID | sed 's/>//' > 00_scaffolds_to_keep.ID

seqtk subseq scaffolds.fasta 00_scaffolds_to_keep.ID > 00_scaffolds.clean.fasta
#4694 scaffolds
```

Get scaffolds that matched with Ectocarpus.

Taxonomic ID: 2880 Ectocarpus siliculosus 867726 Ectocarpus sp. CCAP
1310/34

``` bash
grep -w -f 00_all_scaffolds.ID ../blobtools/masurca.vs.customdb.diamond.blastx.out | grep -e '2880' -e '867726' | cut -f1 | sort -u > 01_all_Ectocarpus_scaffolds.ID

grep -w -f 01_all_noEctocarpus_scaffolds.ID masurca_blobview.blobdb.blobDB.bestsum.table.txt | grep 'Eukaryota' | cut -f1 > 01_eukaryota.ID

#check this IDs in Uniprot
grep -f 01_eukaryota.ID ../blobtools/masurca.vs.uniprot.1e25.diam.taxified.out

#Total 1355 scaffolds matched with Ectocarpus

# Scaffols do not match with Ectocarpus
grep -v -w -f 01_all_Ectocarpus_scaffolds.ID 00_all_scaffolds.ID > 01_noEctocarpus.ID
# 3339 do not match with Ectocarpus
```

**Using transcriptomes to clean the genome**

``` bash
minimap2 -ax splice -C 5 --splice-flank=no --secondary=no -t 16 scaffolds.cleantaxaless1000l.fasta all_transcripts.fasta > mapped_transcripts.sam

#transform in paf format to get CIGAR string info for introns
./paftools.js sam2paf mapped_transcripts.sam > mapped_transcripts.paf
```

CIGAR string explained:

cg:Z:130M11864N224M3939N188M –\> **N** means introns!

2586 scaffold presents introns

**Calculate outliers** Script on R –\> Outrliers.R

**Final genome: scaffolds.clean.fasta**

**3864 scaffolds**

## Assembly Transcriptomes

We used two transcriptomes downloaded from GenBank

*Macrocystis integrifolia* transcriptome:
[PRJNA322132](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA322132)

Reference:
<https://link.springer.com/article/10.1007/s13205-018-1204-4#Sec7>

RNAseq for the brown alga *Macrocystis pyrifera*:
[PRJNA353611](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA353611)

#### Download data from GenBank

``` bash
prefetch <SRA_accession_number>
fastq-dump --split-3 <SRA_accession_number> #Descompress SRA and split in pair-ends 
```

Check quality with FastQC and trim reads with Trimmomatic.

**Two modes were used: de novo and “genome-guided”**
[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)

#### *De novo*

``` bash
Trinity \
        --seqType fq \
        --max_memory 250G \
        --samples_file samples.txt \
        --CPU 20 \
        --output trinity_out_dir
```

#### *Genome guided*

``` bash
module load star/2.7.9a
module load samtools/1.12
module load trinity/2.12.0

#Genome index
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir genome --genomeFastaFiles scaffolds.clean.fasta --genomeSAindexNbases 11

#Mapping reads
STAR --runThreadN 16 --genomeDir genome --readFilesManifest manifest.tsv

#Sam to bam
samtools view -@ 8 -b -t scaffolds.clean.fasta Aligned.out.sam | samtools sort -@ 8 -o Aligned.out.bam

#assembly
Trinity --genome_guided_bam Aligned.out.bam \
        --genome_guided_max_intron 10000 \
        --max_memory 250G \
        --CPU 20
```

|              |    PRJNA322132    |    PRJNA353611    |
|:-------------|:-----------------:|:-----------------:|
| De novo      | Trinity-DNP.fasta | Trinity-DNR.fasta |
|              |      124436       |      310167       |
| Genome Guide | Trinity-GGP.fasta | Trinity-GGR.fasta |
|              |       14586       |       34452       |

Concat all transcriptomes.

## Genome annotation

Workflow based on [Chen et
al. 2019](https://onlinelibrary.wiley.com/doi/full/10.1111/jpy.12947)

<https://github.com/TimothyStephens/Dinoflagellate_Annotation_Workflow>

#### Clean Transcriptomes

[Univec](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/) database

Univec will be used to clean transcriptomes

``` bash
#blast+/2.12.0
#pasa/2.5.1

mkdir UNIVEC/
cd UNIVEC/

wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core
wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/README.uv


makeblastdb  -in UniVec -dbtype nucl
makeblastdb  -in UniVec_Core -dbtype nucl

seqclean all_transcripts.fasta -v UNIVEC/UniVec_Core -c 4

#Cleaned transcript: all_transcripts.fasta.clean
# 483627 transcripts
```

#### Repeat analysis

``` bash

#repeatmasker/4.1.2p1
#rmblast/2.11.0
#repeatmodeler/2.0.2a


### If need TRF ####
#Download Tandem Repeats Finder
#https://tandem.bu.edu/trf/trf409.linux64.download.html
mv trf409.legacylinux64 trf
chmod +x trf
mv trf ~/bin
####################

export GENOME_NAME=scaffolds.clean.fasta
export GENOME_SIZE=51923561

#Running Repeat Modeler

BuildDatabase -name ${GENOME_NAME}_db -engine ncbi ${GENOME_NAME} > BuildDatabase.out
RepeatModeler -engine ncbi -pa ${SLURM_NTASKS} -database ${GENOME_NAME}_db 1> RepeatModeler.sdtout 2> RepeatModeler.sdterr

# Combine repeat modeler and repeat masker libraries.
cp RM_*/consensi.fa* .
cat consensi.fa.classified RepeatMasker.lib > ${GENOME_NAME}_CombinedRepeatLib.lib

#### Run Repeat Masker

RepeatMasker -lib ${GENOME_NAME}_CombinedRepeatLib.lib -e ncbi -gff -x -no_is -a -pa ${SLURM_NTASKS} ${GENOME_NAME} 1> RepeatMasker.sdtout 2> RepeatMasker.sdterr


# SoftMask Genome
maskFastaFromBed -soft -fi ${GENOME_NAME} -fo ${GENOME_NAME}.softmasked -bed ${GENOME_NAME}.out.gff

# Get repeat features in gff3 format (used for EVM later on)
rmOutToGFF3.pl ${GENOME_NAME}.out > ${GENOME_NAME}.out.gff3

# Model repeat landscape (not needed for annotation)
calcDivergenceFromAlign.pl -s ${GENOME_NAME}.divsum ${GENOME_NAME}.align
createRepeatLandscape.pl -div ${GENOME_NAME}.divsum -g ${GENOME_SIZE} > ${GENOME_NAME}.html
```

#### PASA

``` bash
module load perl/5.32.1
module load blat/3.6
module load samtools/1.12
module load pasa/2.5.1

# The job command(s):
export PASA=/apps/pasa/2.5.1/
export GENOME_NAME=scaffolds.clean.fasta
export GENOME_SIZE=51923561
export MAX_INTRON_LENGTH=70000
export SCRIPTS=/scratch1/iha002/WGS_Illumina_AGRF/AGRF_CAGRF21056777_HGHMYDRXY/annotation/Workflow/Dinoflagellate_Annotation_Workflow/scripts

Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g ${GENOME_NAME} --MAX_INTRON_LENGTH ${MAX_INTRON_LENGTH} --ALIGNERS blat --CPU ${SLURM_NTASKS} -T -t all_transcripts.fasta.clean -u all_transcripts.fasta > pasa.log 2>&1

perl $PASA/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta ${GENOME_NAME}_pasadb.sqlite.assemblies.fasta --pasa_transcripts_gff3 ${GENOME_NAME}_pasadb.sqlite.pasa_assemblies.gff3

perl ${SCRIPTS}/GeneStats.pl ${GENOME_NAME}_pasadb.sqlite.assemblies.fasta.transdecoder.cds ${GENOME_NAME}_pasadb.sqlite.assemblies.fasta.transdecoder.genome.gff3 ${GENOME_NAME} > ${GENOME_NAME}_pasadb.sqlite.assemblies.fasta.transdecoder.cds.stats
```

#### Blast PASA output to NR

``` bash
mkdir blast2db

#Link proteins predicted
ln -s scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep .

## Get sequences which are type:complete.
grep '>' ${GENOME_NAME}_pasadb.sqlite.assemblies.fasta.transdecoder.pep | awk '{print $1}' | sed -e 's/>//' > Complete_Sequences.ids

## Get all seq IDs that have only  one CDS.
awk '$3=="CDS"{print $0}' ${GENOME_NAME}_pasadb.sqlite.assemblies.fasta.transdecoder.genome.gff3 | awk '{ print $9}' | sed -e 's/.*Parent=\(.*\)/\1/' | sort | uniq -c | awk '{ if($1==1) print $2 }' > Single_CDS_Genes.ids

## Get seq IDs which have coords on the genome.
awk '$3=="mRNA"{print $0}' ${GENOME_NAME}_pasadb.sqlite.assemblies.fasta.transdecoder.genome.gff3 | awk '{ print $9 }' | sed -e 's/ID=\(.*\);Parent.*/\1/' > Genes_with_genome_coords.ids

## Filter IDs

## Get seq IDs that are NOT Single Exon genes, have genome coords, and are type complete. 
## The later GoldenGenes stage will filter these genes anyway so we might as well filter them out now before the intensive blast stage.

python2 ../Dinoflagellate_Annotation_Workflow/scripts/filter_ids.py -i Complete_Sequences.ids -o Complete_Sequences.ids.filtered -k Genes_with_genome_coords.ids -r Single_CDS_Genes.ids

## Get pep sequences in filtered ID list
xargs samtools faidx ${GENOME_NAME}_pasadb.sqlite.assemblies.fasta.transdecoder.pep < Complete_Sequences.ids.filtered > ${GENOME_NAME}_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered

mkdir fasta_split/

seqkit split2 -s 50 scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered
ls -1 fasta_split/*.filtered > files2run.txt
```

Run BLAST array

``` bash
#SBATCH --array=0-48

# The modules to load:
module load blast+/2.12.0
module load bioref/default

# The job command(s):

if [ -d "$INDIR" ]
then
        SAMPLES=( `ls -1 ${INDIR}/*.filtered | sed -E 's/.+\/(.+)\.filtered/\1/' ` );

        if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
        then
                i=$SLURM_ARRAY_TASK_ID
                blastp -query ${INDIR}/${SAMPLES[$i]}.filtered -db nr  -out ${INDIR}/${SAMPLES[$i]}.blastp.outfmt6 -num_threads ${SLURM_NTASKS} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
        else
                echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
        fi
else
        echo "Error: Missing input file directory as --export env INDIR or doesn't exist"
fi
```

Sorting results

``` bash
cat FASTA_SPLIT/*.outfmt6 > blastp.outfmt6

cat blastp.outfmt6 | awk -F "\t" '$11<1e-20' | awk -F "\t" '{ if (( (($8-$7)+1) / $13) > 0.8) {print $1}}' | sort | uniq | xargs samtools faidx scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered > scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered.top_coverage_onlyquery.faa  #### 1624 sequences
```

#### Transposon-PSI

Download [Transposon-PSI](http://transposonpsi.sourceforge.net/)
Download
[blastall](wget%20ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz).
The program use formatdb.

Create a \~/.ncbirc on home:

    [NCBI]
    data=<full-path-to-TransposonPSI-data>/Transposon-PSI/blast-2.2.26/data

``` bash

perl TransposonPSI_08222010/transposonPSI.pl scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered.top_coverage.faa prot

cat scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered.top_coverage.faa.TPSI.allHits | awk '{ print $5 }' | sort | uniq > TransposonPSI.hit.seq.ids
```

#### CD-HITs

``` bash

mkdir CDHITS
cd CDHITS/

grep '>' scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered.top_coverage.faa | sed -e 's/>//' > Seqs_from_blast.ids

cp ../Dinoflagellate_Annotation_Workflow/scripts/filter_ids.py .

python2 filter_ids.py -i Seqs_from_blast.ids -o Seqs_from_blast.ids.filtered -r ../Transposon-PSI/TransposonPSI.hit.seq.ids

module load samtools/1.12

xargs samtools faidx scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered.top_coverage.faa < Seqs_from_blast.ids.filtered > scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered.top_coverage.final.faa


module load cd-hit/4.8.1
cd-hit \
   -i scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered.top_coverage.final.faa \
   -o scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered.top_coverage.final.faa.cdhit75 \
   -c 0.75 \
   -n 5 \
   -T 16 \
   -M 30000


grep '>' scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered.top_coverage.final.faa.cdhit75 | sed -e 's/>//' > scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered.top_coverage.final.faa.cdhit75.ids
```

#### Golden Genes

``` bash

git clone https://github.com/genomecuration/JAMg.git
date > JAMg/date_retrieved.txt

mkdir Prepare_Golden_Genes

cp JAMg/date_retrieved.txt Prepare_Golden_Genes/
cp JAMg/bin/prepare_golden_genes_for_predictors.pl Prepare_Golden_Genes/prepare_golden_genes_for_predictors.pl
cp JAMg/bin/run_exonerate.pl Prepare_Golden_Genes/run_exonerate.pl
cp -r JAMg/PerlLib Prepare_Golden_Genes/

export PERL5LIB=$PERL5LIB:"/home/iha002/iha002/WGS_Illumina_AGRF/AGRF_CAGRF21056777_HGHMYDRXY/annotation/Prepare_Golden_Genes/PerlLib"

#Install BIO
perl -MCPAN -Mlocal::lib -e shell
cpan[1]> install Bio::Perl

#Install cdbfasta
git clone https://github.com/gpertea/cdbfasta.git
make
cd ~/bin
ln -s ../iha002/WGS_Illumina_AGRF/AGRF_CAGRF21056777_HGHMYDRXY/annotation/Prepare_Golden_Genes/cdbfasta/cdbfasta
ln -s ../iha002/WGS_Illumina_AGRF/AGRF_CAGRF21056777_HGHMYDRXY/annotation/Prepare_Golden_Genes/cdbfasta/cdbyank

## Check it works
module load perl/5.32.1
module load blast+/2.12.0
module load gmap/20210825
module load samtools/1.12
module load exonerate/2.2.0
module load augustus/3.4.0
module load snap-korflab/211201
module load trinity/2.12.0
module load emboss/6.6.0


export PERL5LIB=$PERL5LIB:"<PATH-to-PerlLIB>/Prepare_Golden_Genes/PerlLib"

perl prepare_golden_genes_for_predictors.pl \
        -genome scaffolds.clean.fasta.masked \
        -softmasked scaffolds.clean.fasta.softmasked \
        -no_gmap -same_species -intron 7000 -cpu 16 -norefine -complete -no_single -verbose \
        -pasa_gff *.assemblies.fasta.transdecoder.gff3 \
        -pasa_peptides *.assemblies.fasta.transdecoder.pep \
        -pasa_cds *.assemblies.fasta.transdecoder.cds \
        -pasa_genome *.assemblies.fasta.transdecoder.genome.gff3 \
        -pasa_assembly *.assemblies.fasta
```

Running golden genes script

``` bash
mkdir Golden_genes/

ln -s ../repeat/scaffolds.clean.fasta.masked .
ln -s ../repeat/scaffolds.clean.fasta.softmasked .
ln -s ../CDHITS/scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.pep.complete_only.filtered.top_coverage.final.faa.cdhit75.ids CD-HIT.ids.txt
ln -s ../../Prepare_Golden_Genes/prepare_golden_genes_for_predictors.pl .

## Get PASA output for GoldenSet.
## 
## This script should get just the IDs given in the CD-HIT.ids.txt file in
## the same format as you would expect from PASA. This seems to be required by the 
## prepare_golden_genes script for some reason.
## This script might not work with different versions of PASA if the seq names are in a different format.
python2 ../Dinoflagellate_Annotation_Workflow/scripts/get_PASA_for__prepare_golden_genes.py --ids CD-HIT.ids.txt --pasa_assembly ../PASA/*.assemblies.fasta --pasa_cds ../PASA/*.assemblies.fasta.transdecoder.cds --pasa_peptides ../PASA/*.assemblies.fasta.transdecoder.pep --pasa_gff ../PASA/*.assemblies.fasta.transdecoder.gff3 --pasa_genome ../PASA/*.assemblies.fasta.transdecoder.genome.gff3

Run golden_genes.script
# I had to install Parafly
```

#### GeneMark-ES

``` bash
perl gmes_linux_64/gmes_petap.pl --ES --cores 16 --format GFF3 --v --sequence scaffolds.clean.fasta.masked 1>&2> genemark.log

perl ${MAKER}/genemark_gtf2gff3 genemark.gtf > genemark.gff3
```

#### SNAP

``` bash
mkdir SNAP

ln -s ../repeat/scaffolds.clean.fasta.softmasked .
ln -s ../Golden_genes/final_golden_genes.gff3.nr.golden.zff .

module load samtools/1.12
module load snap-korflab/211201

grep '>' final_golden_genes.gff3.nr.golden.zff | sed 's/>\(.*\)/\1/' | xargs samtools faidx scaffolds.clean.fasta.softmasked > snap.fasta

fathom final_golden_genes.gff3.nr.golden.zff snap.fasta -gene-stats
#118 sequences
#0.506642 avg GC fraction (min=0.465218 max=0.544290)
#142 genes (plus=77 minus=65)
#0 (0.000000) single-exon
#142 (1.000000) multi-exon
#190.731445 mean exon (min=3 max=1576)
#1170.067871 mean intron (min=90 max=6611)

fathom final_golden_genes.gff3.nr.golden.zff snap.fasta -validate > validate.txt

fathom final_golden_genes.gff3.nr.golden.zff snap.fasta -categorize 1000

cat uni.ann wrn.ann > com.ann
cat uni.dna wrn.dna > com.dna

fathom -export 1000 -plus com.ann com.dna

forge export.ann export.dna

hmm-assembler.pl -c 0.000 scaffolds.clean.fasta . > scaffolds.clean.fasta.snap.hmm

mkdir PREDICT; cd PREDICT/

ln -s ../scaffolds.clean.fasta.softmasked .

snap ../scaffolds.clean.fasta.snap.hmm scaffolds.clean.fasta.softmasked -lcmask -quiet 1> stdout.snap 2> stderr.snap

perl ../../Dinoflagellate_Annotation_Workflow/scripts/SNAP_output_to_gff3.pl stdout.snap scaffolds.clean.fasta.softmasked 1> snap.gff3 2> snap.gff3.strerr
```

#### AUGUSTUS

I didn’t used Training mode

``` bash

mkdir PREDICTION

module load augustus/3.4.0

ln -s ../../repeat/scaffolds.clean.fasta.softmasked .

cp ../../Dinoflagellate_Annotation_Workflow/scripts/fasta-splitter.pl .
#Had to install File::Util - perl module

perl fasta-splitter.pl --n-parts 100 --out-dir FASTA_SPLIT ${GENOME_NAME}.softmasked

augustus --softmasking=1 --gff3=on --UTR=off --codingseq=on --protein=on --exonnames=on --species=${SPECIES} ${INDIR}/${SAMPLES[$i]}.softmasked > ${INDIR}/${SAMPLES[$i]}.augustus.out

cat FASTA_SPLIT/*.out | join_aug_pred.pl > augustus.gff3
getAnnoFasta.pl augustus.gff3
#augustus3.aa
#augustus3.codingseq

#Install Evidence modeler
export PERL5LIB=$PERL5LIB:"/home/iha002/sw/EVidenceModeler-1.1.1/PerlLib"
cp ~/sw/EVidenceModeler-1.1.1/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl .

perl augustus_GFF3_to_EVM_GFF3.pl augustus.gff3 | sed '/^$/d' > augustus.cleaned.gff3
perl GeneStats.pl augustus3.codingseq augustus.cleaned.gff3 ../../scaffolds.clean.fasta > augustus.cleaned.gff3.stats
```

#### Evidence Modeler (Finally!!!)

Prepare data:

``` bash

mkdir EVIDENCE_MODELER

cd EVIDENCE_MODELER/

ln -s ../scaffolds.clean.fasta .
ln -s ../GeneMark/genemark.gff3 .
ln -s ../SNAP/PREDICT/snap.gff3 .
ln -s ../AUGUSTUS/augustus.cleaned.gff3 .
ln -s ../PASA/scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.genome.gff3 .
ln -s ../repeat/scaffolds.clean.fasta.out.gff .

grep -v "^#" genemark.gff3 > abinitio_gene_predictions.gff3
grep -v "^#" snap.gff3 | sed '/^$/d' | awk '{print $1 "\tSNAP\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t"$8 "\t" $9}' >> abinitio_gene_predictions.gff3
cat augustus.cleaned.gff3 >> abinitio_gene_predictions.gff3
sed '/^$/d' scaffolds.clean.fasta_pasadb.sqlite.assemblies.fasta.transdecoder.genome.gff3 >> abinitio_gene_predictions.gff3
```

Make weights.txt

    ABINITIO_PREDICTION     GeneMark.hmme 2
    ABINITIO_PREDICTION     SNAP    2
    ABINITIO_PREDICTION     Augustus        6
    OTHER_PREDICTION        transdecoder    10

Run Evidence Modeler

``` bash

mkdir PARTITIONS_EVM
cd PARTITIONS_EVM/

ln -s ../abinitio_gene_predictions.gff3 .

partition_EVM_inputs.pl --genome scaffolds.clean.fasta --gene_predictions abinitio_gene_predictions.gff3 --segmentSize 50000000 --overlapSize 10000 --partition_listing partitions_list.out

write_EVM_commands.pl --genome scaffolds.clean.fasta --gene_predictions abinitio_gene_predictions.gff3 --weights weights.txt --output_file_name evm.out --partitions partitions_list.out > commands.list

ParaFly -v -CPU 16 -c commands.list -failed_cmds commands.list.failed

recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output evm.out --genome /scratch2/iha002/Macrocystis_annotation/scaffolds.clean.fasta

cd ../

find PARTITIONS_EVM/ -name evm.out.gff3 -exec cat '{}' ; > scaffolds.clean.fasta.evm.gff3


# Get CDS and Proteins.

#sed 's@EVM.evm.TU@EVM.evm.TU@' will remove ^\ (file seperator) character from EVM^\evm.TU to EVM.evm.TU

gff3_file_to_proteins.pl scaffolds.clean.fasta.evm.gff3 scaffolds.clean.fasta CDS | sed 's@EVM.evm.TU@EVM.evm.TU@' > scaffolds.clean.fasta.evm.cds.fna

gff3_file_to_proteins.pl scaffolds.clean.fasta.evm.gff3 scaffolds.clean.fasta prot | sed 's@EVM.evm.TU@EVM.evm.TU@' > scaffolds.clean.fasta.evm.protein.faa
```

## DONE!
