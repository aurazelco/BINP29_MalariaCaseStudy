# BINP29_MalariaCaseStudy
Repo containing the scripts used to analyze malaria parasites' genome in the BINP29 course (Lund University):
https://github.com/aurazelco/BINP29_MalariaCaseStudy. 

All analysis is executed on the server.

We start from genomes of malaria parasites, and the endpoint is to build phylogenetic trees using parasites with different hosts, with a particular focus on *Haemoproteus tartakovskyi*. 

The novel python scripts (removeScaffold.py, removeBird.py, GCcontent.py and buscoParser.py) are present in the Github repo. 

## Data Collection

We first create a parent folder, and download the data:
```shell
mkdir -p malaria_project && cd malaria_project  # this will be the parent directory of all analysis
# download all genomes
tar -xvzf plasmodiumGenomes.tgz
```
## *Plasmodium* data - Gene Prediction

We run GeneMark on a genome of choice, in this case *Plasmodium yoelii*. We first check GeneMark version, whihc will be the same throughout the wholoe analysis:
```shell
gmes_petap.pl
#GeneMark-ES Suite version 4.62_lic
```
We run it in a newly created folder:
```shell
mkdir 1_GeneMark && cd 1_GeneMark
nohup gmes_petap.pl --sequence ../Plasmodium_yoelii.genome &
```
After uploading it to the tmp fodler, we proceed to work on the *Haemoproteus tartakovskyi* data. 

## Processing of *Haemoproteus tartakovskyi* (Ht) data
### Clean genome sequence

*Haemoproteus tartakovskyi* is a parasite of birds; therefore, even if the parasite genome was enriched in the experimental settings, we still need to remove the sequences and scaffolds from the hosts before running the phylogenetic analysis. 

From this genome, we need to remove the sequences which have a GC content higher than a certain threshold, and the sequences shorter than 3000 nucleotides. To do so, we run the remoevScaffold.py script. 

But first, we make two new directories, one which contains all python scripts and the other to store the results of removeScaffold.py
```shell
mkdir -p Scripts && cd Scripts
mkdir -p 2_RemoveScaffolds && cd 2_RemoveScaffolds
```
The removeScaffold.py script is run as (see repo for full code):
```shell
python ../Scripts/removeScaffold.py -i ../Haemoproteus_tartakovskyi.genome -o Haemoproteus_tartakovskyi.output
```
### Make a gene prediction
We run again GeneMark, this time on the 'cleaned' Ht genome. 
```shell
mkdir -p 3_GenePredictionHaemoproteus && cd 3_GenePredictionHaemoproteus
nohup gmes_petap.pl --sequence ../2_RemoveScaffolds/Haemoproteus_tartakovskyi.output &
```
Once this is run, we would like to run the gffParse.pl script. However, there is a discrepancy on the scaffold names; therefore, we use sed to remove the extra fields. 
```shell
mv genemark_es.gtf Haemoproteus.gff
cat Haemoproteus.gff | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Ht2.gff
cd ..
```
### gffParse

Once that is done, we retrieve the version of gffParse.pl. 
```shell
mkdir -p 4_GTFtoFASTA && cd 4_GTFtoFASTA
gffParse.pl -v
#gffParse.pl version 1.1
```
The gffParse.pl script is run as:
```shell
gffParse.pl -i ../2_RemoveScaffolds/Haemoproteus_tartakovskyi.output  -g ../3_GenePredictionHaemoproteus/Ht2.gff -c -p
```
### BLAST

Once we obtain the gff files, we would like to BLAST so that we can remove the sequences which align to birds. Theerfore, we run BLASTX on the gffParse.fna file we obtained in the pervious step. 

```shell
cd ..
mkdir -p 5_BLAST && cd 5_BLAST
blastx -query ../4_GTFtoFASTA/gffParse.fna -db SwissProt -out Haemoproteus.blastx -num_threads 20 -evalue 1e-10
# version from BLASTX output: BLASTX 2.11.0+
```
Once the BLASTX analysis is over, we need to remove the sequences which have as first BLASTX hit a sequence belonging to an avian organism. To do so, we use the BLASTX results, the original gffParse FASTA file and two databases for taxonomy and uniprot. 
The databases contain all scientific names (taxonomy.dat) and the genes id (uniprot_sprot.dat), which are used to extract the scaffolds from the gffParse.fna file which had a top BLAST hit belonging to a bird. I used the provided datParser.py script, which runs as:
```shell
python datParser.py Haemoproteus.blastx ../4_GTFtoFASTA/gffParse.fna taxonomy.dat uniprot_sprot.dat > bird_scaffolds.txt
```
The scaffolds to be reomves are saved in bird_scaffolds.txt. 

### Clean again genome from bird scaffolds

Now that we have the scaffolds to be removed, we use a python script called removeBird.py which removes these scaffolds from the Haemoproteus_tartakovskyi.output, and saves the newly cleaned genome in a new output. 

The python script run as (see repo for full code):
```shell
cd ..
mkdir -p 6_FilteredGeneMark && cd 6_FilteredGeneMark
python ../Scripts/removeBird.py -i ../2_RemoveScaffolds/Haemoproteus_tartakovskyi.output -s ../5_BLAST/bird_scaffolds.txt -o Haemoproteus_clean.genome
```
### Gene Prediction and gffParse on the clean genome
We run the GeneMark script again, as before:
```shell
nohup gmes_petap.pl --sequence Haemoproteus_clean.genome &
mv genemark_es.gtf Haemoproteus_clean.gtf # renamed files
```
And the gffParse.pl:
```shell
cd ..
mkdir -p 7_CleanGTFtoFASTA && cd 7_CleanGTFtoFASTA
gffParse.pl -i ../6_FilteredGeneMark/Haemoproteus_clean.genome -g ../6_FilteredGeneMark/Haemoproteus_clean.gtf -c -p
```

#### Question 3 - genome length
```shell
cd ..
cat 6_FilteredGeneMark/Haemoproteus_clean.genome | grep -v "^>" | tr -d "\n" | wc -c
# 22068498
for file in *.genome; do echo $file; grep -v "^>" $file | tr -d "\n" | wc -c; done
#Haemoproteus_tartakovskyi.genome -> this is before filtering (not included in final table)
# 22609283
#Plasmodium_berghei.genome
#17954629
#Plasmodium_cynomolgi.genome
#26181343
#Plasmodium_faciparum.genome
#23270305
#Plasmodium_knowlesi.genome
#23462346
#Plasmodium_vivax.genome
#27007701
#Plasmodium_yoelii.genome
#22222369
#Toxoplasma_gondii.genome
#128105889
```
#### Question 3 - number of genes
```shell
mkdir -p 8_allGTFs && cd 8_allGTFs
# after copying the gtf files
# just changed the name so that all files are .gtf
mv Tg.gff toxoplasma_gondii.gtf
# copy here
cp ../6_FilteredGeneMark/Haemoproteus_clean.gtf .

# number of genes
for file in *.gtf; do echo $file; cut -f9 $file | cut -d ";" -f1 | sort -u | wc -l ; done
# Haemoproteus_clean.gtf
# 4429
# knowlesi.gtf
# 4953
# plasmodium_berghei.gtf
# 7282
# plasmodium_cynomolgi.gtf
# 5787
# plasmodium_faciparum.gtf
# 5207
# plasmodium_vivax.gtf
# 5682
# plasmodium_yoelii.gtf
# 4919
# toxoplasma_gondii.gtf
# 15892
```

#### Question 3 - GC content
To calculate the GC content, I wrote a python script, GCcontent.py, which runs as python Scripts/GCcontent.py -path 9_Genomes/ and prints the result to the console (see repo for full code). 

```shell
mkdir -p 9_Genomes
# moved all genomes to a new folder for better tree organization

python Scripts/GCcontent.py -path 9_Genomes/
#9_Genomes/Plasmodium_berghei.genome
#23.72%
#9_Genomes/Plasmodium_yoelii.genome
#21.77%
#9_Genomes/Plasmodium_vivax.genome
#42.28%
#9_Genomes/Toxoplasma_gondii.genome
#52.35%
#9_Genomes/Plasmodium_knowlesi.genome
#38.83%
#9_Genomes/Plasmodium_cynomolgi.genome
#40.38%
#9_Genomes/Haemoproteus_clean.genome
#26.02%
#9_Genomes/Plasmodium_faciparum.genome
#19.36%
```

All these results can be found summarized in the main.txt file. 

## Phylogenetic trees
### Identify orthologs
Now that we have clean Ht genome (hereafter called Hc), we can run proteinortho and BUSCO to find orthologs. 

For each genome, we need to run gffParse.pl:
```shell
mkdir 10_allGeneMark
# renames file
mv 8_allGTFs/knowlesi.gtf 8_allGTFs/plasmodium_knowlesi.gtf

for file in 8_allGTFs/*; do echo $file; genome=$(echo $file | cut -d '/' -f2 | cut -d '.' -f1); genus=$(echo $file | cut -d '/' -f2 | cut -c 1); species=$(echo $file | cut -d '_' -f3 | cut -c 1); output=(${genus^}$species); gffParse.pl -i 9_Genomes/${genome^}'.genome' -g $file -c -p -b 10_allGeneMark/$output ; done
```

Now all gff results are in the 10_allGeneMark folder,relabelled as Hc,Pb,Pc,Pf,Pk,Pv,Py,Tg. 


### Installation of ProteinOrtho and BUSCO
```shell
conda create -n malariaenv
conda activate malariaenv
conda install proteinortho
conda update --name base conda
conda deactivate
conda create -n malariaBusco -c conda-forge -c bioconda busco=5.3.0
```

#### ProteinOrtho
We create a new directory, we activate the conda environment for proteinortho and then we run it as follows:
```shell
mkdir 11_ProteinOrtho && cd 11_ProteinOrtho
conda activate malariaenv
proteinortho6.pl -h
# PoFF version 6.0.33
# I have Hc instead of Ht because the names are Haemo..._clean
nohup proteinortho6.pl ../10_allGeneMark/{Hc,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa &
conda deactivate
# how many orthologs we have?
grep -v '^#'  11_ProteinOrtho/myproject.proteinortho.tsv | wc -l
# 5126
```

#### BUSCO
We change the environment to the one where we installed BUSCO, and then run BUSCO for all protein files (.faa) obtained from the last gffParse.pl command. 
```shell
conda activate malariaBusco
mkdir -p 12_BUSCO

for file in 10_allGeneMark/*.faa; do echo $file; output=$(echo $file | cut -d '/' -f2 | cut -d '.' -f1); echo $output; busco -i $file -o 12_BUSCO/$output -m prot -l apicomplexa; done
```

The orthologs in BUSCO are classified as:
```shell
cat 12_BUSCO/Tg/run_apicomplexa_odb10/full_table.tsv | grep -v '^#' | cut -f2 | sort -u
Complete
Duplicated
Fragmented
Missing
```
We are therefore interested in creating the trees from orthologs which are duplicated or complete. 

#### Question 7

BUSCO labels the orhtologues found as Complete, Duplicated, Missing or Fragmented. To assess the quality, we can check how many Duplicated and/or Complete hits we have. 

We can first look at both labels:

```shell
for folder in 12_BUSCO/*; do file=$folder$'/run_apicomplexa_odb10/full_table.tsv'; species=$(echo $file | cut -d '/' -f2) ; echo $species; cd_busco=$(grep -v '^#' $file | cut -f1,2 | grep -E 'Duplicated|Complete' | cut -f1 | sort -u| wc -l ) ; echo $cd_busco;  tot=$(cat $file | grep -v '^#' $file | cut -f1 | sort -u| wc -l ); (echo  $cd_busco*100/$tot | bc -l ) ; done
Hc
325
72.86995515695067264573
Pb
372
83.40807174887892376681
Pc
429
96.18834080717488789237
Pf
436
97.75784753363228699551
Pk
323
72.42152466367713004484
Pv
437
97.98206278026905829596
Py
434
97.30941704035874439461
Tg
384
86.09865470852017937219
```

Then we look at the Complete hits alone:

```shell
for folder in 12_BUSCO/*; do file=$folder$'/run_apicomplexa_odb10/full_table.tsv'; species=$(echo $file | cut -d '/' -f2) ; echo $species; cd_busco=$(grep -v '^#' $file | cut -f1,2 | grep -E 'Complete' | cut -f1 | sort -u| wc -l ) ; echo $cd_busco;  tot=$(cat $file | grep -v '^#' $file | cut -f1 | sort -u| wc -l ); (echo  $cd_busco*100/$tot | bc -l ) ; done
Hc
323
72.42152466367713004484
Pb
361
80.94170403587443946188
Pc
429
96.18834080717488789237
Pf
435
97.53363228699551569506
Pk
322
72.19730941704035874439
Pv
435
97.53363228699551569506
Py
433
97.08520179372197309417
Tg
4
.89686098654708520179
```
Except for *Toxoplasma gondii*, the genomes seem to be of good quality. 

#### Question 9-10

```shell
mkdir -p 13_UniqProtID

# creates files to contain the IDs from BUSCO
for folder in 12_BUSCO/*; do file=$folder$'/run_apicomplexa_odb10/full_table.tsv'; species=$(echo $file | cut -d '/' -f2) ; echo $species; (grep -v '^#' $file | cut -f1,2 | grep -E 'Duplicated|Complete' | cut -f1 | sort -u > 13_UniqProtID/$species.txt) ; done

# makes two concatenated files, one with all species, and one without Tg
cat 13_UniqProtID/*.txt >> 13_UniqProtID/all.txt
cat 13_UniqProtID/{Hc,Pb,Pc,Pf,Pk,Pv,Py}.txt >> 13_UniqProtID/noTg.txt

# prints how many BUSCO ids are found in all species
cat 13_UniqProtID/all.txt | sort | uniq -c | grep '8 ' | wc -l
# 185
# prints how many BUSCO ids are found in all species except Tg
cat 13_UniqProtID/noTg.txt | sort | uniq -c | grep '7 ' | wc -l
# 207
```

Interpretation is found in the main.txt file. 

#### BUSCO Parser
Now that we have the results from BUSCO, we want to retrieve the orthologs sequences so that we can run first an alignment, and then build the phylogenetic trees. 

The python script I wrote, buscoParser.py, takes an input the species name (corresponding to the folder in the BUSCO result directory), some input directories and an optional output directory . 

Input files needed for the script:
```shell
# make a list of the BUSCO species folders
for i in 12_BUSCO/*; do sp=$(echo $i | cut -d '/' -f2); echo $sp; done > 13_UniqProtID/species_names.txt

mkdir 14_BUSCOParse_output
cat 13_UniqProtID/all.txt | sort | uniq -c | grep '8 ' | cut -d ' ' -f8 > 13_UniqProtID/BUSCO_uniq_ID.txt 
```
The python script runs as follows (see repo for full code):

```shell
python Scripts/buscoParser.py -busco_id 13_UniqProtID/BUSCO_uniq_ID.txt -species 13_UniqProtID/species_names.txt -busco_path 12_BUSCO/ -faa_path 10_allGeneMark/ -output_path 14_BUSCOParse_output/
```

This script produces 185 files, containing the species name and the ortholog sequence in a FASTA-like format. 

### Alignments and trees

First we create a conda environment:
```shell
conda create -n malariaalign
conda activate malariaalign
conda install -c bioconda clustalo raxml
# Clustal Omega - 1.2.4
# RAxML version 8.2.12
```
Then we run clustalo on all output files from the parser:
```shell
mkdir 15_clustalo
for file in 14_BUSCOParser_output/*.output; do id=$(echo $file | cut -d '/' -f2 | cut -d '.' -f1); output='15_clustalo/'$id'_aligned.faa'; clustalo -i $file -o $output -v ; done
```

And then we run raxml on clustalo outputs:
```shell

```
