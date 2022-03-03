# BINP29_MalariaProject
Repo containing the scripts used to analyze malaria parasites' genome in the BINP29 course (Lund University)
https://github.com/aurazelco/BINP29_MalariaProject.git

All analysis is executed on the server.

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
## Processing of Haemoproteus tartakovskyi data
### Clean genome sequence

mkdir -p Scripts && cd Scripts
scp inf-36-2021@bioinf-biol302436.biol.lu.se:/home/inf-36-2021/Desktop/BINP29/malaria_project/removeScaffold.py .

# insert removeScaffold.py script here

mkdir -p 2_RemoveScaffolds && cd 2_RemoveScaffolds

python ../Scripts/removeScaffold.py -i ../Haemoproteus_tartakovskyi.genome -o Haemoproteus_tartakovskyi.output

mkdir -p 3_GenePredictionHaemoproteus && cd 3_GenePredictionHaemoproteus
nohup gmes_petap.pl --sequence ../2_RemoveScaffolds/Haemoproteus_tartakovskyi.output &

mv genemark_es.gtf Haemoproteus.gff
cat Haemoproteus.gff | sed "s/ GC=.*\tGeneMark.hmm/\tGeneMark.hmm/" > Ht2.gff
cd ..

mkdir -p 4_GTFtoFASTA && cd 4_GTFtoFASTA
gffParse.pl -v
#gffParse.pl version 1.1

gffParse.pl -i ../2_RemoveScaffolds/Haemoproteus_tartakovskyi.output  -g ../3_GenePredictionHaemoproteus/Ht2.gff -c -p

#INFO: Following files were output in the directory
#      /home/inf-36-2021/BINP29/malaria_project/4_GTFtoFASTA:
#      gffParse.fna (fasta file containing the genes)
#      gffParse.log (log file)
#      gffParse.faa (fasta file containing translated genes)


cd ..
mkdir -p 5_BLAST && cd 5_BLAST
blastx -query ../4_GTFtoFASTA/gffParse.fna -db SwissProt -out Haemoproteus.blastx -num_threads 20 -evalue 1e-10
# BLASTX 2.11.0+

ln -s /resources/binp29/Data/malaria/taxonomy.dat taxonomy.dat
ln -s /resources/binp29/Data/malaria/uniprot_sprot.dat uniprot_sprot.dat


scp inf-36-2021@bioinf-biol302436.biol.lu.se:/home/inf-36-2021/Desktop/BINP29/malaria_project/removeScaffold.py ../Scripts

# just to have an input to play with for the second python script
# cp /resources/binp29/Data/malaria/backup_results/Ht.blastout.gz .
# gzip -d Ht.blastout.gz

# decided to use the given script, since it was quite difficult to visualize properly the files without a text editor on the server
python datParser.py Haemoproteus.blastx ../4_GTFtoFASTA/gffParse.fna taxonomy.dat uniprot_sprot.dat > bird_scaffolds.txt

# ignores the lines (-v) which match the given pattern (-f) from bird_scaffolds.txt and saves in a new file
grep -v -f bird_scaffolds.txt ../4_GTFtoFASTA/gffParse.fna > filtered_scaffolds.txt

cd ..

cd Scripts
scp inf-36-2021@bioinf-biol302436.biol.lu.se:/home/inf-36-2021/Desktop/BINP29/malaria_project/removeBird.py .
cd ..
mkdir -p 6_FilteredGeneMark && cd 6_FilteredGeneMark
python ../Scripts/removeBird.py -i ../2_RemoveScaffolds/Haemoproteus_tartakovskyi.output -s ../5_BLAST/bird_scaffolds.txt -o Haemoproteus_clean.genome

nohup gmes_petap.pl --sequence Haemoproteus_clean.genome &

# length of genome
cd 6_FilteredGeneMark/
mv genemark_es.gtf Haemoproteus_clean.gtf
# now we run the gffParse.pl again, so we extract the genes again
cd ..
mkdir -p 7_CleanGTFtoFASTA && cd 7_CleanGTFtoFASTA
gffParse.pl -i ../6_FilteredGeneMark/Haemoproteus_clean.genome -g ../6_FilteredGeneMark/Haemoproteus_clean.gtf -c -p
cd ../6_FilteredGeneMark
cat Haemoproteus_clean.genome | grep -v "^>" | tr -d "\n" | wc -c
# 22068498
cd ..
for file in *.genome; do echo $file; grep -v "^>" $file | tr -d "\n" | wc -c; done
#Haemoproteus_tartakovskyi.genome -> this is before filtering
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

mkdir -p 8_allGTFs && cd 8_allGTFs
cp  /tmp/Prediction/*.gtf .
cp /resources/binp29/Data/malaria/Tg.gff.gz .
gzip -d Tg.gff.gz
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

# Gc content
mkdir -p 9_Genomes
mv *.genome 9_Genomes/
rm 9_Genomes/Haemoproteus_tartakovskyi.genome
cp 6_FilteredGeneMark/Haemoproteus_clean.genome 9_Genomes/
scp inf-36-2021@bioinf-biol302436.biol.lu.se:/home/inf-36-2021/Desktop/BINP29/malaria_project/GCcontent.py ../Scripts
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



#### Phylogenetic trees
mkdir 10_allGeneMark
mv 8_allGTFs/knowlesi.gtf 8_allGTFs/plasmodium_knowlesi.gtf

for file in 8_allGTFs/*; do echo $file; genome=$(echo $file | cut -d '/' -f2 | cut -d '.' -f1); genus=$(echo $file | cut -d '/' -f2 | cut -c 1); species=$(echo $file | cut -d '_' -f3 | cut -c 1); output=(${genus^}$species); gffParse.pl -i 9_Genomes/${genome^}'.genome' -g $file -c -p -b 10_allGeneMark/$output ; done

#------------------------------------------------
8_allGTFs/Haemoproteus_clean.gtf

INFO: The gff or gtf file 8_allGTFs/Haemoproteus_clean.gtf has successfully been parsed.
      There were 1734 scaffolds containing genes.
      The scaffolds contained 4429 genes.
      The genes contained the feature CDS 11953 times.

INFO: The scaffold/genome file 9_Genomes/Haemoproteus_clean.genome was successfully parsed.
      There were 1803 scaffolds. This number may be higher than the one above.

INFO: Following files were output in the directory
      /home/inf-36-2021/BINP29/malaria_project:
      10_allGeneMark/Hc.fna (fasta file containing the genes)
      10_allGeneMark/Hc.log (log file)
      10_allGeneMark/Hc.faa (fasta file containing translated genes)

NOTICE: There are 951 warnings in the log file.
        Reading frames have been adjusted for 951 genes (with the use of -c option).

8_allGTFs/plasmodium_berghei.gtf

INFO: The gff or gtf file 8_allGTFs/plasmodium_berghei.gtf has successfully been parsed.
      There were 5271 scaffolds containing genes.
      The scaffolds contained 7282 genes.
      The genes contained the feature CDS 16453 times.

INFO: The scaffold/genome file 9_Genomes/Plasmodium_berghei.genome was successfully parsed.
      There were 7479 scaffolds. This number may be higher than the one above.

INFO: Following files were output in the directory
      /home/inf-36-2021/BINP29/malaria_project:
      10_allGeneMark/Pb.fna (fasta file containing the genes)
      10_allGeneMark/Pb.log (log file)
      10_allGeneMark/Pb.faa (fasta file containing translated genes)

NOTICE: There are 3087 warnings in the log file.
        Reading frames have been adjusted for 3087 genes (with the use of -c option).

8_allGTFs/plasmodium_cynomolgi.gtf

INFO: The gff or gtf file 8_allGTFs/plasmodium_cynomolgi.gtf has successfully been parsed.
      There were 858 scaffolds containing genes.
      The scaffolds contained 5787 genes.
      The genes contained the feature CDS 17161 times.

INFO: The scaffold/genome file 9_Genomes/Plasmodium_cynomolgi.genome was successfully parsed.
      There were 1663 scaffolds. This number may be higher than the one above.

INFO: Following files were output in the directory
      /home/inf-36-2021/BINP29/malaria_project:
      10_allGeneMark/Pc.fna (fasta file containing the genes)
      10_allGeneMark/Pc.log (log file)
      10_allGeneMark/Pc.faa (fasta file containing translated genes)

NOTICE: There are 259 warnings in the log file.
        Reading frames have been adjusted for 259 genes (with the use of -c option).

8_allGTFs/plasmodium_faciparum.gtf

INFO: The gff or gtf file 8_allGTFs/plasmodium_faciparum.gtf has successfully been parsed.
      There were 15 scaffolds containing genes.
      The scaffolds contained 5207 genes.
      The genes contained the feature CDS 15638 times.

INFO: The scaffold/genome file 9_Genomes/Plasmodium_faciparum.genome was successfully parsed.
      There were 15 scaffolds. This number may be higher than the one above.

INFO: Following files were output in the directory
      /home/inf-36-2021/BINP29/malaria_project:
      10_allGeneMark/Pf.fna (fasta file containing the genes)
      10_allGeneMark/Pf.log (log file)
      10_allGeneMark/Pf.faa (fasta file containing translated genes)

NOTICE: There are 2 warnings in the log file.
        Reading frames have been adjusted for 2 genes (with the use of -c option).

8_allGTFs/plasmodium_knowlesi.gtf

INFO: The gff or gtf file 8_allGTFs/plasmodium_knowlesi.gtf has successfully been parsed.
      There were 14 scaffolds containing genes.
      The scaffolds contained 4953 genes.
      The genes contained the feature CDS 15416 times.

INFO: The scaffold/genome file 9_Genomes/Plasmodium_knowlesi.genome was successfully parsed.
      There were 14 scaffolds. This number may be higher than the one above.

INFO: Following files were output in the directory
      /home/inf-36-2021/BINP29/malaria_project:
      10_allGeneMark/Pk.fna (fasta file containing the genes)
      10_allGeneMark/Pk.log (log file)
      10_allGeneMark/Pk.faa (fasta file containing translated genes)

NOTICE: There are 4024 warnings in the log file.
        Reading frames have been adjusted for 858 genes (with the use of -c option).

        8_allGTFs/plasmodium_vivax.gtf

        INFO: The gff or gtf file 8_allGTFs/plasmodium_vivax.gtf has successfully been parsed.
              There were 786 scaffolds containing genes.
              The scaffolds contained 5682 genes.
              The genes contained the feature CDS 15461 times.

        INFO: The scaffold/genome file 9_Genomes/Plasmodium_vivax.genome was successfully parsed.
              There were 2747 scaffolds. This number may be higher than the one above.

        INFO: Following files were output in the directory
              /home/inf-36-2021/BINP29/malaria_project:
              10_allGeneMark/Pv.fna (fasta file containing the genes)
              10_allGeneMark/Pv.log (log file)
              10_allGeneMark/Pv.faa (fasta file containing translated genes)

        NOTICE: There are 238 warnings in the log file.
                Reading frames have been adjusted for 238 genes (with the use of -c option).

        8_allGTFs/plasmodium_yoelii.gtf

        INFO: The gff or gtf file 8_allGTFs/plasmodium_yoelii.gtf has successfully been parsed.
              There were 116 scaffolds containing genes.
              The scaffolds contained 4919 genes.
              The genes contained the feature CDS 14724 times.

        INFO: The scaffold/genome file 9_Genomes/Plasmodium_yoelii.genome was successfully parsed.
              There were 130 scaffolds. This number may be higher than the one above.

        INFO: Following files were output in the directory
              /home/inf-36-2021/BINP29/malaria_project:
              10_allGeneMark/Py.fna (fasta file containing the genes)
              10_allGeneMark/Py.log (log file)
              10_allGeneMark/Py.faa (fasta file containing translated genes)

        NOTICE: There are 18 warnings in the log file.
                Reading frames have been adjusted for 18 genes (with the use of -c option).

                8_allGTFs/toxoplasma_gondii.gtf

                INFO: The gff or gtf file 8_allGTFs/toxoplasma_gondii.gtf has successfully been parsed.
                      There were 2238 scaffolds containing genes.
                      The scaffolds contained 15892 genes.
                      The genes contained the feature CDS 103760 times.

                INFO: The scaffold/genome file 9_Genomes/Toxoplasma_gondii.genome was successfully parsed.
                      There were 2290 scaffolds. This number may be higher than the one above.

                INFO: Following files were output in the directory
                      /home/inf-36-2021/BINP29/malaria_project:
                      10_allGeneMark/Tg.fna (fasta file containing the genes)
                      10_allGeneMark/Tg.log (log file)
                      10_allGeneMark/Tg.faa (fasta file containing translated genes)

#------------------------------------------------

# installation
conda create -n malariaenv
conda activate malariaenv
conda install proteinortho
conda update --name base conda
conda deactivate
conda create -n malariaBusco -c conda-forge -c bioconda busco=5.3.0


# proteinortho
mkdir 11_ProteinOrtho && cd 11_ProteinOrtho
conda activate malariaenv
# I have Hc instead of Ht because the names are Haemo..._clean
nohup proteinortho6.pl ../10_allGeneMark/{Hc,Pb,Pc,Pf,Pk,Pv,Py,Tg}.faa &
conda deactivate

grep -v '^#'  11_ProteinOrtho/myproject.proteinortho.tsv | wc -l
# 5126


#BUSCO
conda malariaBusco
mkdir -p 12_BUSCO

for file in 10_allGeneMark/*.faa; do echo $file; output=$(echo $file | cut -d '/' -f2 | cut -d '.' -f1); echo $output; busco -i $file -o 12_BUSCO/$output -m prot -l apicomplexa; done



cat 12_BUSCO/Tg/run_apicomplexa_odb10/full_table.tsv | grep -v '^#' | cut -f2 | sort -u
Complete
Duplicated
Fragmented
Missing

for folder in 12_BUSCO/*; do file=$folder$'/run_apicomplexa_odb10/full_table.tsv'; species=$(echo $file | cut -d '/' -f2) ; echo $species; (grep -v '^#' $file | cut -f1,2 | grep -E 'Duplicated|Complete' | cut -f1 | sort -u| wc -l ) ; done
Hc
325
Pb
372
Pc
429
Pf
436
Pk
323
Pv
437
Py
434
Tg
384

mkdir -p 13_UniqProtID

# creates files to contain the IDs from BUSCO
for folder in 12_BUSCO/*; do file=$folder$'/run_apicomplexa_odb10/full_table.tsv'; species=$(echo $file | cut -d '/' -f2) ; echo $species; (grep -v '^#' $file | cut -f1,2 | grep -E 'Duplicated|Complete' | cut -f1 | sort -u > 13_UniqProtID/$species.txt) ; done

# makes one concatenated file
cat 13_UniqProtID/*.txt >> 13_UniqProtID/all.txt
cat 13_UniqProtID/all.txt | sort | uniq -c | grep '8 ' | wc -l
# 185
cat 13_UniqProtID/noTg.txt | sort | uniq -c | grep '7 ' | wc -l
# 207

# made a list of the BUSCO species folders
for i in 12_BUSCO/*; do sp=$(echo $i | cut -d '/' -f2); echo $sp; done > 13_UniqProtID/species_names.txt

# copies the script from local to server
scp inf-36-2021@bioinf-biol302436.biol.lu.se:/home/inf-36-2021/Desktop/BINP29/malaria_project/buscoParser.py Scripts/

mkdir 14_BUSCOParse_output
cat 13_UniqProtID/all.txt | sort | uniq -c | grep '8 ' | cut -d ' ' -f8 > 13_UniqProtID/BUSCO_uniq_ID.txt 
python Scripts/buscoParser.py -busco_id 13_UniqProtID/BUSCO_uniq_ID.txt -species 13_UniqProtID/species_names.txt -busco_path 12_BUSCO/ -faa_path 10_allGeneMark/ -output_path 14_BUSCOParse_output/

