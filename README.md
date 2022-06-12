# Epigenomics_task4_5

Please, find the commands used to complete the task below. You can also find python and R script that were used in task5 in this repository. I also provide archives with regulatory_elements files. 

#entering docker

sudo su

docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course

# TASK4

mkdir ATAC-seq_new
cd ./ATAC-seq_new/
mkdir ./analyses
mkdir ./data
mkdir ./annotation
mkdir ./atac-nf

#downloading metadata file 

../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&assay_title=ATAC-seq&biosample_ontology.term_name=stomach&biosample_ontology.term_name=sigmoid+colon&type=Experiment"

#create directory for fastq files

mkdir ./data/fastq.files

#making directories for bed/bigwig files, inside ATAC-seq_new directory

mkdir data/bigBed.files data/bigWig.files

#getting bigbed peak files

grep -F ATAC-seq metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.peaks.ids.txt

cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done


##check the integrity of bed files
#retrieve original MD5 hash from the metadata
#compute MD5 hash on the downloaded files 
#make sure there are no files for which original and computed MD5 hashes differ

for file_type in bigBed; do
  ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,46 > data/"$file_type".files/md5sum.txt
  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt 
  awk '$2!=$3' data/"$file_type".files/md5sum.txt

done

#no different md5 were found, everything is OK, proceed

##getting gencode annotations

wget -P annotation "https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz"

gunzip annotation/gencode.v24.primary_assembly.annotation.gtf.gz

awk '$3=="gene"' annotation/gencode.v24.primary_assembly.annotation.gtf |\
grep -F "protein_coding" |\
cut -d ";" -f1 |\
awk 'BEGIN{OFS="\t"}{print $1, $4, $5, $10, 0, $7, $10}' |\
sed 's/\"//g' |\
awk 'BEGIN{FS=OFS="\t"}$1!="chrM"{$2=($2-1); print $0}' > annotation/gencode.v24.protein.coding.gene.body.bed

##retrieving peaks 

mkdir data/bed.files

cut -f1 analyses/bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done

#downloading list of promoters

cd ./annotation
wget "https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/gencode.v24.protein.coding.non.redundant.TSS.bed"
cd ..

#get peaks intersecting promoter regions

mkdir analyses/peaks.analysis/

cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a annotation/gencode.v24.protein.coding.non.redundant.TSS.bed -b data/bed.files/"$filename".bed -u |\
  cut -f7 |\
  sort -u > analyses/peaks.analysis/genes.with.peaks."$tissue".ATAC-seq.txt
done

#get peaks NOT intersecting gene body
#using flag -v to find regions that do not overlap


cut -f-2 analyses/bigBed.peaks.ids.txt |\
while read filename tissue; do 
  bedtools intersect -a annotation/gencode.v24.protein.coding.gene.body.bed -b data/bed.files/"$filename".bed -v |\
  cut -f7 |\
  sort -u > analyses/peaks.analysis/genes.no.peaks.gene.body."$tissue".ATAC-seq.txt
done


# TASK5

#Task1

#creating directories

cd ..

mkdir ./regulatory_elements/


#Task2
#copying all files we will need in regulatory_elements
cp ./ATAC-seq_new/analyses/peaks.analysis/genes.no.peaks.gene.body.stomach.ATAC-seq.txt ./regulatory_elements/
cp ./ATAC-seq_new/analyses/peaks.analysis/genes.no.peaks.gene.body.sigmoid_colon.ATAC-seq.txt ./regulatory_elements/
cp ./ATAC-seq_new/annotation/gencode.v24.protein.coding.gene.body.bed ./regulatory_elements/

##inside ChIP-seq folder

cd ./ChIP-seq/

#make directories where we will store files for H3K27ac and H3K4me1

mkdir data/acbigBed.files
mkdir data/me1bigBed.files

###retrieveing H3K27ac "ac" bigbed

grep -F H3K27ac metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/acbigBed.peaks.ids.txt

cut -f1 analyses/acbigBed.peaks.ids.txt |\
while read filename; do
wget -P data/acbigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done

#How many are they?

wc -l ./data/acbigBed.files/*

#8668 ./data/acbigBed.files/ENCFF872UHN.bigBed
#9162 ./data/acbigBed.files/ENCFF977LBD.bigBed
#17830 total



###retrieveing H3K4me1 "me1" bigbed

grep -F H3K4me1 metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $11, $23}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/me1bigBed.peaks.ids.txt


cut -f1 analyses/me1bigBed.peaks.ids.txt |\
while read filename; do
wget -P data/me1bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
done


#How many are they?

wc -l ./data/me1bigBed.files/*

#15883 ./data/me1bigBed.files/ENCFF724ZOF.bigBed
#12362 ./data/me1bigBed.files/ENCFF844XRN.bigBed
#84735 total



##converting to bed files 

cut -f1 analyses/acbigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/acbigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done
cut -f1 analyses/me1bigBed.peaks.ids.txt |\
while read filename; do
  bigBedToBed data/me1bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
done

#copy files to regulatory_elements folder
cd ..
cp ./ChIP-seq/data/bed.files/*.bed ./regulatory_elements/


#move to regulatory_elements folder

cd ../regulatory_elements

#####prepare bed files for ATAC-seq peaks outside gene body

for tissue in stomach sigmoid_colon; do
  ../bin/selectRows.sh genes.no.peaks.gene.body."$tissue".ATAC-seq.txt <(awk 'BEGIN{FS=OFS="\t"}{print $4, $0}' gencode.v24.protein.coding.gene.body.bed) |\
  cut -f2- > genes.no.peaks.gene.body."$tissue".ATAC-seq.bed
done



#######now find intersections for each tissue between ac, me1 and ATAC-seq annotations

#bedtools intersect -a ./genes.no.peaks.gene.body.stomach.ATAC-seq.bed -b ./ENCFF977LBD.bed ./ENCFF844XRN.bed
bedtools intersect -u -a ./genes.no.peaks.gene.body.stomach.ATAC-seq.bed -b ./ENCFF977LBD.bed |sort -u > stomach_reg1.bed
bedtools intersect -u -a ./stomach_reg1.bed -b ./ENCFF844XRN.bed |sort -u > stomach_reg.bed


bedtools intersect -u -a ./genes.no.peaks.gene.body.sigmoid_colon.ATAC-seq.bed -b ./ENCFF872UHN.bed > sigmoid_reg1.bed
bedtools intersect -u -a sigmoid_reg1.bed -b ./ENCFF724ZOF.bed |sort -u > sigmoid_colon_reg.bed


#Task3

#subsetting to the elements located only on chromosome 1

for tissue in stomach sigmoid_colon; do
 awk '$1="chr1"' "$tissue"_reg.bed > ./"$tissue"_reg_chr1.bed 
 awk '{if ($6=="+"){start=$2} else {start=$3}; print $4, start}' FS=' ' OFS='\t' "$tissue"_reg_chr1.bed > ./"$tissue"_reg_chr1_start.tsv ;
done


#Task4

awk '$1="chr1"' ./gencode.v24.protein.coding.gene.body.bed > ./gene.body.chr1.bed

awk '{if ($6=="+"){start=$2} else {start=$3}; print $4, start}' FS=' ' OFS='\t' ./gene.body.chr1.bed > ./task4.gene.starts.tsv 


#Task5 

#run nano and insert there a code provided, modify it accordingly, and save file as get.distance.py in ../bin/.

nano 

####create, modify, save the script. The final version of the script can be found in this repository.

python ./bin/get.distance.py --input ./regulatory_elements/task4.gene.starts.tsv --start 980000


#Task6

cat ./stomach_reg_chr1_start.tsv | while read element start; do  python ../bin/get.distance.py --input ./task4.gene.starts.tsv --start $start  | cut -f1,3 ;  done > ./stomach_regulatoryElements.genes.distances.tsv

cat ./sigmoid_colon_reg_chr1_start.tsv | while read element start; do  python ../bin/get.distance.py --input ./task4.gene.starts.tsv --start $start  | cut -f1,3 ;  done > ./sigmoid_colon_regulatoryElements.genes.distances.tsv


Task7

#inside regulatory_elements directory
#a script created in named R_epigenomics.R - run nano and write there an R code. the script can be found in this repository.

nano

Rscript ./R_epigenomics.R > task7mean_median_out.txt

#looking at the output

cat task7mean_median_out.txt

#[1] "the median distance between the gene and the regulatory element in sigmoid colon is: 2020 bp; the mean distance is: 5457.99481865285 bp"
#[1] "the median distance between the gene and the regulatory element in stomach is: 2168 bp; the mean distance is: 6731.58851674641 bp"

