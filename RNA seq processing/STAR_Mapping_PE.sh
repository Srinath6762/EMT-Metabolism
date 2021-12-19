## Mapping by STAR-aligner of paired end reads

#ulimit -n 999999
#ulimit -f unlimited

echo "Mapping Started"
date

for i in $(ls *.fastq.gz | sed -e 's/_1.fastq.gz//' -e 's/_2.fastq.gz//' | sort -u)

do

	STAR \
	--runThreadN 20\
	--genomeDir /media/csb/'New Volume'/susmita/hg38_STARIndexed\
	--readFilesCommand zcat\
	--readFilesIn ${i}_1.fastq.gz ${i}_2.fastq.gz\
	--outSAMattributes All\
	--outSAMtype BAM Unsorted\
	--quantMode GeneCounts\
	--outFileNamePrefix ./BAM/${i}_

done

date
echo "Mapping Done"
