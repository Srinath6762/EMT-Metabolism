## Mapping by STAR-aligner of single end reads

#ulimit -n 999999
#ulimit -f unlimited

echo "Mapping Started"
date

for i in $(ls *.fastq.gz | sed -e 's/.fastq.gz//' | sort -u)

do

	STAR \
	--runThreadN 20\
	--genomeDir /media/csb/'New Volume'/susmita/hg38_STARIndexed\
	--readFilesCommand zcat\
	--readFilesIn ${i}.fastq.gz\
	--outSAMattributes All\
	--outSAMtype BAM Unsorted\
	--quantMode GeneCounts\
	--outFileNamePrefix ./BAM/${i}_

done

date
echo "Mapping Done"
