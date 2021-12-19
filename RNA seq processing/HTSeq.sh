#Sorting the BAM file by name (htseq default: name)

date
echo "Sorting Started"

for i in $(ls *.bam | sed -e 's/_Aligned.out.bam//' | sort -u)

do

	samtools sort\
	-@ 20\
	-n\
	${i}_Aligned.out.bam\
	-o ${i}.bam

	rm -rf ${i}_Aligned.out.bam

done


echo "Sorting Done"
date

# Change directory & reference genome accordingly!
echo "HTSeq Started"

date
parallel -j 20 'htseq-count --additional-attr=gene_name -s reverse -f bam {} /media/csb/New\ Volume/susmita/hg38/Homo_sapiens.GRCh38.99.gtf > /media/csb/Backup\ Plus/tanishq/GSE146532/htseq/{}_htseq.txt' ::: *.bam
date

echo "HTSeq Done"
