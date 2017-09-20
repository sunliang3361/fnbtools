1.Produce Bed File With Gaps In Coverage From A Bam File
https://www.biostars.org/p/76809/
2.Question: how to visualize the bedgraph and bed file
https://www.biostars.org/p/102929/



command lines:
1.extract all gaps
genomeCoverageBed -ibam accepted_hits.bam -bga | awk '$4 == 0' > intervals.with.0.cov.bedg


Dec 7th
1. 3156 and 3164 should merge 3 fastq files
or7t3d ---- or7       3156
22rtd7 ---- 22r       3164

2. 3196 and 6507 should merge 2 fastq files
r2dax3d                 3196
fgdetr66                6507


----------
Data analysis for Jiangqi

1.FNBtool6.pl
nohup perl FNBscan6.pl -g mt4.fa -1 S1-FN3156_S1_R1.fastq.gz,S2-FN3164_S2_R1.fastq.gz,S3-FN3196_S3_R1.fastq.gz,S4-FN6507_S4_R1.fastq.gz -2 S1-FN3156_S1_R2.fastq.gz,S2-FN3164_S2_R2.fastq.gz,S3-FN3196_S3_R2.fastq.gz,S4-FN6507_S4_R2.fastq.gz > log.out 2>&1
2.VarDiff.py
python VarDiff.py -m result/AF_S1-FN3156_S1_R1.fastq.gz.out -c result/AF_S2-FN3164_S2_R1.fastq.gz.out result/AF_S3-FN3196_S3_R1.fastq.gz.out result/AF_S4-FN6507_S4_R1.fastq.gz.out -o result/AF_unique_S1.txt
python VarDiff.py -m result/AF_S2-FN3164_S2_R1.fastq.gz.out -c result/AF_S1-FN3156_S1_R1.fastq.gz.out result/AF_S3-FN3196_S3_R1.fastq.gz.out result/AF_S4-FN6507_S4_R1.fastq.gz.out -o result/AF_unique_S2.txt
python VarDiff.py -m result/AF_S3-FN3196_S3_R1.fastq.gz.out -c result/AF_S1-FN3156_S1_R1.fastq.gz.out result/AF_S2-FN3164_S2_R1.fastq.gz.out result/AF_S4-FN6507_S4_R1.fastq.gz.out -o result/AF_unique_S3.txt
python VarDiff.py -m result/AF_S4-FN6507_S4_R1.fastq.gz.out -c result/AF_S1-FN3156_S1_R1.fastq.gz.out result/AF_S2-FN3164_S2_R1.fastq.gz.out result/AF_S3-FN3196_S3_R1.fastq.gz.out -o result/AF_unique_S4.txt
 

3.python DelDiff.py -m result/fnb.S1-FN3156_S1_R1.fastq.gz_gap.bedg -c result/fnb.S2-FN3164_S2_R1.fastq.gz_gap.bedg result/fnb.S3-FN3196_S3_R1.fastq.gz_gap.bedg result/fnb.S4-FN6507_S4_R1.fastq.gz_gap.bedg -o result/fnb_big_deletion_S1.txt
python DelDiff.py -m result/fnb.S2-FN3164_S2_R1.fastq.gz_gap.bedg -c result/fnb.S1-FN3156_S1_R1.fastq.gz_gap.bedg result/fnb.S3-FN3196_S3_R1.fastq.gz_gap.bedg result/fnb.S4-FN6507_S4_R1.fastq.gz_gap.bedg -o result/fnb_big_deletion_S2.txt
python DelDiff.py -m result/fnb.S3-FN3196_S3_R1.fastq.gz_gap.bedg -c result/fnb.S1-FN3156_S1_R1.fastq.gz_gap.bedg result/fnb.S2-FN3164_S2_R1.fastq.gz_gap.bedg result/fnb.S4-FN6507_S4_R1.fastq.gz_gap.bedg -o result/fnb_big_deletion_S3.txt
 python DelDiff.py -m result/fnb.S4-FN6507_S4_R1.fastq.gz_gap.bedg -c result/fnb.S1-FN3156_S1_R1.fastq.gz_gap.bedg result/fnb.S2-FN3164_S2_R1.fastq.gz_gap.bedg result/fnb.S3-FN3196_S3_R1.fastq.gz_gap.bedg -o result/fnb_big_deletion_S4.txt



4. python VarAnnot.py -i result/fnb_big_deletion_S1.txt -f Mtruncatula_285_Mt4.0v1.gene.gff3 -o result/fnb_big_deletion_S1_annot.txt
python VarAnnot.py -i result/fnb_big_deletion_S2.txt -f Mtruncatula_285_Mt4.0v1.gene.gff3 -o result/fnb_big_deletion_S2_annot.txt
python VarAnnot.py -i result/fnb_big_deletion_S3.txt -f Mtruncatula_285_Mt4.0v1.gene.gff3 -o result/fnb_big_deletion_S3_annot.txt
python VarAnnot.py -i result/fnb_big_deletion_S4.txt -f Mtruncatula_285_Mt4.0v1.gene.gff3 -o result/fnb_big_deletion_S4_annot.txt

python VarAnnot.py -i result/AF_unique_S2.txt -f Mtruncatula_285_Mt4.0v1.gene.gff3 -o result/AF_unique_S2_annot.txt
python VarAnnot.py -i result/AF_unique_S3.txt -f Mtruncatula_285_Mt4.0v1.gene.gff3 -o result/AF_unique_S3_annot.txt
python VarAnnot.py -i result/AF_unique_S4.txt -f Mtruncatula_285_Mt4.0v1.gene.gff3 -o result/AF_unique_S4_annot.txt 

5. python CircosVis.py -s result/AF_unique_S1.txt  -l result/fnb_big_deletion_S1_annot.txt -o S1.png
python CircosVis.py -s result/AF_unique_S2.txt  -l result/fnb_big_deletion_S2_annot.txt -o S2.png
python CircosVis.py -s result/AF_unique_S3.txt  -l result/fnb_big_deletion_S3_annot.txt -o S3.png
python CircosVis.py -s result/AF_unique_S4.txt  -l result/fnb_big_deletion_S4_annot.txt -o S4.png

-------------
Filter bad mapping reads from bam file
-------------
http://seqanswers.com/forums/archive/index.php/t-64879.html
samtools view -q 30 inputfile.bam >good.bam


