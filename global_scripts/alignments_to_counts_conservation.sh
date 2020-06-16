
home="/home/fcalvet/Desktop/uORFs";

conservation_home="${home}/REFERENCE_DATA/conservation";
riboseq_software="${home}/REFERENCE_DATA/riboseq/software";

#### download the alignments data

wget -P $conservation_home http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.bw

$riboseq_software/bigWigToBedGraph ${conservation_home}/hg38.phastCons100way.bw ${conservation_home}/hg38.phastCons100way.bg;

rm ${conservation_home}/hg38.phastCons100way.bw;




wget -P $conservation_home http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw

$riboseq_software/bigWigToBedGraph ${conservation_home}/hg38.phyloP100way.bw ${conservation_home}/hg38.phyloP100way.bg;

rm ${conservation_home}/hg38.phyloP100way.bw;

