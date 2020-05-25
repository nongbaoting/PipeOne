


## apa
qapa fasta -f $hg38_fa  qapa_3utrs.gencode_V31.hg38.bed qapa_3utrs.gencodeV31.hg38.fa
python3 /dsk2/who/nbt/python/binPy3/fasta.py shorten_id qapa_3utrs.gencodeV31.hg38.fa qapa_3utrs.gencodeV31.hg38.shortedID.fa
salmon index -t qapa_3utrs.gencodeV31.hg38.shortedID.fa -i 3utr -p 8