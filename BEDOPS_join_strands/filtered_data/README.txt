[ ] blastn_blaster
    * It has the data of the blastn dictionary from the positive elements and the fasta for the ingi element.

[ ] blaster_neg
    * It has the data of the balstn dictionary from the negative elements and the fasta for the ingi element.

[ ] gff
    * It has the data of the positive elements in a GFF file format.

negatives_testing_elements.csv
    * elements negative to the filtering of at least in 5 chromosomes and evalue > 1E-09

positives_testing_elements.csv
    * elements negative to the filtering of at least in 5 chromosomes and evalue > 1E-09

to_fasta_neg.ipynb
    * jupyter notebook file to create a fasta and BLASTn analysis from the `negatives.fasta` elements and the `ingi.fasta`


to_fasta.ipynb
    * jupyter notebook file to create a fasta and BLASTn analysis from the `positives.fasta` elements and the `ingi.fasta`

to_gff.ipynb
    * jupyter notebook file to create a GFF file from the positive elements.