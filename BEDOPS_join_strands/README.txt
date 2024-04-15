[] filtered_data
    * Folder to place the results from the filtering methods.

 --- sequences.csv
        * Output from `extract_seqs`.ipynb

 --- negatives_testing_elements.csv
        * data filtered through `filter_test_1.py`
 --- positives_testing_elements.csv
        * data filtered through `filter_test_1.py`

data_bedops.gff ==>
    * Last output GFF from `testing_1.ipynb`

data_gff_sorted.bed ==>
    * GFF file used in `testing_1.ipynb`

extract_seqs.ipynb ==>
    * reads `data_bedops.gff`
    * Uses BLASTn local `blastdbcmd` to get sequences.
    * Prepares data frame.
    * Export data into `sequences.csv`

filter_test_1.ipynb ==>
    * jupyter notebook to test code for `filter_test_1`

filter_test_1 ==>
    * Script to filter elements by being in >= 5 chromosomes and e-value <= 1.0E-09

testing_1.ipynb ==>
    * Imported GFF data
    * created BED file
    * joined BED file with data_bedops
    * Created GFF file after that