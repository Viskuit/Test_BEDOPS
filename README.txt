Description of elements.
================================

BEDOPS_join_strands:
    [3]
    * Problem: the 6k elements are in both strands and sometimes are overlaped.
    * Need to remove the overlaps in that way.
    * Creation: 15/04/24

BEDOPS_Test_Final: 
    [2]
    * Created to experiment and check how BEDOPS works.
    * Code results implemented in the main Leishmania Project.

Compare_OrigPositvSIDER_vs_LastFiltered
    [4]
    * Let's compare the output from BEDOPS_join_strands vs the SIDER positive values.

Compare_run20_input_chr1_vs_input_chr32
    [1]
    * Scripts made to compare the data from the SOFTWARE with the input from chr1 and the input from chr20

Data:
    * Will have the genome data, analysis data, CSV, and such.
 --- diff_formats:
        * At the moment it has the GFF format.
     --- to_gff.ipynb: to transform the CSV file from my project to GFF

Filter_Test_1
    * Created to filter the RUN20.csv data with different methods.
 --- filter_test_1.py ==> filtered data by 1.0E-09
 --- filtering_draft1.ipynb ===> test jupyter notebook to define `filter_test_1.py`

