Description of elements.
================================

[]Data: [0]
    * Will have the genome data, analysis data, CSV, and such.
 --- diff_formats:
        * At the moment it has the GFF format.
     --- to_gff.ipynb: to transform the CSV file from my project to GFF

[]Compare_run20_input_chr1_vs_input_chr32: [1]
    * Scripts made to compare the data from the SOFTWARE with the input from chr1 and the input from chr20

[]Manual analysis: [2]
 --- Created to compare the original TRUE POSITIVES data vs my owm.

[]BEDOPS_Test_Final: [3]
    * Created to experiment and check how BEDOPS works.
    * Code results implemented in the main Leishmania Project.

[]BEDOPS_join_strands: [4]
    * Problem: the 6k elements are in both strands and sometimes are overlaped.
    * Need to remove the overlaps in that way.
    * Creation: 15/04/24

[]Compare_OrigPositvSIDER_vs_LastFiltered: [5]
    * Let's compare the output from BEDOPS_join_strands vs the SIDER positive values.

[]Subfamilies_test_1: [6]
    * Try to make a subfamily classification with elements >85% and alignment length of > 100 nt.

[]Filter_Test_1
    * Created to filter the RUN20.csv data with different methods.
 --- filter_test_1.py ==> filtered data by 1.0E-09
 --- filtering_draft1.ipynb ===> test jupyter notebook to define `filter_test_1.py`


