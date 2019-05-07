# GateID
Cell type purification by single-cell transcriptome-trained sorting

File description:

# Functions_gate_design.R: functions to source to run the code.

# Code_gate_design.R contains the code to perform GateID gate design. 

To test the code, you can download the following test dataset: "gateID_training_dataset_test.csv" 
The code designs gates for cluster 1.

The following files are the gate solutions for the test dataset if the code ran correctly:
# gateID_solutions_test.csv : this csv file is the equivalent of the gate_sol dataframe that should be obtain after running the code.
# gateID_solutions_test.pdf: are the plots of the gate solutions if the code ran correctly. They are rank by decreasing purity (first page of the pdf= higher purity solution).

You can input your own dataset instead of the test dataset. Your dataset should have:
# single cells as rows
#the first column with a cell type identification: ideally a cluster number
#all subsequent column are FACS index values for all available channels

Happy gate design!
