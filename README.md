# FM_Partitioning
Program which implements Fiduccia–Mattheyses partitioning algorithm
PROGRAMMING ASSIGNMENT PA3
Mark Sears and Navya Sreeram
10/15/2018

==================================================
PROGRAM EXECUTION INSTRUCTIONS

The code assumes that the example design data can be found in a directory "FPGA Benchmarks" (this is set in source code by string "file_path").

Compiling the code can be done with the included makefile. 

To run the code, type in the terminal: "./fmpart.out #" where "#" is the example design number you wish to use. The code will look for the design in "./FPGA Benchmarks/FPGA-example#"

==================================================
EXAMPLES

In terminal, navigate to a working directory which includes the source code "fmpart.cpp", the makefile, and a directory "FPGA Benchmarks".

Type "make" to compile.

Type "./fmpart.out 1" to run the partitioning on example 1.
Type "./fmpart.out 2" to run the partitioning on example 2.


==================================================
DESCRIPTION OF ALGORITHM

The Fidducia-Mattheyses method of bi-partitioning is implemented. Classes for cells and nets are used to organize the design information. The cells and nets files are first read in and stored to memory. As they are processed, different data structures are created to better organize the information for quick reference and manipulation.

	1) All cells and nets are stored in a vector. This is the primary method of accessing all cells in an arbitrary order.	
	2) Cells are stored in a "gain bucket" multimap with the gains used as keys. This allows the algorithm to quickly pick a cell which has the highest gain.

The general program flow can be seen in the main function:

	1) The design files are read.
	2) Initial partitioning. Cells are randomly assigned to partitions until the area in each is about 50%.
	3) Initial cutsize. All nets are looped through to see if it is cut by the partitioning.
	4) Cell gains are computed.
	5) Perform FM passes until a local minimum is found that cannot be "climbed" out of.

Each FM pass consists of moving cells between partitions until either all cells are locked or no legal moves remain (defined by the area constraint). Cells with highest gain are selected first.

The algorithm continues performing passes until the result of a pass does not yield any cutsize improvement. In practice, this usually occurs after about 15 passes. A maximum pass count of 20 was also added to put a cap on execution time.

==================================================
RESULTS

Outputs of running the algorithm are tabulated in "FM partition results.xls". The 1st tab shows results of running passes until a local minimum is found (i.e. a pass does not yield an improvement in cutsize). The 2nd tab shows results of stopping after a specified number of passes.

During a run, the current pass and cutsize are printed to the console. After finishing, the results of the partitioning, including a list of which partition each cell should be in, is printed to a file "fmpartition.nodes" in the example's design directory.

Overall, results appear good. The percent change from initial random partitioning to final is between 90-95% for most runs. 

Execution time was kept low by ensuring the code does not do any unnecessary work. For example, the "best_partition" is not updated until all non-negative gains are exhausted during a pass. This helps cut down on unnecessary operations without impacting the final results.

The Ratio Cut is a measure of partition balance and should be minimized. Rc = cutsize / |A|*|B|
This was accomplished primarily because the next cell to move was selected from alternating partitions, so on average the number of cells in each partition was not changed very much.


