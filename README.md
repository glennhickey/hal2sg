# hal2sg
Prototype code for converting HAL to Side Graph SQL

## algorithm

Iteratatively add genomes to side graph.  Each genome is aligned to nearest genome already in side graph.  This alignment is used to thread the genome onto the graph.

### CAMEL Input

The following logic is used to support CAMEL output:  If the root DNA sequence is all N's, then it will be inferred from the sequences of its children. 

## instructions

**Dependencies:**  HAL and SonLib.  Expected to be in sister directories to hal2sg but can be changed in include.mk.  Tested only with latest version from development branch for both. 

**Note** Best to compile in debug mode (cppflags_dbg in sonLib/include.mk).  asserts are used to check if the output is honest.  At the very least run the unit tests:

	  make test

To run the converter:

	  hal2sg input.hal output.fa output.sql > output.out

`input.hal` Input alignment to convert

`output.fa` Output fasta file of all Side Graph sequences

`output.sql` Output text file listing INSERT commands for Sequences, Joins and Paths (for each input sequence) in the graph.

## sorry state of affairs (Update May 13)

* **NOTE:** The only working option is `--rootGenome`.  This means that hal2sg will only process complete clades (including ancestors) in a top-down manner. 

## next updates

* Support bottom-up mappings in order to activate remaining genome selection options.  




