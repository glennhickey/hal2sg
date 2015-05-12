# hal2sg
Prototype code for converting HAL to Side graph

## algorithm

Iteratatively add genomes to side graph.  Each genome is aligned to nearest genome already in side graph.  This alignment is used to thread the genome onto the graph.

## instructions

**Dependencies:**  HAL and SonLib.  Expected to be in sister directories to hal2sg but can be changed in include.mk.  Tested only with latest version from development branch for both. 

**Note** Best to compile in debug mode (cppflags_dbg in sonLib/include.mk).  asserts are used to check if the output is honest.  At the very least run the unit tests:

	  make test

To run the converter:

	  hal2sg input.hal output.fa output.sql > output.out

`input.hal` Input alignment to convert

`output.fa` Output fasta file of all Side Graph sequences

`output.sql` Output text file listing INSERT commands for Sequences, Joins and Paths (for each input sequence) in the graph.

* **None of the options except `--targetGenomes` (and the standard HAL caching options) work.**
* Please ignore voluminous garbage printed to stdout.

## sorry state of affairs (Update May 12)

`hal2sg` only works on input of this form:
* **Root reference:** root genome used as initial side graph sequence.
* **Star Tree** all other genomes children of root. 

## next updates

* Non-star trees.  Only missing a unit test
* Clean up debug output on stdout
* Non-root reference.   There is a paralogy case that is not yet implemented that can arise when mapping up then down tree.  





