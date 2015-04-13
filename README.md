# hal2sg
Prototype code for converting HAL to Side graph

## algorithm

Iteratatively add genomes to side graph.  Each genome is aligned to nearest genome already in side graph.  This alignment is used to thread the genome onto the graph.

## instructions

**Dependencies:**  HAL and SonLib.  Expected to be in sister directories to hal2sg but can be changed in include.mk.  Tested only with latest version from development branch for both. 

**Note** Best to compile in debug mode (cppflags_dbg in sonLib/include.mk).  asserts are used to check if the output is honest. 

	  hal2sg input.hal output.fa output.sql 

`input.hal` Input alignment to convert

`output.fa` Output fasta file of all Side Graph sequences

`output.sql` Output text file listing INSERT commands for Sequences, Joins and Paths (for each input sequence) in the graph.

* **None of the options except `--targetGenomes` (and the standard HAL caching options) work.**
* Please ignore voluminous garbage printed to stdout.
* SQL output needs work.  Obvious things are directly updating db; figuring out what References are and creating them; using the proper Path format.. 

## sory state of affairs (Update April 13)

`hal2sg` only works on input of this form:
* **Mutation free:** ie as produced by Adam's pipeline
* **Root reference:** root genome used as initial side graph sequence.
* **Star Tree**
* **One sequence per genome:**

## next updates

* Fixes to output format?
* Clean code and tests
* Multiple sequences per genome (easy).  Already supported but not tested
* Cactus graphs.  Need to finish mutation support
* Arbitrary topology and references and genome subselection.  Much is done, just need to generalize path code, and handle paralogies in insertions, and of course tests.






