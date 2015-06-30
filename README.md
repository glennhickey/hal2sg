# hal2sg
Prototype code for converting [HAL](https://github.com/glennhickey/hal) to [Side Graph SQL](https://github.com/ga4gh/schemas/wiki/Human-Genome-Variation-Reference-(HGVR)-Pilot-Project#graph-format)

(c) 2015 Glenn Hickey. See [LICENSE](https://github.com/glennhickey/hal2sg/blob/development/LICENSE) for details.

## Algorithm

Iteratatively add genomes to side graph.   The HAL alignment is used to thread each successive genome onto the nearest genome already in the graph.   Genomes are added in breadth-first search order starting from the reference (specifiable by option, root by default).

### CAMEL Input

The following special logic is used to support HAL files generated by [CAMEL](https://github.com/adamnovak/sequence-graphs):  If the root DNA sequence is all N's, then it will be automatically inferred from the sequences of its children. 

## Instructions

**Dependencies:**   [HAL](https://github.com/glennhickey/hal)  and [SonLib](https://github.com/benedictpaten/sonLib).  Expected to be in sister directories to hal2sg but can be changed in include.mk.  Tested only with latest version from development branch for both. 

**Note** Useful (but very slow) to compile in debug mode (cppflags_dbg in sonLib/include.mk), as asserts are used to check if the output is honest.  At the very least run the unit tests:

	  make test

To run the converter:

	  hal2sg input.hal output.fa output.sql

`input.hal` Input alignment to convert

`output.fa` Output fasta file of all Side Graph sequences

`output.sql` Output text file listing INSERT commands for Sequences, Joins and Paths (for each input sequence) in the graph.

To see all the options, run with no args or use `--help`.




