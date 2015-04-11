# hal2sg
Prototype code for converting HAL to Sequence/String/Side/? Graph 

## idea
Want a quick way to get Cactus genome alignments into the latest Human Variation graph structure.  This is mostly for internal testing purposes, as I understand. 

This is being written in a bigtime rush.   Functinality that is not working in this version. 

1. Paralogy in reference will not be written (not issue if reference is root)
2. Paralogy in inserted region (vs reference) not written (not issue if star tree with reference as root)
3. Point mutations (not issue if using Adam's hals as they have no subs)
4. Traceback ie paths for each genome.  These are computed just not written yet.
5. Join correction cases.  There are some situations where pairs of trivial or invalid joins need to be merged.  Shouldn't be hard to add.

