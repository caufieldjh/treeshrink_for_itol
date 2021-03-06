# treeshrink_for_itol.py
**A small utility for summarizing protein presence/absence in summary phylogenetic trees.**

Aids in creating annotated trees for the Interactive Tree of Life
(iTOL, http://itol.embl.de/). 

Uses ETE (http://etetoolkit.org/) for tree analysis and NCBI Taxonomy handling.

ETE will already do many of the things iTOL does...
but iTOL is friendlier to programming laymen.

This utility also relies on PhyloT (http://phylot.biobyte.de/) to ensure
that trees match NCBI taxonomy and to convert trees to iTOL output.

**INPUT**: 
The name of an eggNOG v4 NOG/COG or a Newick tree as output from and eggNOG NOG. 

The script assumes that each leaf node in the tree indicates a single genome containing a member of an orthologous group.

**OUTPUT**:
Two files: a list of taxids for assembly into a new tree and an iTOL annotation file.
