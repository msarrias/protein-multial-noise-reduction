# Reducing noise in protein multialignments

In this research a noise-reduction method is implemented in protein multialigments
to evaluate its impact on phylogenic inference. The data used in
this project is a reduced data set from the original data used by the creators
of TrimAI.

To test the data the program fastprot was used to obtain the
distance matrices and fnj to infer a phylogenic tree for each alignment. The
Python package DendroPy was used to measure the symmetric distance between
the inferred tree (of the original and noise reduced alignment) and the
reference tree. In most of the cases the method
applied did not make any change between the noise reduced inferred tree
and the original inferred tree. This study focuses in determining whether the
noise reduction makes any effect in the symmetric distance calculation.
