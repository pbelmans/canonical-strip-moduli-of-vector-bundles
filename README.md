# Counterexamples to the canonical strip hypothesis

This code checks the canonical strip hypothesis (and its variations) for the moduli space of rank 2 bundles with fixed determinant on a smooth projective curve of genus 2. It follows that this hypothesis fails for genus at least 10.

# Introduction

The Hilbert polynomial of the homogeneous coordinate ring associated to the anticanonical bundle is subject to [Golyshev's hypothesis](https://mathscinet.ams.org/mathscinet-getitem?mr=2503098) which states for any Fano variety and anticanonical polarisation that

> The real parts of roots of the Hilbert polynomial are all negative.

Let C be a smooth projective curve of genus g at least 2. Then the moduli space of rank 2 bundles with fixed determinant of odd degree is a smooth projective Fano variety of dimension 3g-3 and Picard rank 1. The code in this repository accompanies a paper that shows that this hypothesis is *em* satisfied for these Fano varieties, provided that the genus is at least 10.

# How to use this code

There are two files, which implement the same computation:

* `canonical-strip.sage` is written in [Sage](https://sagemath.org)
* `canonical-strip.gp` is written in [Pari/GP](https://pari.math.u-bordeaux.fr)

The code checks the canonical strip hypothesis for the moduli space of vector bundles by computing sufficiently many terms in the Hilbert series using the Verlinde formula, computing the Hilbert polynomial from this, and then checking the location of its roots. The paper also discusses variations on the canonical strip hypothesis for varieties which can be constructed from the moduli space of vector bundles, and in all these cases the corresponding hypothesis is violated, starting around genus 10.
