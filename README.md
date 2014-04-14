# fullsniplings

A package to infer full-sibling relationships from genetic data.
Geared toward datasets with many individuals and about 100 SNPs which are
treated as if they are unlinked. (While that is clearly false, it might be
reasonable).  The name derives from the fact that the idea for the approach
taken here came to me when I was working on sibling inference issues.  It
became clear that being able to compute the posterior predictive
distribution of genotypes conditional on observed sibling groups would be a
powerful way forward in the sibling inference problem. I will eventually generalize
that to more general pedigree inference, but for now just want to implement
the easy full-sibling case as a prototypical exercise.
