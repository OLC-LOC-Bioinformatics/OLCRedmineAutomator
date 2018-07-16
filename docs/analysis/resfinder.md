# ResFinder

### What does it do?

ResFinder is a program developed by the Danish Center for Genomic Epidemiology for detection of antibiotic resistance
in draft genome assemblies. It is very important to note that the Redmine version will only look for acquired antibiotic
resistance genes (generally plasmid-borne) and not chromosomally encoded AMR genes that are caused by point mutations.

If you're interested in chromosomally encoded AMR genes, you can [External Retrieve](data/external_retrieve.md) your
assemblies of interest and submit them to an alternate AMR predictor, such as McMaster's [CARD](https://card.mcmaster.ca/analyze/rgi).

### How do I use it?

#### Subject

In the `Subject` field, put `ResFinder`. Spelling counts, but case sensitivity doesn't.

#### Description

All you need to put in the description is a list of SEQIDs you want to detect AMR in, one per line.

#### Example

For an example ResFinder, see [issue 12854](https://redmine.biodiversity.agr.gc.ca/issues/12854).

### How long does it take?

ResFinder is very fast - it should only take a few seconds to analyze each SEQID requested.

### What can go wrong?

1) Requested SEQIDs are not available. If we can't find some of the SEQIDs that you request, you will get a warning
message informing you of it.

