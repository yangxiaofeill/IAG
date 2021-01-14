# IAG for three Papaver species
To figure out what genome structural changes occured to shape the present-day karyotypes of three Papaver species. We developed a novel bottom-up framework IAG (inferring ancestor genome) to reconstruct pre- and post-WGD-1 ancestors of Papaver genomes. This is a pipline only for Papaver species. And a more user-friendly version is being developed.

## Dependencies
python 3.6.4
```Bash
conda create -n environment python=3.6.4
```
```Bash
conda install packge=version
```
packges  | version|
--------- | --------|
numpy  | 1.16.4 |
pandas  | 0.20.3 |
networkx | 2.1 |

[Gurobi 9.0.2](https://www.gurobi.com)

## Introductions

inferringAncestorGenomeStructure
Containing the source code of three inferring models, including:

* GMP(genome median problem): Given some ordinary genome(with no WGD, one copy), finding the ancestor genome by minimizing the genomic distance.

* GGAP(guided genome aliquoting problem): Given a duplicated genome(with a WGD or WGT, two or three copies) and a ordinary genome, finding pre-duplicated ancestor genome by minimizing the genomic distance.

* MO(matcing optimaztion): Given two multi-duplicated genome, finding the block matcing relation by minimaizing the genomic distance.

We built three integer programming model to solve above problems. All above optimization instances were solved with the GUROBI solver.

inferringGeneContent

