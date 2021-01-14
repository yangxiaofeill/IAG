# IAG for three Papaver species
To figure out what genome structural changes occured to shape the present-day karyotypes of three Papaver species. We developed a novel bottom-up framework IAG (inferring ancestor genome) to reconstruct pre- and post-WGD-1 ancestors of Papaver genomes. This is a pipline only for Papaver species. And a more user-friendly version and article are in developing.

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

### inferringAncestorGenomeStructure
Containing the source code of three models for inferring ancestor genome structures which are focusing on the block connection, including:

* GMP(genome median problem): Given some ordinary genome(with no WGD, one copy), finding the ancestor genome by minimizing the genomic distance.

* GGAP(guided genome aliquoting problem): Given a duplicated genome(with a WGD or WGT, two or three copies) and a ordinary genome, finding pre-duplicated ancestor genome by minimizing the genomic distance.

* MO(matcing optimaztion): Given two multi-duplicated genome, finding the block matcing relation by minimaizing the genomic distance.

We built three integer programming model to solve above problems. All above optimization instances were solved with the GUROBI solver.

### inferringGeneContent
Containing the source code for inferring ancestor genome gene order. Each block contains multi-copies in sub current genomes and with some different gene indels and micro rearrangements. This step is intend to infer the possible ancestral gene order in ancestors. We used topological sorting with a greedy strategy. When the topological sorting removed nodes, we recorded their out-degree edges weight and used this weight as decision condition. When there is not node without in-degree edge (the number of in-degree is 0), we removed the node with maximal decision weight to break cycle and continue sorting to get the final order.

### inputFiles

* block files which is the block order and direction of three species. P. rhoeas(PR), P. somniferum(PS) and P. setigerum(PT).

* multi-synteny files which is the block gene content in multi-copies with the homo gene ID and the gene name.





