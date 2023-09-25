# ML-FWRF (Multi-label feature weighted random forest)
Code for correcting for lineage-dependency in random forest prediction of drug-resistant phenotypes from *Mycobacterium tuberculosis* genome sequences using multi-label feaature weighted ranom forest.
Code is supplied in pipelines which can be run using the following examples. 
There are a few requirements:
* A csv file containing the train dataset. Samples are named 'id'. Phenotypes are name of outcome to be used in the model.
* A csv file containing the test dataset. Samples are named 'id'. Phenotype columns are used to assess performance of the model.
* Phenotypes (character vector for phenotypes)
* Phenotype oder (character vector describing the order of the chain)
* File containing feature weights
* Output id (tag to add to results file)
* Drug (for interaction analysis only)

## Examples of how to run code on command line
### MLFWRF
```
Rscript MLFWRF.R --train "train_dataset.csv" --test "train_dataset.csv"  --phenotypes c("rifampicin", "isoniazid","ethambutol", "pyrazinamide", "streptomycin", "kanamycin","capreomycin", "amikacin", "ofloxacin", "moxifloxacin") --pheno_order c("rifampicin", "isoniazid","ethambutol", "pyrazinamide", "streptomycin", "ofloxacin", "moxifloxacin","amikacin","capreomycin", "kanamycin") --fw_file "feature_weights.txt" --output_id "order1"
```
