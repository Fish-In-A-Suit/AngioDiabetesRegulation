# Gene Ontology Reverse Lookup Tool

[![DOI](/doc/images/DOI.svg)](not yet)
[![Latest PyPI version](https://img.shields.io/pypi/v/goatools.svg)](not yet)
[![bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](not yet)
[![Github Actions](https://github.com/tanghaibao/goatools/workflows/build/badge.svg)](not yet)
[![Downloads](https://pepy.tech/badge/goatools)](not yet)

|         |                                                                       |
| ------- | --------------------------------------------------------------------- |
| Authors | Vladimir Smrkolj ([ladismrkolj](http://github.com/ladismrkolj))       |
|         | Aljoša Škorjanc ([Fish-In-A-Suit](https://github.com/Fish-In-A-Suit)) |
| Email   | <??@gmail.com>                                                        |
| License | cc-by-nc                                                              |

## Instalation TODO
```
cd goreverselookup
python setup.py sdist
pip install dist/goreverselookup-{VERSION}.tar.gz
```
For developers:
```
cd goreverselookup
pip install -e .
```
TODO: how to publish

## Usage
see: TestReverseLookupLibrary.py
Currently I have to sort out the ugly imports!

```python
import ReverseLookup

# Define model from input file
model = ReverseLookup.ReverseLookup.from_input_file("diabetes_angio_1/input.txt")

# Fetch all GO term names and descriptions
model.fetch_all_go_term_names_descriptions()

# Fetch all GO term products
model.fetch_all_go_term_products()

# Create products from GO terms
model.create_products_from_goterms()

# Fetch UniProt ID products
model.fetch_UniprotID_products()

# Prune products
model.prune_products()

# Fetch UniProt information
model.fetch_Uniprot_infos()

# Score products
model.score_products()

# Optional: Fetch mRNA sequences
model.fetch_mRNA_sequences()

# Predict miRNAs
model.predict_miRNAs()

# Score miRNAs
model.score_miRNAs()

# Generate report
report = ReverseLookup.ReportGenerator(model, verbosity=3)
report.general_report("diabetes_angio_1/general.txt")

# Save model
model.save_model("diabetes_angio_1/data.json") # save to a specific filepath
model.save_model() # save to program_data_files/models/self.model_name/data.json

# Alternatively, you can use model.compute_all(str:report_filepath) to achieve all of the above steps:
model.compute_all("diabetes_angio1/report.txt")
```

In future runs, you can import the model with:

```python
model = ReverseLookup.ReverseLookup.load_model("diabetes_angio_1/data.json")
# or
model = ReverseLookup.ReverseLookup.load_model(model_name="V1") # if a model with the same model_name has been saved before
```

### TODO: Future additions:
  - [ ] GAF (offline master file) parsing for uniprotkb orthologs
mogoče bi se dal sparsat iz GAF (Gene Association File-a) podobno kar parsava na internetu. Ker trenutno se mi zdi, da je parsing za GAF napisan tko, da samo uposteva prvi in drugi column (kar predstavlja uniprotkb kodo) in pol isce katermu GOTermu to pripada. Sm pa lihkar opazu, da je v GAF tut dodatn column, kjer zgleda so zapisani tudi vsi analogi, kot npr. zfin in drugi: 
![todo1 image](https://i.ibb.co/M8ks35R/todo1.png)

  - [ ] input file parse for homo sapiens only has no logic
Classmethod _from_input_file_ in ReverseLookup class has no logic specified for parsing HomoSapiens_only = True from input.txt. Fix this!


### Legacy Usage
This section explains how to use the legacy, non-class workflow.

1. Create an empty folder and set constants.py -> TARGET_FOLDER to point to it's filepath. This is the folder where the results will be saved.
2. Populate lists in constants.py with the terms you wish to analyse
3. Run dev_go_download.py: it downloads a list of products associated with each GO term. Resulting folder is saved as terms_direct_products.json into TARGET_FOLDER. By default, it processes only homosapiens products and excludes non-homosapiens related products.
4. Run dev_scoring_products_basic.py using useLegacy=True -> result = product_scores.json saved to TARGET_FOLDER
5. Run dev_mrna_download.py. Results are product_mRNA.json and product_mRNA_refseq.json saved to TARGET_FOLDER
6. Run dev_predict_miRNA_v2.py. It expects product_mRNA_refseq.json to exist in TARGET FOLDER, result is product_mRNA_miRDB-predicted-miRNA-matching.json
7. Run dev_scoring_miRNA_basic.py; result is mirna_scores.json.
8. Run dev_generate_report.py; result is report.txt
