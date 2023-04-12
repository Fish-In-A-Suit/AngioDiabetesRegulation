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
pip install dist/goreverselookup-version.tar.gz
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
