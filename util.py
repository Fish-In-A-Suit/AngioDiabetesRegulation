# Utility functions file
import requests
import constants
import json

#x = requests.get("http://amigo.geneontology.org/amigo/search/ontology?q=angiogenesis")
#print(x.text)

def get_array_terms(array_name):
    """
    Returns the specified array terms, possible options:
      - 'ALL': all arrays combined
      - 'ANGIOGENESIS': all angiogenesis arrays combined
      - 'ANGIOGENESIS-POSITIVE'
      - 'ANGIOGENESIS-NEGATIVE'
      - 'ANGIOGENESIS-GENERAL'
      - 'DIABETES': all diabetes arrays combined
      - 'DIABETES-POSITIVE'
      - 'DIABETES-NEGATIVE'
      - 'DIABETES-GENERAL'
    """
    if array_name == 'ALL':
        return constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY + constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY
    elif array_name == 'ANGIOGENESIS':
        return constants.TERMS_ANGIOGENESIS_GENERAL + constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY + constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY
    elif array_name == 'ANGIOGENESIS-NEGATIVE':
        return constants.TERMS_ANGIOGENESIS_NEGATIVE_ARRAY
    elif array_name == 'ANGIOGENESIS-POSITIVE':
        return constants.TERMS_ANGIOGENESIS_POSITIVE_ARRAY
    elif array_name == 'ANGIOGENESIS-GENERAL':
        return constants.TERMS_ANGIOGENESIS_GENERAL
    elif array_name == 'DIABETES':
        return constants.TERMS_DIABETES_GENERAL + constants.TERMS_DIABETES_NEGATIVE_ARRAY + constants.TERMS_DIABETES_POSITIVE_ARRAY
    elif array_name == 'DIABETES-NEGATIVE':
        return constants.TERMS_DIABETES_NEGATIVE_ARRAY
    elif array_name == 'DIABETES-POSITIVE':
        return constants.TERMS_DIABETES_POSITIVE_ARRAY
    elif array_name == 'DIABETES-GENERAL':
        return constants.TERMS_DIABETES_GENERAL
    else:
        print(array_name + " could not be found! Returning empty array.")
        empty = []
        return empty

def read_file_as_json(filepath):
    """
    Reads the file into json
    """
    with open(filepath, "r") as read_content:
        return json.load(read_content)


# print(read_file_as_json("demofile2.json")) #this now works but breaks the console xd