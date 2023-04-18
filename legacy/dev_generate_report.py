import util as util
import constants as constants
import os
from tabulate import tabulate, SEPARATING_LINE
import requests
import re

import logging
logger = logging.getLogger(__name__)

def generate_report(top_n):

    string = "Term Enrichment System, Version 1\n"
    string += "Authors: Vladimir Smrkolj (SI), Aljosa Skorjanc (SI)\n"
    string += "March 2023\n"
    string += "-"*50+"\n"
    string += f"NOTE: Where applicable, displaying top {top_n} relevant results!\n"
    #string += "INPUT:"
    string += "-"*50+"\n"
    string += "PRODUCT SCORES".center(50)+"\n"
    string += "-"*50+"\n"

    scores = util.read_file_as_json(os.path.join(constants.TARGET_FOLDER, "product_scores.json"))
    top_products=[]
    for i in range(top_n):
        gene = scores[i]["product"]

        gene_id = gene.split(":")[1]
        query = util._uniprot_query_API(gene_id, type="prot") #Get protein name from UniProtAPI
        index = next((index for (index, d) in enumerate(query["results"]) if d["primaryAccession"] == gene_id), None)
        protein_name = query["results"][index]["proteinDescription"]["recommendedName"]["fullName"]["value"]
        scores[i]["protein_name"] = protein_name

        top_products.append([scores[i]["product"], f"{scores[i]['al_corr_score']:.2f}", protein_name])
    bottom_products=[]
    for i in range(len(scores)-top_n-1, len(scores)):
        gene = scores[i]["product"]

        gene_id = gene.split(":")[1]
        query = util._uniprot_query_API(gene_id, type="prot")
        index = next((index for (index, d) in enumerate(query["results"]) if d["primaryAccession"] == gene_id), None)
        protein_name = query["results"][index]["proteinDescription"]["recommendedName"]["fullName"]["value"]
        scores[i]["protein_name"] = protein_name

        bottom_products.append([scores[i]["product"], f"{scores[i]['al_corr_score']:.2f}", protein_name])
    
    table1 = [["UniProtID", "Score", "Protein name"]] + top_products + [["----", "----", "----"]] + bottom_products
    string += tabulate(table1, headers="firstrow", tablefmt="pretty") + "\n\n"

    string += "DETAILED OUTPUT".center(50)+"\n\n"

    #now lets describe each entry in detail
    for i in range(top_n):
        string += " "*10+f"{scores[i]['product']} - {scores[i]['protein_name']} - {scores[i]['al_corr_score']:.2f}".center(100)+"\n"
        tables_go = []
        tables_go.append(["GO term", "GO label", "GO description"])
        for term in scores[i]["terms:"]:
            response = requests.get(f"http://api.geneontology.org/api/ontology/term/{term}", params={}).json() #"goid", "label", "definition"
            #temp_text = f"{term} - {response['label']} - {response['definition']}"
            #string += re.sub("(.{64})", "\\1\n", temp_text, 0, re.DOTALL)+"\n"
            tables_go.append([term, response["label"], response["definition"]])
        string += tabulate(tables_go, headers="firstrow", tablefmt="grid", maxcolwidths = 100) + "\n\n"

    string += ("-"*30)+"\n\n"

    for i in range(len(scores)-top_n-1, len(scores)):
        string += " "*10+f"{scores[i]['product']} - {scores[i]['protein_name']} - {scores[i]['al_corr_score']:.2f}"+"\n"
        tables_go = []
        tables_go.append(["GO term", "GO label", "GO description"])
        for term in scores[i]["terms:"]:
            response = requests.get(f"http://api.geneontology.org/api/ontology/term/{term}", params={}).json() #"goid", "label", "definition"
            #temp_text = f"{term} - {response['label']} - {response['definition']}"
            #string += re.sub("(.{64})", "\\1\n", temp_text, 0, re.DOTALL)+"\n"
            tables_go.append([term, response["label"], response["definition"]])
        string += tabulate(tables_go, headers="firstrow", tablefmt="grid", maxcolwidths = 100) + "\n\n"
    

    string += "-"*50+"\n"
    string += "miRNA SCORES".center(50)+"\n"
    string += "-"*50+"\n"

    mirna_scores = util.read_file_as_json(os.path.join(constants.TARGET_FOLDER, "mirna_scores.json"))
    table_mirna = []
    table_mirna.append(["miRNA", "score", "suppressed products"])
    for i in range(top_n):
        processed_overlaps = ""
        temp_list=[]
        for j in range(len(mirna_scores[i]["products_inhibited"])):
            product = mirna_scores[i]["products_inhibited"][j]
            if any(product in sub for sub in top_products) or any(product in sub for sub in bottom_products):
                temp_list.insert(0,f"*{product}")
            else:
                temp_list.append(f"{product}")
        for j in range(len(temp_list)): 
            processed_overlaps += temp_list[j]
            if j != len(mirna_scores[i]["products_inhibited"])-1:
                processed_overlaps += "\n"
        table_mirna.append([mirna_scores[i]["mirna"], mirna_scores[i]["basic_score"], processed_overlaps])
    string += tabulate(table_mirna, headers="firstrow", tablefmt="grid") + "\n\n"
    


    logger.debug(string)
    return string

def main():
    out = generate_report(top_n=5)

    util.save_txt(out, os.path.join(constants.TARGET_FOLDER, "report.txt"))    

if __name__ == '__main__':
    import logging.config
    import logging_config as cfg
    logging.config.dictConfig(cfg.log_dict)
    main()