import logging
import networkx as nx
from .GOTerm import GOTerm

logger = logging.getLogger(__name__)

class OboParser:
    def __init__(obo_filepath:str):
        """
        Parses the Gene Ontology OBO file.

        Params:
          - (str) obo_filepath: the filepath to the obo file
        """
        dag = nx.DiGraph()

        def _reset_term_data():
            """
            Used during file parsing
            """
            return {
                'id': None, # obo: id
                'name': None, # obo: name
                'category': None, # obo: namespace
                'description': None, # obo: definition
                'parent_term_ids': [] # obo: is_a 
            }
        
        all_goterms = {} # mapping of all go ids to GOTerm objects
        with open(obo_filepath, 'r') as obo_file:
            term_data = {
                'id': None, # obo: id
                'name': None, # obo: name
                'category': None, # obo: namespace
                'description': None, # obo: def
                'parent_term_ids': [] # obo: is_a 
            }
            is_obsolete = False
            for line in obo_file:
                if line == "":
                    continue
                if "[Term]" in line:
                    term_data = _reset_term_data()
                    if is_obsolete == False:
                        current_goterm = GOTerm(
                            id=term_data['id'],
                            name=term_data['name'],
                            category=term_data['category'],
                            description=term_data['description'],
                            parent_term_ids=term_data['parent_term_ids']
						)
                        all_goterms[current_goterm.id] = current_goterm
                    else:
                        logger.debug(f"Term {term_data['id']} is obsolete! It wasn't added to DAG.")
                        is_obsolete = False # reset so other terms aren't faultily missed
                    
                line = line.strip()
                chunks = line.split(' ', 1) # split only first spacebar
                line_identifier = chunks[0]
                line_value = chunks[1]
                match line_identifier:
                    case 'id':
                        term_data['id'] = line_value
                    case 'name':
                        term_data['name'] = line_value
                    case 'def':
                        line_value = line_value.strip("\"") # definition line value contains double quotes in obo, strip them
                        term_data['description'] = line_value
                    case 'namespace':
                        term_data['category'] = line_value
                    case 'is_a':
                        line_value = line_value.split(' ')[0] # GO:0000090 ! mitotic anaphase -> split into GO:0000090
                        term_data['parent_term_ids'].append(line_value)
                    case 'is_obsolete':
                        if line_value == "true":
                            is_obsolete = True

        # all go terms from OBO are now constructed as GOTerm objects in all_goterms dictionary
        for goid,goterm in all_goterms.items():
            assert isinstance(goterm, GOTerm)
            dag.add_node(goterm) # add this goterm as a node
            for parent_id in goterm.parent_term_ids:
                dag.add_edge(goterm, all_goterms[parent_id]) # add all goterm parents to this goterm node

        # TODO: based on this, create the function to get all parents
        #def get_parent_terms(dag, term_id):
        #    parents = set()
        #    for parent in nx.ancestors(dag, term_id):
        #        parents.add(parent)
        #    return parents

        # Find all parent terms for GO:0001938 (positive regulation of endothelial cell proliferation)
        #term_id_to_find = 'GO:0001938'
        #parent_terms = get_parent_terms(dag, term_id_to_find)