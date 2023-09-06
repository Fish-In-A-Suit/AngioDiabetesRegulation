import logging
import networkx as nx
from .GOTerm import GOTerm

logger = logging.getLogger(__name__)

class OboParser:
    def __init__(self, obo_filepath:str="src_data_files/go.obo"):
        """
        Parses the Gene Ontology OBO file.

        Params:
          - (str) obo_filepath: the filepath to the obo file
        """
        dag = nx.MultiDiGraph() # DiGraph cannot store multiple edges !!!

        def _reset_term_data():
            """
            Used during file parsing
            """
            return {
                'id': None, # obo: id
                'name': None, # obo: name
                'category': None, # obo: namespace
                'description': None, # obo: definition
                'parent_term_ids': [], # obo: is_a 
                'is_obsolete': False
            }
        
        # read all GO terms from the OBO file
        all_goterms = {} # mapping of all go ids to GOTerm objects
        with open(obo_filepath, 'r') as obo_file:
            term_data = {
                'id': None, # obo: id
                'name': None, # obo: name
                'category': None, # obo: namespace
                'description': None, # obo: def
                'parent_term_ids': [], # obo: is_a
                'is_obsolete': False
            }
            is_obsolete = False
            for line in obo_file: # also strips \n etc from lines
                line = line.strip()
                if line == "":
                    continue
                if "[Typedef]" in line: # typedefs are at the end of the file, no more go term data is expected
                    break
                if "[Term]" in line:
                    if term_data['id'] != None: # term_data != _reset_term_data check is used so that if all values in JSON are none, the goterm creation block isn't executed
                        # 'is_obsolete' is not present in all GO Terms. If it isn't present, set 'is_obsolete' to false
                        if 'is_obsolete' not in term_data:
                            term_data['is_obsolete'] = False

                        current_goterm = GOTerm(
                            id=term_data['id'],
                            name=term_data['name'],
                            category=term_data['category'],
                            description=term_data['description'],
                            parent_term_ids=term_data['parent_term_ids'],
                            is_obsolete=term_data['is_obsolete']
						)
                        all_goterms[current_goterm.id] = current_goterm
                    term_data = _reset_term_data() # reset term data for a new goterm
                else: # Term is not in line -> line is GO Term data -> process term data in this block
                    chunks = line.split(': ', 1) # split only first ": " element
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
                            is_obsolete = True if line_value == "true" else False
                            term_data['is_obsolete'] = is_obsolete
                            
        # all go terms from OBO are now constructed as GOTerm objects in all_goterms dictionary
        # create a Direcected Acyclic Graph from the created GO Terms
        for goid,goterm in all_goterms.items():
            assert isinstance(goterm, GOTerm)
            if dag.has_node(goterm.id) == False:
                dag.add_node(goterm.id) # add this goterm id as a node
            for parent_id in goterm.parent_term_ids:
                # nodes are automatically added if they are not yet in the graph when using add_edge
                dag.add_edge(all_goterms[parent_id].id, goterm.id) # add_edge(PARENT, CHILD) = "from parent to child"

        self.dag = dag
        self.all_goterms = all_goterms
        self.previously_computed_parents_cache = {} # cache dictionary between already computed goterms and their parents

    def get_parent_terms(self, term_id:str, return_as_class:bool = False, ordered:bool = True):
        """
        Gets all of GO Term parents of 'term_id'.

        Parameters:
          - (str) term_id: The GO Term whose parents you wish to obtain
          - (bool) return_as_class: If False, will return a list of string ids of parent GO Terms.
                                    If True, will return a list of GO Term parent classes.
          - (bool) ordered: If True, parents will be returned topologically (closest parents will be listed first in the returned list)
        
        Returns: A list of parent GO Terms (either ids or classes)
        """
        # attempt to cache old data
        if term_id in self.previously_computed_parents_cache:
            return self.previously_computed_parents_cache[term_id]
        
        parents = [] # WARNRING!! Using set() DESTROYS THE ORDER OF PARENTS !!!
        ancestors = nx.ancestors(self.dag, term_id)

        if ordered == True:
            # Calculate the distance from each ancestor to the given node
            distances = {ancestor: nx.shortest_path_length(self.dag, source=ancestor, target=term_id) for ancestor in ancestors}
            # Sort ancestors by distance in ascending order
            sorted_distances = dict(sorted(distances.items(), key=lambda item: item[1]))
            ancestors = sorted_distances.keys()

        for parent_id in ancestors:
            if return_as_class == True:
                parents.append(self.all_goterms[parent_id])
            else:
                parents.append(parent_id)
        
        # cache
        self.previously_computed_parents_cache[term_id] = parents
        
        return parents