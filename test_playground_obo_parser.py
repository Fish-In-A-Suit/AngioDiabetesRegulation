# Displays the OboParser's GO Term parent query functionality

import logging
from goreverselookuplib .OboParser import OboParser

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s [%(levelname)s] %(funcName)s: %(message)s')
logger = logging.getLogger(__name__)

# test obo file parser
go_id = "GO:0060253" # negative regulation of glial cell proliferation
obo_parser = OboParser()
parent_terms = obo_parser.get_parent_terms(go_id)
logger.info(f"parent terms: {len(parent_terms)}")
for term in parent_terms:
    logger.info(f"  - {term}")
