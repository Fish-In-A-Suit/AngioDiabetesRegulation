# This test file demonstrates how to query GO Terms associated with a gene id.
# In the background, the following link is used: http://api.geneontology.org/api/bioentity/gene/{gene_id}/function
from goreverselookuplib.AnnotationProcessor import GOApi
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s [%(levelname)s] %(funcName)s: %(message)s')
# Create a file handler
file_handler = logging.FileHandler('./log_output/test_json_dump.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(logging.Formatter('%(asctime)s [%(filename)s:%(lineno)s] - [%(funcName)s] %(message)s'))
logger = logging.getLogger(__name__)

gene_id = "UniProtKB:P15692"
goapi = GOApi()
goterms_associated_with_geneid = goapi.get_goterms(gene_id, go_categories=['molecular_activity', 'biological_process', 'cellular_component'])
logger.info(f"len goterms = {len(goterms_associated_with_geneid)}")
logger.info(f"goterms = {goterms_associated_with_geneid}")