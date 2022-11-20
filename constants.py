# MISC:
#   File tree in term_genes:
#     - homosapiens_only = use Ladi's new get_go_genes_API that only returns homosapiens taxon without any test-animal orthologs
#
#
#
#

# Positive regulatory terms for search "angiogenesis"
TERMS_ANGIOGENESIS_POSITIVE_ARRAY = [
    "positive regulation of blood vessel endothelial cell proliferation involved in sprouting angiogenesis", "GO:1903589",
    "positive regulation of sprouting angiogenesis", "GO:1903672",
    "positive regulation of angiogenesis", "GO:0045766",
    "positive regulation of cell adhesion involved in sprouting angiogenesis", "GO:0106090",
    "positive regulation of cell migration involved in sprouting angiogenesis", "GO:0090050",
    "endothelial cell proliferation","GO:0001935"
    ]

# Negative regulatory terms for search "angiogenesis"
TERMS_ANGIOGENESIS_NEGATIVE_ARRAY = [
    "negative regulation of blood vessel endothelial cell proliferation involved in sprouting angiogenesis", "GO:1903588",
    "negative regulation of sprouting angiogenesis", "GO:1903671",
    "negative regulation of angiogenesis", "GO:0016525",
    "negative regulation of cell adhesion involved in sprouting angiogenesis", "GO:0106089",
    "negative regulation of cell migration involved in sprouting angiogenesis", "GO:0090051",
    "angiostatin binding", "GO:0043532",
    "angiogenin-PRI complex", "GO:0032311",
    "negative regulation of endothelial cell proliferation","GO:0001937",
    "negative regulation of vascular endothelial cell proliferation","GO:1905563",
    ]

# Positive regulatory terms for search "diabetes"
TERMS_DIABETES_POSITIVE_ARRAY = [
    "negative regulation of cellular response to insulin stimulus","GO:1900077",
    "negative regulation of insulin receptor signaling pathway","GO:0046627",
    ]

# Negative regulatory terms for search "diabetes"
TERMS_DIABETES_NEGATIVE_ARRAY = [
    "positive regulation of cellular response to insulin stimulus","GO:1900078",
    "positive regulation of insulin receptor signaling pathway","GO:0046628",
    ]

# Terms that are neither positive or negative
TERMS_ANGIOGENESIS_GENERAL = [
    "angiogenesis", "GO:0001525",
    "regulation of blood vessel endothelial cell proliferation involved in sprouting angiogenesis", "GO:1903587",
    "regulation of sprouting angiogenesis", "GO:1903670",
    "regulation of angiogenesis", "GO:0045765",
    "angiogenesis involved in wound healing", "GO:0060055",
    "sprouting angiogenesis", "GO:0002040",
    "intussusceptive angiogenesis", "GO:0002041",
    "cell migration involved in sprouting angiogenesis", "GO:0002042",
    "blood vessel endothelial cell proliferation involved in sprouting angiogenesis", "GO:0002043",
    "blood vessel endothelial cell migration involved in intussusceptive angiogenesis", "GO:0002044",
    "regulation of cell adhesion involved in sprouting angiogenesis", "GO:0106088",
    "regulation of cell adhesion involved in intussusceptive angiogenesis", "GO:0002045",
    "cell adhesion involved in sprouting angiogenesis", "GO:0120078",
    "regulation of cell migration involved in sprouting angiogenesis", "GO:0090049",
    "alpha6-beta1 integrin-CYR61 complex", "GO:0071116",
    "angiogenic sprout fusion", "GO:0120077",
    "blood vessel endothelial cell migration", "GO:0043534",
    "endothelium development","GO:0003158",
    "positive regulation of endothelial cell proliferation","GO:0001938",
    "positive regulation of vascular endothelial cell proliferation","GO:1905564",
    "regulation of endothelial cell proliferation","GO:0001936",
    "regulation of vascular endothelial cell proliferation","GO:1905562",
    "vascular endothelial cell proliferation","GO:0101023",
    ]

TERMS_DIABETES_GENERAL = [
    "cellular ketone body metabolic process", "GO:0046950",
    "cellular response to insulin stimulus","GO:0032869",
    "glucose import in response to insulin stimulus","GO:0044381",
    "glucose metabolic process","GO:0006006",
    "insulin binding","GO:0043559",
    "insulin catabolic process","GO:1901143",
    "insulin metabolic process","GO:1901142",
    "insulin receptor complex","GO:0005899",
    "insulin receptor internalization","GO:0038016",
    "insulin receptor recycling","GO:0038020",
    "insulin receptor signaling pathway via phosphatidylinositol 3-kinase","GO:0038028",
    "insulin receptor signaling pathway","GO:0008286",
    "insulin receptor substrate binding","GO:0043560",
    "insulin-activated receptor activity","GO:0005009",
    "insulin-like growth factor binding","GO:0005520",
    "insulin-responsive compartment","GO:0032593",
    "regulation of insulin receptor signaling pathway","GO:0046626",
    "negative regulation of skeletal muscle cell proliferation","GO:0014859",
    "negative regulation of skeletal muscle hypertrophy","GO:1904205",
    "negative regulation of skeletal muscle tissue growth","GO:0048632",
    "positive regulation of skeletal muscle cell proliferation","GO:0014858",
    "positive regulation of skeletal muscle hypertrophy","GO:1904206",
    "positive regulation of skeletal muscle tissue growth","GO:0048633",
    "regulation of skeletal muscle adaptation","GO:0014733",
    "regulation of skeletal muscle cell proliferation","GO:0014857",
    "regulation of skeletal muscle hypertrophy","GO:1904204",
    "regulation of skeletal muscle tissue growth","GO:0048631",
    "skeletal muscle adaptation","GO:0043501",
    "skeletal muscle atrophy","GO:0014732",
    "skeletal muscle cell proliferation","GO:0014856",
    "skeletal muscle fiber differentiation","GO:0098528",
    "skeletal muscle hypertrophy","GO:0014734",
    "skeletal muscle satellite cell proliferation","GO:0014841",
    "skeletal muscle tissue development","GO:0007519",
    "skeletal muscle tissue growth","GO:0048630",
    "skeletal muscle tissue regeneration","GO:0043403",
    ]

# Terms that came up when searching "angiogenesis", but are irrelevant for us
TERMS_EXCLUDED = [
    "cell migration involved in coronary angiogenesis", "GO:0060981",
    "angiogenesis involved in coronary vascular morphogenesis", "GO:0060978",
    ]

# array of genes (like the zebrafish), that have been confirmed to be orthologs by the user
TRUSTED_GENES = []

# array of empty terms (get_GO_genes_API returned empty list), loaded at startup from appropriate terms_empty.txt file
TERMS_EMPTY = []

GENE_DESCRIPTION_ANALYSIS_KEYWORDS = [
    "vessel", "angiogenesis", "capillary", "artery", "arteriole"
]

# single array enums in use to make the cross section between angiogenesis and diabetes terms
TERM_ARRAYS_ENUMS_SINGLE = {
    "terms_angiogenesis_positive": TERMS_ANGIOGENESIS_POSITIVE_ARRAY,
    "terms_angiogenesis_negative": TERMS_ANGIOGENESIS_NEGATIVE_ARRAY,
    "terms_angiogenesis_general": TERMS_ANGIOGENESIS_GENERAL,
    "terms_diabetes_positive": TERMS_DIABETES_POSITIVE_ARRAY,
    "terms_diabetes_negative": TERMS_DIABETES_NEGATIVE_ARRAY,
    "terms_diabetes_general": TERMS_DIABETES_GENERAL,
}

# multiple array enums in use to make the cross section between angiogenesis and diabetes terms
TERM_ARRAY_ENUMS_MULTIPLE = {
    "terms_diabetes_all": TERMS_DIABETES_GENERAL + TERMS_DIABETES_NEGATIVE_ARRAY + TERMS_DIABETES_POSITIVE_ARRAY,
    "terms_angiogenesis_all": TERMS_ANGIOGENESIS_GENERAL + TERMS_ANGIOGENESIS_POSITIVE_ARRAY +  TERMS_ANGIOGENESIS_NEGATIVE_ARRAY
}