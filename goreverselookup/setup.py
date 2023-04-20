from setuptools import setup, find_packages

setup(
    name="goreverselookup",
    version="{{VERSION_PLACEHOLDER}}",
    author="Vladimir Smrkolj, Aljoša Škorjanc",
    author_email="vladimir.smrkolj@gmail.com",
    description="Library for Gene Ontology Reverse Lookup",
    url="https://github.com/SurgeonsStrikeBack/GeneOntologyAnalysis",
    packages=['goreverselookuplib'],
    include_package_data=False,
    install_requires=[
        'requests',
        'tqdm',
        'tabulate',
    ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows"
    ]
)