from setuptools import setup, find_packages

setup(
    name="goreverselookup", #must match the root foldername
    version="0.0.1",
    author="Vladimir Smrkolj, Aljoša Škorjanc",
    author_email="your.email@example.com",
    description="A short description of your library",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'requests',
        'tqdm',
        'tabulate',
    ],
)