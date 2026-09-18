from setuptools import find_packages, setup


setup(
    name="IGenotyper",
    version="1.2.0",
    packages=find_packages(),
    python_requires=">=3.9",
    url="https://github.com/oscarlr/IGenotyper",
    author="Oscar Rodriguez",
    author_email="oscar.rodriguez@icahn.mssm.edu",
    description="Long-read genotyping of immunoglobulin loci",
    package_data={
        "IGenotyper": [
            "data/*.bed",
            "data/alleles.fasta",
            "data/alleles.fasta.fai",
            "data/assembly.sh",
            "data/reference.fasta.fai",
            "data/rhesus/*",
            "data/immune_receptor_genomics/*",
            "data/immune_receptor_genomics/*/*",
            "templates/*",
            "scripts/*",
            "ancestry/*",
        ]
    },
    entry_points={"console_scripts": ["IG=IGenotyper.main:main"]},
)
