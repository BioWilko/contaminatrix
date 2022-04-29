import os
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), "contaminatrix", "version.py")
version = open(version_py).read().strip().split("=")[-1].replace('"', "").strip()
# long_description = """
# ``align_trim`` is a standalone version of the align_trim step of the artic fieldbioinformatics pipeline for the purposes of primer-trimming and normalising amplicon-based viral sequencing data
# """

HERE = os.path.dirname(__file__)

with open(os.path.join(HERE, "requirements.txt"), "r") as f:
    install_requires = [x.strip() for x in f.readlines()]

setup(
    name="contaminatrix",
    version=version,
    install_requires=install_requires,
    requires=["python (>=3.5)"],
    packages=["contaminatrix"],
    author="Sam Wilkinson",
    description="A tool for spotting contamination in amplicon-based sequencing data by exploting read-group overlaps",
    # long_description=long_description,
    url="",  # Add later
    package_dir={"contaminatrix": "contaminatrix"},
    package_data={"contaminatrix": []},
    zip_safe=False,
    include_package_data=True,
    entry_points={"console_scripts": ["contaminatrix=contaminatrix.cli:main"],},
    author_email="s.a.j.wilkinson@bham.ac.uk",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
