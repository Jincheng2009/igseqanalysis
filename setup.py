import re
from setuptools import setup, find_packages

setup(
    name = "igseqanalysis",
    packages = ['igseqanalysis'],
    entry_points = {
        "console_scripts": ['parseigblast = igseqanalysis.parseIgBlast:main',
                            'pairbyid = igseqanalysis.pairByID:main',
                            'csv2fasta = igseqanalysis.csv2fasta:main',
                            'formatcluster = igseqanalysis.formatCluster:main',
                            'translatetable = igseqanalysis.translateTable:main',
                            'countunique = igseqanalysis.countUnique:main'
                            ]
        },
    version = "0.0.1",
    description = "Python command line application to process IgBlast alignments",
    author = "Jincheng Wu",
    author_email = "wuji@medimmune.com",
    )