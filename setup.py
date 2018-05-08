from setuptools import setup

setup(
    name = "igseqanalysis",
    packages = ['igseqanalysis'],
    entry_points = {
        "console_scripts": ['igblastIG =  igseqanalysis.igblastIG:main',
                            'parseIgBlast = igseqanalysis.parseIgBlast:main',
                            'pairById = igseqanalysis.pairByID:main',
                            'csv2fasta = igseqanalysis.csv2fasta:main',
                            'formatCluster = igseqanalysis.formatCluster:main',
                            'translateTable = igseqanalysis.translateTable:main',
                            'countUnique =  igseqanalysis.countUnique:main'
                            ]
        },
    package_data = {'igseqanalysis': ['imgt/optional_file/*.aux', 'imgt/ig/*', 'imgt/tcr/*']},
    include_package_data=True,
    version = "0.2.0",
    description = "Python command line application to process IgBlast alignments",
    author = "Jincheng Wu",
    author_email = "wuji@medimmune.com",
    )
