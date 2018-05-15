from setuptools import setup

setup(
    name = "igseqanalysis",
    packages = ['igseqanalysis'],
    install_requires=[
        'pandas',
        'biopython'
    ],
    entry_points = {
        "console_scripts": ['igblast_IG =  igseqanalysis.igblast_IG:main',
                            'igblast_TCR =  igseqanalysis.igblast_TCR:main',
                            'parse_igblast = igseqanalysis.parse_igblast:main',
                            'pair_by_id = igseqanalysis.pair_by_id:main',
                            'csv2fasta = igseqanalysis.csv2fasta:main',
                            'format_cluster = igseqanalysis.format_cluster:main',
                            'translate_table = igseqanalysis.translate_table:main',
                            'count_unique =  igseqanalysis.count_unique:main'
                            ]
        },
    package_data = {'igseqanalysis': ['imgt/optional_file/*.aux', 'imgt/ig/*', 'imgt/tcr/*']},
    include_package_data=True,
    version = "0.2.1",
    description = "Python command line application to process IgBlast alignments",
    author = "Jincheng Wu",
    author_email = "wuji@medimmune.com",
    )
