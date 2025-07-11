from setuptools import setup

setup(
    name="evodock",
    version="1.0.0",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    description="EvoDOCK: protein-protein docking with evolutionary algorithm",
    author="Daniel Varela, Mads Jeppesen",
    author_email="daniel.varela@biochemistry.lu.se, mads.jeppesen@biochemistry.lu.se",
    url="https://github.com/Andre-lab/evodock",
    keywords="protein-docking,evolutionary-algorithms",
    install_requires=[
        "cubicsym @ git+https://github.com/Andre-lab/cubicsym@20ef715b4da5c9f0c710c7f20e41d17e66b627f2",
        "cloudcontactscore @ git+https://github.com/Andre-lab/cloudcontactscore@9f2b8429e35dea94d8e20f0bff49dbfa6171de64",
        "symmetryhandler @ git+https://github.com/Andre-lab/symmetryhandler@cdfd3c017f93101bc47be6412deaace2728e2413"
        "numpy>=1.21.0",
        "pandas>=1.3.4",
        "pillow>=9.1.0",
        "scipy>=1.7.1",
        "seaborn>=0.11.2",
        "setuptools>=44.0.0",
        "imageio>=2.10.1",
        "matplotlib>=3.4.3",
        "scikit-learn>=1.0.2",
        "biopython<=1.79",
    ],
    scripts=[
        "scripts/make_evolution_plot.py",
        "scripts/super_cleaner.py",
        "scripts/make_scatter_plot.py",
        "scripts/experiments_creator.py",
        "scripts/recover_pdbs.py",
        "scripts/scatter_plot_from_pdbs.py",
        "scripts/fake_scatter.py",
        "scripts/make_scatter_vs_rosetta.py",
        "scripts/prepacking.py",
        "scripts/symmetric_relax.py",
        "scripts/af_to_evodock.py",
    ],
    packages=[
        "src",
    ],
)
