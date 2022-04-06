from setuptools import setup

setup(
    name="evodock",
    version="1.0.0",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    description="EvoDOCK: protein-protein docking with evolutionary algorithm",
    author="Daniel Varela",
    author_email="daniel.varela@biochemistry.lu.se",
    url="https://github.com/Andre-lab/evodock",
    keywords="protein-docking,evolutionary-algorithms",
    install_requires=[
        "imageio==2.10.1",
        "matplotlib==3.4.3",
        "numpy>=1.21.0",
        "pandas==1.3.4",
        "scipy==1.7.1",
        "seaborn==0.11.2",
        "setuptools==44.0.0",
        "vector3d==1.1.1",
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
    ],
    packages=[
        "src",
    ],
)
