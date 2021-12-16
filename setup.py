from setuptools import setup, find_packages

setup(
    name='evodock',
    version='1.0.0',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        'imageio==2.11.1',
        'matplotlib==3.3.3',
        'numpy==1.19.4',
        'pandas==1.1.4',
        'scipy==1.6.1',
        'seaborn==0.11.1',
        'vector3d==1.1.1',
    ],
    scripts=[
        'scripts/make_evolution_plot.py',
        'scripts/make_scatter_plot.py',
        'scripts/scatter_plot_from_pdbs.py',
        'scripts/make_scatter_vs_rosetta.py',
        'scripts/prepacking.py',
    ],
    packages=find_packages(include=['src', 'src.*']),
)
