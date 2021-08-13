
from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name="spec_map_analysis",
    version="1.1.1",
    description="An astronomical spectral map analysis package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kjdoore/spec_map_analysis",
    author="Keith Doore",
    author_email="kjdoore@gmail.com",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy"],
    license='MIT',
    packages=find_packages(),
    python_requires=">=3.6",
    install_requires=[
        'numpy>=1.20',
        'astropy>=4.2',
        'scipy>=1.6',
        'reproject>=0.7',
        'PyNeb>=1.1'],
    extras_require={
        'all_dep': ['matplotlib>=3.3',
                    'photutils>=1.1',
                    'emcee>=3.0'],
    }
)
