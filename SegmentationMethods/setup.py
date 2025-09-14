from setuptools import setup, find_packages

setup(
    name="crp8_sereias",
    version="0.1.0",
    description="Segmentation methods: Sereias module",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="The Cosmostatistics Initiative",
    author_email="contact@cosmostatistics-initiative.org",
    url="https://github.com/COINtoolbox/crp8_segmentation",
    license="",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "astropy"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
)

