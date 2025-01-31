from setuptools import setup, find_packages

setup(
    name="engine",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        # Add any required dependencies here, for example:
        "rdkit",
        "pytest",
        "biopython",
        "matplotlib",
        "pypdb",
        "redis",
        "google-genai",
        "google-cloud-aiplatform",
        "google-cloud-storage",
        "sentence-transformers",
        "pandas"
    ],
    author="TRY-ER",
    author_email="",
    description="",
)