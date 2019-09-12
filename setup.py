import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='geostatsmodels',
    version='0.3.3',
    author='Connor Johnson',
    author_email='connor.labaume.johnson@gmail.com',
    description='Kriging and geostatistics tools.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/cjohnson318/geostatsmodels',
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)