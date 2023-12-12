from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='dolfen',
    version='0.0.1',
    description='Downsampling Likelihood Function Estimation',
    py_modules=['dolfen'],
    package_dir={'':'src'},
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Jethro Linley',
    author_email='jethro.linley@hotmail.co.uk',
    url='http://dolfen.readthedocs.org/',
    download_url='http://pypi.python.org/pypi/dolfen',
    install_requires=['numpy>=1.25.2','scipy>=1.11.2'],
)
