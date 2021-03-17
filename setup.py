import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PythonFastCycles", 
    version="1.0.4",
    author="Roxane Ferry",
    author_email='roxane.ferry@ens.fr',
    description="Package to interact with FastCycles",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rferry/PythonFastCycles",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
