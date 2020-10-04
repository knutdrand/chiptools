import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name='chiptools',  
    version='0.1',
    # scripts=['chiptools/main.py'] ,
    author="Knut Rand",
    author_email="knutdrand@gmail.com",
    description="Package for working with chip-seq data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/knutdrand/chiptools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
     ],
    entry_points = {
        'console_scripts': ['chiptools=chiptools.main:main',
                            'chipplots=chiptools.plot_cli:main',
        ],
    },
    install_requires=['trackhub', 'numpy', 'cigar', 'matplotlib']
 )
