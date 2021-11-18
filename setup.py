import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
	long_description = fh.read()

setuptools.setup(
	name="GeneGrouper",
	version="1.0.2",
	author="Alexander G. McFarland",
	author_email="alexandermcfarland2022@u.northwestern.edu",
	description="Find and cluster genomic regions containing a seed gene",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/agmcfarland/GeneGrouper",
	project_urls={
		"Bug Tracker": "https://github.com/agmcfarland/GeneGrouper/issues",
	},
	classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
	],
	package_dir={"": "src"},
	packages=setuptools.find_packages(where="src", exclude=['docs','test_data']), #exlude = ['*.egg-info', ]
	entry_points={
		'console_scripts': ['GeneGrouper = GeneGrouper.__main__:main'], },
	package_data={'GeneGrouper' : ['Rscripts/*']},
	python_requires=">=3.6")

# https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments
# https://packaging.python.org/tutorials/packaging-projects/
# https://trstringer.com/easy-and-nice-python-cli/
# https://github.com/pypa/sampleproject/blob/main/setup.py
# https://chriswarrick.com/blog/2014/09/15/python-apps-the-right-way-entry_points-and-scripts/