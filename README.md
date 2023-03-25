# Dog Breed Project
This is a package designed to compare an input DNA sequence, compare its similarity to a database of dog breed DNA sequences, and generate the following results:

1. The most similar sequence found in the dog breed database and a measure of how similar the sequence is.

2. A list of all sequences in the database with a p-value estimating probability that the input sequence is an exact match to it.

3. An adaptable phylogeny tree containing your input sequence and all sequences within the dog breed database.

The dog breed database taken from https://doi.org/10.1016/j.celrep.2017.03.079, original data available at https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000002285.5/.
This project was developed as part of the Birkbeck University Biocomputing module (Spring 2023).
For a list of all packages and software versions used in developing this module, please see the requirements.txt and pyproject.toml files.

# Project File Structure

1. **/data**: Holds the **/reference** dog breed database files, and your **/seq** mystery test sequence.
2. **/src**: Holds source code for the project
3. **/results**: Stores the outputs from this project after running (split into **/temp** and **/final**, where /temp stores temporary backups of files during running of the project, such as extracted versions of the input file, and /final holds the final outputs as described above)
4. **/tests**: A folder for unit and integration testing of the program. Holds a **/test_scripts** file which holds a script to test each function of the code separately, and several dummy folders representing different test cases of the /data folder. **To ensure proper running of the test scripts, it is vital these dummy folders are left untouched**. When running the .py file, progress of each test and each outcome will also be printed to the terminal.
5. **/docs**: A folder containing documentation for each function (with default and non-default arguments).

# Using this package
This package can be either called in its entirely using the command line, which will generate all outputs using all default parameters, or each function can be called individually for customisation to the defaults (including specification of where the unknown sequence is saved, and where the outputs should be saved).

Instructions are provided for producing the default 3 main outputs using either method, but please see the **/docs** folder for further information on customisation.

## Run all functions using default settings (*for generating all coursework outputs)

1. Clone this repository and open it as below:
*please ensure you have the packages listed in requirement.txt already installed*

```bash
git clone https://github.com/sparrow-rose/DogBreed_Project.git
cd ./DogBreed_Project
```

2. Copy your unknown sequence into the **/data/seq** folder - *please ensure your mystery sequence is in .fasta format!*

3. Run the similarity_checker.py file within the **/src** folder:

```bash
python ./src/similarity_checker.py
```

4. Once the program has run, locate your results in the **/results/final** folder. 

## Calling specific functions (please see documentation for further information)

1. Clone this repository and open it as below:
*please ensure you have the packages listed in requirement.txt already installed*

```bash
git clone https://github.com/sparrow-rose/DogBreed_Project.git
cd ./DogBreed_Project/
```

2. In your python environment (or from the command line), import all functions from similarity_checker

```python
from src.similarity_checker import *
```

3. Run the following functions to generate the main outputs for your unknown sequence (please see the documentation for further information, the below uses the default parameters)

```python
unknown = MysterySequence('file path to your unknown sequence')
output1 = unknown.alignment()
output2 = unknown.karlin_altschul(output1)
output3 = unknown.phylogeny()
```

# To run tests

1. Run the unit test files using pytest as below (also see requirements.txt file and ensure all packages are installed in your environment)

```bash
git clone https://github.com/sparrow-rose/DogBreed_Project.git
cd ./DogBreed_Project/
pytest ./tests/test_scripts/test_similarity_checker_parent_class.py
pytest ./tests/test_scripts/test_similarity_checker_init.py
pytest ./tests/test_scripts/test_similarity_checker_alignment.py
pytest ./tests/test_scripts/test_similarity_checker_karlin_altschul.py
pytest ./tests/test_scripts/test_similarity_checker_phylogeny.py
```