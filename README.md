# Dog Breed Project
This is a package designed to compare an input DNA sequence, compare its similarity to a database of dog breed DNA sequences, and output the following results:

1. The most similar sequence found in the dog breed database and a measure of how similar the sequence is.

2. A list of all sequences in the database with a p-value estimating probability that the input sequence is an exact match to it.

3. An adaptable phylogeny tree containing your input sequence and all sequences within the dog breed database.

Dog breed database taken from https://doi.org/10.1016/j.celrep.2017.03.079, original data available at https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000002285.5/.
This project was developed as part of the Birkbeck University Biocomputing module (Spring 2023).

# Project File Structure

1. **/data**: Holds the dog breed database files, and your mystery test sequence. Also holds a test mystery sequence for the purpose of running tests (see **/tests**)
2. **/src**: Holds source code for the project
3. **/results**: Stores the outputs from this project after running (split into **/temp** and **/final**, where /temp stores temporary backups of files during running of the project, such as extracted versions of the input file, and /final holds the final outputs as described above)
4. **/tests**: A folder for unit and integration testing of the program. Holds a .py file to run a battery of tests, and an output file recording the outcome of each test. When running the .py file, progress of each test and each outcome will also be printed to the terminal.
5. **/docs**: A folder containing documentation for each function (default and non-default arguments).

# Using this package
This package can be either called in its entirely using the command line, which will generate all outputs using all default parameters, or each function can be called individually for customisation to the defaults (including specification of where the unknown sequence is saved, and where the outputs should be saved).

Instructions are provided for producing the default 3 main outputs using either method, but please see the **/docs** folder for further information on customisation.

## (Automatic) Using this package to automatically perform all functions - alignment, statistical testing, and phylogeny outputs

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

## (Manual) Importing this package to your python environment and manually calling specific functions (please see documentation)

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

1. Navigate to the tests folder and run the unit test files using pytest (see requirements.txt file)

```bash
git clone https://github.com/sparrow-rose/DogBreed_Project.git
cd ./DogBreed_Project/
pytest ./tests/test_similarity_checker_parent_class.py
pytest ./tests/test_similarity_checker_init.py
pytest ./tests/test_similarity_checker_alignment.py
pytest ./tests/test_similarity_checker_karlin_altschul.py
pytest ./tests/test_similarity_checker_phylogeny.py
```