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
4. **/tests**: A file for unit and integration testing of the program. Holds a .py file to run a battery of tests, and an output file recording the outcome of each test. When running the .py file, progress of each test and each outcome will also be printed to the terminal.

# Running the project

1. Clone this repository and open it as below:

```bash
git clone https://github.com/sparrow-rose/DogBreed_Project.git
cd ./DogBreed_Project
```

2. Copy your mystery sequence into the **/data/seq** folder - *please ensure your mystery sequence is in .fasta format!*

3. Run the .py file within the **/src** folder

```bash
cd ./DogBreed_Project/src/
python similarity_checker
```

4. Once the program has run, locate your results in the **/results/final** folder. 