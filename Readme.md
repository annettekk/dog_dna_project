Welcome to https://github.com/annettekk/dog_dna_project

This project allows you to find the nearest DNA sequence from the a list of DNA sequences in the database file.
Both the sequence-of-interest and the database file should be in fasta format.
The output is description of the best match taken from its fasta record, best alighnment score calculated using default settings of the local PairwiseAligner from biopython package, and the alignment itself.
Attention: this tool does not align to the reverse strands of sequences in the database.
Files can be DNA sequences from any organism, not only Canis lupus familiaris, nevermind the name of the project!

This is an open source code. Enjoy!

## Installation

Before running this project, ensure that you have Python installed. Then, install the required dependencies:

```sh
pip install -r requirements.txt

```
Don't forget to dowload dog_dna_project folder to your computer.

## Input files

Replace mystery.fa and dog_breeds.fa in the 'data' folder with your sequence-of-interest and database files.
Both files should be in fasta format.

## Running the code

In terminal, run code1.py from dog_dna_project folder:

```sh
python ./src/code1.py

```

## Results

The best alignment results are saved as txt file in 'results' folder.

## Testing and debugging

You can test whether the input files that you provided are in the correct format.
In terminal, run test_IO.py from dog_dna_project folder:

```sh
python ./tests/test_IO.py

```
The test file also helps you delete old results file, in case you run code1.py multiple times.