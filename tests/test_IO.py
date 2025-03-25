from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

def check_mystery_fasta_format(mystery_fasta: str):
    if not os.path.exists(mystery_fasta):
        raise FileNotFoundError(f"The mystery file does not exist: {mystery_fasta}")
    
    try:
        mystery_record = SeqIO.read(mystery_fasta, "fasta")
    except ValueError as e:
        raise Exception(f"Error reading mystery file: {str(e)}")
    
    assert isinstance(mystery_record, SeqRecord), "Mystery file should contain a single sequence"
    assert len(mystery_record.seq) > 0, "Mystery sequence should not be empty"
    
    print(f"Mystery file {mystery_fasta} loaded successfully with sequence ID: {mystery_record.id}")

def check_database_fasta_format(database_fasta: str):
    if not os.path.exists(database_fasta):
        raise FileNotFoundError(f"The database file does not exist: {database_fasta}")
    
    try:
        database_records = list(SeqIO.parse(database_fasta, "fasta"))
    except ValueError as e:
        raise Exception(f"Error reading database file: {str(e)}")
    
    assert len(database_records) > 1, "Database file should contain multiple sequences"
    
    for record in database_records:
        assert isinstance(record, SeqRecord), "Each record in the database should be a SeqRecord"
        assert len(record.seq) > 0, "Each sequence in the database should not be empty"
    
    print(f"Database file {database_fasta} loaded successfully with {len(database_records)} sequences")

def check_best_alignment_file():
    alignment_file = os.path.join("results", "best_alignment.txt")
    if os.path.exists(alignment_file):
        response = input(f"The file {alignment_file} already exists. Do you want to remove it? (yes/no): ")
        if response.lower() in ["yes", "y"]:
            os.remove(alignment_file)
            print(f"Old alignment file {alignment_file} removed.")
        else:
            print("Please remove or rename the old file before proceeding.")
            exit(1)

if __name__ == "__main__":
    mystery_fasta = os.path.join('data', 'mystery.fa')
    database_fasta = os.path.join('data', 'dog_breeds.fa')
    
    check_mystery_fasta_format(mystery_fasta)
    check_database_fasta_format(database_fasta)
    check_best_alignment_file()
