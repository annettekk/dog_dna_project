import pytest
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from typing import List
import os

# Define the fixture for the mystery FASTA file
@pytest.fixture
def mystery_fasta():
    return os.path.join('data', 'mystery.fa')

# Define the fixture for the database FASTA file
@pytest.fixture
def database_fasta():
    return os.path.join('data', 'dog_breeds.fa')

# Test that the mystery file is a valid FASTA file and contains a single sequence
def test_mystery_fasta_format(mystery_fasta: str):
    
    from src.code import load_fasta_records
 
    # Ensure the file exists
    assert os.path.exists(mystery_fasta), f"The mystery file does not exist: {mystery_fasta}"

    # Ensure it's a valid FASTA file by attempting to parse it with SeqIO
    try:
        mystery_record = SeqIO.read(mystery_fasta, "fasta")
    except ValueError as e:
        pytest.fail(f"Error reading mystery file: {str(e)}")

    # Test that it contains only one sequence
    assert isinstance(mystery_record, SeqRecord), "Mystery file should contain a single sequence"
    assert len(mystery_record.seq) > 0, "Mystery sequence should not be empty"

    print(f"Mystery file {mystery_fasta} loaded successfully with sequence ID: {mystery_record.id}")

# Test that the database file is a valid FASTA file and contains multiple sequences
def test_database_fasta_format(database_fasta: str):
    from src.code import load_fasta_records  

    # Ensure the file exists
    assert os.path.exists(database_fasta), f"The database file does not exist: {database_fasta}"

    # Ensure it's a valid FASTA file by attempting to parse it with SeqIO
    try:
        database_records = list(SeqIO.parse(database_fasta, "fasta"))
    except ValueError as e:
        pytest.fail(f"Error reading database file: {str(e)}")

    # Test that the database contains multiple sequences
    assert len(database_records) > 0, "Database file should contain multiple sequences"
    
    # Ensure each sequence in the database is properly formatted
    for record in database_records:
        assert isinstance(record, SeqRecord), "Each record in the database should be a SeqRecord"
        assert len(record.seq) > 0, "Each sequence in the database should not be empty"
    
    print(f"Database file {database_fasta} loaded successfully with {len(database_records)} sequences")

# Test load_fasta_records function
def test_load_fasta_records(mystery_fasta: str, database_fasta: str):
    from src.code import load_fasta_records  

    # Call function with actual user files
    mystery_record, database_records = load_fasta_records(mystery_fasta, database_fasta)

    # Test if the mystery record is correctly loaded
    assert isinstance(mystery_record, SeqRecord), "Mystery record should be a SeqRecord object"
    assert len(mystery_record.seq) > 0, "Mystery sequence should not be empty"
    assert isinstance(database_records, list), "Database records should be a list"
    assert len(database_records) > 0, "Database should contain sequences"
    
    # Check individual database sequence format
    for record in database_records:
        assert isinstance(record, SeqRecord), "Each record in the database should be a SeqRecord"
        assert len(record.seq) > 0, "Each sequence in the database should not be empty"
    
    print(f"Mystery and database sequences loaded successfully.")

# Test find_best_alignment function (based on properly loaded files)
def test_find_best_alignment(mystery_fasta: str, database_fasta: str):
    from src.code import find_best_alignment  
    from src.code import load_fasta_records  

    # Load records using load_fasta_records
    mystery_record, database_records = load_fasta_records(mystery_fasta, database_fasta)

    # Call function to find the best alignment
    best_alignment, best_score, best_description = find_best_alignment(mystery_record, database_records)

    # Test that the best alignment is correct
    assert best_description in [record.description for record in database_records], "Best alignment description should be from the database"
    assert best_score >= 0, "Best alignment score should be non-negative"
    
    print(f"Best alignment found with score: {best_score} and description: {best_description}")
