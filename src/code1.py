from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Align
from Bio.Align import Alignment
from typing import List, Tuple

def load_fasta_records(mystery_file: str, database_file: str) -> Tuple[SeqRecord, List[SeqRecord]]:
    """
    Loads a single sequence (mystery_record) and a list of sequences (database_records) from FASTA files.

    Parameters:
    - mystery_file (str): Path to the FASTA file containing the mystery sequence.
    - database_file (str): Path to the FASTA file containing the database of target sequences.

    Returns:
    - Tuple[SeqRecord, List[SeqRecord]]:
        - A SeqRecord object for the mystery sequence.
        - A list of SeqRecord objects for the database sequences.
    """
    mystery_record = SeqIO.read(mystery_file, 'fasta')
    database_records = list(SeqIO.parse(database_file, 'fasta'))

    return mystery_record, database_records

def find_best_alignment(mystery_record: SeqRecord, database_records: List[SeqRecord]) -> Tuple[Align.Alignment, float, str]:
    """
    Performs local pairwise alignment between a query sequence and multiple target sequences,
    returning the best local alignment based on the highest score. 
    Uses local PairwiseAligner from biopython package. 
    Best alighnment score is calculated using default settings.

    Parameters:
    - mystery_record (SeqRecord): A SeqRecord object containing the query sequence.
    - database_records (List[SeqRecord]): A list of SeqRecord objects containing target sequences.

    Returns:
    - Tuple[Align.Alignment, float, str]: 
        - Best alignment object.
        - Best alignment score.
        - Best alignment description (from the target sequence).
    """
    
    # Initialize the aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"

    query = mystery_record.seq
    best_alignments: List[Tuple[Align.Alignment, str]] = []

    # Iterate over target sequences and compute alignments
    for target in database_records:
        target_seq = target.seq
        alignments = aligner.align(target_seq, query)
        best_alignment = alignments[0]
        best_alignments.append((best_alignment, target.description))

    # Find the overall best alignment across all targets    
    search_result = max(best_alignments, key=lambda x: x[0].score)
    search_result_alignment = search_result[0]
    search_result_score = search_result_alignment.score
    search_result_description = search_result[1]

    return search_result_alignment, search_result_score, search_result_description

def save_best_alignment_result(search_result_alignment: Alignment, search_result_score: float, search_result_description: str) -> None:
    """
    Generates and saves the best alignment result to a text file.

    Parameters:
    - search_result_alignment (Alignment): The best alignment object.
    - search_result_score (float): The score of the best alignment.
    - search_result_description (str): The description of the best alignment (e.g., target sequence description).

    Returns:
    - None
    """
    
    # Generate the results string
    search_results = (f"Best Alignment Score: {search_result_score}\n"
                      f"Best Alignment Description: {search_result_description}\n"
                      f"{str(search_result_alignment)}")
        
    # Define the file path and name
    new_file = 'results\\' + 'best_alignment.txt'

    # Write the results to the file
    with open(new_file, 'w') as f_out:
        f_out.write(search_results)

    print(f"Best alignment results saved to: {new_file}")

if __name__ == "__main__":
    mystery_file = 'data\\mystery.fa'
    database_file = 'data\\dog_breeds.fa'
    
    mystery_record, database_records = load_fasta_records(mystery_file, database_file)
    search_result_alignment, search_result_score, search_result_description = find_best_alignment(mystery_record, database_records)
    save_best_alignment_result(search_result_alignment, search_result_score, search_result_description)
