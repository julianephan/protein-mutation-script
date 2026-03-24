from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import product

def generate_mutations(input_file_name, input_positions, output_file_name, max_mutations=None, position_mutation_options=None):
    """
    Parameters (inputs):
    - input_file_name: The name of the FASTA file containing the scaffold protein sequence.
    - input_positions: A list of positions (each indexed from 1) in the protein sequence where the residue will be mutated
    - output_file_name: The name of the output FASTA file that will contain all the mutated sequences.
    - max_mutations (optional): The maximum number of positions allowed to be mutated. If None, all possible mutations will be generated.
    - position_mutation_options (optional): A dictionary mapping positions to lists of allowed amino acid options.

    Outputs:
    - A FASTA file containing all the possible mutations, where each mutated sequence has a unique name.
    """
    positions = [position - 1 for position in input_positions]  # Convert to 0-based indexing
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

    # Read the scaffold protein sequence
    scaffold_record = SeqIO.read(input_file_name, "fasta")
    scaffold_sequence = str(scaffold_record.seq)

    # Generate mutation options per position
    mutation_options = []
    for pos_1based, pos_0based in zip(input_positions, positions):
        if scaffold_sequence[pos_0based] == "-":  # Account for gaps in the scaffold sequence (no mutations at gaps)
            mutation_options.append("-")
        elif position_mutation_options is not None and pos_1based in position_mutation_options:
            mutation_options.append(position_mutation_options[pos_1based])
        else:
            mutation_options.append(amino_acids)

    # List of all possible combinations of mutations at the specified positions
    all_combinations = product(*mutation_options)  
    
    # Apply mutations to the scaffold sequence and create SeqRecords for each mutant
    mutant_records = []
    mutant_index = 0  # Index for valid mutants only
    for combo in all_combinations:
        mutated_sequence = list(scaffold_sequence)
        mutation_count = 0

        # Count the number of mutations in the current combination
        for position, amino_acid in zip(positions, combo):
            if scaffold_sequence[position] == "-":
                continue  # Skip mutations at gaps
            
            if amino_acid != scaffold_sequence[position]:
                mutation_count += 1
            mutated_sequence[position] = amino_acid  # Mutate the specified position
                
        # Skip this combination if it exceeds the maximum allowed mutations
        if max_mutations is not None and mutation_count > max_mutations:
            continue  

        mutated_sequence = "".join(mutated_sequence)
        
        # Create a SeqRecord for the mutant sequence
        mutant_record = SeqRecord(
            Seq(mutated_sequence),
            id=f"scaffold_0_mutant_{mutant_index}",
            description=""
        )
        mutant_records.append(mutant_record)
        mutant_index += 1  # Increment index for valid mutants only
    
    # Write all mutated sequences to a FASTA file
    SeqIO.write(mutant_records, output_file_name, "fasta")


if __name__ == "__main__":
    # ----- CHANGE INPUTS BELOW IF NEEDED -----
    input_file_name = "./test_inputs/test.fasta"
    output_file_name = "mutants_test.fasta"

    positions = [1, 2, 3]  # List of positions to mutate (1-based indexing)
    max_mutations = None  # Optional: Maximum number of positions allowed to be mutated (set to None to allow all mutations)

    # Optional: specify amino acid options for specific positions (set to None to use all 20 amino acids for all positions)
    # For any position not specified in position_mutation_options, the default amino acid options (all 20 amino acids) will be used.
    position_mutation_options = { 
        1: list("AC"),
        2: list("DE"),
        3: list("EF")
    }  
    # -----------------------------------------

    generate_mutations(input_file_name, positions, output_file_name, max_mutations, position_mutation_options)

    # Print the number of mutants generated
    records = list(SeqIO.parse(output_file_name, "fasta"))
    print(f"Number of mutants generated: {len(records)}")