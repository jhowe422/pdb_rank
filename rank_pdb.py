import os
import numpy as np
import math
from Bio import PDB
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
import argparse

# Define favorable and allowed phi/psi regions
favored_phi = (-70, -50)
favored_psi = (-60, -30)
allowed_phi = (50, 70)
allowed_psi = (50, 70)

def extract_sequence_from_pdb(pdb_file):
    """Extracts the amino acid sequence from a PDB file using Bio.PDB."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    sequence = []

    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            for poly_index, poly_index_seq in enumerate(polypeptides):
                sequence.append(''.join([seq1(res.get_resname()) for res in poly_index_seq]))

    return ''.join(sequence)

def calculate_intensity(sequence):
    """Calculates the intensity of a sequence."""
    N = len(sequence)  # Length of the protein sequence
    
    # Frequency of each amino acid
    freq = {aa: sequence.count(aa) / N for aa in set(sequence)}
    
    # Calculate A_x(i) for each position and total intensity
    total_intensity = 0
    for i, aa in enumerate(sequence, start=1):
        f_x = freq[aa]
        Yx_i = f_x * i
        
        if Yx_i > 0:
            Ax_i = (10 * math.log(Yx_i)) / N
        else:
            Ax_i = 0
        
        total_intensity += Ax_i
    
    return total_intensity

def percent_difference(value1, value2):
    """Calculates the percent difference between two values."""
    return abs(value1 - value2) / ((value1 + value2) / 2) * 100

def calculate_phi_psi(structure):
    """Calculate phi and psi angles from a PDB structure."""
    phi_psi = []
    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            for poly in polypeptides:
                angles = poly.get_phi_psi_list()
                for i, (phi, psi) in enumerate(angles):
                    residue = poly[i].get_resname() + str(poly[i].get_id()[1])
                    if phi is not None and psi is not None:
                        phi_psi.append((residue, np.degrees(phi), np.degrees(psi)))
                    elif phi is None or psi is None:
                        print(f"Warning: Residue {residue} in chain {chain.id} has incomplete phi/psi angles.")
    return phi_psi

def evaluate_stability(angles):
    """Evaluate the stability of residues based on phi and psi angles."""
    scores = {}
    
    for residue, phi, psi in angles:
        if phi is None or psi is None:
            continue
        # Calculate distance from favored/allowed region
        distance = min(
            abs(phi - favored_phi[0]), abs(phi - favored_phi[1]),
            abs(psi - favored_psi[0]), abs(psi - favored_psi[1]),
            abs(phi - allowed_phi[0]), abs(phi - allowed_phi[1]),
            abs(psi - allowed_psi[0]), abs(psi - allowed_psi[1])
        )
        
        if (favored_phi[0] <= phi <= favored_phi[1] and
            favored_psi[0] <= psi <= favored_psi[1]):
            score = 1 - (distance / 180.0)  # 1 is the maximum score
        elif (allowed_phi[0] <= phi <= allowed_phi[1] and
              allowed_psi[0] <= psi <= allowed_psi[1]):
            score = 0.25 - (distance / 180.0)  # 0.25 is the maximum allowed score
        else:
            score = -0.5 - (distance / 180.0)  # Negative score for disallowed regions

        scores[residue] = max(score, -1)  # Ensure the score doesn't go below -1
    
    return scores

def process_pdb_files(input_fasta, pdb_directory, output_file):
    """Process PDB files to calculate combined similarity and stability scores."""
    # Read the input FASTA file
    with open(input_fasta, 'r') as f:
        lines = f.readlines()
    
    input_sequence = ''.join([line.strip() for line in lines if not line.startswith('>')])
    input_intensity = calculate_intensity(input_sequence)
    
    pdb_files = [f for f in os.listdir(pdb_directory) if f.endswith('.pdb')]
    all_scores = []
    
    for pdb_file in pdb_files:
        pdb_path = os.path.join(pdb_directory, pdb_file)
        print(f"Processing {pdb_file}...")
        
        # Extract sequence and calculate intensity
        sequence = extract_sequence_from_pdb(pdb_path)
        intensity = calculate_intensity(sequence)
        similarity_score = 100 - percent_difference(input_intensity, intensity)
        
        # Calculate phi/psi angles and evaluate stability
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_file, pdb_path)
        angles = calculate_phi_psi(structure)
        stability_scores = evaluate_stability(angles)
        stability_score = 2*sum(stability_scores.values())
        
        # Combine similarity and stability scores
        final_score = stability_score + similarity_score
        all_scores.append((pdb_file, stability_score, similarity_score, final_score))
        
        print(f"Processed {pdb_file} with final combined score {final_score}")

    # Sort scores by final combined score (greater is better)
    all_scores.sort(key=lambda x: x[3], reverse=True)
    
    # Write final scores to the output file with better spacing
    with open(output_file, 'w') as f:
        header = f"{'PDB File':<30}{'Stability Score':<30}{'Similarity Score':<30}{'Final Combined Score':<30}\n"
        f.write(header)
        f.write("="*120 + "\n")
        for pdb_file, stability, similarity, final_score in all_scores:
            line = f"{pdb_file:<30}{stability:<30}{similarity:<30}{final_score:<30}\n"
            f.write(line)

    print(f"Scoring complete. Results saved to {output_file}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate stability and similarity of PDB files.')
    parser.add_argument('input_fasta', help='Path to the input FASTA file.')
    parser.add_argument('pdb_directory', type=str, help='Directory containing PDB files')
    parser.add_argument('output_file', type=str, help='Path to the output file to save the final scores.')
    args = parser.parse_args()
    
    process_pdb_files(args.input_fasta, args.pdb_directory, args.output_file)
