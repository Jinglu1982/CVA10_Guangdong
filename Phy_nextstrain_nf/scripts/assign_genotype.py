import argparse
from Bio import Phylo, AlignIO
import pandas as pd

def assign_genotype(sequence_file, tree_file, dataframe_file, output_file):
    # Read the dataframe file
    df = pd.read_csv(dataframe_file,sep='\t')

    # Create two dictionaries, one for storing the names and genotypes of reference sequences,
    # and another for storing the names of query sequences
    reference_genotypes = df[df['datatype'] == 'ref'].set_index('strain')['clade'].to_dict()
    query_sequences = df[df['datatype'] == 'query']['strain'].tolist()

    # Read the reference and query sequences
    alignment = AlignIO.read(sequence_file, 'fasta')

    # Construct the phylogenetic tree
    tree = Phylo.read(tree_file, 'newick')

    # For each query sequence
    for query_name in query_sequences:
        query = next(rec for rec in alignment if rec.id == query_name)
        closest_distance = float('inf')
        closest_genotype = None

        # Iterate through each leaf node in the tree
        for clade in tree.get_terminals():
            # If the leaf node is a reference sequence, we use the predefined reference genotype
            if clade.name in reference_genotypes:
                genotype = reference_genotypes[clade.name]
            else:
                continue

            # Calculate the distance between the query sequence and the leaf node
            distance = tree.distance(query.name, clade)
            if distance < closest_distance:
                closest_distance = distance
                closest_genotype = genotype

        # Update the dataframe
        df.loc[(df['strain'] == query.id) & (df['datatype'] == 'query'), 'clade'] = closest_genotype

    # Save the updated dataframe
    df.to_csv(output_file,sep='\t',index=False)
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sequence_file', required=True)
    parser.add_argument('--tree_file', required=True)
    parser.add_argument('--dataframe_file', required=True)
    parser.add_argument('--output_file', required=True)
    args = parser.parse_args()
    assign_genotype(args.sequence_file, args.tree_file, args.dataframe_file, args.output_file)
