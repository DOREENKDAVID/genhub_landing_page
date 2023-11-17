
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sqlite3
from tabulate import tabulate
import os

def read_sequences_from_file(file_path):
    """function that reads a file from user input gets the sequence ans stores in db"""
    if os.path.isfile(file_path):
        #try extracting from file and getting the accession code and sequence
        try:
            sequences = list(SeqIO.parse(file_path, "fasta"))
            if sequences:
                return sequences
            else:
                print("No sequences found in the file.")
        except Exception as e:
            print("Error:", str(e))
            return None
    else:
        print(f"File not found: {file_path}")
        return None
    
def read_sequence_from_database(db, accession_code):
    """function to get the file from a database and modify it"""
    Entrez.email = "A.N.Other@example.com"
    try:            
            with Entrez.efetch(db=db, rettype="fasta", retmode="text", id=accession_code) as handle:
                seq_record = SeqIO.read(handle, "fasta")
                return {accession_code, seq_record.description, str(seq_record.seq)}
    except Exception as e:            
            print("Error during online query:", str(e))
            return None
    
    
# Define functions for DNA/RNA/protein operations
def complement(sequence):
    seq_complement = sequence.complement()
    return seq_complement

def reverse_complement(sequence):
    reverse_complement_dna = sequence.reverse_complement()
    return reverse_complement_dna

def transcribe(sequence):
    mRNA = sequence.transcribe()
    return mRNA

def reverse_transcribe(mRNA):
    coding_dna = mRNA.back_transcribe()
    return coding_dna

def translate_mRNA(mRNA):
    codon_remainder = len(mRNA) % 3
    if codon_remainder > 0:
        mRNA += "N" * (3 - codon_remainder)

    try:
        protein_seq_mRNA = mRNA.translate(table=1, to_stop=True)
        return str(protein_seq_mRNA)
    except Exception as e:
        print(f"Translation Error: {str(e)}")
        return "Translation Error"

def translate_dna(sequence):
    codon_remainder = len(sequence) % 3
    if codon_remainder > 0:
        sequence += "N" * (3 - codon_remainder)

    try:
        protein_seq_dna = sequence.translate(table=1, to_stop=True)
        return str(protein_seq_dna)
    except Exception as e:
        print(f"Translation Error: {str(e)}")
        return "Translation Error"
    

def save_sequence_to_file(sequence, accession_code, step, filename):
    """function to save a sequence and modified sequence to a file"""

    # Check if the input is a SeqRecord or a Seq object
    if isinstance(sequence, SeqRecord):
        header = f">{accession_code} {step} {sequence.description}\n"
        sequence = str(sequence)
    else:
        header = f">{accession_code} {step}\n"
        sequence = str(sequence)
        
    with open(filename, "w") as file:
        file.write(header)
        file.write(sequence)


# Database logic
def create_database(database_name):
    """function that created the sequence db"""
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()
    # Commit and close the connection to create the database file
    conn.commit()
    conn.close()

def create_tables(database_name):
    """functions to create the tables to th db"""
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS sequences (
            accession_code TEXT PRIMARY KEY NOT NULL,
            seq_record TEXT NOT NULL
        )
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS modified_sequences (
            accession_code TEXT PRIMARY KEY NOT NULL,
            seq_complement TEXT NOT NULL,
            seq_reverse_complement TEXT NOT NULL,
            mRNA_seq TEXT NOT NULL,
            Protein_seq TEXT NOT NULL
        )
    ''')

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS users (
            user_email TEXT PRIMARY KEY,
            user_password TEXT NOT NULL,
            user_name TEXT
        )
    ''')
    
    conn.commit()
    conn.close()

def insert_or_update_sequences(database_name, accession_code, sequence, table_name):
    """function to insert the original sequence to sequences db"""
    
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()
    cursor.execute('''
    INSERT OR REPLACE INTO sequences (accession_code, seq_record)
    VALUES (?, ?)
    ''', (accession_code, sequence))
    conn.commit()
    conn.close()
    
def insert_or_update_modified_sequences(database_name, accession_code, seq_complement, reverse_complement_dna, mRNA, Protein_seq_mRNA):
    """function to insert the modified seruences to modified table"""
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()
    cursor.execute('''
    INSERT OR REPLACE INTO modified_sequences (accession_code, seq_complement, seq_reverse_complement, mRNA_seq, protein_seq)
    VALUES (?, ?, ?, ?, ?)
    ''', (accession_code, str(seq_complement), str(reverse_complement_dna), str(mRNA), str(Protein_seq_mRNA))
    )
    conn.commit()
    conn.close()

def query_database(database_name, accession_code,table_name):
    """querys the database to get sequences bassed on the accession_code"""
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()
    query = f'SELECT seq_record FROM {table_name} WHERE accession_code = ?'
    cursor.execute (query, (accession_code,))
    results = cursor.fetchone()
    conn.close()
    return results[0] if results else None

def query_modified_seq_table(database_name, accession_code):
    """querys the database to get sequences bassed on the accession_code"""
    conn = sqlite3.connect(database_name)
    cursor = conn.cursor()

    # Query the modified_sequences table
    query = f'SELECT seq_complement, seq_reverse_complement, mRNA_seq, protein_seq FROM modified_sequences WHERE accession_code = ?'
    cursor.execute (query, (accession_code,))
    modified_results = cursor.fetchone()

    # Query the sequences table to get the original sequence
    query = f'SELECT seq_record FROM sequences WHERE accession_code = ?'
    cursor.execute (query, (accession_code,))
    original_results = cursor.fetchone()
    conn.close()
    # Combine the results, the original  and modified sequence
    if modified_results and original_results:
        seq_record = original_results[0]
        seq_complement, seq_reverse_complement, mRNA_seq, protein_seq = modified_results

        return ( seq_record, seq_complement, seq_reverse_complement, mRNA_seq, protein_seq)
    else:
        None


if __name__ == "__main__":
    create_database('sequences_data.db')
    create_tables('sequences_data.db')
    file_path = input("Enter the path to a local file (or press Enter to skip this option): ")

    if file_path:
        sequences = read_sequences_from_file(file_path)
        if sequences:
            for sequence in sequences:
                accession_code = sequence.description.split('|')[1]
                seq = str(sequence.seq)
                insert_or_update_sequences('sequences_data.db', accession_code, seq, 'sequences')


    else:
        sequence = None
    #if there is sequence in file perform the central dogma
    if sequence is not None:
        seq_complement = complement(sequence.seq)
        reverse_complement_dna = reverse_complement(sequence.seq)
        mRNA = transcribe(sequence.seq)
        coding_dna = reverse_transcribe(mRNA)
        protein_seq_mRNA = translate_mRNA(mRNA)
        protein_seq_dna = translate_dna(sequence.seq)

        #save the modifies sequences to file for users to download
        save_sequence_to_file(sequence, accession_code, "Original Sequence", "original.fasta")
        save_sequence_to_file(seq_complement, accession_code, "Complement Sequence", "complement.fasta")
        save_sequence_to_file(reverse_complement_dna, accession_code, "Reverse Complement Sequence", "reverse_complement_dna.fasta")
        save_sequence_to_file(mRNA, accession_code, "mRNA Sequence", "mRNA.fasta")
        save_sequence_to_file(coding_dna, accession_code, "Coding DNA Sequence", "coding_dna.fasta")
        save_sequence_to_file(protein_seq_mRNA, accession_code, "Protein Sequence (from mRNA)", "protein_seq_mRNA.fasta")
        insert_or_update_modified_sequences('sequences_data.db', accession_code, seq_complement, reverse_complement_dna, mRNA, protein_seq_mRNA)
        accession_code = input("Enter the ID of the query record: ")
                    #query the modified database
        result = query_modified_seq_table('sequences_data.db', accession_code)

        if result:
                seq_record, seq_complement, seq_reverse_complement, mRNA_seq, protein_seq = result
                table = [["Record_ID", "sequence"], [accession_code, seq_record]]
                        
                print("Original Sequence:")
                print(tabulate(table, headers="firstrow", tablefmt="grid"))
                        
                print("Modified Sequences:")
                table = [
                        ["Sequence Type", "Sequence"],
                        ["Complement", seq_complement],
                        ["Reverse Complement", seq_reverse_complement],
                        ["mRNA", mRNA_seq], ["Protein", protein_seq],
                        ]
                        
                print(tabulate(table, headers="firstrow", tablefmt="grid"))
        else:
            print("No result found.")




    #query the online databases to get the sequences
    else:
        db = input("Enter the database (e.g., 'nucleotide'): ")
        accession_code = input("Enter the ID of the record: ")

        if db and accession_code:
            sequence = read_sequence_from_database(db, accession_code)
            if sequence:
                sequence = Seq(sequence)
                #if the sequence is found add to sequences table
                insert_or_update_sequences('sequences_data.db', accession_code, str(sequence), 'sequences')
                accession_code = input("Enter the ID of the query record: ")

                #user can query the db to get the original seq vs its modified sequeces
                result = query_database('sequences_data.db', accession_code, 'sequences')
                if result:
                    table = [["Record_ID", "sequence"], [accession_code, result]]
                    print("Result:")
                    print(tabulate(table, headers="firstrow", tablefmt="grid"))
                else:
                    print("No result found.")
                seq_complement = complement(sequence)
                reverse_complement_dna = reverse_complement(sequence)
                mRNA = transcribe(sequence)
                coding_dna = reverse_transcribe(mRNA)
                protein_seq_mRNA = translate_mRNA(mRNA)

                #save the sequences searched from the  online db to files for users to download
                save_sequence_to_file(sequence, accession_code, "Original Sequence", "original.fasta")
                save_sequence_to_file(seq_complement, accession_code, "Complement Sequence", "complement.fasta")
                save_sequence_to_file(reverse_complement_dna, accession_code, "Reverse Complement Sequence", "reverse_complement_dna.fasta")
                save_sequence_to_file(mRNA, accession_code, "mRNA Sequence", "mRNA.fasta")
                save_sequence_to_file(coding_dna, accession_code, "Coding DNA Sequence", "coding_dna.fasta")
                save_sequence_to_file(protein_seq_mRNA, accession_code, "Protein Sequence (from mRNA)", "protein_seq_mRNA.fasta")
                insert_or_update_modified_sequences('sequences_data.db', accession_code, seq_complement, reverse_complement_dna, mRNA, protein_seq_mRNA)
                accession_code = input("Enter the ID of the query record: ")
                #query the modified database
                result = query_modified_seq_table('sequences_data.db', accession_code)

                if result:
                    seq_record, seq_complement, seq_reverse_complement, mRNA_seq, protein_seq = result
                    table = [["Record_ID", "sequence"], [accession_code, seq_record]]
                    
                    print("Original Sequence:")
                    print(tabulate(table, headers="firstrow", tablefmt="grid"))
                    
                    print("Modified Sequences:")
                    table = [
                        ["Sequence Type", "Sequence"],
                        ["Complement", seq_complement],
                        ["Reverse Complement", seq_reverse_complement],
                        ["mRNA", mRNA_seq],
                        ["Protein", protein_seq],
                    ]
                    
                    print(tabulate(table, headers="firstrow", tablefmt="grid"))
                else:
                    print("No result found.")


            else:
                print("Error retrieving sequence from the database.")
        else:
            print("Invalid input. Database query requires both 'db' and 'id'.")

