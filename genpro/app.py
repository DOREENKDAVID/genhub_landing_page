#import the required modules
from flask import Flask, render_template, url_for, redirect, request, send_file
from Bio.Seq import Seq
import sqlite3
from tabulate import tabulate
from genpro.working import read_sequences_from_file, read_sequence_from_database, \
    complement, reverse_complement, transcribe, reverse_transcribe, translate_mRNA, translate_dna, \
    insert_or_update_sequences, insert_or_update_modified_sequences, \
    query_database, query_modified_seq_table, save_sequence_to_file

app = Flask(__name__)

user_sequences = None
online_sequences = None
modified_sequences = None

#this toute takes you to the  landing page
@app.route('/')
def landing_page():
    return render_template('index.html', user_sequences=user_sequences, online_sequences=online_sequences,
                           modified_sequences=modified_sequences)


#this route takes you to landing page  and if you press get started it will start the logic
@app.route('/home', methods=['POST', 'GET'])
def upload_file():
    """takes user input of file to perform central dogma"""
    global user_sequences
    if request.method == 'POST':
        file = request.files['file']
        if file:
            file_path = 'original.fasta'
            file.save(file_path)
            user_sequences = read_sequences_from_file(file_path)
            if user_sequences:
                for sequence in user_sequences:
                    accession_code = sequence.id
                    seq = str(sequence.seq)
                    insert_or_update_sequences('sequences_data.db', accession_code, seq, 'sequences')

                    seq_complement = complement(sequence.seq)
                    reverse_complement_dna = reverse_complement(sequence.seq)
                    mRNA = transcribe(sequence.seq)
                    coding_dna = reverse_transcribe(mRNA)
                    protein_seq_mRNA = translate_mRNA(mRNA)

                    # Save the modified sequences to files
                    save_sequence_to_file(sequence, accession_code, "Original Sequence", f"original_{accession_code}.fasta")
                    save_sequence_to_file(seq_complement, accession_code, "Complement Sequence", f"complement_{accession_code}.fasta")
                    save_sequence_to_file(reverse_complement_dna, accession_code, "Reverse Complement Sequence", f"reverse_complement_{accession_code}.fasta")
                    save_sequence_to_file(mRNA, accession_code, "mRNA Sequence", f"mRNA_{accession_code}.fasta")
                    save_sequence_to_file(coding_dna, accession_code, "Coding DNA Sequence", f"coding_dna_{accession_code}.fasta")
                    save_sequence_to_file(protein_seq_mRNA, accession_code, "Protein Sequence (from mRNA)", f"protein_seq_MRNA_{accession_code}.fasta")

                    insert_or_update_modified_sequences('sequences_data.db', accession_code, seq_complement, reverse_complement_dna, mRNA, protein_seq_mRNA)
    return render_template("upload.html", user_sequences=user_sequences, modified_sequences=modified_sequences)


@app.route('/query', methods=['POST'])
def query_online_db():
    """queries an online database if user does not have an input file"""

    db = request.form.get('db')
    accession_code = request.form.get('accession_code')
    sequence = None
    download_links = []
    if db and accession_code:
            sequence = read_sequence_from_database(db, accession_code)
            #if the sequence is found add to sequences table 
            insert_or_update_sequences('sequences_data.db', accession_code, str(sequence), 'sequences')
            if sequence:
                if isinstance(sequence, str) and len(sequence) > 0:
                    try:

                        sequence = Seq(sequence)
                        # Perform central dogma
                        seq_complement = complement(sequence.seq)
                        reverse_complement_dna = reverse_complement(sequence.seq)
                        mRNA = transcribe(sequence.seq)
                        coding_dna = reverse_transcribe(mRNA)
                        protein_seq_mRNA = translate_mRNA(mRNA)

                        # Save the modified sequences to files
                        save_sequence_to_file(sequence, accession_code, "Original Sequence", f"original_{accession_code}.fasta")
                        save_sequence_to_file(seq_complement, accession_code, "Complement Sequence", f"complement_{accession_code}.fasta")
                        save_sequence_to_file(reverse_complement_dna, accession_code, "Reverse Complement Sequence", f"reverse_complement_{accession_code}.fasta")
                        save_sequence_to_file(mRNA, accession_code, "mRNA Sequence", f"mRNA_{accession_code}.fasta")
                        save_sequence_to_file(coding_dna, accession_code, "Coding DNA Sequence", f"coding_dna_{accession_code}.fasta")
                        save_sequence_to_file(protein_seq_mRNA, accession_code, "Protein Sequence (from mRNA)", f"protein_seq_MRNA_{accession_code}.fasta")

                        # Insert or update modified sequences in the database
                        insert_or_update_modified_sequences('sequences_data.db', accession_code, seq_complement, reverse_complement_dna, mRNA, protein_seq_mRNA)

                        # Generate download links for users
                        download_links = {
                        "Original Sequence": f"/download/original_{accession_code}.fasta",
                        "Complement Sequence": f"/download/complement_{accession_code}.fasta",
                        "Reverse Complement Sequence": f"/download/reverse_complement_{accession_code}.fasta",
                        "mRNA Sequence": f"/download/mRNA_{accession_code}.fasta",
                        "Coding DNA Sequence": f"/download/coding_dna_{accession_code}.fasta",
                        "Protein Sequence": f"/download/protein_seq_MRNA_{accession_code}.fasta"
                    }
                        # Generate download links for users
                        for seq_type, download_link in download_links.items():
                            print(f"Download {seq_type}: <a href='{download_link}'>{seq_type}</a>")

                    except Exception as e:
                        # Log the exception or handle the error appropriately
                        print(f"Error in sequence manipulation: {e}")
                    
                else:
                    # Handle empty sequence case here
                     print("Sequence is empty or invalid")
            else:
                # Handle case where sequence is not found in the database
                print("Sequence not found in the database")
                

    return render_template("upload.html", sequence=sequence, download_links=download_links)

#this route performs central dogma on files


@app.route('/view_modified', methods=['GET', 'POST'])
def display_modified_seq():
    """this route displays the mdified sequences"""
    if request.method == 'POST':
        accession_code = request.form.get('accession_code')
        modified_sequences = query_modified_seq_table('sequences_data.db', accession_code)

        if modified_sequences:
                    
            return render_template('upload.html', modified_sequences=modified_sequences)

        else:
            return render_template('upload.html', message='No result found.')
    else:
            print("No result found.")    

   

#this route directs the user to modified files download
@app.route('/download/<filename>', methods=['GET', 'POST'])
def download(filename):
    file_name = f"{filename}"
    return send_file(file_name, as_attachment=True)


if __name__ == "__main__":
    app.run(debug=True)

