#import the required modules
from flask import Flask, render_template, url_for, redirect, request, send_file
from Bio.Seq import Seq
import sqlite3
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
            file_path = 'temp.fasta'
            file.save(file_path)
            user_sequences = read_sequences_from_file(file_path)
    return render_template("upload.html", user_sequences=user_sequences, modified_sequences=modified_sequences)


@app.route('/query', methods=['POST'])
def query_online_db():
    """queries an online database if user does not have an input file"""
    global online_sequences
    db = request.form.get('db')
    accession_code = request.form.get('accession_code')
    online_sequences = read_sequence_from_database(db, accession_code)
    return render_template("upload.html", online_sequences=online_sequences, modified_sequences=modified_sequences)

#this route performs central dogma on files
@app.route('/perform_central_dogma', methods=['POST'])
def perform_central_dogma():
    global user_sequences, online_sequences, modified_sequences

    if request.method == 'POST':
        if 'file' in request.files:
            file = request.files['file']
            sequences = read_sequences_from_file(file)

            if sequences:
                modified_sequences = []

                for sequence in sequences:
                    accession_code = sequence.description.split('|')[1]
                    seq = str(sequence.seq)
                    insert_or_update_sequences('sequences_data.db', accession_code, seq, 'sequences')

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
                    insert_or_update_modified_sequences('sequences_data.db', accession_code, seq_complement, reverse_complement_dna, mRNA, protein_seq_mRNA)

                    # Append the modified sequences for display
                    modified_sequences.append((accession_code, seq_complement, reverse_complement_dna, mRNA, protein_seq_mRNA))

                return render_template('upload.html', user_sequences=user_sequences, online_sequences=online_sequences, modified_sequences=modified_sequences)
            else:
                return render_template('upload.html', error="No sequences found in the file.")

        elif 'db' in request.form and 'accession_code' in request.form:
            db = request.form['db']
            accession_code = request.form['accession_code']

            sequence = read_sequence_from_database(db, accession_code)

            if sequence:
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
                insert_or_update_modified_sequences('sequences_data.db', accession_code, seq_complement, reverse_complement_dna, mRNA, protein_seq_mRNA)

                modified_sequences = [(accession_code, seq_complement, reverse_complement_dna, mRNA, protein_seq_mRNA)]

                return render_template('upload.html', user_sequences=user_sequences, online_sequences=online_sequences, modified_sequences=modified_sequences)
            else:
                return render_template('upload.html', error="Error retrieving sequence from the database.")

    return render_template('upload.html', user_sequences=user_sequences, online_sequences=online_sequences, modified_sequences=modified_sequences)



#this route directs the user to modified files download
@app.route('/download/<filename>', methods=['GET', 'POST'])
def download(filename):
    """download the modified files"""
    download_file = save_sequence_to_file()
    return send_file(download_file, as_attachment=True)


if __name__ == "__main__":
    app.run(debug=True)
