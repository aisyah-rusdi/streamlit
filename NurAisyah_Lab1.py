import streamlit as st
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# 1.1 Retrieve Protein Sequence Data
Entrez.email = "aisyah-@graduate.utm.my"  

def retrieve_data(protein_id):
    handle = Entrez.efetch(db='protein', id=protein_id, rettype='fasta', retmode='text')
    record = SeqIO.read(handle, 'fasta')
    handle.close()
    return record

# 1.2 Basic analyses of protein sequence
def basic_analyses(sequence):
    seq_analysis = ProteinAnalysis(record.seq)
    length = len(record.seq)
    amino_acid_comp = seq_analysis.count_amino_acids()
    molecular_weight = seq_analysis.molecular_weight()
    isoelectric_point = seq_analysis.isoelectric_point()
    
    return {
        "length": length,
        "amino_acid_composition": amino_acid_comp,
        "molecular_weight": molecular_weight,
        "isoelectric_point": isoelectric_point
    }

# 1.4 Streamlit interface
st.title('LAB 1 - NUR AISYAH BINTI RUSDI')  
protein_id = st.text_input('Enter Uniprot ID')
retrieve = st.button('Retrieve')

# your function calls and execution of processes at HERE
if retrieve:
    if protein_id:
        record = retrieve_data(protein_id)
        
        st.write(f"Raw record returned: \n{record}")
        st.write(f"Protein Name: {record.name}")
        st.write(f"Sequence: {record.seq}")
        st.write(f"Description: {record.description}")        
        
        seq_analysis = basic_analyses(record.seq)
        st.write(f"Sequence Length: {seq_analysis['length']}")
        st.write(f"Amino Acid Composition: {seq_analysis['amino_acid_composition']}")
        st.write(f"Molecular Weight: {seq_analysis['molecular_weight']}")
        st.write(f"Isoelectric Point: {seq_analysis['isoelectric_point']}")
    else:
        st.warning('Please enter a Uniprot ID')