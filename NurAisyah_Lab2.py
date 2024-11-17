import streamlit as st
import pandas as pd
import networkx as nx
from io import StringIO
import requests
import matplotlib.pyplot as plt

def retrieve_ppi_biogrid(target_protein):
    biogrid_url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": "476017c234f7e7131318e91e394ff7c6",
        "format": "json",
        "searchNames": True,
        "geneList": target_protein,
        "organism": 9606,
        "includeInteractors": True
    }
    response = requests.get(biogrid_url, params=params)
    if response.status_code == 200:
        data = response.json()
        # Convert JSON response to DataFrame
        records = [(interaction['OFFICIAL_SYMBOL_A'], interaction['OFFICIAL_SYMBOL_B']) for interaction in data]
        df = pd.DataFrame(records, columns=['proteinA', 'proteinB'])
        return df
    else:
        st.error("Failed to retrieve data from BioGRID")
        return pd.DataFrame()

def retrieve_ppi_string(target_protein):
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": target_protein,
        "species": 9606
    }
    response = requests.get(string_url, params=params)
    if response.status_code == 200:
        data = response.json()
        # Convert JSON response to DataFrame
        records = [(interaction['preferredName_A'], interaction['preferredName_B']) for interaction in data]
        df = pd.DataFrame(records, columns=['proteinA', 'proteinB'])
        return df
    else:
        st.error("Failed to retrieve data from STRING")
        return pd.DataFrame()

def generate_network(dataframe):
    G = nx.Graph()
    for _, row in dataframe.iterrows():
        protein_a = row['proteinA']
        protein_b = row['proteinB']
        G.add_edge(protein_a, protein_b)
    return G

def get_centralities(network_graph):
    try:
        centralities = {
            "Degree Centrality": nx.degree_centrality(network_graph),
            "Betweenness Centrality": nx.betweenness_centrality(network_graph),
            "Closeness Centrality": nx.closeness_centrality(network_graph),
            "Eigenvector Centrality": nx.eigenvector_centrality(network_graph, max_iter=1000),
            "PageRank": nx.pagerank(network_graph)
        }
        return centralities
    except Exception as e:
        st.error(f"Error calculating centralities: {e}")
        return {}

st.title("Protein-Protein Interaction (PPI) Network Analysis")

protein_id = st.text_input("Enter Protein ID")
database = st.selectbox("Select Database", ["BioGRID", "STRING"])
retrieve = st.button("Retrieve PPI Data")

if retrieve and protein_id:
    if database == "BioGRID":
        ppi_data = retrieve_ppi_biogrid(protein_id)
    else:
        ppi_data = retrieve_ppi_string(protein_id)
    
    if not ppi_data.empty:
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("PPI Data Information")
            st.dataframe(ppi_data)
            network_graph = generate_network(ppi_data)
            
            st.write(f"**Number of Nodes**: {network_graph.number_of_nodes()}")
            st.write(f"**Number of Edges**: {network_graph.number_of_edges()}")
            
            fig, ax = plt.subplots()
            slayout = nx.spring_layout(network_graph, seed=123)
            nx.draw(network_graph, slayout, with_labels=True, node_size=20, font_size=8, ax=ax)
            plt.tight_layout()
            st.pyplot(fig)

        with col2:
            st.subheader("Centrality Measures")
            centralities = get_centralities(network_graph)
            
            for measure, values in centralities.items():
                st.write(f"**{measure}:**")
                st.write(values)
    else:
        st.warning("No data found. Please check the Protein ID.")