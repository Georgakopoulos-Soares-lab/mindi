import streamlit as st
import os
import pandas as pd
import plotly.graph_objects as go
from coverage.density_plot import extract_density

# Set up layout to have side-by-side plots
st.set_page_config(layout="wide")

window_size = st.number_input("Window Size", value=500, min_value=50)

# Get list of datasets from 'extractions' directory
dataset_dir = "extractions"
datasets = [f for f in os.listdir(dataset_dir) if os.path.isfile(os.path.join(dataset_dir, f))]

# Dropdown to select a dataset
selected_dataset = st.selectbox("Choose a dataset", datasets)

gff_files = {file.name.split('.gff')[0]: file for file in Path("accessions").glob("*.gff")}
extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2])

if selected_dataset:
    extracted_file = [file for file in Path("extractions").glob("*.tsv")][0]

    chosen_file_id = extract_id(extracted_file)
    associated_gff = gff_files[chosen_file_id]
    densities = dict()

    # if the name is the same it will overwrite the class above watch out?
    for eDensity in extract_density(
                              extraction=extracted_file,
                              gff_file=associated_gff,
                              window_size=window_size
                              ):

        densities.update({eDensity.loci: eDensity.density / np.mean(eDensity.density)})

    
    # Create two line plots side by side using Plotly
    xrange = range(-window_size, window_size+1)
    fig1 = go.Figure()
    fig1.add_trace(go.Scatter(x=list(xrange), 
                        y=densities['tss'], 
                        mode='lines', 
                        line=dict(color='black')
                              )
                    )
    fig1.update_layout(title="Plot 1", xaxis_title="X-axis", yaxis_title="Y-axis")

    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(x=list(xrange), 
                              y=densities['tes'], 
                              mode='lines', 
                              line=dict(color='black')
                              )
                   )
    fig2.update_layout(title="Plot 2", xaxis_title="X-axis", yaxis_title="Y-axis")

    # Display plots side by side
    col1, col2 = st.columns(2)
    col1.plotly_chart(fig1, use_container_width=True)
    col2.plotly_chart(fig2, use_container_width=True)
