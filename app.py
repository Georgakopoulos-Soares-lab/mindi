from pathlib import Path
from collections import defaultdict
import os
import streamlit as st
import seaborn as sns
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from mindi.coverage.dna_stats import find_spacer_distribution, find_top_10
from mindi.coverage.density_plot import extract_density
from mindi.coverage.coverage_utils import extract_coverage
from mindi.coverage.pwm_density import strand_evaluators
import base64

color_palette = sns.color_palette("Set2").as_hex() + sns.color_palette("Set3").as_hex()[2:3]

x = 1
# Function to load local image and encode it to base64
def get_base64_of_image(file_path):
    with open(file_path, "rb") as f:
        return base64.b64encode(f.read()).decode()

vertical_line_style = dict(color="lightblue", 
                                width=3, 
                                dash="dot")


st.set_page_config(layout="wide")

# background_image = get_base64_of_image("dna.png")
background_image = get_base64_of_image("dna2.jpg")

# Inject custom CSS for the background image
page_bg_img = f"""
    <style>
    [data-testid="stAppViewContainer"] {{
        background-image: url("data:image/png;base64,{background_image}");
        background-size: cover;
        background-repeat: no-repeat;
        background-attachment: fixed;
        background-position: center;
    }}
    </style>
    """

extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2])
st.markdown(page_bg_img, unsafe_allow_html=True)
# Your Streamlit app content
st.title("Mindi: Your favourite NON B-DNA Visualizer")
col1, col2, col3, col4, col5 = st.columns(5)

with col1:
    window_size = st.number_input("Window Size", value=500, min_value=50)

dataset_dir = Path("mindi/extractions")
datasets = [f for f in dataset_dir.glob("**/*.tsv")]

repeats = {
        "MR": {extract_id(file): file for file in datasets if "MR" in file.name},
        "IR": {extract_id(file): file for file in datasets if "IR" in file.name},
        "RE": {extract_id(file): file for file in datasets if "RE" in file.name},
           }
modes = ["MR", "IR"] #, "RE"]
mode_colors = {"MR": {
                      "general": "#6d24d4",
                      "template": "#a586d1",
                      "non_template": "#6d24d4",
                      },
               "IR": "#3c820a", 
               "RE": "pink"}

density_types = {"MR": "template",
                 "IR": "density",
                 "RE": "density"
                 }
with col1:
    selected_dataset = st.selectbox("Choose a dataset", list(repeats["MR"].keys()))
    mode = st.selectbox("Select mode:", ["template", "density", "pwm"])

gff_files = {file.name.split('.gff')[0]: file for file in Path("mindi/accessions").glob("*.gff")}
if selected_dataset and selected_dataset in gff_files:
    chosen_accession_id = selected_dataset
    associated_gff = gff_files[chosen_accession_id]
    densities = defaultdict(dict)
    coverage_df_all = []
    spacers = {}
    top_10 = {}

    for mode in modes:
        extraction_file = repeats[mode][chosen_accession_id]
        df = pd.read_table(extraction_file)

        if mode == "IR" or mode == "MR":
            spacers[mode] = find_spacer_distribution(df)
            top_10[mode] = find_top_10(df)

            # extract coverage with genomic compartments
            coverage_df = extract_coverage(extraction_file, 
                                       associated_gff, 
                                       list(map(str, range(9))),
                                       partition_col="spacerLength",
                                       )
        else:
            # for regex case (g4)
            coverage_df = extract_coverage(extraction_file, associated_gff, ["all"], partition_col=None)

        coverage_df["mode"] = mode
        coverage_df_all.append(coverage_df)
        for eDensity in extract_density(
                                extraction=extraction_file,
                                gff_file=associated_gff,
                                density_type=density_types[mode],
                                window_size=window_size,
                                determine_strand=strand_evaluators["HDNA"]
                              ):
            densities[mode].update({(eDensity.loci, eDensity.category): eDensity.density / np.mean(eDensity.density)})

    xrange = range(-window_size, window_size+1)
    coverage_df_all = pd.concat(coverage_df_all)
    coverage_df_all["spacerLength"] = coverage_df_all["spacerLength"].astype(str)

    spacers_barplot = go.Figure()
    compartment_barplot = {}

    fig_tss = go.Figure()
    fig_tes = go.Figure()

    y1_max_tss = -float('inf')
    y1_max_tes = -float('inf')

    for mode in modes:
        # spacer barplot
        if mode == "IR" or mode == "MR":
            if mode == "IR":
                color = mode_colors[mode]
            elif mode == "MR":
                color = mode_colors[mode]["template"]
            spacers_temp = spacers[mode]
            total = np.sum(spacers_temp.values)
            spacers_barplot.add_trace(go.Bar(x=spacers_temp.index, 
                                         y=spacers_temp.values, 
                                         marker=dict(color=color),
                                         text=[f"{(val / total) * 100:.1f}%" for val in spacers_temp.values],
                                         textposition='outside', 
                                         name=mode))
        # transcription start site 
        if mode == "MR":
            tss_density_template = densities[mode]['tss', 'template']
            tss_density_nontemplate = densities[mode]['tss', 'non_template']
            fig_tss.add_trace(go.Scatter(x=list(xrange), 
                                y=tss_density_template,
                                mode='lines', 
                                name=f"{mode} Template",
                                line=dict(color=mode_colors[mode]["template"], 
                                          width=2)
                        )
                    )
            fig_tss.add_trace(go.Scatter(x=list(xrange), 
                                y=tss_density_nontemplate,
                                mode='lines', 
                                name=f"{mode} Non Template",
                                line=dict(color=mode_colors[mode]["non_template"], 
                                          width=2)
                            )
                        )
            # transcription end site
            tes_density_template = densities[mode]['tes', 'template']
            tes_density_nontemplate = densities[mode]['tes', 'non_template']
            fig_tes.add_trace(
                    go.Scatter(x=list(xrange), 
                               y=tes_density_template,
                               name=f"{mode} Template",
                               mode='lines', 
                               line=dict(color=mode_colors[mode]["template"], width=2)
                            )
                )
            fig_tes.add_trace(
                    go.Scatter(x=list(xrange), 
                               y=tes_density_nontemplate,
                               name=f"{mode} Non Template",
                               mode='lines', 
                               line=dict(color=mode_colors[mode]["non_template"], 
                                         width=2)
                            )
                    )
            y1_max_tss = max(y1_max_tss, np.max(densities[mode]['tss', 'template']), np.max(densities[mode]['tss', 'non_template']))
            y1_max_tes = max(y1_max_tes, np.max(densities[mode]['tes', 'non_template']), np.max(densities[mode]['tes', 'template']))
        else:
            tss_density = densities[mode]['tss', None]
            tes_density = densities[mode]['tes', None]
            fig_tss.add_trace(go.Scatter(x=list(xrange), 
                                y=tss_density,
                                mode='lines', 
                                name=mode,
                                line=dict(color=mode_colors[mode], width=2)
                            )
                        )
            fig_tes.add_trace(
                    go.Scatter(x=list(xrange), 
                               y=tes_density,
                               name=mode,
                               mode='lines', 
                               line=dict(color=mode_colors[mode], width=2)
                            )
                    )
            y1_max_tss = max(y1_max_tss, np.max(densities[mode]['tss', None]))
            y1_max_tes = max(y1_max_tes, np.max(densities[mode]['tes', None]))
        compartment_barplot[mode] = px.bar(coverage_df_all.query(f"mode == '{mode}'"),
                     x='compartment', 
                     y='totalCoverage', 
                     color='spacerLength', 
                     labels={"spacerLength": "Spacer Length"},
                      color_discrete_sequence=color_palette,
                     title="Subcompartment Coverage", 
                     facet_col="mode",
                     barmode='group')

    # spacers barplot
    spacers_barplot.update_layout(barmode='group', 
                                    font=dict(family="Comic Sans MS", size=14), 
                                    xaxis=dict(tickfont=dict(size=14)),
                                    xaxis_title_font=dict(size=18), 
                                    yaxis_title_font=dict(size=18),  
                                    xaxis_title='Spacer Length', 
                                    yaxis_title='Occurrences (%)'
                                  )
    spacers_barplot.update_xaxes(tickmode='linear')

    # subcompartments barplot
    for mode in modes:
        compartment_barplot[mode].update_layout(barmode='group', 
                                      xaxis_title='', 
                                      title_x=0.4,
                                      title=mode,
                                      font=dict(family="Comic Sans MS", size=14), 
                                      xaxis=dict(tickfont=dict(size=14)),
                                      xaxis_title_font=dict(size=18), 
                                      yaxis_title_font=dict(size=18),  
                                      legend_title_font=dict(size=16),  
                                      yaxis_title='Density per Mb')
        compartment_barplot[mode].update_xaxes(tickmode='linear')

    col1, col2 = st.columns(2)
    with col1:
        col1.plotly_chart(spacers_barplot) #use_container_width=True)
    
    cols = st.columns(len(modes))
    combinations = {mode: cols[i] for i, mode in enumerate(modes)}
    for mode in modes:
        with combinations[mode]:
            st.plotly_chart(compartment_barplot[mode])

    # TSS & TES
    fig_tss.add_shape(type='line', x0=0, y0=0, y1=y1_max_tss, line=vertical_line_style)
    fig_tes.add_shape(type='line', x0=0, y0=0, y1=y1_max_tes, line=vertical_line_style)

    fig_tss.update_layout(title="Transcription Start Site", 
                              xaxis=dict(tickmode='linear', 
                                         tick0=-500, 
                                         dtick=100),
                              title_x=0.4,
                              xaxis_title="", 
                              font=dict(family="Comic Sans MS", size=16), 
                              xaxis_title_font=dict(size=18), 
                              yaxis_title_font=dict(size=18),  
                              legend_title_font=dict(size=16),  
                              yaxis=dict(range=[0, None]),
                              yaxis_title="Enrichment",
                              )

    fig_tes.update_layout(
        xaxis_title='',   # No label for x-axis
        yaxis_title='Enrichment',  # Y-axis label
        xaxis=dict(tickmode='linear', tick0=-500, dtick=100),  # Increase x-ticks
        title_text='Transcription End Site',
        yaxis=dict(range=[0, None]),
        title_x=0.4, 
    )
    
    col1, col2 = st.columns(2)
    col1.plotly_chart(fig_tss, use_container_width=True, height=800)
    col2.plotly_chart(fig_tes, use_container_width=True)
