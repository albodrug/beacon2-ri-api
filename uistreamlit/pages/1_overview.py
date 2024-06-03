#!/usr/bin/python3
# abodrug
# 2024.04.09

import streamlit as st
import requests
import plotly.graph_objects as go
import json
from pathlib import Path

# storage
storagefile = open("storage/savedvalues.json", 'r+', encoding='utf-8')
storagejson =  json.load(storagefile)
storagedata = {}

# config
apptitle = '# Overview of data in this Beacon'

st.set_page_config(
    page_title="Overview",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded")

st.markdown(apptitle)

# general
try:
    nb_vars = storagejson[Path(__file__).stem]['nb_vars']
    nb_biosamples = storagejson[Path(__file__).stem]['nb_biosamples']
    nb_individuals = storagejson[Path(__file__).stem]['nb_individuals']
except:
    response_gvar = requests.get("http://localhost:5050/api/g_variants")
    response_biosamples = requests.get("http://localhost:5050/api/biosamples")
    response_individuals = requests.get("http://localhost:5050/api/individuals")
    nb_vars = response_gvar.json()['responseSummary']['numTotalResults']
    nb_biosamples = response_biosamples.json()['responseSummary']['numTotalResults']
    nb_individuals = response_individuals.json()['responseSummary']['numTotalResults']
    #
    storagedata[Path(__file__).stem] = {}
    storagedata[Path(__file__).stem]['nb_vars'] = nb_vars
    storagedata[Path(__file__).stem]['nb_biosamples'] = nb_biosamples
    storagedata[Path(__file__).stem]['nb_individuals'] = nb_individuals
    #json.dump(storagedata, storagefile, ensure_ascii=False, indent=4)
    #st.write(storagejson)


st.write("This Beacon contains ", nb_vars, "variations", 
        ", ", nb_biosamples, " biosamples", 
        " and ", nb_individuals, " individuals.")
storagefile.close()

col = st.columns((2, 2, 2), gap='medium')

with col[1]:
    st.write("Variations")
    # variation type
    callparam = {'meta': {'apiVersion': '2.0'}, 
                'query': { 
                        'requestParameters': {
                            'geneId' : 'ANGPTL6',
                            'variantType' : 'SNP'
                            }, 
                        'filters': [], 
                        'includeResultsetResponses': 'HIT', 
                        'pagination': {'skip': 0, 'limit': 10000}, 
                        'testMode': False, 
                        'requestedGranularity': 'count'
                            }
                }
    # snp
    response_snp = requests.post("http://localhost:5050/api/g_variants", json=callparam)
    nbvar_snp = response_snp.json()['responseSummary']['numTotalResults']
    # indel
    callparam['query']['requestParameters']['variantType'] = 'INDEL'
    response_ins = requests.post("http://localhost:5050/api/g_variants", json=callparam)
    nbvar_ins = response_ins.json()['responseSummary']['numTotalResults']
    #
    labels = ["SNP", "INDEL"]
    values = [nbvar_snp, nbvar_ins]
    colors = ['#d8e5ab', '#4aa1db']
    fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.5)])
    fig.update_traces(hoverinfo='value+percent', textinfo='label', textfont_size=20,
                    marker=dict(colors=colors, line=dict(color='#000000', width=2)))
    fig.update_layout(title_text="Variant type")
    st.plotly_chart(fig, use_container_width=True) 

    # variation molecular effect
    labels = ['5_prime_UTR_variant','3_prime_UTR_variant',  'downstream_gene_variant', \
              'synonymous_variant', 'stop_gained', \
              'missense_variant', 'frameshift_variant', 'intron_variant' ]
    values = []
    for effect in labels:
        response = requests.get('http://localhost:5050/api/g_variants?filters=molecularAttributes.molecularEffects:'+effect)
        count = response.json()['responseSummary']['numTotalResults']
        values.append(count)
    
    import plotly.graph_objects as go

    fig = go.Figure(go.Bar(
            x=values,
            y=labels,
            orientation='h', 
            marker=dict(
                color='rgba(58, 71, 80, 0.6)',
                line=dict(color='rgba(58, 71, 80, 1.0)', width=3)
                ),
            ))
    fig.update_layout(title_text="Variant effect")

    st.plotly_chart(fig, use_container_width=True) 

with col[0]:
    st.write("Individuals")
    # sex distribution
    response_individuals_female = requests.get("http://localhost:5050/api/individuals?filters=NCIT:C16576")
    response_individuals_male = requests.get("http://localhost:5050/api/individuals?filters=NCIT:C20197")
    females = response_individuals_female.json()['responseSummary']['numTotalResults']
    males = response_individuals_male.json()['responseSummary']['numTotalResults']

    labels = ['females', 'males']
    values = [females, males]
    colors = ['#bca3bc', '#9a9e5e']
    fig = go.Figure(data=[go.Pie(labels=labels, values=values, hole=.5)])
    fig.update_traces(hoverinfo='label+percent', textinfo='value', textfont_size=20,
                    marker=dict(colors=colors, line=dict(color='#000000', width=2)))
    fig.update_layout(title_text="Sex distribution")
    st.plotly_chart(fig, use_container_width=True)

    # ethnicity
    #ethnicity_id = ['NCIT:C128457', 'NCIT:C161419', 'NCIT:C41263', 'NCIT:C42331', 'NCIT:C43851']
    ethnicity_label = ['American', 'East%20Asian', 'South%20Asian', 'African', 'European']
    pretty_ethnicity_label = ['American', 'East Asian', 'South Asian', 'African', 'European']
    colors = ['#55bcb5', '#e3827c', '#f48e50', '#f2f0ee', '#fac591']
    values = []
    for label in ethnicity_label:
        response = requests.get("http://localhost:5050/api/individuals?filters=ethnicity:"+label)
        count = response.json()['responseSummary']['numTotalResults']
        values.append(count)

    fig = go.Figure(data=[go.Pie(labels=pretty_ethnicity_label, values=values, hole=.5)])
    fig.update_traces(hoverinfo='value+percent', textinfo='label', textfont_size=20,
                    marker=dict(colors=colors, line=dict(color='#000000', width=2)))
    fig.update_layout(title_text="Ethnicity distribution")
    st.plotly_chart(fig, use_container_width=True)

with col[2]:
    st.write("Diseases")
