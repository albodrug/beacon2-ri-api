#!/usr/bin/python3
# abodrug
# 2024.04.09

import streamlit as st
import requests
import utils.uifunc as uu

st.set_page_config(
    page_title="Position search",
    page_icon="🔎",
    layout="wide",
    initial_sidebar_state="expanded")


st.markdown("# Search by allele position 🔎")
st.sidebar.markdown("# Search by allele position 🔎")

# functions

#######################################################
# Do you have biosamples with variation at position ? #
#######################################################


with st.sidebar:
    # https://discuss.streamlit.io/t/\
    # how-to-use-multiple-columns-in-forms-so-that-the-input-is-side-by-side-instead-of-below-each/21035/4
    # @ yanissi

    with st.form(key='columns_in_form'):
        c1, c2, c3, c4, c5 = st.columns(5)
        with c1:
            chromosome = st.text_input("chr",value = 19, placeholder = 19 )
        with c2:
            position_start = st.text_input("start", value = 11131631, placeholder = 11131631)
        with c3:
            position_end = st.text_input("end",placeholder = 11231631)
        with c4:
            reference = st.text_input("ref", placeholder = 'G')
        with c5:
            alternate = st.text_input("alt", placeholder = 'C')

        submitButton = st.form_submit_button(label = 'Search')


callparam = {'meta': {'apiVersion': '2.0'}, 
            'query': { 
                    'requestParameters': {
                        'variantType' : "SNP", 
                        }, 
                    'filters': [], 
                    'includeResultsetResponses': 'HIT', 
                    'pagination': {'skip': 0, 'limit': 100}, 
                    'testMode': False, 
                    'requestedGranularity': 'faketyfake'
                        }
            }

if position_start:
    callparam['query']['requestParameters']['start'] = [ int(position_start) -1 ]
    callparam['query']['requestParameters']['end'] =[ int(position_start) ]
if position_end:
    callparam['query']['requestParameters']['end'] =[ int(position_end) ]
if reference:
    callparam['query']['requestParameters']['referenceBases'] = reference
if alternate:
    callparam['query']['requestParameters']['alternateBases'] = alternate    


if chromosome and position_start: 
    response_bool = uu.get_position_search_bool_response(callparam)
    response_count = uu.get_position_search_count_response(callparam)
    #response_gvar = requests.post("http://localhost:5050/api/g_variants", json=callparam)
    #jgvar = response_gvar.json()
    st.write("Variation at location exists: ", response_bool)
    st.write("Total number of variant positions: ", response_count)
    
else:
    st.write("To send an API call I need at least a chromosome number and a start position.")
