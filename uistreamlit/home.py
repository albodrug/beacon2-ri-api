#!/usr/bin/python3
# abodrug
# 2024.04.09

'''
UI: 1) side bar for user input 2) first column for api returns and
image rendering 3) second column for external links
'''

import streamlit as st
import requests

# background

# config
apptitle = '# Beacon2 API for iCAN-like data'

st.set_page_config(
    page_title="Home",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded")

st.markdown(apptitle)

# sidebar
#st.sidebar.markdown(apptitle)

# column api return
col = st.columns((4, 2), gap='medium')
with col[0]:
    #st.markdown('## API returns')
    response = requests.get("http://localhost:5050/api")
    j = response.json()
    r = j['response']
    #st.json(r)
    beaconname = beacondescription = r['name']
    st.markdown("## "+beaconname)
    orgname = r['organization']['name']
    orgdescription = r['organization']['description']
    orgaddress = r['organization']['address']
    beacondescription = r['description']
    st.write(beacondescription)
    st.write('@', orgname)
    st.write(orgdescription)
    st.write(orgaddress)
    
# column external info
with col[1]:
    #st.markdown('## External links')
    welcomeurl = r['organization']['welcomeUrl']
    contacturl = r['organization']['contactUrl']
    logo = r['organization']['logoUrl']
    st.write(welcomeurl)
    st.write(contacturl)
    st.image(logo)