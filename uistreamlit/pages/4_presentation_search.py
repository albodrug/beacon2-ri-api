import streamlit as st
import utils.uifunc as uu
import pandas as pd

st.set_page_config(layout="wide")
st.subheader("Beacon API requests")

column = st.columns((2,4,4), gap="medium")

genelist = ["ANGPTL6", "SHFL", "PPAN", "SPC24", "CTSO"]

with column[0]:
    st.markdown("**Gene selection**")
    geneid = st.selectbox("", options=genelist, 
                        index=None,
                        placeholder="Select geneId...")
with column[1]:
    st.markdown("**Response by Granularity**")
    st.markdown("> **genomicVariants** endpoint")
    tab = st.tabs(["Boolean", "Count", "Record", 'json'])
    if geneid:
        with tab[0]:
            exists, response_api_boolgene = uu.get_bool_info(geneid)
            st.markdown(f":green[{exists}]")
        with tab[1]:
            count, response_api_countgene = uu.get_count_info(geneid)
            st.markdown(f":green[{count}]")
        with tab[2]:
            if exists:
                snp_allele_count = uu.get_snp_frequency(geneid)
                positions = [[int(i["position"])] for i in snp_allele_count]
                acs = [i["allele_count"] for i in snp_allele_count]
                labels = [i["variantInternalId"] for i in snp_allele_count]
                af = [i/(3202*2) for i in acs]
                data = []
                for i in range(len(positions)):
                    data.append([float(af[i]), labels[i]])
                df = pd.DataFrame(data, columns=['allele_frequency', 'variant'])
                st.write(df)
        with tab[3]:
            st.json(response_api_boolgene.json())
with column[2]:
    st.markdown("**Filters in requests**")
    st.markdown("> **individuals** endpoint: **count** granularity")
    tab = st.tabs(["_", "BMI", "Sex", "Aneurysm size", "json"])
    if geneid and exists:
        with tab[1]:
            uu.build_bmi_pie(geneid)
        with tab[2]:
            uu.build_sex_pie(geneid)
        with tab[3]:
            response = uu.build_aneurysmsize_pie(geneid)
        with tab[4]:
             st.json(response.json())