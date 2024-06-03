#!/usr/bin/python3
# abodrug
# 2024.04.09

import streamlit as st
import requests
import time
import os
import copy

from dna_features_viewer import annotate_biopython_record, BiopythonTranslator, GraphicRecord, GraphicFeature
from Bio import Entrez, SeqIO
Entrez.email = 'alexandrina.bodrug@univ-nantes.fr'
import dask.dataframe as dd
import pandas as pd
import numpy as np
import plotly.figure_factory as ff
import matplotlib.pyplot as plt

import utils.uifunc as uu

todo = ''' 
Modify the image record to really reflect the cds/mrna arrangements.
https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer/issues/85
Adjust the translation frame, as the coordinates in the translations are
set out of convenience and do not reflect the biology at all. 
Some gene names such as PLPP2 appear in beacon annotations but not in
genbank files. Solve this issue by making a representation surrounding the
variants, not taking into account the whole gene, or find a solution to
find gene synonyms (also PLPP2 is not a synonym in any genbank file genes)
Add the image rendering and the api search possible for all chromosomes.
'''

# setting up
st.set_page_config(
    page_title="Gene search",
    page_icon="🔎",
    layout="wide",
    initial_sidebar_state="expanded")

# beacon container deployed and dataset loaded

st.markdown("# Search by gene ID 🔎")
st.sidebar.markdown("# Search by gene ID 🔎")

# side bar questions and variable input
genelist = uu.get_geneIds_list()

with st.sidebar:
    st.write("Do you have biosamples with a variation in specified gene?")
    #displaygeneview = st.checkbox("Display gene structure ?")
    geneId_select = st.selectbox("", options=genelist, 
                        index=None,
                        placeholder="Select geneId...")
geneid = geneId_select #if geneId_select is not None else st.session_state.geneid_write 

with st.sidebar:
    if geneid:
        checknbindividuals = st.checkbox("Check how many individuals have a mutation in gene "+geneid+" ?")
        if checknbindividuals:
            groupsex = st.checkbox("Group by male/female ?")
        #    groupbmi = st.checkbox("Group by BMI ?")
        #    groupethno = st.checkbox("Group by ethnicity ?")
        checkvariantdistrib = st.checkbox("Check variant distribution within "+geneid+" ?")

# building response
if geneid:
    with st.spinner("Looking for samples with variation in gene "+geneid):
        time.sleep(1)
        exists, response_api_boolgene = uu.get_bool_info(geneid)
    if exists:
        st.markdown(f"🗒️Found variation in [**{geneid}**](https://www.genecards.org/cgi-bin/carddisp.pl?gene={geneid}) : :green[{exists}]")
    else:
        st.markdown(f"🗒️Found variation in [**{geneid}**](https://www.genecards.org/cgi-bin/carddisp.pl?gene={geneid}) : :red[{exists}]")
    
    expander = st.expander("See API request")
    expander.write(response_api_boolgene.json())
   
    if checknbindividuals:
        exists, nbvar, nbind, nbbio, response_api_countvar, response_api_countind = uu.get_record_info(geneid)
        st.markdown(f"🗒️:green[{nbvar}] variations found within \
                    :green[{nbind}] individuals.") 
                    #:green[{nbbio}] biosamples.")
        expander = st.expander("See API request for variants count")
        expander.write(response_api_countvar.json())
        expander = st.expander("See API request for individuals count")
        expander.write(response_api_countind.json())

        if groupsex:
            nbmale, nbfemale, response_sexcount = uu.get_individuals_by_sex(geneid)
            st.markdown(f"🗒️:green[{nbmale}] males, \
                        :green[{nbfemale}] females.")
            import plotly.express as px
            import pandas as pd 
            d = {'sex': ['female', 'male'], 'count': [nbfemale, nbmale]}
            df = pd.DataFrame(data=d)
            fig = px.pie(df, values='count', names="sex", hole=.4)
            st.plotly_chart(fig)
            expander = st.expander("See API request for female individuals with variants in "+geneid)
            expander.write(response_sexcount.json())
        #if groupbmi:
        #    data, response_bmicount = uu.get_individuals_by_bmi(geneid)
        #    st.write(data)
        #    expander = st.expander("See API request for BMI 15+ individuals with variants in "+geneid)
        #    expander.write(response_bmicount.json())

        #
        if checkvariantdistrib:
            snp_allele_count = uu.get_snp_frequency(geneid)
            #positions_2d = [[i["position"]]*i["allele_count"] for i in snp_allele_count]
            #positions = [item for sub_list in positions_2d for item in sub_list]
            positions = [[int(i["position"])] for i in snp_allele_count]
            acs = [i["allele_count"] for i in snp_allele_count]
            labels = [i["variantInternalId"] for i in snp_allele_count]
            #hist_data = [positions]
            #group_labels = ["SNP"]
            #fig = ff.create_distplot(
            #        hist_data, group_labels, bin_size=[1])
            #st.write("SNP distribution")     
            #df = pd.DataFrame.from_dict({'position' : positions, 'allele_count': acs, 'labels': labels})
            af = [i/(3202*2) for i in acs]
            #chart_data =  pd.DataFrame({'position': positions, 'allele_frequency': af, 'labels' :  labels})  

            #st.scatter_chart(
            #    chart_data,
            #    x='position',
            #    y='allele_frequency'
            #    )
            data = []
            for i in range(len(positions)):
                data.append([int(positions[i][0]), float(af[i]), int(acs[i]), labels[i]])
            df = pd.DataFrame(data, columns=['positions', 'allele_frequency', 'allele_count', 'variant'])
            st.markdown(f"🗒️Check [known variants in **{geneid}**](https://www.genecards.org/cgi-bin/carddisp.pl?gene={geneid}#snp)")
            st.write(df)

            import plotly.graph_objects as go

            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=df['positions'],
                y=df['allele_frequency'],
                marker=dict(color="gray", size=5),
                mode="markers",
                name="allele frequency",
            ))
            st.plotly_chart(fig)

            chromosome = 'hg38.chr19'
            chromosome_gb_id = uu.get_chromosome_gb_id(chromosome)
            # get gene coordinates
            genefeature = uu.get_gene_genbank_info(geneid)
            #st.write(genefeature)
            gene_start_coordinate = genefeature.location.start
            gene_end_coordinate = genefeature.location.end
            varcoord = [gene_start_coordinate-50, gene_end_coordinate+50]
            try:
                genegbfile = "temp/tmp_"+geneid+".gb"
                if True == True : #not os.path.exists(genegbfile):
                    # efetch -db nuccore -id NC_000962.3 -format gb -seq_start 1 -seq_stop 99000 > seq.gb
                    gbstart, gbend = uu.get_gbcoordinates_to_efetch(gene_start_coordinate, gene_end_coordinate, varcoord)
                    handle = Entrez.efetch(db="nuccore", id=chromosome_gb_id, rettype="gb", retmode="text", \
                                        seq_start=gbstart, seq_stop=gbend)
                    record = SeqIO.read(handle, "genbank")
                    # loop over variants and variant labels
                    #for i in range(len(varcoord)):
                    #    coordinateduo = varcoord[i]
                    #    varlabel = varmeta[i]
                    #    coordinate_in_gene_start = coordinateduo[0] - gene_start_coordinate
                    #    coordinate_in_gene_end = coordinateduo[1] - gene_start_coordinate
                        # annotate genbank file with labels of variant at variant location
                        #annotate_biopython_record(
                        #    record, location=(coordinate_in_gene_start, coordinate_in_gene_end), 
                        #    feature_type='variant', label=varlabel.replace("chrchr19_", "").replace("_", ":", 1).replace("_", ">")
                        #)
                    # loop over features
                    comment = '''
                    c = 0
                    for feature in record.features:
                        if len(feature.location.parts) != 1:  # implicit check, compound locations will have multiple parts
                            for i in range(0, len(feature.location.parts)):
                                part = feature.location.parts[i]
                                new_feature = copy.deepcopy(feature)  # otherwise overwrites the original feature
                                new_feature.location = copy.deepcopy(part)
                                new_feature.type = 'my-'+str(c)+"_"+feature.type+"_"+str(i)
                                record.features += [new_feature]
                            c = c + 1
                    '''
                    # delete certain features
                    record.features = [
                        f for f in record.features if f.type not in ["CDS", "mRNA"] and "mRNA" not in f.type
                    ]
        
                    SeqIO.write(record, genegbfile, "genbank")
                    handle.close()
                    # add variants as gb features
                graphic_record = uu.MyCustomTranslator().translate_record(genegbfile)    
                ax, _ = graphic_record.plot(figure_width=10, strand_in_label_threshold=7)
                #
                genegbimage = "images/"+geneid+"_genbank.png"
                ax.set_title("", loc='left', weight='bold')
                ax.figure.savefig(genegbimage)
                st.write("Gene structure")
                st.image(genegbimage)
            except Exception as error:
                st.write("Gene view rendering impossible.")
                st.write(error)
            
    exit()


        # propose to display gene view only if checked nb individuals previously
    if True == False: #displaygeneview:
        # get chromosome id for genbank
        chromosome = 'hg38.chr19'
        chromosome_gb_id = uu.get_chromosome_gb_id(chromosome)
        # get gene coordinates
        genefeature = uu.get_gene_genbank_info(geneid)
        #st.write(genefeature)
        gene_start_coordinate = genefeature.location.start
        gene_end_coordinate = genefeature.location.end
        try:
            genegbfile = "temp/tmp_"+geneid+".gb"
            if True == True : #not os.path.exists(genegbfile):
                # efetch -db nuccore -id NC_000962.3 -format gb -seq_start 1 -seq_stop 99000 > seq.gb
                gbstart, gbend = uu.get_gbcoordinates_to_efetch(gene_start_coordinate, gene_end_coordinate, varcoord)
                handle = Entrez.efetch(db="nuccore", id=chromosome_gb_id, rettype="gb", retmode="text", \
                                    seq_start=gbstart, seq_stop=gbend)
                record = SeqIO.read(handle, "genbank")
                # loop over variants and variant labels
                #for i in range(len(varcoord)):
                #    coordinateduo = varcoord[i]
                #    varlabel = varmeta[i]
                #    coordinate_in_gene_start = coordinateduo[0] - gene_start_coordinate
                #    coordinate_in_gene_end = coordinateduo[1] - gene_start_coordinate
                    # annotate genbank file with labels of variant at variant location
                    #annotate_biopython_record(
                    #    record, location=(coordinate_in_gene_start, coordinate_in_gene_end), 
                    #    feature_type='variant', label=varlabel.replace("chrchr19_", "").replace("_", ":", 1).replace("_", ">")
                    #)
                # loop over features
                comment = '''
                c = 0
                for feature in record.features:
                    if len(feature.location.parts) != 1:  # implicit check, compound locations will have multiple parts
                        for i in range(0, len(feature.location.parts)):
                            part = feature.location.parts[i]
                            new_feature = copy.deepcopy(feature)  # otherwise overwrites the original feature
                            new_feature.location = copy.deepcopy(part)
                            new_feature.type = 'my-'+str(c)+"_"+feature.type+"_"+str(i)
                            record.features += [new_feature]
                        c = c + 1
                '''
                # delete certain features
                record.features = [
                    f for f in record.features if f.type not in ["CDS", "mRNA"] and "mRNA" not in f.type
                ]
    
                SeqIO.write(record, genegbfile, "genbank")
                handle.close()
                # add variants as gb features
            graphic_record = uu.MyCustomTranslator().translate_record(genegbfile)    
            ax, _ = graphic_record.plot(figure_width=10, strand_in_label_threshold=7)
            #
            genegbimage = "images/"+geneid+"_genbank.png"
            ax.set_title("", loc='left', weight='bold')
            ax.figure.savefig(genegbimage)
            st.write("Gene structure")
            st.image(genegbimage)
        except Exception as error:
            st.write("Gene view rendering impossible.")
            st.write(error)
                
comment = """
                record = uu.MyCustomTranslator().translate_record(genegbfile)
                from faker import Factory
                fake = Factory.create()
                for i in range(len(varcoord)):
                    coordinateduo = varcoord[i]
                    coordinate_in_gene = coordinateduo[1] - gene_start_coordinate
                    s = coordinate_in_gene - 25
                    e = coordinate_in_gene + 25
                    tsvet = fake.hex_color()
                    #st.write(coordinateduo, coordinate_in_gene, s, e, gbstart, gbend)
                    ax.fill_between((s, e), +1000, -1000, alpha=0.25, color = tsvet)
                    zoom_start, zoom_end = s, e  # coordinates of the "detail"
                    #st.write(zoom_start, zoom_end)
                    try:
                        cropped_record = record.crop((zoom_start, zoom_end))
                        ax_cropped, _ = cropped_record.plot(figure_width=10, strand_in_label_threshold=7)
                        # TRANSLATION COORDINATES MUST START AT CDS
                        cropped_record.plot_translation(ax=ax_cropped, location=(zoom_start, zoom_end),
                                fontdict={'weight': 'bold'})
                        ax_cropped.fill_between((coordinate_in_gene-1, coordinate_in_gene+1), 
                                                +1000, -1000, alpha=0.25, color=tsvet)
                        ax_cropped.set_title("Variation detail "+varmeta[i], loc='left', weight='bold')
                        variationdetailimage = "images/"+geneid+"_"+str(coordinate_in_gene)+"_"+str(count)+"_genbank.png"
                        ax_cropped.figure.savefig(variationdetailimage)
                        st.image(variationdetailimage)
                    except Exception as error:
                        st.write("variant view rendering impossible: ", coordinate_in_gene)
                        st.write(error)
                    count = count + 1
                """
            
