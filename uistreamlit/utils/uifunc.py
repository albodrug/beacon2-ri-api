
from Bio import SeqIO, Entrez, GenBank
from dna_features_viewer import BiopythonTranslator
import requests
import streamlit as st
import os

verbose = 0

# classes
class MyCustomTranslator(BiopythonTranslator):
    """Custom translator implementing the following theme:

    - Color terminators in green, CDS in blue, all other features in gold.
    - Do not display features that are restriction sites unless they are BamHI
    - Do not display labels for restriction sites.
    - For CDS labels just write "CDS here" instead of the name of the gene.

    @ dna_feature_viewer examples

    """

    def compute_feature_color(self, feature):
        if "CDS" in feature.type:
            return "#e377c2" # cooked asparagus green
        elif "mRNA" in feature.type:
            return "#2ca02c" # raspberry yohurt pink
        elif feature.type == "variant":
            return "gray"
        elif feature.type == "gene":
            return "#17becf"  # blue teal
        else:
            return "gray"

    def compute_feature_label(self, feature):
        if feature.type == 'CDS' or feature.type == 'mRNA':
            return None
        elif feature.type == 'gene':
            return BiopythonTranslator.compute_feature_label(self, feature)
        elif "my-" in feature.type:
            if "CDS" in feature.type:
                return None
            elif "mRNA" in feature.type:
                return None
            else:
                return None #BiopythonTranslator.compute_feature_label(self, feature)
        else:
            return None #BiopythonTranslator.compute_feature_label(self, feature)


# functions
def get_geneIds_list():
    ''' Get geneId list from API?
    This is for test purposes. On a real whole dataset, this would take
    a very long time. Parsing all the geneIds where variation is found is
    not that smart. I'm extracting it from the raw json file. Maybe a smarter
    way to get the list will come up.
    '''
    itxgenes = []
    '''
    callparam = {'meta': {'apiVersion': '2.0'}, 
                'query': { 
                        'requestParameters': {
                            "geneId": "CPAMD8"
                            }, 
                        'filters': [], 
                        'includeResultsetResponses': 'HIT', 
                        'pagination': {'skip': 0, 'limit': 10000}, 
                        'testMode': False, 
                        'requestedGranularity': 'record'
                            }
                }
    response_gvar = requests.post("http://localhost:5050/api/g_variants", json=callparam)
    try:
        jgvar = response_gvar.json()
        st.json(jgvar)
    except:
        itxgenes = ['LDLR']
    '''

    #bodrug-a ▶ pp-irs1-4071ylt ▶~/Devlopment/beacon2-ri-api/deploy 
    #(master)> sed 's,.*geneIds,geneIds,' data/1kGP-hg38-chr19-exome/genomicVariations.json | \
    #sed -e 's/],.*//' -e 's,]}.*,,' | sed -e 's/geneIds":\["//' -e 's/[",]/ /g' | uniq | xargs -n 1 \
    # | uniq | sort | uniq |awk '{printf "\""$1"\", "}' 
    genesinjson = ["AC008764.1", "AC020931.1", "AC026803.2", "AC092069.1", "ANGPTL6", "BNIP3P8", "C19orf44", "CALR", "CALR3", "CDC42EP3P1", "CDC42EP5", "COLGALT1", "CPAMD8", "DYRK1B", "F2RL3", "FARSA", "FBL", "FTL", "GYS1", "LDLR", "LENG8", "LENG9", "MIR6515", "MIR6719", "MIR6886", "NIBAN3", "PPAN", "PPAN-P2RY11", "RAD23A", "RUVBL2", "SHFL", "SNORD105", "SPC24", "ZNF91"]

    genelist = list(set([*itxgenes, *genesinjson]))
    genelist.sort()
    return genelist
#
def get_record_info(geneid):
    ''' Getting counts of variation, individuals and biosamples.
    '''
    callparam = {'meta': {'apiVersion': '2.0'}, 
                'query': { 
                        'requestParameters': {
                            'geneId': geneid,
                            'variantType' : 'SNP'
                            }, 
                        'filters': [], 
                        'includeResultsetResponses': 'HIT', 
                        'pagination': {'skip': 0, 'limit': 10000}, 
                        'testMode': False, 
                        'requestedGranularity': 'count'
                            }
                }
    exists = False
    nbvar = -1
    nbind = -1
    nbbio = -1
    response_gvar = requests.post("http://localhost:5050/api/g_variants", json=callparam)
    jgvar = response_gvar.json()
    #st.write(jgvar)
    nbvar = jgvar["responseSummary"]["numTotalResults"]
    if verbose: st.json(jgvar) 
    # getting nb of biosamples from g_variations end point
    # because the biosamples endpoint count does not work (returns 15)
    #nbbio = get_nb_biosamples_from_g_variants(jgvar)
    #variant_coordinates, varmeta = get_variant_coordinates(callparam)
    try:
        #
        callparam_ind = callparam
        callparam_ind["query"]["filters"] = []
        response_individuals = requests.post("http://localhost:5050/api/individuals", json=callparam_ind)
        jind = response_individuals.json()
        #st.json(jind)
        nbind = jind["responseSummary"]["numTotalResults"]
        #
        response_biosamples = requests.post("http://localhost:5050/api/biosamples", json=callparam)
        jbio = response_biosamples.json()
        nbbio = jbio["responseSummary"]["numTotalResults"]
    except:
        pass
    return exists, nbvar, nbind, nbbio, response_gvar, response_individuals #, variant_coordinates, varmeta
#
def get_variant_coordinates(callparam):
    ''' Get list of coordinates at which variation occurs.
    '''
    list = []
    meta = []
    callparam["query"]["requestedGranularity"] =  "record"
    response_gvar = requests.post("http://localhost:5050/api/g_variants", json=callparam)
    jgvar = response_gvar.json()
    variants = jgvar['response']['resultSets'][0]['results']
    for el in variants:
        #st.json(el)
        position = el['_position']
        start = position['startInteger']
        end = position['endInteger']
        list.append([start, end]) # list of lists
        #
        identifier = el['variantInternalId']
        meta.append(identifier)
    #st.write(meta)
    return list, meta
#
def get_bool_info(geneid):
    ''' Just return boolean info about existence of variation in a gene.
    '''
    callparam = {'meta': {'apiVersion': '2.0'}, 
            'query': { 
                    'requestParameters': {
                        'geneId': geneid,
                        'variantType' : 'SNP'
                        }, 
                    'filters': [], 
                    'includeResultsetResponses': 'HIT', 
                    'pagination': {'skip': 0, 'limit': 10000}, 
                    'testMode': False, 
                    'requestedGranularity': 'boolean'
                        }
            }
    response_gvar = requests.post("http://localhost:5050/api/g_variants", json=callparam)
    jgvar = response_gvar.json()
    response = jgvar["responseSummary"]["exists"]
    return response, response_gvar
#
def get_chromosome_gb_id(chrom):
    '''
    '''
    dict = {"hg19.chr19" : "NC_000019.9", 
            "hg38.chr19" : "NC_000019.10"}
    return dict[chrom]
#
def get_gbcoordinates_to_efetch(c1, c2, l):
    ''' among gene coordiantes extracted from gtf
    and variation coordinates found with api
    retrieve smallest and largest coordinates
    to pass to entrez efetch for gb file extraction
    '''
    all = []
    for row in l:
        all.extend(row)
    all.append(c1)
    all.append(c2)
    all.sort()
    start = all[0] 
    end = all[-1] 
    return start, end
#
def get_nb_biosamples_from_g_variants(j):
    ''' Obtaining number of distinct biosamples where 
    a variation is present in a speficied gene.
    Beware: all biosamples will be counted, meaning those
    that belong to the dataset and those that do not but that
    were present in the vcf.
    '''
    biosamples = []
    try:
        results = j['response']['resultSets'][0]['results']
    except:
        return 0
    for el in results:
        for info in el['caseLevelData']:
            biosampleId = info["biosampleId"]
            biosamples.append(biosampleId)
    ubiosamples = list(set(biosamples))
    nbubiosamples = len(ubiosamples)
    #st.write(j)
    return nbubiosamples
#
def get_gene_genbank_info(geneid):
    ''' Get gene information such as start and end from
    genbank files.
    '''
    chromosome_gb_id = get_chromosome_gb_id('hg38.chr19')
    chromgbfile = "temp/tmp_"+chromosome_gb_id+".gb"

    if not os.path.exists(chromgbfile): # download the chromosome gb file if it does not exist already
        handle = Entrez.efetch(db="nuccore", id=chromosome_gb_id, rettype="gb", retmode="text", \
                                        seq_start=1, seq_stop=249250621)
        feature_table = handle.read()
        with open(chromgbfile, "w") as f:
                print(feature_table, file=f)
        handle.close()

    for rec in SeqIO.parse(chromgbfile, "gb"):
        for feature in rec.features:
            if feature.type == 'gene':
                for key, val in feature.qualifiers.items():
                    if key == 'gene' or key == 'gene_synonym':
                        for id in val:
                            if id == geneid:                               
                                return feature
            #elif feature.type == 'CDS':
            #    st.write(feature.type, feature.location)
            #    st.write('feature start: ', feature.location.start, feature.location.end)
    return None
#
def get_position_search_bool_response(callparam):
    '''
    '''
    callparam['query']['requestedGranularity'] = 'boolean'
    response = requests.post("http://localhost:5050/api/g_variants", json=callparam)
    jgvar = response.json()
    exists = jgvar["responseSummary"]["exists"]
    return exists
#
def get_position_search_count_response(callparam):
    '''
    '''
    callparam['query']['requestedGranularity'] = 'count'
    response = requests.post("http://localhost:5050/api/g_variants", json=callparam)
    jgvar = response.json()
    numTotalResults = jgvar["responseSummary"]["numTotalResults"]
    return numTotalResults
#
def get_position_search_reconrd_response(callparam):
    '''
    '''
    return 0
#
def get_individuals_by_sex(geneid):
    ''' An example of request with filters
    '''
    callparam = {'meta': {'apiVersion': '2.0'}, 
                'query': { 
                    'requestParameters': {
                        'geneId': geneid,
                        'variantType' : 'SNP'
                        }, 
                    'filters': [], 
                    'includeResultsetResponses': 'HIT', 
                    'pagination': {'skip': 0, 'limit': 10000}, 
                    'testMode': False, 
                    'requestedGranularity': 'count'
                        }
            }
    # males
    callparam["query"]["filters"] = [{"id": "NCIT:C20197", "label" : "male"}]
    response_individuals_males = requests.post("http://localhost:5050/api/individuals", json=callparam)
    jind_males = response_individuals_males.json()
    nbmales = jind_males["responseSummary"]["numTotalResults"]
    # females
    callparam["query"]["filters"] = [{"id": "NCIT:C16576", "label" : "female"}]
    response_individuals_females = requests.post("http://localhost:5050/api/individuals", json=callparam)
    jind_females = response_individuals_females.json()
    nbfemales = jind_females["responseSummary"]["numTotalResults"]
    return nbmales, nbfemales, response_individuals_females
#
def get_individuals_by_bmi(geneid):
    ''' An example of request with filters
    '''
    callparam = {'meta': {'apiVersion': '2.0'}, 
                'query': { 
                    'requestParameters': {
                        'geneId': geneid,
                        'variantType' : 'SNP'
                        }, 
                    'filters': [], 
                    'includeResultsetResponses': 'HIT', 
                    'pagination': {'skip': 0, 'limit': 10000}, 
                    'testMode': False, 
                    'requestedGranularity': 'record'
                        }
            }
    # bmi > 15
    callparam["query"]["filters"] = [{"id": "LOINC:35925-4", "label" : "BMI"}]
    response_bmi15 = requests.post("http://localhost:5050/api/individuals", json=callparam)
    #
    jind_bmi15 = response_bmi15.json()
    bmi15 = jind_bmi15["responseSummary"]["numTotalResults"]
    return bmi15, response_bmi15
#
def get_individuals_by_enthinicty(geneid):
    ''' individuals endpoint
    '''
    callparam = {'meta': {'apiVersion': '2.0'}, 
            'query': { 
                'requestParameters': {
                    'geneId': geneid,
                    'variantType' : 'SNP'
                    }, 
                'filters': [], 
                'includeResultsetResponses': 'HIT', 
                'pagination': {'skip': 0, 'limit': 10000}, 
                'testMode': False, 
                'requestedGranularity': 'count'
                    }
        }
    callparam["query"]["filters"] = [{"id": "NCIT:C43851", "label" : "European" }]
    return {}
#
def get_snp_frequency(geneid):
    '''
    '''
    callparam = {'meta': {'apiVersion': '2.0'}, 
            'query': { 
                'requestParameters': {
                    'geneId': geneid,
                    'variantType' : 'SNP'
                    }, 
                'filters': [], 
                'includeResultsetResponses': 'HIT', 
                'pagination': {'skip': 0, 'limit': 10000}, 
                'testMode': False, 
                'requestedGranularity': 'record'
                    }
        }
    response = requests.post("http://localhost:5050/api/g_variants", json=callparam)
    jgvar = response.json()
    results = jgvar['response']['resultSets'][0]['results']
    AF = []
    for el in results:
        position = el["_position"]['endInteger']
        internalId = el["variantInternalId"]
        allele_count = 0
        for biosampleinfo in el['caseLevelData']:
            biosampleId = biosampleinfo["biosampleId"]
            zygosity = biosampleinfo["zygosity"]['id']
            if zygosity == "GENO:GENO_0000458":
                allele_count = allele_count + 1
            elif zygosity == "GENO:GENO_0000136":
                allele_count = allele_count + 2
        AF.append({'position' : position, 'allele_count' : allele_count, 'variantInternalId' : internalId.replace("chr", "")})
              
    return AF
