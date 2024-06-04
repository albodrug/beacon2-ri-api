import yaml

with open("beacon/api_version.yml") as api_version_file:
    api_version = yaml.safe_load(api_version_file)

"""Beacon Configuration."""


#
# Beacon general info
#
beacon_id = 'itx.french.beacon'  # ID of the Beacon
beacon_name = 'ITX French Beacon'  # Name of the Beacon service
api_version = api_version['api_version'] # Version of the Beacon implementation
uri = 'https://umr1087.univ-nantes.fr/'

#
# Beacon granularity
#
default_beacon_granularity = "record"
max_beacon_granularity = "record"

#
#  Organization info
#
org_id = 'ITX'  # Id of the organization
org_name = 'L`Institut du Thorax'  # Full name
org_description = ('L`institut du thorax is a joint structure in translational research '
                   'dedicated to cardiac, vascular, metabolic and respiratory diseases.  '
                   'Recognized by the INSERM, the CNRS, Nantes Université '
                   'and Nantes University-Hospital, it involves around 800 collaborators. ')
org_adress = ('8 quai Moncousu'
              '44000, Nantes, France')
org_welcome_url = 'https://umr1087.univ-nantes.fr'
org_contact_url = 'u1087@univ-nantes.fr'
org_logo_url = 'https://umr1087.univ-nantes.fr/medias/photo/itx_1564495949407-png'
org_info = ''

#
# Project info
#
#description = (r"This <a href='https://beacon-project.io/'>Beacon</a> "
#               r"is based on the GA4GH Beacon "
#               r"<a href='https://github.com/ga4gh-beacon/specification-v2/blob/master/beacon.yaml'>v2.0</a>")
description = r"This Beacon is based on synthetic data"
version = 'v2.0'
welcome_url = 'https://umr1087.univ-nantes.fr'
alternative_url = 'https://10-54-1-83.gcp.glicid.fr/'
create_datetime = '2023-11-29T12:00:00.000000'
update_datetime = ''
# update_datetime will be created when initializing the beacon, using the ISO 8601 format

#
# Service
#
service_type = 'org.ga4gh:beacon:1.0.0'  # service type
service_url = 'https://10-54-1-83.gcp.glicid.fr/'
entry_point = False
is_open = True
documentation_url = 'https://10-54-1-83.gcp.glicid.fr/doc'  # Documentation of the service
environment = 'test'  # Environment (production, development or testing/staging deployments)

# GA4GH
ga4gh_service_type_group = 'org.ga4gh'
ga4gh_service_type_artifact = 'beacon'
ga4gh_service_type_version = '1.0'

# Beacon handovers
beacon_handovers ={
        'handoverType': {
            'id': 'CUSTOM:000001',
            'label': 'Project description'
        },
        'note': 'Project description',
        'url': 'https://www.nist.gov/programs-projects/genome-bottle'
    }

#
# Database connection
#
database_host = 'mongo'
database_port = 27017
database_user = 'root'
database_password = 'example'
database_name = 'beacon'
database_auth_source = 'admin'
# database_schema = 'public' # comma-separated list of schemas
# database_app_name = 'beacon-appname' # Useful to track connections

#
# Web server configuration
# Note: a Unix Socket path is used when behind a server, not host:port
#
beacon_host = '0.0.0.0'
beacon_port = 5050
beacon_tls_enabled = False
beacon_tls_client = False
beacon_cert = '/etc/ega/server.cert'
beacon_key = '/etc/ega/server.key'
CA_cert = '/etc/ega/CA.cert'

#
# Permissions server configuration
#
permissions_url = 'http://beacon-permissions:5051/'

#
# IdP endpoints (OpenID Connect/Oauth2)
#
# or use Elixir AAI (see https://elixir-europe.org/services/compute/aai)
#

idp_url = 'http://idp:8080/'
#idp_url = 'http://localhost:8080/'

#
# UI
#
autocomplete_limit = 16
autocomplete_ellipsis = '...'

#
# Ontologies
#
ontologies_folder = "ontologies"

alphanumeric_terms = ['libraryStrategy', 'molecularAttributes.geneIds', 'diseases.ageOfOnset.iso8601duration']