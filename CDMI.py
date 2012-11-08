try:
    import json
except ImportError:
    import sys
    sys.path.append('simplejson-2.3.3')
    import simplejson as json
    
import urllib



class CDMI_API:

    def __init__(self, url):
        if url != None:
            self.url = url

    def fids_to_annotations(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_annotations',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        print resp_str
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_functions(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_functions',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_literature(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_literature',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_protein_families(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_protein_families',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_roles(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_roles',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_subsystems(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_subsystems',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_co_occurring_fids(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_co_occurring_fids',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_locations(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_locations',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def locations_to_fids(self, region_of_dna_strings):

        arg_hash = { 'method': 'CDMI_API.locations_to_fids',
                     'params': [region_of_dna_strings],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def locations_to_dna_sequences(self, locations):

        arg_hash = { 'method': 'CDMI_API.locations_to_dna_sequences',
                     'params': [locations],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def proteins_to_fids(self, proteins):

        arg_hash = { 'method': 'CDMI_API.proteins_to_fids',
                     'params': [proteins],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def proteins_to_protein_families(self, proteins):

        arg_hash = { 'method': 'CDMI_API.proteins_to_protein_families',
                     'params': [proteins],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def proteins_to_literature(self, proteins):

        arg_hash = { 'method': 'CDMI_API.proteins_to_literature',
                     'params': [proteins],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def proteins_to_functions(self, proteins):

        arg_hash = { 'method': 'CDMI_API.proteins_to_functions',
                     'params': [proteins],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def proteins_to_roles(self, proteins):

        arg_hash = { 'method': 'CDMI_API.proteins_to_roles',
                     'params': [proteins],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def roles_to_proteins(self, roles):

        arg_hash = { 'method': 'CDMI_API.roles_to_proteins',
                     'params': [roles],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def roles_to_subsystems(self, roles):

        arg_hash = { 'method': 'CDMI_API.roles_to_subsystems',
                     'params': [roles],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def roles_to_protein_families(self, roles):

        arg_hash = { 'method': 'CDMI_API.roles_to_protein_families',
                     'params': [roles],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_coexpressed_fids(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_coexpressed_fids',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def protein_families_to_fids(self, protein_families):

        arg_hash = { 'method': 'CDMI_API.protein_families_to_fids',
                     'params': [protein_families],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def protein_families_to_proteins(self, protein_families):

        arg_hash = { 'method': 'CDMI_API.protein_families_to_proteins',
                     'params': [protein_families],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def protein_families_to_functions(self, protein_families):

        arg_hash = { 'method': 'CDMI_API.protein_families_to_functions',
                     'params': [protein_families],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def protein_families_to_co_occurring_families(self, protein_families):

        arg_hash = { 'method': 'CDMI_API.protein_families_to_co_occurring_families',
                     'params': [protein_families],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def co_occurrence_evidence(self, pairs_of_fids):

        arg_hash = { 'method': 'CDMI_API.co_occurrence_evidence',
                     'params': [pairs_of_fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def contigs_to_sequences(self, contigs):

        arg_hash = { 'method': 'CDMI_API.contigs_to_sequences',
                     'params': [contigs],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def contigs_to_lengths(self, contigs):

        arg_hash = { 'method': 'CDMI_API.contigs_to_lengths',
                     'params': [contigs],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def contigs_to_md5s(self, contigs):

        arg_hash = { 'method': 'CDMI_API.contigs_to_md5s',
                     'params': [contigs],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def md5s_to_genomes(self, md5s):

        arg_hash = { 'method': 'CDMI_API.md5s_to_genomes',
                     'params': [md5s],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def genomes_to_md5s(self, genomes):

        arg_hash = { 'method': 'CDMI_API.genomes_to_md5s',
                     'params': [genomes],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def genomes_to_contigs(self, genomes):

        arg_hash = { 'method': 'CDMI_API.genomes_to_contigs',
                     'params': [genomes],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def genomes_to_fids(self, genomes, types_of_fids):

        arg_hash = { 'method': 'CDMI_API.genomes_to_fids',
                     'params': [genomes, types_of_fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def genomes_to_taxonomies(self, genomes):

        arg_hash = { 'method': 'CDMI_API.genomes_to_taxonomies',
                     'params': [genomes],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def genomes_to_subsystems(self, genomes):

        arg_hash = { 'method': 'CDMI_API.genomes_to_subsystems',
                     'params': [genomes],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def subsystems_to_genomes(self, subsystems):

        arg_hash = { 'method': 'CDMI_API.subsystems_to_genomes',
                     'params': [subsystems],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def subsystems_to_fids(self, subsystems, genomes):

        arg_hash = { 'method': 'CDMI_API.subsystems_to_fids',
                     'params': [subsystems, genomes],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def subsystems_to_roles(self, subsystems, aux):

        arg_hash = { 'method': 'CDMI_API.subsystems_to_roles',
                     'params': [subsystems, aux],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def subsystems_to_spreadsheets(self, subsystems, genomes):

        arg_hash = { 'method': 'CDMI_API.subsystems_to_spreadsheets',
                     'params': [subsystems, genomes],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_roles_used_in_models(self, ):

        arg_hash = { 'method': 'CDMI_API.all_roles_used_in_models',
                     'params': [],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def complexes_to_complex_data(self, complexes):

        arg_hash = { 'method': 'CDMI_API.complexes_to_complex_data',
                     'params': [complexes],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def genomes_to_genome_data(self, genomes):

        arg_hash = { 'method': 'CDMI_API.genomes_to_genome_data',
                     'params': [genomes],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_regulon_data(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_regulon_data',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def regulons_to_fids(self, regulons):

        arg_hash = { 'method': 'CDMI_API.regulons_to_fids',
                     'params': [regulons],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_feature_data(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_feature_data',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def equiv_sequence_assertions(self, proteins):

        arg_hash = { 'method': 'CDMI_API.equiv_sequence_assertions',
                     'params': [proteins],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_atomic_regulons(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_atomic_regulons',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def atomic_regulons_to_fids(self, atomic_regulons):

        arg_hash = { 'method': 'CDMI_API.atomic_regulons_to_fids',
                     'params': [atomic_regulons],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_protein_sequences(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_protein_sequences',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_proteins(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_proteins',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_dna_sequences(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_dna_sequences',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def roles_to_fids(self, roles, genomes):

        arg_hash = { 'method': 'CDMI_API.roles_to_fids',
                     'params': [roles, genomes],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def reactions_to_complexes(self, reactions):

        arg_hash = { 'method': 'CDMI_API.reactions_to_complexes',
                     'params': [reactions],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def reaction_strings(self, reactions, name_parameter):

        arg_hash = { 'method': 'CDMI_API.reaction_strings',
                     'params': [reactions, name_parameter],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def roles_to_complexes(self, roles):

        arg_hash = { 'method': 'CDMI_API.roles_to_complexes',
                     'params': [roles],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def complexes_to_roles(self, complexes):

        arg_hash = { 'method': 'CDMI_API.complexes_to_roles',
                     'params': [complexes],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_subsystem_data(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_subsystem_data',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def representative(self, genomes):

        arg_hash = { 'method': 'CDMI_API.representative',
                     'params': [genomes],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def otu_members(self, genomes):

        arg_hash = { 'method': 'CDMI_API.otu_members',
                     'params': [genomes],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def fids_to_genomes(self, fids):

        arg_hash = { 'method': 'CDMI_API.fids_to_genomes',
                     'params': [fids],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def text_search(self, input, start, count, entities):

        arg_hash = { 'method': 'CDMI_API.text_search',
                     'params': [input, start, count, entities],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def corresponds(self, fids, genome):

        arg_hash = { 'method': 'CDMI_API.corresponds',
                     'params': [fids, genome],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def corresponds_from_sequences(self, g1_sequences, g1_locations, g2_sequences, g2_locations):

        arg_hash = { 'method': 'CDMI_API.corresponds_from_sequences',
                     'params': [g1_sequences, g1_locations, g2_sequences, g2_locations],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def close_genomes(self, genomes, how, n):

        arg_hash = { 'method': 'CDMI_API.close_genomes',
                     'params': [genomes, how, n],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def representative_sequences(self, seq_set, rep_seq_parms):

        arg_hash = { 'method': 'CDMI_API.representative_sequences',
                     'params': [seq_set, rep_seq_parms],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result']
        else:
            return None

    def align_sequences(self, seq_set, align_seq_parms):

        arg_hash = { 'method': 'CDMI_API.align_sequences',
                     'params': [seq_set, align_seq_parms],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None




class CDMI_EntityAPI:

    def __init__(self, url):
        if url != None:
            self.url = url

    def get_entity_AlignmentTree(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_AlignmentTree',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_AlignmentTree(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_AlignmentTree',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_AlignmentTree(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_AlignmentTree',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_AlleleFrequency(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_AlleleFrequency',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_AlleleFrequency(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_AlleleFrequency',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_AlleleFrequency(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_AlleleFrequency',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Annotation(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Annotation',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Annotation(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Annotation',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Annotation(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Annotation',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Assay(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Assay',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Assay(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Assay',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Assay(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Assay',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_AtomicRegulon(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_AtomicRegulon',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_AtomicRegulon(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_AtomicRegulon',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_AtomicRegulon(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_AtomicRegulon',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Attribute(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Attribute',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Attribute(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Attribute',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Attribute(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Attribute',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Biomass(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Biomass',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Biomass(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Biomass',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Biomass(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Biomass',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_BiomassCompound(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_BiomassCompound',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_BiomassCompound(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_BiomassCompound',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_BiomassCompound(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_BiomassCompound',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Compartment(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Compartment',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Compartment(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Compartment',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Compartment(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Compartment',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Complex(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Complex',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Complex(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Complex',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Complex(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Complex',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Compound(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Compound',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Compound(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Compound',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Compound(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Compound',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Contig(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Contig',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Contig(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Contig',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Contig(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Contig',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_ContigChunk(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_ContigChunk',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_ContigChunk(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_ContigChunk',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_ContigChunk(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_ContigChunk',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_ContigSequence(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_ContigSequence',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_ContigSequence(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_ContigSequence',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_ContigSequence(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_ContigSequence',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_CoregulatedSet(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_CoregulatedSet',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_CoregulatedSet(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_CoregulatedSet',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_CoregulatedSet(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_CoregulatedSet',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Diagram(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Diagram',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Diagram(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Diagram',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Diagram(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Diagram',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_EcNumber(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_EcNumber',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_EcNumber(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_EcNumber',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_EcNumber(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_EcNumber',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Experiment(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Experiment',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Experiment(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Experiment',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Experiment(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Experiment',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Family(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Family',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Family(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Family',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Family(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Family',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Feature(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Feature',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Feature(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Feature',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Feature(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Feature',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Genome(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Genome',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Genome(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Genome',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Genome(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Genome',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Locality(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Locality',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Locality(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Locality',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Locality(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Locality',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Media(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Media',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Media(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Media',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Media(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Media',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Model(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Model',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Model(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Model',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Model(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Model',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_ModelCompartment(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_ModelCompartment',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_ModelCompartment(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_ModelCompartment',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_ModelCompartment(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_ModelCompartment',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_OTU(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_OTU',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_OTU(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_OTU',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_OTU(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_OTU',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_ObservationalUnit(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_ObservationalUnit',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_ObservationalUnit(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_ObservationalUnit',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_ObservationalUnit(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_ObservationalUnit',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_PairSet(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_PairSet',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_PairSet(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_PairSet',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_PairSet(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_PairSet',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Pairing(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Pairing',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Pairing(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Pairing',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Pairing(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Pairing',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_ProbeSet(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_ProbeSet',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_ProbeSet(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_ProbeSet',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_ProbeSet(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_ProbeSet',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_ProteinSequence(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_ProteinSequence',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_ProteinSequence(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_ProteinSequence',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_ProteinSequence(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_ProteinSequence',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Publication(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Publication',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Publication(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Publication',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Publication(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Publication',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Reaction(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Reaction',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Reaction(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Reaction',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Reaction(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Reaction',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_ReactionRule(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_ReactionRule',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_ReactionRule(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_ReactionRule',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_ReactionRule(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_ReactionRule',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Reagent(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Reagent',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Reagent(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Reagent',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Reagent(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Reagent',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Requirement(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Requirement',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Requirement(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Requirement',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Requirement(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Requirement',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Role(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Role',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Role(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Role',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Role(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Role',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_SSCell(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_SSCell',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_SSCell(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_SSCell',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_SSCell(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_SSCell',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_SSRow(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_SSRow',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_SSRow(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_SSRow',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_SSRow(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_SSRow',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Scenario(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Scenario',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Scenario(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Scenario',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Scenario(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Scenario',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Source(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Source',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Source(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Source',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Source(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Source',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_StudyExperiment(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_StudyExperiment',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_StudyExperiment(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_StudyExperiment',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_StudyExperiment(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_StudyExperiment',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Subsystem(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Subsystem',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Subsystem(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Subsystem',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Subsystem(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Subsystem',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_SubsystemClass(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_SubsystemClass',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_SubsystemClass(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_SubsystemClass',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_SubsystemClass(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_SubsystemClass',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_TaxonomicGrouping(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_TaxonomicGrouping',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_TaxonomicGrouping(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_TaxonomicGrouping',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_TaxonomicGrouping(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_TaxonomicGrouping',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Trait(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Trait',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Trait(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Trait',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Trait(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Trait',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_entity_Variant(self, ids, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_entity_Variant',
                     'params': [ids, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def query_entity_Variant(self, qry, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.query_entity_Variant',
                     'params': [qry, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def all_entities_Variant(self, start, count, fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.all_entities_Variant',
                     'params': [start, count, fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_AffectsLevelOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_AffectsLevelOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsAffectedIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsAffectedIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Aligns(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Aligns',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsAlignedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsAlignedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_AssertsFunctionFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_AssertsFunctionFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasAssertedFunctionFrom(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasAssertedFunctionFrom',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Concerns(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Concerns',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsATopicOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsATopicOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Contains(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Contains',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsContainedIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsContainedIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Controls(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Controls',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsControlledUsing(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsControlledUsing',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Describes(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Describes',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsDescribedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsDescribedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Displays(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Displays',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsDisplayedOn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsDisplayedOn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Encompasses(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Encompasses',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsEncompassedIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsEncompassedIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Formulated(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Formulated',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_WasFormulatedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_WasFormulatedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_GeneratedLevelsFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_GeneratedLevelsFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_WasGeneratedFrom(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_WasGeneratedFrom',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasCompoundAliasFrom(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasCompoundAliasFrom',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_UsesAliasForCompound(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_UsesAliasForCompound',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasIndicatedSignalFrom(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasIndicatedSignalFrom',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IndicatesSignalFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IndicatesSignalFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasMember(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasMember',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsMemberOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsMemberOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasParticipant(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasParticipant',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_ParticipatesIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_ParticipatesIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasPresenceOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasPresenceOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsPresentIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsPresentIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasProteinMember(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasProteinMember',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsProteinMemberOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsProteinMemberOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasReactionAliasFrom(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasReactionAliasFrom',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_UsesAliasForReaction(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_UsesAliasForReaction',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasRepresentativeOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasRepresentativeOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsRepresentedIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsRepresentedIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasResultsIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasResultsIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasResultsFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasResultsFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasSection(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasSection',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsSectionOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsSectionOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasStep(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasStep',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsStepOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsStepOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasTrait(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasTrait',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Measures(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Measures',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasUnits(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasUnits',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsLocated(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsLocated',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasUsage(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasUsage',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsUsageOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsUsageOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasValueFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasValueFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasValueIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasValueIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasVariationIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasVariationIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsVariedIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsVariedIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Impacts(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Impacts',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsImpactedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsImpactedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Includes(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Includes',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsIncludedIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsIncludedIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IncludesPart(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IncludesPart',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsPartOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsPartOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IndicatedLevelsFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IndicatedLevelsFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasLevelsFrom(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasLevelsFrom',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Involves(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Involves',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsInvolvedIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsInvolvedIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsARequirementIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsARequirementIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsARequirementOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsARequirementOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsAnnotatedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsAnnotatedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Annotates(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Annotates',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsAssayedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsAssayedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsLocatedOn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsLocatedOn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsClassFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsClassFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsInClass(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsInClass',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsCollectionOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsCollectionOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsCollectedInto(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsCollectedInto',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsComposedOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsComposedOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsComponentOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsComponentOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsComprisedOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsComprisedOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Comprises(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Comprises',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsConfiguredBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsConfiguredBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_ReflectsStateOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_ReflectsStateOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsConsistentWith(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsConsistentWith',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsConsistentTo(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsConsistentTo',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsCoregulatedWith(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsCoregulatedWith',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasCoregulationWith(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasCoregulationWith',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsCoupledTo(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsCoupledTo',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsCoupledWith(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsCoupledWith',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsDefaultFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsDefaultFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_RunsByDefaultIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_RunsByDefaultIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsDefaultLocationOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsDefaultLocationOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasDefaultLocation(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasDefaultLocation',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsDeterminedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsDeterminedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Determines(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Determines',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsDividedInto(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsDividedInto',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsDivisionOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsDivisionOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsExemplarOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsExemplarOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasAsExemplar(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasAsExemplar',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsFamilyFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsFamilyFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_DeterminesFunctionOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_DeterminesFunctionOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsFormedOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsFormedOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsFormedInto(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsFormedInto',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsFunctionalIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsFunctionalIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasFunctional(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasFunctional',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsGroupFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsGroupFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsInGroup(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsInGroup',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsImplementedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsImplementedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Implements(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Implements',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsInPair(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsInPair',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsPairOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsPairOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsInstantiatedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsInstantiatedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsInstanceOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsInstanceOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsLocatedIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsLocatedIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsLocusFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsLocusFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsModeledBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsModeledBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Models(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Models',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsOwnerOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsOwnerOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsOwnedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsOwnedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsProposedLocationOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsProposedLocationOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasProposedLocationIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasProposedLocationIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsProteinFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsProteinFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Produces(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Produces',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsRealLocationOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsRealLocationOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasRealLocationIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasRealLocationIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsReferencedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsReferencedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_UsesReference(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_UsesReference',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsRegulatedIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsRegulatedIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsRegulatedSetOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsRegulatedSetOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsRelevantFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsRelevantFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsRelevantTo(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsRelevantTo',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsRepresentedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsRepresentedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_DefinedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_DefinedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsRequiredBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsRequiredBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Requires(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Requires',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsRoleOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsRoleOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasRole(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasRole',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsRowOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsRowOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsRoleFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsRoleFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsSequenceOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsSequenceOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasAsSequence(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasAsSequence',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsSubInstanceOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsSubInstanceOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Validates(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Validates',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsSuperclassOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsSuperclassOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsSubclassOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsSubclassOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsTargetOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsTargetOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Targets(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Targets',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsTaxonomyOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsTaxonomyOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsInTaxa(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsInTaxa',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsTerminusFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsTerminusFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HasAsTerminus(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HasAsTerminus',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsTriggeredBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsTriggeredBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Triggers(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Triggers',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsUsedAs(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsUsedAs',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsUseOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsUseOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Manages(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Manages',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsManagedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsManagedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_OperatesIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_OperatesIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsUtilizedIn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsUtilizedIn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Overlaps(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Overlaps',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IncludesPartOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IncludesPartOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_ParticipatesAs(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_ParticipatesAs',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsParticipationOf(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsParticipationOf',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_ProducedResultsFor(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_ProducedResultsFor',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_HadResultsProducedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_HadResultsProducedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_ProjectsOnto(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_ProjectsOnto',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsProjectedOnto(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsProjectedOnto',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Provided(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Provided',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_WasProvidedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_WasProvidedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Shows(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Shows',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsShownOn(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsShownOn',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Submitted(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Submitted',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_WasSubmittedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_WasSubmittedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Summarizes(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Summarizes',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_SummarizedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_SummarizedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_Uses(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_Uses',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def get_relationship_IsUsedBy(self, ids, from_fields, rel_fields, to_fields):

        arg_hash = { 'method': 'CDMI_EntityAPI.get_relationship_IsUsedBy',
                     'params': [ids, from_fields, rel_fields, to_fields],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None




        
