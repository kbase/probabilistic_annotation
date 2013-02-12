try:
    import json
except ImportError:
    import sys
    sys.path.append('simplejson-2.3.3')
    import simplejson as json
    
import urllib



class ProbabilisticAnnotation:

    def __init__(self, url):
        if url != None:
            self.url = url

    def annotation_probabilities(self, genomeTO):

        arg_hash = { 'method': 'ProbabilisticAnnotation.annotation_probabilities',
                     'params': [genomeTO],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None

    def annotation_probabilities_id(self, genome_id):

        arg_hash = { 'method': 'ProbabilisticAnnotation.annotation_probabilities_id',
                     'params': [genome_id],
                     'version': '1.1'
                     }

        body = json.dumps(arg_hash)
        resp_str = urllib.urlopen(self.url, body).read()
        resp = json.loads(resp_str)

        if 'result' in resp:
            return resp['result'][0]
        else:
            return None




        
