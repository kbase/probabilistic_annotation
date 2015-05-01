
import requests
from biokbase.probabilistic_annotation.kegg.KEGGReaction import KEGGReaction
from biokbase.probabilistic_annotation.kegg.KEGGEnzyme import KEGGEnzyme

class QueryError(Exception):
    pass

''' Query Kyoto Encyclopedia of Genes and Genomes using KEGG API '''

class QueryKEGG:
    
    def __init__(self, url=None):
        ''' Initialize object.

            @param url: Base URL of KEGG web service endpoint
        '''

        if url is None:
            self.keggUrl = 'http://rest.kegg.jp'
        else:
            self.keggUrl = url
        self.result = list()

    def _request(self, url):
        ''' Send a request to the KEGG web service.
        
            @param url: KEGG API request URL
            @return Nothing
        '''

        # Send the request and receive the response from the web service.
        try:
            rget = requests.get(url)
        except Exception as e:
            raise QueryError(u'Unable to connect to KEGG server %s: %s' %(url, e))
        
#         if not (rget.ok and rget.text):
#             raise Exception(u'Unable to connect to KEGG server %s: %s' %(url, rget.raise_for_status()))

        # When the status code indicates a failure, raise an exception with the details.
        if rget.status_code != 200:
            if rget.status_code == 400:
                reason = 'bad request (syntax error, wrong database name, etc.)'
            elif rget.status_code == 404:
                reason = 'not found'
            else:
                reason = 'unknown (%d)' %(rget.status_code)
            raise QueryError('Query %s failed: %s' %(url, reason))
        
        # Add the response to the current results.
        self.result.extend(rget.text.split('\n')[:-1]) # Remove the empty last line
        return

    def get(self, dbname, idList, option=None, handle=None):
        ''' Get records from a KEGG database.
        
            @param dbname: Name of database to get records from
            @param idList: List of IDs of records to get
            @param option: Option value for request
            @param filename: Path to file to save retrieved records to
            @return Nothing
        '''

        # Delete the previous result and start fresh with an empty list.        
        del self.result
        self.result = list()
        
        # The web service limits the size of the query to 10 items.
        increment = 10
        start = 0
        if len(idList) > increment:
            end = start + increment
        else:
            end = len(idList)
        counter = len(idList)
        
        # Run through the list of IDs until all have been processed.
        while counter > 0:
            # Build the query string.
            query = ''
            for index in range(start, end):
                query += dbname+':'+idList[index]+'+'

            # Adjust for the next query.
            start += increment
            end += increment
            if end > len(idList):
                end = len(idList)
            counter -= increment
            
            # Send the request.
            url = self.keggUrl+'/get/'+query[:-1] # Take off trailing +
            if option is not None:
                url = url+'/'+option
            try:
                self._request(url)
            except QueryError as e:
                print e
        
        # Save the results to a file if specified.
        if handle is not None:
            for index in range(len(self.result)):
                handle.write(self.result[index]+'\n')

        return

    def list(self, database):
        ''' Return a list of entry identifiers and associated definitions for a given database.
        
            @param database: Database name or set of database entries
            @return List of something
        '''
        
        # Need to limit to 100 database entries.
        
        # Delete the previous result and start fresh with an empty list.        
        del self.result
        self.result = list()
        
        # The web service limits the size of the query to 100 items.
        increment = 100
        start = 0
        if len(idList) > increment:
            end = start + increment
        else:
            end = len(idList)
        counter = len(idList)
        
        url = self.keggUrl+'/list/'+database
        try:
            self._request(url)
        except QueryError as e:
            print e
        return self.result

    def getReactions(self, idList, handle=None):
        ''' Get reactions from the reaction database.
        
            @param idList: List of reaction IDs
            @param filename: Path to file to save returned records to
            @return List of KEGGReaction objects for returned reactions
        '''

        # Get the records from the reaction database.
        self.get('rn', idList, handle=handle)

        # Convert the returned reaction records to KEGGReaction objects.  
        reactionList = list()
        start = 0
        for index in range(len(self.result)):
            if self.result[index] == '///':
                reaction = KEGGReaction()
                reaction.parse(self.result[start:index+1])
                reactionList.append(reaction)
                start = index + 1
    
        return reactionList

    def getEnzymes(self, idList, handle=None):
        ''' Get enzymes from the enzyme database.

            @param idList: List of enzyme IDs
            @param filename: Path to file to save returned records to
            @return List of KEGGEnzyme objects for returned enzymes
        '''
        
        # Get the records from the enzyme database.
        self.get('ec', idList, handle=handle)
        
        # Convert the returned enzyme records to KEGGEnzyme objects.
        enzymeList = list()
        start = 0
        for index in range(len(self.result)):
            if self.result[index] == '///':
                enzyme = KEGGEnzyme()
                enzyme.parse(self.result[start:index+1])
                enzymeList.append(enzyme)
                start = index + 1
    
        return enzymeList
    
    def getAminoAcidSeq(self, code, geneList, handle):
        ''' Get amino acid sequences for a list of genes from an organism.
        
            @param code Organism code (3 or 4 character string)
            @param geneList: List of gene IDs
            @param handle: Handle to file to save returned records to
            @return Nothing
        '''
        
        self.get(code, geneList, option='aaseq', handle=handle)
        return
        # Delete the previous result and start fresh with an empty list.        
        del self.result
        self.result = list()
        
        # Can replace all of this with call to get()
        # After building list of strings.
        
        # The service endpoint limits the size of the query to 10 items.
        increment = 10
        start = 0
        if len(geneList) > increment:
            end = start + increment
        else:
            end = len(geneList)
        counter = len(geneList)
        
        # Run through the list of reaction IDs until all have been processed.
        while counter > 0:
            # Build the query string.
            query = ''
            for index in range(start, end):
                query += geneList[index][0]+':'+geneList[index][1]+'+'
            start += increment
            end += increment
            if end > len(geneList):
                end = len(geneList)
            counter -= increment
            
            # Send the request.
            url = self.keggUrl+'/get/'+query[:-1]+'/aaseq' # Take off trailing +
            try:
                self._request(url)
            except QueryError as e:
                print e
        
        for index in range(len(self.result)):
            handle.write(self.result[index]+'\n')
        return #aminoAcidList
