
from biokbase.probabilistic_annotation.kegg.KEGGDatabase import KEGGDatabase
from biokbase.probabilistic_annotation.kegg.KEGGOrganism import KEGGOrganism
from biokbase.probabilistic_annotation.kegg.QueryKEGG import QueryKEGG
from biokbase.probabilistic_annotation.CDMExtractor import CDMExtractor
from biokbase.cdmi.client import CDMI_EntityAPI

''' Manage a KEGG organism flat file database. '''

class KEGGOrganismDatabase(KEGGDatabase):
    
    def __init__(self, filename):
        ''' Initialize object.
        
            @param filename: Path to organism database file
        '''

        self.filename = filename
        self.organisms = dict() # Keyed by T number
        self.codeToId = dict() # Keyed by organism code
        return
    
    def download(self, config):
        ''' Download the organism database from KEGG web service.

            @param config: Dictionary of configuration variables
            @return Nothing
        '''

        # Get the current organism database from KEGG.        
        query = QueryKEGG()
        organismList = query.list('organism')

        # Get the list of all genomes available in the central data model.
        cdmi_entity = CDMI_EntityAPI(config['cdmi_url'])
        genomes = cdmi_entity.all_entities_Genome(0, 50000, ['id', 'scientific_name'])
        nameDict = dict()
        for key in genomes:
            # Remove 'substr' and 'str' from scientific name for better matching to KEGG genome names.
            name = genomes[key]['scientific_name'].replace('substr. ', '')
            name = name.replace('str. ', '')
            nameDict[name] = key

        # Get the list of representative OTU genome IDs from the central data model.
        extractor = CDMExtractor(config)
        otuGenomes = extractor.getOtuGenomeDictionary(1000)
        
        # For every prokaryote in organism database, see if there is a match on name.
        numNoMatch = 0
        numMatch = 0
        numProkaryotes = 0
        numRepresentatives = 0
        for index in range(len(organismList)):
            # Add the organism to the database.
            record = organismList[index]
            record +='\t0' # Mark as not a representative by default
            organism = KEGGOrganism()
            organism.parse(record)
            self.organisms[organism.id] = organism
            
            # Only need to check prokaryote organisms.
            if organism.isProkaryote():
                numProkaryotes += 1
                matches = list()
                
                # Check for exact match on name to CDM genome.
                if organism.name in nameDict:
                    matches.append(nameDict[organism.name])
                else:
                    # Strip words after a parenthesis and check again for exact match on
                    # name to CDM genome.
                    name = organism.name.split(' (')[0]
                    if name in nameDict:
                        matches.append(nameDict[name])
                    else:
                        # Extract the first two words in the name and check again for CDM
                        # genomes that start with the two words.
                        parts = name.split()
                        m = parts[0]
                        if len(parts) > 1:
                            m += ' '+parts[1]
                        for key in nameDict:
                            if key.startswith(m):
                                matches.append(nameDict[key])

                # If there was a match to just one CDM genome, see if the genome is the
                # OTU representative.  If so, mark the KEGG organism as the representative.
                if len(matches) == 1:
                    numMatch += 1
                    if matches[0] in otuGenomes:
                        organism.otuRepresentative = 1
                        numRepresentatives += 1
                else:
                    # When there are no matches or multiple matches, we cannot tell if the
                    # KEGG organism is the same as the CDM genome.  We are conservative to
                    # avoid having lots of E. coli genomes marked as representatives.
                    numNoMatch += 1
                        
        # Save the organism to the flat file database.     
        self.store()
        
        print 'prokaryotes: %d, matches: %d, no match: %d, representatives: %d' %(numProkaryotes, numMatch, numNoMatch, numRepresentatives)
        return
            
    def load(self):
        ''' Load the organism database.
        
            @return Nothing
        '''
        
        # Read all of the records form the flat file database and build a dictionary
        # of KEGGOrganism objects.
        with open(self.filename, 'r') as handle:
            for record in handle:
                organism = KEGGOrganism()
                organism.parse(record)
                self.organisms[organism.id] = organism
                self.codeToId[organism.code] = organism.id
        return
    
    def store(self, filename=None):
        ''' Save the organism database to a flat file.
        
            @param filename: Path to organism database file
            @return Nothing
        '''
        
        if filename is None:
            filename = self.filename
        
        # Convert all of the KEGGOrganism objects to flat file database records
        # and write the records to the file.  Each line describes an organism and
        # includes the ID, abbreviation, name, taxonomy, and OTU representative flag.
        with open(filename, 'w') as handle:
            for key in sorted(self.organisms):
                record = self.organisms[key].makeRecord()
                handle.write(record+'\n')
        return

    def size(self):
        ''' Get the number of records in the organism database.
        
            @return Number of records
        '''

        return len(self.organisms)
    
    def numOtuRepresentatives(self):
        ''' Get the number of OTU representative organisms in the database.

            @return Number of OTU representative organisms
        '''

        count = 0
        for key in self.organisms:
            if self.organisms[key].isRepresentative():
                count += 1
        return count

    def numProkaryotes(self):
        ''' Get the number of prokaryote organisms in the database.

            @return Number of prokaryote organisms
        '''

        count = 0
        for key in self.organisms:
            if self.organisms[key].isProkaryote():
                count += 1
        return count

    def numEukaryotes(self):
        ''' Get the number of eukaryote organisms in the database.

            @return Number of eukaryote organisms
        '''

        count = 0
        for key in self.organisms:
            if self.organisms[key].isEukaryote():
                count += 1
        return count

    def exists(self, id):
        ''' Check if an ID exists in the organism database.
        
            @param id: ID to check
            @return True when ID exists, otherwise False
        '''

        if id in self.organisms:
            return True
        return False
    
    def get(self, id):
        ''' Get an organism with the specified ID.
        
            @param id: Organism ID to get from database
            @return KEGGOrganism object for specified ID or None when not found
        '''

        if id in self.organisms:
            return self.organisms[id]
        return None
#        raise RecordNotFound('Reaction ID %s was not found' %(id))

    def getByCode(self, code):
        ''' Get an organism with the specified code.
        
            @param code: Code of organism to get from database
            @return KEGGOrganism object for specified code or None when not found
        '''

        if code in self.codeToId:
            return self.organisms[self.codeToId[code]]
        return None

    def update(self, organism):
        ''' Update an organism in the database (add new or replace existing organism).
        
            @param organism: KEGGOrganism object for organism to add or replace
            @return Nothing
        '''

        # Delete the current KEGGOrganism object if it exists in the dictionary.
        if organism.id in self.organisms:
            del self.organisms[organism.id]
            
        # Add the new KEGGOrganism object to the dictionary.
        self.organisms[organism.id] = organism
        return