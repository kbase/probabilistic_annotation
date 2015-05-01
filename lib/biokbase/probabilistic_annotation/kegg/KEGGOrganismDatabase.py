
from biokbase.probabilistic_annotation.kegg.KEGGDatabase import KEGGDatabase
from biokbase.probabilistic_annotation.kegg.KEGGOrganism import KEGGOrganism
from biokbase.probabilistic_annotation.kegg.QueryKEGG import QueryKEGG

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
    
    def download(self):
        ''' Download the organism database from KEGG web service.
        
            @return Nothing
        '''

        # Get the current organism database from KEGG.        
        query = QueryKEGG()
        organism = query.list('organism')

        # What if I also got the list of genomes using all_entities_Genome with id, name, and domain.
        # Then look for a match on name and link to kbase id.  And then figure out if OTU representative.
        # Then limit the amino acid download to OUT representatives.  Just like with probanno today
        cdmi_entity = CDMI_EntityAPI(config["cdmi_url"])
        genomes = cdmi_entity.all_entities_Genome(0, 50000, ['id', 'scientific_name'])
        nameDict = dict()
        for key in genomes:
            name = genomes[key]['scientific_name'].replace('substr. ', '')
            name = name.replace('str. ', '')
            nameDict[name] = key
        print len(nameDict)

        otuGenomes = getOtuGenomeDictionary(1000, config)
        
        # For every prokaryote in organism database, see if there is a match on name.
        numNoMatch = 0
        numProkaryotes = 0
        for id in organismDB.organisms:
            organism = organismDB.get(id)
            if organism.isProkaryote():
                numProkaryotes += 1
                matches = list()
                if organism.name in nameDict:
                    matches.append(nameDict[organism.name])
                else:
                    name = organism.name.split(' (')[0]
                    if name in nameDict:
                        matches.append(nameDict[name])
                    else:
                        parts = name.split()
                        m = parts[0]
                        if len(parts) > 1:
                            m += ' '+parts[1]
                        for key in nameDict:
                            if key.startswith(m):
                                matches.append(nameDict[key])
    #                    print '%s: %s' %(organism.name, matches)
                if len(matches) == 0:
                    print organism.name
                    numNoMatch += 1
                
                for index in range(len(matches)):
                    if matches[index] in otuGenomes:
                        organism.otuRepresentative = 1
                        organismDB.update(organism)  
                        # go through the list of matches, if any kbase id is a representative for the OTU, mark it
                    # go through the keys. find all that start with the organism.name
                    # if only one, pick it
                    # if more than one, search in name for each remaining word
    #                print 'No match for %s %s' %(organism.name, organism.taxonomy)
        print '%d %d' %(numProkaryotes, numNoMatch)
        organismDB.store()
        # Save the organism to the flat file database.  Each line describes an
        # organism and includes the ID, abbreviation, name, and taxonomy.   
        with open(self.filename, 'w') as handle:
            for line in organism:
                handle.write(line+'\t0\n')
            
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
        # and write the records to the file.    
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