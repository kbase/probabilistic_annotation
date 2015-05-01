
''' KEGG organism '''

class KEGGOrganism:
    
    ''' A KEGG organism record contains the following fields:

        ID      Each entry is identified by the unique identifier called the T number ('T' followed
                by five-digit number).
        Code    The code for the organism is a three or four character string based on the name
                and is used for identifying genes in complete genomes.
        Name    The scientific name of the organism (common name in parenthesis).
        Taxonomy    The taxonomy of the organism as a string with levels separated by semicolons.
        OtuRep  Flag indicating if the organism is representative of the OTU.

        Each record is a single line with fields separated by tabs.
    '''
    
    def __init__(self):
        self.id = None
        self.code = None
        self.name = None
        self.taxonomy = tuple()
        self.otuRepresentative = 0
        return
    
    def asDict(self):
        ''' Return the Organism object as a dictionary.
        
            @return Dictionary representing Organism object
        '''

        organism = dict()
        if self.id is not None:
            organism['id'] = self.id
        if self.code is not None:
            organism['code'] = self.code
        if self.name is not None:
            organism['name'] = self.name
        organism['taxonomy'] = self.taxonomy
        organism['otuRepresentative'] = self.otuRepresentative
        return organism
    
    def parse(self, record):
        ''' Parse a record from the flat file database.

            @param record: Line with record
            @return Nothing
        '''

        # Each record is a single line with the fields separated by tabs. 
        fields = record.split('\t')
        self.id = fields[0]
        self.code = fields[1]
        self.name = fields[2]
        self.taxonomy = tuple(fields[3].split(';')) # Taxonomy levels are separated by semicolons
        self.otuRepresentative = int(fields[4])
        return

    def makeRecord(self):
        ''' Make a record for the flat file database.
        
            @return Line with organism record
        '''
        
        taxonomy = ';'.join(self.taxonomy)
        record = '%s\t%s\t%s\t%s\t%d' %(self.id, self.code, self.name, taxonomy, self.otuRepresentative)
        return record

    def isProkaryote(self):
        ''' Check if the organism is a prokaryote.

            @return True if the organism is a prokaryote.
        '''

        if self.taxonomy[0] == 'Prokaryotes':
            return True
        return False

    def isEukaryote(self):
        ''' Check if the organism is a eukaryote.

            @return True if the organism is a eukaryote.
        '''
        
        if self.taxonomy[0] == 'Eukaryotes':
            return True
        return False

    def isBacteria(self):
        ''' Check if the organism is a bacteria.

            @return True if the organism is a bacteria.
        '''
        
        if self.taxonomy[1] == 'Bacteria':
            return True
        return False
    
    def isArchaea(self):
        ''' Check if the organism is an archaea.

            @return True if the organism is an archaea.
        '''
        
        if self.taxonomy[1] == 'Archaea':
            return True
        return False
    
    def isRepresentative(self):
        ''' Check if the organism is a representative of the OTU.

            @return True if organism is an OTU representative
        '''

        if self.otuRepresentative:
            return True
        return False
                