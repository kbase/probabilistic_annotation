
from biokbase.probabilistic_annotation.kegg.KEGGDatabase import KEGGDatabase
from biokbase.probabilistic_annotation.kegg.KEGGEnzyme import KEGGEnzyme

''' Manage a KEGG enzyme flat file database. '''

class KEGGEnzymeDatabase(KEGGDatabase):
    
    def __init__(self, filename):
        ''' Initialize object.
        
            @param filename: Path to enzyme database file
        '''

        self.filename = filename
        self.enzymes = dict()
        return
    
    def load(self):
        ''' Load the enzyme database.
        
            @return Nothing
        '''
        
        # Read all of the records form the flat file database and build a dictionary
        # of KEGGEnzyme objects.
        with open(self.filename, 'r') as handle:
            for record in self.getRecord(handle):
                enzyme = KEGGEnzyme()
                enzyme.parse(record)
                self.enzymes[enzyme.id] = enzyme
        
        return
    
    def store(self, filename=None):
        ''' Save the enzyme database to a flat file.
        
            @param filename: Path to enzyme database file
            @return Nothing
        '''
        
        if filename is None:
            filename = self.filename
        
        # Convert all of the KEGGEnzyme objects to flat file database records
        # and write the records to the file.    
        with open(filename, 'w') as handle:
            for key in sorted(self.enzymes):
                record = self.enzymes[key].makeRecord()
                for line in record:
                    handle.write(line+'\n')
        return

    def size(self):
        ''' Get the number of records in the enzyme database.
        
            @return Number of records
        '''

        return len(self.enzymes)
    
    def exists(self, id):
        ''' Check if an ID exists in the enzyme database.
        
            @param id: ID to check
            @return True when ID exists, otherwise False
        '''

        if id in self.enzymes:
            return True
        return False
    
    def get(self, id):
        ''' Get an enzyme with the specified ID.
        
            @param id: Enzyme ID to get from database
            @return KEGGEnzyme object for specified ID or None when not found
        '''

        if id in self.enzymes:
            return self.enzymes[id]
        return None
#        raise RecordNotFound('Reaction ID %s was not found' %(id))

    def update(self, enzyme):
        ''' Update an enzyme in the database (add new or replace existing enzyme).
        
            @param enzyme: KEGGEnzyme object for enzyme to add or replace
            @return Nothing
        '''

        # Delete the current KEGGEnzyme object if it exists in the dictionary.
        if enzyme.id in self.enzymes:
            del self.enzymes[enzyme.id]
            
        # Add the new KEGGEnzyme object to the dictionary.
        self.enzymes[enzyme.id] = enzyme
        return