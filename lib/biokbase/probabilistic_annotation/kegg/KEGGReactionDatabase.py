
from biokbase.probabilistic_annotation.kegg.KEGGDatabase import KEGGDatabase
from biokbase.probabilistic_annotation.kegg.KEGGReaction import KEGGReaction

''' Manage a KEGG reaction flat file database. '''

class KEGGReactionDatabase(KEGGDatabase):
    
    def __init__(self, filename):
        ''' Initialize object.
        
            @param filename: Path to reaction database file
        '''

        self.filename = filename
        self.reactions = dict()
        return
    
    def load(self):
        ''' Load the reaction database from a flat file.
        
            @return Nothing
        '''
        
        # Read all of the records form the flat file database and build a dictionary
        # of KEGGReaction objects.
        with open(self.filename, 'r') as handle:
            for record in self.getRecord(handle):
                rxn = KEGGReaction()
                rxn.parse(record)
                self.reactions[rxn.id] = rxn
        
        return
    
    def store(self, filename=None):
        ''' Save the reaction database to a flat file.
        
            @param filename: Path to reaction database file
            @return Nothing
        '''
        
        if filename is None:
            filename = self.filename
            
        # Convert all of the KEGGReaction objects to flat file database records
        # and write the records to the file.    
        with open(filename, 'w') as handle:
            for key in self.reactions:
                record = self.reactions[key].makeRecord()
                for line in record:
                    handle.write(line+'\n')
        return

    def size(self):
        ''' Get the number of records in the reaction database.
        
            @return Number of records
        '''

        return len(self.reactions)
    
    def get(self, id):
        ''' Get a reaction with the specified ID.
        
            @param id: Reaction ID of reaction to get from database
            @return KEGGReaction object for specified ID or None when not found
        '''

        if id in self.reactions:
            return self.reactions[id]
        return None
#        raise RecordNotFound('Reaction ID %s was not found' %(id))

    def exists(self, id):
        ''' Check if an ID exists in the reaction database.
        
            @param id: ID to check
            @return True when ID exists, otherwise False
        '''

        if id in self.reactions:
            return True
        return False
    
    def update(self, reaction):
        ''' Update a reaction in the database (add new or replace existing reaction).
        
            @param reaction: KEGGReaction object for reaction to add or replace
            @return Nothing
        '''

        # Delete the current KEGGReaction object if it exists in the dictionary.
        if reaction.id in self.reactions:
            del self.reactions[reaction.id]
            
        # Add the new KEGGReaction object to the dictionary.
        self.reactions[reaction.id] = reaction
        return