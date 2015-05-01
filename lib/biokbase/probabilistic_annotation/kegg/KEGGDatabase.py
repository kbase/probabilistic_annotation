
''' Exception when a record is not found. '''
class RecordNotFound(Exception):
    pass

''' Base class for managing a KEGG flat file database '''

class KEGGDatabase:
    
    def __init__(self, filename):
        return
    
    def getRecord(self, handle):
        ''' Get a record from a database file.
        
            @param handle: File handle of database file
            @return List of lines in record
        '''
        
        record = list()
        for line in handle:
            record.append(line.strip('\n'))
            if line[:3] == '///':
                yield record
                record = list()
                continue
        