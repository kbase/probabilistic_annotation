''' KEGG enzyme '''

class KEGGEnzyme:
    
    ''' A KEGG enzyme record contains the following fields:

        Entry    KEGG ENZYME contains the information about Enzyme Nomenclature obtained from the
                 ExplorEnz database. Additional information is included both computationally and
                 manually. Manually added information includes the KEGG reaction data with parent-child
                 (general to more specific) relationship and the source organism and protein sequence
                 information, whenever available, in each reference. Links to various data in KEGG and
                 other databases are computationally generated. The accession number of an ENZYME entry
                 is the EC (Enzyme Commission) number given by the Nomenclature Committee of the IUBMB
                 and the IUPAC. It is preceded by EC here, simply because of a historical reason.
        Name    The accepted name and alternative names of the enzyme.
        Class   The class, subclass, and sub-subclass of the enzyme, which correspond to the first
                three figures of the EC number.
        Sysname    The systematic name of the enzyme given by the Nomenclature Committee. It represents
                the nature of the chemical reaction.
        Reaction(KEGG)    This field contains all the reactions (R numbers) with this EC number in KEGG.
                Additional reactions are either examples of the general reaction formula by IUBMB,
                designated by ">" for the parent-child relationship, or unrelated to the IUBMB reaction
                formula designated by "(other)".
        Substrate  The chemical compounds that appear on the left side of the reaction equation.
        Product    The chemical compounds that appear on the right side of the reaction equation.
        Comment    Text information commenting on the enzyme.
        Pathway    Links to the KEGG pathway maps, where the corresponding enzyme is marked in red.
        Orthology  Link to the K number entry in the KEGG Orthology (KO) database, which corresponds to
                the ortholog group for the enzyme.
        Genes   Links to the GENES database entries with the assignment (through the KO system) of the
                corresponding EC number.
        Reaction(IUBMB)    The chemical reaction catalyzed by the enzyme is shown in the form of an
                equation or in a text description as defined by the IUBMB/IUPAC Nomenclature Committee.
                The corresponding R number in the KEGG REACTION database is also shown.
        History    Update history given by the IUBMB Enzyme Nomenclature Committee
        Reference  References describing the enzyme.

    '''
    
    def __init__(self):
        ''' Initialize object.
        '''

        self.id = None
        self.name = list()
        self.obsolete = 0
        self.classname = list()
        self.sysname = None
        self.reaction = list()
        self.rxnid = list()
        self.allreactions = list()
        self.substrate = list()
        self.product = list()
        self.comment = list()
        self.pathway = list()
        self.orthology =  list()
        self.genes = dict()
        return
    
    def asDict(self):
        ''' Return the Enzyme object as a dictionary.
        
            @return Dictionary representing Enzyme object
        '''

        enzyme = dict()
        if self.id is not None:
            enzyme['id'] = self.id
        enzyme['obsolete'] = self.obsolete
        if len(self.name) > 0 :
            enzyme['name'] = self.name
        if len(self.classname) > 0:
            enzyme['classname'] = self.classname
        if self.sysname is not None:
            enzyme['sysname'] = self.sysname
        if len(self.reaction) > 0:
            enzyme['reaction'] = self.reaction
        if len(self.allreactions) > 0:
            enzyme['allreactions'] = self.allreactions
        if len(self.substrate) > 0:
            enzyme['substrate'] = self.substrate
        if len(self.product) > 0:
            enzyme['product'] = self.product
        if len(self.comment) > 0:
            enzyme['comment'] = self.comment
        if len(self.pathway) > 0:
            enzyme['pathway'] = self.pathway
        if len(self.orthology) > 0 :
            enzyme['orthology'] = self.orthology
        if len(self.genes) > 0:
            enzyme['genes'] = self.genes
        
        return enzyme
    
    def parse(self, record):
        ''' Parse a record from the flat file database to complete this object.

            @param record: Array of lines with enzyme record from flat file database
            @return Nothing 
        '''
        
        for line in record:
            if line[:3] == '///': # End of record delimiter
                return
            
            # A field in the record has the name in the first 12 characters of the line.
            # If the 12 character prefix is blank, this line is a part of the current field.
            if line[:12] != '            ':
                fieldName = line[:12].strip()
            value = line[12:].strip()
            
            # Process each field in the record.
            if fieldName == 'ENTRY':
                # Entry is in the format EC n.n.n.n
                parts = value.split()
                self.id = parts[1]
                if parts[2] == 'Obsolete':
                    self.obsolete = 1
            elif fieldName == 'NAME':
                # A continuation of the field is delimited by a semicolon at the end of the line.
                self.name.append(value.strip(';'))
            elif fieldName == 'CLASS':
                self.classname.append(value.strip(';'))
            elif fieldName == 'SYSNAME':
                self.sysname = value
            elif fieldName == 'REACTION':
                # Extract the reaction numbers for easy cross referencing.
                val = value.strip(';')
                startPos = val.find('[RN:')
                if startPos >= 0:
                    idList = val[startPos+4:-1]
                    self.rxnid.extend(idList.split())
                self.reaction.append(val)
            elif fieldName == 'ALL_REAC':
                self.allreactions.append(value.strip(';'))
            elif fieldName == 'SUBSTRATE':
                self.substrate.append(value.strip(';'))
            elif fieldName == 'PRODUCT':
                self.product.append(value.strip(';'))
            elif fieldName == 'COMMENT':
                self.comment.append(value)
            elif fieldName == 'PATHWAY':
                self.pathway.append( [ value[:7], value[8:].strip() ] )
            elif fieldName == 'ORTHOLOGY':
                self.orthology.append( [ value[:6], value[7:].strip() ] )
            elif fieldName == 'GENES':
                # Field has organism code, a colon, and a space separated list of genes
                pos = value.find(':')
                organism = value[:pos].lower()
                self.genes[organism] = value[pos+2:].split()
            
        return   

    def makeRecord(self):
        ''' Make an enzyme record for the flat file database.
        
            @return List of lines with enzyme record
        '''
        record = list()
        entry = 'ENTRY       EC %s' %(self.id)
        if self.obsolete:
            pad = 30 - len(entry)
            for index in range(pad):
                entry += ' '
            entry += 'Obsolete  Enzyme'
        else:
            pad = 40 - len(entry)
            for index in range(pad):
                entry += ' '
            entry += 'Enzyme'
        record.append(entry)
        if len(self.name) > 0:
            line = 'NAME        %s' %(self.name[0])
            if len(self.name) == 1:
                record.append(line)
            else:
                record.append(line+';')
                for index in range(1, len(self.name)-1):
                    record.append('            %s;' %(self.name[index]))
                record.append('            %s' %(self.name[-1]))
        if len(self.classname) > 0:
            line = 'CLASS       %s' %(self.classname[0])
            if len(self.classname) == 1:
                record.append(line)
            else:
                record.append(line+';')
                for index in range(1, len(self.classname)-1):
                    record.append('            %s;' %(self.classname[index]))
                record.append('            %s' %(self.classname[-1]))
        if self.sysname is not None:
            record.append('SYSNAME     %s' %(self.sysname))
        if len(self.reaction) > 0:
            line = 'REACTION    %s' %(self.reaction[0])
            if len(self.reaction) == 1:
                record.append(line)
            else:
                record.append(line+';')
                for index in range(1, len(self.reaction)-1):
                    record.append('            %s;' %(self.reaction[index]))
                record.append('            %s' %(self.reaction[-1]))
        if len(self.allreactions) > 0:
            line = 'ALL_REAC    %s' %(self.allreactions[0])
            if len(self.allreactions) == 1:
                record.append(line)
            else:
                record.append(line+';')
                for index in range(1, len(self.allreactions)-1):
                    record.append('            %s;' %(self.allreactions[index]))
                record.append('            %s' %(self.allreactions[-1]))
        if len(self.substrate) > 0:
            line = 'SUBSTRATE   %s' %(self.substrate[0])
            if len(self.substrate) == 1:
                record.append(line)
            else:
                record.append(line+';')
                for index in range(1, len(self.substrate)-1):
                    record.append('            %s;' %(self.substrate[index]))
                record.append('            %s' %(self.substrate[-1]))
        if len(self.product) > 0:
            line = 'PRODUCT     %s' %(self.product[0])
            if len(self.product) == 1:
                record.append(line)
            else:
                record.append(line+';')
                for index in range(1, len(self.product)-1):
                    record.append('            %s;' %(self.product[index]))
                record.append('            %s' %(self.product[-1]))
        if len(self.comment) > 0:
            record.append('COMMENT     %s' %(self.comment[0]))
            for index in range(1, len(self.comment)):
                record.append('            %s' %(self.comment[index]))
        if len(self.pathway) > 0:
            record.append('PATHWAY     %s  %s' %(self.pathway[0][0], self.pathway[0][1]))
            for index in range(1, len(self.pathway)):
                record.append('            %s  %s' %(self.pathway[index][0], self.pathway[index][1]))
        if len(self.orthology) > 0:
            record.append('ORTHOLOGY   %s  %s' %(self.orthology[0][0], self.orthology[0][1]))
            for index in range(1, len(self.orthology)):
                record.append('            %s  %s' %(self.orthology[index][0], self.orthology[index][1]))
        if len(self.genes) > 0:
            first = True
            for key in self.genes:
                if first:
                    prefix = 'GENES       '
                    first = False
                else:
                    prefix = '            '
                record.append('%s%s: %s' %(prefix, key, ' '.join(self.genes[key])))
        record.append('///')
        return record
