
''' KEGG reaction '''

class KEGGReaction:
    
    ''' A KEGG reaction record contains the following fields:

        Entry   The KEGG REACTION database is a manually curated collection of biochemical reactions,
                mostly enzymatic reactions but including some spontaneous reactions, which are derived
                from KEGG ENZYME (Enzyme Nomenclature) or KEGG PATHWAY (KEGG metabolic pathway maps).
                Each entry is identified by the unique identifier called the R number ('R' followed
                by five-digit number).
        Name    The name of the reaction, usually the systematic name of the enzyme that catalyzes the
                reaction (can be not specified).
        Definition    The chemical reaction in the form of an equation. The reaction is assumed to be
                reversible and reactants (substrates and products) are separated by '<=>'. Each compound
                in the left or the right side is separated by ' + '. There may be a coefficient before
                the compound name.
        Equation    The C number representation of the reaction equation with links to the COMPOUND
                database entries.
        Remark  The same reaction entries, if any, are given.
        Comment Text information commenting on the reaction.
        RPair   Links to the corresponding KEGG RPAIR database entries, which contain chemical structure
                transformation patterns of reactant pairs (substrate-product pairs).
        Enzyme  Links to the corresponding KEGG ENZYME database entries.
        Pathway Links to the KEGG pathway maps, where the corresponding reaction is marked in red.
        Orthology    Links to the corresponding KEGG ORTHOLOGY (KO) database entries.
        Module  Links to the corresponding KEGG MODULE database entries.
        Reference    References for the reaction.
    '''
    
    def __init__(self):
        ''' Initialize object.
        '''

        self.id = None
        self.name = list()
        self.definition = None
        self.equation = None
        self.remark = None
        self.comment = list()
        self.rpair = list()
        self.enzyme = list()
        self.pathway = list()
        self.orthology =  list()
        self.reference = list()
        self.module = list()
        return
    
    def asDict(self):
        ''' Return the Reaction object as a dictionary.
        
            @return Dictionary representing Reaction object
        '''

        rxn = dict()
        if self.id is not None:
            rxn['id'] = self.id
        if len(self.name) > 0 :
            rxn['name'] = self.name
        if self.definition is not None:
            rxn['definition'] = self.definition
        if self.equation is not None:
            rxn['equation'] = self.equation
        if self.remark is not None:
            rxn['remark'] = self.remark
        if len(self.comment) > 0:
            rxn['comment'] = self.comment
        if len(self.rpair) > 0:
            rxn['rpair'] = self.rpair
        if len(self.enzyme) > 0:
            rxn['enzyme'] = self.enzyme
        if len(self.pathway) > 0:
            rxn['pathway'] = self.pathway
        if len(self.orthology) > 0 :
            rxn['orthology'] = self.orthology
        if len(self.reference) > 0:
            rxn['reference'] = self.reference
        if len(self.module) > 0:
            rxn['module'] = self.module
        
        return rxn
    
    def parse(self, record):
        ''' Parse a record from the flat file database to complete this object.

            @param record: List of lines with reaction record from flat file database
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
                parts = value.split()
                self.id = parts[0]
            elif fieldName == 'NAME':
                # A continuation of the field is delimited by a semicolon at the end of the line.
                self.name.append(value.strip(';'))
            elif fieldName == 'DEFINITION':
                self.definition = value
            elif fieldName == 'EQUATION':
                self.equation = value
            elif fieldName == 'REMARK':
                self.remark = value
            elif fieldName == 'COMMENT':
                self.comment.append(value)
            elif fieldName == 'RPAIR':
                parts = value.split()
                self.rpair.append(tuple(parts))
            elif fieldName == 'ENZYME':
                # Multiple enzymes are all listed on the same line.
                self.enzyme.extend(value.split())
            elif fieldName == 'PATHWAY':
                self.pathway.append( [ value[:7], value[8:].strip() ] )
            elif fieldName == 'ORTHOLOGY':
                self.orthology.append( [ value[:6], value[7:].strip() ] )
            elif fieldName == 'REFERENCE':
                self.reference.append(value)
            elif fieldName == 'MODULE':
                self.module.append( [ value[:7], value[9:].strip() ] )
            else:
                print 'Skipping field '+fieldName+' with value '+value
        return

    def makeRecord(self):
        ''' Make a reaction record for the flat file database.
        
            @return List of lines with reaction record
        '''

        record = list()
        record.append('ENTRY       %s                      Reaction' %(self.id))
        if len(self.name) > 0:
            line = 'NAME        %s' %(self.name[0])
            if len(self.name) == 1:
                record.append(line)
            else:
                record.append(line+';')
                for index in range(1, len(self.name)-1):
                    record.append('            %s;' %(self.name[index]))
                record.append('            %s' %(self.name[-1]))
        if self.definition is not None:
            record.append('DEFINITION  %s' %(self.definition))
        if self.equation is not None:
            record.append('EQUATION    %s' %(self.equation))
        if self.remark is not None:
            record.append('REMARK      %s' %(self.remark))
        if len(self.comment) > 0:
            record.append('COMMENT     %s' %(self.comment[0]))
            for index in range(1, len(self.comment)):
                record.append('            %s' %(self.comment[index]))
        if len(self.rpair) > 0:
            line = 'RPAIR       %s  ' %(self.rpair[0][0])
            line += ' '.join(self.rpair[0][1:])
            record.append(line)
            for index in range(1, len(self.rpair)):
                line = '            %s  ' %(self.rpair[index][0])
                line += ' '.join(self.rpair[index][1:])
                record.append(line)
        if len(self.enzyme) > 0:
            enzymes = '       '.join(self.enzyme)
            record.append('ENZYME      '+enzymes)
        if len(self.pathway) > 0:
            record.append('PATHWAY     %s  %s' %(self.pathway[0][0], self.pathway[0][1]))
            for index in range(1, len(self.pathway)):
                record.append('            %s  %s' %(self.pathway[index][0], self.pathway[index][1]))
        if len(self.orthology) > 0:
            record.append('ORTHOLOGY   %s  %s' %(self.orthology[0][0], self.orthology[0][1]))
            for index in range(1, len(self.orthology)):
                record.append('            %s  %s' %(self.orthology[index][0], self.orthology[index][1]))
        if len(self.reference) > 0:
            record.append('REFERENCE   %s' %(self.reference[0]))
            for index in range(1, len(self.reference)):
                record.append('            %s' %(self.reference[index]))
        if len(self.module) > 0:
            record.append('MODULE      %s  %s' %(self.module[0][0], self.module[0][1]))
            for index in range(1, len(self.module)):
                record.append('            %s  %s' %(self.module[index][0], self.module[index][1]))
        record.append('///')
        return record
            
            
                