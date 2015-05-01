
import re
import time
from biokbase.fbaModelServices.Client import fbaModelServices, ServerError as FbaModelServerError
from biokbase.probabilistic_annotation.kegg.KEGGReactionDatabase import KEGGReactionDatabase
from biokbase.probabilistic_annotation.kegg.KEGGEnzymeDatabase import KEGGEnzymeDatabase
from biokbase.probabilistic_annotation.kegg.KEGGOrganismDatabase import KEGGOrganismDatabase
from biokbase.probabilistic_annotation.kegg.QueryKEGG import QueryKEGG
from biokbase.probabilistic_annotation.Helpers import now

''' Extract data to build static database files from KEGG databases '''

class KEGGExtractor:
    
    def __init__(self, config, reactionFile, enzymeFile, organismFile):
        ''' Initialize object.

            @param config: Dictionary of configuration variables
            @param reactionFile: Path to reaction database file
            @param enzymeFile: Path to enzyme database file
            @param organismFile: Path to organism database file
        '''

        # Save the configuration variables.
        self.config = config
        
        # Load the reaction database.
        self.reactionDB = KEGGReactionDatabase(reactionFile)
        self.reactionDB.load()
        
        # Load the enzyme database.    
        self.enzymeDB = KEGGEnzymeDatabase(enzymeFile)
        self.enzymeDB.load()
        
        # Assign a complex ID to every enzyme in the enzyme database.
        self.nextComplexId = 10000
        self.enzymeComplexDict = dict()
        for key in sorted(self.enzymeDB.enzymes):
            self.enzymeComplexDict[key] = 'kb|cpx.%d' %(self.nextComplexId)
            self.nextComplexId += 1

        # Create the organism database.
        self.organismDB = KEGGOrganismDatabase(organismFile)

        # Create a query object.
        self.query = QueryKEGG(url=self.config['kegg_url'])

        return

    def getKeggReactionList(self, kbaseRxnIdList, increment, biochem, biochemws):
        ''' Get the KEGG reaction IDs for the corresponding KBase reaction IDs.

            The returned list contains tuples with three elements: (1) KEGG reaction
            ID, (2) KBase reaction ID, (3) KBase reaction name.

            @param kbaseRxnIdList: List of KBase reaction IDs (any of the aliases)
            @param increment: Size for breaking list into smaller parts
            @param biochem: ID of Biochemistry object to use for lookups
            @param biochemws: ID of workspace containing Biochemistry object
            @return List of tuples as described above
        '''

        # Create client object.
        fbaClient = fbaModelServices(self.config['fbamodelservices_url'])
        
        # Get the reaction details from the fba modeling service.  For every reaction that
        # has an alias in KEGG format, add it to a list of tuples with the KEGG reaction
        # ID, KBase reaction ID, and reaction name.
        # Maybe keep a set to prune duplicates
        keggReactionList = list()
        
        # Setup get_reactions() input parameters that are the same for every call.
        getReactionsParams = dict()
        getReactionsParams['biochemistry'] = biochem
        getReactionsParams['biochemistry_workspace'] = biochemws
        
        # Break up the list into smaller parts to avoid timeouts.  
        start = 0
        if len(kbaseRxnIdList) < increment:
            end = start + len(kbaseRxnIdList)
        else:
            end = start + increment
        counter = len(kbaseRxnIdList)
        while counter > 0:
            # Get the reaction details for this part of the list.
            getReactionsParams['reactions'] = kbaseRxnIdList[start:end]
            reactionList = fbaClient.get_reactions(getReactionsParams)
    
            # Look through the aliases in each reaction.  If an alias matches the KEGG
            # reaction ID format, add it to the set.  A KEGG reaction ID is in the
            # format Rddddd.  Note that a KBase reaction can have an alias to more than
            # one KEGG reaction and sometimes there are duplicates.
            for rindex in range(len(reactionList)):
                if reactionList[rindex] is None:
                    print 'Reaction %s was not found in biochemistry' %(kbaseRxnIdList[start+rindex])
                    continue
                if 'aliases' in reactionList[rindex]:
                    rxnAliasSet = set() # Use a set to remove duplicates
                    for aindex in range(len(reactionList[rindex]['aliases'])):
                        alias = reactionList[rindex]['aliases'][aindex]
                        if re.match(r'^R\d+$', alias):
                            rxnAliasSet.add(alias)
                    for alias in rxnAliasSet:
                        keggReactionList.append( [ alias, reactionList[rindex]['id'], reactionList[rindex]['name'] ] )

            # Update for the next loop iteration.
            print 'Processed %d reactions at %s' %(end, now())
            counter -= len(reactionList)
            start += increment
            end += increment
            if end > len(kbaseRxnIdList):
                end = len(kbaseRxnIdList)
    
        # Sort the list by KEGG reaction ID.
        keggReactionList.sort()
        return keggReactionList
    
    def updateLocalDatabases(self, keggReactionList):
        ''' Update the local reaction and enzyme databases from the current databases.

            The input list contains tuples with three elements: (1) KEGG reaction
            ID, (2) KBase reaction ID, (3) KBase reaction name.  If a reaction is
            not in the local database, query the web service and update the local
            database with more current data.  For each reaction, make sure the
            linked enzymes are found in the enzyme database.  If an enzyme is not
            in the local database, query the web service and update the local
            database with more current data.

            @param keggReactionList: List of tuples with KEGG reaction IDs
            @return Nothing
        '''

        # Build a list of the KEGG reaction IDs that are not found in the local database.
        # The first element in the keggReactionList tuple is the KEGG reaction ID.
        queryList = list()
        for index in range(len(keggReactionList)):
            if self.reactionDB.exists(keggReactionList[index][0]) == False:
                queryList.append(keggReactionList[index][0])
        print '%d reactions not found in local KEGG reaction database' %(len(queryList))
    
        # Query the current KEGG database through the web service for the missing reactions.
        if len(queryList) > 0:
            rxnList = self.query.getReactions(queryList)
            print 'Found %d reactions in current KEGG database' %(len(rxnList))
            
            if len(rxnList) > 0:
                # Add the reactions returned from the query to the local database.
                for index in range(len(rxnList)):
                    self.reactionDB.update(rxnList[index])
            
                # Save the updated database.
                self.reactionDB.store()
                print 'There are %d reactions in updated local reaction database' %(self.reactionDB.size())

        # Go through the list of KEGG reactions and find the links to enzymes.
        numEnzymeFound = 0
        numEnzymeNotFound = 0
        numObsolete = 0
        reactionNoEnzyme = list() # Reactions that do not have a link to an enzyme
        enzymeNotFound = list() # Enzymes that are not found in the database
        for rindex in range(len(keggReactionList)):
            # Get the reaction record.
            reaction = self.reactionDB.get(keggReactionList[rindex][0])
            if reaction is not None:
                if len(reaction.enzyme) > 0:
                    # There can be more than one enzyme so check each one.
                    found = False
                    for eindex in range(len(reaction.enzyme)):
                        enzyme = self.enzymeDB.get(reaction.enzyme[eindex])
                        if enzyme is not None:
                            found = True
                            if enzyme.obsolete:
                                print 'Enyzme %s in reaction %s is obsolete' %(reaction.enzyme[eindex], reaction.id)
                                numObsolete += 1
                        else:
                            enzymeNotFound.append(reaction.enzyme[eindex]) 
                    if found:
                        numEnzymeFound += 1
                    else:
                        numEnzymeNotFound += 1
                else:
                    reactionNoEnzyme.append(reaction.id)
        print 'Found %d reactions with links to enzymes in local enzyme database' %(numEnzymeFound)
        print 'Found %d enzymes that are obsolete' %(numObsolete)
        print 'Found %d reactions where enzymes were not found in local enzyme database' %(numEnzymeNotFound)
        
        # For the reactions that have no enyzme, query the current database to see if the
        # reaction has been updated.
        print 'Found %d reactions with no link to an enzyme' %(len(reactionNoEnzyme))
        rxnList = self.query.getReactions(reactionNoEnzyme)
        
        # Update the reactions returned from the query in the local database.
        if len(rxnList) > 0:
            numNoEnzyme = 0
            for rindex in range(len(rxnList)):                
                # If the updated reaction has at least one enzyme, update the local reaction database.
                reaction = rxnList[rindex]
                if len(reaction.enzyme) > 0:
                    self.reactionDB.update(reaction)
                    
                    # Look for the new enzymes in the local enzyme database.
                    for eindex in range(len(reaction.enzyme)):
                        if not self.enzymeDB.exists(reaction.enzyme[eindex]):
                            enzymeNotFound.append(reaction.enzyme[eindex])
                else:
                    numNoEnzyme += 1

            # Save the updated database if reactions were updated.
            print 'Found %d reactions with updated enzymes in current reaction database' %(len(rxnList)-numNoEnzyme)
            print 'Found %d reactions with no enzyme in current reaction database' %(numNoEnzyme)
            if len(rxnList)-numNoEnzyme > 0:
                self.reactionDB.store()
                print 'There are %d reactions in updated local reaction database' %(self.reactionDB.size())

        # Query the current database for the missing enzymes.
        if len(enzymeNotFound) > 0:
            print '%d enzymes not found in local enzyme database' %(len(enzymeNotFound))
            enzymeList = self.query.getEnzymes(enzymeNotFound)
            print 'Found %d enzymes in the current database' %(len(enzymeList))
            
            if len(enzymeList) > 0:
                # Update the enzymes returned from the query to the local database.
                for index in range(len(enzymeList)):
                    self.enzymeDB.update(enzymeList[index])
                    self.enzymeComplexDict[enzymeList[index].id] = self.nextComplexId
                    self.nextComplexId += 1
                
                # Save the updated database.
                self.enzymeDB.store()
                print 'There are %d enzymes in updated local enzyme database' %(self.enzymeDB.size())

        return
    
    def mapReactionToComplex(self, keggReactionList):
        ''' Make a mapping of reactions to a list of complexes for the reaction.

            The input list contains tuples with three elements: (1) KEGG reaction
            ID, (2) KBase reaction ID, (3) KBase reaction name.  For each reaction,
            the list of complexes is generated from the list of enzymes.

            @param keggReactionList: List of tuples with KEGG reaction IDs
            @return Dictionary keyed by reaction ID with list of complexes for reaction
        '''

        # Build a dictionary keyed by reaction ID with a list of complexes.
        # We call the enzyme the complex and use our assigned complex IDs.
        numNoEnzyme = 0
        reactionToComplexes = dict()
        print 'Building mapping with %d reactions' %(len(keggReactionList))
        for rindex in range(len(keggReactionList)):
            # Get the reaction record.
            reaction = self.reactionDB.get(keggReactionList[rindex][0])
            if reaction is not None:
                if len(reaction.enzyme) > 0:
                    complexList = list()
                    for eindex in range(len(reaction.enzyme)):
                        try:
                            complexList.append(self.enzymeComplexDict[reaction.enzyme[eindex]])
                        except KeyError: # Just ignore the enzyme if it is not found
                            pass
                    if len(complexList) > 0:
                        rxnid = re.sub(r'rxn0*(\d+)', r'kb|rxn.\1', keggReactionList[rindex][1]) # Convert ModelSEED format to KBase format
                        if rxnid not in reactionToComplexes:
                            reactionToComplexes[rxnid] = complexList
                        else:    
                            reactionToComplexes[rxnid].extend(complexList)
                    else:
                        numNoEnzyme += 1
                else:
                    numNoEnzyme += 1
        
        print '%d reactions have no valid enzyme links' %(numNoEnzyme)
        return reactionToComplexes
    
    def mapComplexToRole(self):
        ''' Make a mapping of complexes to a list of roles for the complex.
        
            @return Dictionary keyed by complex ID with list of roles for complex
        '''

        # Build a dictionary keyed by complex ID with a list of roles.
        # We call the enzyme the complex and the roles are the names of the enzyme.
        complexToRoles = dict()
        for key in self.enzymeComplexDict:
            enzyme = self.enzymeDB.get(key)
            complexToRoles[self.enzymeComplexDict[key]] = list()
            if len(enzyme.name) > 0:
                for nindex in range(len(enzyme.name)):
                    complexToRoles[self.enzymeComplexDict[key]].append('%s (%s)' %(enzyme.name[nindex], key))
            else:
                complexToRoles[self.enzymeComplexDict[key]].append('Unknown (%s)' %(key))

        return complexToRoles
    
    def mapFeatureToRole(self, fastaPath):
        ''' Make a mapping of features to roles and download amino acid sequences
            for the features from organisms that are the representative for the OTU.

            @note If fastaPath is None, sequence download is skipped.
            @param fastaPath: Path to fasta file for storing amino acid sequences
            @return fidsToRoles dictionary
        '''

        # Load the organism database.
        self.organismDB.load()
    
        # Build a dictionary keyed by feature ID with a list of roles and get
        # the amino acid sequences for the feature IDs.
        fidsToRoles = dict()
        
        # Open the fasta file for storing the amino acid sequences.
        if fastaPath is not None:
            handle = open(fastaPath, 'w')
        
        # For every enzyme, look at each organism with genes.  If the organism
        # is a prokaryote and is the representative for the OTU, add the genes
        # to the dictionary and download the amino acid sequences for the genes. 
        numSkipped = 0
        for key in self.enzymeComplexDict:
            enzyme = self.enzymeDB.get(key)
            
            # The genes in the enzyme are stored in a dictionary keyed by organism
            # code.  Build a second dictionary keyed by organism code for just
            # the prokaryote organisms that tracks if the organism is a
            # representative for its OTU.
            prokaryotes = dict()
            numRepresentatives = 0
            for code in enzyme.genes:
                organism = self.organismDB.getByCode(code)
                if organism.isProkaryote(): 
                    prokaryotes[code] = organism.isRepresentative()
                    if prokaryotes[code]:
                        numRepresentatives += 1
            
            # If there are no genes from prokaryote organisms, skip this enzyme.
            if len(prokaryotes) == 0:
                numSkipped += 1
                continue

            # There are no prokaryote organisms that are a representative for
            # their OTU.  For this enzyme, get the amino acid sequences for
            # every gene.
            if numRepresentatives == 0:
                for key in prokaryotes:
                    prokaryotes[key] = True
            
            # For the prokaryote organisms that are OTU representatives, get the
            # amino acid sequence of the gene and add it to the fids to roles mapping.
            for code in prokaryotes:
                if prokaryotes[code]:
                    geneList = list()
                    for gene in enzyme.genes[code]:
                        # Remove the common name of the gene if it is specified.
                        pos = gene.find('(')
                        if pos >= 0:
                            gene = gene[:pos]
                        geneList.append(gene)
                        featureId = '%s:%s' %(code, gene)
                        fidsToRoles[featureId] = list()
                        if len(enzyme.name) > 0:
                            for nindex in range(len(enzyme.name)):
                                fidsToRoles[featureId].append('%s (%s)' %(enzyme.name[nindex], enzyme.id))
                        else:
                            fidsToRoles[featureId].append('Unknown (%s)' %(enzyme.id))
                    if fastaPath is not None:
                        self.query.getAminoAcidSeq(code, geneList, handle)
                        time.sleep(1) # Wait a second between queries to be nice to KEGG

        # Close the amino acid sequence fasta file.
        if fastaPath is not None:
            handle.close()
        
        print '%d enzymes skipped because there were no genes from prokaryotes' %(numSkipped)
        return fidsToRoles
    
    def updateOrganism(self):
        organismDB = KEGGOrganismDatabase('organism.db')
    #    organismDB.download()
        organismDB.load()
        print '%d %d' %(organismDB.numOtuRepresentatives(), organismDB.numProkaryotes())
        otuGenomes = getOtuGenomeDictionary(1000, config)
        print len(otuGenomes)
        exit(0)
        # Get the list of all genomes with scientific names.
        # Move this to organismDB.download()
        cdmi_entity = CDMI_EntityAPI(config["cdmi_url"])
        genomes = cdmi_entity.all_entities_Genome(0, 50000, ['id', 'scientific_name'])
        nameDict = dict()
        for key in genomes:
            name = genomes[key]['scientific_name'].replace('substr. ', '')
            name = name.replace('str. ', '')
            nameDict[name] = key
        print len(nameDict)
    #     for key in nameDict:
    #         print '%s: %s' %(key, nameDict[key])
       
    #     dataParser = DataParser(config) 
    #     otus, prokotus = dataParser.readOtuData()
        
        
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

