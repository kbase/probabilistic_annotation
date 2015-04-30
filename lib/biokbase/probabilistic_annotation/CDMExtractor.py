#!/usr/bin/python

from biokbase.cdmi.client import CDMI_API, CDMI_EntityAPI
from biokbase.probabilistic_annotation.Helpers import now
from urllib2 import URLError, HTTPError
import urllib
import sys
import operator #for itemgetter
try:
    import json
except ImportError:
    sys.path.append('simplejson-2.3.3')
    import simplejson as json

''' Extract data to build static database files from Central Data Model '''

class CDMExtractor():

    def __init__(self, config):
        ''' Initialize object.
        
            @param config: Dictionary of configuration variables
        '''

        # Save the configuration variables.
        self.config = config

        # Create client objects.
        self.cdmi = CDMI_API(self.config['cdmi_url'])
        self.cdmiEntity = CDMI_EntityAPI(self.config['cdmi_url'])
        return

    def getFieldFromEntity(self, seedEntity, fieldName):
        ''' Get a field from a entity returned by a get_entity_XXX() or all_entities_XXX() function.
    
            @param seedEntity: Dictionary keyed by ID to a dictionary of key-value pairs
            @param fieldName: Name of key in dictionary of key-value pairs
            @return List of values for the specified field (key) in all of the entities
            @raise ValueError: Value of input argument is invalid
        '''
    
        if seedEntity is None:
            msg = 'Value of seedEntity is None which usually means you searched for something that does not exist in the CDM'
            raise ValueError(msg)

        # Check for an error I seem to make all the time and yell at me in a USEFUL way
        if not isinstance(seedEntity, dict):
            msg = 'Value of seedEntity is not a dictionary.  Perhaps you meant to call getFieldFromRelationship?'
            raise ValueError(msg)

        f = list()
        for entry in seedEntity:
            if fieldName not in seedEntity[entry]:
                msg = 'Field name %s not found in seedEntity' %(fieldName)
                raise ValueError(msg)
            f.append(seedEntity[entry][fieldName])
        return f
    
    def getFieldFromRelationship(self, seedRelationship, fieldName, objtype):
        ''' Get a field from an object returned by a get_relationship_XXX() function.
    
            The get_relationship_XXX() functions return lists of lists.
            The first list is a list of all the links
            The second list has three dictionaries in it: the TO dictionary, the REL dictionary
            and the FROM dictionary describing properties on either end of the link and of the link itself.
    
            If you want to  maintain duplicate relationships (many-to-one, one-to-many, many-to-many),
            this function should be called at least twice (once on each end of the relationship, or once
            on an end and once in the middle).
    
            @param seedRelationship: Output from a get_relationship_XXX() function
            @param fieldName: Field to extract from the object
            @param objtype: Type of object, "TO", "REL", or "FROM"
            @return List (in the same order as the list from the get_relationship function)
                of the values with the specified field name
            @raise ValueError: Value of input argument is invalid
        '''
    
        if seedRelationship is None:
            msg = 'Value of seedRelationship is None which usually means you were searching for something that does not exist in the CDM'
            raise ValueError(msg)
        objIndex = None
        if objtype.lower() == 'from':
            objIndex = 0
        elif objtype.lower() == 'rel':
            objIndex = 1
        elif objtype.lower() == 'to':
            objIndex = 2
        else:
            msg = 'Value of objtype is not TO, REL, or FROM'
            raise ValueError(msg)
        if not isinstance(seedRelationship, list):
            msg = 'Value of seedRelationship is not a list.  Perhaps you meant to call getFieldFromEntity?'
            raise ValueError(msg)

        # Unravel
        f = list()
        for entry in seedRelationship:
            # TO CHECK: Is it possible to have one of the links lead to nothing?
            # Check field name validity - if it links to something there has to be the data request there
            # or else something is wrong.
            if fieldName not in entry[objIndex]:
                msg = 'Field name %s not found in seedRelationship' %(fieldName)
                raise ValueError(msg)
            f.append(entry[objIndex][fieldName])
        return f
    
    def getSubsystemFids(self, count):
        ''' Query the CDMI for a list of feature IDs in the subsystems.
    
            @param count: Number of entities to retrieve in each function call
            @return List of subsystem feature IDs
        '''
    
        # Get the genes that are in subsystems and in OTUs.
        ssdict = dict()
        start = 0
        done = False
        while not done:
            subdict = self.cdmiEntity.all_entities_Subsystem(start, count, ["id"])
            ssdict.update(subdict)
            start += count
            if len(subdict) < count:
                done = True
        ssids = self.getFieldFromEntity(ssdict, "id")
        print 'Found %d subsystems' %(len(ssids))
    
        # Now lets get a list of FIDs within those subsystems
        # Break the complete list into smaller sub-lists to avoid timeouts
        start = 0
        increment = 10
        end = start + increment
        counter = len(ssids)
        ssfids = []
        while counter > 0:
            try:
                ssfiddict = self.cdmi.subsystems_to_fids(ssids[start:end], [])
            except HTTPError as e:
                if increment > 1:
                    increment = increment / 2
                    end = start + increment
                sys.stderr.write("caught '%s' error, increment is now %d\n" %(e.reason, increment))
                continue
            for key in ssfiddict:
                for ssfid in ssfiddict[key]:
                    ls = ssfiddict[key][ssfid]
                    for arr in ls:
                        if len(arr) > 1:
                            gl = arr[1]
                            for l in gl:
                                ssfids.append(l)
                                
            # Move to next sub-list
            start += increment
            end += increment
            if end >= len(ssids):
                end = len(ssids)
            counter -= increment
    
        # Uniquify!
        return list(set(ssfids))
    
    def getDlitFids(self, count):
        ''' Query the CDMI for a list of feature IDs with direct literature evidence (dlits).
    
            @param count: Number of entities to retrieve in each function call
            @return List of literature feature IDs
        '''
    
        pubdict = dict()
        start = 0
        done = False
        while not done:
            subdict = self.cdmiEntity.all_entities_Publication(start, count, ["id"])
            pubdict.update(subdict)
            start += count
            if len(subdict) < count:
                done = True
    
        pubids = self.getFieldFromEntity(pubdict, "id")
        print 'Found %d publication IDs' %(len(pubids))
        pub2seq = self.cdmiEntity.get_relationship_Concerns(pubids, [], [], ["id"])
        pubseqs = self.getFieldFromRelationship(pub2seq, "id", "to")
        print 'Found %d protein sequences from publications' %(len(pubseqs))
        seq2fids = self.cdmiEntity.get_relationship_IsProteinFor(pubseqs, [], [], ["id"])
        fids = self.getFieldFromRelationship(seq2fids, "id", "to")
        return fids
    
    def filterFidsByOtus(fidlist, otus, config):
        '''
        Obsolete (I think this isn't used any more)
    
        Given a list of representative organism IDs (OTUs) and a list of
        FIDs, returns only those FIDs found in an OTU.'''
    
        # Identify the organism belonging to each fid
        # If this fails to find an organism we don't want it anyway...
        orgdict = self.cdmiEntity.get_relationship_IsOwnedBy(fidlist, [], [], ["id"])
        flist = self.getFieldFromRelationship(orgdict, "from_link", "rel")
        olist = self.getFieldFromRelationship(orgdict, "id", "to")
    
        fids = []
        for ii in range(len(olist)):
            if olist[ii] in otus:
                fids.append(flist[ii])
        return fids
    
    def filterFidsByOtusBetter(fidsToRoles, rolesToFids, oturepsToMembers, config):
        '''Attempt to do a more intelligent filtering of FIDs by OTU.
    
        Given all FIDs attached to a role in the unfiltered set we do the following:
        
        Initialize KEEP
        For each OTU and each role:
           If role is found in the representative, add to KEEP and continue;
           Otherwise, iterate over other genomes.
               If role is found in one other genome, add to KEEP and continue;
    
        This process should make our calculation less sensitive to the choice of OTUs...
    
        '''
    
        # Identify the organism belonging to each fid
        # If this fails to find an organism we don't want it anyway...
        fidlist = fidsToRoles.keys()
        orgdict = []
         # Break the complete list into smaller sub-lists to avoid timeouts
        start = 0
        increment = 5000
        end = start + increment
        counter = len(fidlist)
        while counter > 0:
            try:
                od = self.cdmiEntity.get_relationship_IsOwnedBy(fidlist[start:end], [], [], ["id"])
            except HTTPError as e:
                if increment > 1:
                    increment = increment / 2
                    end = start + increment
                sys.stderr.write("caught '%s' error, increment is now %d\n" %(e.reason, increment))
                continue
            orgdict.extend(od)
            start += increment
            end += increment
            if end >= len(fidlist):
                end = len(fidlist)
            counter -= increment
        fidlist = self.getFieldFromRelationship(orgdict, "from_link", "rel")
        orglist = self.getFieldFromRelationship(orgdict, "id", "to")
        fidToOrg = {}
        for ii in range(len(fidlist)):
            fidToOrg[fidlist[ii]] = orglist[ii]
        
        keptFidsToRoles = {}
        keptRolesToFids = {}
        # If the OTUs are comprehensive this should be empty.
        missingRoles = []
    
        # For each OTU
        for oturep in oturepsToMembers:
            # for each role
            for role in rolesToFids:
                fidlist = rolesToFids[role]
                keepFid = None
                keepRole = None
                for fid in fidlist:
                    # This can happen due to MOL issues
                    if fid not in fidToOrg:
                        continue
                    org = fidToOrg[fid]
                    # If the organism is the representative we keep it and go to the next role
                    if org == oturep:
                        keepFid = fid
                        keepRole = role
                        break
                    # Otherwise look at the rest of the list (note that I just pick one without really paying
                    # attention to WHICH one...). We save them in case there are no examples of the role in the
                    # representative organism, but continue on anyway.
                    if org in oturepsToMembers[oturep]:
                        keepFid = fid
                        keepRole = role
                if keepFid is not None:
                    if keepFid in keptFidsToRoles:
                        keptFidsToRoles[keepFid].append(keepRole)
                    else:
                        keptFidsToRoles[keepFid] = [ keepRole ]
                    if keepRole in keptRolesToFids:
                        keptRolesToFids[keepRole].append(keepFid)
                    else:
                        keptRolesToFids[keepRole] = [ keepFid ]
    
        missingRoles = list(set(rolesToFids.keys()) - set(keptRolesToFids.keys()))
    
    #    print oturepsToMembers
    #    print missingRoles
    #    print keptRolesToFids
    
        return keptFidsToRoles, keptRolesToFids, missingRoles
    
    def filterFidsByOtusOptimized(self, featureIdList, rolesToFids, otuRepsToMembers):
        ''' Filter feature IDs by OTU (optimized version).
    
            To minimize the amount of redundancy in the list of target proteins, filter
            the feature IDs so there is at most one protein from each OTU for each
            functional role.
    
            @param featureIdList: List of unfiltered feature IDs
            @param rolesToFids: Dictionary keyed by role of list of feature IDs
            @param otuRepsToMembers: Dictionary keyed by OTU representative to list of OTU members
            @return Dictionary keyed by feature ID of list of roles, dictionary keyed by role
                of list of feature IDs
        '''
    
        # Identify the organism belonging to each feature ID.
        # If this fails to find an organism we don't want it anyway...
        fidToOrganism = dict() # Map feature IDs to organisms
    
         # Break the complete list into smaller sub-lists to avoid timeouts
        start = 0
        increment = 100000
        end = start + increment
        counter = len(featureIdList)
        while counter > 0:
            try:
                ownedBy = self.cdmiEntity.get_relationship_IsOwnedBy(featureIdList[start:end], [], ['from_link'], ['id'])
            except HTTPError as e:
                if increment > 1:
                    increment = increment / 2
                    end = start + increment
                sys.stderr.write("caught '%s' error, increment is now %d\n" %(e.reason, increment))
                continue
            # just build the dictionary here, run the list of ob, extracting fid from from_link and organism from id
            fidList = self.getFieldFromRelationship(ownedBy, "from_link", "rel")
            organismList = self.getFieldFromRelationship(ownedBy, "id", "to")
            for index in range(len(fidList)):
                fidToOrganism[fidList[index]] = organismList[index]
    
            start += increment
            end += increment
            if end >= len(featureIdList):
                end = len(featureIdList)
            counter -= increment
    
        # Add all possible keys to the dictionaries and initialize the value.
        # Then we don't have to check if the key exists in the main loop below.
        keptFidsToRoles = dict()
        for index in range(len(featureIdList)):
            keptFidsToRoles[featureIdList[index]] = list()
        keptRolesToFids = dict()
        for role in rolesToFids:
            keptRolesToFids[role] = list()
    
        # Find the feature ID (protein) from each OTU for each functional role.
        otuCounter = 0
        for otuRepresentative in otuRepsToMembers:
            # This loop takes a very long time so print a message every so often
            # to track progress.
            otuCounter += 1
            if otuCounter % 10 == 0:
                print 'Processed %d OTUs at %s' %(otuCounter, now())
    
            # Check every functional role.
            for role in rolesToFids:
                keepFid = None
                keepRole = None
                for fid in rolesToFids[role]:
                    # This can happen due to MOL issues
                    if fid not in fidToOrganism:
                        continue
                    organism = fidToOrganism[fid]
    
                    # If the organism is the representative we keep it and go to the next role
                    if organism == otuRepresentative:
                        keepFid = fid
                        keepRole = role
                        break
    
                    # Otherwise look at the rest of the list (note that I just pick one without really paying
                    # attention to WHICH one...). We save them in case there are no examples of the role in the
                    # representative organism, but continue on anyway.
                    if organism in otuRepsToMembers[otuRepresentative]:
                        keepFid = fid
                        keepRole = role
    
                # Add to the dictionaries if we are keeping the feature ID.
                if keepFid is not None:
                    keptFidsToRoles[keepFid].append(keepRole)
                    keptRolesToFids[keepRole].append(keepFid)
    
        # Look for any empty lists and remove them.
        keysToRemove = list()
        for fid in keptFidsToRoles:
            if len(keptFidsToRoles[fid]) == 0:
                keysToRemove.append(fid)
        for key in keysToRemove:
            del keptFidsToRoles[key]
        keysToRemove = list()
        for role in keptRolesToFids:
            if len(keptRolesToFids[role]) == 0:
                keysToRemove.append(role)
        for key in keysToRemove:
            del keptRolesToFids[key]
    
        return keptFidsToRoles, keptRolesToFids
    
    def getOtuGenomeDictionary(self, count):
        ''' Obtain a dictionary from OTU representatives to all genomes in the OTU.
    
            @param count: Number of entities to retrieve in each function call
            @return Dictionary keyed by OTU representative of list of OTU members
        '''

        # Get list of OTUs
        otulist = self.getOtuGenomeIds(count)
        otudict = self.cdmi.otu_members(otulist[0])
        return otudict
    
    def mapFidsToRoles(self, fidlist):
        ''' Given a list of feature IDs return a dictionary keyed by feature ID to
            the list of roles the encoding gene performs and a dictionary keyed by
            roles to the feature IDs performing them.
    
            @param fidlist: List of feature IDs
            @return Dictionary keyed by feature ID of list of roles encoding gene performs, dictionary
                keyed by role of list of feature IDs performing the role
        '''
    
        # Break the complete list into smaller sub-lists to avoid timeouts
        start = 0
        increment = 1000
        end = start + increment
        counter = len(fidlist)
        fidsToRoles = {}
        rolesToFids = {}
        while counter > 0:
            try:
                roledict = self.cdmiEntity.get_relationship_HasFunctional(fidlist[start:end], [], [], ["id"])
            except HTTPError as e:
                if increment > 1:
                    increment = increment / 2
                    end = start + increment
                sys.stderr.write("caught '%s' error, increment is now %d\n" %(e.reason, increment))
                continue
            flist = self.getFieldFromRelationship(roledict, "from_link", "rel")
            rolelist = self.getFieldFromRelationship(roledict, "id", "to")
            for ii in range(len(flist)):
                # We have to use sets here because a bug(?) in get_relationship_HasFunctional allows multiple identical
                # links between fids and roles.
                # See for example what happens when you call it on g.9647.peg.2332
                if flist[ii] in fidsToRoles:
                    fidsToRoles[flist[ii]].add(rolelist[ii])
                else:
                    fidsToRoles[flist[ii]] = set([rolelist[ii]])
                if rolelist[ii] in rolesToFids:
                    rolesToFids[rolelist[ii]].add(flist[ii])
                else:
                    rolesToFids[rolelist[ii]] = set([flist[ii]])
                    
            # Move to next sub-list
            start += increment
            end += increment
            if end >= len(fidlist):
                end = len(fidlist)
            counter -= increment
            
        # Convert back to lists to not break other functions.
        for f in fidsToRoles:
            fidsToRoles[f] = list(fidsToRoles[f])
        for r in rolesToFids:
            rolesToFids[r] = list(rolesToFids[r])
        return fidsToRoles, rolesToFids
    
    def getAminoAcidSequences(self, fidlist):
        ''' Given a list of feature IDs, returns a dictionary keyed by feature ID to
            its amino acid sequence.
    
            @note Features with no amino acid sequence are discarded.
            @param fidlist: List of feature IDs
            @return Dictionary keyed by feature ID of amino acid sequence for feature
        '''
    
        fidlist = list(set(fidlist))
        start = 0
        increment = 5000
        end = start + increment
        counter = len(fidlist)
        seqs = {}
        while counter > 0:
            try:
                ps = self.cdmi.fids_to_protein_sequences(fidlist[start:end])
            except HTTPError as e:
                if increment > 1:
                    increment = increment / 2
                    end = start + increment
                sys.stderr.write("caught '%s' error, increment is now %d\n" %(e.reason, increment))
                continue
            seqs.update(ps)
            
            # Move to next sub-list
            start += increment
            end += increment
            if end >= len(fidlist):
                end = len(fidlist)
            counter -= increment
        
        return seqs
    
    def getPegs(self, genomes):
        ''' Given a list of genome IDs, returns a list of feature IDs for
            protein-encoding genes in the specified genomes.
    
            @param genomes: List of genome IDs
            @return List of feature IDs for protein-encoding genes in specified genomes
        '''
    
        fiddict = self.cdmiEntity.get_relationship_IsOwnerOf(genomes, [], [], ["id", "feature_type"])
        fidlist = self.getFieldFromRelationship(fiddict, "id", "to")
        typelist = self.getFieldFromRelationship(fiddict, "feature_type", "to")
        # We want protein-encoding genes only (not e.g. operons, operators, etc...)
        # The type of protein-encoding genes is CDS now but will possibly be changed to peg later...
        pegs = []
        for ii in range(len(fidlist)):
            if typelist[ii] == "peg" or typelist[ii] == "CDS":
                pegs.append(fidlist[ii])    
        return pegs
    
    def getOtuGenomeIds(self, count):
        ''' Query the CDMI for a list of OTU genome IDs.
    
            @param count: Number of entities to retrieve in each function call
            @return List of all OTU genome IDs, list of only prokaryote OTUs
        '''
    
        # Get the complete list of OTUs.
        otudict = dict()
        start = 0
        done = False
        while not done:
            subdict = self.cdmiEntity.all_entities_OTU(start, count, ["id"])
            otudict.update(subdict)
            start += count
            if len(subdict) < count:
                done = True
    
        # Find out if a OTU is marked as representative and if it is prokaryotic.
        otuids = self.getFieldFromEntity(otudict, "id")
        gendict = self.cdmiEntity.get_relationship_IsCollectionOf(otuids, [], ["representative"], ["id", "prokaryotic"])
        isrep = self.getFieldFromRelationship(gendict, "representative", "rel")
        isprok = self.getFieldFromRelationship(gendict, "prokaryotic", "to")
        genomeid = self.getFieldFromRelationship(gendict, "id", "to")
        prokotus = []
        otus = []
        for ii in range(len(genomeid)):
            if int(isrep[ii]) == 1 and int(isprok[ii]) == 1:
                prokotus.append(genomeid[ii])
            if int(isrep[ii]) == 1:
                otus.append(genomeid[ii])
        return otus, prokotus
    
    ################
    # Input: genomes - List of genome IDs
    # Output 1: Tuples of gene neighborhood information
    #   (contig_id, feature_id, start_location, strand)
    #   for each input OTU
    #   The output is sorted by contig id, then by start location
    #
    # Output 2: A dictionary from fid to roles
    ################
    def getGenomeNeighborhoodsAndRoles(genomes, config):
        self.cdmiEntity = CDMI_EntityAPI(config["cdmi_url"])
    
        pegs = genomesToPegs(genomes)
        # Get contigs
        fidlocdict = cdmiEntity.get_relationship_IsLocatedIn(pegs, [], ["begin", "dir"], ["id"])
        fids = self.getFieldFromRelationship(fidlocdict, "from_link", "rel")
        begins = self.getFieldFromRelationship(fidlocdict, "begin", "rel")
        dirs = self.getFieldFromRelationship(fidlocdict, "dir", "rel")
        cids = self.getFieldFromRelationship(fidlocdict, "id", "to")
    
        tuplist = []
        for ii in range(len(cids)):
            tuplist.append( (cids[ii], fids[ii], int(begins[ii]), dirs[ii]) )
        # Sort by contig first, then by start location.
        tuplist = sorted(tuplist, key=operator.itemgetter(0,2))
    
        # Now lets get the role for all of these IDs
        # Note that a single protein can have multiple roles.
        roledict = cdmiEntity.get_relationship_HasFunctional(fids, [], [], ["id"])
        fids = self.getFieldFromRelationship(roledict, "from_link", "rel")
        roles = self.getFieldFromRelationship(roledict, "id", "to")
        fidToRoles = {}
        rolesToFids = {}
        for ii in range(len(fids)):
            if fids[ii] in fidToRoles:
                fidToRoles[fids[ii]].append(roles[ii])
            else:
                fidToRoles[fids[ii]] = [ roles[ii] ]
            if roles[ii] in rolesToFids:
                rolesToFids[roles[ii]].append(fids[ii])
            else:
                rolesToFids[roles[ii]] = [ fids[ii] ]
        return tuplist, fidToRoles
    
    def mapComplexToRole(self, count):
        ''' Query the CDM for a list of links from complexes to roles.
    
            Only roles listed as "required" are included in the links.
    
            @note Someday the modeling service will provide this data from Mapping objects.
            @param count: Number of entities to retrieve in each function call
            @return Dictionary keyed by role of a list of complex IDs, dictionary keyed by
                complex ID to a list of roles.
        '''
    
        # Get a list of complexes
        cplxdict = dict()
        start = 0
        done = False
        while not done:
            subdict = self.cdmiEntity.all_entities_Complex(start, count, ["id"])
            cplxdict.update(subdict)
            start += count
            if len(subdict) < count:
                done = True
        cplxlist = self.getFieldFromEntity(cplxdict, "id")
    
        # Get a list of roles linked to those complexes
        roledict = self.cdmiEntity.get_relationship_IsTriggeredBy(cplxlist, [], ["optional"], ["id"])
        cplx = self.getFieldFromRelationship(roledict, "from_link", "rel")
        opt = self.getFieldFromRelationship(roledict, "optional", "rel")
        role = self.getFieldFromRelationship(roledict, "id", "to")
        complexToRequiredRoles = {}
        requiredRolesToComplex = {}
        for ii in range(len(cplx)):
            # For now - we don't want to deal with the "optional" components. I'm not sure how I'd incorporate them into a likelihood calculation anyway.
            if int(opt[ii]) == 1:
                continue
            # Note - this becomes an all-AND GPR - (role1 AND role2 AND ... )
            if cplx[ii] in complexToRequiredRoles:
                complexToRequiredRoles[cplx[ii]].append(role[ii])
            else:
                complexToRequiredRoles[cplx[ii]] = [ role[ii] ]
            if role[ii] in requiredRolesToComplex:
                requiredRolesToComplex[role[ii]].append(cplx[ii])
            else:
                requiredRolesToComplex[role[ii]] = [ cplx[ii] ]
        return complexToRequiredRoles, requiredRolesToComplex
    
    def mapReactionToComplex(self, count):
        ''' Query the CDM for a list of links from reactions to complexes.
    
            @note Someday the modeling service will provide this data from Mapping objects.
            @param count: Number of entities to retrieve in each function call
            @return Dictionary keyed by reaction ID to lists of complexes performing them,
                dictionary keyed by complex ID to list of reactions they perform.
        '''
    
        # The API was recently changed to use model IDs and to not use the reactions_to_complexes
        # but use the ER model instead.
        # I reflect that here...
        rxndict = dict()
        start = 0
        done = False
        while not done:
            subdict = self.cdmiEntity.all_entities_Reaction(start, count, ['id'])
            rxndict.update(subdict)
            start += count
            if len(subdict) < count:
                done = True
        rxns = self.getFieldFromEntity(rxndict, "id")
        cplxdict = self.cdmiEntity.get_relationship_IsStepOf(rxns, [], [], ["id"])
        rxnlist = self.getFieldFromRelationship(cplxdict, "from_link", "rel")
        cplxlist = self.getFieldFromRelationship(cplxdict, "id", "to")
        
        rxnToComplex = {}
        complexToRxn = {}
        for ii in range(len(rxnlist)):
            if rxnlist[ii] in rxnToComplex:
                rxnToComplex[rxnlist[ii]].append(cplxlist[ii])
            else:
                rxnToComplex[rxnlist[ii]] = [ cplxlist[ii] ]
            if cplxlist[ii] in complexToRxn:
                complexToRxn[cplxlist[ii]].append(rxnlist[ii])
            else:
                complexToRxn[cplxlist[ii]] = [ rxnlist[ii] ]
    
        return rxnToComplex, complexToRxn
