#!/usr/bin/python
#
# Call the CDMI API to get data needed to run the autorecon algorithm
# Needed data includes:
#
# For likelihood computations:
# - List of substems and the sequences of their members
# - List of OTU organisms
# - Gene neighborhoods for OTU organisms
#
# For optimization:
# - Reactions [only want mass/charge-balanced ones]
# - Metabolites
# -[Growth data? ]

# The CDMI_API is for "well-trodden paths" functions
# CDMI_EntityAPI is for ER functions (all_entities_..., get_Relationship_....)
from CDMI import CDMI_API, CDMI_EntityAPI
import urllib
import sys
import simplejson as json
import operator #for itemgetter
from PYTHON_GLOBALS import *

# Get field "fieldName" from the entity seedEntity, the result of a 
# get_entity_xxx or an all_entities_xxx call...
#
# The entities are dictionaries from ID to a dictionary of key-value pairs
# where keys are whatever you tell it you want.
def getFieldFromEntity(seedEntity, fieldName):
    if seedEntity is None:
        sys.stderr.write("INTERNAL ERROR: Provided seedEntity was None - usually this means you were searching for something that doesnt exist in the database\n")
        raise ValueError
    # Check for an error I seem to make all the time and yell at me in a USEFUL way
    if not isinstance(seedEntity, dict):
        sys.stderr.write("INTERNAL ERROR: getFieldFromEntity expects a dictionary - perhaps you meant to call getFieldFromRelationship?\n")
        raise ValueError
    f = []
    for entry in seedEntity:
        if fieldName not in seedEntity[entry]:
            sys.stderr.write("INTERNAL ERROR: Field name %s not found in provided entity\n" %(fieldName))
            raise ValueError
        f.append(seedEntity[entry][fieldName])
    return f

################
# INPUTS:
# seedRelationship: The result of one of the get_relationship_xxx functions
# fieldName: The field you want to extract from the object.
# objtype: "TO", "REL", or "FROM"
#
# OUTPUTS:
# A list (in the same order as the list from the get_relationship function)
# of the values with the specified field name.
#
# The get_relationship_xxx functions return lists of lists.
# The first list is a list of all the links
# The second list has three dictionaries in it: the TO dictionary, the REL dictionary and the FROM dictionary
#     describing properties on either end of the link and of the link itself...
#
#
# If you want to  maintain duplicicate relationships (many-to-one, one-to-many, many-to-many), this function should be called at
# least twice (once on each end of the relationship, or once on an end and once in the middle)..
#
##################
def getFieldFromRelationship(seedRelationship, fieldName, objtype):
    if seedRelationship is None:
        sys.stderr.write("INTERNAL ERROR: The provided relationship was None - usually this means you were searching for something that doesn't exist in the database.\n")
        raise ValueError
    objidx = None
    if objtype.lower() == "from":
        objidx = 0
    elif objtype.lower() == "rel":
        objidx = 1
    elif objtype.lower() == "to":
        objidx = 2
    else:
        sys.stderr.write("INTERNAL ERROR: In getFieldFromRelationship - objtype must be TO, REL, or FROM\n")
        raise ValueError
    if not isinstance(seedRelationship, list):
        sys.stderr.write("INTERNAL ERROR: getFieldFromRelationship expects a list - perhaps you meant to call getFieldFromEntity?\n")
        raise ValueError
    # Unravel
    f = []
    for entry in seedRelationship:
        # TO CHECK: Is it possible to have one of the links lead to nothing?
        # Check field name validity - if it links to something there has to be the data request there
        # or else something is wrong.
        if fieldName not in entry[objidx]:
            sys.stderr.write("INTERNAL ERROR: Field name %s not found in provided relationship\n" %(fieldName))
            raise ValueError
        f.append(entry[objidx][fieldName])
    return f

######
# Input: Minimum index (MINN) and count (COUNT) of items to extract
# Output: A list of fids in the specified subsystems
######
def subsystemFids(MINN, COUNT):
    cdmi = CDMI_API(URL)
    cdmi_entity = CDMI_EntityAPI(URL)
    # Get Genes that are in Subsystems and in Otus.
    ssdict = cdmi_entity.all_entities_Subsystem(MINN,COUNT,["id"])
    ssids = getFieldFromEntity(ssdict, "id")

    # Now lets get a list of FIDs within those subsystems
    ssfiddict = cdmi.subsystems_to_fids(ssids, [])
    ssfids = []
    for key in ssfiddict:
        for ssfid in ssfiddict[key]:
            ls = ssfiddict[key][ssfid]
            for arr in ls:
                if len(arr) > 1:
                    gl = arr[1]
                    for l in gl:
                        ssfids.append(l)
    # Uniquify!
    return list(set(ssfids))

#######
# NOTE - this is the same as the previous one BUT
# it tries to work around a bug in the current KBase implementation of
# subsystems_to_fids
# Namely, that although we can go from fids to subsystems, we cannot
# go back the other way to get the fids in some cases.
#
######
def subsystemFids_WORKAROUND(MINN, COUNT):
    cdmi = CDMI_API(URL)
    cdmi_entity = CDMI_EntityAPI(URL)
    # Our work-around is to go from subsystems to roles to fids and then back
    # to subsystems.
    # All of these links appear to be intact.
    ssdict = cdmi_entity.all_entities_Subsystem(MINN,COUNT,["id"])
    ssids = getFieldFromEntity(ssdict, "id")
    ssroledict = cdmi.subsystems_to_roles(ssids, [])
    print ssids
    print ssroledict
    ssroles = []
    exit(2)

#######
# Unput: Minimum index (MINN) and count (COUNT) of items to extract
# Output: List of FIDs with direct literature (dlit) evidence for their function.
#######
def getDlitFids(MINN, COUNT):
    cdmi = CDMI_API(URL)
    cdmi_entity = CDMI_EntityAPI(URL)
    pubdict = cdmi_entity.all_entities_Publication(MINN, COUNT, ["id"])
    pubids = getFieldFromEntity(pubdict, "id")
    pub2seq = cdmi_entity.get_relationship_Concerns(pubids, [], [], ["id"])
    pubseqs = getFieldFromRelationship(pub2seq, "id", "to")
    seq2fids = cdmi_entity.get_relationship_IsProteinFor(pubseqs, [], [], ["id"])
    fids = getFieldFromRelationship(seq2fids, "id", "to")
    return fids

######
#
# Input: 1: A list of fids (fidlist) 
#        2: A list of representative organisms (otus)
# Output: A list of FIDs containing only those organisms in the OTUs
# 
# fids with no organism are thrown away as are any that don't match an organism in orgidlist
######
def filterFidsByOtus(fidlist, otus):
    cdmi_entity = CDMI_EntityAPI(URL)

    # Identify the organism belonging to each fid
    # If this fails to find an organism we don't want it anyway...
    orgdict = cdmi_entity.get_relationship_IsOwnedBy(fidlist, [], [], ["id"])
    flist = getFieldFromRelationship(orgdict, "from_link", "rel")
    olist = getFieldFromRelationship(orgdict, "id", "to")

    fids = []
    for ii in range(len(olist)):
        if olist[ii] in otus:
            fids.append(flist[ii])
    return fids

######
# input: fidlist: A list of FIDs
# output 1: fidsToRoles: A dictionary from FID to a list of roles
# output 2: rolesToFids: A dictionary from a role to a list of FIDS (from the original list)
#
# fids with no roles are thrown away.
######
def fidsToRoles(fidlist):
    cdmi = CDMI_API(URL)
    cdmi_entity = CDMI_EntityAPI(URL)
    roledict = cdmi_entity.get_relationship_HasFunctional(fidlist, [], [], ["id"])

    flist = getFieldFromRelationship(roledict, "from_link", "rel")
    rolelist = getFieldFromRelationship(roledict, "id", "to")
    fidsToRoles = {}
    rolesToFids = {}
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
    # Convert back to lists to not break other functions.
    for f in fidsToRoles:
        fidsToRoles[f] = list(fidsToRoles[f])
    for r in rolesToFids:
        rolesToFids[r] = list(rolesToFids[r])
    return fidsToRoles, rolesToFids

#############
# Input: A list of feature IDs
# Output: A dictionary fid --> sequence
#
# features with no sequences are discarded
#############
def fidsToSequences(fidlist):
    cdmi = CDMI_API(URL)
    fidlist = list(set(fidlist))
    seqs = cdmi.fids_to_protein_sequences(fidlist)
    return seqs

def genomesToPegs(genomes):
    cdmi_entity = CDMI_EntityAPI(URL)
    fiddict = cdmi_entity.get_relationship_IsOwnerOf(genomes, [], [], ["id", "feature_type"])
    fidlist = getFieldFromRelationship(fiddict, "id", "to")
    typelist = getFieldFromRelationship(fiddict, "feature_type", "to")
    # We want protein-encoding genes only (not e.g. operons, operators, etc...)
    # The type of protein-encoding genes is CDS now but will possibly be changed to peg later...
    pegs = []
    for ii in range(len(fidlist)):
        if typelist[ii] == "peg" or typelist[ii] == "CDS":
            pegs.append(fidlist[ii])    
    return pegs

############
# Input: MINN - Minimum index and COUNT - number of items to extract
# Output: A dictionary from OTU ID to a list of genome IDs within the OTU...
############
def getOtuGenomeIds(MINN, COUNT):
    cdmi_entity = CDMI_EntityAPI(URL)

    # Get the requested number of OTU
    otudict = cdmi_entity.all_entities_OTU(MINN, COUNT, ["id"])
    otuids = getFieldFromEntity(otudict, "id")
    gendict = cdmi_entity.get_relationship_IsCollectionOf(otuids, [], ["representative"], ["id", "prokaryotic"])
    isrep = getFieldFromRelationship(gendict, "representative", "rel")
    isprok = getFieldFromRelationship(gendict, "prokaryotic", "to")
    genomeid = getFieldFromRelationship(gendict, "id", "to")
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
def getGenomeNeighborhoodsAndRoles(genomes):
    cdmi_entity = CDMI_EntityAPI(URL)

    pegs = genomesToPegs(genomes)
    # Get contigs
    fidlocdict = cdmi_entity.get_relationship_IsLocatedIn(pegs, [], ["begin", "dir"], ["id"])
    fids = getFieldFromRelationship(fidlocdict, "from_link", "rel")
    begins = getFieldFromRelationship(fidlocdict, "begin", "rel")
    dirs = getFieldFromRelationship(fidlocdict, "dir", "rel")
    cids = getFieldFromRelationship(fidlocdict, "id", "to")

    tuplist = []
    for ii in range(len(cids)):
        tuplist.append( (cids[ii], fids[ii], int(begins[ii]), dirs[ii]) )
    # Sort by contig first, then by start location.
    tuplist = sorted(tuplist, key=operator.itemgetter(0,2))

    # Now lets get the role for all of these IDs
    # Note that a single protein can have multiple roles.
    roledict = cdmi_entity.get_relationship_HasFunctional(fids, [], [], ["id"])
    fids = getFieldFromRelationship(roledict, "from_link", "rel")
    roles = getFieldFromRelationship(roledict, "id", "to")
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

#########
# Get a dictionary from roles to complexes and back
# Any functions that all form a complex are assumed to be AND-relationships
#
# Input: Minimum index (MINN) and count (COUNT)
# Output: Two dictionaries:
#    role --> complex
#    complex --> roles
#########
def complexRoleLinks(MINN, COUNT):
    cdmi_entity = CDMI_EntityAPI(URL)
    # Get a list of complexes
    cplxdict = cdmi_entity.all_entities_Complex(MINN, COUNT, ["id"])
    cplxlist = getFieldFromEntity(cplxdict, "id")
    # Get a list of roles linked to those complexes
    roledict = cdmi_entity.get_relationship_IsTriggeredBy(cplxlist, [], ["optional"], ["id"])
    cplx = getFieldFromRelationship(roledict, "from_link", "rel")
    opt = getFieldFromRelationship(roledict, "optional", "rel")
    role = getFieldFromRelationship(roledict, "id", "to")
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

############
# Get dictionaries from reactions to complexes
############
def reactionComplexLinks(MINN, COUNT):
    cdmi_entity = CDMI_EntityAPI(URL)

    # The API was recently changed to use model IDs and to not use the reactions_to_complexes
    # but use the ER model instead.
    # I reflect that here...
    rxndict = cdmi_entity.all_entities_Reaction(MINN, COUNT, ["id"])
    rxns = getFieldFromEntity(rxndict, "id")
    cplxdict = cdmi_entity.get_relationship_IsStepOf(rxns, [], [], ["id"])
    rxnlist = getFieldFromRelationship(cplxdict, "from_link", "rel")
    cplxlist = getFieldFromRelationship(cplxdict, "id", "to")
    
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
