#!/usr/bin/python

# Read and write data files
from PYTHON_GLOBALS import *
import os
import hashlib
import sys
import math

PSEUDOCOUNT=40
MIN_EVALUE = 1E-200 #E values of less than 1E-200 are treated as 1E-200 to avoid log of 0 issues.

###########
#  OTUs   #
###########

def readOtuData(folder):
    fid = open(os.path.join(folder, OTU_ID_FILE), "r")
    otus = []
    prokotus = []
    for line in fid:
        spl = line.strip("\r\n").split("\t")
        otus.append(spl[0])
        if int(spl[1]) == 1:
            prokotus.append(spl[0])
    fid.close()
    return otus, prokotus

def writeOtuData(otus, prokotus, folder):
    fid = open(os.path.join(folder, OTU_ID_FILE), "w")
    for otu in otus:
        if otu in prokotus:
            fid.write("%s\t%d\n" %(otu, 1))
        else:
            fid.write("%s\t%d\n" %(otu, 0))
    fid.close()
    return

##################
# Subsystem FIDs #
##################

def readSubsystemFids(folder):
    fid = open(os.path.join(folder, SUBSYSTEM_FID_FILE), "r")
    sub_fids = []
    for line in fid:
        spl = line.strip("\r\n")
        sub_fids.append(spl)
    fid.close()
    return sub_fids

def writeSubsystemFids(sub_fids, folder):
    fid = open(os.path.join(folder, SUBSYSTEM_FID_FILE), "w")
    for f in sub_fids:
        fid.write("%s\n" %(f))
    fid.close()
    return

##################
# OTU FIDs       #
##################

def readDlitFids(folder):
    fid = open(os.path.join(folder, DLIT_FID_FILE), "r")
    otu_fids = []
    for line in fid:
        spl = line.strip("\r\n")
        otu_fids.append(spl)
    fid.close()
    return otu_fids

def writeDlitFids(otu_fids, folder):
    fid = open(os.path.join(folder, DLIT_FID_FILE), "w")
    for f in otu_fids:
        fid.write("%s\n" %(f))
    fid.close()
    return

###########################
# All FID roles           #
###########################

def readAllFidRoles(folder):
    fid = open(os.path.join(folder, CONCATINATED_FID_ROLE_FILE), "r")
    all_fidsToRoles = {}
    for line in fid:
        spl = line.strip("\r\n").split("\t")
        roles = spl[1].split(SEPARATOR)
        if spl[0] in all_fidsToRoles:
            all_fidsToRoles[spl[0]] += roles
        else:
            all_fidsToRoles[spl[0]] = roles
    fid.close()

    all_rolesToFids = {}
    for fid in all_fidsToRoles:
        roles = all_fidsToRoles[fid]
        for role in roles:
            if role in all_rolesToFids:
                all_rolesToFids[role].append(fid)
            else:
                all_rolesToFids[role] = [ fid ]

    return all_fidsToRoles, all_rolesToFids

def writeAllFidRoles(otu_fidsToRoles, folder):
    fid = open(os.path.join(folder, CONCATINATED_FID_ROLE_FILE), "w")
    for f in otu_fidsToRoles:
        fid.write("%s\t%s\n" %(f, SEPARATOR.join(otu_fidsToRoles[f])))
    fid.close()
    return

######################
# Filtered OTU roles #
######################

def readFilteredOtuRoles(folder):
    fid = open(os.path.join(folder, SUBSYSTEM_OTU_FID_ROLES_FILE), "r")
    otu_fidsToRoles = {}
    for line in fid:
        spl = line.strip("\r\n").split("\t")
        roles = spl[1].split(SEPARATOR)
        if spl[0] in otu_fidsToRoles:
            otu_fidsToRoles[spl[0]] += roles
        else:
            otu_fidsToRoles[spl[0]] = roles
    fid.close()

    otu_rolesToFids = {}
    for fid in otu_fidsToRoles:
        roles = otu_fidsToRoles[fid]
        for role in roles:
            if role in otu_rolesToFids:
                otu_rolesToFids[role].append(fid)
            else:
                otu_rolesToFids[role] = [ fid ]

    return otu_fidsToRoles, otu_rolesToFids

def writeFilteredOtuRoles(otu_fidsToRoles, folder):
    fid = open(os.path.join(folder, SUBSYSTEM_OTU_FID_ROLES_FILE), "w")
    for f in otu_fidsToRoles:
        fid.write("%s\t%s\n" %(f, SEPARATOR.join(otu_fidsToRoles[f])))
    fid.close()
    return

########################
# Subsystem FASTA file #
########################

def readSubsystemFasta(folder):
    fid = open(os.path.join(folder, SUBSYSTEM_OTU_FASTA_FILE), "r")
    fid.close()
    return

def writeSubsystemFasta(fidsToSeqs, folder):
    '''
    Write a FASTA for all the subsystem OTU genes, sort them,
    and compute an MD5. Then compile them (using makeblastdb) into a
    BLAST database.

    The location of the FASTA file is folder/SUBSYSTEM_OTU_FASTA_FILE
    where SUBSYSTEM_OTU_FASTA_FILE is defined in the config file.

    Returns the MD5 of the FASTA file.
    '''
    fid = open(os.path.join(folder, SUBSYSTEM_OTU_FASTA_FILE), "w")

    # Sort the proteins by ID so that if we have the same list multiple times 
    # we get the same hash for the FASTA file output.
    for fids in sorted(fidsToSeqs.keys()):
        fid.write(">%s\n%s\n" %(fids, fidsToSeqs[fids]))
    fid.close()
    # Compute the MD5 hash of the FASTA file
    f = open(os.path.join(folder, SUBSYSTEM_OTU_FASTA_FILE), "r")
    md5 = hashlib.md5()
    block_size = 65536
    while True:
        data = f.read(block_size)
        if not data:
            break
        md5.update(data)
    fasta_md5 = md5.hexdigest()

    # Compile the BLAST database for the fasta file
    status = os.system("makeblastdb -in %s -dbtype prot > /dev/null" %(os.path.join(folder, SUBSYSTEM_OTU_FASTA_FILE)))
    # TODO: Throw exceptions for errors
    if os.WIFEXITED(status):
        if os.WEXITSTATUS(status) != 0:
            print("makeblastdb failed with status %d\n" %(os.WEXITSTATUS(status)))
    if os.WIFSIGNALED(status):
        print("makeblastdb ended by signal %d\n" %(os.WTERMSIG(status)))

    return fasta_md5

#####################
# OTU neighborhoods #
#####################

#def readOtuNeighborhoods(folder):
#    fid = open(os.path.join(folder, OTU_NEIGHBORHOOD_FILE), "r")##

#    tuplist = []
#    fidToRoles = {}
#    for line in fid:
#        spl = line.strip("\r\n").split("\t")
#        roles = spl[4].split(SEPARATOR)
#        if spl[1] in fidToRoles:
#            fidToRoles[spl[1]] += roles
#        else:
#            fidToRoles[spl[1]] = roles
#        tuplist.append( (spl[0], spl[1], spl[2], spl[3],) )
#    fid.close()
#    return tuplist, fidToRoles

#def writeOtuNeighborhoods(tuplist, fidToRoles, verbose, fname):
#    fid = open(fname, "w")
#    for f in tuplist:
#        if f[1] in fidToRoles:
#            roles = fidToRoles[f[1]]
#        else:
#            if verbose:
#                sys.stderr.write("WARNING: Fid %s has no role despite being a neighbor of an OTU gene!\n" %(f[1]) )
#            roles = ""
#        try:
#            fid.write("%s\t%s\t%s\t%s\t%s\n" %(f[0], f[1], f[2], f[3], SEPARATOR.join(roles)))
#        except UnicodeEncodeError:
#            sys.stderr.write("ERROR: encountered roles that contain non-ASCII characters?\n")
#            sys.stderr.write("In gene ID %s\n" %(f[1]))
#            sys.stderr.write("Skipping...\n")
#            continue
#    fid.close()
#    return

#####################
# Complex --> roles #
#####################

def readComplexRoles(folder):
    fid = open(os.path.join(folder, COMPLEXES_ROLES_FILE), "r")
    complexToRequiredRoles = {}
    for line in fid:
        spl = line.strip("\r\n").split("\t")
        complexes = spl[0]
        roles = spl[1].split(SEPARATOR)
        # This shouldn't be necessary but just to be safe...
        if complexes in complexToRequiredRoles:
            complexToRequiredRoles[complexes] += roles
        else:
            complexToRequiredRoles[complexes]  = roles
    fid.close()
    return complexToRequiredRoles

def writeComplexRoles(complexToRequiredRoles, folder):
    fid = open(os.path.join(folder, COMPLEXES_ROLES_FILE), "w")
    for complexes in complexToRequiredRoles:
        fid.write("%s\t%s\n" %(complexes, SEPARATOR.join(complexToRequiredRoles[complexes])))
    fid.close()
    return

#########################
# Reaction --> complex  #
#########################

def readReactionComplex(folder):
    fid = open(os.path.join(folder, REACTION_COMPLEXES_FILE), "r")
    rxnToComplexes = {}
    for line in fid:
        spl = line.strip("\r\n").split("\t")
        rxn = spl[0]
        cplxlist = spl[1].split(SEPARATOR)
        # This shouldn't be necessary but just to be safe...
        if rxn in rxnToComplexes:
            rxnToComplexes[rxn] += cplxlist
        else:
            rxnToComplexes[rxn] = cplxlist
    fid.close()
    return rxnToComplexes

def writeReactionComplex(rxnToComplexes, folder):
    fid = open(os.path.join(folder, REACTION_COMPLEXES_FILE), "w")
    for rxn in rxnToComplexes:
        fid.write("%s\t%s\n" %(rxn, SEPARATOR.join(rxnToComplexes[rxn])))
    fid.close()
    return

# Read in the BLAST results file and store the results in a convenient structure
# Query ID --> [ (target ID, score) ]
#
# Score is the negative log-E value
def parseBlastOutput(blast_result_file):
    idToTargetList = {}
    for line in open(blast_result_file, "r"):
        spl = line.strip("\r\n").split("\t")
        queryid = spl[0]
        targetid = spl[1]
        logeval = -1.0 * math.log10(float(spl[10]) + MIN_EVALUE)
        tup = ( targetid, logeval )
        if queryid in idToTargetList:
            idToTargetList[queryid].append( tup )
        else:
            idToTargetList[queryid] = [ tup ]
    return idToTargetList

# Read the roleset probability file and returns a dictionary
# Query --> list of (rolelist, probability)
def readRolesetProbabilityFile(roleset_probability_file):
    queryToTuplist = {}
    for line in open(roleset_probability_file, "r"):
        spl = line.strip("\r\n").split("\t")
        if spl[0] in queryToTuplist:
            queryToTuplist[spl[0]].append( (spl[1], float(spl[2])) )
        else:
            queryToTuplist[spl[0]] = [ (spl[1], float(spl[2])) ]
    return queryToTuplist
