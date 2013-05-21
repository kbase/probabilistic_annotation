#!/usr/bin/python

# Read and write data files
import os, sys, math, subprocess
from ConfigParser import ConfigParser

# E values of less than 1E-200 are treated as 1E-200 to avoid log of 0 issues.
MIN_EVALUE = 1E-200

# Exception thrown when makeblastdb command failed
class MakeblastdbError(Exception):
    pass

###########
#  OTUs   #
###########

def readOtuData(config):
    fid = open(os.path.join(config["data_folder_path"], config["otu_id_file"]), "r")
    otus = []
    prokotus = []
    for line in fid:
        spl = line.strip("\r\n").split("\t")
        otus.append(spl[0])
        if int(spl[1]) == 1:
            prokotus.append(spl[0])
    fid.close()
    return otus, prokotus

def writeOtuData(otus, prokotus, config):
    fid = open(os.path.join(config["data_folder_path"], config["otu_id_file"]), "w")
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

def readSubsystemFids(config):
    fid = open(os.path.join(config["data_folder_path"], config["subsystem_fid_file"]), "r")
    sub_fids = []
    for line in fid:
        spl = line.strip("\r\n")
        sub_fids.append(spl)
    fid.close()
    return sub_fids

def writeSubsystemFids(sub_fids, config):
    fid = open(os.path.join(config["data_folder_path"], config["subsystem_fid_file"]), "w")
    for f in sub_fids:
        fid.write("%s\n" %(f))
    fid.close()
    return

##################
# OTU FIDs       #
##################

def readDlitFids(config):
    fid = open(os.path.join(config["data_folder_path"], config["dlit_fid_file"]), "r")
    otu_fids = []
    for line in fid:
        spl = line.strip("\r\n")
        otu_fids.append(spl)
    fid.close()
    return otu_fids

def writeDlitFids(otu_fids, config):
    fid = open(os.path.join(config["data_folder_path"], config["dlit_fid_file"]), "w")
    for f in otu_fids:
        fid.write("%s\n" %(f))
    fid.close()
    return

###########################
# All FID roles           #
###########################

def readAllFidRoles(config):
    fid = open(os.path.join(config["data_folder_path"], config["concatenated_fid_role_file"]), "r")
    all_fidsToRoles = {}
    for line in fid:
        spl = line.strip("\r\n").split("\t")
        roles = spl[1].split(config["separator"])
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

def writeAllFidRoles(otu_fidsToRoles, config):
    fid = open(os.path.join(config["data_folder_path"], config["concatenated_fid_role_file"]), "w")
    for f in otu_fidsToRoles:
        fid.write("%s\t%s\n" %(f, config["separator"].join(otu_fidsToRoles[f])))
    fid.close()
    return

######################
# Filtered OTU roles #
######################

def readFilteredOtuRoles(config):
    fid = open(os.path.join(config["data_folder_path"], config["subsystem_otu_fid_roles_file"]), "r")
    otu_fidsToRoles = {}
    for line in fid:
        spl = line.strip("\r\n").split("\t")
        roles = spl[1].split(config["separator"])
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

def writeFilteredOtuRoles(otu_fidsToRoles, config):
    fid = open(os.path.join(config["data_folder_path"], config["subsystem_otu_fid_roles_file"]), "w")
    for f in otu_fidsToRoles:
        fid.write("%s\t%s\n" %(f, config["separator"].join(otu_fidsToRoles[f])))
    fid.close()
    return

########################
# Subsystem FASTA file #
########################

def readSubsystemFasta(config):
    fid = open(os.path.join(config["data_folder_path"], config["subsystem_otu_fasta_file"]), "r")
    fid.close()
    return

def writeSubsystemFasta(fidsToSeqs, config):
    filepath = os.path.join(config["data_folder_path"], config["subsystem_otu_fasta_file"])
    fid = open(filepath, "w")
    for fids in fidsToSeqs:
        fid.write(">%s\n%s\n" %(fids, fidsToSeqs[fids]))
    fid.close()
    # Compile the BLAST database for the fasta file
    args = ["/usr/bin/makeblastdb", "-in", filepath, "-dbtype", "prot"]
    try:
        retcode = subprocess.call(args)
        if retcode < 0:
            cmd = ' '.join(args)
            raise MakeblastdbError("'%s' was terminated by signal %d" %(cmd, -retcode))
        else:
            if retcode > 0:
                cmd = ' '.join(args)
                raise MakeblastdbError("'%s' failed with status %d" %(cmd, retcode))
    except OSError as e:
        cmd = ' '.join(args)
        raise MakeblastdbError("Failed to run '%s': %s" %(cmd, e.strerror))
    return

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

def readComplexRoles(config):
    fid = open(os.path.join(config["data_folder_path"], config["complexes_roles_file"]), "r")
    complexToRequiredRoles = {}
    for line in fid:
        spl = line.strip("\r\n").split("\t")
        complexes = spl[0]
        roles = spl[1].split(config["separator"])
        # This shouldn't be necessary but just to be safe...
        if complexes in complexToRequiredRoles:
            complexToRequiredRoles[complexes] += roles
        else:
            complexToRequiredRoles[complexes]  = roles
    fid.close()
    return complexToRequiredRoles

def writeComplexRoles(complexToRequiredRoles, config):
    fid = open(os.path.join(config["data_folder_path"], config["complexes_roles_file"]), "w")
    for complexes in complexToRequiredRoles:
        fid.write("%s\t%s\n" %(complexes, config["separator"].join(complexToRequiredRoles[complexes])))
    fid.close()
    return

#########################
# Reaction --> complex  #
#########################

def readReactionComplex(config):
    fid = open(os.path.join(config["data_folder_path"], config["reaction_complexes_file"]), "r")
    rxnToComplexes = {}
    for line in fid:
        spl = line.strip("\r\n").split("\t")
        rxn = spl[0]
        cplxlist = spl[1].split(config["separator"])
        # This shouldn't be necessary but just to be safe...
        if rxn in rxnToComplexes:
            rxnToComplexes[rxn] += cplxlist
        else:
            rxnToComplexes[rxn] = cplxlist
    fid.close()
    return rxnToComplexes

def writeReactionComplex(rxnToComplexes, config):
    fid = open(os.path.join(config["data_folder_path"], config["reaction_complexes_file"]), "w")
    for rxn in rxnToComplexes:
        fid.write("%s\t%s\n" %(rxn, config["separator"].join(rxnToComplexes[rxn])))
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

def getConfig(filename):
    retconfig = {}
    config = ConfigParser()
    config.read(filename)
    for nameval in config.items("probabilistic_annotation"):
        retconfig[nameval[0]] = nameval[1]
    return retconfig
