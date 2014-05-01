#!/usr/bin/python

# Read and write data files
import os, sys, math, subprocess, time
from ConfigParser import ConfigParser

# E values of less than 1E-200 are treated as 1E-200 to avoid log of 0 issues.
MIN_EVALUE = 1E-200

# File names for static database files
DatabaseFiles = {
    "otu_id_file": "OTU_GENOME_IDS",
    "subsystem_fid_file": "SUBSYSTEM_FIDS",
    "dlit_fid_file": "DLIT_FIDS",
    "concatenated_fid_file": "ALL_FIDS",
    "concatenated_fid_role_file": "ALL_FID_ROLES",
    "subsystem_otu_fid_roles_file": "SUBSYSTEM_OTU_FID_ROLES",
    "subsystem_otu_fasta_file": "SUBSYSTEM_FASTA",
    "subsystem_otu_index_file": "SUBSYSTEM_FASTA.pin",
    "subsystem_otu_sequence_file": "SUBSYSTEM_FASTA.psq",
    "subsystem_otu_header_file": "SUBSYSTEM_FASTA.phr",
    "subsystem_udb_file": "SUBSYSTEM.udb",
    "complexes_roles_file": "COMPLEXES_ROLES",
    "reaction_complexes_file": "REACTIONS_COMPLEXES"
}

# File names for tracking static database files
StatusFiles = {
    "status_file": "staticdata.status",
    "cache_file": "staticdata.cache"    
}
# Exception thrown when makeblastdb command failed
class MakeblastdbError(Exception):
    pass

# Exception thrown when static database files are not ready
class NotReadyError(Exception):
    pass

###########
#  OTUs   #
###########

# The OTU ID file is a list of representative OTU genome IDs.  Each line has these fields:
#   1. Genome ID in KBase format (e.g. kb|g.0)
#   2. Flag indicating if the genome is a prokaryote (1 means yes, 0 means no)

''' Read data from the representative OTU genome ID file.

    @param config Dictionary of configuration variables
    @returns List of all OTU genome IDs, list of prokaryote OTU genome IDs
'''

def readOtuData(config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["otu_id_file"]), "r")
    otus = []
    prokotus = []
    for line in fid:
        spl = line.strip("\r\n").split("\t")
        otus.append(spl[0])
        if int(spl[1]) == 1:
            prokotus.append(spl[0])
    fid.close()
    return otus, prokotus

''' Write data to the representative OTU genome ID file.

    @param otus List of all OTU genome IDs
    @param prokouts List of prokaryote OTU genome IDs
    @param config Dictionary of configuration variables
    @returns Nothing
'''

def writeOtuData(otus, prokotus, config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["otu_id_file"]), "w")
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

# The subsystem feature ID file is a list of feature IDs from SEED subsystems.  Each line has
# one field that is the feature ID in KBase format (e.g. kb|g.3.peg.541).

''' Read data from the subsystem feature ID file.

    @param config Dictionary of configuration variables
    @returns List of feature IDs from SEED subystems
'''

def readSubsystemFids(config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_fid_file"]), "r")
    sub_fids = []
    for line in fid:
        spl = line.strip("\r\n")
        sub_fids.append(spl)
    fid.close()
    return sub_fids

''' Write data to the subsystem feature ID file.

    @param sub_fids List of feature IDs from SEED subsystems
    @param config Dictionary of configuration variables
    @returns Nothing
'''

def writeSubsystemFids(sub_fids, config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_fid_file"]), "w")
    for f in sub_fids:
        fid.write("%s\n" %(f))
    fid.close()
    return

##################
# DLIT FIDs      #
##################

# The direct literature-supported feature ID file is a list of feature IDs identified in the
# literature.  Each line has one field that is the feature ID in KBase format (e.g. kb|g.428.peg.6254).

''' Read data from the direct literature-supported feature ID file.

    @param config Dictionary of configuration variables
    @returns List of feature IDs from literature
'''

def readDlitFids(config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["dlit_fid_file"]), "r")
    otu_fids = []
    for line in fid:
        spl = line.strip("\r\n")
        otu_fids.append(spl)
    fid.close()
    return otu_fids

''' Write data to the direct literature-supported feature ID file.

    @param otu_fids List of feature IDs from literature
    @param config Dictionary of configuration variables
    @returns Nothing
'''

def writeDlitFids(otu_fids, config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["dlit_fid_file"]), "w")
    for f in otu_fids:
        fid.write("%s\n" %(f))
    fid.close()
    return

###########################
# All FID roles           #
###########################

# The concatenated feature ID to role file is a mapping of feature IDs to functional roles.
# Each line has these fields:
#   1. Feature ID in KBase format (e.g. kb|g.0.peg.2094)
#   2. List of names of functional roles (e.g. Conserved ATP-binding protein YghS)
#
# Note that functional roles must be separated by a string that does not occur in any role.

''' Read data from the concatenated feature ID to role file.

    @param config Dictionary of configuration variables
    @returns Dictionary mapping a feature ID to list of names of roles, dictionary mapping a role to feature ID
'''

def readAllFidRoles(config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["concatenated_fid_role_file"]), "r")
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

''' Write data to the concatenated feature ID to role file.

    @param otu_fidsToRoles Dictionary mapping a feature ID to list of names of roles
    @param config Dictionary of configuration variables
    @returns Nothing
'''

def writeAllFidRoles(otu_fidsToRoles, config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["concatenated_fid_role_file"]), "w")
    for f in otu_fidsToRoles:
        fid.write("%s\t%s\n" %(f, config["separator"].join(otu_fidsToRoles[f])))
    fid.close()
    return

######################
# Filtered OTU roles #
######################

# The filtered feature ID to roles file is a mapping of feature IDs to functional roles where
# there is one protein from each OTU for each functional role.  Each line has these fields:
#   1. Feature ID in KBase format (e.g. kb|g.0.peg.2094)
#   2. List of names of functional roles (e.g. Conserved ATP-binding protein YghS)
#
# Note that functional roles must be separated by a string that does not occur in any role.

''' Read data from the filtered feature ID to roles file.

    @param config Dictionary of configuration variables
    @returns Dictionary mapping a feature ID to list of names of roles, dictionary mapping a role to feature ID
'''

def readFilteredOtuRoles(config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_otu_fid_roles_file"]), "r")
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

''' Write data to the filtered feature ID to roles file.

    @param otu_fidsToRoles Dictionary mapping a feature ID to list of names of roles
    @param config Dictionary of configuration variables
    @returns Nothing
'''

def writeFilteredOtuRoles(otu_fidsToRoles, config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_otu_fid_roles_file"]), "w")
    for f in otu_fidsToRoles:
        fid.write("%s\t%s\n" %(f, config["separator"].join(otu_fidsToRoles[f])))
    fid.close()
    return

########################
# Subsystem FASTA file #
########################

# The subsystem FASTA file contains the amino acid sequences for a set of feature IDs.  When the file
# is written a BLAST database is automatically generated.

''' Read data from the subsystem FASTA file.

    @note This function does not return data (not sure it is really needed except for symmetry).
    @param config Dictionary of configuration variables
    @returns Nothing
'''

def readSubsystemFasta(config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_otu_fasta_file"]), "r")
    fid.close()
    return

''' Write data to the subsystem FASTA file and build BLAST database.

    @param fidsToSeqs Dictionary mapping a feature ID to amino acid sequence
    @param config Dictionary of configuration variables
    @returns Nothing
'''

def writeSubsystemFasta(fidsToSeqs, config):
    filepath = os.path.join(config["data_folder_path"], DatabaseFiles["subsystem_otu_fasta_file"])
    fid = open(filepath, "w")
    # Sort the fids so that fasta files containing the same proteins hash to the same MD5 (for
    # data provenance purposes)
    for fids in sorted(fidsToSeqs.keys()):
        fid.write(">%s\n%s\n" %(fids, fidsToSeqs[fids]))
    fid.close()

    # Build the command based on the configured search program.
    if config['search_program'] == 'usearch':
        args = [ config['search_program_path'], '-makeudb_ublast', filepath, '-output', os.path.join(config['data_folder_path'], DatabaseFiles['subsystem_udb_file']) ]
    else:
        args = [ "/usr/bin/makeblastdb", "-in", filepath, "-dbtype", "prot" ]

    # Run the command to compile the database from the subsystem fasta file.
    try:
        proc = subprocess.Popen(args, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        (stdout, stderr) = proc.communicate()
        if proc.returncode < 0:
            cmd = ' '.join(args)
            raise MakeblastdbError("'%s' was terminated by signal %d" %(cmd, -proc.returncode))
        else:
            if proc.returncode > 0:
                cmd = ' '.join(args)
                details = "'%s' failed with return code %d:\nCommand: '%s'\nStdout: '%s'\nStderr: '%s'" \
                    %(args[0], proc.returncode, cmd, stdout, stderr)
                raise MakeblastdbError(details)
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

# The complex to roles file contains a mapping of complex IDs to functional roles.
# Each line has these fields:
#   1. Complex ID in KBase format (e.g. kb|cpx.51240)
#   2. List of names of functional roles (e.g. V-type ATP synthase subunit H (EC 3.6.3.14))
#
# Note that functional roles must be separated by a string that does not occur in any role.

''' Read data from the complex to roles file.

    @param config Dictionary of configuration variables
    @returns Dictionary mapping a complex ID to list of names of functional roles
'''

def readComplexRoles(config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["complexes_roles_file"]), "r")
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

''' Write data to the complex to roles file.

    @param complexToRequiredRoles Dictionary mapping a complex ID to list of names of functional roles
    @param config Dictionary of configuration variables
    @returns Nothing
'''

def writeComplexRoles(complexToRequiredRoles, config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["complexes_roles_file"]), "w")
    for complexes in complexToRequiredRoles:
        fid.write("%s\t%s\n" %(complexes, config["separator"].join(complexToRequiredRoles[complexes])))
    fid.close()
    return

#########################
# Reaction --> complex  #
#########################

# The reaction to complexes file contains a mapping of reaction IDs to complex IDs.
# Each line has these fields:
#   1. Reaction ID in KBase format (e.g. kb|rxn.5682)
#   2. List of complex IDs in KBase format (e.g. kb|cpx.1507///kb|cpx.1813)
#
# Note that complex IDs must be separated by a string that does not occur in any complex ID.

''' Read data from the reaction to complexes file.

    @param config Dictionary of configuration variables
    @returns Dictionary mapping a reaction ID to list of complex IDs
'''

def readReactionComplex(config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["reaction_complexes_file"]), "r")
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

''' Write data to reaction to complexes file.

    @param rxnToComplexes Dictionary mapping a reaction ID to list of complex IDs
    @param config Dictionary of configuration variables
    @returns Nothing
'''

def writeReactionComplex(rxnToComplexes, config):
    fid = open(os.path.join(config["data_folder_path"], DatabaseFiles["reaction_complexes_file"]), "w")
    for rxn in rxnToComplexes:
        fid.write("%s\t%s\n" %(rxn, config["separator"].join(rxnToComplexes[rxn])))
    fid.close()
    return

''' Read BLAST results file and store in a convenient structure.

    @note Score is the negative log E-value
    @param blast_result_file Path to BLAST results file
    @returns Dictionary mapping query ID to tuple of target ID and score
'''

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

''' Read the roleset probability file.

    @param roleset_probability_file Path to roleset probability file
    @returns Dictionary mapping query ID to list of tuples with list of roles and probability
'''

def readRolesetProbabilityFile(roleset_probability_file):
    queryToTuplist = {}
    for line in open(roleset_probability_file, "r"):
        spl = line.strip("\r\n").split("\t")
        if spl[0] in queryToTuplist:
            queryToTuplist[spl[0]].append( (spl[1], float(spl[2])) )
        else:
            queryToTuplist[spl[0]] = [ (spl[1], float(spl[2])) ]
    return queryToTuplist

''' Read a configuration file.

    If filename is None, the default configuration file is read.  If the
    probabilistic_annotation section is not present in the configuration
    file, it is automatically added.

    @param filename Path to configuration file.
    @returns Config object
'''

def readConfig(filename=None):
    # Use default config file if one is not specified.
    if filename == None:
        filename = os.path.join(os.environ['KB_TOP'], 'deployment.cfg')

    # Read the config file.
    config = ConfigParser()
    try:
        config.read(filename)
    except Exception as e:
        print "Error while reading config file %s: %s" % (filename, e)

    # Make sure there is a probabilistic_annotation section in the config file.
    if not config.has_section('probabilistic_annotation'):
        config.add_section('probabilistic_annotation')
        with open(filename, 'w') as configfile:
            config.write(configfile)

    return config

''' Get the probabilistic_annotation section from a configuration file.

    If filename is None, the default configuration file is read.

    @param filename Path to configuration file.
    @returns Dictionary mapping configuration variables to values
'''

def getConfig(filename=None):
    # Read the config file.
    config = readConfig(filename)

    # Extract the probabilistic annotation section.
    sectionConfig = dict()
    for nameval in config.items("probabilistic_annotation"):
        sectionConfig[nameval[0]] = nameval[1]
    return sectionConfig

# The status file is used to track the status of setting up the static database files when
# the server starts.  The first line of the file contains the status which is one of
# these values:
#   1. 'building' when the pa-gendata command is building the files
#   2. 'running' when the server initialization is in progress
#   3. 'ready' when the server initialization is complete or a build is complete
#   4. 'failed' when there was an error building or loading the files
#
# The second line has the timestamp of when the status was last changed.

''' Read the current status value from the status file.

    @param config Dictionary of configuration variables
    @returns Current status string
'''

def readStatusFile(config):
    fid = open(os.path.join(config["data_folder_path"], StatusFiles["status_file"]), "r")
    statusLine = fid.readline()
    fid.close()
    return statusLine.strip("\r\n")

''' Write new status value to the status file.

    @param config Dictionary of configuration variables
    @param status New status value
    @returns Nothing
'''

def writeStatusFile(config, status):
    fid = open(os.path.join(config["data_folder_path"], StatusFiles["status_file"]), "w")
    fid.write("%s\ncompleted at %s\n" %(status, time.strftime("%a %b %d %Y %H:%M:%S %Z", time.localtime())))
    fid.close()
    return

''' Check for existence of all of the database files.

    @param data_folder_path Path to directory containing static database files
    @raise NotReadyError: A database file does not exist
    @returns Nothing
'''

def checkIfDatabaseFilesExist(data_folder_path):
    for key in DatabaseFiles:
        localPath = os.path.join(data_folder_path, DatabaseFiles[key])
        if not os.path.exists(localPath):
            raise NotReadyError("Static database file '%s' does not exist" %(localPath))
    return


