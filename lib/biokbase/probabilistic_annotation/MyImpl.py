import os, sys
from biokbase.probabilistic_annotation.DataExtractor import *
from biokbase.probabilistic_annotation.DataParser import *
from biokbase.probabilistic_annotation.PYTHON_GLOBALS import *


def generate_data(params):
    
    # When regenerating or deleting the data, remove all of the files first.
    if params.regenerate or params.delete_only:
        safeRemove(OTU_ID_FILE, params.folder)
        safeRemove(SUBSYSTEM_FID_FILE, params.folder)
        safeRemove(DLIT_FID_FILE, params.folder)
        safeRemove(CONCATINATED_FID_FILE, params.folder)
        safeRemove(SUBSYSTEM_OTU_FID_ROLES_FILE, params.folder)
        safeRemove(SUBSYSTEM_OTU_FASTA_FILE, params.folder)
        safeRemove(SUBSYSTEM_OTU_FASTA_FILE + ".psq", params.folder) 
        safeRemove(SUBSYSTEM_OTU_FASTA_FILE + ".pin", params.folder)
        safeRemove(SUBSYSTEM_OTU_FASTA_FILE + ".phr", params.folder)
        safeRemove(COMPLEXES_ROLES_FILE, params.folder)
        safeRemove(REACTION_COMPLEXES_FILE, params.folder)
    
    # Our job is done if all we want to do is delete files.
    if params.delete_only:
        return(1)
    
    sys.stderr.write("Generating requested data:....\n")
    
    # Get lists of OTUs
    sys.stderr.write("OTU data...")
    try:
        if params.verbose:
            sys.stderr.write("reading from file...")
        otus, prokotus = readOtuData(params.folder)
    except IOError:
        if params.verbose:
            sys.stderr.write("failed...generating file...")
        otus, prokotus = getOtuGenomeIds(MINN, COUNT)
        writeOtuData(otus, prokotus, params.folder)
    sys.stderr.write("done\n")
    
    # Get a list of subsystem FIDs
    sys.stderr.write("List of subsystem FIDS...")
    try:
        if params.verbose:
            sys.stderr.write("reading from file...")
        sub_fids = readSubsystemFids(params.folder)
    except IOError:
        if params.verbose:
            sys.stderr.write("failed...generating file...")
        sub_fids = subsystemFids(MINN, COUNT)
        writeSubsystemFids(sub_fids, params.folder)
    sys.stderr.write("done\n")
    
    # Get a list of Dlit FIDSs
    # We include these because having them greatly expands the
    # number of roles for which we have representatives.
    sys.stderr.write("Getting a list of DLit FIDs...")
    try:
        if params.verbose:
            sys.stderr.write("reading from file...")
        dlit_fids = readDlitFids(params.folder)
    except IOError:
        if params.verbose:
            sys.stderr.write("failed...generating file...")
        dlit_fids = getDlitFids(MINN, COUNT)
        writeDlitFids(dlit_fids, params.folder)
    sys.stderr.write("done\n")
    
    # Concatenate the two FID lists before filtering
    # (Note - doing so after would be possible as well but
    # can lead to the same kinds of biases as not filtering
    # the subsystems... I'm not sure the problem would
    # be as bad for these though)
    sys.stderr.write("Combining lists of subsystem and DLit FIDS...")
    fn = os.path.join(params.folder, CONCATINATED_FID_FILE)
    try:
        if params.verbose:
            sys.stderr.write("reading from file...")
        all_fids = set()
        for line in open(fn, "r"):
            all_fids.add(line.strip("\r\n"))
        all_fids = list(all_fids)
    except IOError:
        if params.verbose:
            sys.stderr.write("failed...generating file...")
        all_fids = list(set(sub_fids + dlit_fids))
        f = open(fn, "w")
        for fid in all_fids:
            f.write("%s\n" %(fid))
        f.close()
    sys.stderr.write("done\n")
    
    # Identify roles for the OTU genes
    sys.stderr.write("Roles for un-filtered list...")
    try:
        if params.verbose:
            sys.stderr.write("reading from file...")
        all_fidsToRoles, all_rolesToFids = readAllFidRoles(params.folder)
    except IOError:
        if params.verbose:
            sys.stderr.write("failed...generating file...")
        all_fidsToRoles, all_rolesToFids = fidsToRoles(all_fids)
        writeAllFidRoles(all_fidsToRoles, params.folder)
    sys.stderr.write("done\n")
    
    # Filter the subsystem FIDs by organism. We only want OTU genes.
    # Unlike the neighborhood analysis, we don't want to include only 
    # prokaryotes here.
    sys.stderr.write("Filtered list by OTUs...")
    try:
        if params.verbose:
            sys.stderr.write("reading from file...")
        otu_fidsToRoles, otu_rolesToFids = readFilteredOtuRoles(params.folder)
    except IOError:
        if params.verbose:
            sys.stderr.write("failed...generating file...")
        otudict = getOtuGenomeDictionary(MINN, COUNT)
        otu_fidsToRoles, otu_rolesToFids, missing_roles = filterFidsByOtusBetter(all_fidsToRoles, all_rolesToFids, otudict)
        writeFilteredOtuRoles(otu_fidsToRoles, params.folder)
    sys.stderr.write("done\n")
    
    # Generate a FASTA file for the fids in fidsToRoles
    sys.stderr.write("Subsystem FASTA file...")
    try:
        if params.verbose:
            sys.stderr.write("reading from file...")
        readSubsystemFasta(params.folder)
    except IOError:
        if params.verbose:
            sys.stderr.write("failed...generating file...")
        fidsToSeqs = fidsToSequences(otu_fidsToRoles.keys())
        writeSubsystemFasta(fidsToSeqs, params.folder)
    sys.stderr.write("done\n")
    
    # Complexes --> Roles
    # Needed to go from annotation likelihoods
    # to reaction likelihoods
    # Note that it is easier to go in this direction 
    #    Because we need all the roles in a complex to get the probability of that complex.
    sys.stderr.write("Complexes to roles...")
    try:
        if params.verbose:
            sys.stderr.write("reading from file...")
        complexToRequiredRoles = readComplexRoles(params.folder)
    except IOError:
        if params.verbose:
            sys.stderr.write("failed...generating file...")
        complexToRequiredRoles, requiredRolesToComplexes = complexRoleLinks(MINN, COUNT)
        writeComplexRoles(complexToRequiredRoles, params.folder)
    sys.stderr.write("done\n")
    
    # reaction --> complex
    # Again it is easier to go in this direction since we'll be filtering multiple complexes down to a single reaction.    
    sys.stderr.write("Reactions to complexes...")
    try:
        if params.verbose:
            sys.stderr.write("reading from file...")
        rxnToComplexes = readReactionComplex(params.folder)
    except IOError:
        if params.verbose:
            sys.stderr.write("failed...generating file...")
        rxnToComplexes, complexesToReactions = reactionComplexLinks(MINN, COUNT)
        writeReactionComplex(rxnToComplexes, params.folder)
    sys.stderr.write("done\n")
    
    sys.stderr.write("Data gathering done...\n")
    return(1)
    
def safeRemove(fname, dirname):
    totalfname = os.path.join(dirname, fname)
    try:
        # Check for file existence
        fid = open(totalfname, "r")
        fid.close()
        os.remove(totalfname)
    # If there is still an OSError despite the file existing we want to raise that, it will probably
    # cause problems trying to write to the files anyway. but an IOError just means the file isn't there.
    except IOError:
        pass

