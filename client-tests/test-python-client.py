import sys
import unittest
import subprocess
import traceback
import time
import os
import json
from biokbase.workspace.client import Workspace
from biokbase.workspace.client import ServerError as WorkspaceServerError
from biokbase.fbaModelServices.Client import fbaModelServices
from biokbase.userandjobstate.client import UserAndJobState, ServerError as JobStateServerError
from biokbase.probabilistic_annotation.Client import ProbabilisticAnnotation
from biokbase.probabilistic_annotation.DataParser import getConfig

class TestPythonClient(unittest.TestCase):    

    def setUp(self):
        # Set configuration variables.
        self._config = getConfig(os.environ["KB_TEST_CONFIG"])

        # Get an authorization token for the test user.
        wsClient = Workspace(self._config["workspace_url"], user_id=self._config["test_user"], password=self._config["test_pwd"])
        self._token = wsClient._headers['AUTHORIZATION']

    def test_loadGenome(self):
        ''' Load a test Genome object into the test workspace. '''
        
        # Create the test workspace.
        wsClient = Workspace(self._config["workspace_url"], token=self._token)
        try:
            # When the workspace exists, make sure all of the objects are deleted so there is a clean slate for the test.
            objectList = wsClient.list_objects( { 'workspaces': [ self._config['test_ws'] ] } )
            deleteList = list()
            for object in objectList:
                deleteList.append( { 'wsid': object[6], 'objid': object[0], 'ver': object[4] })
            if len(deleteList) > 0:
                wsClient.delete_objects(deleteList)

        except WorkspaceServerError as e:
            # Hopefully this means the workspace does not exist. (It could also mean someone messed up setting up the URLs)
            traceback.print_exc(file=sys.stderr)
            wsInfo = wsClient.create_workspace( { "workspace": self._config["test_ws"] } )

        # We also need to put in a mapping and a biochemistry object somewhere.
        # To do this, I just create a "dependency workspace" and pull them from there.
        self._config["dependency_ws"]
        try:
            # When the workspace exists, delete it so there is a clean slate for the test.
            wsInfo = wsClient.get_workspace_info( { "workspace": self._config["dependency_ws"] } )
        except WorkspaceServerError as e:
            # Hopefully this means the workspace does not exist. (It could also mean someone messed up setting up the URLs)
#            traceback.print_exc(file=sys.stderr)
            depWsInfo = wsClient.create_workspace( { "workspace": self._config["dependency_ws"] } )

        # Load the mapping and biochemsitry objects
        # (Note - this won't be allowed with the new version of the workspace but it agrees with what is in there now)
        #mappingObject = json.load(open(self._config["test_mapping_object"], "r"))
        #mappingObjectData = { 'type': 'Mapping', 'name': 'default', 'data': mappingObject }

#         biochemistryObject = json.load(open(self._config["test_biochemistry_object"], "r"))
#         biochemistryObjectData = { 'type': 'Biochemistry', 'name': 'default', 'data': biochemistryObject }
#         
#         objectList = list()
        #objectList.append(mappingObjectData)
#         objectList.append(biochemistryObjectData)
#         wsClient.save_objects( { 'workspace': self._config["dependency_ws"], 'objects': objectList } )

        # Load the test Genome object into the workspace.
#         fbaClient = fbaModelServices(self._config["fbamodelservices_url"])
#         testGenome = json.load(open(self._config["genome_file"], "r"))
#         genomeMetadata = fbaClient.genome_object_to_workspace( { 
#             "genomeobj": testGenome,
#             "workspace": self._config["test_ws"],
#             "auth": self._token
#             "mapping_workspace": self._config["dependency_ws"]
#             } )
        testContigSet = json.load(open(self._config['contigset_file'], 'r'))
        contigSetSaveData = dict()
        contigSetSaveData['type'] = 'KBaseGenomes.ContigSet'
        contigSetSaveData['name'] = self._config['contigsetid']
        contigSetSaveData['data'] = testContigSet        
        wsClient.save_objects( { 'workspace': self._config['test_ws'], 'objects': [ contigSetSaveData ] } )
        testGenome = json.load(open(self._config["genome_file"], "r"))
        genomeSaveData = dict()
        genomeSaveData['type'] = 'KBaseGenomes.Genome'
        genomeSaveData['name'] = self._config['genomeid']
        genomeSaveData['data'] = testGenome
        wsClient.save_objects( { 'workspace': self._config['test_ws'], 'objects': [ genomeSaveData ] } )
        
    def test_annotate(self):
        ''' Run pa-annotate on a valid Genome object and verify that the job runs and returns a valid ProbAnno object in the expected time.'''

        # Run the annotate() function to generate a ProbAnno object.
        paClient = ProbabilisticAnnotation(self._config["probanno_url"], token=self._token)
        jobid = paClient.annotate( {
            "genome": self._config["genomeid"],
            "genome_workspace": self._config["test_ws"],
            "probanno": self._config["probannoid"],
            "probanno_workspace": self._config["test_ws"] } )
        
        # Allow time for the command to run.
        time.sleep(float(self._config["runtime"]))
        
        # Make sure the job has completed.
        ujsClient = UserAndJobState(self._config['ujs_url'], token=self._token)
        jobList = ujsClient.list_jobs([ self._config['test_user'] ], 'C')
        jobCompleted = False
        for job in jobList:
            if jobid == job[0]:
                jobCompleted = True
        self.assertEqual(jobCompleted, True)
        
        # Look for the ProbAnno object in the test workspace.
        wsClient = Workspace(self._config["workspace_url"], token=self._token)
        try:
            probannoObjectId = { 'workspace': self._config['test_ws'], 'name': self._config['probannoid'] }
            objectList = wsClient.get_objects( [ probannoObjectId ] )
            probannoObject = objectList[0]
            self.assertEqual(probannoObject['info'][1], self._config['probannoid'])
            # TODO Could add some checking of the object data here
        except WorkspaceServerError as e:
            traceback.print_exc(file=sys.stderr)
            self.fail(msg = "The expected object %s did not get created in the workspace %s!\n" %(self._config["probannoid"], self._config["test_ws"]))
        
    def test_calculate(self):
        ''' Run pa-calculate on a valid ProbAnno object and verify that the job runs and returns a valid RxnProbs object.'''
        
        # Run the calculate() function to generate a RxnProbs object.
        paClient = ProbabilisticAnnotation(self._config["probanno_url"], token=self._token)
        rxnprobsMetadata = paClient.calculate( {
            "probanno":           self._config["probannoid"],
            "probanno_workspace": self._config["test_ws"],
            "rxnprobs":           self._config["rxnprobsid"],
            "rxnprobs_workspace": self._config["test_ws"] 
            } )
         
        # Look for the RxnProbs object in the test workspace.
        wsClient = Workspace(self._config["workspace_url"], token=self._token)
        try:
            rxnprobsObjectId = { 'workspace': self._config['test_ws'], 'name': self._config['rxnprobsid'] }
            objectList = wsClient.get_objects( [ rxnprobsObjectId ] )
            rxnprobsObject = objectList[0]
            self.assertEqual(rxnprobsObject['info'][1], self._config['rxnprobsid'])
            # TODO Could add some checking of the object data here
        except WorkspaceServerError as e:
            traceback.print_exc(file=sys.stderr)
            self.fail(msg = "The expected object %s did not get created in the workspace %s!\n" %(self._config["rxnprobsid"], self._config["test_ws"]))

    def test_get_rxnprobs(self):
        ''' Verify that we can successfully get a list of rxnprobs data from a valid RxnProbs object.'''
        paClient = ProbabilisticAnnotation(self._config["probanno_url"], token=self._token)
        rxnProbsData = paClient.get_rxnprobs( {
                "rxnprobs":           self._config["rxnprobsid"],
                "rxnprobs_workspace": self._config["test_ws"]
                })

    def test_get_probanno(self):
        return
    
    def test_cleanup(self):
        ''' Cleanup objects created by tests. '''
        
        # Delete all of the objects in the test workspace.
        wsClient = Workspace(self._config["workspace_url"], token=self._token)
        listObjectsParams = dict()
        listObjectsParams['workspaces'] = [ self._config['test_ws'] ]
        objectList = wsClient.list_objects(listObjectsParams)
        deleteList = list()
        for object in objectList:
            deleteList.append( { 'wsid': object[6], 'objid': object[0], 'ver': object[4] })
        wsClient.delete_objects(deleteList)
        objectList = wsClient.list_objects(listObjectsParams)
        for object in objectList:
            print 'After delete Object %s %s' %(object[1], object[4])
        
        
if __name__ == '__main__':
    # Create a suite, add tests to the suite, run the suite.
    suite = unittest.TestSuite()
    suite.addTest(TestPythonClient('test_loadGenome'))
    suite.addTest(TestPythonClient('test_annotate'))
    suite.addTest(TestPythonClient('test_calculate'))
    suite.addTest(TestPythonClient('test_get_rxnprobs'))
    suite.addTest(TestPythonClient('test_cleanup'))
    unittest.TextTestRunner().run(suite)
    
