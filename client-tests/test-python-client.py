import sys
import unittest
import subprocess
import traceback
import time
import os
import json
from biokbase.workspaceService.Client import workspaceService
from biokbase.workspaceService.Client import ServerError as WorkspaceServerError
from biokbase.fbaModelServices.Client import fbaModelServices
from biokbase.probabilistic_annotation.Client import ProbabilisticAnnotation
from biokbase.probabilistic_annotation.DataParser import getConfig

class TestPythonClient(unittest.TestCase):    

    def setUp(self):
        # Set configuration variables.
        self._config = getConfig(os.environ["KB_TEST_CONFIG"])

        # Get an authorization token for the test user.
        ws = workspaceService(self._config["workspace_url"], user_id=self._config["test_user"], password=self._config["test_pwd"])
        self._token = ws._headers['AUTHORIZATION']

        # In the test environment (as of 08-09-2013) the RxnProbs object isn't defined in the workspace so it isn't allowed to save the results.
        # This is a patch until we fix that more permanently.
        try:
            ws.add_type( { "type" : "RxnProbs",
                           "auth" : self._token 
                           } )
        except:
            # Type already exists
            pass

    def test_loadGenome(self):
        ''' Load a test Genome object into the test workspace. '''
        
        # Create the test workspace.
        wsClient = workspaceService(self._config["workspace_url"])
        try:
            # When the workspace exists, delete it so there is a clean slate for the test.
            wsMetadata = wsClient.get_workspacemeta( { "workspace": self._config["wsid"], "auth": self._token } )
            wsClient.delete_workspace( { "workspace": self._config["wsid"], "auth": self._token } )
        except WorkspaceServerError as e:
            # Hopefully this means the workspace does not exist. (It could also mean someone messed up setting up the URLs)
#            traceback.print_exc(file=sys.stderr)
            pass
        
        wsMetadata = wsClient.create_workspace( { "workspace": self._config["wsid"], "default_permission": "n", "auth": self._token } )

        # We also need to put in a mapping and a biochemistry object somewhere.
        # To do this, I just create a "depednecy workspace" and pull them from there.
        self._config["dependency_ws"]
        try:
            # When the workspace exists, delete it so there is a clean slate for the test.
            wsMetadata = wsClient.get_workspacemeta( { "workspace": self._config["dependency_ws"], "auth": self._token } )
            wsClient.delete_workspace( { "workspace": self._config["dependency_ws"], "auth": self._token } )
        except WorkspaceServerError as e:
            # Hopefully this means the workspace does not exist. (It could also mean someone messed up setting up the URLs)
#            traceback.print_exc(file=sys.stderr)
            pass

        depWsMetadata = wsClient.create_workspace( { "workspace": self._config["dependency_ws"], "default_permission": "n", "auth": self._token } )

        # Load the mapping and biochemsitry objects
        # (Note - this won't be allowed with the new version of the workspace but it agrees with what is in there now)
        mappingObject = json.load(open(self._config["test_mapping_object"], "r"))
        wsClient.save_object( {
                "id" : "default",
                "type" : "Mapping",
                "data" : mappingObject,
                "auth" : self._token,
                "workspace" : self._config["dependency_ws"]
               } )

        biochemistryObject = json.load(open(self._config["test_biochemistry_object"], "r"))
        wsClient.save_object( {
                "id" : "default",
                "type" : "Biochemistry",
                "data" : biochemistryObject,
                "auth" : self._token,
                "workspace" : self._config["dependency_ws"]
               } )


        # Load the test Genome object into the workspace.
        fbaClient = fbaModelServices(self._config["fbamodelservices_url"])
        testGenome = json.load(open(self._config["genome_file"], "r"))
        genomeMetadata = fbaClient.genome_object_to_workspace( { 
            "genomeobj": testGenome,
            "workspace": self._config["wsid"],
            "auth": self._token,
            "mapping_workspace": self._config["dependency_ws"]
            } )
        
    def test_annotate(self):
        ''' Run pa-annotate on a valid Genome object and verify that the job runs and returns a valid ProbAnno object in the expected time.'''

        # Run the annotate() function to generate a ProbAnno object.
        paClient = ProbabilisticAnnotation(self._config["probanno_url"], token=self._token)
        jobid = paClient.annotate( {
            "genome": self._config["genomeid"],
            "genome_workspace": self._config["wsid"],
            "probanno": self._config["probannoid"],
            "probanno_workspace": self._config["wsid"] } )
        
        # Allow time for the command to run.
        time.sleep(float(self._config["runtime"]))
        
        # Look for the ProbAnno object in the test workspace.
        wsClient = workspaceService(self._config["workspace_url"])
        try:
            output = wsClient.get_object( {
                "workspace": self._config["wsid"],
                "type": "ProbAnno",
                "id": self._config["probannoid"],
                "auth": self._token } )
            self.assertEqual(output["metadata"][4], "pa-annotate")
            # TODO Could add some checking of the object data here
        except WorkspaceServerError as e:
            self.fail(msg = "The expected object %s did not get created in the workspace %s!\n" %(self._config["probannoid"], self._config["wsid"]))
        
    def test_calculate(self):
        ''' Run pa-calculate on a valid ProbAnno object and verify that the job runs and returns a valid RxnProbs object.'''
        
        # Run the calculate() function to generate a RxnProbs object.
        paClient = ProbabilisticAnnotation(self._config["probanno_url"], token=self._token)
        rxnprobsMetadata = paClient.calculate( {
            "probanno":           self._config["probannoid"],
            "probanno_workspace": self._config["wsid"],
            "rxnprobs":           self._config["rxnprobsid"],
            "rxnprobs_workspace": self._config["wsid"] 
            } )
         
        # Look for the RxnProbs object in the test workspace.
        wsClient = workspaceService(self._config["workspace_url"])
        try:
            output = wsClient.get_object( {
                "workspace": rxnprobsMetadata[7],
                "type": rxnprobsMetadata[1],
                "id": rxnprobsMetadata[0],
                "auth": self._token } )
            self.assertEqual(output["metadata"][4], "pa-calculate")
            # TODO Could add some checking of the object data here
        except WorkspaceServerError as e:
            self.fail(msg = "The expected object %s did not get created in the workspace %s!\n" %(self._config["rxnprobsid"], self._config["wsid"]))

    def test_get_rxnprobs(self):
        ''' Verify that we can successfully get a list of rxnprobs data from a valid RxnProbs object.'''
        paClient = ProbabilisticAnnotation(self._config["probanno_url"], token=self._token)
        rxnProbsData = paClient.get_rxnprobs( {
                "rxnprobs":           self._config["rxnprobsid"],
                "rxnprobs_workspace": self._config["wsid"]
                })

    def test_cleanup(self):
        ''' Cleanup objects created by tests. '''
        
        # Delete the test workspace.
        wsClient = workspaceService(self._config["workspace_url"])
        wsClient.delete_workspace( { "workspace": self._config["wsid"], "auth": self._token } )
        wsClient.delete_workspace( { "workspace": self._config["dependency_ws"], "auth": self._token } )
        
if __name__ == '__main__':
    # Create a suite, add tests to the suite, run the suite.
    suite = unittest.TestSuite()
    suite.addTest(TestPythonClient('test_loadGenome'))
    suite.addTest(TestPythonClient('test_annotate'))
    suite.addTest(TestPythonClient('test_calculate'))
    suite.addTest(TestPythonClient('test_get_rxnprobs'))
    suite.addTest(TestPythonClient('test_cleanup'))
    unittest.TextTestRunner().run(suite)
    
