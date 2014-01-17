#! /bin/bash

GENOMELIST="158879.1 171101.1 208964.1 224308.1 243273.1 243277.1 272635.1 401614.5 62977.3 71421.1 83333.1 85962.1 99287.1"
#GENOMELIST="171101.1"

SCRIPT_DIR=~/kb/dev_container/modules/probabilistic_annotation/internalScripts/analysis-scripts
WORKSPACE=RecalcultedSeedModelGapfills_FINAL_January_2013

for GENOME in $GENOMELIST; do
    MODEL1=${GENOME}.model.std.int.minimal.int
    KOSIM1=${GENOME}.model.std.int.minimal.int.knockoutsim
    MODEL2=${GENOME}.model.pa.int.minimal.int
    KOSIM2=${GENOME}.model.pa.int.minimal.int.knockoutsim
    MODEL3=${GENOME}.model.std.iterative.int.minimal.int
    KOSIM3=${GENOME}.model.std.iterative.int.minimal.int.knockoutsim
    MODEL4=${GENOME}.model.pa.iterative.int.minimal.int
    KOSIM4=${GENOME}.model.pa.iterative.int.minimal.int.knockoutsim

    ${SCRIPT_DIR}/GeneratePhenotypeROCCurve.py --script-dir ${SCRIPT_DIR} --fba-url http://localhost:7036 ${GENOME} ${MODEL1} ${KOSIM1} ${WORKSPACE}
    ${SCRIPT_DIR}/GeneratePhenotypeROCCurve.py --script-dir ${SCRIPT_DIR} --fba-url http://localhost:7036 ${GENOME} ${MODEL2} ${KOSIM2} ${WORKSPACE} 
    ${SCRIPT_DIR}/GeneratePhenotypeROCCurve.py --script-dir ${SCRIPT_DIR} --fba-url http://localhost:7036 ${GENOME} ${MODEL3} ${KOSIM3} ${WORKSPACE}
    ${SCRIPT_DIR}/GeneratePhenotypeROCCurve.py --script-dir ${SCRIPT_DIR} --fba-url http://localhost:7036 ${GENOME} ${MODEL4} ${KOSIM4} ${WORKSPACE} 
done
