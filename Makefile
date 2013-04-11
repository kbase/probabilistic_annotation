# Set default locations for runtime and deployment
# if the directories are not already set:
DEPLOY_RUNTIME ?= /kb/runtime
TARGET         ?= /kb/deployment
TOOLS_DIR      ?= /kb/dev_container/tools

# Include standard makefile
TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

INTERNAL_PYTHON = $(wildcard internalScripts/*.py)

SERVICE_NAME=probabilistic_annotation
SERV_SERVER_SPEC 	= ProbabilisticAnnotation.spec
SERV_SERVER_MODULE 	= ${SERVICE_NAME}
SERV_SERVICE 		= ${SERVICE_NAME}
# the lib/ prefix is added automatically to this.
SERV_PSGI_PATH 		= ${SERVICE_NAME}.psgi
SERV_SERVICE_PORT 	= 7073
SERV_SERVICE_DIR = $(TARGET)/services/$(SERV_SERVICE)
SERV_TPAGE = $(KB_RUNTIME)/bin/perl $(KB_RUNTIME)/bin/tpage
# but the lib/ is NOT added automatically here.
SERV_TPAGE_ARGS = --define kb_top=$(TARGET) --define kb_runtime=$(KB_RUNTIME) --define kb_service_name=$(SERV_SERVICE) \
	--define kb_service_port=$(SERV_SERVICE_PORT) --define kb_service_psgi=lib/$(SERV_PSGI_PATH)

all: compile-typespec

# TESTS
# Note I don't have any test scripts yet but when I make them they'll go in these locations
CLIENT_TESTS = $(wildcard client-tests/*.t)
SCRIPT_TESTS = $(wildcard script-tests/*.sh)
SERVER_TESTS = $(wildcard server-tests/*.t)

test: test-service test-client test-scripts
	@echo "running server, script and client tests"

test-service:
	for t in $(SERVER_TESTS) ; do \
		if [ -f $$t ] ; then \
			$(DEPLOY_RUNTIME)/bin/prove $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

test-scripts:
	for t in $(SCRIPT_TESTS) ; do \
		if [ -f $$t ] ; then \
			/bin/sh $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

test-client:
	for t in $(CLIENT_TESTS) ; do \
		if [ -f $$t ] ; then \
			$(DEPLOY_RUNTIME)/bin/prove $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

# DEPLOYMENT
deploy: deploy-dir deploy-client deploy-service

deploy-service: deploy-libs deploy-scripts deploy-service-files deploy-config
deploy-client: deploy-libs deploy-scripts deploy-docs deploy-config

deploy-scripts: deploy-perlscripts deploy-pythonscripts

deploy-config:
	# We need to create a config file if it doesn't exist because auth requires it.
	# This will evantually become more standardized but for now I want it to "just work".
	if [ ! -e $(TARGET)/deployment.cfg ]; then touch $(TARGET)/deployment.cfg; fi

deploy-dir:
	if [ ! -d $(SERV_SERVICE_DIR) ] ; then mkdir -p $(SERV_SERVICE_DIR) ; fi
#        if [ ! -d $(SERV_SERVICE_DIR)/webroot ] ; then mkdir -p $(SERV_SERVICE_DIR)/webroot ; fi

deploy-service-files:
	tpage $(SERV_TPAGE_ARGS) service/start_service.tt > $(SERV_SERVICE_DIR)/start_service; \
	chmod +x $(SERV_SERVICE_DIR)/start_service; \
        tpage $(SERV_TPAGE_ARGS) service/stop_service.tt > $(SERV_SERVICE_DIR)/stop_service; \
        chmod +x $(SERV_SERVICE_DIR)/stop_service; \
        tpage $(SERV_TPAGE_ARGS) service/process.tt > $(SERV_SERVICE_DIR)/process.$(SERV_SERVICE); \
        chmod +x $(SERV_SERVICE_DIR)/process.$(SERV_SERVICE);

deploy-perlscripts:
	# These three are needed to make these variables appear in the wrapped script
	export KB_TOP=$(TARGET); \
	export KB_RUNTIME=$(DEPLOY_RUNTIME); \
	export KB_PERL_PATH=$(TARGET)/lib bash ; \
	for src in $(SRC_PERL) ; do \
		basefile=`basename $$src`; \
		base=`basename $$src .pl`; \
		cp $$src $(TARGET)/plbin ; \
		bash $(TOOLS_DIR)/wrap_perl.sh "$(TARGET)/plbin/$$basefile" $(TARGET)/bin/$$base ; \
	done

deploy-pythonscripts:
	# These three are needed to make these variables appear in the wrapped script
	export KB_TOP=$(TARGET); \
	export KB_RUNTIME=$(DEPLOY_RUNTIME); \
	export KB_PYTHON_PATH=$(TARGET)/lib bash ; \
	for src in $(INTERNAL_PYTHON) $(SRC_PYTHON) ; do \
		basefile=`basename $$src`; \
		base=`basename $$src .py`; \
		cp $$src $(TARGET)/pybin ; \
		bash $(TOOLS_DIR)/wrap_python.sh "$(TARGET)/pybin/$$basefile" $(TARGET)/bin/$$base ; \
	done

deploy-libs:
	# lib/ contains basically all of the guts of the code. 
	# Reference a python module within here with "import biokbase.probabilistic_annotation.[module_name]"
	rsync -arv lib/. $(TARGET)/lib/.

deploy-docs:
	# I have nothing here yet.
	# The python code doesn't use pod but would need something different to make the nice HTML... I'll deal with this later.
	if [ ! -d docs ] ; then mkdir -p docs ; fi

compile-typespec:
	mkdir -p lib/biokbase/${SERVICE_NAME}
	touch lib/biokbase/__init__.py
	touch lib/biokbase/${SERVICE_NAME}/__init__.py
	mkdir -p lib/javascript/${SERVICE_NAME}
	compile_typespec \
	-impl Bio::KBase::${SERVICE_NAME}::Impl \
	-service Bio::KBase::${SERVICE_NAME}::Server \
	-psgi $(SERV_PSGI_PATH) \
	-client Bio::KBase::${SERVICE_NAME}::Client \
	-js javascript/${SERVICE_NAME}/Client \
	-py biokbase/${SERVICE_NAME}/client \
	${SERV_SERVER_SPEC} lib
	rm -f lib/${SERVICE_NAME}Impl.py
	rm -f lib/${SERVICE_NAME}Server.py
	rm -rf Bio
