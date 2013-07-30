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
CLIENT_TESTS_PYTHON = $(wildcard client-tests/*.py)
SCRIPT_TESTS = $(wildcard script-tests/*.py)
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
			PYTHONPATH="$(TARGET)/lib:" TARGET=$(TARGET) KB_TEST_CONFIG=test.cfg python $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

test-client:
	for t in $(CLIENT_TESTS_PYTHON) ; do \
		if [ -f $$t ] ; then \
			PYTHONPATH="$(TARGET)/lib:" KB_TEST_CONFIG=test.cfg python $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

# DEPLOYMENT
deploy: deploy-dir deploy-client deploy-service

deploy-service: deploy-libs deploy-scripts deploy-service-files deploy-cfg
deploy-client: deploy-libs deploy-scripts deploy-docs

deploy-scripts: deploy-perlscripts deploy-pythonscripts

deploy-dir:
	if [ ! -d $(SERV_SERVICE_DIR) ] ; then mkdir -p $(SERV_SERVICE_DIR) ; fi
	if [ ! -d docs ] ; then mkdir -p docs ; fi
	if [ ! -d ${SERV_SERVICE_DIR}/webroot/ ]; then mkdir -p ${SERV_SERVICE_DIR}/webroot/ ; fi

deploy-service-files:
	tpage $(SERV_TPAGE_ARGS) service/process.tt > $(SERV_SERVICE_DIR)/process.$(SERV_SERVICE); \
	chmod +x $(SERV_SERVICE_DIR)/process.$(SERV_SERVICE); \
	tpage $(SERV_TPAGE_ARGS) service/start_service.tt > $(SERV_SERVICE_DIR)/start_service; \
	chmod +x $(SERV_SERVICE_DIR)/start_service; \
	tpage $(SERV_TPAGE_ARGS) service/stop_service.tt > $(SERV_SERVICE_DIR)/stop_service; \
	chmod +x $(SERV_SERVICE_DIR)/stop_service;

deploy-perlscripts:
	# These three are needed to make these variables appear in the wrapped script
	export KB_TOP=$(TARGET); \
	export KB_RUNTIME=$(DEPLOY_RUNTIME); \
	export KB_PERL_PATH=$(TARGET)/lib bash ; \
	for src in $(SRC_PERL) ; do \
		basefile=`basename $$src`; \
		base=`basename $$src .pl`; \
		cp $$src $(TARGET)/plbin ; \
		$(WRAP_PERL_SCRIPT) "$(TARGET)/plbin/$$basefile" $(TARGET)/bin/$$base ; \
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
		$(WRAP_PYTHON_SCRIPT) "$(TARGET)/pybin/$$basefile" $(TARGET)/bin/$$base ; \
	done

deploy-docs:
	# pod2html doesn't work on the Python client (whcih has no docstrings anyway) but it does work on the Perl client.
	# Note - we can also run pydoc -w on our Impl.py file to get some documentation for those functions
	# but the formatting isn't consistent with what the other services use...
	$(KB_RUNTIME)/bin/pod2html -t "${SERVICE_NAME}" lib/Bio/KBase/${SERVICE_NAME}/Client.pm > docs/${SERVICE_NAME}.html
	cp docs/*html $(SERV_SERVICE_DIR)/webroot/

compile-typespec:
	mkdir -p lib/biokbase/${SERVICE_NAME}
	touch lib/biokbase/__init__.py
	touch lib/biokbase/${SERVICE_NAME}/__init__.py
	mkdir -p lib/javascript/${SERVICE_NAME}
	# The PSGI file that gets created is in Perl and it creates incorrect "use" statements
	# for the impl file and the service if we dont' specify these.
	compile_typespec \
	--pyimpl biokbase.${SERVICE_NAME}.Impl \
	--pyserver biokbase.${SERVICE_NAME}.Server \
	--psgi $(SERV_PSGI_PATH) \
	--client Bio::KBase::${SERVICE_NAME}::Client \
	--js javascript/${SERVICE_NAME}/Client \
	--py biokbase/${SERVICE_NAME}/Client \
	${SERV_SERVER_SPEC} lib
	rm -f lib/${SERVICE_NAME}Impl.py
	rm -f lib/${SERVICE_NAME}Server.py
	rm -rf Bio

include $(TOP_DIR)/tools/Makefile.common.rules
