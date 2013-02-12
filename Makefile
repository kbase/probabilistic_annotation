# Set default locations for runtime and deployment
# if the directories are not already set:
DEPLOY_RUNTIME ?= /kb/runtime
TARGET         ?= /kb/deployment

# Include standard makefile
TOP_DIR = ../..
#include $(TOP_DIR)/tools/Makefile.common
SRC_PYTHON = $(wildcard scripts/*.py)

SERVICE_NAME=probabilistic_annotation
SERV_SERVER_SPEC 	= ProbabilisticAnnotation.spec
SERV_SERVER_MODULE 	= ${SERVICE_NAME}
SERV_SERVICE 		= ${SERVICE_NAME}
SERV_PSGI_PATH 		= lib/${SERVICE_NAME}.psgi
SERV_SERVICE_PORT 	= 7073
SERV_SERVICE_DIR = $(TARGET)/services/$(SERV_SERVICE)
SERV_TPAGE = $(KB_RUNTIME)/bin/perl $(KB_RUNTIME)/bin/tpage
SERV_TPAGE_ARGS = --define kb_top=$(TARGET) --define kb_runtime=$(KB_RUNTIME) --define kb_service_name=$(SERV_SERVICE) \
	--define kb_service_port=$(SERV_SERVICE_PORT) --define kb_service_psgi=$(SERV_PSGI_PATH)

all: bin compile-typespec


$(BIN_DIR)/%: scripts/%.pl 
	$(TOOLS_DIR)/wrap_perl '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

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

deploy: deploy-client deploy-service
deploy-all: deploy-client deploy-service

deploy-service: deploy-dir deploy-libs deploy-scripts deploy-services
deploy-client: deploy-dir deploy-libs deploy-scripts deploy-docs

deploy-dir:
	if [ ! -d $(SERV_SERVICE_DIR) ] ; then mkdir -p $(SERV_SERVICE_DIR) ; fi
	if [ ! -d $(SERV_SERVICE_DIR)/webroot ] ; then mkdir -p $(SERV_SERVICE_DIR)/webroot ; fi

deploy-scripts:
	export KB_TOP=$(TARGET); \
	export KB_RUNTIME=$(KB_RUNTIME); \
	export KB_PYTHON_PATH=$(TARGET)/lib bash ; \
	for src in $(SRC_PYTHON) ; do \
		# ??
	done 

deploy-libs: compile-typespec
	rsync -arv lib/. $(TARGET)/lib/.

deploy-services: deploy-basic-service


deploy-docs:
	if [ ! -d docs ] ; then mkdir -p docs ; fi
	# The python code doesn't use pod but would need something different to make the nice HTML... I'll deal with this later.
	#	$(KB_RUNTIME)/bin/pod2html -t "workspaceService" lib/Bio/KBase/workspaceService/Client.pm > docs/workspaceService.html
	#	cp docs/*html $(SERV_SERVICE_DIR)/webroot/.

compile-typespec:
	mkdir -p lib/biokbase/${SERVICE_NAME}
	touch lib/biokbase/__init__.py
	touch lib/biokbase/${SERVICE_NAME}/__init__.py
	mkdir -p lib/javascript/${SERVICE_NAME}
	compile_typespec \
	-impl Bio::KBase::${SERVICE_NAME}::Impl \
	-service Bio::KBase::${SERVICE_NAME}::Server \
	-psgi lib/${SERVICE_NAME}/psgi \
	-client Bio::KBase::${SERVICE_NAME}::Client \
	-js javascript/${SERVICE_NAME}/Client \
	-py biokbase/${SERVICE_NAME}/client \
	${SERV_SERVER_SPEC} lib
	rm -f lib/${SERVICE_NAME}Impl.py
	rm -f lib/${SERVICE_NAME}Server.py
	rm -rf Bio
