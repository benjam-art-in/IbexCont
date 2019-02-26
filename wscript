#! /usr/bin/env python
# encoding: utf-8

from waflib import Logs, ConfigSet

######################
###### options #######
######################
def options (opt):
	opt.add_option ("--with-cont1", action="store_true", dest="WITH_CONT1",
			help = "install IbexSolve plugin")

######################
##### configure ######
######################
def configure (conf):
	conf.env.WITH_CONT1 = conf.options.WITH_CONT1

	conf.start_msg ("plugin One-Parameter Continuation (cont1)")
	if not conf.env.WITH_CONT1:
		conf.end_msg ("not used")
		return

	conf.end_msg ("enabled")

	conf.env.append_unique ("IBEX_PLUGIN_USE_LIST", "CONT1")

	# check that conf.env.LP_LIB is not set to NONE (i.e., we need a LP library
	# for the cont1 plugin).
	#if conf.env.LP_LIB == "NONE":
	#	conf.fatal ("cont1 plugin needs a LP library")

	# We need -std=c++11 to compile ibexcont1
	#conf.check_cxx(cxxflags = "-std=c++11", uselib_store="IBEXCONT1")

	# To fix Windows compilation problem (strdup with std=c++11, see issue #287)
	#conf.check_cxx(cxxflags = "-U__STRICT_ANSI__", uselib_store="IBEXCONT1")

	# Add information in ibex_Setting
	conf.setting_define ("WITH_CONT1", 1)

	# add CONT1 plugin include directory
	for f in conf.path.ant_glob ("src/** src", dir = True, src = False):
		conf.env.append_unique("INCLUDES_CONT1", f.abspath())

	# The build and install steps will be run from the main src/wscript script so
	# we need to give path relative to the main src directory
	mainsrc = conf.srcnode.make_node ("src")

	# add CONT1 headers
	for f in conf.path.ant_glob ("src/**/ibex_*.h"):
		conf.env.append_unique ("IBEX_HDR", f.path_from (mainsrc))

	# add CONT1 source files
	for f in conf.path.ant_glob ("src/**/ibex_*.cpp"):
		conf.env.append_unique ("IBEX_SRC", f.path_from (mainsrc))


	# ------- Tests in this plugin are not using CPPUNIT so far ------------

	# The utest step will be run from the main tests/wscript script so we need to
	# give path relative to the main tests directory
	#maintests = conf.srcnode.make_node ("tests")

	# add CONT1 test files
	#for f in conf.path.ant_glob ("tests/**/Test*.cpp"):
	#	conf.env.append_unique ('TEST_SRC', f.path_from (maintests))

	# Add cont1/tests directory to list of INCLUDES for TESTS
	#testsnode = conf.path.make_node ("tests")
	#conf.env.append_unique ("INCLUDES_TESTS", testsnode.abspath ())

	# ----------------------------------------------------------------------
######################
####### build ########
######################
def build (bld):
	
	if bld.env.WITH_CONT1:
		# build cont1 binary
		bld.program (
		target = "ibexcont1",
		use = [ "ibex", "IBEXCONT1" ], # add dependency on ibex library
		source = bld.path.ant_glob ("main/**/*.cpp"),
		install_path = bld.env.BINDIR,
		)
