#!/usr/bin/env python

from __future__ import (absolute_import, division, print_function)
from multiprocessing import Process, Array, Manager
import optparse
import os
import sys

#########################################################################
# Set up the command line options
#########################################################################

VERSION = "1.1"
THIS_MODULE_DIR = os.path.dirname(os.path.realpath(__file__))
DEFAULT_FRAMEWORK_LOC = os.path.realpath(os.path.join(THIS_MODULE_DIR, "..","lib","systemtests"))
DATA_DIRS_LIST_PATH = os.path.join(THIS_MODULE_DIR, "datasearch-directories.txt")
SAVE_DIR_LIST_PATH = os.path.join(THIS_MODULE_DIR, "defaultsave-directory.txt")

info = []
info.append("This program will configure mantid run all of the system tests located in")
info.append("the 'tests/analysis' directory.")
info.append("This program will create a temporary 'Mantid.user.properties' file which")
info.append("it will rename to 'Mantid.user.properties.systest' upon completion. The")
info.append("current version of the code does not print to stdout while the test is")
info.append("running, so the impatient user may ^C to kill the process. In this case")
info.append("all of the tests that haven't been run will be marked as skipped in the")
info.append("full logs.")

parser = optparse.OptionParser("Usage: %prog [options]", None,
                               optparse.Option, VERSION, 'error', ' '.join(info))
parser.add_option("", "--email", action="store_true",
                  help="send an email with test status.")
parser.add_option("-x", "--executable", dest="executable",
                  help="The executable path used to run each test. Default is the sys.executable")
parser.add_option("-a", "--exec-args", dest="execargs",
                  help="Arguments passed to executable for each test Default=[]")
parser.add_option("", "--frameworkLoc",
                  help="location of the stress test framework (default=%s)" % DEFAULT_FRAMEWORK_LOC)
parser.add_option("", "--disablepropmake", action="store_false", dest="makeprop",
                  help="By default this will move your properties file out of the "
                  + "way and create a new one. This option turns off this behavior.")
parser.add_option("-R", "--tests-regex", dest="testsInclude",
                  help="String specifying which tests to run. Simply uses 'string in testname'.")
parser.add_option("-E", "--excluderegex", dest="testsExclude",
                  help="String specifying which tests to not run. Simply uses 'string in testname'.")
loglevelChoices=["error", "warning", "notice", "information", "debug"]
parser.add_option("-l", "--loglevel", dest="loglevel",
                  choices=loglevelChoices,
                  help="Set the log level for test running: [" + ', '.join(loglevelChoices) + "]")
parser.add_option("-j", "--parallel", dest="ncores", action="store", type="int",
                  help="The number of instances to run in parallel, like the -j option in ctest. Default is 1.")
parser.add_option("-q", "--quiet", dest="quiet", action="store_true",
                  help="Prints detailed log to terminal.")
parser.add_option("-c", "--clean", dest="clean", action="store_true",
                  help="Performs a cleanup of the data generated by the test suite (does not run the tests).")
parser.add_option("", "--output-on-failure", dest="output_on_failure", action="store_true",
                  help="Print full log for failed tests.")
parser.add_option("", "--showskipped", dest="showskipped", action="store_true",
                  help="List the skipped tests.")
parser.add_option("-d", "--datapaths", dest="datapaths",
                  help="A semicolon-separated list of directories to search for data")
parser.add_option("-s", "--savedir", dest="savedir",
                  help="A directory to use for the Mantid save path")
parser.add_option("", "--archivesearch", dest="archivesearch", action="store_true",
                  help="Turn on archive search for file finder.")
parser.add_option("", "--exclude-in-pull-requests", dest="exclude_in_pr_builds",action="store_true",
                  help="Skip tests that are not run in pull request builds")
parser.set_defaults(frameworkLoc=DEFAULT_FRAMEWORK_LOC, executable=sys.executable, makeprop=True,
                    loglevel="information", ncores=1, quiet=False, output_on_failure=False, clean=False)
(options, args) = parser.parse_args()

# import the stress testing framework
sys.path.append(options.frameworkLoc)
import stresstesting

#########################################################################
# Configure mantid
#########################################################################

# Parse files containing the search and save directories, unless otherwise given
data_paths = options.datapaths
if data_paths is None or data_paths == "":
    with open(DATA_DIRS_LIST_PATH, 'r') as f_handle:
        data_paths = f_handle.read().strip()

save_dir = options.savedir
if save_dir is None or save_dir == "":
    with open(SAVE_DIR_LIST_PATH, 'r') as f_handle:
        save_dir = f_handle.read().strip()

# Configure properties file
mtdconf = stresstesting.MantidFrameworkConfig(loglevel=options.loglevel,
                                              data_dirs=data_paths, save_dir=save_dir,
                                              archivesearch=options.archivesearch)
if options.makeprop:
    mtdconf.config()

#########################################################################
# Run the tests
#########################################################################

# Multi-core processes
processes = [] # an array to hold the processes
results_array = Array('i', 4*options.ncores) # shared array to hold skipped, failed and total number of tests + status
test_counter = Array('i', [0, 0, 0]) # shared array for number of executed tests, total number and max name length
manager = Manager() # a manager to create a shared dict to store names of skipped and failed tests
status_dict = manager.dict() # a shared dict to store names of skipped and failed tests

# Prepare ncores processes
for ip in range(options.ncores):
    processes.append(Process(target=stresstesting.testProcess,args=(mtdconf.testDir, mtdconf.saveDir,
                     options, results_array, status_dict, test_counter, ip)))
# Start and join processes
for p in processes:
    p.start()
for p in processes:
    p.join()


# Gather results
skippedTests = sum(results_array[:options.ncores])
failedTests = sum(results_array[options.ncores:2*options.ncores])
totalTests = sum(results_array[2*options.ncores:3*options.ncores])
# Find minimum of status: if min == 0, then success is False
success = bool(min(results_array[3*options.ncores:4*options.ncores]))

#########################################################################
# Cleanup
#########################################################################

# Put the configuration back to its original state
if options.makeprop:
    mtdconf.restoreconfig()

#########################################################################
# Output summary to terminal (skip if this was a cleanup run)
#########################################################################

if (not options.clean):
    nwidth = 80
    banner = "#" * nwidth
    print()
    print(banner)
    print("#"+(" "*(int(nwidth/2)-4))+"SUMMARY")
    print(banner)

    if (skippedTests > 0) and options.showskipped:
        print("\nSKIPPED:")
        for key in status_dict.keys():
            if status_dict[key] == 'skipped':
                print(key)
    if failedTests > 0:
        print("\nFAILED:")
        for key in status_dict.keys():
            if status_dict[key] == 'failed':
                print(key)

    # Report global statistics on tests
    print()
    if skippedTests == totalTests:
        print("All tests were skipped")
        success = False # fail if everything was skipped
    else:
        percent = 1.-float(failedTests)/float(totalTests-skippedTests)
        percent = int(100. * percent)
        print("%d%s tests passed, %d tests failed out of %d (%d skipped)" %
              (percent, '%', failedTests, (totalTests-skippedTests), skippedTests))
    print('All tests passed? ' + str(success))
    if not success:
        sys.exit(1)
