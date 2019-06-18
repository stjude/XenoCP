#!/bin/sh
# Provides functions for use by QC scripts.
#
# IMPORTANT: This script is intended to be sourced.
#
# The overall paradigm here is similar to that of a unit test.  Every step you
# want to test will get a "test case" which is simply a shell script that runs
# tests.  This script provides helper functionality to the test case scripts.
# The contract of the test case scripts is that they print to stdout lines
# stating the success/failure of tests, and nothing else (as far as possible),
# concluding with summary information.
#
# So, test case scripts should use this script as follows:
# 1. Source it: source qclib.sh
# 2. Call starttestcase <test case name>
# 3. For each test in the test case:
#    a. Call starttest <test name>
#    b. Determine success or failure
#    c. Call pass, fail, or warn
# 4. At the end, call summarize
#
# The fastest way to interpret results is to look only at the last line; in
# fact, the first word of the last line will be PASS, WARN, or FAIL.  If you
# want to be really safe, you should also check to make sure that the last
# line really is the summary line.  The easiest way to do this is to check
# the first word of the second-to-last line; it should be SUMMARY.

# Prints a message
# $1 = pass, fail, warn, etc.
# $2 = test or test case name
# $3... = message (optional)
function printmessage {
  COL1=$1
  COL2=$2
  shift 2
  echo -e "$COL1\t$COL2\t$*"
}

# Prints a message for the current test
# $1 = pass, fail, warn, etc.
# $2... = message (optional)
function printtestmessage {
  COL1=$1
  shift 1
  printmessage $COL1 $CURTEST $*
}

# Starts a test case
# $1 = test case name
# $2... = message (optional)
function starttestcase {
  # Get test case name
  CURTESTCASE=$1
  # Print header to stdout
  shift 1
  printmessage HEADER $CURTESTCASE $*
  # Print a header to stderr
  echo "==================== $CURTESTCASE" 1>&2
  # Init variables
  TESTS_STARTED=0
  TESTS_PASSED=0
  TESTS_FAILED=0
  TESTS_WARNED=0
}

# Starts a test
# $1 = test name
function starttest {
  # Print a header to stderr
  echo "--------------- $1" 1>&2
  # Init variables
  CURTEST=$1
  let "++TESTS_STARTED"
}

# Passes the current test, if there is a current test
# $1... = message (optional)
function passtestbydefault {
  if [ -n "$CURTEST" ]; then passtest $*; fi
}

# Passes the current test
# $1... = message (optional)
function passtest {
  let "++TESTS_PASSED"
  printtestmessage PASS $*
  CURTEST=
}

# Warns for the current test, if there is a current test
# $1... = message (optional)
function warntestbydefault {
  if [ -n "$CURTEST" ]; then warntest $*; fi
}

# Warns for the current test
# $1... = message (optional)
function warntest {
  let "++TESTS_WARNED"
  printtestmessage WARN $*
  CURTEST=
}

# Fails the current test, if there is a current test
# $1... = message (optional)
function failtestifactive {
  if [ -n "$CURTEST" ]; then failtest $*; fi
}

# Fails the current test
# $1... = message (optional)
function failtest {
  let "++TESTS_FAILED"
  printtestmessage FAIL $*
  CURTEST=
}

# Fails the current test, summarizes the test case and exits the script
# $1... = message (optional)
function aborttestcase {
  if [ -n "$CURTEST" ]; then failtest $*; fi
  summarize
  exit
}

# Prints summary information for the test case
function summarize {
  # Check counts
  BAD_TESTS_STARTED=
  let "CHECK = TESTS_PASSED + TESTS_FAILED + TESTS_WARNED"
  if [ "$CHECK" -ne "$TESTS_STARTED" ]
  then
    BAD_TESTS_STARTED=$TESTS_STARTED
  fi
  starttest CheckCounts
  if [ -n "$BAD_TESTS_STARTED" ];
  then
    failtest "Started $BAD_TESTS_STARTED tests, but finished $CHECK tests"
  else
    passtest "Test counts are consistent"
  fi
  
  # Print counts
  printmessage SUMMARY $CURTESTCASE "Started: $TESTS_STARTED  Passed: $TESTS_PASSED  Failed: $TESTS_FAILED  Warnings: $TESTS_WARNED"
  
  # Compute what will be the return value before the final pseudo-test
  retval=$TESTS_FAILED
  
  # Final message reflects worst case (do a pseudo-test using the test case 
  # name)
  starttest $CURTESTCASE
  if [ "$TESTS_FAILED" -gt 0 ]
  then
    failtest "At least one test failed"
  elif [ "$TESTS_WARNED" -gt 0 ]
  then
    warntest "No tests failed, but there was at least one warning"
  else
    passtest "All tests passed"
  fi
  
  return $retval
}

# Helper functions

# Get filesize; alias for stat -c '%s'
# $1 = filename
function filesize {
  stat -L -c '%s' $*
}
