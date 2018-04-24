#! /usr/bin/env python

import sys
from time import time
from sofia import __version_full__ as sofia_version_full


# ==================================
# SETTINGS: Terminal colour settings
# ==================================

TERM_COL_RESET = "\x1B[0m"
TERM_COL_WARN  = "\x1B[33m"
TERM_COL_ERROR = "\x1B[31m"
TERM_COL_FATAL = "\x1B[35m"


# =====================================
# FUNCTION: Print informational message
# =====================================

def message(message, verbose=True):
	if verbose:
		sys.stdout.write(str(message) + "\n")
		sys.stdout.flush()
	return


# ===============================
# FUNCTION: Print warning message
# ===============================

def warning(message, fatal=False, frame=False):
	message = message.replace("\n", "\n         ")
	sys.stderr.write(TERM_COL_WARN + separator(frame) + "WARNING: " + str(message) + "\n" + separator(frame) + TERM_COL_RESET)
	if fatal: sys.exit(0)
	return


# =============================
# FUNCTION: Print error message
# -============================

def error(message, fatal=True, frame=False):
	if fatal:
		message = message.replace("\n", "\n             ")
		sys.stderr.write(TERM_COL_FATAL + separator(frame) + "FATAL ERROR: " + str(message) + "\n" + separator(frame) + TERM_COL_RESET)
		sys.exit(1)
	else:
		message = message.replace("\n", "\n       ")
		sys.stderr.write(TERM_COL_ERROR + separator(frame) + "ERROR: " + str(message) + "\n" + separator(frame) + TERM_COL_RESET)
	return


# =========================================================
# FUNCTION: Check condition and exit with signal 1 if false
# =========================================================

def ensure(condition, message, fatal=True, frame=False):
	if not condition: error(message, fatal=fatal, frame=frame)
	return


# ===============================
# FUNCTION: Return separator line
# ===============================

def separator(frame = True):
	if frame: return("______________________________________________________________________________\n\n")
	return ""


# ==========================
# FUNCTION: Print line break
# ==========================

def linebreak():
	sys.stdout.write("\n")
	return


# ==========================
# FUNCTION: Print time stamp
# ==========================

def print_progress_time(t0):
	message("    {0:.2f} seconds since start".format(time() - t0))
	return


# ================================
# FUNCTION: Print progress message
# ================================

def print_progress_message(msg, t0=None):
	msg = "--- {0:}: {1:} ".format(sofia_version_full, msg)
	msg = msg.ljust(78, "-")
	message("\n" + msg)
	if t0 is not None: print_progress_time(t0)
	linebreak()
	return