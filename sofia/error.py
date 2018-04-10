import sys


# -------------------------------------
# FUNCTION: Print informational message
# -------------------------------------

def message(message, verbose=True):
	if verbose:
		sys.stdout.write(str(message) + "\n")
		sys.stdout.flush()
	return


# -------------------------------
# FUNCTION: Print warning message
# -------------------------------

def warning(message, fatal=False, frame=False):
	message = message.replace("\n", "\n         ")
	separator(frame)
	sys.stderr.write("\033[33mWARNING: " + str(message) + "\033[0m\n")
	separator(frame)
	if fatal: sys.exit(0)
	return


# -----------------------------
# FUNCTION: Print error message
# -----------------------------

def error(message, fatal=True, frame=False):
	if fatal:
		message = message.replace("\n", "\n             ")
		separator(frame)
		sys.stderr.write("\033[35mFATAL ERROR: " + str(message) + "\033[0m\n")
		separator(frame)
		sys.exit(1)
	else:
		message = message.replace("\n", "\n       ")
		separator(frame)
		sys.stderr.write("\033[31mERROR: " + str(message) + "\033[0m\n")
		separator(frame)
	return


# ---------------------------------------------------------
# FUNCTION: Check condition and exit with signal 1 if false
# ---------------------------------------------------------

def ensure(condition, message, fatal=True, frame=False):
	if not condition: error(message, fatal=fatal, frame=frame)
	return


# ------------------------------
# FUNCTION: Print separator line
# ------------------------------

def separator(frame = True):
	if frame: sys.stderr.write("______________________________________________________________________________\n\n")
	return
