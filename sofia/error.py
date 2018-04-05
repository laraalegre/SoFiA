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

def warning(message, border=False):
	message = message.replace("\n", "\n         ")
	separator(border)
	sys.stderr.write("WARNING: " + str(message) + "\n")
	separator(border)
	return


# -----------------------------
# FUNCTION: Print error message
# -----------------------------

def error(message, fatal=True, border=False):
	if fatal:
		message = message.replace("\n", "\n             ")
		separator(border)
		sys.stderr.write("FATAL ERROR: " + str(message) + "\n")
		separator(border)
		sys.exit(1)
	else:
		message = message.replace("\n", "\n       ")
		separator(border)
		sys.stderr.write("ERROR: " + str(message) + "\n")
		separator(border)
	return


# ---------------------------------------------------------
# FUNCTION: Check condition and exit with signal 1 if false
# ---------------------------------------------------------

def ensure(condition, message, fatal=True, border=False):
	if not condition: print_error(message, fatal=fatal, border=border)
	return


# ------------------------------
# FUNCTION: Print separator line
# ------------------------------

def separator(border = True):
	if border: sys.stderr.write("______________________________________________________________________________\n\n")
	return
