import sys

# Print informational message
def print_message(message, verbose=True):
	if verbose:
		sys.stdout.write(str(message) + "\n")
		sys.stdout.flush()
	return

# Print warning message
def print_warning(message):
	message = message.replace("\n", "\n         ")
	sys.stderr.write("WARNING: " + str(message) + "\n")
	return

# Print error message
def print_error(message, fatal=True):
	if fatal:
		message = message.replace("\n", "\n             ")
		sys.stderr.write("FATAL ERROR: " + str(message) + "\n")
		sys.exit(1)
	else:
		message = message.replace("\n", "\n       ")
		sys.stderr.write("ERROR: " + str(message) + "\n")
	return

# Check condition and exit with signal 1 if false
def ensure(condition, message):
	if not condition: print_error(message, fatal=True)
	return
