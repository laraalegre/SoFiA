import sys

# Check condition and exit with signal 1 if false
def ensure(condition, message):
	if not condition:
		sys.stderr.write("FATAL ERROR: " + str(message) + "\n")
		sys.exit(1)
	return

# Print informational message
def print_info(message, verbose=True):
	if verbose:
		sys.stdout.write(str(message) + "\n")
		sys.stdout.flush()
	return
