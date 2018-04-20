#! /usr/bin/env python

# The purpose of this file is to define global constants
# that are needed across different modules.


# FITS header keywords and values

KEYWORDS_VELO = ["VOPT", "VRAD", "VELO", "FELO"]  # NOTE: "felo" is not a valid FITS coordinate code!
KEYWORDS_FREQ = ["FREQ"]

def check_values(values, keyword):
	for value in values:
		if value in keyword.upper(): return True
	return False
