#! /usr/bin/env python

# The purpose of this file is to define global constants and functions
# that are needed across modules.


# FITS header keywords and values

KEYWORDS_VELO = ["VOPT", "VRAD", "VELO", "FELO"]  # NOTE: "felo" is not a valid FITS coordinate code!
KEYWORDS_FREQ = ["FREQ"]

def check_header_keywords(values, keyword):
	for value in values:
		if value in keyword.upper(): return True
	return False
