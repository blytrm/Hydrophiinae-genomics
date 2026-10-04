#!/usr/bin/env python3
"""Run compleasm with its offline fallback repaired.

"""
import sys
import urllib.error

import compleasm

compleasm.URLError = urllib.error.URLError

if __name__ == "__main__":
    sys.exit(compleasm.main())
