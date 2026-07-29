#!/usr/bin/env python3
"""Compatibility entry point for the portable GIMs Portal build stage."""

from __future__ import annotations

import sys

from gims_portal import main


if __name__ == "__main__":
    raise SystemExit(main(["build", *sys.argv[1:]]))
