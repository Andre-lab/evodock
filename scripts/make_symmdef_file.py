#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Wrapper for creating symmetry files
@Author: Mads Jeppesen
@Date: 7/15/21
"""
from pathlib import Path
import subprocess

ROOT = Path("..").resolve()
BENCHMARK = ROOT.joinpath("benchmark")

script_path = Path("~/Rosetta/main/source/src/apps/public/symmetry/make_symmdef_file.pl").resolve()




