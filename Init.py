# -*- coding: utf-8 -*-
# SubD Modeler - application (non-GUI) start-up. Only makes sure the modules in this folder can be imported;
# saved SubD objects refer to the module "subd_modeler", so documents open even without the GUI part.
import os
import sys

try:
    _dir = os.path.dirname(os.path.abspath(__file__))
except NameError:                       # some FreeCAD versions do not set __file__; Mod/ folders are on sys.path
    _dir = None
if _dir and os.path.isfile(os.path.join(_dir, "subd_modeler.py")) and _dir not in sys.path:
    sys.path.insert(0, _dir)
