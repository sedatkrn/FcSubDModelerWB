# -*- coding: utf-8 -*-
# SubD Modeler - application (non-GUI) start-up. Only makes sure the modules in this folder can be imported;
# saved SubD objects refer to the module "subd_modeler", so documents open even without the GUI part.
import os
import sys

import FreeCAD

try:
    _dir = os.path.dirname(os.path.abspath(__file__))
except NameError:                       # some FreeCAD versions do not set __file__; Mod/ folders are on sys.path
    _dir = None


def _same(a, b):
    return os.path.normcase(os.path.realpath(a)) == os.path.normcase(os.path.realpath(b))


if _dir and os.path.isfile(os.path.join(_dir, "subd_modeler.py")):
    # Klasor Mod'a baglanti (junction / symlink) ile kuruluysa modul Mod yolundan yuklensin: FreeCAD 26.3+ belge
    # acilirken yalniz FreeCAD.__ModDirs__ altindaki modullere izin verir; gercek yoldan (or. Documents\...)
    # yuklenen modul engellenir ve SubD nesneleri Proxy'siz acilir.
    _dir = next((d for d in getattr(FreeCAD, "__ModDirs__", []) if _same(d, _dir)), _dir)
    sys.path[:] = [p for p in sys.path if not (p and p != _dir and os.path.isdir(p) and _same(p, _dir))]
    if _dir in sys.path:
        sys.path.remove(_dir)
    sys.path.insert(0, _dir)
