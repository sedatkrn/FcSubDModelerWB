# -*- coding: utf-8 -*-
# SubD Modeler - FreeCAD commands (toolbar / menu). Copyright (c) 2026 Sedat Kuran - MIT License
import os

import FreeCAD as App
import FreeCADGui as Gui

from subd_i18n import T

ICONS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "icons")
ICON = os.path.join(ICONS, "subd_modeler.svg")


def _selected_subd():
    import subd_modeler
    sel = Gui.Selection.getSelection()
    return sel[0] if len(sel) == 1 and subd_modeler.is_subd(sel[0]) else None


class _OpenModeler:
    def GetResources(self):
        return {"Pixmap": ICON,
                "MenuText": T("SubD Modeler"),
                "ToolTip": T("Open the SubD modeler (edits the selected SubD object, or starts with a box)")}

    def IsActive(self):
        return True

    def Activated(self):
        import subd_modeler
        subd_modeler.start()


class _MakeSolid:
    def GetResources(self):
        return {"Pixmap": os.path.join(ICONS, "subd_solid.svg"),
                "MenuText": T("SubD to solid"),
                "ToolTip": T("Make a parametric B-rep solid (STEP) from the selected SubD object")}

    def IsActive(self):
        return App.ActiveDocument is not None

    def Activated(self):
        import subd_modeler
        obj = _selected_subd()
        if obj is None:
            App.Console.PrintWarning(T("Select a SubD object first.\n"))
            return
        obj.Document.openTransaction("SubD solid")
        subd_modeler.make_solid(obj)
        obj.Document.commitTransaction()


class _Zebra:
    """Zebra seritleri (basili kalan dugme): modelleyicinin puruzsuz yuzeyinde yuzey surekliligi analizi."""

    def GetResources(self):
        return {"Pixmap": os.path.join(ICONS, "subd_zebra.svg"),
                "MenuText": T("Zebra stripes"),
                "ToolTip": T("Reflection stripes on the surface: kinks show as broken stripes (G0), corners in the "
                             "stripes as G1, smooth stripes as G2"),
                "Checkable": False}

    def IsActive(self):
        import sys
        m = sys.modules.get("subd_modeler")
        return m is not None and m._session is not None

    def Activated(self, checked=None):
        import subd_modeler
        s = subd_modeler._session
        if s is not None and (checked is None or bool(checked) != s.zebra):
            s.toggle_zebra()


COMMANDS = ["SubD_Modeler", "SubD_MakeSolid", "SubD_Zebra"]

Gui.addCommand("SubD_Modeler", _OpenModeler())
Gui.addCommand("SubD_MakeSolid", _MakeSolid())
Gui.addCommand("SubD_Zebra", _Zebra())
