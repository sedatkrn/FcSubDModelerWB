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


COMMANDS = ["SubD_Modeler", "SubD_MakeSolid"]

Gui.addCommand("SubD_Modeler", _OpenModeler())
Gui.addCommand("SubD_MakeSolid", _MakeSolid())
