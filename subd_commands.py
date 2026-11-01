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


class _RailPipe:
    """PartDesign Body'de Rail Pipe (ekler) / Rail Groove (oyar): secili bir ya da iki duz yuz arasinda kilavuz egri
    boyunca supurme (subd_railpipe)."""

    def __init__(self, subtractive):
        self.sub = subtractive

    def GetResources(self):
        if self.sub:
            return {"Pixmap": os.path.join(ICONS, "subd_railgroove.svg"), "MenuText": T("Rail Groove"),
                    "ToolTip": T("PartDesign: cut a groove along a guide curve from a flat face or closed sketch (and "
                                 "to a second one). Edit GuidePoints, Tension and EndScale in the property view")}
        return {"Pixmap": os.path.join(ICONS, "subd_railpipe.svg"), "MenuText": T("Rail Pipe"),
                "ToolTip": T("PartDesign: add a pipe along a guide curve from a flat face or closed sketch to a "
                             "second one or a free end. Edit GuidePoints, Tension and EndScale in the property "
                             "view")}

    def IsActive(self):
        return App.ActiveDocument is not None

    def Activated(self):
        import subd_railpipe
        prof = []                                          # (nesne, "FaceN") ya da sketch icin (sketch, None)
        for so in Gui.Selection.getSelectionEx():
            if subd_railpipe._is_sketch(so.Object):       # sketch: kenarina tiklanmis ya da agacta secili
                if all(p[0] is not so.Object for p in prof):
                    prof.append((so.Object, None))
            elif so.SubElementNames:
                prof += [(so.Object, sub) for sub in so.SubElementNames if sub.startswith("Face")]
            elif hasattr(so.Object, "Shape") and len(so.Object.Shape.Faces) == 1:
                prof.append((so.Object, "Face1"))           # agacta secilen tek yuzlu nesne (or. Binder)
        if not 1 <= len(prof) <= 2:
            App.Console.PrintWarning(T("Rail Pipe: select one or two flat faces or closed sketches of a "
                                       "PartDesign Body.\n"))
            return
        obj = prof[0][0]
        body = obj.getParentGeoFeatureGroup() if hasattr(obj, "getParentGeoFeatureGroup") else None
        if body is None or body.TypeId != "PartDesign::Body":
            App.Console.PrintWarning(T("Rail Pipe: the faces must belong to a PartDesign Body.\n"))
            return
        doc = obj.Document
        doc.openTransaction(T("Rail Groove") if self.sub else T("Rail Pipe"))
        new = subd_railpipe.make(body, prof[0], prof[1] if len(prof) > 1 else None, self.sub)
        if isinstance(new, str):
            doc.abortTransaction()
            App.Console.PrintWarning(T(new) + "\n")
            return
        doc.recompute()
        doc.commitTransaction()
        for p in prof:                                     # PartDesign gibi: kullanilan sketch'ler gizlenir
            if p[1] is None:
                p[0].ViewObject.Visibility = False
        if "Invalid" in new.State:
            App.Console.PrintWarning(T("Rail Pipe: the first sweep failed; adjust GuidePoints in the property "
                                       "view.\n"))
        elif new.BaseFeature is not None:
            new.BaseFeature.ViewObject.Visibility = False
            new.ViewObject.Visibility = True
        Gui.Selection.clearSelection()
        Gui.Selection.addSelection(new)


COMMANDS = ["SubD_Modeler", "SubD_MakeSolid", "SubD_Zebra", "Separator", "SubD_RailPipe", "SubD_RailGroove"]

Gui.addCommand("SubD_Modeler", _OpenModeler())
Gui.addCommand("SubD_MakeSolid", _MakeSolid())
Gui.addCommand("SubD_Zebra", _Zebra())
Gui.addCommand("SubD_RailPipe", _RailPipe(False))
Gui.addCommand("SubD_RailGroove", _RailPipe(True))
