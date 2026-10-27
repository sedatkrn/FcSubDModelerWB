# -*- coding: utf-8 -*-
# SubD Modeler - FreeCAD workbench registration. Copyright (c) 2026 Sedat Kuran - MIT License
import FreeCADGui as Gui


class SubDModelerWorkbench(Gui.Workbench):
    MenuText = "SubD Modeler"
    ToolTip = "LightWave-style subdivision (SubD) modeling with B-rep / STEP output"

    def __init__(self):
        import subd_commands
        self.__class__.Icon = subd_commands.ICON

    def Initialize(self):
        import subd_commands
        self.appendToolbar("SubD", subd_commands.COMMANDS)
        self.appendMenu("SubD", subd_commands.COMMANDS)

    def Activated(self):
        pass

    def Deactivated(self):
        pass

    def GetClassName(self):
        return "Gui::PythonWorkbench"


class _RestoreObserver:
    """Belge acilinca: Proxy'si kopmus SubD nesnelerini onar (FreeCAD modulu acilista ice aktaramadiysa yuzey ve
    katilar guncellenmez). Modul yalniz belgede SubD nesnesi varsa yuklenir."""

    def slotFinishRestoreDocument(self, doc):
        if any(("CageFaces" in o.PropertiesList or "Samples" in o.PropertiesList)
               and getattr(o, "Proxy", None) is None for o in doc.Objects):
            import subd_modeler
            subd_modeler.repair_proxies(doc)


import FreeCAD as App
App.addDocumentObserver(_RestoreObserver())
Gui.addWorkbench(SubDModelerWorkbench())
