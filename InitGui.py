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
        if any("GuidePoints" in o.PropertiesList and type(getattr(o, "Proxy", None)).__name__ != "RailPipe"
               for o in doc.Objects):                     # Rail Pipe / Groove (PartDesign)
            import subd_railpipe
            subd_railpipe.repair_proxies(doc)


class _PartDesignTools:
    """Rail Pipe / Rail Groove dugmeleri PartDesign'in "Part Design Modeling Features" arac cubugunda da olsun
    (Additive / Subtractive Pipe'in yaninda); FreeCAD 1.0+ workbench manipulator. Pimp My PartDesign kuruluysa
    dugmeleri o ekler (ayni arac iki kez gorunmesin); SubD_RailPipe komutlari eski belgeler icin yine kayitli."""

    def modifyToolBars(self):
        try:
            if Gui.activeWorkbench().name() != "PartDesignWorkbench":
                return []
            if "PDA_RailPipe" in Gui.listCommands():
                return []
        except Exception:
            return []
        return [{"append": c, "toolBar": "Part Design Modeling Features"} for c in ("SubD_RailPipe", "SubD_RailGroove")]


import FreeCAD as App
App.addDocumentObserver(_RestoreObserver())
Gui.addWorkbench(SubDModelerWorkbench())
if hasattr(Gui, "addWorkbenchManipulator"):
    Gui.addWorkbenchManipulator(_PartDesignTools())
