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


Gui.addWorkbench(SubDModelerWorkbench())
