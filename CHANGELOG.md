# Changelog

## 0.2.0 (2026-10-27)

First tagged release. Requires FreeCAD 1.0 or newer (developed on 26.3).

### New tools
- **From face**: turn selected faces or closed edge loops of Part / PartDesign objects into one SubD cage.
  Four-sided faces (flat or curved: cylinders, B-splines) become quad grids, flat faces become fans or quad
  strips, neighbouring faces are welded along shared edges, sharp joins stay creased. *Exact Edge* solves the
  cage so the SubD surface passes through the B-rep (about 0.08 mm on a filleted 4 mm box).
- **Snap** (bottom bar): free moves snap to the vertices and edges of visible Part / PartDesign shapes.
- **Make polygon** (P): from 3+ points, a closed edge loop, or two edges.
- **Bridge**: connect two edge loops or two face regions, with cuts, twist, smoothness and profile.
- **Hole**: a round hole in the middle of the selected faces, with automatic re-meshing, radius, segments and a
  support loop.
- **Mirror X / Y / Z**: mirror polygons (welded on the plane) or whole objects.
- **Copy / Cut / Paste** polygons (Ctrl+C / X / V), paste as a new object.
- **Box panel**: width / depth / height and X / Y / Z segments (default 100 mm).
- **Extrude** open edges, with segments (mouse wheel) and X / Y / Z axis lock (keys or panel).

### Interface
- LightWave-style **Create / Modify** tabs in the FreeCAD toolbar; tool menus that do not fit flow into columns.
- Knife (loop), Knife (line) and Select loop stay on (yellow) until turned off.
- T opens the move gizmo, G is the free mouse move; panel fields change with the mouse wheel.
- The Numeric panel is now the **Options** panel.
- Loop select follows open boundaries and grows from the selection (L).

### Fixes
- SubD objects lost their link to the code when FreeCAD 26.3 blocked the import while opening a document (the
  object mode showed an old shape); the module now loads from the Mod folder and broken objects are repaired.
- Closing the document closes the modeler (the tab toolbar no longer stays behind and the modeler opens again).
- Ctrl+S saves the current cage.
- Catmull-Clark: patch corners and creased corners on open boundaries stay sharp.

## 0.1.0

Initial public version: SubD cage modeling with extrude, bevel / smooth shift, edge bevel, loop cut, line knife,
crease, symmetry, gizmos, numeric input and a parametric B-rep (STEP) solid.
