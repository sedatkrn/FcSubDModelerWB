# Changelog

## Unreleased

### New tools
- **From face** is now a drop-down with three commands. *Face*: as before. *Edge*: every chain of selected
  edges becomes one polygon, flat or not; ends that meet within 0.01 mm (also across objects) are joined, an
  open chain is closed with a straight edge (split by *Max Length*), separate chains stay separate polygons.
  *Point*: the selected points become one polygon in pick order (kept across objects). Edge and Point polygons
  face the camera.
- **Rail bridge** (Polygon): a tube along a live guide curve. Select two polygon regions / edge loops (or one
  region to extend from); the curve leaves each end perpendicular to its polygon and every ring in between
  stays perpendicular to the curve, carried along without twist, its shape blending from one end to the
  other. While editing, click a point and move it with the gizmo or the panel's X / Y / Z, set its ring
  *Scale*, *Segments* re-spaces the rings along the curve, I adds a ring in the clicked span, Del removes the
  selected one, Enter finishes. The guide is stored with the cage (undo, saved in the document): double-click
  the curve later to edit it again. If the tube is cut or rebuilt by hand, the guide detaches with a message.
  *Move neighbors* makes the guide act like a cable: the whole guide follows the moved point, bends smoothly
  (leaving the end rings perpendicular), keeps its segment lengths, pulls straight and stretches evenly when
  pulled further than its length. A free tail lies like a cable on a table: pulling a middle point raises a
  bulge there while the free end barely moves, until the tail is pulled taut and drags along. *Elasticity*
  (default 0.1) sets how easily the
  segments stretch: 0 chain, 1 rubber band. The end rings never move. Bridge and Rail bridge share one
  drop-down button.
- **Align** (two small icon buttons next to Snap in the bottom bar, also in the Options panel and the Rail
  bridge panel): lines up two or more selected points horizontally or vertically on screen, along the world
  axis nearest to the screen's up / right direction, at their mean. With polygons or edges selected their
  centers are aligned instead: touching selected elements move together as one piece, unchanged.
  While a Rail bridge guide is edited it aligns the guide points instead: the start (A) and end (B) points now
  show as handles too, Shift-click adds points to the selection, and if A or B is selected the other points
  line up with it (A and B themselves never move).
- Rail bridge: the guide no longer bulges after a long gap followed by a short one (tangents scale with each
  span's own length), and a new *Tension* field (0 smooth .. 1 straight between points, saved per guide) pulls
  the curve tighter where it still rounds out too much.
- Rail bridge: the panel's X / Y / Z fields of the selected guide point now work (the method was shadowed by a drawing node).
- Tool panel: optional button row under the fields.
- **Flip** (Modify, Shift+F; LightWave Flip): reverses the selected polygons (all when nothing is selected),
  so their normals point the other way. F stays Save.
- **Bridge > Auto direction** (on by default): an end made of free (pasted, unattached) polygons is turned to
  face the other end, so a smooth bridge no longer wraps around behind a copied polygon. Ends embedded in a
  surface keep their direction (handles over a surface still work); the strip takes its winding from the
  embedded end.
- **Reduce points** (Modify, LightWave Construct > Reduce): removes in-between points along polygon edges whose
  bend is below *Degrees* (default 5°), on the selection or the whole cage. Only points on a single edge line
  are taken, so neighbouring polygons stay joined; polygons never drop below triangles; creases carry over to
  the merged edge. *Delete points* removes the freed points from the cage.
- **Extrude > Individual faces (own normals)**: each selected face is extruded on its own, with its own side
  walls, along its own normal (Blender's Extrude Individual Faces), e.g. for knurled / ribbed rings on a
  cylinder. Off: the region moves along the averaged normal as before. The choice is remembered; an axis lock
  turns it off. It sits on its own row under the axis checks.
- **Shell** (Modify): object mode shells the whole selected object, edit mode the selected polygons. Thickness,
  inside / outside / centered, *Detail filter* (details narrower than this width are skipped on the inner side:
  short inner edges such as support loops and narrow steps are merged, and a robust local quadric fit flattens
  small dents without losing the main curvature), *Fill* joins inner and outer at open edges (holes, open rims),
  *Sharp* creases those rims, *Flat rim* keeps each filled rim in the plane of its opening (inner and outer edge
  line up, also on sloped walls). Solidify-style cage offset with even thickness and a crossing guard (no inner
  cage poking through at sharp, narrow corners); an excessive thickness (more than a quarter of the part's
  smallest size) is refused with a message. Selection borders stay open and are selected for weaving with P /
  Bridge.
- **Cylinder / Topo cylinder** (Create > Cylinder menu): diameter / height / segments / height segments in the
  panel. Cylinder has triangle-fan caps; Topo cylinder has all-quad caps (grid plus a ring, segments a multiple of 4).
- **SubD surface options** (Options panel): separate edit / final subdivision level, open edge rule (Edge and
  corner / Edge only), smooth triangles (Loop-like weights, smoother triangle fans such as cylinder caps).
- **Solid (STEP) precision and continuity**: Draft / Normal / Fine presets and G0 / G1 / G2. G1 and G2 convert
  regular patches exactly to bicubic Bezier patches and fit the rest to the limit surface (G2 refines irregular
  patches first); deviation from the SubD surface drops from about 0.03-0.06 mm to below 0.0001 mm, and it is faster.
  New solids use G2 / Normal; existing solids keep G0 until changed.
- **Zebra stripes** (toolbar button next to SubD to solid, or Options): reflection stripes on the SubD surface to
  check continuity.
- **Knife (line) options**: *Split in two* separates the model along the cut (open edges on both sides),
  *Split into two objects* moves the smaller side to a new SubD object. A cut limited to a selection still splits
  the whole model; a cut that does not go all the way around is not split into objects. Ctrl+Z also removes the
  new object.
- The edited object is marked (bold, yellow) in the model tree.
- **From face: flat faces with several holes** (e.g. a motor mounting plate): a quad ring around each hole, the
  area between the outer boundary and the rings triangulated (ear clipping with hole bridges, Delaunay flips) and
  merged into quads.
- **Gizmo size** (Options > Gizmo): length of the gizmo arrows on screen in px (30-400, default 90), remembered.
- Info and error messages are shown inside the 3D view, right of the tool column above the bottom bar (they were
  cut off in the bottom bar); long messages wrap.

### Fixes
- From face: hole walls (full cylinders with a seam) folded into a broken cage, alone or together with the plate:
  closed edges (circles) now take their direction from each face's own wire, and the periodic parameter runs
  continuously along each edge.
- From face *Exact Edge*: only boundary points (and curved four-sided faces) are solved; interior points of flat
  faces stay put, so irregular flat meshes no longer fold.
- Tool panel checkboxes are shown when the panel is built (the Shell panel was squashed).
- The tree marker no longer sticks to a second object after Knife splits rebuild the tree.
- Undo: autosaves no longer become FreeCAD undo steps; a whole edit session is one step ("SubD edit"), so Ctrl+Z
  in object mode undoes the last edit instead of an intermediate autosave. FreeCAD Undo / Redo while editing
  reloads the cage.
- Hole is centred on the true area centroid and aligned with the main axis (it was rotated a few degrees).
- Select loop in polygon mode follows selected polygons that are in the same strip but not adjacent.
- Tool column no longer splits into two columns when the modeler first opens; the Options panel scrolls instead
  of pushing the FreeCAD window off screen.

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
