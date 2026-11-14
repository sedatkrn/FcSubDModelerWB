# Changelog

## Unreleased

- **New: Sphere and Topo sphere (Create tab).** *Sphere* is the classic UV sphere (poles on Z, triangle fans at the
  poles; panel: Diameter, Segments, Rings). *Topo sphere* is a cube-sphere made only of quads (panel: Diameter,
  Segments 2-16), which gives an even SubD surface without pole pinching.
- **New: Torus (Create tab).** Ring / tube diameters and both segment counts in the panel; all quads.
- **New: Cage B-rep (File tab).** Adds the control cage itself to the document as a parametric B-rep with flat
  faces (a Solid when the cage is closed, otherwise a Shell). Polygons that are not flat are split into triangles
  (`PlanarTolerance`), neighbouring coplanar faces are merged (`MergeFlatFaces`). Meant for sending a hand-built cage
  straight into FreeCAD, e.g. for jewellery or mechanical parts. Updates when the cage changes.

## 0.5.0 (2026-11-14)

- **Changed: parametric groups no longer use FreeCAD datum planes.** Existing groups keep their points and colour;
  their datum is removed. A group can get a **reference point** (panel: *Set reference point*, select 1, 2 or 4
  points; the middle is used) shown as a sphere.
- **New: PD (Parametric Dimension)** in the top bar. Shows the reference spheres and the dimension types
  Horizontal / Vertical (along the world axis nearest to the screen direction), Aligned and Angle (3 spheres).
  Click the spheres in order: the first is fixed (red), the next ones (green) move when you type a new value into the
  dimension's text. Dimensions are saved with the object.
- **PD: stretch.** A group can have several reference spheres (each *Set reference point* adds one). A dimension
  between two spheres of the *same* group stretches the group along the measured direction from the first (red)
  sphere, so a width can be changed; between different groups the second group still moves.
- **PD dimensions:** the second sphere's extension line goes straight to the dimension line; drag a dimension's text
  to move the dimension line sideways (undoable); double click the text to type a new value. PD dimensions also
  appear, in orange, in the Distances list.
- **Fix: Undo right after adding a Box / Cylinder removes the new object.**
- Topo cylinder segments: 8 to 32 in steps of 4.
- **New: Circle (Modify tab).** Arranges the selected points, an edge loop / chain, or the boundary of a polygon
  region (one region without holes) on a circle in their best-fit plane; the region's inner points follow
  proportionally. Panel: *Radius* (type a number; *Auto radius* turns off), *Strength* %, *Auto radius*,
  *Even spacing*. No default shortcut.

## 0.4.3 (2026-11-10)

- **Fix: FreeCAD crashed when a document was closed while SubD mode was open** (closing the tab or the file,
  especially after an object was created). The modeler no longer touches the already destroyed 3D view, its widgets
  or its scene when the document goes away; the SubD tool panel is still removed.
- **Fix: parametric groups could get mixed up after Knife (line).** Unmerged points at the same place now belong
  to the group together and move together; a point can no longer end up in two groups.

## 0.4.2 (2026-11-10)

- **New: parametric groups.** Select polygons, edges or points and press *Parametric group*: the group gets a FreeCAD
  datum plane (`SubD_PG_NNNN`) in the tree and its own colour (also shown in the tree's Description column). Moving
  the datum (Placement, Attachment, expressions) moves the group's points.
- **Distances between groups.** Attach a group's datum to another group's datum, an LCS, a plane or a face: a
  dimension with arrows is drawn (offset from the centre, alternating sides), its value is bound to a `Distance`
  property through a FreeCAD expression and can be edited in the 3D view or in the Distances list.
- **Locked groups.** Groups are locked by default. Selecting a group's points opens a panel with *Unlock / Lock*,
  *Add to group* and *Remove from group*. Moving or rotating one member moves the whole group and its datum;
  scaling is blocked until the group is unlocked. Undo / Redo restore the datum too.
- Attaching a group datum keeps its current place instead of jumping onto the reference.
- Group datums are shown under their SubD object in the tree. The SubD does not depend on its datums in the FreeCAD
  graph, so a datum can reference the SubD itself, its faces or an LCS attached to it without a dependency cycle.

## 0.4.1 (2026-11-07)

- **Fixed: two documents open.** Each document now has its own modeler session. Switching to another document tab
  suspends the modeler (panel, bars, capture hide) and switching back resumes it with its Undo / Redo history;
  starting the modeler in a second document no longer closes the first one's session.
- **Fixed:** a closed session could come back to life when switching tabs, and a new session deleted the tool panel of
  the other document's session. A suspended session no longer autosaves into the wrong document.
- **Fixed: bars overlapping FreeCAD's overlay panels.** When the Tree / Property view is used in Overlay mode (floating
  over the 3D view), the Object / Edit bar, the bottom bar, the tool options panel and messages now move aside and
  follow the overlay panels as they open, close or resize.
- **Fixed:** Qt warning `objectName not set for QDockWidget` in the Report View (the Options panel had no name).
- **New: right-click a SubD object** in the tree (SubD Modeler workbench) and choose *SubD mode* to edit it.

## 0.4.0 (2026-11-07)

- **Fixed: switching to another workbench left the SubD tools on screen.** The *SubD tools* panel, the bars in the 3D
  view and the mouse / key / Undo capture stayed active in the other workbench and got in its way. Now the modeler
  is **suspended** when another workbench becomes active: the panel and bars hide, the capture is released and an unfinished
  tool is cancelled, but the session and the cage stay (sketch links keep following, autosave keeps running). It comes
  back as it was when you return to the SubD workbench. The modeler is shown in the SubD workbench and in the workbench
  it was started from (for example with the macro).
- **Distance** (top bar next to Object / Edit, key **U**; assignable with `subd_modeler.set_user_key("dim", "...")`):
  measures the distance between two points in **Edit and Object mode**. The snap is always on (independent of Q):
  cage points, edge middles, polygon centers and edges, corners / edges of other objects, and finally the surface of
  the model (SubPatch surface or other objects). The panel shows the distance and dX / dY / dZ.
  - **Keep (fixed)** off (default): the measurement is temporary and goes away with the next one or when you leave the tool.
    In Object mode measurements are always temporary.
  - **Keep (fixed)** on (Edit mode): the measurement stays on screen, is **bound to the cage points** (it updates live
    when you move, rotate, twist or array them), is saved with the object, works with Undo / Redo, and is deleted with the
    **×** on the left of its text (works with or without the tool; hover turns it red). *Clear all* deletes them all. Points on the surface (not on the cage) are temporary only.
  - The measurements are drawn inside the 3D scene, so they stay in sync while you orbit / zoom. The text size follows the gizmo size setting. Esc / right click cancels the first point, then leaves the tool.
- **Distances list** (the **≡** button next to Distance in the top bar): a window with all fixed distances of the
  object. Give each distance a **name** (shown in front of the value in the 3D view, switch with *Show names*), an optional
  **target** and **tolerance** (the value turns red when it is out of tolerance, the *Dev.* column shows the deviation),
  **show / hide** it (checkbox, *Show all* / *Hide all*) or delete it. The selected row is drawn blue in the view.
  Names and targets are saved with the object, so the list is there in the next session. **Export CSV...** writes the
  list to a file. In Object mode the list shows the distances of all SubD objects (read-only). Edits are undoable.
- **Twist** (Modify > Twist; LightWave Twist): the selected polygons are twisted around an axis. Press Twist, **move the
  mouse** for the angle (right = clockwise, Ctrl = steps), click to accept (Esc cancels), then fine-tune in the panel.
  Edit mode only, symmetry is ignored.
  - **Falloff** along the axis: the angle grows from 0 at *Start* to the full *Angle* at *End*. *Linear*, *Smooth*,
    *Ease in*, *Ease out* keep one end fixed (*Reverse* picks the other end); *Peak* and *Dome* keep both ends fixed and
    twist the middle the most. *Start* / *End* are mm along the axis from the axis point; *Auto range* (default) uses the
    ends of the selection, typing a Start or End switches it off.
  - **Cuts**: the selected polygons are first cut by *Cuts* evenly spaced planes between Start and End (perpendicular
    to the axis), so the twist looks smooth. Neighbouring polygons take the new points, the cage stays closed.
  - **The axis** works as in Array: *Draw axis* (drag a line with snapping, a click gives an axis along the view),
    *Axis X / Y / Z*, or the axis point and direction as numbers. It starts through the centre of the selection along
    its longest world axis.
- **Assignable shortcuts**: Array and Twist have no default key, but you can give them one. In FreeCAD's Python console:
  `import subd_modeler; subd_modeler.set_user_key("twist", "Shift+T")` (an empty string removes it). It works at once
  in an open modeler, the button shows the key the next time the modeler opens, and it wins over a built-in key that
  uses the same combination.

- **Array** (Polygon > Array, no shortcut; LightWave Array): the selected polygons are copied as **separate parts**
  (the original stays where it is; *Merge points* can weld them afterwards). Edit mode only, symmetry is ignored.
  One button: it starts as Radial around the world Z axis through the origin, and *Radial / Rectangular* is switched
  in the panel.
  - **Radial**: *Count* copies between *Start angle* and *End angle* around an axis; a full circle (360°) skips the
    last copy so it does not land on the first. *Offset along axis* lifts each copy (spiral). *Rotate copies* turns
    the copies with the circle (off: they only move).
  - **Rectangular**: *Count X / Y / Z* and *Offset X / Y / Z* along the world axes.
  - **Offset type** (both modes): *Between* (distance between copies), *Total* (distance from the first to the last),
    *Multiply* (the gap is multiplied by *Factor* at each step).
  - **The axis is a line, which is the difference from LightWave** (X / Y / Z only there): *Draw axis* in the panel
    lets you drag a line in the view, with the same snapping as Mirror (line) (points, edge middles, edges, Part /
    PartDesign corners; Ctrl = 15° steps, Ctrl also turns snapping off). A short click gives an axis along the view
    direction through the clicked point (Top view: a vertical axis). Or type the axis point and direction in the
    panel, or press *Axis X / Y / Z*. Every value updates the result live and Esc / Undo bring it back.

## 0.3.3 (2026-11-03)

- **Add points** (Polygon): the selected edges are split into equal parts (*Count* in the panel). With nothing selected
  the tool stays on: a click on an edge adds a point there, snapping to the middle (exact *Position* % in the panel).
  Neighbouring polygons take the new points, and creased edges stay creased. (#3)
- **Connect** (Ctrl+L, LightWave Split): adds an edge between two selected points of one polygon, or between the
  middles of two selected edges of one polygon. The polygon is split in two. (#3)
- **Subdivide** (Shift+D, LightWave Subdivide): the selected polygons are split into quads (n sides become n quads).
  Neighbouring polygons take the edge middles, so the cage stays closed. (#6)
- **From edge polygons follow their sketch**: when a Sketcher sketch (or other source edges) changes size or
  position, a cage that is still just the flat polygon is rebuilt on the new curve. In an extruded or bevelled cage
  the linked points and the points straight above them move with the sketch, and other points stay where they are.
  Cancelling a sketch edit brings the cage back. Options > *Follow sketch* (or the object's *SketchLink* property)
  turns the link off.
- From edge: straight edges stay one polygon edge, and curves are sampled by *Arc Segments*.
- The modeler **opens empty** in object mode, with no default box, and shows "SubD Modeler is ready". (#7)
- **Object / Edit and Point / Edge / Polygon** buttons moved to the top-left corner of the view, as in Blender, so
  the bottom bar fits in one row. On narrow views the bottom bar hides labels step by step instead of cutting them off.
- Knife (loop) and Select loop now pass through polygons that were merged with Merge polygons. (#5)

## 0.3.2 (2026-11-02)

- Matcap and environment reflection in the modeler: section planes now cut the surface too.

## 0.3.1 (2026-11-02)

- **Undo goes step by step everywhere while editing a cage**: FreeCAD's Undo / Redo buttons, the Edit menu and Ctrl+Z with the focus outside the 3D view now undo one modeling step. They used to take back the whole edit and lose the steps. The arrow on FreeCAD's toolbar button still lists the document's own steps.
- **Redo**: Ctrl+Y or Ctrl+Shift+Z. **Undo / Redo buttons** are on the bottom bar.
- The steps survive a trip to object mode and back, as long as the object did not change in between.
- Ghost mode is about 3x faster while the view or a drag is moving: the transparent surface is sorted triangle by triangle
  only when everything stands still, so the still image stays clean.
- FreeCAD 27: tool and bottom bar buttons no longer get cut off (the theme forced every button to the same width).
- Magnet: the right mouse button no longer opens FreeCAD's context menu.
- **Tool tips** on all tool, menu and bottom bar buttons (English and Turkish).
- **Circle on face** (Polygon, under Hole): works like Hole, but the circle stays filled and comes out selected, ready for
  Extrude or Bevel. It takes several separate patches too.
- Shell *Detail filter*: merging short inner edges could leave an edge shared by four polygons, which can make the solid
  invalid. Points that would do this are no longer merged.
- Shell *Centered* with *Flat rim*: the outer rim edge now stays in the plane of the opening too. It used to dip with the
  surface normal (0.5 mm on a mouse housing with a 2 mm wall), so the rim was not flat.
- For the FreeCAD Addon Index rules: the workbench no longer changes `sys.path`, and it no longer adds Rail Pipe / Rail
  Groove buttons to PartDesign's toolbar.
- The Rail Pipe / Rail Groove buttons are gone from the SubD toolbar and menu too: install the **Pimp My PartDesign**
  add-in to create them. Old documents with these features still open and recompute.

## 0.3.0 (2026-11-02)

- **Snap** can be switched on and off with the **Q** key.
- **Hole** now works on several separate polygon patches at once, one hole in each (for example opposite walls). The
  radius and segments are shared, so the holes match; if one patch cannot take the hole, nothing changes.
- **Magnet** (Modify): drag on the cage to pull or push it softly. Points inside the radius (shown as a ring) move by
  their distance from the pressed spot; a selection limits it to the selected points. Falloff: Smooth, Linear, Peak,
  Dome. **Spherical** measures in 3D, **Radial** ignores depth, so it affects everything inside the ring through the
  whole model along the view direction. Dragging with the right mouse button (or the wheel during a pull) sets the
  radius; a right click without moving closes the tool. Ctrl while pulling keeps the move straight, horizontal or
  vertical on screen. With Snap on, the pressed spot and the end of the pull catch the cage's own points and edges (not the
  ones being moved) and Part / PartDesign corners and edges; the target and the ring show before you press.
  The panel re-applies Radius, offset and falloff afterwards.
  Works with symmetry; it only moves points, so no polygons are added.
- **Align** now lines everything up with the first clicked point (polygon, edge, guide point), which stays in place,
  instead of the mean. A selection made only with a box still aligns at the mean, since it has no first click.
- Rail Pipe / Rail Groove now live in the **Pimp My PartDesign** add-in (github.com/sedatkrn/partdesign_addins).
  When it is installed, this workbench no longer adds its own Rail Pipe buttons to PartDesign's toolbar (no
  duplicates); old documents keep opening, and the add-in links their Rail Pipe objects to itself.
- Deleting every polygon of an object now removes the object (with its solid) from the tree when you leave its cage:
  it used to stay there with its old shape. Ctrl+Z still works inside the edit, and FreeCAD's Undo brings the object
  back afterwards.
- The tool column and the Create / Modify tabs now live in a FreeCAD panel (**SubD tools**) instead of on top of the 3D
  view: it opens as a tab next to Model / Tasks and can be dragged to a side or floated, and its place is remembered.
  The bottom bar, numeric panel and messages stay in the view.
### New tools
- **Cage color** (Options): color of the cage faces when Subpatch is off; remembered.
- **Environment reflection** (Options): the smooth preview in Edit mode reflects a loaded panorama, read per pixel and
  wrapped on a sphere around the model (fixed to the world, Z up), so flat faces show an image instead of streaks. *Roughness* blurs the reflection and *Metal* goes
  from plastic (surface color lit by the image, with a soft reflection stronger at the edges; *Surface color* sets it) to metal. Use an equirectangular (2:1) image: Radiance .hdr is best, JPG and PNG also work. HDR is tone
  mapped to 8 bits. Image and setting are remembered. The surface does not reflect itself or cast shadows.
- **Matcap** (bottom bar): surface looks for the smooth preview in Edit mode, with light, highlights and darkened edges
  painted into a texture: Clay, Plastic, Black plastic, Ceramic, Metal. The texture is read per pixel, on smooth normals
  of the surface (creases and very sharp breaks stay hard), so no blocks appear in curved areas. The choice is
  remembered. Zebra takes over while it is on.
- **Merge polygons** (Polygon, Shift+Z, LightWave Merge Polygons): selected adjacent polygons become one polygon;
  inner edges and points go, creases on the outline stay. Polygons that are not adjacent, or a selection with no single
  outline (a ring with a hole, a whole closed shape, mixed directions), give a warning and leave the cage untouched.
- **Knife (line)** snaps to the cage when Snap is on: both ends of the line catch points, edge midpoints and edges
  (points and midpoints win over edges). Over a target the cursor turns into a cross and a mark appears (square:
  point, triangle: midpoint, X: edge), and the cut passes exactly through it. Ctrl (15° steps) turns snapping off for
  the end of the line.
- **Merge points** (Modify, key M, LightWave Detail > Points > Merge Points): joins points on the selection or the
  whole cage. *Automatic* takes points in the same place; *Fixed distance* takes points within *Distance*. A point
  moves onto the other (no averaging). Polygons that collapse are removed, a polygon that touches itself is split,
  and two polygons left back to back on the same points are removed, so two closed pieces become one.
- **Weld** (Modify, Ctrl+W, LightWave Detail > Points > Weld): selected points join at the last one clicked, which
  stays in place, with the same polygon clean-up as Merge points. With Symmetry on, the mirror points weld too. A box
  selection has no order, so Shift+click the target last.
- Rail Pipe: a profile face that belongs to no solid (SubShapeBinder of a face, a lone face) no longer sends
  the pipe out of its back side: like a sketch, it leaves towards the other profile. A binder or sketch that
  lies on the Body's surface is treated like the Body's own face (outward, end carried into the body, so the
  fuse no longer splits or comes out empty). A single-face object picked in the tree is accepted.
- **Rail Pipe / Rail Groove** (toolbar, also in PartDesign's Modeling Features toolbar): PartDesign features in
  a Body. From one flat face or closed sketch to a second one or a free end (an empty Body works too), the start
  profile is swept (OCCT pipe shell) along a guide curve built like the modeler's
  Rail bridge guide (perpendicular to the faces, *Tension*, no bulge on uneven spacing) and fused into or cut
  out of the Body. Properties: *GuidePoints*, *Tension*, *EndScale* (free end), *Reversed*. Double-click edits the guide on
  screen: the selected point gets the modeler's Move gizmo (axis arrows, plane squares, free cube) or is
  dragged directly (live curve, the pipe recomputes as you drag; a missed button release no longer leaves the
  point stuck to the mouse), while editing, the Body shows as with Pad / Pocket: the previous solid as usual and the added (green) or cut
  (red) shape as FreeCAD's own translucent preview; click the curve to pick a span (*Add point* inserts exactly
  where you clicked), delete points, *Move neighbors*
  (cable, with *Elasticity*), X / Y / Z of the selected point; OK is one undo step, Cancel reverts. The ends are
  carried a little past the faces so no coincident faces are left for the boolean, and the sweep falls back
  through several frame / continuity modes until it gives a valid solid. The cross-section is oriented by an
  auxiliary spine built from our own rotation-minimizing frame (OCCT's CorrectedFrenet let the seam and edge
  lines wobble a few degrees around the tube); while dragging on screen only the pipe is swept (fast mode, no
  boolean / refine / orientation test, coarse preview mesh: ~0.05 s per step) and the smooth one is computed on
  release. The boolean result is checked by volume (fuse: body + pipe - common, cut: body - common;
  a sweep that looks clean could lose the pipe in the boolean), falling back to the next sweep mode; the slow
  BOP check (up to 25 s on sharp bends) is not used. Works without the GUI (freecadcmd).
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
  axis nearest to the screen's up / right direction. With polygons or edges selected their
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
