# SubD Modeler for FreeCAD

LightWave / Blender-style **subdivision surface (SubD) modeling inside FreeCAD's 3D view** — and a
**real B-rep solid** at the end that exports to STEP and works with Part / PartDesign.

![SubD Modeler user interface](docs/ui_overview.png)

You model a low-poly control cage with polygon tools (extrude, bevel, knife, loop cut...) and see the smooth
Catmull-Clark surface live. Everything has exact numeric input, so it stays a CAD tool: after every operation a
small panel in the viewport shows its values — type a number, press Enter, and the operation is redone with it.

> **Status:** early version (0.4.0, see the [changelog](CHANGELOG.md)). It works and is fun to use, but expect
> rough edges. Feedback and bug reports are very welcome.

## Is the result a mesh or a B-rep?

**A B-rep.** The modeler stores the cage in a `SubD` document object (its Catmull-Clark surface is shown as a
mesh while you work). **Solid (STEP)** turns it into a parametric `Part` solid made of bicubic B-spline faces:

- after one Catmull-Clark step the cage is all quads; every quad becomes one untrimmed B-spline patch
  (9 x 9 interpolated points by default),
- neighbouring patches share their boundary curves exactly, so the faces sew into a closed, valid solid,
- the solid is linked to the `SubD` object and updates when you edit the cage,
- creased edges stay sharp, so a fully creased 20 mm cube becomes exactly 20 x 20 x 20 mm (8000 mm³).

![B-rep solid made of B-spline patches](docs/step_solid.png)

On a 95 mm test model the B-rep deviates at most ~0.02 mm from the subdivision surface. Evaluating the exact
Catmull-Clark limit surface is on the roadmap.

## Features

**Cage modeling (in the 3D view)**
- Point / edge / polygon selection, box select, loop select, ghost mode (see and select hidden elements)
- Extrude faces or open edges (edge mode: the surface continues outward), with segments (mouse wheel) and an
  X / Y / Z axis lock in the panel or with the X / Y / Z keys while dragging (no axis = along the normal)
- Bevel / Smooth Shift and Inset (mouse, or exact Inset / Shift values)
- **Edge bevel** (LightWave Rounder / Blender Ctrl+B): width + segments (mouse wheel), true circular profile,
  clean all-quad corners where three bevelled edges meet, straight 45° miters on edge loops
- **Knife (loop)**: loop cut that slides along the ring, exact position in % or mm
- **Knife (line)**: LightWave knife — draw a line and the model is cut **live** by the plane through it;
  then set the plane exactly (normal + offset)
- **Add points** (LightWave Add Points): selected edges are split into equal parts (*Count* in the panel); with
  nothing selected the tool stays on and a click on an edge adds a point there (it snaps to the middle, exact
  *Position* % in the panel). Neighbouring polygons take the new points, creased edges stay creased
- **Connect** (Ctrl+L, LightWave Split): an edge between two selected points of one polygon, or between the
  middles of two selected edges of one polygon; the polygon is split in two
- **Subdivide** (Shift+D, LightWave Subdivide): the selected polygons are split into quads (n sides → n quads);
  neighbouring polygons take the edge middles, so the cage stays closed
- **Array** (Polygon > Array, LightWave Array): the selected polygons are copied as separate parts, either around an
  axis (*Radial*: count, start / end angle, offset along the axis for a spiral) or along the world axes
  (*Rectangular*: count and offset for X / Y / Z). Offset type: Between, Total or Multiply. The axis is a drawn line
  (*Draw axis*, with the same snapping as Mirror (line); a click gives an axis along the view) or numbers in the
  panel.
- **Twist** (Modify > Twist, LightWave Twist): the selected polygons are twisted around an axis; the angle grows along
  the axis by a falloff curve (Linear, Smooth, Ease in, Ease out, Peak, Dome). Move the mouse for the angle, click to
  accept, then edit the values in the panel. *Cuts* splits the selection along the axis first so the twist is smooth.
  The axis is a drawn line (snapping), Axis X / Y / Z or numbers, as in Array.
- **Distance** (top bar, key U): measures the distance between two points (cage points, edge middles, polygon centers,
  other objects or the surface) in Edit and Object mode. With *Keep (fixed)* the measurement stays, follows the cage
  points when they move and is deleted with the × next to its text; otherwise it is temporary.
  The **≡** button opens a list of the distances: names, target / tolerance, show / hide, delete, CSV export.
- **Shortcuts for Array and Twist**: none by default; assign one with `subd_modeler.set_user_key("twist", "Shift+T")`.
- Crease (sharp or semi-sharp edges) — makes the surface hit exact dimensions
- Copy / cut / paste polygons (with their creases): paste in place into the same cage, into another SubD
  object, or as a new object (e.g. copy the neck of a bottle to start its cap)
- **From face** (drop-down: Face / Edge / Point), from a Part / PartDesign object (in object mode FreeCAD's own
  selection works):
  - *Face*: select faces (or closed flat edge loops) and turn them into one cage: four-sided faces (flat or
    curved: cylinders, B-splines...) become quad grids, flat faces with one or two boundaries become fans / quad
    strips, flat faces with several holes get a quad ring around each hole and a quad-dominant fill between
    them, neighbouring faces are welded along their shared edges, sharp joins and corners stay creased. With
    *Exact Edge* the cage is solved so the SubD surface passes through all sample points (about 0.08 mm on a
    4 mm thick filleted box)
  - *Edge*: each chain of selected edges becomes one polygon (flat or not, e.g. a mating line). Edges whose ends
    meet (within 0.01 mm, also from different objects) form one chain; an open chain is closed with a straight
    edge; separate chains give separate polygons. Straight edges stay one polygon edge; curves are sampled by
    *Arc Segments* (per 90°) / *Max Length*. The polygon stays **linked to its source** (e.g. a Sketcher
    sketch); when the sketch changes size or position:
    - a cage that is still just that flat polygon is rebuilt on the new curve, with a point count to suit it;
    - in an edited cage (extruded, bevelled...) the linked points move to the same samples on the new curve and
      the points straight above / below them (extruded copies) move with them. Other points (bevel and inset
      rings) stay where they are: nothing is recomputed, the modeler has no history. A linked point you moved by
      hand leaves the link; when none is left, the link ends.

    Cancelling a sketch edit brings the cage back. Options > *Follow sketch* (or the object's *SketchLink*
    property) turns the link off.
  - *Point*: the selected points become one polygon, in the order they were picked.
- **Make polygon** (P, LightWave Make Polygon / Blender F): a polygon from 3 or more selected points, a closed
  edge loop (fills a hole), or two edges / an open chain (bridges their ends); faces its neighbours
- **Snap** (bottom bar): free moves (G, the gizmo's centre cube) snap the point nearest to the cursor to the
  vertices and edges of visible Part / PartDesign shapes
- **Hole**: a round hole in the middle of the selected faces (one, two, four... in one flat patch); the patch is
  re-meshed around it automatically (quads where the counts match), with radius, segments and a support loop
- **Bridge** (Blender Bridge Edge Loops): connects two open edge loops / chains, or two face regions (their
  faces are removed), with number of cuts, twist, smoothness (tangent to the surfaces) and profile factor
- Mirror X / Y / Z (LightWave Mirror): mirrors the selected polygons (or the whole cage) across a world axis
  plane (exact center in the panel) and welds the points on the plane, so a half model becomes one piece
- The tools are in the **SubD tools** panel, a FreeCAD panel you can dock or float; its LightWave-style **Create /
  Modify** tabs choose which tool menu it shows
- X symmetry for every tool
- Move / rotate / stretch / size gizmo (ported from FC YBT Raytracer), Ctrl = mm / degree steps
- Exact numeric input for every operation, plus selection center / size in mm

**Object mode (O)**
- Work with several SubD objects: click / Shift+click / box select them (or select them in the tree)
- Move and rotate with the gizmo (stored in FreeCAD's standard `Placement`), stretch / size the cage
- Exact position, rotation and size in the viewport panel
- Shell (Modify): whole object or selected polygons, inside / outside / centered, detail filter for a smooth inner side, fill / sharp / flat rims
- SubD options: edit / final level, open edge rule, smooth triangles, zebra stripes; solid precision presets and G0 / G1 / G2 continuity (regular patches exported exactly)
- New box with a Box panel (width / depth / height and X / Y / Z segments, rebuilt as you type), new cylinder or all-quad topo cylinder (diameter / height / segments), duplicate (D), mirror across X / Y / Z, delete, solids for all selected objects
- O (or Tab) on a selected object edits its cage; the other objects stay visible as reference

**CAD side**
- The cage is saved in the document (`SubD` object), auto-saved every 20 s
- Parametric B-rep solid (`SubD to solid`), STEP export, booleans, PartDesign
- Shows both the cage size and the real surface size
- **Rail Pipe / Rail Groove** (the buttons now come with the [Pimp My PartDesign](https://github.com/sedatkrn/partdesign_addins)
  add-in; documents that already use them keep working here): the exact B-rep counterpart of Rail bridge as a step in a
  PartDesign Body. Select one or two flat faces or closed sketches of the Body (a sketch: click its edge or pick it in the
  tree) and press Rail Pipe (adds) or Rail Groove (cuts): the start profile is swept along a guide curve to the
  second profile, or to a free end; in an empty Body the pipe is the first solid. The
  curve leaves (and enters) the faces perpendicular, passes through *GuidePoints* and follows the same rules as
  the modeler's guide (*Tension*, no bulge on uneven spacing); *EndScale* sizes a free end, *Reversed* flips
  the start direction. Double-click the feature to edit the guide on screen (shown like Pad / Pocket editing: the previous solid
  plus a green / red translucent preview): click a yellow point and move it with
  the modeler's Move gizmo (arrows: along X / Y / Z, squares: in a plane, cube: in the screen plane) or drag
  the point itself (the white end points belong to the profiles); click the curve to pick a span and *Add
  point* puts a point exactly there; delete points, *Move neighbors* makes the guide follow like a cable
  (as in the modeler); the task panel also has the selected point's X / Y / Z, *Tension*, *End scale* and
  *Reversed*. OK keeps the changes as one undo step, Cancel reverts them. The feature recomputes like any
  PartDesign feature and later steps (fillets, pockets) build on it.

![Edge bevel with 4 segments](docs/edge_bevel_smooth.png)

![Line knife with exact plane values](docs/knife_panel.png)

## Installation

Requires **FreeCAD 1.0 or newer** (developed on 26.3). No extra Python packages — only `numpy`, which ships
with FreeCAD.

**Addon Manager (custom repository)**
1. *Tools → Addon manager → ⚙ (settings) → Custom repositories*, add
   `https://github.com/sedatkrn/FcSubDModelerWB` with branch `main`.
2. Find **SubD Modeler** in the list, install, restart FreeCAD.

**Manual**
Clone or unzip this repository into FreeCAD's user `Mod` folder and restart FreeCAD:
- Windows: `%APPDATA%\FreeCAD\<version>\Mod\FcSubDModelerWB`
- Linux: `~/.local/share/FreeCAD/<version>/Mod/FcSubDModelerWB` (older versions: `~/.FreeCAD/Mod`)
- macOS: `~/Library/Application Support/FreeCAD/<version>/Mod/FcSubDModelerWB`

## Usage

Switch to the **SubD Modeler** workbench and press **SubD Modeler** (or run `macros/SubDModeler.FCMacro` from
any workbench). With one `SubD` object selected its cage is opened for editing; otherwise you start in object
mode, and an empty document stays empty until you create something (**Create > Box**, **Cylinder**, ...).
**O** switches between object mode and editing the selected object (**Tab** toggles the SubPatch preview, as
in LightWave). **Close** writes the cage back to the
document; save the document as usual (Ctrl+S).

| Key | Action | Key | Action |
|---|---|---|---|
| O | object mode / edit the selected object | D | duplicate (object mode) |
| 1 / 2 / 3 | point / edge / polygon mode | Tab / S | SubPatch preview on/off |
| Click, Shift+click | select, add/remove | Drag | box select (Shift add, Ctrl remove) |
| A | select all | L | loop from the selection (edges, points, two faces), or under the cursor |
| E | extrude (faces / open edges, wheel = segments) | B | bevel (polygon mode) / edge bevel (edge mode) |
| I | inset | C | crease / uncrease selected edges |
| K | loop cut + slide (stays on) | J | line knife, live plane cut (stays on) |
| T / R / H | move / rotate / stretch-size gizmo | G | free move with the mouse |
| Y | X symmetry | X | ghost mode |
| Q | Snap on/off | | |
| Del | delete faces | Ctrl+Z / Ctrl+Y | undo / redo one step (object mode: one whole edit) |
| P | make polygon from selected points / edges | M | merge points (selection or whole cage) |
| Ctrl+W | weld selected points at the last one clicked | Shift+Z | merge selected adjacent polygons |
| Ctrl+C / Ctrl+X | copy / cut polygons | Ctrl+V | paste in place (object mode: as a new SubD object) |
| F | save to the document | N | Options panel |
| Ctrl (while dragging) | mm / degree steps | Shift (while dragging) | fine control |
| Esc / right click | cancel the running tool / turn the knife off | Mouse wheel (edge bevel) | segments |

## Language

The interface is English by default and **Turkish** is built in. The **EN / TR** button at the right end of the
bottom bar switches it (the modeler reopens with the same model). The setting is stored in
`Preferences/Mod/SubDModeler/Language` (`auto` follows FreeCAD's language).

![Turkish interface](docs/turkish_ui.png)

## Roadmap

- Exact Catmull-Clark limit surface for the B-rep (Stam evaluation)
- Edge slide, more selection tools (face strips, grow / shrink)
- Export / import of the raw cage with creases (OBJ / USD) for Blender round-trips
- Odd segment counts for all-quad edge-bevel corners
- From face: curved faces that are not four-sided
- Rail Pipe: scale per guide point, sketches with several profiles (holes)

## License

MIT — see [LICENSE](LICENSE). Code comments are mostly in Turkish (the author's language); the user
interface and documentation are in English.

---

### Türkçe

FreeCAD için LightWave / Blender tarzı SubD (alt bölümlemeli yüzey) modelleyici: 3B görünümün içinde kafes
modelleme ve sonunda STEP'e aktarılabilen gerçek bir B-rep katı. Arayüz dili alt çubuktaki **EN / TR**
düğmesiyle Türkçeye çevrilebilir. Kurulum ve kullanım yukarıdaki İngilizce bölümle aynıdır.
