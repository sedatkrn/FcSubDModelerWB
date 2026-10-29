# -*- coding: utf-8 -*-
# SubD Modeler - Rail Pipe / Rail Groove: kilavuz egriyi ekranda duzenleme. Copyright (c) 2026 Sedat Kuran - MIT
# License
#
# RailPipe'a cift tiklayinca (setEdit): kilavuz egri sari cizgi, noktalari tutamak (A / B uclari sabit, beyaz).
# Secili noktada modelleyicinin Tasi gizmosu (oklar: eksen, kareler: duzlem, kup: ekran duzlemi); noktanin kendisi
# de surukulenebilir (ekran duzleminde). Egri hemen, boru kisa araliklarla yeniden hesaplanir. Gorev paneli: secili noktanin X / Y / Z'si, nokta ekle / sil, Move
# neighbors (modelleyicideki kablo davranisi, subd_modeler.rail_cable) ve Elasticity, Tension, EndScale, Reversed.
# OK degisiklikleri tek geri alma adimi olarak birakir, Cancel geri alir.
import numpy as np

import FreeCAD as App
import FreeCADGui as Gui
from pivy import coin
from PySide import QtCore, QtWidgets

import subd_railpipe as rp
from subd_i18n import T

V = App.Vector
_editor = None


def start(obj):
    global _editor
    stop()
    _editor = Editor(obj)


def stop():
    global _editor
    if _editor is not None:
        e, _editor = _editor, None
        e.close()


def _points(color, size, coord, node):
    lm = coin.SoLightModel()
    lm.model = coin.SoLightModel.BASE_COLOR
    bc = coin.SoBaseColor()
    bc.rgb = color
    ds = coin.SoDrawStyle()
    ds.pointSize = size
    ds.lineWidth = size
    sep = coin.SoSeparator()
    for n in (lm, bc, ds, coord, node):
        sep.addChild(n)
    return sep


def _set(coord, node, P, points=True):
    n = len(P)
    P = [tuple(p) for p in P] or [(0.0, 0.0, 0.0)]
    coord.point.setValues(0, len(P), P)
    coord.point.setNum(len(P))
    if points:
        node.numPoints = n
    elif n:
        node.numVertices.setValues(0, 1, [n])
    else:
        node.numVertices.setNum(0)


class Editor(object):
    def __init__(self, obj):
        import subd_modeler
        self.obj, self.doc = obj, obj.Document
        self.view = Gui.ActiveDocument.ActiveView
        self.gp = obj.getGlobalPlacement() if hasattr(obj, "getGlobalPlacement") else App.Placement()
        self.inv = self.gp.inverse()
        self.sel = None                                    # secili dugum (dugum listesinde sira)
        self.span = None                                   # secili egri parcasi: (parca sirasi, tiklanan parametre)
        self.drag = None
        self.cable = False
        self.elastic = 0.1
        self.doc.openTransaction(T("Edit Rail Pipe"))
        # cizim: egri, sabit uclar, tasinabilir noktalar, secili nokta (derinlik testi olmadan: govdenin icinde de)
        self.c_coord, self.c_line = coin.SoCoordinate3(), coin.SoLineSet()
        self.f_coord, self.f_pts = coin.SoCoordinate3(), coin.SoPointSet()
        self.m_coord, self.m_pts = coin.SoCoordinate3(), coin.SoPointSet()
        self.s_coord, self.s_pts = coin.SoCoordinate3(), coin.SoPointSet()
        self.p_coord, self.p_line = coin.SoCoordinate3(), coin.SoLineSet()
        self.root = coin.SoAnnotation()
        for n in (_points((0.98, 0.80, 0.15), 2.5, self.c_coord, self.c_line),
                  _points((1.0, 0.45, 0.0), 5, self.p_coord, self.p_line),
                  _points((0.95, 0.95, 0.95), 10, self.f_coord, self.f_pts),
                  _points((0.98, 0.80, 0.15), 10, self.m_coord, self.m_pts),
                  _points((1.0, 0.30, 0.0), 14, self.s_coord, self.s_pts)):
            self.root.addChild(n)
        self.parent = subd_modeler._scene_parent(self.view)
        self.parent.addChild(self.root)
        self.gizmo = subd_modeler.MoveGizmo(self.view)     # modelleyicideki Tasi gizmosu: oklar, duzlemler, kup
        self.parent.addChild(self.gizmo.root)
        self.zoom = QtCore.QTimer()                        # yakinlastirinca gizmo ekranda ayni boyda kalsin
        self.zoom.timeout.connect(self.gizmo.update_scale)
        self.zoom.start(50)
        self.transparency = obj.ViewObject.Transparency    # duzenlerken boru yari saydam: egri ve noktalar gorunsun
        obj.ViewObject.Transparency = max(self.transparency, 60)
        self.cb_btn = self.view.addEventCallbackPivy(coin.SoMouseButtonEvent.getClassTypeId(), self.on_btn)
        self.cb_mov = self.view.addEventCallbackPivy(coin.SoLocation2Event.getClassTypeId(), self.on_move)
        self.timer = QtCore.QTimer()                       # surukleme sirasinda boru en fazla bu aralikla hesaplanir
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self._recompute)
        self.panel = Panel(self)
        Gui.Control.closeDialog()
        Gui.Control.showDialog(self.panel)
        self.redraw()

    # --- koordinatlar ---
    def world(self, p):
        return self.gp.multVec(V(p))

    def local_dir(self, d):
        return self.inv.Rotation.multVec(V(d))

    def nodes(self):
        res = rp.guide_nodes(self.obj)
        return None if isinstance(res, str) else res

    def movable(self, j, N):
        return 0 < j < (len(N) if self.obj.EndFace is None or not self.obj.EndFace[0] else len(N) - 1)

    def free_end(self):
        return self.obj.EndFace is None or not self.obj.EndFace[0]

    # --- cizim ---
    def redraw(self):
        res = self.nodes()
        if res is None:
            for c, n in ((self.c_coord, self.c_line), (self.f_coord, self.f_pts), (self.m_coord, self.m_pts),
                         (self.s_coord, self.s_pts)):
                _set(c, n, [], n is not self.c_line)
            return
        N, ta, tb = res
        bs = self.bs = rp.guide_spline(N, ta, tb, self.obj.Tension)
        a, b = bs.FirstParameter, bs.LastParameter
        self.curve = [(a + (b - a) * k / 200.0,) for k in range(201)]
        self.curve = [(u, self.world(bs.value(u))) for (u,) in self.curve]
        _set(self.c_coord, self.c_line, [w for u, w in self.curve], False)
        K = bs.getKnots()                                  # K[i] -> dugum i (uzantisiz egri)
        if self.span is not None and self.span[0] < len(K) - 1:
            i = self.span[0]
            _set(self.p_coord, self.p_line, [self.world(bs.value(K[i] + (K[i + 1] - K[i]) * k / 40.0))
                                             for k in range(41)], False)
        else:
            _set(self.p_coord, self.p_line, [], False)
        W = [self.world(p) for p in N]
        _set(self.f_coord, self.f_pts, [w for j, w in enumerate(W) if not self.movable(j, N)])
        _set(self.m_coord, self.m_pts, [w for j, w in enumerate(W) if self.movable(j, N)])
        _set(self.s_coord, self.s_pts, [W[self.sel]] if self.sel is not None and self.sel < len(W) else [])
        on = self.sel is not None and self.sel < len(W) and self.movable(self.sel, N)
        if self.drag is None or self.drag["key"][0] == "point":
            self.gizmo.place(W[self.sel] if on else None)

    def _recompute(self):
        rp.PREVIEW = True                                  # surukleme onizlemesi: hizli supurme
        try:
            self.obj.recompute()
        except Exception:
            pass
        finally:
            rp.PREVIEW = False
        Gui.updateGui()

    # --- fare ---
    def pick(self, pos, N, radius=12.0):
        best = None
        for j, p in enumerate(N):
            s = self.view.getPointOnViewport(self.world(p))
            d = ((s[0] - pos[0]) ** 2 + (s[1] - pos[1]) ** 2) ** 0.5
            if d < radius and (best is None or d < best[0]):
                best = (d, j)
        return None if best is None else best[1]

    def pick_span(self, pos, radius=8.0):
        """Imlecin altindaki egri parcasi: (parca sirasi, en yakin noktanin parametresi) ya da None."""
        if not getattr(self, "curve", None):
            return None
        S = [self.view.getPointOnViewport(w) for u, w in self.curve]
        best = None
        for k in range(len(S) - 1):
            ax, ay = S[k][0], S[k][1]
            dx, dy = S[k + 1][0] - ax, S[k + 1][1] - ay
            L2 = dx * dx + dy * dy
            t = 0.0 if L2 < 1e-9 else max(0.0, min(1.0, ((pos[0] - ax) * dx + (pos[1] - ay) * dy) / L2))
            d = ((pos[0] - ax - t * dx) ** 2 + (pos[1] - ay - t * dy) ** 2) ** 0.5
            if d < radius and (best is None or d < best[0]):
                best = (d, self.curve[k][0] + t * (self.curve[k + 1][0] - self.curve[k][0]))
        if best is None:
            return None
        K = self.bs.getKnots()
        i = max(k for k in range(len(K) - 1) if K[k] <= best[1] + 1e-12) if best[1] > K[0] else 0
        return min(i, len(K) - 2), best[1]

    def ray_plane(self, pos, p0):
        """Imlecten isinin p0'dan gecen, bakisa dik duzlemle kesisimi (dunya)."""
        p = V(self.view.getPoint(int(pos[0]), int(pos[1])))
        n = V(self.view.getViewDirection())
        d = V(n)
        if self.view.getCameraType() == "Perspective":
            o = V(*self.view.getCameraNode().position.getValue().getValue())
            d = p - o
        den = d.dot(n)
        if abs(den) < 1e-12:
            return p0
        return p + d * ((p0 - p).dot(n) / den)

    def on_btn(self, cb):
        ev = cb.getEvent()
        if ev.getButton() != coin.SoMouseButtonEvent.BUTTON1:
            return
        pos = ev.getPosition().getValue()
        down = ev.getState() == coin.SoButtonEvent.DOWN
        if not down:
            if self.drag is not None:
                cb.setHandled()
                self.end_drag()
            return
        if self.drag is not None:                          # birakma kacirildiysa yeni tiklama once onu bitirir
            self.end_drag()
        res = self.nodes()
        if res is None:
            return
        N = res[0]
        # once baska bir nokta (gizmonun oklari / kareleri yakindaki noktalarin ustune biniyordu: tiklama noktayi
        # secmek yerine gizmoyu tutuyordu), sonra secili noktadaki gizmo
        j = self.pick(pos, N)
        key = None
        if j is None or j == self.sel:
            key = self.gizmo.pick(pos)
            if key is not None:
                j = self.sel
        if j is None:
            sp = self.pick_span(pos)
            if sp is None:
                return                                     # bosluk: gezinme / secim FreeCAD'e kalsin
            cb.setHandled()
            self.sel, self.span = None, sp
            self.redraw()
            self.panel.sync()
            return
        cb.setHandled()
        self.sel, self.span = j, None
        if self.movable(j, N):
            w0 = self.world(N[j])
            self.drag = dict(key=key or ("point", 0), j=j, N=[V(p) for p in N],
                             pts=[V(p) for p in self.obj.GuidePoints], w0=w0, ta=res[1], tb=res[2])
            self.drag["ref"] = self.measure(pos)
        self.redraw()
        self.panel.sync()

    def measure(self, pos):
        """Surukleme turune gore olcum (dunya): ok -> eksen parametresi, duzlem karesi -> duzlemdeki nokta, kup /
        noktanin kendisi -> bakisa dik duzlemdeki nokta (modelleyicideki drag_measure gibi)."""
        import subd_modeler as sm
        g = self.drag
        kind, i = g["key"]
        o, d = sm._ray(self.view, pos)
        if kind == "move":
            return sm._line_param(o, d, g["w0"], sm.AXES[i])
        if kind == "plane":
            return sm._plane_hit(o, d, g["w0"], sm.AXES[i])
        return self.ray_plane(pos, g["w0"])

    def end_drag(self):
        self.drag = None
        self.timer.stop()
        self.obj.touch()                                   # onizleme hizli yontemle hesaplanmisti: kaliteli hal
        self.doc.recompute()
        self.redraw()
        self.panel.sync()

    def on_move(self, cb):
        import subd_modeler as sm
        g = self.drag
        pos = cb.getEvent().getPosition().getValue()
        if g is None:
            self.gizmo.set_hover(self.gizmo.pick(pos))
            return
        if not sm._left_held():                            # birakma olayi gelmedi (gorunum disinda / FreeCAD yuttu):
            self.end_drag()                                # nokta fareye yapisik kalmasin
            return
        cb.setHandled()
        cur = self.measure(pos)
        if cur is None or g["ref"] is None:
            return
        kind, i = g["key"]
        dw = sm.AXES[i] * (cur - g["ref"]) if kind == "move" else cur - g["ref"]
        self.gizmo.place(g["w0"] + dw)
        self.move_to(g, self.local_dir(dw))

    def move_to(self, g, d):
        """Surukleme baslangicina gore d (yerel) tasima: yalniz nokta ya da Move neighbors ile kablo gibi."""
        pts = [V(p) for p in g["pts"]]
        j = g["j"]
        if self.cable:
            import subd_modeler
            N = np.array([tuple(p) for p in g["N"]], float)
            tb = None if g["tb"] is None else np.array(tuple(g["tb"]), float)
            X = subd_modeler.rail_cable(N, j, np.array(tuple(d), float), not self.free_end(), self.elastic,
                                        np.array(tuple(g["ta"]), float), tb)
            pts = [V(*X[k + 1]) for k in range(len(pts))]
        else:
            pts[j - 1] = pts[j - 1] + d
        self.obj.GuidePoints = pts
        self.redraw()
        self.panel.sync()
        if not self.timer.isActive():
            self.timer.start(120)

    # --- panel islemleri ---
    def set_coord(self, axis, val):
        res = self.nodes()
        if res is None or self.sel is None or not self.movable(self.sel, res[0]):
            return
        N, ta, tb = res
        w = self.world(N[self.sel])
        t = V(w)
        setattr(t, "xyz"[axis], val)
        if (t - w).Length < 1e-9:
            return
        g = dict(j=self.sel, N=[V(p) for p in N], pts=[V(p) for p in self.obj.GuidePoints], w0=w, ta=ta, tb=tb)
        self.move_to(g, self.local_dir(t - w))
        self.doc.recompute()

    def add_point(self):
        """Secili egri parcasinda tiklanan yere; yoksa secili noktadan sonraki parcanin (secim yoksa en uzun
        parcanin) ortasina yeni nokta."""
        res = self.nodes()
        if res is None:
            return
        N, ta, tb = res
        n = len(N)
        if self.sel is not None:
            i = min(self.sel, n - 2)
        else:
            i = max(range(n - 1), key=lambda k: (N[k + 1] - N[k]).Length)
        bs = rp.guide_spline(N, ta, tb, self.obj.Tension)
        K = bs.getKnots()
        u = (K[i] + K[i + 1]) / 2
        if self.span is not None and self.span[0] < len(K) - 1:   # secili parca: tiklanan yere
            i, u = self.span
            u = min(max(u, K[i] + 1e-6 * (K[i + 1] - K[i])), K[i + 1] - 1e-6 * (K[i + 1] - K[i]))
        mid = bs.value(u)
        pts = list(self.obj.GuidePoints)
        pts.insert(i, mid)                                 # dugum i + 1 olur
        self.obj.GuidePoints = pts
        self.sel, self.span = i + 1, None
        self.doc.recompute()
        self.redraw()
        self.panel.sync()

    def delete_point(self):
        res = self.nodes()
        if res is None or self.sel is None or not self.movable(self.sel, res[0]):
            return
        pts = list(self.obj.GuidePoints)
        if self.free_end() and len(pts) <= 1:
            return                                         # serbest uc tek noktadan az olamaz
        del pts[self.sel - 1]
        self.obj.GuidePoints = pts
        self.sel = None
        self.doc.recompute()
        self.redraw()
        self.panel.sync()

    def set_prop(self, name, val):
        if getattr(self.obj, name) != val:
            setattr(self.obj, name, val)
            self.doc.recompute()
            self.redraw()

    # --- bitis ---
    def finish(self, ok):
        self.finished = True
        self.close()
        if ok:
            self.doc.commitTransaction()
        else:
            self.doc.abortTransaction()
        self.doc.recompute()
        global _editor
        _editor = None
        Gui.ActiveDocument.resetEdit()

    def close(self):
        if getattr(self, "closed", False):
            return
        self.closed = True
        self.timer.stop()
        if not getattr(self, "finished", False):            # duzenleme baska yoldan bitti: degisiklikler kalsin
            self.doc.commitTransaction()
        try:
            self.view.removeEventCallbackPivy(coin.SoMouseButtonEvent.getClassTypeId(), self.cb_btn)
            self.view.removeEventCallbackPivy(coin.SoLocation2Event.getClassTypeId(), self.cb_mov)
        except Exception:
            pass
        self.zoom.stop()
        try:
            self.parent.removeChild(self.gizmo.root)
        except Exception:
            pass
        try:
            self.parent.removeChild(self.root)
            self.obj.ViewObject.Transparency = self.transparency
        except Exception:
            pass
        Gui.Control.closeDialog()


class Panel(object):
    """Gorev paneli: secili nokta X / Y / Z, nokta ekle / sil, Move neighbors + Elasticity, Tension, EndScale,
    Reversed."""

    def __init__(self, ed):
        self.ed = ed
        self.form = w = QtWidgets.QWidget()
        w.setWindowTitle(T("Rail Pipe"))
        lay = QtWidgets.QVBoxLayout(w)
        tip = QtWidgets.QLabel(T("Click a guide point in the 3D view and move it with the gizmo (arrows: along an "
                                 "axis, squares: in a plane, cube: free) or drag the point itself. The white end "
                                 "points belong to the profiles and stay put. Click the curve to pick a span: Add "
                                 "point puts the new point where you clicked."))
        tip.setWordWrap(True)
        lay.addWidget(tip)
        g = QtWidgets.QGroupBox(T("Selected point"))
        grid = QtWidgets.QGridLayout(g)
        self.xyz = []
        for k, a in enumerate("XYZ"):
            grid.addWidget(QtWidgets.QLabel(a), k, 0)
            b = self._spin(-1e7, 1e7, 3, " mm")
            b.editingFinished.connect(lambda k=k, b=b: self.ed.set_coord(k, b.value()))
            grid.addWidget(b, k, 1)
            self.xyz.append(b)
        row = QtWidgets.QHBoxLayout()
        for text, fn in ((T("Add point"), ed.add_point), (T("Delete point"), ed.delete_point)):
            b = QtWidgets.QPushButton(text)
            b.clicked.connect(lambda checked=False, fn=fn: fn())
            row.addWidget(b)
        grid.addLayout(row, 3, 0, 1, 2)
        lay.addWidget(g)
        g = QtWidgets.QGroupBox(T("Guide"))
        grid = QtWidgets.QGridLayout(g)
        self.cable = QtWidgets.QCheckBox(T("Move neighbors"))
        self.cable.setToolTip(T("The guide follows the dragged point like a cable"))
        self.cable.toggled.connect(lambda on: setattr(self.ed, "cable", bool(on)))
        grid.addWidget(self.cable, 0, 0, 1, 2)
        grid.addWidget(QtWidgets.QLabel(T("Elasticity")), 1, 0)
        self.elastic = self._spin(0.0, 1.0, 2, "", 0.05)
        self.elastic.setValue(ed.elastic)
        self.elastic.valueChanged.connect(lambda v: setattr(self.ed, "elastic", float(v)))
        grid.addWidget(self.elastic, 1, 1)
        grid.addWidget(QtWidgets.QLabel(T("Tension")), 2, 0)
        self.tension = self._spin(0.0, 1.0, 2, "", 0.05)
        self.tension.setValue(ed.obj.Tension)
        self.tension.valueChanged.connect(lambda v: self.ed.set_prop("Tension", float(v)))
        grid.addWidget(self.tension, 2, 1)
        grid.addWidget(QtWidgets.QLabel(T("End scale")), 3, 0)
        self.scale = self._spin(0.01, 100.0, 3, " x", 0.05)
        self.scale.setValue(ed.obj.EndScale)
        self.scale.setEnabled(ed.free_end())
        self.scale.valueChanged.connect(lambda v: self.ed.set_prop("EndScale", float(v)))
        grid.addWidget(self.scale, 3, 1)
        self.rev = QtWidgets.QCheckBox(T("Reversed"))
        self.rev.setChecked(bool(ed.obj.Reversed))
        self.rev.toggled.connect(lambda on: self.ed.set_prop("Reversed", bool(on)))
        grid.addWidget(self.rev, 4, 0, 1, 2)
        lay.addWidget(g)
        lay.addStretch()
        self.sync()

    @staticmethod
    def _spin(lo, hi, dec, suffix, step=1.0):
        b = QtWidgets.QDoubleSpinBox()
        b.setRange(lo, hi)
        b.setDecimals(dec)
        b.setSuffix(suffix)
        b.setSingleStep(step)
        b.setKeyboardTracking(False)
        return b

    def sync(self):
        res = self.ed.nodes()
        j = self.ed.sel
        ok = res is not None and j is not None and j < len(res[0]) and self.ed.movable(j, res[0])
        for k, b in enumerate(self.xyz):
            b.setEnabled(ok)
            if ok and not b.hasFocus():
                b.blockSignals(True)
                b.setValue(getattr(self.ed.world(res[0][j]), "xyz"[k]))
                b.blockSignals(False)

    def getStandardButtons(self):
        return int(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)

    def accept(self):
        self.ed.finish(True)
        return True

    def reject(self):
        self.ed.finish(False)
        return True
