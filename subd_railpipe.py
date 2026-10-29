# -*- coding: utf-8 -*-
# SubD Modeler - Rail Pipe / Rail Groove: PartDesign features (Body icinde Pad / Pocket gibi bir adim). Copyright (c)
# 2026 Sedat Kuran - MIT License
#
# Modelleyicideki Rail bridge'in kesin B-rep karsiligi. Baslangic yuzu (A, duz) ile bitis yuzu (B) arasinda ya da A'dan
# serbest uca, kilavuz egri boyunca OCCT MakePipeShell supurmesi: A'nin dis teli A'dan B'nin teline (tek uclu:
# EndScale ile olceklenmis A teline) kademeli gecer. Egri modelleyicideki gibi Hermite parcalaridir (rail_nodes /
# rail_seg ile ayni kurallar): A'dan (ve B'ye) yuzlere dik, ara noktalarda Catmull-Rom tegetleri, araliklar esit
# degilse parca boyuyla olceklenir, Tension tegetleri kisaltir. Tek bir kubik B-spline olarak kurulur (Bezier
# parcalari, dugumlerde G1). Rail Pipe Body'ye ekler (FeatureAdditivePython), Rail Groove oyar
# (FeatureSubtractivePython; egri yuzden iceri girer). GUI gerektirmez: belge freecadcmd ile de yeniden hesaplanir.
import math

import FreeCAD as App
import Part

from subd_i18n import T

V = App.Vector
PROPS = "GuidePoints"     # bu ozellik varsa nesne Rail Pipe / Groove'dur (onarim, komutlar)


def _face(link):
    """PropertyLinkSub (nesne, [alt eleman]) -> Part yuzu ya da None."""
    if not link or not link[1]:
        return None
    try:
        sh = link[0].Shape.getElement(link[1][0])
    except Exception:
        return None
    return sh if sh.ShapeType == "Face" else None


def face_frame(f):
    """Yuzun merkezi ve (katinin disina bakan) normali."""
    c = f.CenterOfMass
    u, v = f.Surface.parameter(c)
    return c, f.normalAt(u, v)


def _hermite(P0, P1, m0, m1, t):
    t2, t3 = t * t, t * t * t
    return P0 * (2 * t3 - 3 * t2 + 1) + m0 * (t3 - 2 * t2 + t) + P1 * (-2 * t3 + 3 * t2) + m1 * (t3 - t2)


def guide_spline(N, ta, tb=None, tension=0.0, ext_a=0.0, ext_b=0.0):
    """Dugumlerden (V listesi, ilki A merkezi; tb varsa sonuncusu B merkezi) kubik B-spline. ta: A'dan cikis, tb: B'ye
    giris yonu (birim); tb None: serbest uc, son parcanin yonu. ext_a / ext_b > 0: egri A'nin gerisinden (B'nin
    otesinden) o kadar duz bir parcayla baslar (biter); dugum araliklari boylarla orantili oldugundan C1 baglanir."""
    n = len(N)
    L = [(N[i + 1] - N[i]).Length for i in range(n - 1)]
    M = [V() for _ in N]
    for i in range(1, n - 1):
        M[i] = (N[i + 1] - N[i - 1]) * 0.5
    M[0] = ta * L[0]
    M[-1] = tb * L[-1] if tb is not None else N[-1] - N[-2]
    k = max(1.0 - min(max(float(tension), 0.0), 1.0), 0.02)
    poles = [N[0]]
    for i in range(n - 1):
        m0, m1 = M[i] * k, M[i + 1] * k
        if i > 0:                                          # esit olmayan araliklar: parcanin kendi boyuyla (gobek yok)
            m0 = m0 * (2 * L[i] / max(L[i - 1] + L[i], 1e-12))
        if i + 1 < n - 1:
            m1 = m1 * (2 * L[i] / max(L[i] + L[i + 1], 1e-12))
        poles += [N[i] + m0 * (1 / 3.0), N[i + 1] - m1 * (1 / 3.0), N[i + 1]]   # Hermite = Bezier
    # dugum araliklari parca boylariyla orantili: boyla olceklenen tegetler o zaman parametrede de C1 olur; katlar
    # 3'ten 2'ye indirilir (sekil ayni). C0 kalsaydi supurulen yuzeyler dugumlerde kirik sayilir, OCCT boolean'lari
    # (Body'ye ekleme / oyma) gecersiz sonuc verebiliyordu.
    knots = [0.0]
    for l in L:
        knots.append(knots[-1] + l)
    if ext_a > 0:                                          # duz on parca: tegeti A'daki ile ayni (k ile), C1
        P = N[0] - ta * ext_a
        poles = [P, P + ta * (ext_a * k / 3.0), N[0] - ta * (ext_a * k / 3.0)] + poles
        knots = [0.0] + [x + ext_a for x in knots]
    if ext_b > 0 and tb is not None:
        P = N[-1] + tb * ext_b
        poles += [N[-1] + tb * (ext_b * k / 3.0), P - tb * (ext_b * k / 3.0), P]
        knots.append(knots[-1] + ext_b)
    m = len(knots)
    bs = Part.BSplineCurve()
    bs.buildFromPolesMultsKnots(poles, [4] + [3] * (m - 2) + [4], knots, False, 3)
    for i in range(2, m):
        bs.removeKnot(i, 2, 1e-6 * max(knots[-1], 1.0))
    return bs


def _transport(bs, t0, steps=64):
    """Egri boyunca burulmasiz tasima: baslangic yonu t0'dan sona kadar en kisa donmelerin carpimi (App.Rotation)."""
    a, b = bs.FirstParameter, bs.LastParameter
    R, prev = App.Rotation(), V(t0)
    for k in range(1, steps + 1):
        t = bs.tangent(a + (b - a) * k / steps)[0]
        R = App.Rotation(prev, t).multiply(R)
        prev = t
    return R


def _plane_normal(P, along, tol=1e-6):
    """Noktalar bir duzlemdeyse o duzlemin normali (dogru uzerindeyse along'a dik herhangi bir yon), degilse None."""
    o = P[0]
    size = max((p - o).Length for p in P) or 1.0
    n = V()
    for p in P[1:]:                                        # en buyuk capraz carpim: en iyi kosullu normal
        for q in P[1:]:
            c = (p - o).cross(q - o)
            if c.Length > n.Length:
                n = c
    if n.Length < 1e-9 * size * size:                     # hepsi bir dogru uzerinde
        n = along.cross(V(1, 0, 0)) if abs(along.x) < 0.9 else along.cross(V(0, 1, 0))
        n.normalize()
        return n
    n.normalize()
    return n if all(abs((p - o).dot(n)) <= tol * size for p in P) else None


def _oriented(w, along):
    """Telin donus yonu (sag el kurali) along ile ayni olsun: iki uc teli ayni yonde donmezse supurme burulur."""
    try:
        n = Part.Face(w).normalAt(0, 0)
    except Exception:
        return w
    if n.dot(along) < 0:
        w = w.copy()
        w.reverse()
    return w


def _start_near(w, p):
    """Teli, baslangic noktasi p'ye en yakin kenardan baslayacak sekilde yeniden diz (supurmede koseler eslesir)."""
    E = w.OrderedEdges
    if len(E) < 2:
        return w
    k = min(range(len(E)), key=lambda i: (E[i].Vertexes[0].Point - p).Length)
    return Part.Wire(E[k:] + E[:k])


def build_tool(obj):
    """Supurulen kati (Body koordinatinda) ya da hata metni."""
    fa = _face(obj.StartFace)
    if fa is None:
        return "Rail Pipe: the start face is missing"
    if not (fa.isPlanar() if hasattr(fa, "isPlanar") else fa.Surface.TypeId == "Part::GeomPlane"):
        return "Rail Pipe: the start face must be flat"
    sub = obj.TypeId == "PartDesign::FeatureSubtractivePython"
    out = -1.0 if (sub != bool(obj.Reversed)) else 1.0    # Pipe disari, Groove iceri (Reversed cevirir)
    ca, na = face_frame(fa)
    ta = na * out
    fb = _face(obj.EndFace)
    pts = list(obj.GuidePoints)
    if fb is None and not pts:
        return "Rail Pipe: a free end needs at least one guide point"
    tb = None
    N = [ca] + pts
    if fb is not None:
        cb, nb = face_frame(fb)
        tb = nb * -out                                     # B'ye disaridan (Groove: iceriden) girer
        N.append(cb)
    for i in range(len(N) - 1):
        if (N[i + 1] - N[i]).Length < 1e-6:
            return "Rail Pipe: two guide points are on top of each other"
    # uclar yuzun obur yanindan baslar / biter (Pipe: govdenin icinden, Groove: disindan): supurmenin yaklasik ucu
    # yuzle tam cakisinca boolean cok kucuk kenar / kendini kesen kose birakiyor, sonraki adimlar bozuluyordu
    ea = min(0.5, 0.05 * fa.BoundBox.DiagonalLength)
    eb = min(0.5, 0.05 * fb.BoundBox.DiagonalLength) if fb is not None else 0.0
    bs = guide_spline(N, ta, tb, obj.Tension, ea, eb)
    spine = Part.Wire(bs.toShape())
    wa = _oriented(fa.OuterWire, ta).copy()
    wa.translate(-ta * ea)
    ca = ca - ta * ea
    if fb is not None:
        wb = _oriented(fb.OuterWire, tb).copy()
        wb.translate(tb * eb)
    else:                                                  # serbest uc: A'nin teli egri boyunca tasinip olceklenir
        R = _transport(bs, ta)
        wb = wa.copy()
        wb.translate(-ca)
        wb.rotate(V(), R.Axis, math.degrees(R.Angle))
        wb.scale(max(obj.EndScale, 0.01))
        wb.translate(N[-1])
    R = _transport(bs, ta)                                 # A'nin ilk kosesi sonda nereye gelir: B teli oradan baslar
    first = wa.OrderedEdges[0].Vertexes[0].Point
    end = N[-1] + (tb * eb if fb is not None else V())
    wb = _start_near(wb, R.multVec(first - ca) * (obj.EndScale if fb is None else 1.0) + end)
    single = fb is None and abs(obj.EndScale - 1.0) < 1e-9   # olceksiz serbest uc: tek profil (en saglami)
    box = spine.BoundBox
    box.enlarge(2.0 * max(wa.BoundBox.DiagonalLength, wb.BoundBox.DiagonalLength))
    mid = bs.value((bs.FirstParameter + bs.LastParameter) / 2)
    binormal = _plane_normal([ca] + N[1:] + [end], ta)
    err = "Rail Pipe: the sweep gave no solid"
    # sirayla denenir, gecerli ve egrinin cevresinde kalan ilk sonuc alinir: CorrectedFrenet (burulmasiz; egrilige
    # dayanir, duz cizgide tanimsiz kalip 'No geometry' veriyordu), sonra duzlemsel egride sabit binormal (duz
    # cizgide de calisir, ama iki profilli kivrimda bazen gecersiz). Her birinde once C1 yaklasim (C0 kenar / yuz
    # kalirsa sonraki boolean'lar bozuluyor), sonra C1'siz (iki profilde C1 yaklasim yuzeyi savurabiliyor).
    modes = [None] + ([binormal] if binormal is not None else [])
    for bn, c1 in [(m, c) for m in modes for c in (True, False)]:
        ps = Part.BRepOffsetAPI.MakePipeShell(spine)
        if bn is not None:
            ps.setBiNormalMode(bn)
        else:
            ps.setFrenetMode(False)
        ps.setForceApproxC1(c1)
        ps.add(wa, False, False)
        if not single:
            ps.add(wb, False, False)
        try:
            ps.build()
            ps.makeSolid()
            tool = ps.shape()
        except Exception as e:
            err = T("Rail Pipe: the sweep failed (%s)") % e
            continue
        if tool.isNull() or not tool.Solids:
            continue
        tool = tool.Solids[0]
        if not tool.isInside(mid, 1e-7, True):            # ters yuzlu (ici disina donuk) cikabiliyor: cevir
            tool.reverse()
        if tool.isValid() and box.isInside(tool.BoundBox):
            return tool
        err = "Rail Pipe: the sweep gave an invalid solid"
    return err


class RailPipe:
    """Rail Pipe (FeatureAdditivePython) / Rail Groove (FeatureSubtractivePython) nesnesinin Proxy'si."""

    def __init__(self, obj):
        obj.addProperty("App::PropertyLinkSub", "StartFace", "Rail", "Flat face the pipe starts from (A)")
        obj.addProperty("App::PropertyLinkSub", "EndFace", "Rail", "Face the pipe ends on (B); empty: free end")
        obj.addProperty("App::PropertyVectorList", PROPS, "Rail",
                        "Guide curve points between A and B (with a free end: the last one is the end)")
        obj.addProperty("App::PropertyFloat", "Tension", "Rail", "0 smooth curve .. 1 straight between points")
        obj.addProperty("App::PropertyFloat", "EndScale", "Rail", "Size of the free end relative to A").EndScale = 1.0
        obj.addProperty("App::PropertyBool", "Reversed", "Rail", "Leave A to the other side")
        if "Refine" not in obj.PropertiesList:             # PartDesign'in kendi Refine ozelligi (yoksa ekle)
            obj.addProperty("App::PropertyBool", "Refine", "Base", "Merge coplanar faces of the result")
        obj.Refine = True
        obj.Proxy = self

    def execute(self, obj):
        tool = build_tool(obj)
        if isinstance(tool, str):
            raise RuntimeError(T(tool))
        base = obj.BaseFeature.Shape if obj.BaseFeature is not None else None
        obj.AddSubShape = tool
        if base is None or base.isNull():
            res = tool
        elif obj.TypeId == "PartDesign::FeatureSubtractivePython":
            res = base.cut(tool)
        else:
            res = base.fuse(tool)
        if obj.Refine:
            res = res.removeSplitter()
        obj.Shape = res

    def dumps(self):
        return None

    def loads(self, state):
        return None


class ViewProviderRailPipe:
    def __init__(self, vobj, icon):
        vobj.Proxy = self
        self.icon = icon

    def attach(self, vobj):
        self.Object = vobj.Object

    def getIcon(self):
        import os
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "icons")
        groove = getattr(self, "Object", None) is not None and \
            self.Object.TypeId == "PartDesign::FeatureSubtractivePython"
        return os.path.join(d, "subd_railgroove.svg" if groove else "subd_railpipe.svg")

    def dumps(self):
        return None

    def loads(self, state):
        return None


def initial_points(fa, fb, subtractive, k=3):
    """Baslangic kilavuz noktalari: iki yuz arasinda A'ya ve B'ye dik yumusak egri uzerinde esit aralikli k nokta;
    tek yuzde normal yonunde yuz boyunun iki kati uzunlukta duz cizgi (son nokta serbest uc)."""
    out = -1.0 if subtractive else 1.0
    ca, na = face_frame(fa)
    if fb is None:
        L = 2.0 * fa.BoundBox.DiagonalLength
        return [ca + na * out * (L * j / k) for j in range(1, k + 1)]
    cb, nb = face_frame(fb)
    L = (cb - ca).Length
    curve = [_hermite(ca, cb, na * out * L, nb * -out * L, j / 128.0) for j in range(129)]
    acc = [0.0]
    for p, q in zip(curve, curve[1:]):
        acc.append(acc[-1] + (q - p).Length)
    pts = []
    for j in range(1, k + 1):
        d = acc[-1] * j / (k + 1)
        m = next(i for i, a in enumerate(acc) if a >= d)
        pts.append(curve[m])
    return pts


def make(body, start, end=None, subtractive=False):
    """Body'de yeni Rail Pipe / Rail Groove: start / end (nesne, "FaceN") ya da None. Nesneyi dondurur."""
    obj = body.newObject("PartDesign::FeatureSubtractivePython" if subtractive else "PartDesign::FeatureAdditivePython",
                         "RailGroove" if subtractive else "RailPipe")
    RailPipe(obj)
    obj.StartFace = (start[0], [start[1]])
    if end is not None:
        obj.EndFace = (end[0], [end[1]])
    fa = start[0].Shape.getElement(start[1])
    fb = end[0].Shape.getElement(end[1]) if end is not None else None
    obj.GuidePoints = initial_points(fa, fb, subtractive)
    if App.GuiUp:
        ViewProviderRailPipe(obj.ViewObject, None)
    return obj


def repair_proxies(doc):
    """Proxy'si kopmus Rail Pipe / Groove nesnelerini onar (bkz. subd_modeler.repair_proxies)."""
    fixed = []
    for o in (doc.Objects if doc is not None else []):
        if PROPS in getattr(o, "PropertiesList", []) and type(getattr(o, "Proxy", None)).__name__ != "RailPipe":
            o.Proxy = RailPipe.__new__(RailPipe)
            if App.GuiUp and type(getattr(o.ViewObject, "Proxy", None)).__name__ != "ViewProviderRailPipe":
                vp = ViewProviderRailPipe.__new__(ViewProviderRailPipe)
                o.ViewObject.Proxy = vp
            o.touch()
            fixed.append(o)
    for o in fixed:
        o.recompute()
    return fixed
