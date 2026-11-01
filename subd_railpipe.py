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
PREVIEW = False           # ekranda surukleme sirasinda True: yavas (kaliteli) supurme yontemi atlanir
LAST_MODE = None          # son supurmenin yontemi (tani)


def _is_sketch(o):
    return o is not None and (o.TypeId == "Sketcher::SketchObject" or o.isDerivedFrom("Part::Part2DObject"))


def profile(link):
    """PropertyLinkSub'tan profil: (dis tel, merkez, normal, katinin yuzu mu) ya da None / hata metni. Katinin duz
    yuzu: normali disari bakar, boru oradan cikar. Sketch (alt eleman yok ya da kenar): en buyuk kapali teli, normali
    sketch duzleminin Z'si. Katisi olmayan yuz (SubShapeBinder, tek basina yuz / kabuk; alt eleman yoksa tek yuzu):
    normalinin yonu rastgele (Binder'da boru ters taraftan cikiyordu), sketch gibi isareti egriye gore secilir."""
    if not link or link[0] is None:
        return None
    o, subs = link[0], list(link[1] or [])
    if _is_sketch(o):
        ws = [w for w in o.Shape.Wires if w.isClosed()]
        if not ws:
            return "Rail Pipe: the sketch has no closed profile"
        w = max(ws, key=lambda w: w.BoundBox.DiagonalLength)
        return w, Part.Face(w).CenterOfMass, o.Placement.Rotation.multVec(V(0, 0, 1)), False
    f = None
    if subs:
        try:
            f = o.Shape.getElement(subs[0])
        except Exception:
            f = None
    elif len(o.Shape.Faces) == 1:                          # agacta secilen tek yuzlu nesne (or. Binder)
        f = o.Shape.Faces[0]
    if f is None or f.ShapeType != "Face":
        return "Rail Pipe: the start face is missing"
    if not (f.isPlanar() if hasattr(f, "isPlanar") else f.Surface.TypeId == "Part::GeomPlane"):
        return "Rail Pipe: the start face must be flat"
    c, n = face_frame(f)
    return f.OuterWire, c, n, bool(o.Shape.Solids)


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


def _aux_spine(bs, t0, ref, size, n=80):
    """Yardimci omurga: egri boyunca burulmasiz tasinan yan yonun (ref, t0'a dik) size kadar ucundan gecen egri."""
    r = V(ref) - t0 * V(ref).dot(t0)
    if r.Length < 1e-9:
        r = t0.cross(V(1, 0, 0)) if abs(t0.x) < 0.9 else t0.cross(V(0, 1, 0))
    r.normalize()
    a, b = bs.FirstParameter, bs.LastParameter
    steps = 8 * n
    R, prev, pts = App.Rotation(), V(t0), []
    for k in range(steps + 1):
        u = a + (b - a) * k / steps
        t = bs.tangent(u)[0]
        R = App.Rotation(prev, t).multiply(R)
        prev = t
        if k % 8 == 0:
            pts.append(bs.value(u) + R.multVec(r) * size)
    try:
        c = Part.BSplineCurve()
        c.interpolate(pts)
        return Part.Wire(c.toShape())
    except Exception:
        return None


def _tight_box(sh):
    """Gercek sinir kutusu: BoundBox B-spline yuzeyin kontrol noktalarina bakiyor, C1 yaklasimda kutu sekilden cok
    tasiyor ve dogru supurme 'savrulmus' sayilip reddediliyordu."""
    try:
        return sh.optimalBoundingBox()
    except Exception:
        return sh.BoundBox


def _volume_ok(base, tool, res, sub, rel=0.005):
    """Boolean sonucunun hacmi beklenenle tutuyor mu (ekleme: govde + arac - kesisim, oyma: govde - kesisim)."""
    try:
        common = base.common(tool).Volume
        want = base.Volume - common if sub else base.Volume + tool.Volume - common
        return abs(res.Volume - want) <= rel * max(base.Volume, tool.Volume, 1e-9)
    except Exception:
        return False


def _bop_clean(sh):
    try:
        sh.check(True)
        return True
    except Exception:
        return False


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


def _on_body(p, base):
    """Govdenin yuzeyinde duran profil (govdenin kendi yuzunun Binder'i, yuze cizilmis sketch): katinin yuzu gibi
    davranilir; normali govdenin disina cevrilir, uc govdenin icine tasirilir. Yoksa cakisan yuzler boolean'i
    bozuyordu (sonuc bos ya da iki ayri kati) ve her basarisiz deneme saniyeler suruyordu."""
    if p is None or p[3] or base is None or base.isNull() or not base.Solids:
        return p
    w, c, n, _ = p
    size = max(w.BoundBox.DiagonalLength, 1.0)
    try:
        # yuzun icinden ornek noktalar: cogu govdenin yuzeyindeyse profil govdenin uzerinde. Yalniz kose noktalari
        # yaniltiyordu: delik agzini kapatan disk (Binder) kenariyla govdeye degiyor ama uzerinde degil.
        f = Part.Face(w)
        u0, u1, v0, v1 = f.ParameterRange
        on, inner, e = [], 0, 1e-3 * size
        for i in range(1, 6):
            for j in range(1, 6):
                q = f.Surface.value(u0 + (u1 - u0) * i / 6.0, v0 + (v1 - v0) * j / 6.0)
                if not f.isInside(q, 1e-6, True):
                    continue
                inner += 1
                if base.distToShape(Part.Vertex(q))[0] <= 1e-4 * size:
                    on.append(q)
        if not inner or len(on) < 0.6 * inner:
            return p
        vote = sum(int(base.isInside(q + n * e, 1e-7, False)) - int(base.isInside(q - n * e, 1e-7, False))
                   for q in on)                            # hangi yan govdenin ici (cogunluk)
        if vote > 0:                                       # normal govdenin icine bakiyor: disa cevir
            n = -n
    except Exception:
        return p
    return w, c, n, True


def _directions(pa, pb, pts, sub, reversed_):
    """A'dan cikis (ta) ve B'ye giris (tb) yonu. Katinin yuzu: Pipe disari, Groove iceri (Reversed cevirir); sketch
    ya da katisiz yuz (Binder): normalin isareti obur profile (tek profilde ilk kilavuz noktasina) dogru; B'de son
    noktadan B'ye. Reversed cevirir."""
    out = -1.0 if (sub != bool(reversed_)) else 1.0
    wa, ca, na, face_a = pa
    if face_a:
        ta = na * out
    else:
        to = pb[1] if pb is not None else (pts[0] if pts else ca + na)   # iki profilde obur profile dogru
        ta = na if na.dot(to - ca) >= 0 else -na
        if reversed_:
            ta = -ta
    tb = None
    if pb is not None:
        wb, cb, nb, face_b = pb
        if face_b:
            tb = nb * -out                                 # B'ye disaridan (Groove: iceriden) girer
        else:
            frm = pts[-1] if pts else ca
            tb = nb if nb.dot(cb - frm) >= 0 else -nb
    return ta, tb


def guide_nodes(obj):
    """Ekranda duzenleme icin: (dugumler [A merkezi, kilavuz noktalari, (B merkezi)], ta, tb) ya da hata metni."""
    pa = profile(obj.StartFace)
    if pa is None or isinstance(pa, str):
        return pa or "Rail Pipe: the start face is missing"
    pb = profile(obj.EndFace)
    if isinstance(pb, str):
        return pb
    pts = list(obj.GuidePoints)
    base = obj.BaseFeature.Shape if obj.BaseFeature is not None else None
    pa, pb = _on_body(pa, base), _on_body(pb, base)
    ta, tb = _directions(pa, pb, pts, obj.TypeId == "PartDesign::FeatureSubtractivePython", obj.Reversed)
    return [pa[1]] + pts + ([pb[1]] if pb is not None else []), ta, tb


def build_tool(obj):
    """Supurulen ilk uygun kati (Body koordinatinda) ya da hata metni."""
    for t in build_tools(obj):
        return t
    return "Rail Pipe: the sweep gave no solid"


def build_tools(obj):
    """Supurulen katilar, denenecek sirayla (bkz. asagidaki yontem sirasi); hicbiri olmazsa tek bir hata metni."""
    pa = profile(obj.StartFace)
    if pa is None:
        yield "Rail Pipe: the start face is missing"
        return
    if isinstance(pa, str):
        yield pa
        return
    pb = profile(obj.EndFace)
    if isinstance(pb, str):
        yield pb
        return
    pts = list(obj.GuidePoints)
    if pb is None and not pts:
        yield "Rail Pipe: a free end needs at least one guide point"
        return
    base = obj.BaseFeature.Shape if obj.BaseFeature is not None else None
    pa, pb = _on_body(pa, base), _on_body(pb, base)
    ta, tb = _directions(pa, pb, pts, obj.TypeId == "PartDesign::FeatureSubtractivePython", obj.Reversed)
    ca = pa[1]
    N = [ca] + pts + ([pb[1]] if pb is not None else [])
    for i in range(len(N) - 1):
        if (N[i + 1] - N[i]).Length < 1e-6:
            yield "Rail Pipe: two guide points are on top of each other"
            return
    # yuz uclari yuzun obur yanindan baslar / biter (Pipe: govdenin icinden, Groove: disindan): supurmenin yaklasik
    # ucu yuzle tam cakisinca boolean cok kucuk kenar / kendini kesen kose birakiyor, sonraki adimlar bozuluyordu.
    # Sketch uclari tam yerinde (cakisan yuz yok).
    ea = min(0.5, 0.05 * pa[0].BoundBox.DiagonalLength) if pa[3] else 0.0
    eb = min(0.5, 0.05 * pb[0].BoundBox.DiagonalLength) if pb is not None and pb[3] else 0.0
    fb = pb
    bs = guide_spline(N, ta, tb, obj.Tension, ea, eb)
    spine = Part.Wire(bs.toShape())
    wa = _oriented(pa[0], ta).copy()
    wa.translate(-ta * ea)
    ca = ca - ta * ea
    if fb is not None:
        wb = _oriented(pb[0], tb).copy()
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
    # sirayla denenir, gecerli ve egrinin cevresinde kalan ilk sonuc alinir:
    # 1. yardimci omurga: kendi burulmasiz cercevemizin (egri boyunca en kisa donmelerle tasinan yan yon) ucundan
    #    gecen egri kesitin yonunu verir. CorrectedFrenet kesiti egri boyunca birkac derece salindiriyordu (dikis
    #    ve kenar cizgileri borunun cevresinde dalgali, ~2 mm r=40'ta); bununla ~0.01 derece.
    # 2. CorrectedFrenet (egrilige dayanir, duz cizgide tanimsiz kalip 'No geometry' veriyordu).
    # 3. duzlemsel egride sabit binormal (duz cizgide de calisir, ama iki profilli kivrimda bazen gecersiz).
    # Her birinde once C1 yaklasim (C0 kenar / yuz kalirsa sonraki boolean'lar bozuluyor), sonra C1'siz (iki
    # profilde C1 yaklasim yuzeyi savurabiliyor).
    size = max(wa.BoundBox.DiagonalLength, 1.0)
    auxs = [] if PREVIEW else [a for a in (_aux_spine(bs, ta, first - ca, size * k) for k in (1.0, 2.0)) if a]
    # Yardimci omurga OCCT'de kirilgan: bazi ayarlarda isValid gecip ici bozuk kati veriyor, sonraki boolean boruyu
    # kaybediyordu; egrisel eslesme kapaliyken bazen dakikalarca takiliyor. Olcumlerde egrisel eslesme + C1 hic
    # takilmadi: yalniz o denenir, yardimci egri profil boyu ve iki kati uzakta (~0.3-1 s), sonra CorrectedFrenet.
    # Dogruluk boolean sonucunun hacmiyle denetlenir (RailPipe.execute); BOP kontrolu sert kivrimda 25 s
    # surebiliyordu. Ekranda surukleme sirasinda (PREVIEW) hic denenmez: birakinca kaliteli hali hesaplanir.
    tries = [(("aux", a, True), True) for a in auxs]
    tries += [((m, a, None), c1) for m, a in [("frenet", None)] + ([("binormal", binormal)] if binormal is not None
                                                                    else []) for c1 in ((False, True) if PREVIEW
                                                                                        else (True, False))]
    given = False
    for (mode, arg, cureq), c1 in tries:
        ps = Part.BRepOffsetAPI.MakePipeShell(spine)
        if mode == "aux":
            ps.setAuxiliarySpine(arg, cureq, 0)
        elif mode == "binormal":
            ps.setBiNormalMode(arg)
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
        if not PREVIEW and not tool.isInside(mid, 1e-7, True):   # ters yuzlu (ici disina donuk) cikabiliyor:
            # cevir (onizlemede gerekmez: yalniz boolean icin onemli, nokta testi surukleme adiminin 3/4'u)
            tool.reverse()
        if tool.isValid() and (PREVIEW or box.isInside(_tight_box(tool))):   # onizleme: kutu denetimi yok (hiz)
            global LAST_MODE
            LAST_MODE = (mode, cureq, c1)                  # tani icin: son kullanilan supurme yontemi
            yield tool
            given = True
            continue
        err = "Rail Pipe: the sweep gave an invalid solid"
    if not given:
        yield err


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
        base = obj.BaseFeature.Shape if obj.BaseFeature is not None else None
        if base is not None and base.isNull():
            base = None
        sub = obj.TypeId == "PartDesign::FeatureSubtractivePython"
        first = None
        # boolean sonucu da denetlenir: gecersizse, govde tek katiyken bolunduyse ya da hacmi tutmuyorsa (ekleme:
        # govde + boru - kesisim, oyma: govde - kesisim; supurme tek basina temiz gorunup boolean'da boruyu
        # kaybedebiliyordu) siradaki supurme denenir; hicbiri olmazsa ilk gecerli sonuc
        for tool in build_tools(obj):
            if isinstance(tool, str):
                raise RuntimeError(T(tool))
            res = tool if base is None else base.cut(tool) if sub else base.fuse(tool)
            if obj.Refine:
                res = res.removeSplitter()
            if not res.isValid() or not res.Solids:
                continue
            if base is not None and not _volume_ok(base, tool, res, sub):
                continue
            if first is None:
                first = (tool, res)
            if base is None or len(res.Solids) <= max(len(base.Solids), 1):
                first = (tool, res)
                break
        if first is None:
            raise RuntimeError(T("Rail Pipe: the sweep gave an invalid solid"))
        obj.AddSubShape, obj.Shape = first

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

    def doubleClicked(self, vobj):
        import FreeCADGui as Gui
        Gui.ActiveDocument.setEdit(vobj.Object, 0)
        return True

    def setEdit(self, vobj, mode=0):
        if mode != 0:
            return None
        import subd_railpipe_gui                           # ekranda kilavuz duzenleme (tutamaklar + gorev paneli)
        subd_railpipe_gui.start(vobj.Object)
        return True

    def unsetEdit(self, vobj, mode=0):
        if mode != 0:
            return None
        import subd_railpipe_gui
        subd_railpipe_gui.stop()
        return True

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


def initial_points(pa, pb, subtractive, k=3):
    """Baslangic kilavuz noktalari (pa / pb: profile()): iki profil arasinda A'ya ve B'ye dik yumusak egri uzerinde
    esit aralikli k nokta; tek profilde cikis yonunde profil boyunun iki kati uzunlukta duz cizgi (son nokta serbest
    uc)."""
    ta, tb = _directions(pa, pb, [], subtractive, False)
    ca = pa[1]
    if pb is None:
        L = 2.0 * pa[0].BoundBox.DiagonalLength
        return [ca + ta * (L * j / k) for j in range(1, k + 1)]
    cb = pb[1]
    L = (cb - ca).Length
    curve = [_hermite(ca, cb, ta * L, tb * L, j / 128.0) for j in range(129)]
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
    """Body'de yeni Rail Pipe / Rail Groove: start / end (nesne, "FaceN"), sketch icin (sketch, None); end None:
    serbest uc. Nesneyi dondurur; profil kullanilamiyorsa hata metni (nesne olusturulmaz)."""
    link = lambda p: (p[0], [p[1]] if p[1] else [])
    pa = profile(link(start))
    pb = profile(link(end)) if end is not None else None
    for p in (pa, pb):
        if isinstance(p, str):
            return p
    base = body.Tip.Shape if getattr(body, "Tip", None) is not None else None   # yeni adimin BaseFeature'i
    pa, pb = _on_body(pa, base), _on_body(pb, base)
    obj = body.newObject("PartDesign::FeatureSubtractivePython" if subtractive else "PartDesign::FeatureAdditivePython",
                         "RailGroove" if subtractive else "RailPipe")
    RailPipe(obj)
    obj.StartFace = link(start)
    if end is not None:
        obj.EndFace = link(end)
    obj.GuidePoints = initial_points(pa, pb, subtractive)
    if App.GuiUp:
        ViewProviderRailPipe(obj.ViewObject, None)
        obj.ViewObject.Deviation = 0.1                     # varsayilan 0.5: kavisli boru ekranda dalgali gorunuyordu
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
