# -*- coding: utf-8 -*-
# SubD Modeler - LightWave-style cage + SubPatch (Catmull-Clark) modeling for FreeCAD
# Copyright (c) 2026 Sedat Kuran - MIT License (see LICENSE)
#
# While modeling nothing is written to the document: the cage and its Catmull-Clark preview are drawn
# directly in the Coin3D scene, so every operation is instant. The cage is stored in a "SubD" document
# object (every 20 s and on close); select it and open the modeler again to keep editing.
#
# Shortcuts (in the 3D view):
#   1 / 2 / 3: point / edge / polygon mode     Ctrl+G / Ctrl+H: point / polygon mode (LightWave)
#   Click: select   Shift+click: add/remove   Drag: box select   A: select all   L: select loop   X: ghost
#   E: extrude (edge mode: open edges; wheel = segments; X / Y / Z key or panel = axis lock, none = normal)   B: bevel / Smooth Shift (edge mode: edge bevel, wheel = segments)   I: inset
#   K: loop cut + slide   J: LightWave knife (draw a line, live plane cut)   Shift+V: mirror across a drawn line
#   C: crease   Del: delete
#   Knife (K / J) and Select loop stay on (yellow button) until K / J / the button / Esc / right click
#   T: move gizmo   G: free move (mouse)   R: rotate gizmo   H: stretch / size gizmo   Y: X symmetry   Tab / S: SubPatch
#   Ctrl+C / Ctrl+X / Ctrl+V: copy / cut / paste polygons (object mode Ctrl+V: as a new SubD object)
#   Ctrl+Z: undo   Ctrl+Y: redo   F: save to the document   N: Options panel   Esc / right click: cancel
#   After every operation the small panel in the viewport shows its values; type a value + Enter and the
#   operation is redone with it. Ctrl while dragging = steps (mm / degrees).
#
# Code comments are mostly in Turkish (the author's language); the user interface is English by default
# with a Turkish option (subd_i18n.py).
import colorsys
import copy
import json
import math
import os
import time

import numpy as np
import FreeCAD as App
import FreeCADGui as Gui
import Mesh
from pivy import coin
from PySide import QtCore, QtGui, QtWidgets

from subd_i18n import T
import subd_i18n


# ---------- Catmull-Clark (numpy, vektorize) ----------

CREASE_HARD = 10.0       # bu keskinlik her bolme seviyesinden buyuk: kenar tamamen keskin kalir
BOUNDARY_RULES = ("Edge and corner", "Edge only")   # acik kenar kurali (OpenSubdiv adlari)
RULES = ("Edge and corner", False)                  # varsayilan kurallar: (acik kenar kurali, ucgen yumusatma)


def catmull_clark(V, faces, crease=None, rules=RULES):
    """Bir Catmull-Clark adimi. crease: {frozenset((a, b)): keskinlik}. Keskinlik >= 1 kenar bu adimda keskin
    kurallarla bolunur, 0-1 arasi yumusakla harmanlanir; cocuk kenarlara keskinlik - 1 gecer (yari keskin kenar).
    rules: (acik kenar kurali, ucgen yumusatma). "Edge and corner": tek yuzlu acik kenar koseleri sivri kalir,
    "Edge only": onlar da yuvarlanir. Ucgen yumusatma: ucgenlere komsu kenar noktalari Loop'a yakin agirlikla
    (OpenSubdiv "smooth triangle" gibi; ucgen yelpaze kapaklarindaki cukurlasmayi azaltir).
    Dondurur: (noktalar, dortgenler, cocuk crease sozlugu)."""
    if isinstance(faces, np.ndarray):
        deg = np.full(len(faces), faces.shape[1])
        fv = faces.ravel()
    else:
        deg = np.array([len(f) for f in faces])
        fv = np.fromiter((v for f in faces for v in f), dtype=np.int64, count=int(deg.sum()))
    nF, N = len(deg), len(V)
    fidx = np.repeat(np.arange(nF), deg)
    start = np.concatenate([[0], np.cumsum(deg)[:-1]])
    pos = np.arange(len(fv)) - start[fidx]
    nxt = start[fidx] + (pos + 1) % deg[fidx]
    prv = start[fidx] + (pos - 1) % deg[fidx]

    fp = np.zeros((nF, 3))
    np.add.at(fp, fidx, V[fv])
    fp /= deg[:, None]

    a, b = fv, fv[nxt]
    keys = np.minimum(a, b) * N + np.maximum(a, b)
    uniq, einv = np.unique(keys, return_inverse=True)
    einv = einv.ravel()
    E = len(uniq)
    ea, eb = uniq // N, uniq % N
    ecnt = np.bincount(einv, minlength=E)
    efs = np.zeros((E, 3))
    np.add.at(efs, einv, fp[fidx])
    mid = (V[ea] + V[eb]) / 2
    ep = mid.copy()
    inner = ecnt >= 2
    ep[inner] = mid[inner] / 2 + efs[inner] / (2 * ecnt[inner][:, None])
    if rules[1] and (deg == 3).any():
        # iki yuzlu kenar: ep = wv (a + b) + toplam(wf * yuz noktasi); ucgen wf = 3/8, dortgen 1/4. Iki ucgen
        # arasinda tam Loop kenar kurali (3/8 uclar, 1/8 karsi koseler) cikar; yalniz dortgenlerde degismez.
        wf = np.where(deg[fidx] == 3, 0.375, 0.25)
        ws, wfs = np.zeros(E), np.zeros((E, 3))
        np.add.at(ws, einv, wf)
        np.add.at(wfs, einv, wf[:, None] * fp[fidx])
        two = ecnt == 2
        ep[two] = (1 - ws[two])[:, None] * mid[two] + wfs[two]

    n = np.bincount(fv, minlength=N).astype(float)
    fsum = np.zeros((N, 3))
    np.add.at(fsum, fv, fp[fidx])
    rsum = np.zeros((N, 3))
    np.add.at(rsum, ea, mid)
    np.add.at(rsum, eb, mid)
    m = np.bincount(np.concatenate([ea, eb]), minlength=N).astype(float)
    Vn = V.copy()
    used = n > 0
    nn = n[used][:, None]
    Vn[used] = (fsum[used] / nn + 2 * rsum[used] / m[used][:, None] + (nn - 3) * V[used]) / nn

    # acik kenarlar: B-spline egrisi kurali (limitte nokta = (onceki + 4 kendisi + sonraki) / 6). Keskin
    # kenarlardan ONCE: 3+ keskin kenarli acik kenar noktasi kose olarak yerinde kalsin (eskiden bu kural
    # sonra geldigi icin acik kenardaki keskin koseler de yuvarlaniyordu).
    bd = ecnt == 1
    if bd.any():
        bsum, bcnt = np.zeros((N, 3)), np.zeros(N)
        np.add.at(bsum, ea[bd], V[eb[bd]])
        np.add.at(bsum, eb[bd], V[ea[bd]])
        np.add.at(bcnt, ea[bd], 1)
        np.add.at(bcnt, eb[bd], 1)
        bv = bcnt == 2
        Vn[bv] = 0.75 * V[bv] + 0.125 * bsum[bv]
        corner = bv & (n == 1) & (rules[0] == BOUNDARY_RULES[0])   # tek yuzlu kose: yerinde kalir
    else:                                              # (OpenSubdiv varsayilani "edge and corner"); keskin kenar
        corner = np.zeros(N, bool)                     # kurallarindan SONRA uygulanir (yoksa egri kurali ezerdi)

    # keskin kenarlar (crease)
    child = {}
    if crease:
        ck = np.array([min(e) * N + max(e) for e in crease], dtype=np.int64)
        cs = np.array(list(crease.values()), float)
        at = np.minimum(np.searchsorted(uniq, ck), E - 1)
        ok = uniq[at] == ck
        s = np.zeros(E)
        s[at[ok]] = cs[ok]
        w = np.clip(s, 0, 1)[:, None]
        ep = ep * (1 - w) + mid * w                    # keskin kenar noktasi: kenarin ortasi
        cr = s > 0
        if cr.any():
            nc = np.bincount(np.concatenate([ea[cr], eb[cr]]), minlength=N)
            nb, ssum = np.zeros((N, 3)), np.zeros(N)
            np.add.at(nb, ea[cr], V[eb[cr]])
            np.add.at(nb, eb[cr], V[ea[cr]])
            np.add.at(ssum, ea[cr], s[cr])
            np.add.at(ssum, eb[cr], s[cr])
            Vc = V.copy()                              # 3+ keskin kenar: kose, nokta yerinde kalir
            two = nc == 2                              # 2 keskin kenar: egri kurali
            Vc[two] = 0.75 * V[two] + 0.125 * nb[two]
            many = nc >= 2
            wv = np.clip(ssum[many] / nc[many], 0, 1)[:, None]
            Vn[many] = Vn[many] * (1 - wv) + Vc[many] * wv
        for e in np.nonzero(s > 1)[0]:
            child[_edge(int(ea[e]), N + int(e))] = s[e] - 1
            child[_edge(int(eb[e]), N + int(e))] = s[e] - 1
    Vn[corner] = V[corner]

    V2 = np.concatenate([Vn, ep, fp])
    quads = np.stack([fv, N + einv, N + E + fidx, N + einv[prv]], 1)
    return V2, quads, child


def subdivide(V, F, levels, crease=None, rules=RULES):
    for _ in range(levels):
        V, F, crease = catmull_clark(V, F, crease, rules)
    return V, F


SHARP_COS = math.cos(math.radians(60.0))   # yuzey onizlemesi: komsu yuzler arasinda bundan buyuk aci sert kose gibi golgelenir


def surface_normals(S, Q, crease=None):
    """Dortgen yuzeyin kose normalleri (dortgen basina 4; n x 4 x 3): her nokta, cevresindeki yuzlerin alana gore
    agirlikli ortalamasini alir (matcap piksel basina bunlardan okur). Keskin kenar (crease) ve 60 dereceden sert kirik
    iki yana ayri golgelenir."""
    fn = np.cross(S[Q[:, 2]] - S[Q[:, 0]], S[Q[:, 3]] - S[Q[:, 1]])
    unit = lambda a: a / np.maximum(np.linalg.norm(a, axis=-1, keepdims=True), 1e-30)
    vn = np.zeros_like(S)
    np.add.at(vn, Q.ravel(), np.repeat(fn, 4, 0))
    out = unit(vn)[Q]
    N = len(S)
    a, b = Q.ravel(), np.roll(Q, -1, 1).ravel()
    keys = np.minimum(a, b) * N + np.maximum(a, b)
    order = np.argsort(keys, kind="stable")
    k = keys[order]
    pair = np.nonzero(k[1:] == k[:-1])[0]                  # ayni kenari paylasan iki kose
    fu = unit(fn)
    fi, fj = order[pair] // 4, order[pair + 1] // 4
    sharp = set(k[pair][(fu[fi] * fu[fj]).sum(1) < SHARP_COS].tolist())
    if crease:
        sharp.update(min(e) * N + max(e) for e in crease)
    if not sharp:
        return out
    ks = np.array(sorted(sharp), dtype=np.int64)
    corners = {}                                           # yalniz sert kenara degen noktalar: kose gruplari ayrilir
    for v in set((ks // N).tolist()) | set((ks % N).tolist()):
        corners[v] = []
    for ci, v in enumerate(a.tolist()):
        if v in corners:
            corners[v].append(ci)
    for v, cs in corners.items():
        group = {c: c for c in cs}

        def root(c):
            while group[c] != c:
                group[c] = group[group[c]]
                c = group[c]
            return c
        side = {}
        for c in cs:                                       # bu kosenin iki kenari; sert olmayanla komsu yuzu birlesir
            q, i = c // 4, c % 4
            for w in (Q[q, (i + 1) % 4], Q[q, i - 1]):
                w = int(w)
                if min(v, w) * N + max(v, w) in sharp:
                    continue
                if w in side:
                    group[root(c)] = root(side[w])
                else:
                    side[w] = c
        total = {}
        for c in cs:
            r = root(c)
            total[r] = total.get(r, 0) + fn[c // 4]
        for c in cs:
            out[c // 4, c % 4] = unit(total[root(c)])
    return out


# ---------- Kafes (poligon mesh) ve duzenleme islemleri ----------

def _edge(a, b):
    return frozenset((a, b))


class Cage:
    def __init__(self, V=None, F=None, C=None, W=None):
        self.V = np.zeros((0, 3)) if V is None else V
        self.F = [] if F is None else F
        self.C = {} if C is None else C   # keskin kenarlar: frozenset((a, b)) -> keskinlik
        # serbest kenarlar (yuzsuz, LightWave 2 noktali cokgen gibi): kopyalanan kenarlar, uzerine yuzey kurulur.
        # Bir yuz ayni kenari kullanmaya baslarsa burada kalmasi zararsiz (kayitta ayiklanir: free_edges).
        # Hicbir yuz / kenarda olmayan noktalar da kafeste kalabilir (serbest noktalar).
        self.W = set() if W is None else W
        self.rails = []                   # Rail bridge kilavuzlari (rail_* fonksiyonlari; nokta indeksleriyle)
        self.links = []                   # Sketch baglari (link_update; nokta konumlariyla, From face > Edge)
        self.dims = []                    # Distance olculeri (kalici; _dim_* ve DIM_KEEP; nokta indeksleriyle)
        self.groups = []                  # Parametrik gruplar (group_update; datum tutamagina bagli noktalar)

    @staticmethod
    def cube(size=20.0):
        """size: kenar ya da (en, derinlik, yukseklik) mm; merkez orijinde."""
        hx, hy, hz = np.broadcast_to(np.asarray(size, float), 3) / 2
        V = np.array([[x, y, z] for z in (-hz, hz) for y in (-hy, hy) for x in (-hx, hx)], float)
        F = [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]]
        return Cage(V, F)

    @staticmethod
    def box(size=20.0, segs=1):
        """Segmentli kutu (LightWave Box): size (en, derinlik, yukseklik) mm, segs (X, Y, Z) dilim; merkez
        orijinde. Kenar ve kose noktalari komsu yuzlerde ortak; yuzler disari bakar."""
        n = np.maximum(np.broadcast_to(np.asarray(segs, int), 3), 1)
        h = np.broadcast_to(np.asarray(size, float), 3) / 2
        index, V, F = {}, [], []

        def vid(p):
            if p not in index:
                index[p] = len(V)
                V.append([-h[a] + 2 * h[a] * p[a] / n[a] for a in range(3)])
            return index[p]
        for ax in range(3):
            for side in (0, 1):
                u, v = ((1, 2), (2, 0), (0, 1))[ax]        # u x v = +ax; eksi tarafta ters
                if not side:
                    u, v = v, u
                for a in range(n[u]):
                    for b in range(n[v]):
                        quad = []
                        for da, db in ((0, 0), (1, 0), (1, 1), (0, 1)):
                            p = [0, 0, 0]
                            p[ax], p[u], p[v] = n[ax] * side, a + da, b + db
                            quad.append(vid(tuple(p)))
                        F.append(quad)
        return Cage(np.array(V, float), F)

    @staticmethod
    def cylinder(diameter=100.0, height=100.0, segs=16, hsegs=1, topo=False):
        """Silindir: Z ekseni boyunca, merkez orijinde; yan duvar segs x hsegs dortgen. topo=False: kapaklar
        merkeze ucgen yelpaze (segs ucgen). topo=True: kapaklar tamamen dortgen (SubD'ye uygun): ortada k x k
        izgara (k = segs / 4, segs 4'un kati olmali) ve izgara siniri ile cevre arasinda bir dortgen halkasi."""
        R, h = diameter / 2.0, height / 2.0
        V, F = [], []
        a0 = -math.pi / 4 if topo else 0.0                 # topo: ilk cevre noktasi izgaranin sag alt kosesine bakar
        for j in range(hsegs + 1):                         # yan halkalar: nokta j * segs + i
            z = -h + height * j / hsegs
            V += [[R * math.cos(a0 + 2 * math.pi * i / segs), R * math.sin(a0 + 2 * math.pi * i / segs), z]
                  for i in range(segs)]
        ring = lambda j: [j * segs + i for i in range(segs)]
        for j in range(hsegs):
            lo, hi = ring(j), ring(j + 1)
            F += [[lo[i], lo[(i + 1) % segs], hi[(i + 1) % segs], hi[i]] for i in range(segs)]
        for z, rim, top in ((-h, ring(0), False), (h, ring(hsegs), True)):
            faces = []                                     # ust kapakta +Z'den bakinca saat yonu tersi
            if not topo:
                c = len(V)
                V.append([0.0, 0.0, z])
                faces = [[c, rim[i], rim[(i + 1) % segs]] for i in range(segs)]
            else:
                k = segs // 4
                r_in = R * k / (k + 2.0)                   # halka genisligi ~ izgara hucresi
                base = len(V)
                for b in range(k + 1):
                    for a in range(k + 1):
                        u, v = 2.0 * a / k - 1, 2.0 * b / k - 1
                        du, dv = u * math.sqrt(1 - v * v / 2), v * math.sqrt(1 - u * u / 2)   # kare -> daire
                        V.append([r_in * (u + du) / 2, r_in * (v + dv) / 2, z])               # yari yuvarlak
                g = lambda a, b: base + b * (k + 1) + a
                faces = [[g(a, b), g(a + 1, b), g(a + 1, b + 1), g(a, b + 1)] for b in range(k) for a in range(k)]
                edge = [g(k, b) for b in range(k)] + [g(k - a, k) for a in range(k)] + \
                    [g(0, k - b) for b in range(k)] + [g(a, 0) for a in range(k)]           # sinir, saat yonu tersi
                faces += [[edge[i], rim[i], rim[(i + 1) % segs], edge[(i + 1) % segs]] for i in range(segs)]
            F += faces if top else [f[::-1] for f in faces]
        return Cage(np.array(V, float), F)

    def copy(self):
        c = Cage(self.V.copy(), [list(f) for f in self.F], dict(self.C), set(self.W))
        c.rails = copy.deepcopy(self.rails)
        c.links = copy.deepcopy(self.links)
        c.dims = copy.deepcopy(self.dims)
        c.groups = copy.deepcopy(self.groups)
        return c

    def face_edges(self):
        return {_edge(f[i], f[(i + 1) % len(f)]) for f in self.F for i in range(len(f))}

    def edges(self):
        """Butun kenarlar: yuz kenarlari ve serbest kenarlar (cizim, secim, yakalama icin)."""
        return self.face_edges() | self.W if self.W else self.face_edges()

    def free_edges(self):
        """Hicbir yuzun kullanmadigi serbest kenarlar."""
        return self.W - self.face_edges() if self.W else set()

    def loose_points(self):
        """Hicbir yuz ya da serbest kenarda olmayan noktalar."""
        used = {v for f in self.F for v in f} | {v for e in self.W for v in e}
        return [v for v in range(len(self.V)) if v not in used]

    def delete_parts(self, edges=(), points=(), faces=()):
        """Kenar / nokta silme (Blender Delete Edges / Vertices gibi): verilen kenarlar (serbest olanlar kalkar,
        yuz kenariysa o kenari kullanan yuzler de silinir), noktalar (onlari kullanan yuzler, serbest kenarlar
        da) ve yuzler. Bosa dusen noktalar gider; onceden serbest olanlar kalir (silinenler haric)."""
        edges, points = set(edges), set(points)
        before = {v for f in self.F for v in f} | {v for e in self.W for v in e}
        self.W = {e for e in self.W if e not in edges and not e & points}
        dead = {fi for fi, f in enumerate(self.F) if fi in faces or points.intersection(f)
                or any(_edge(f[i], f[(i + 1) % len(f)]) in edges for i in range(len(f)))}
        self.F = [f for fi, f in enumerate(self.F) if fi not in dead]
        after = {v for f in self.F for v in f}
        self._compact([v for v in range(len(self.V)) if v not in points and (v in after or v not in before)])

    def edge_faces(self):
        ef = {}
        for fi, f in enumerate(self.F):
            for i in range(len(f)):
                ef.setdefault(_edge(f[i], f[(i + 1) % len(f)]), []).append(fi)
        return ef

    def normal(self, fi):
        p = self.V[self.F[fi]]
        n = np.cross(p, np.roll(p, -1, 0)).sum(0)
        ln = np.linalg.norm(n)
        return n / ln if ln else n

    def _layers(self, verts, segs):
        """verts noktalarinin segs kat kopyasi (hepsi ayni yerde): {nokta: [kendisi, kat 1, ..., kat segs]}."""
        L = {v: [v] for v in verts}
        for _ in range(segs):
            n = len(self.V)
            self.V = np.vstack([self.V, self.V[verts]])
            for j, v in enumerate(verts):
                L[v].append(n + j)
        return L

    def extrude(self, sel, segs=1, layers=None):
        """Secili yuzleri bolge olarak extrude eder (yan duvarlar segs dilim), tasinacak en ust noktalari
        dondurur (sirali noktalar sirasiyla). layers verilirse {taban noktasi: [taban, kat 1..segs]} doldurulur."""
        count = {}
        for fi in sel:
            f = self.F[fi]
            for i in range(len(f)):
                e = _edge(f[i], f[(i + 1) % len(f)])
                count[e] = count.get(e, 0) + 1
        verts = sorted({v for fi in sel for v in self.F[fi]})
        L = self._layers(verts, segs)
        top = {v: L[v][-1] for v in verts}
        sides = []
        for fi in sel:
            f = self.F[fi]
            for i in range(len(f)):
                a, b = f[i], f[(i + 1) % len(f)]
                e = _edge(a, b)
                if count[e] == 1:
                    sides += [[L[a][k], L[b][k], L[b][k + 1], L[a][k + 1]] for k in range(segs)]
                if e in self.C:                    # keskinlik yeni (tasinan) kenara gecer;
                    s = self.C[e]                  # sinirdakiler alt kenarda da kalir
                    self.C[_edge(top[a], top[b])] = s
                    if count[e] > 1:
                        del self.C[e]
            self.F[fi] = [top[v] for v in f]
        self.F += sides
        if layers is not None:
            layers.update(L)
        return [top[v] for v in verts]

    def extrude_edges(self, edges, segs=1, layers=None, wire_dir=None):
        """Acik (tek yuzlu) kenarlari yuzeyin devami gibi disari uzatir: her kenardan segs dortgen. Ic kenarlar
        atlanir (SubD'de yuzgec / manifold olmayan kenar olurdu). (sirali taban noktalari, birim disari yonleri)
        ya da hic acik kenar yoksa None dondurur; layers extrude'daki gibi doldurulur. Serbest kenarlar (yuzsuz)
        da uzar: duzlemsel kapali zincir kendi duzleminin normalinde, digerleri ekranda zincire dik (wire_dir:
        bakis yonu, yerel; normal ona dogru cevrilir)."""
        directed = self._open_directed(edges)
        ef = self.edge_faces()
        wires = [e for e in edges if e in self.W and e not in ef]
        if not directed and not wires:
            return None
        out = self._outward(directed)
        wired = []                                         # serbest kenarlar zincir sirasiyla: (a, b, None)
        for chain, closed in self.edge_chains(wires) or [(list(e), False) for e in wires]:
            pairs = list(zip(chain, chain[1:])) + ([(chain[-1], chain[0])] if closed else [])
            wired += [(a, b, None) for a, b in pairs]
            for v, d in self._wire_dirs(chain, closed, wire_dir).items():
                out[v] = d if v not in out else out[v] + d
        verts = sorted({v for a, b, _ in directed + wired for v in (a, b)})
        dirs = np.array([out[v] / max(np.linalg.norm(out[v]), 1e-12) for v in verts])
        L = self._layers(verts, segs)
        for a, b, fi in directed + wired:                  # yeni yuz b -> a sirasiyla: komsu yuzle ayni yonde
            self.F += [[L[b][k], L[a][k], L[a][k + 1], L[b][k + 1]] for k in range(segs)]
            e = _edge(a, b)
            self.W.discard(e)                              # artik yuz kenari
            if e in self.C:                                # keskinlik yeni acik kenara gecer
                self.C[_edge(L[a][-1], L[b][-1])] = self.C.pop(e)
        if layers is not None:
            layers.update(L)
        return verts, dirs

    def _wire_dirs(self, chain, closed, view=None):
        """Serbest kenar zincirinin uzama yonleri {nokta: birim yon}: kapali ve duzlemsel zincirde duzlem normali
        (hepsi ayni yon, bakisa dogru); yoksa her noktada zincire ve bakisa dik yon (ekranda yana)."""
        P = self.V[chain]
        view = np.array([0.0, 0.0, -1.0]) if view is None else np.asarray(view, float)
        if closed and len(chain) >= 3:
            c = P.mean(0)
            _, s, vt = np.linalg.svd(P - c)
            if s[1] > 1e-9 and s[2] < 0.05 * s[1]:         # duzleme yakin halka
                n = vt[2] * (-1.0 if vt[2].dot(view) > 0 else 1.0)   # bakan tarafa (bakis yonunun tersi)
                return {v: n for v in chain}
        out = {}
        for i, v in enumerate(chain):
            a = chain[i - 1] if i > 0 or closed else v
            b = chain[(i + 1) % len(chain)] if i < len(chain) - 1 or closed else v
            t = self.V[b] - self.V[a]
            d = np.cross(t, view)
            ln = np.linalg.norm(d)
            if ln < 1e-12:                                 # zincir bakis yonunde: herhangi bir dik yon
                d = np.cross(t, [1.0, 0, 0] if abs(t[0]) < 0.9 * np.linalg.norm(t) else [0, 1.0, 0])
                ln = np.linalg.norm(d)
            out[v] = d / max(ln, 1e-12)
        return out

    def _open_directed(self, edges):
        """Acik (tek yuzlu) kenarlar, yuzdeki sirayla: [(a, b, yuz)] (yuzde a -> b). Ic kenarlar atlanir."""
        ef = self.edge_faces()
        directed = []
        for e in edges:
            fs = ef.get(e, ())
            if len(fs) != 1:
                continue
            f = self.F[fs[0]]
            for i in range(len(f)):
                a, b = f[i], f[(i + 1) % len(f)]
                if _edge(a, b) == e:
                    directed.append((a, b, fs[0]))
                    break
        return directed

    def _outward(self, directed):
        """Acik kenar noktalarinda yuzeyin disari (devam) yonu: yuz duzleminde, kenara dik, yuzden uzaga. Yuz
        normaline gore saat yonu tersi sirada ic taraf a -> b'nin solunda (cross(n, t)); disari = cross(t, n).
        Noktada komsu kenarlarin ortalamasi, birim. {nokta: yon}"""
        acc = {}
        for a, b, fi in directed:
            o = np.cross(self.V[b] - self.V[a], self.normal(fi))
            ln = np.linalg.norm(o)
            for v in (a, b):
                acc[v] = acc.get(v, np.zeros(3)) + (o / ln if ln > 1e-12 else 0.0)
        return {v: d / max(np.linalg.norm(d), 1e-12) for v, d in acc.items()}

    def edge_chains(self, edges):
        """Kenarlari bagli zincirlere ayirir: [(sirali noktalar, kapali mi)]; dallanan nokta (3+ kenar) varsa None."""
        adj = {}
        for e in edges:
            a, b = tuple(e)
            adj.setdefault(a, []).append(b)
            adj.setdefault(b, []).append(a)
        if any(len(x) > 2 for x in adj.values()):
            return None
        seen, out = set(), []
        for start in sorted(adj, key=lambda v: (len(adj[v]) != 1, v)):   # acik zincirler ucundan baslasin
            if start in seen:
                continue
            chain, prev, cur = [start], None, start
            seen.add(start)
            while True:
                nxt = [w for w in adj[cur] if w != prev and w not in seen]
                if not nxt:
                    break
                prev, cur = cur, nxt[0]
                chain.append(cur)
                seen.add(cur)
            out.append((chain, len(chain) > 2 and len(adj[start]) == 2))
        return out

    def bridge(self, A, B, closed, cuts=0, smooth=0.0, profile=0.0, dA=None, dB=None, flip=False):
        """Iki eslenmis nokta dizisini (A[i] <-> B[i]) dortgen seritle baglar (Blender Bridge Edge Loops), arada
        cuts ara halka. smooth > 0: ara noktalar Hermite egrisiyle, iki uctaki yuzeye teget cikar (dA / dB: yuzeyden
        disari birim yonler, (n, 3)); profile: ortadan disari siskinlik (eksene gore, oransal). flip: yuz sirasi
        ters (komsu yuzlerle ayni yon). Yeni yuz indeksleri."""
        n = len(A)
        PA, PB = self.V[A], self.V[B]
        L = np.linalg.norm(PB - PA, axis=1)[:, None]
        cA, cB = PA.mean(0), PB.mean(0)
        ax = cB - cA
        ax = ax / np.linalg.norm(ax) if np.linalg.norm(ax) > 1e-9 else None
        rows = [list(A)]
        for j in range(1, cuts + 1):
            t = j / (cuts + 1.0)
            if smooth > 0 and dA is not None:
                h00, h10 = 2 * t ** 3 - 3 * t ** 2 + 1, t ** 3 - 2 * t ** 2 + t
                h01, h11 = -2 * t ** 3 + 3 * t ** 2, t ** 3 - t ** 2
                P = h00 * PA + h10 * dA * smooth * L + h01 * PB - h11 * dB * smooth * L
            else:
                P = PA + t * (PB - PA)
            if profile:
                r = P - (cA + t * (cB - cA))
                if ax is not None:
                    r = r - np.outer(r @ ax, ax)
                P = P + profile * math.sin(math.pi * t) * r
            rows.append(list(range(len(self.V), len(self.V) + n)))
            self.V = np.vstack([self.V, P])
        rows.append(list(B))
        nf = len(self.F)
        for j in range(len(rows) - 1):
            for i in range(n if closed else n - 1):
                k = (i + 1) % n
                q = [rows[j][i], rows[j][k], rows[j + 1][k], rows[j + 1][i]]
                self.F.append(q[::-1] if flip else q)
        return list(range(nf, len(self.F)))

    def hole(self, region, radius=0.0, segs=0, support=-1.0, cap=False):
        """Secili yuzlerin birbirine kenardan bagli her parcasina bir delik (karsilikli delikler icin). Otomatik
        degerler parcalarin hepsine sigan ortak degerlerdir: en kucuk yaricap ve destek, en cok segment. Bir parca
        olmazsa kafes degismez. cap: delik yerine dolu daire (Circle on face). Donus: _hole_patch gibi."""
        ef = self.edge_faces()
        patches, left = [], set(region)
        while left:
            todo, patch = [left.pop()], set()
            while todo:
                fi = todo.pop()
                patch.add(fi)
                f = self.F[fi]
                for i in range(len(f)):
                    for o in ef[_edge(f[i], f[(i + 1) % len(f)])]:
                        if o in left:
                            left.remove(o)
                            todo.append(o)
            patches.append(patch)
        if len(patches) < 2:
            return self._hole_patch(region, radius, segs, support, cap)
        res = self._hole_patches(patches, radius, segs, support, cap)
        if isinstance(res, str):
            return res
        if radius <= 0 or segs < 3 or support < 0:
            used = res[2]
            res = self._hole_patches(patches, radius if radius > 0 else min(u[0] for u in used),
                                     segs if segs >= 3 else max(u[1] for u in used),
                                     support if support >= 0 else min(u[2] for u in used), cap)
            if isinstance(res, str):
                return res
        self.__dict__.update(res[0].__dict__)
        return (res[1],) + res[2][0]

    def _hole_patches(self, patches, radius, segs, support, cap=False):
        """Parcalari bir kopyada sirayla deler: (kopya, yeni yuzler, [(radius, segs, support)]) ya da hata metni."""
        c, new, used = self.copy(), [], []
        patches = [set(p) for p in patches]
        for i, p in enumerate(patches):
            res = c._hole_patch(p, radius, segs, support, cap)
            if isinstance(res, str):
                return res
            gone = sorted(p)                               # silinenler kadar kayan indeksler
            shift = lambda fi: fi - int(np.searchsorted(gone, fi))
            new = [shift(f) for f in new] + list(res[0])
            patches[i + 1:] = [{shift(f) for f in q} for q in patches[i + 1:]]
            used.append(tuple(res[1:]))
        return c, new, used

    def _hole_patch(self, region, radius=0.0, segs=0, support=-1.0, cap=False):
        """Secili yuz bolgesinin ortasina daire delik acar (bolgenin yuzleri yeniden orulur). Bolgenin dis siniri
        (N nokta) delik cevresindeki M noktaya aci sirasiyla baglanir: M = N ise hepsi dortgen; M > N besgen,
        M < N ucgen olur. support > 0: delik kenarina destek halkasi (keskin, SubD'de daire kalir). radius / segs /
        support <= 0 ise otomatik. Daire bolge yuzeyine normal boyunca izdusurulur (hafif egri yuzeyde de olur).
        cap: delik acilmaz, dairenin ici tek cokgenle kapanir. Donus: (yeni yuz indeksleri, cap ise yalniz daire yuzu;
        kullanilan radius, segs, support) ya da hata metni (kafes degismez)."""
        name = "Circle on face" if cap else "Hole"         # hata metinleri
        if not region:
            return T(name + ": select faces")
        ef = self.edge_faces()
        edges = {e for e, fs in ef.items() if sum(f in region for f in fs) == 1}
        chains = self.edge_chains(edges)
        if not chains or len(chains) != 1 or not chains[0][1]:
            return T(name + ": each selected patch must be closed, without holes inside")
        B = chains[0][0]
        # bolgenin normali ve gercek alan merkezi (ucgen yelpazeleri; kose ortalamasi simetrik bolgede bile
        # merkezi eksenden kaydiriyordu)
        nsum, csum, asum = np.zeros(3), np.zeros(3), 0.0
        for fi in region:
            P = self.V[self.F[fi]]
            nsum += np.cross(P, np.roll(P, -1, 0)).sum(0) / 2
            for k in range(1, len(P) - 1):
                a = np.linalg.norm(np.cross(P[k] - P[0], P[k + 1] - P[0])) / 2
                csum += (P[0] + P[k] + P[k + 1]) / 3 * a
                asum += a
        if asum < 1e-12 or np.linalg.norm(nsum) < 1e-9:
            return T(name + ": the selected faces are not flat")
        n = nsum / np.linalg.norm(nsum)
        c = csum / asum
        if any(self.normal(fi).dot(n) < 0.8 for fi in region):     # ortalamadan ~37 dereceden fazla sapan yuz
            return T(name + ": the selected faces are not flat")
        # sinir sirasi bolge yuzleriyle ayni yonde olsun (yuzde B[0] -> B[1])
        f0 = next(self.F[f] for f in ef[_edge(B[0], B[1])] if f in region)
        if not any(f0[i] == B[0] and f0[(i + 1) % len(f0)] == B[1] for i in range(len(f0))):
            B = B[::-1]
        # daire kafesin ana eksenine hizali: duzleme en yatik X / Y / Z ekseni u, bir kose u uzerinde
        ax = np.eye(3)[int(np.argmin(np.abs(n)))]
        u = ax - ax.dot(n) * n
        u /= np.linalg.norm(u)
        w = np.cross(n, u)
        D = self.V[B] - c
        D2 = np.stack([D @ u, D @ w], 1)
        # merkezden sinira en kisa uzaklik (duzlemde): delik + destek halkasi bunun icinde kalmali
        dmin = min(_seg_dist(0.0, 0.0, D2[i], D2[(i + 1) % len(B)]) for i in range(len(B)))
        N = len(B)
        segs = int(segs) if segs and segs >= 3 else max(8, N + N % 2)
        radius = radius if radius > 0 else 0.5 * dmin
        auto = support < 0                                 # otomatik destek yaricapla buyur: 0.2 * radius
        support = support if support >= 0 else 0.2 * radius
        if radius + support > 0.95 * dmin:
            return T(name + ": radius too large for the selected faces (max %.2f mm)") % (
                0.95 * dmin / 1.2 if auto else 0.95 * dmin - support)
        th = np.unwrap(np.arctan2(D2[:, 1], D2[:, 0]))
        if th[-1] < th[0]:                                 # sinir normal etrafinda saat yonunde: aci azalir
            th = -th
            w = -w
        step = 2 * math.pi / segs
        a0 = round(th[0] / step) * step                    # u'ya hizali kose acilarindan ilk sinir noktasina en yakini
        j = [int(round((t - a0) / step)) for t in th] + [segs]
        j = list(np.maximum.accumulate(np.clip(j, 0, segs)))

        def ring(r):
            """Bolge yuzeyine izdusurulmus daire noktalari (yeni nokta indeksleri)."""
            pts = []
            for k in range(segs):
                a = a0 + k * step
                p = c + r * (math.cos(a) * u + math.sin(a) * w)
                pts.append(self._project_on(region, p, n))
            i0 = len(self.V)
            self.V = np.vstack([self.V, np.array(pts)])
            return list(range(i0, i0 + segs))
        S = ring(radius + support) if support > 0 else None
        H = ring(radius)
        outer = S if S is not None else H
        nf = len(self.F)
        for i in range(N):
            ja, jb = j[i], j[i + 1]
            self.F.append([B[i], B[(i + 1) % N]] + [outer[k % segs] for k in range(jb, ja - 1, -1)])
        if S is not None:
            for k in range(segs):
                k2 = (k + 1) % segs
                self.F.append([S[k], S[k2], H[k2], H[k]])
        if cap:
            self.F.append(list(H))                         # bolgeyle ayni yonde
        self.delete(region)                                # yeni yuzler sonda: indeksleri silinenler kadar kayar
        new = [len(self.F) - 1] if cap else list(range(nf - len(region), len(self.F)))
        return new, radius, segs, support

    def _project_on(self, region, p, n):
        """p noktasini n dogrultusunda bolge yuzeyine (ucgen yelpazeleri) izdusurur; kesmezse p."""
        best = None
        for fi in region:
            f = self.F[fi]
            for k in range(1, len(f) - 1):
                a, b, cc = self.V[f[0]], self.V[f[k]], self.V[f[k + 1]]
                e1, e2 = b - a, cc - a
                q = np.cross(n, e2)
                det = e1.dot(q)
                if abs(det) < 1e-12:
                    continue
                s = p - a
                uu = s.dot(q) / det
                r = np.cross(s, e1)
                vv = n.dot(r) / det
                if uu < -1e-6 or vv < -1e-6 or uu + vv > 1 + 1e-6:
                    continue
                t = e2.dot(r) / det
                if best is None or abs(t) < abs(best):
                    best = t
        return p + n * best if best is not None else p

    def bevel_dirs(self, region):
        """Bevel / Smooth Shift icin, bolgenin her noktasi (sirali, extrude ile ayni sira) icin iki yon:
        kaydirma (komsu secili yuzlerin ortalama normali) ve iceri alma (sinir noktalarinda yuz duzleminde
        iceri, kose acisina gore olceklenir ki kenarlar esit mesafede iceri girsin; ic noktalarda 0)."""
        verts = sorted({v for fi in region for v in self.F[fi]})
        idx = {v: k for k, v in enumerate(verts)}
        count = {}
        for fi in region:
            f = self.F[fi]
            for i in range(len(f)):
                e = _edge(f[i], f[(i + 1) % len(f)])
                count[e] = count.get(e, 0) + 1
        Nv = np.zeros((len(verts), 3))
        Ov = np.zeros((len(verts), 3))
        inward = [[] for _ in verts]
        for fi in region:
            f = self.F[fi]
            n = self.normal(fi)
            for v in f:
                Nv[idx[v]] += n
            for i in range(len(f)):
                a, b = f[i], f[(i + 1) % len(f)]
                if count[_edge(a, b)] != 1:
                    continue
                d = self.V[b] - self.V[a]
                ln = np.linalg.norm(d)
                if ln < 1e-12:
                    continue
                w = np.cross(n, d / ln)                    # saat yonu tersi yuzde kenarin solu = ic
                for v in (a, b):
                    Ov[idx[v]] += w
                    inward[idx[v]].append(w)
        for k in range(len(verts)):
            ln = np.linalg.norm(Nv[k])
            if ln > 1e-12:
                Nv[k] /= ln
            ln = np.linalg.norm(Ov[k])
            if ln > 1e-12:
                m = Ov[k] / ln
                cos = min(m.dot(w) for w in inward[k])
                Ov[k] = m / max(cos, 0.25)
        return Nv, Ov

    def _ring_exit(self, r):
        """Halkanin n-gon'dan (n >= 5) cikis kenari: r[0]r[1] giris kenari; giris ortasindan, yuz duzleminde
        girise dik ve ice dogru bir isin hangi kenara varirsa (r[j], r[j+1]) o. j dondurur (None = cikis yok).
        Birlestirilmis dortgenler (Merge polygons) dogrusal noktali n-gon olur; isin karsi kenardaki dogru
        yere varir (topolojik karsi kenar orada yanlis olurdu)."""
        n = len(r)
        P = self.V[r]
        nrm = np.zeros(3)
        for i in range(n):                                 # Newell normali
            a, b = P[i], P[(i + 1) % n]
            nrm += np.cross(a, b)
        d = P[1] - P[0]
        ln, dl = np.linalg.norm(nrm), np.linalg.norm(d)
        if ln < 1e-12 or dl < 1e-12:
            return n // 2 if n % 2 == 0 else None
        x = d / dl
        y = np.cross(nrm / ln, x)
        m = (P[0] + P[1]) / 2
        if (P.mean(0) - m).dot(y) < 0:
            y = -y
        Q = np.column_stack([(P - m).dot(x), (P - m).dot(y)])   # isin: (0, 0)'dan +y yonunde
        best, bt = None, None
        for j in range(2, n - 1):
            (ax, ay), (bx, by) = Q[j], Q[(j + 1) % n]
            if abs(bx - ax) < 1e-12:
                continue
            s = -ax / (bx - ax)
            if -1e-6 <= s <= 1 + 1e-6:
                t = ay + s * (by - ay)
                if t > 1e-9 and (bt is None or t < bt):
                    best, bt = j, t
        return best

    def face_ring(self, fi, u, v):
        """(u, v) kenarindan gecen yuz halkasi: {yuz: (dondurulmus yuz r, j)}; r[0]r[1] giris, r[j]r[j+1] cikis
        kenari. Dortgende j = 2 (karsi kenar); 5+ kenarli yuzde cikis _ring_exit ile bulunur; ucgende halka biter."""
        ef = self.edge_faces()

        def rotated(f, e):
            for i in range(len(f)):
                if {f[i], f[(i + 1) % len(f)]} == e:
                    return list(f[i:]) + list(f[:i])

        ring = {}

        def walk(f, e):
            while f is not None and f not in ring and len(self.F[f]) >= 4:
                q = rotated(self.F[f], e)
                j = 2 if len(q) == 4 else self._ring_exit(q)
                if j is None:
                    break
                ring[f] = (q, j)
                e = {q[j], q[(j + 1) % len(q)]}
                nb = [g for g in ef[_edge(*e)] if g != f]
                f = nb[0] if nb else None

        walk(fi, {u, v})
        nb = [g for g in ef[_edge(u, v)] if g != fi]
        if nb:
            walk(nb[0], {u, v})
        return ring

    def edge_loop(self, e0):
        """Kenar halkasi: 4 kenarli noktalarda, mevcut kenarin yuzlerini paylasmayan karsi kenarla devam eder.
        Acik (sinir) kenar: sinir boyunca devam eder (delik agzi, acik kenar halkasi; Blender gibi); oradaki
        noktalar 3 kenarli oldugundan eskiden ilk noktada duruyordu."""
        ef = self.edge_faces()
        ve = {}
        for e in ef:
            for v in e:
                ve.setdefault(v, []).append(e)
        loop = {e0}
        for start in e0:
            e, v = e0, start
            while True:
                faces = set(ef[e])
                if len(faces) == 1:                        # sinir: bu noktadaki oteki sinir kenari
                    nxt = [c for c in ve[v] if c != e and len(ef[c]) == 1]
                elif len(ve.get(v, ())) == 4:
                    nxt = [c for c in ve[v] if c != e and not faces & set(ef[c])]
                else:                                      # kutup (3 / 5+ kenarli nokta): halka biter
                    break
                if len(nxt) != 1 or nxt[0] in loop:
                    break
                e = nxt[0]
                loop.add(e)
                v = next(iter(e - {v}))
        return loop

    def loop_cut(self, fi, u, v):
        """(u, v) kenarini kesen, dortgen halka boyunca yeni bir kenar halkasi ekler. Kaydirma icin
        [(yeni nokta, a, b)] dondurur: a hep u ile ayni taraftaki uc (bos liste = kesilemedi)."""
        ring = self.face_ring(fi, u, v)
        if not ring:
            return []
        return self._split_ring(ring, u)

    def _split_ring(self, ring, seed):
        """Halkadaki her yuzu giris / cikis kenarlarinin ortasindan ikiye boler."""
        # taraf: q0 ile q[j+1], q1 ile q[j] ayni tarafta; halka sirasiyla seed'in tarafi yayilir
        side = {seed}
        for q, j in ring.values():
            k = (j + 1) % len(q)
            if q[1] in side or q[j] in side:
                side |= {q[1], q[j]}
            else:
                side |= {q[0], q[k]}
        mids, extra = {}, []

        def mid(a, b):
            k = _edge(a, b)
            if k not in mids:
                mids[k] = len(self.V) + len(extra)
                extra.append((self.V[a] + self.V[b]) / 2)
            return mids[k]

        split = []
        for i, f in enumerate(self.F):
            if i in ring:
                q, j = ring[i]
                m, n = mid(q[0], q[1]), mid(q[j], q[(j + 1) % len(q)])
                split += [[q[0], m, n] + q[j + 1:], [m] + q[1:j + 1] + [n]]
            else:
                split.append(f)
        self.V = np.vstack([self.V, extra])
        # halkanin bittigi dortgen olmayan komsu yuzlere yeni noktalari ekle
        out = []
        for f in split:
            g = []
            for i in range(len(f)):
                g.append(f[i])
                k = _edge(f[i], f[(i + 1) % len(f)])
                if k in mids:
                    g.append(mids[k])
            out.append(g)
        self.F = out
        for k, m in mids.items():                  # bolunen keskin kenarin iki yarisi da keskin
            if k in self.C:
                s = self.C.pop(k)
                a, b = tuple(k)
                self.C[_edge(a, m)] = s
                self.C[_edge(m, b)] = s
        out = []
        for k, m in mids.items():
            a, b = tuple(k)
            out.append((m, a, b) if a in side else (m, b, a))
        return out

    def bevel_edges(self, edges, w, segs=1):
        """Kenar bevel'i (LightWave Rounder / Blender Ctrl+B): secili kenarlar w genisliginde, segs segmentli
        yuvarlak bir seritle degisir. Her uc noktanin yuz yelpazesi, secili kenarlarla dilimlere ayrilir; her
        dilime bir yeni nokta. Yeni nokta dilimin icindeki secili olmayan kenar boyunca kayar (ki secili kenara
        uzakligi w olsun); dilim tek yuzse o yuzde aciortay boyunca. Iki yeni nokta arasindaki profil, kose
        noktasini kontrol noktasi alan dairesel bir egridir (90 derecede tam ceyrek daire). Zincirin ucu 3 kenarli
        bir koseyse kose kalkar ve yan yuz profili alir (Blender gibi); daha cok kenarliysa kose kalir ve profille
        arasi bir yuzle kapanir. 3+ secili kenarin bulustugu noktada bosluk bir kapak yuzuyle kapanir.
        Acik kenarli noktalar desteklenmez (False doner)."""
        segs = max(int(segs), 1)
        ef = self.edge_faces()
        edges = {e for e in edges if len(ef.get(e, ())) == 2}
        if not edges or w <= 0:
            return False
        V, F = self.V, self.F
        verts = {v for e in edges for v in e}
        vf = {}
        for fi, f in enumerate(F):
            for v in f:
                if v in verts:
                    vf.setdefault(v, []).append(fi)
        extra = []

        def add(p):
            extra.append(p)
            return len(V) + len(extra) - 1

        def other(c, v):
            return next(iter(c - {v}))

        def unit(v, x):
            d = V[x] - V[v]
            return d / max(np.linalg.norm(d), 1e-12)

        def slide(v, along, bounds):
            """v'den 'along' kenarlari yonunde, 'bounds' (secili) kenarlarinin dogrularina uzaklik w olacak kadar."""
            u = sum(unit(v, other(c, v)) for c in along)
            u = u / max(np.linalg.norm(u), 1e-12)
            dens = [np.linalg.norm(u - u.dot(s) * s) for s in (unit(v, other(b, v)) for b in bounds)]
            t = w / max(float(np.mean(dens)), 1e-6)
            lim = 0.95 * min(np.linalg.norm(V[other(c, v)] - V[v]) for c in along)
            return V[v] + u * min(t, lim)

        repl, inserts, caps, ends, nsel = {}, [], [], {}, {}
        kept, fans = set(), []                         # koseyi koruyan noktalar ve yelpaze kapaklari
        for v in verts:
            f0, fan, cross = vf[v][0], [], []
            f = f0
            for _ in range(len(vf[v]) + 1):            # yuz yelpazesi: onceki-kenardan komsu yuze gec
                g = F[f]
                prev = g[g.index(v) - 1]
                c = _edge(v, prev)
                fan.append(f)
                cross.append(c)                        # cross[i]: fan[i] ile fan[i+1] arasindaki kenar
                nb = [h for h in ef[c] if h != f]
                if not nb:
                    return False                       # acik kenar
                f = nb[0]
                if f == f0:
                    break
            else:
                return False
            k = len(fan)
            S = [i for i in range(k) if cross[i] in edges]
            nsel[v] = len(S)
            if len(S) >= 2:
                ring = []                              # dilim basina (ilk, son) yeni nokta
                for j, s0 in enumerate(S):
                    s1 = S[(j + 1) % len(S)]
                    n = (s1 - s0) % k or k
                    faces = [fan[(s0 + 1 + q) % k] for q in range(n)]
                    inner = [cross[(s0 + 1 + q) % k] for q in range(n - 1)]
                    if len(inner) >= 2:
                        # dilimde birden cok secili olmayan kenar: tek ortalama nokta kenarlardan kacip sivri
                        # bir uc yapiyordu. Blender gibi: her pahli kenar kendi komsu kenari boyunca kayar,
                        # aradaki yuzler koseyi korur; bosluk koseden yelpazeyle kapanir (kalan kose).
                        p1 = add(slide(v, [inner[0]], [cross[s0]]))
                        p2 = add(slide(v, [inner[-1]], [cross[s1]]))
                        repl[(faces[0], v)] = p1
                        repl[(faces[-1], v)] = p2
                        inserts.append((faces[1], v, other(inner[0], v), p1))
                        inserts.append((faces[-2], v, other(inner[-1], v), p2))
                        ring.append((p1, p2))
                        kept.add(v)
                        continue
                    if inner:
                        p = slide(v, inner, [cross[s0], cross[s1]])
                    else:                              # dilim tek yuz: yuzde aciortay
                        g = F[faces[0]]
                        i = g.index(v)
                        a, b = g[i - 1], g[(i + 1) % len(g)]
                        e1, e2 = unit(v, a), unit(v, b)
                        m = e1 + e2
                        m = m / max(np.linalg.norm(m), 1e-12)
                        t = w / max(np.linalg.norm(m - m.dot(e1) * e1), 1e-6)
                        lim = 0.95 * min(np.linalg.norm(V[a] - V[v]), np.linalg.norm(V[b] - V[v]))
                        p = V[v] + m * min(t, lim)
                    nv = add(p)
                    for fi in faces:
                        repl[(fi, v)] = nv
                    ring.append((nv, nv))
                if v in kept:
                    fans.append((v, ring))
                elif len(S) >= 3:
                    caps.append((v, [a for a, _ in ring]))
            else:                                      # zincirin ucu: nokta kalir, iki yanina nokta acilir
                s = S[0]
                fs, fs1 = fan[s], fan[(s + 1) % k]
                ca, cb = cross[(s - 1) % k], cross[(s + 1) % k]
                v1 = add(slide(v, [ca], [cross[s]]))
                v2 = add(slide(v, [cb], [cross[s]]))
                repl[(fs, v)] = v1
                repl[(fs1, v)] = v2
                inserts.append((fan[(s - 1) % k], v, other(ca, v), v1))
                inserts.append((fan[(s + 2) % k], v, other(cb, v), v2))
                ends[v] = (k, fan[(s - 1) % k])            # 3 kenarli kosede bu yan yuz koseyi kaybeder

        Vn = np.vstack([V, extra]) if extra else V
        # zincir uclari: komsu yuzun kenarina yeni noktalari ekle. Ayni kenara iki uctan nokta gelebilir
        # (iki ayri kenar birlikte pahlaninca): hepsini birlikte, kenar boyunca sirali ekle
        groups = {}
        for (ff, a, b, m) in inserts:
            groups.setdefault((ff, _edge(a, b)), []).append(m)
        out = []
        for fi, f in enumerate(F):
            g = [repl.get((fi, x), x) for x in f]
            for (ff, e), ms in groups.items():
                if ff != fi:
                    continue
                a, b = tuple(e)                        # kenarin bir ucu bu yuzde baska bir bevelle kaymis olabilir
                e2 = _edge(repl.get((fi, a), a), repl.get((fi, b), b))
                for i in range(len(g)):
                    if _edge(g[i], g[(i + 1) % len(g)]) in (e, e2):
                        p0 = Vn[g[i]]
                        ms = sorted(ms, key=lambda m: np.linalg.norm(Vn[m] - p0))
                        g[i + 1:i + 1] = ms
                        break
            out.append(g)
        prof, ctrl = {}, {}

        def pos(i):
            return V[i] if i < len(V) else extra[i - len(V)]

        def profile(v, x, y, e=None):
            """v kosesinde x'ten y'ye profil noktalari [x, ..., y]; ayni uc cifti tum seritlerde paylasilir.
            Agirlikli (rasyonel) ikinci derece Bezier, kontrol noktasi iki teget dogrunun kesisimi; kose
            acisinda tam daire yayi."""
            if (v, x, y) in prof:
                return prof[(v, x, y)]
            if (v, y, x) in prof:
                return prof[(v, y, x)][::-1]
            A, B, C = Vn[x], Vn[y], V[v]
            # iki pahli kenarin bulustugu noktada (zincir ortasi / halka kosesi) x ve y zincirin iki yaninda:
            # tegetler kosede kesisir, kontrol kose olur ve profil gonye duzleminde kalir. Zincir ucunda ve
            # 3+ kenarli kosede x-y tek kenarla ayrilir: tegetler o kenarin dogrusunda kesisir.
            if e is not None and (nsel.get(v) != 2 or v in kept):
                d = V[other(e, v)] - V[v]
                d = d / max(np.linalg.norm(d), 1e-12)
                C = V[v] + d * ((A + B) / 2 - V[v]).dot(d)
            ctrl[(v, x, y)] = C
            pts = [x]
            if segs > 1:
                da, db = A - C, B - C
                cosang = da.dot(db) / max(np.linalg.norm(da) * np.linalg.norm(db), 1e-12)
                ang = math.acos(min(max(cosang, -1.0), 1.0))  # kosedeki aci; yay acisi = pi - ang
                wt = max(math.cos((math.pi - ang) / 2.0), 0.05)
                for q in range(1, segs):
                    t = q / segs
                    b0, b1, b2 = (1 - t) ** 2, 2 * t * (1 - t) * wt, t * t
                    extra.append((b0 * A + b1 * C + b2 * B) / (b0 + b1 + b2))
                    pts.append(len(V) + len(extra) - 1)
            pts.append(y)
            prof[(v, x, y)] = pts
            return pts

        drop = {}                                      # 3 kenarli zincir ucu: yan yuzde kose yerine profil
        for e in edges:                                # serit yuzleri: segs adet dortgen
            a, b = tuple(e)
            f1, f2 = ef[e]
            g1 = F[f1]
            if g1[(g1.index(a) + 1) % len(g1)] != b:   # f1'de a -> b olsun
                f1, f2 = f2, f1
            a1, b1 = repl.get((f1, a), a), repl.get((f1, b), b)
            a2, b2 = repl.get((f2, a), a), repl.get((f2, b), b)
            Pa, Pb = profile(a, a1, a2, e), profile(b, b1, b2, e)
            for q in range(segs):
                out.append([Pa[q], Pa[q + 1], Pb[q + 1], Pb[q]])
            for x, P in ((a, Pa), (b, Pb[::-1])):
                if x not in ends:
                    continue
                k, side = ends[x]
                if k == 3:                             # kose kalkar: yan yuz profili alir
                    drop[(side, x)] = P
                else:                                  # kose kalir: profille arasi koseden yelpaze (Blender gibi);
                    for q in range(len(P) - 1):        # tek buyuk n-gen icbukey olur, cizimde ve SubD'de bozulurdu
                        out.append([x, P[q + 1], P[q]])
        for v, ring in fans:                           # koseyi koruyan nokta: seritlerin uclari arasi koseden yelpaze
            n = len(ring)
            for j in range(n):
                if ring[j][0] == ring[j][1]:
                    continue
                path = [ring[j][1]]                    # korunan dilimin son noktasindan basla,
                q = j
                while True:                            # sonraki seritlerin profilleri boyunca ilerle,
                    nq = (q + 1) % n
                    path += profile(v, ring[q][1], ring[nq][0])[1:]
                    q = nq
                    if ring[q][0] != ring[q][1]:       # bir sonraki korunan dilimin ilk noktasinda dur
                        break
                Vt = np.vstack([V, extra])
                P = Vt[[v] + path]
                nrm = np.cross(P, np.roll(P, -1, 0)).sum(0)
                if nrm.dot(sum(self.normal(fi) for fi in vf[v])) < 0:
                    path = path[::-1]
                for i in range(len(path) - 1):
                    out.append([v, path[i], path[i + 1]])
        for (side, x), P in drop.items():
            g = out[side]
            i = g.index(x)
            inner = P[1:-1]
            g[i:i + 1] = inner[::-1] if g[i - 1] == P[-1] else inner
        corner_n = lambda v: sum(self.normal(fi) for fi in vf[v])
        for v, ring in caps:
            if len(ring) == 3 and segs >= 2 and segs % 2 == 0:
                out += self._corner_patch(v, ring, segs, profile, ctrl, pos, extra, len(V), corner_n(v))
                continue
            # 3+ kenarin bulustugu kose (tek segment ya da tek sayi): kapak yuzu (profillerle)
            poly = []
            for i in range(len(ring)):
                poly += profile(v, ring[i], ring[(i + 1) % len(ring)])[:-1]
            Vt = np.vstack([V, extra])
            P = Vt[poly]
            nrm = np.cross(P, np.roll(P, -1, 0)).sum(0)
            if nrm.dot(sum(self.normal(fi) for fi in vf[v])) < 0:
                poly = poly[::-1]
            out.append(poly)
        Vn = np.vstack([V, extra]) if extra else V
        self.V, self.F = Vn, out
        # kullanilmayan (bolunen) noktalari at, indeksleri sikistir
        used = sorted({x for f in self.F for x in f} | {x for e in self.W for x in e})
        remap = {x: i for i, x in enumerate(used)}
        self.V = self.V[used]
        self.F = [[remap[x] for x in f] for f in self.F]
        self.W = {_edge(*(remap[x] for x in e)) for e in self.W}
        self.dims = [d for d in (_dim_remap(d, remap) for d in self.dims) if d is not None]
        alive = self.face_edges()
        C = {}
        for e, sh in self.C.items():
            if e in edges or not all(x in remap for x in e):
                continue
            ne = _edge(*(remap[x] for x in e))
            if ne in alive:
                C[ne] = sh
        self.C = C
        return True

    @staticmethod
    def _corner_patch(v, ring, segs, profile, ctrl, pos, extra, nV, outward):
        """3 pahli kenarin bulustugu kose icin dortgen yama (Blender 'patch' kosesi): kose ortada bir noktada
        bulusan 3 bolgeye ayrilir, her bolge (segs/2 x segs/2) dortgen izgara. Noktalar kenar yuvarlaklarinin
        olusturdugu kureye oturtulur: yuvarlaklar kosede kesintisiz birlesir. Yuz listesini dondurur."""
        m = segs // 2

        def add(p):
            extra.append(p)
            return nV + len(extra) - 1

        profs = [profile(v, ring[i], ring[(i + 1) % 3]) for i in range(3)]
        # kure: her profil yayinin merkezi (A + B - kontrol); hepsinin ortalamasi, yaricap = sinir ortalamasi
        cents = []
        for i in range(3):
            x, y = ring[i], ring[(i + 1) % 3]
            C = ctrl.get((v, x, y))
            if C is None:
                C = ctrl[(v, y, x)]
            cents.append(pos(x) + pos(y) - C)
        O = np.mean(cents, 0)
        R = float(np.mean([np.linalg.norm(pos(p) - O) for pr in profs for p in pr]))
        proj = lambda p: O + (p - O) * (R / max(np.linalg.norm(p - O), 1e-12))
        Cn = add(proj(np.mean([pos(r) for r in ring], 0)))
        spokes = []                                    # profil ortasindan merkeze; komsu bolgelerle paylasilir
        for i in range(3):
            M = profs[i][m]
            pts = [M]
            for q in range(1, m):
                t = q / m
                pts.append(add(proj((1 - t) * pos(M) + t * pos(Cn))))
            spokes.append(pts + [Cn])
        faces = []
        for i in range(3):
            a = profs[i][:m + 1]                       # ring[i]'den kendi profilinin ortasina
            b = profs[i - 1][::-1][:m + 1]             # ring[i]'den onceki profilin ortasina
            G = [[None] * (m + 1) for _ in range(m + 1)]
            for u in range(m + 1):
                G[u][0] = a[u]
                G[u][m] = spokes[i - 1][u]
            for w in range(m + 1):
                G[0][w] = b[w]
                G[m][w] = spokes[i][w]
            P = lambda u, w: pos(G[u][w])
            for u in range(1, m):
                for w in range(1, m):                  # Coons yamasi, sonra kureye
                    s, t = u / m, w / m
                    p = ((1 - s) * P(0, w) + s * P(m, w) + (1 - t) * P(u, 0) + t * P(u, m)
                         - ((1 - s) * (1 - t) * P(0, 0) + s * (1 - t) * P(m, 0)
                            + (1 - s) * t * P(0, m) + s * t * P(m, m)))
                    G[u][w] = add(proj(p))
            quads = [[G[u][w], G[u + 1][w], G[u + 1][w + 1], G[u][w + 1]] for u in range(m) for w in range(m)]
            q0 = np.array([pos(x) for x in quads[0]])
            nrm = np.cross(q0, np.roll(q0, -1, 0)).sum(0)
            if nrm.dot(outward) < 0:
                quads = [q[::-1] for q in quads]
            faces += quads
        return faces

    def plane_cut(self, n, d, verts=None):
        """n.p = d duzlemiyle keser (LightWave Knife): duzlemin gectigi kenarlar bolunur, duzlemin iki yanina
        tasan yuzler kesim cizgisi boyunca ikiye ayrilir. verts: yalniz butun noktalari bu kumede olan yuzler
        (None = hepsi). Yeni noktalari dondurur (nokta indeksleri sabit kalir, yeniler sona eklenir)."""
        n = np.asarray(n, float)
        eps = 1e-7 * max(np.ptp(self.V, 0).max(), 1.0)
        s = self.V @ n - d
        s[np.abs(s) < eps] = 0.0                       # duzleme cok yakin noktalar duzlemde sayilir
        target = {fi for fi, f in enumerate(self.F) if verts is None or all(v in verts for v in f)}
        mids, extra = {}, []
        for fi in target:
            f = self.F[fi]
            for i in range(len(f)):
                a, b = f[i], f[(i + 1) % len(f)]
                k = _edge(a, b)
                if s[a] * s[b] < 0 and k not in mids:
                    t = s[a] / (s[a] - s[b])
                    mids[k] = len(self.V) + len(extra)
                    extra.append(self.V[a] + t * (self.V[b] - self.V[a]))
        if extra:
            self.V = np.vstack([self.V, extra])
        on = {int(v) for v in np.nonzero(s == 0)[0]} | set(mids.values())
        out = []
        for fi, f in enumerate(self.F):
            g = []
            for i in range(len(f)):
                g.append(f[i])
                k = _edge(f[i], f[(i + 1) % len(f)])
                if k in mids:
                    g.append(mids[k])
            if fi in target and (s[f] > 0).any() and (s[f] < 0).any():
                cut = [j for j, v in enumerate(g) if v in on]
                if len(cut) == 2 and (cut[1] - cut[0]) % len(g) not in (1, len(g) - 1):
                    i0, i1 = cut                       # kesim kirisi boyunca iki cokgen (yon korunur)
                    out.append(g[i0:i1 + 1])
                    out.append(g[i1:] + g[:i0 + 1])
                    continue
            out.append(g)
        self.F = out
        for k, m in mids.items():                      # bolunen keskin kenarin iki yarisi da keskin
            if k in self.C:
                sh = self.C.pop(k)
                a, b = tuple(k)
                self.C[_edge(a, m)] = sh
                self.C[_edge(m, b)] = sh
        return list(mids.values())

    def split_plane(self, n, d, verts=None):
        """plane_cut'tan sonra modeli kesim cizgisinden ayirir. Kesilen yuzlerin (verts: yalniz butun noktalari bu
        kumede olanlar, None = hepsi) pozitif yanindakilerinden baslayip kesim cizgisini gecmeden ulasilan BUTUN
        yuzler bir yan olur (secimle kesilse de model ikiye bolunur). Kesim cizgisindeki noktalar o yan icin
        kopyalanir; iki parca da kesim boyunca acik kenarli olur. Donus: (bir yanin yuzleri, ayrildi mi). Kesim
        modelin cevresini tam dolasmiyorsa ayrilamaz (False): o zaman yalniz kesim boyunca bir yarik acilir."""
        n = np.asarray(n, float)
        eps = 1e-7 * max(np.ptp(self.V, 0).max(), 1.0)
        sd = self.V @ n - d
        sd[np.abs(sd) < eps] = 0.0
        target = [fi for fi, f in enumerate(self.F) if verts is None or all(v in verts for v in f)]
        seeds = [fi for fi in target if (sd[self.F[fi]] > 0).any() and not (sd[self.F[fi]] < 0).any()]
        negs = {fi for fi in target if (sd[self.F[fi]] < 0).any()}
        cut = {_edge(f[i], f[(i + 1) % len(f)]) for fi in target for f in [self.F[fi]] for i in range(len(f))
               if sd[f[i]] == 0 and sd[f[(i + 1) % len(f)]] == 0}
        ef = self.edge_faces()
        side, stack = set(seeds), list(seeds)
        while stack:                                       # kesim kenarlarini gecmeden yayil
            f = self.F[stack.pop()]
            for i in range(len(f)):
                e = _edge(f[i], f[(i + 1) % len(f)])
                if e in cut:
                    continue
                for g in ef[e]:
                    if g not in side:
                        side.add(g)
                        stack.append(g)
        closed = bool(side) and not (side & negs) and len(side) < len(self.F)
        if not closed:
            side = set(seeds)                              # ayrilamiyor: yalniz kesilen yuzler (yarik)
        if not side:
            return [], False
        in_side = {v for fi in side for v in self.F[fi]}
        in_other = {v for fi, f in enumerate(self.F) if fi not in side for v in f}
        split = sorted(v for v in in_side & in_other if sd[v] == 0)
        if split:
            dup = {v: len(self.V) + i for i, v in enumerate(split)}
            self.V = np.vstack([self.V, self.V[split]])
            e_side = {_edge(f[i], f[(i + 1) % len(f)]) for fi in side for f in [self.F[fi]] for i in range(len(f))}
            e_other = {_edge(f[i], f[(i + 1) % len(f)]) for fi, f in enumerate(self.F) if fi not in side
                       for i in range(len(f))}
            for e, sh in list(self.C.items()):             # o yandaki keskin kenarlar kopyalara gecer
                if e in e_side:
                    self.C[_edge(*(dup.get(v, v) for v in e))] = sh
                    if e not in e_other:
                        del self.C[e]
            for fi in side:
                self.F[fi] = [dup.get(v, v) for v in self.F[fi]]
        return sorted(side), closed

    def split_edges(self, edges, count=1, at=None):
        """Kenarlara nokta ekler (LightWave Add Points): edges'teki her kenar count + 1 esit parcaya bolunur; at
        {kenar: (a, t)}: o kenara a ucundan t oraninda tek nokta. Kenari paylasan butun cokgenler yeni noktalari
        alir (bosluk kalmaz), keskin kenarin parcalari da keskin. Nokta indeksleri sabit, yeniler sona eklenir.
        Dondurur: {kenar: (a, [a'dan diger uca yeni noktalar])}."""
        spec = {e: (min(e), [k / (count + 1.0) for k in range(1, count + 1)]) for e in edges}
        for e, (a, t) in (at or {}).items():
            spec[e] = (a, [t])
        new, extra = {}, []
        for e, (a, ts) in spec.items():
            b = next(iter(e - {a}))
            ids = list(range(len(self.V) + len(extra), len(self.V) + len(extra) + len(ts)))
            extra += [self.V[a] + t * (self.V[b] - self.V[a]) for t in ts]
            new[e] = (a, ids)
        if not extra:
            return {}
        self.V = np.vstack([self.V, extra])
        out = []
        for f in self.F:
            g = []
            for i in range(len(f)):
                u, w = f[i], f[(i + 1) % len(f)]
                g.append(u)
                s = new.get(_edge(u, w))
                if s:
                    g += s[1] if u == s[0] else s[1][::-1]
            out.append(g)
        self.F = out
        for e, (a, ids) in new.items():
            if e in self.C:
                sh = self.C.pop(e)
                chain = [a] + ids + [next(iter(e - {a}))]
                for p, q in zip(chain, chain[1:]):
                    self.C[_edge(p, q)] = sh
        return new

    def connect(self, pairs):
        """Ayni cokgendeki iki nokta arasina kenar (LightWave Split): cokgen ikiye bolunur, yon korunur. pairs
        sirayla uygulanir; komsu ya da ayni cokgende olmayan cift atlanir. Dondurur: eklenen kenarlar."""
        done = []
        for u, v in pairs:
            for fi, f in enumerate(self.F):
                if u not in f or v not in f:
                    continue
                i, j = sorted((f.index(u), f.index(v)))
                if j - i in (0, 1, len(f) - 1):
                    continue
                self.F[fi] = f[i:j + 1]
                self.F.append(f[j:] + f[:i + 1])
                done.append(_edge(u, v))
                break
        return done

    def subdivide_faces(self, faces):
        """Cokgenleri dortgenlere boler (LightWave Subdivide, faceted): her kenarin ortasina ve cokgenin ortasina
        nokta; n kenarli cokgen n dortgen olur. Komsu cokgenler kenar ortalarini alir. Dondurur: yeni cokgenler."""
        faces = sorted(set(faces))
        old = {fi: list(self.F[fi]) for fi in faces}
        E = {_edge(f[i], f[(i + 1) % len(f)]) for f in old.values() for i in range(len(f))}
        mid = {e: ids[0] for e, (a, ids) in self.split_edges(E, 1).items()}
        centers = [self.V[f].mean(0) for f in old.values()]
        c0 = len(self.V)
        self.V = np.vstack([self.V, centers])
        out = []
        for k, (fi, f) in enumerate(old.items()):
            n = len(f)
            quads = [[f[i], mid[_edge(f[i], f[(i + 1) % n])], c0 + k, mid[_edge(f[i - 1], f[i])]] for i in range(n)]
            self.F[fi] = quads[0]
            out.append(fi)
            for q in quads[1:]:
                out.append(len(self.F))
                self.F.append(q)
        return out

    def limit_frame(self, rules=RULES, levels=2):
        """Kafes noktalarinin SubD limit yuzeyindeki yerleri ve birim normalleri (levels adim bolup limit
        kurali; Catmull-Clark'ta bir noktanin cocugu hep ayni indekstedir)."""
        N = len(self.V)
        V, F, C = self.V, self.F, self.C
        for _ in range(levels):
            V, F, C = catmull_clark(V, F, C, rules)
        L = _limit(V, F, C, rules)
        P = V[F]
        fn = np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 1])
        Nv = np.zeros_like(V)
        for j in range(4):
            np.add.at(Nv, F[:, j], fn)
        ln = np.linalg.norm(Nv, axis=1, keepdims=True)
        Nv = Nv / np.where(ln > 1e-12, ln, 1.0)
        return L[:N], Nv[:N], (V, F, C, L)

    def shell(self, region, thickness, mode=0, radius=0.0, fill=True, sharp=True, rules=RULES, flat=True):
        """Kabuk (Blender Solidify gibi). region: kabuklanacak yuzler. Bolgenin kopyasi ters cevrilir ve her nokta
        komsu yuzlerin acıyla agirlikli kafes normali boyunca kaydirilir (mode 0: iceri, 1: disari, 2: ortada = asil
        yuzey yarim disari, kopya yarim iceri). Esit kalinlik: koselerde mesafe yuzlerin kalinligi tutacagi kadar
        artar (en fazla 2 kat). Kafes ofseti kullanilir: limit yuzeyine oturtma turlari keskin / dar koselerde kafes
        noktalarini dis kafesin ote yuzune itiyordu.
        radius > 0: akilli suzgec (radius: atlanacak detay genisligi, mm). Her noktanin cevresindeki limit yuzeyine
        oturtulan kuadrikten sapmasi (dar oluk, kucuk tumsek) kopyada duzeltilir; ana formun egriligi kalir. Modelin
        acik kenarlari (delik agzi, acik alt kenar) suzgecten muaftir.
        Kesisme korumasi: kopyanin bir kenari dis kafesi / kendini deliyorsa (ya da dis kafesin kenari kopyayi)
        o noktalarin kaydirmasi adim adim kisaltilir.
        fill: modelin acik kenarlarinda ic ve dis kabuk dortgen seritle birlesir (sharp: serit kenarlari keskin).
        flat: acik kenar noktalari duvar boyunca degil yalniz duvara dik kaydirilir (tam kalinlik): Fill seridi
        duvara dik ve duz olur, ic ve dis kenar yandan bakinca ayni hizada kalir (egri yuzeydeki delikte de).
        Secim siniri (secili / secili olmayan yuz arasi) her zaman acik kalir: orayi kullanici orer (P, Bridge).
        Donus: (kopya yuzleri, kopyanin secim sinirindaki acik kenarlari, kaydirmasi kisaltilan nokta sayisi) ya
        da hata."""
        region = sorted(set(region))
        if not region:
            return T("Shell: nothing to shell")
        if thickness <= 0:
            return T("Shell: thickness must be greater than 0")
        N = len(self.V)
        ef = self.edge_faces()
        rset = set(region)
        R = sorted({v for fi in region for v in self.F[fi]})
        # asiri kalinlik (or. yanlis yazilan deger): parcanin en kucuk boyutunun dortte birinden kalin kabuk kafesi
        # birbirine gecirir. Duz bolgede (secili duz yuzler) sifira yakin boyut sayilmaz.
        ext = np.ptp(self.V[R], 0)
        lim = 0.25 * ext[ext > 0.1 * ext.max()].min()
        if thickness > lim:
            return T("Shell: thickness %.2f mm is too large for this part (max %.2f mm)") % (thickness, lim)
        ri = {v: k for k, v in enumerate(R)}
        open_e = [e for e, fs in ef.items() if len(fs) == 1 and fs[0] in rset]          # modelin acik kenarlari
        border = [e for e, fs in ef.items() if len(fs) == 2 and (fs[0] in rset) != (fs[1] in rset)]
        # acik kenar noktalari suzgecten muaf; bir halka icerisi yarim
        pinned = {v for e in open_e for v in e}
        ring1 = {w for fi in region for f in [self.F[fi]] if pinned & set(f) for w in f} - pinned
        pin = np.array([0.0 if v in pinned else 0.5 if v in ring1 else 1.0 for v in R])[:, None]
        # kafes normalleri (yuz normalleri, kosedeki aciyla agirlikli) ve esit kalinlik carpani
        Nv, cw, ws = np.zeros((len(R), 3)), np.zeros(len(R)), np.zeros(len(R))
        fnorm = {fi: self.normal(fi) for fi in region}
        for fi in region:
            f = self.F[fi]
            for i, v in enumerate(f):
                e1, e2 = self.V[f[i - 1]] - self.V[v], self.V[f[(i + 1) % len(f)]] - self.V[v]
                ang = math.acos(np.clip(e1 @ e2 / max(np.linalg.norm(e1) * np.linalg.norm(e2), 1e-12), -1, 1))
                Nv[ri[v]] += ang * fnorm[fi]
        Nv /= np.maximum(np.linalg.norm(Nv, axis=1, keepdims=True), 1e-12)
        for fi in region:
            f = self.F[fi]
            for i, v in enumerate(f):
                e1, e2 = self.V[f[i - 1]] - self.V[v], self.V[f[(i + 1) % len(f)]] - self.V[v]
                ang = math.acos(np.clip(e1 @ e2 / max(np.linalg.norm(e1) * np.linalg.norm(e2), 1e-12), -1, 1))
                cw[ri[v]] += ang * (Nv[ri[v]] @ fnorm[fi])
                ws[ri[v]] += ang
        k_even = 1.0 / np.clip(cw / np.maximum(ws, 1e-12), 0.5, 1.0)
        delta = np.zeros((len(R), 3))                      # suzgecin kopyaya duzeltmesi
        if radius > 1e-9:
            L, Nrm, _ = self.limit_frame(rules)
            Lr, Nr = L[R], Nrm[R]
            Ls = Lr.copy()
            # kafes seyrek olabilir (oluk noktalari radius icinde yalniz birbirini gorur): yuzeyden yogun ornek al,
            # her noktanin radius cevresindeki yuzeye agirlikli duzlem oturt, nokta bu duzleme gore yerlessin
            S, Sn = self._region_samples(region, rules, self._sample_levels(region, radius), normals=True)
            grid = {}
            for k, c in enumerate(map(tuple, np.floor(S / (3 * radius)).astype(int))):
                grid.setdefault(c, []).append(k)
            for k, (p, n) in enumerate(zip(Lr, Nr)):
                c = np.floor(p / (3 * radius)).astype(int)
                cand = [j for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
                        for j in grid.get((c[0] + dx, c[1] + dy, c[2] + dz), ())]
                if not cand:
                    continue
                P, Q = S[cand], Sn[cand]
                d = np.linalg.norm(P - p, axis=1)
                w0 = np.exp(-(d / (1.5 * radius)) ** 2) * (d < 3 * radius) * (Q @ n > 0.5)   # ote yuz degil
                if w0.sum() < 1e-12:
                    continue
                # yerel cerceve (n normal): z = ax^2 + bxy + cy^2 + dx + ey + f oturt. Kuadrik ana formun egriligini
                # tasir (yuvarlatilmis kenar duzlesip duvar kalinlasmaz / incelmez); aykiri ornekler (dar oluk dibi,
                # kucuk cikinti) Tukey agirligiyla ayiklanir. Sonuc: noktanin kuadrik uzerindeki yeri ve normali.
                u = np.cross(n, [1.0, 0, 0] if abs(n[0]) < 0.9 else [0, 1.0, 0])
                u /= np.linalg.norm(u)
                v = np.cross(n, u)
                D3 = P - p
                X, Y, Z = D3 @ u, D3 @ v, D3 @ n
                A = np.stack([X * X, X * Y, Y * Y, X, Y, np.ones_like(X)], 1)
                w = w0
                coef = None
                for _ in range(5):
                    sw = np.sqrt(w)
                    if (w > 1e-6).sum() < 8:
                        break
                    coef = np.linalg.lstsq(A * sw[:, None], Z * sw, rcond=None)[0]
                    res = np.abs(Z - A @ coef)
                    sc = 1.4826 * np.median(res[w0 > 0.05 * w0.max()]) + 1e-4 * radius
                    w = w0 * np.clip(1 - (res / (3 * sc)) ** 2, 0, 1) ** 2   # Tukey agirligi
                if coef is None:
                    continue
                Ls[k] = p + coef[5] * n                    # kuadrigin noktadaki yeri: kucuk cukur / tumsek atlanir
            delta = pin * (Ls - Lr)
            dn = (delta * Nv).sum(1, keepdims=True)        # duzeltme yalniz normal yonunde: kafes kaymasin
            delta = np.clip(dn, -0.5 * thickness, 0.5 * thickness) * Nv
        t = thickness / 2.0 if mode == 2 else thickness
        sgn = 1.0 if mode == 1 else -1.0
        off = sgn * t * k_even[:, None] * Nv               # kopyanin kaydirmasi
        out = t * k_even[:, None] * Nv                     # ortada: asil yuzeyin kaydirmasi
        if flat and open_e:
            # duvar yonu d: acik kenar noktasindan duvara (acik olmayan kenarlarla komsulara) giden ortalama yon.
            # Once kaydirmanin d bileseni atilir (kenar seridi duvara dik); halka duzlemselse ic kenar noktasi ic duvar
            # boyunca halkanin duzlemine kaydirilir: egimli duvarda da ic ve dis kenar ayni hizada (ayni duzlemde)
            oset = set(open_e)
            walld = {}
            for v in pinned:
                nbs = {w for w in {x for fi in region for f in [self.F[fi]] if v in f
                                   for x in (f[f.index(v) - 1], f[(f.index(v) + 1) % len(f)])}
                       if _edge(v, w) not in oset}
                if not nbs:
                    continue
                d = sum((self.V[w] - self.V[v]) / max(np.linalg.norm(self.V[w] - self.V[v]), 1e-12) for w in nbs)
                if np.linalg.norm(d) < 1e-9:
                    continue
                d /= np.linalg.norm(d)
                k = ri[v]
                walld[v] = d
                n_perp = Nv[k] - (Nv[k] @ d) * d
                if np.linalg.norm(n_perp) > 1e-6:
                    off[k] = sgn * t * n_perp / np.linalg.norm(n_perp)
                    out[k] = t * n_perp / np.linalg.norm(n_perp)
            for chain, closed in self.edge_chains(oset):
                P = self.V[chain]
                if len(chain) < 3:
                    continue
                ctr = P.mean(0)
                pn = np.linalg.svd(P - ctr)[2][2]          # halkanin en uygun duzlemi
                if np.abs((P - ctr) @ pn).max() > 0.05 * max(np.ptp(P, 0).max(), 1e-9):
                    continue                               # duzlemsel degil: duvara dik kenar kalsin
                for v in chain:
                    d = walld.get(v)
                    if d is None or abs(d @ pn) < 0.2:
                        continue
                    k = ri[v]
                    for o in (off, out):                   # ic duvar boyunca duzleme (ortada dis duvar da)
                        sl = -((self.V[v] + o[k] - ctr) @ pn) / (d @ pn)
                        if abs(sl) <= 2 * t:
                            o[k] = o[k] + sl * d
        base = self.V[R].copy()
        if mode == 2:                                      # ortada: asil yuzey de yarim disari
            self.V[R] = base + out
        # kopyada kisa kenarlari birlestir (akilli suzgec): destek halkasi seritleri, kalinliktan dar basamaklar ic
        # tarafta kopyalanmaz. Yoksa keskin katlanmadaki kisa kenarlar kalinlik kadar kayinca kafesler kesisir.
        # Acik kenar (delik agzi, acik alt kenar) ve secim siniri noktalari birlesmez: Fill ve orme bozulmasin.
        pos0, off0 = base + delta, off
        bverts = {x for e in border for x in e}
        fixed = np.array([v in pinned or v in bverts for v in R])
        thr = max(0.5 * radius, thickness) if radius > 1e-9 else 0.0
        P0 = pos0 + off0
        cand = sorted((np.linalg.norm(P0[ri[x]] - P0[ri[y]]), ri[x], ri[y]) for x, y in (tuple(e) for e in redges_of(self, region))) if thr > 0 else []

        def clusters(keep):
            """Her nokta icin kume temsilcisi; keep'teki noktalar birlesmez."""
            par = list(range(len(R)))

            def find(x):
                while par[x] != x:
                    par[x] = par[par[x]]
                    x = par[x]
                return x
            lo, hi = P0.copy(), P0.copy()
            for ln, x, y in cand:
                if ln >= thr:
                    break
                if x in keep or y in keep:
                    continue
                rx, ry = find(x), find(y)
                if rx == ry or (fixed[rx] and fixed[ry]):
                    continue
                l2, h2 = np.minimum(lo[rx], lo[ry]), np.maximum(hi[rx], hi[ry])
                if np.linalg.norm(h2 - l2) > thr:          # kume buyumesin (sik bolge bir noktaya cokmesin)
                    continue
                if fixed[ry]:
                    rx, ry = ry, rx
                par[ry] = rx                               # sabit nokta kume temsilcisi olur
                lo[rx], hi[rx] = l2, h2
            return [find(k) for k in range(len(R))]
        keep = set()
        while True:                                        # birlesme bir kenari ikiden cok yuze ortak ediyorsa o kumeler
            root = clusters(keep)                          # birlesmesin: kafes manifold kalsin (Solid gecerli cikar)
            count = {}
            for fi in region:
                g = []                                     # asagidaki kopya yuzle ayni: ust uste gelenler tek nokta,
                for v in self.F[fi]:                       # ikiden az kenari kalan yuz yok sayilir
                    w = root[ri[v]]
                    if not g or g[-1] != w:
                        g.append(w)
                if len(g) > 1 and g[0] == g[-1]:
                    g.pop()
                if len(g) >= 3:
                    for a, b in zip(g, g[1:] + g[:1]):
                        count[_edge(a, b)] = count.get(_edge(a, b), 0) + 1
            over = {r for e, n in count.items() if n > 2 for r in e}
            add = {k for k in range(len(R)) if root[k] in over and k not in keep}
            if not add:
                break
            keep |= add
        roots = sorted(set(root))
        rid = {r: j for j, r in enumerate(roots)}
        grp = np.array([rid[r] for r in root])
        cnt = np.bincount(grp, minlength=len(roots)).astype(float)[:, None]
        cpos, coff = np.zeros((len(roots), 3)), np.zeros((len(roots), 3))
        np.add.at(cpos, grp, pos0)
        np.add.at(coff, grp, off0)
        cpos, coff = cpos / cnt, coff / cnt
        for r in roots:                                    # sabit uyeli kume: onun yerinde
            if fixed[r]:
                cpos[rid[r]], coff[rid[r]] = pos0[r], off0[r]
        cp = {v: N + int(grp[ri[v]]) for v in R}
        self.V = np.vstack([self.V, cpos + coff])
        nf = len(self.F)
        for fi in region:
            g = []
            for v in self.F[fi][::-1]:
                w = cp[v]
                if not g or g[-1] != w:
                    g.append(w)
            if len(g) > 1 and g[0] == g[-1]:
                g.pop()
            if len(g) >= 3:
                self.F.append(g)
        inner = list(range(nf, len(self.F)))
        redges = redges_of(self, region)
        for e, sh in list(self.C.items()):
            if e in redges:
                a2, b2 = (cp[v] for v in e)
                if a2 != b2:
                    self.C[_edge(a2, b2)] = sh
        if fill:
            for e in open_e:
                f = self.F[ef[e][0]]
                a, b = next((f[i], f[(i + 1) % len(f)]) for i in range(len(f)) if _edge(f[i], f[(i + 1) % len(f)]) == e)
                strip = [b, a, cp[a], cp[b]] if cp[a] != cp[b] else [b, a, cp[a]]
                self.F.append(strip)                       # dis yuzde a -> b: serit b -> a yonunde
                if sharp:
                    self.C[_edge(a, b)] = CREASE_HARD
                    if cp[a] != cp[b]:
                        self.C[_edge(cp[a], cp[b])] = CREASE_HARD
        idx = np.arange(N, N + len(roots))
        thin = self._shell_guard(region, inner, idx, cpos, coff)
        thin += len(R) - len(roots)
        open_inner = [_edge(cp[a], cp[b]) for a, b in (tuple(e) for e in border) if cp[a] != cp[b]]
        if not fill:
            open_inner += [_edge(cp[a], cp[b]) for a, b in (tuple(e) for e in open_e) if cp[a] != cp[b]]
        return inner, open_inner, thin

    def _shell_guard(self, region, inner, idx, base, off, rounds=12):
        """Kopyanin kenari dis kafesi ya da kendini (komsu olmayan yuzlerini) deliyorsa, ya da dis kafesin kenari
        kopyayi deliyorsa, ilgili kopya noktalarinin kaydirmasini kisalt (her turda %70, en az %15).
        Kisaltilan nokta sayisini dondurur."""
        fac = np.ones(len(idx))
        pos = {int(v): k for k, v in enumerate(idx)}

        def tris(faces, copy=False):
            # kopya yuzler ters sirali: dis yuzle ayni kosegenden bolunsun diye son noktadan (dis f[0]'in kopyasi)
            # baslat; yoksa duzlemsel olmayan dortgende iki kafes sifira yakin ofsette bile "kesisir"
            T = []
            for fi in faces:
                f = self.F[fi]
                if len(f) > 4:                         # icbukey n-gen yelpazesi sahte kesisme verir: atla
                    continue
                if copy:
                    f = f[-1:] + f[:-1]
                for k in range(1, len(f) - 1):
                    T.append((f[0], f[k], f[k + 1]))
            return np.array(T, int)

        def edges(faces):
            return np.array(sorted({tuple(sorted((f[i], f[(i + 1) % len(f)]))) for fi in faces for f in [self.F[fi]]
                                    for i in range(len(f))}), int)
        To, Ti = tris(region), tris(inner, True)
        Eo, Ei = edges(region), edges(inner)
        touched = set()
        for _ in range(rounds):
            bad = set()
            for E, Tt, own in ((Ei, To, False), (Eo, Ti, True), (Ei, Ti, None)):
                hit = _seg_tri_hits(self.V, E, Tt)
                for (a, b), tri in hit:
                    if own is False:                   # kopya kenari dis kafeste
                        bad |= {a, b}
                    elif own:                          # dis kenar kopya ucgeninde
                        bad |= set(tri)
                    else:                              # kopya kendini kesiyor
                        bad |= {a, b} | set(tri)
            bad = [pos[v] for v in bad if v in pos]
            if not bad:
                break
            fac[bad] = np.maximum(fac[bad] * 0.7, 0.15)
            touched |= set(bad)
            self.V[idx] = base + off * fac[:, None]
        return len(touched)

    def _sample_levels(self, region, radius):
        """Ornek araligi suzgec yaricapinin yarisindan kucuk olsun diye bolme seviyesi (2-4)."""
        ln = np.mean([np.linalg.norm(self.V[f[i]] - self.V[f[i - 1]]) for fi in region for f in [self.F[fi]]
                      for i in range(len(f))])
        return int(min(max(np.ceil(np.log2(max(ln / (0.5 * radius), 1.0))), 2), 4))

    def _region_samples(self, region, rules, levels=2, normals=False):
        """Bolge yuzlerinin limit yuzeyinden yogun ornek noktalari (levels adim bolunmus kafesin limit yerleri);
        normals: noktalarin birim normalleri de."""
        V, F, C = self.V, self.F, self.C
        deg = np.array([len(f) for f in F])
        start = np.concatenate([[0], np.cumsum(deg)[:-1]])
        ids = np.concatenate([start[fi] + np.arange(deg[fi]) for fi in region])   # 1. seviye cocuklari
        V, F, C = catmull_clark(V, F, C, rules)
        for _ in range(levels - 1):
            V, F, C = catmull_clark(V, F, C, rules)
            ids = (4 * ids[:, None] + np.arange(4)).ravel()
        L = _limit(V, F, C, rules)
        pts = np.unique(F[ids])
        if not normals:
            return L[pts]
        P = V[F]
        fn = np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 1])
        Nv = np.zeros_like(V)
        for j in range(4):
            np.add.at(Nv, F[:, j], fn)
        Nv = Nv[pts] / np.maximum(np.linalg.norm(Nv[pts], axis=1, keepdims=True), 1e-12)
        return L[pts], Nv

    def subset(self, sel, wires=(), points=()):
        """Secili yuzlerden yeni (bagimsiz) kafes: kullanilan noktalar yeniden numaralanir, iki ucu da kalan
        keskin kenarlar keskin kalir. wires: serbest kenar olarak alinacak kenarlar, points: serbest noktalar.
        Kopyala / yapistir icin."""
        faces = [self.F[fi] for fi in sorted(sel)]
        wires = {e for e in wires if len(e) == 2}
        used = sorted({v for f in faces for v in f} | {v for e in wires for v in e} | set(points))
        remap = {v: i for i, v in enumerate(used)}
        F = [[remap[v] for v in f] for f in faces]
        sub = Cage(self.V[used].copy() if used else np.zeros((0, 3)), F,
                   W={_edge(*(remap[v] for v in e)) for e in wires})
        alive = sub.face_edges()
        for e, s in self.C.items():
            if all(v in remap for v in e):
                me = _edge(*(remap[v] for v in e))
                if me in alive:
                    sub.C[me] = s
        return sub

    def mirror(self, sel, Vm, keep):
        """sel yuzlerinin aynasini ayri parca olarak ekler. Vm: her noktanin aynadaki konumu; keep[v]: nokta ayna
        duzleminde -> yeni parca onu paylasir (kaynak; iki yari birlesir). Butun noktalari duzlemde olan yuz
        (kesim duvari) silinir: iki yari birlesince ic duvar olur, kenarlari uc yuzlu (manifold disi) kalirdi.
        Aynalama yonu cevirdiginden yeni yuzler ters sirali. Yeni yuz indeksleri."""
        walls = {fi for fi in sel if all(keep[v] for v in self.F[fi])}
        faces = [self.F[fi] for fi in sorted(sel) if fi not in walls]
        used = sorted({v for f in faces for v in f})
        add = [v for v in used if not keep[v]]
        n0, nf = len(self.V), len(self.F)
        if add:
            self.V = np.vstack([self.V, Vm[add]])
        new = {v: n0 + i for i, v in enumerate(add)}
        new.update({v: v for v in used if keep[v]})
        self.F += [[new[v] for v in f[::-1]] for f in faces]
        alive = {_edge(f[i], f[(i + 1) % len(f)]) for f in faces for i in range(len(f))}
        for e, s in list(self.C.items()):                  # keskin kenarlar aynada da keskin
            if e in alive:
                self.C[_edge(*(new[v] for v in e))] = s
        if walls:
            self.delete(walls)                             # yeni yuzler sonda: indeksleri silinenler kadar kayar
        return list(range(nf - len(walls), len(self.F)))

    def append(self, other):
        """Baska bir kafesi (ayri parca olarak) ekler; eklenen yuzlerin indekslerini dondurur."""
        n, nf = len(self.V), len(self.F)
        self.V = np.vstack([self.V, other.V]) if len(self.V) else other.V.copy()
        self.F += [[v + n for v in f] for f in other.F]
        for e, s in other.C.items():
            self.C[_edge(*(v + n for v in e))] = s
        self.W |= {_edge(*(v + n for v in e)) for e in other.W}
        return list(range(nf, len(self.F)))

    def reduce_points(self, verts, degrees):
        """LightWave Reduce Points: cokgen kenarlari boyunca neredeyse duz giden ara noktalari (gelen ve giden
        kenar arasindaki sapma < degrees) cokgenlerden cikarir. Yalniz iki kenarli noktalar (bir kenarin
        ortasindaki nokta) alinir; boylece komsu cokgenler ayni kenari paylasmaya devam eder (catlak yok). Ucgen
        altina inen cokgen olmaz. Kenar keskinligi yeni kenara gecer. Cikarilan nokta sayisini dondurur."""
        limit = math.cos(math.radians(max(degrees, 0.0)))
        vf = {}
        for fi, f in enumerate(self.F):
            for v in f:
                vf.setdefault(v, []).append(fi)
        removed = 0
        for v in sorted(verts):
            faces = vf.get(v, [])
            if not faces:
                continue
            nb, ok = set(), True
            for fi in faces:
                f = self.F[fi]
                if len(f) <= 3 or f.count(v) != 1:
                    ok = False
                    break
                i = f.index(v)
                nb.update((f[i - 1], f[(i + 1) % len(f)]))
            if not ok or len(nb) != 2:
                continue
            a, b = nb
            d1, d2 = self.V[v] - self.V[a], self.V[b] - self.V[v]
            l1, l2 = np.linalg.norm(d1), np.linalg.norm(d2)
            if l1 < 1e-12 or l2 < 1e-12 or d1.dot(d2) / (l1 * l2) < limit:
                continue
            for fi in faces:
                self.F[fi].remove(v)
            ea, eb = _edge(a, v), _edge(v, b)
            if ea in self.C and eb in self.C:              # iki yarisi da keskinse birlesik kenar da keskin
                self.C[_edge(a, b)] = max(self.C[ea], self.C[eb])
            self.C.pop(ea, None)
            self.C.pop(eb, None)
            vf[v] = []
            removed += 1
        return removed

    def merge_points(self, verts, dist):
        """LightWave Merge Points: verts icinde birbirine dist mesafesindeki noktalar birlesir. Yuksek numarali nokta
        dusuk numaralinin yerine gider (ortalama alinmaz, cokgen esner). Zincirleme birlesme yok: mesafe hayatta kalan
        noktaya olculur. Tek noktaya / kenara inen cokgen silinir (kafeste 1 noktali cokgen yok); kendine degen cokgen
        ikiye ayrilir. Birlesince ayni noktalara oturan cokgenler: zit yonlu cift ic duvardir, ikisi de silinir (iki
        kapali parca tek parca olur, bkz. mirror); ayni yonlulerden biri kalir. Birlesen nokta sayisini dondurur."""
        verts = sorted(verts)
        d = max(float(dist), 1e-9)
        keys = [tuple(c) for c in np.floor(self.V[verts] / d).astype(np.int64).tolist()]   # 3B izgara: yalniz komsular
        cells = {}
        for v, c in zip(verts, keys):
            cells.setdefault(c, []).append(v)
        to = {}
        for v, (x, y, z) in zip(verts, keys):
            if v in to:
                continue
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for w in cells.get((x + dx, y + dy, z + dz), ()):
                            if w > v and w not in to and np.linalg.norm(self.V[w] - self.V[v]) <= d:
                                to[w] = v
        return self.merge_map(to)

    def merge_map(self, to):
        """to: {nokta: yerine gececegi nokta} (hedefler to'da anahtar degil). Temizlik kurallari merge_points'te.
        Birlesen nokta sayisini dondurur."""
        if not to:
            return 0
        before = {v for f in self.F for v in f}
        faces, touched = [], set()
        for f in self.F:
            path, loops = [], []
            for v in (to.get(v, v) for v in f):
                if v in path:                              # kendine degdi: aradaki parca ayri cokgen
                    i = path.index(v)
                    loops.append(path[i:])
                    del path[i + 1:]
                else:
                    path.append(v)
            loops.append(path)
            hit = any(v in to for v in f)
            for L in loops:
                if len(L) >= 3:
                    if hit:
                        touched.add(len(faces))
                    faces.append(L)

        def canon(f):
            i = f.index(min(f))
            return tuple(f[i:] + f[:i])
        groups = {}
        for fi, f in enumerate(faces):
            groups.setdefault(frozenset(f), []).append(fi)
        dead = set()
        for g in groups.values():
            if len(g) < 2 or touched.isdisjoint(g):
                continue
            kept = {}
            for fi in g:
                c, r = canon(faces[fi]), canon(faces[fi][::-1])
                if r in kept:
                    dead.update((fi, kept.pop(r)))
                elif c in kept:
                    dead.add(fi)
                else:
                    kept[c] = fi
        self.F = [f for fi, f in enumerate(faces) if fi not in dead]
        C = {}
        for e, s in self.C.items():
            a, b = (to.get(v, v) for v in e)
            if a != b:
                C[_edge(a, b)] = max(s, C.get(_edge(a, b), s))
        self.C = C
        after = {v for f in self.F for v in f}             # birlesmeyle bosa dusenler gider, onceden bostakiler kalir
        self._compact([v for v in range(len(self.V)) if v not in to and (v in after or v not in before)])
        return len(to)

    def face_regions(self, sel):
        """sel yuzlerini kenar komsulugu ile parcalara ayirir (yalniz noktada degmek komsuluk sayilmaz)."""
        ef = {}
        for fi in sel:
            f = self.F[fi]
            for i in range(len(f)):
                ef.setdefault(_edge(f[i], f[(i + 1) % len(f)]), []).append(fi)
        left, out = set(sel), []
        while left:
            stack, region = [left.pop()], set()
            while stack:
                fi = stack.pop()
                region.add(fi)
                f = self.F[fi]
                for i in range(len(f)):
                    for g in ef[_edge(f[i], f[(i + 1) % len(f)])]:
                        if g in left:
                            left.remove(g)
                            stack.append(g)
            out.append(region)
        return out

    def face_outline(self, region):
        """Yuz bolgesinin dis siniri, yuzlerin yonunde tek kapali nokta listesi; delikli / kendine degen bolge, tamamen
        kapali kabuk ya da karisik yonlu yuzlerde None."""
        count, step = {}, []
        for fi in region:
            f = self.F[fi]
            for i in range(len(f)):
                a, b = f[i], f[(i + 1) % len(f)]
                count[_edge(a, b)] = count.get(_edge(a, b), 0) + 1
                step.append((a, b))
        nxt = {}
        for a, b in step:
            if count[_edge(a, b)] == 1:
                if a in nxt:
                    return None
                nxt[a] = b
        if len(nxt) < 3:
            return None
        start = min(nxt)
        loop, v = [start], nxt[start]
        while v != start:
            if v not in nxt or len(loop) >= len(nxt):
                return None
            loop.append(v)
            v = nxt[v]
        return loop if len(loop) == len(nxt) else None

    def merge_faces(self, regions, loops):
        """LightWave Merge Polygons: her bolge sinir halkasindan tek cokgen olur; ic kenarlar ve yalniz onlara bagli
        ic noktalar gider, sinirdaki keskin kenarlar kalir. Yeni cokgenlerin indekslerini dondurur."""
        before = {v for f in self.F for v in f}
        repl = {min(r): L for r, L in zip(regions, loops)}
        dead = set().union(*regions) - set(repl)
        F, new = [], []
        for fi, f in enumerate(self.F):
            if fi in repl:
                new.append(len(F))
                F.append(list(repl[fi]))
            elif fi not in dead:
                F.append(f)
        self.F = F
        after = {v for f in self.F for v in f}
        self._compact([v for v in range(len(self.V)) if v in after or v not in before])
        return new

    def delete(self, sel):
        before = {v for f in self.F for v in f}
        self.F = [f for i, f in enumerate(self.F) if i not in sel]
        after = {v for f in self.F for v in f}             # bosa dusenler gider, onceden serbest noktalar kalir
        self._compact([v for v in range(len(self.V)) if v in after or v not in before])

    def _compact(self, used):
        """Yalniz used noktalar (ve serbest kenarlarin noktalari) kalir; yuzler, keskin kenarlar, serbest kenarlar
        ve kilavuzlar yeni numaralara gecer."""
        if self.W:
            used = sorted(set(used) | {v for e in self.W for v in e})
        remap = {v: i for i, v in enumerate(used)}
        self.V = self.V[used] if used else np.zeros((0, 3))
        self.F = [[remap[v] for v in f] for f in self.F]
        self.W = {_edge(*(remap[v] for v in e)) for e in self.W}
        alive = self.face_edges()
        self.C = {_edge(*(remap[v] for v in e)): s for e, s in self.C.items()
                  if all(v in remap for v in e) and _edge(*(remap[v] for v in e)) in alive}
        self.rails = [r for r in (_rail_remap(r, remap) for r in self.rails) if r is not None]
        self.dims = [d for d in (_dim_remap(d, remap) for d in self.dims) if d is not None]


# ---------- Rail bridge (kilavuzlu kopru) ----------
# Iki uc halka (A, B; tek uclu kullanimda yalniz A) arasinda bir kilavuz egri: A merkezi, ara noktalar (her
# birinde bir halka) ve B merkezi / serbest uc. Egri dugumlerden gecen Hermite parcalaridir; A'dan (ve B'ye)
# halka normali yonunde, yani poligona dik cikar. Her ara halka egriye diktir: cerceve A'dan baslayip egri boyunca
# burulmadan tasinir; halka sekli A'dan B'ye kademeli gecer (tek uclu: A'nin sekli). Kilavuz kafesle saklanir
# (Cage.rails, belgede Rails ozelligi): sonradan cift tiklayip yeniden duzenlenir. Kayit (sozluk, yerel koordinat):
#   a, b: uc halka nokta indeksleri (b None: tek uclu) | rings: ara halkalar (tek uclu: son halka serbest uc)
#   pts, sc: her halkanin merkezi ve olcegi | ta, tb: A'dan cikis, B'ye giris yonu | closed: halkalar kapali mi
#   cap, capq: tek uclu ucun butun noktalari ve A cercevesindeki konumlari (uc kapagi kati tasinir)

def _rail_remap(r, remap):
    """Cage.delete sonrasi nokta indekslerini yeni numaralara cevir (kayit yerinde guncellenir: elinde tutan
    duzenleme ayni kaydi gormeye devam eder); noktasi silinen kilavuz kopar (None)."""
    try:
        m = lambda L: [remap[v] for v in L]
        new = dict(a=m(r["a"]), b=m(r["b"]) if r.get("b") is not None else None, rings=[m(x) for x in r["rings"]],
                   cap=m(r["cap"]) if r.get("cap") is not None else None)
    except KeyError:
        return None
    r.update(new)
    return r


def rail_full(r):
    """Uctan uca halka listesi: A, ara halkalar, B."""
    return [r["a"]] + r["rings"] + ([r["b"]] if r.get("b") is not None else [])


def rail_ok(cage, r):
    """Kilavuzun halkalari kafeste hala ardisik kenarlarla bagli mi (loop cut, silme vb. bozmadiysa)."""
    full = rail_full(r)
    n = len(r["a"])
    if any(len(x) != n for x in full) or any(v >= len(cage.V) for x in full for v in x):
        return False
    if len(r["pts"]) != len(r["rings"]) or len(r["sc"]) != len(r["rings"]):
        return False
    E = cage.edges()
    return all(_edge(X[j], Y[j]) in E for X, Y in zip(full, full[1:]) for j in range(n))


def _axis_rot(d, ang):
    """Birim d ekseni etrafinda ang radyan donme matrisi (Rodrigues)."""
    x, y, z = d
    K = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) + math.sin(ang) * K + (1.0 - math.cos(ang)) * (K @ K)


def _rotate_to(a, b):
    """a birim yonunu b'ye en kisa yoldan ceviren donme matrisi."""
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if c < -0.999999:                                      # ters yon: a'ya dik bir eksen etrafinda 180 derece
        ax = np.cross(a, [1.0, 0, 0])
        if np.linalg.norm(ax) < 1e-6:
            ax = np.cross(a, [0, 1.0, 0])
        ax /= np.linalg.norm(ax)
        return 2 * np.outer(ax, ax) - np.eye(3)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K / (1.0 + c)


def _ring_frame(P, hint):
    """Halka cercevesi: merkez, u (merkezden ilk noktaya), v, n (hint yonune bakan halka normali)."""
    c = P.mean(0)
    n = np.zeros(3)
    for i in range(len(P)):                                # Newell normali
        n += np.cross(P[i] - c, P[(i + 1) % len(P)] - c)
    hint = np.asarray(hint, float)
    if np.linalg.norm(n) < 1e-9:
        n = hint
    n = n / np.linalg.norm(n)
    if n.dot(hint) < 0:
        n = -n
    u = P[0] - c
    u = u - u.dot(n) * n
    if np.linalg.norm(u) < 1e-9:
        u = np.cross(n, [1.0, 0, 0]) if abs(n[0]) < 0.9 else np.cross(n, [0, 1.0, 0])
    u = u / np.linalg.norm(u)
    return c, u, np.cross(n, u), n


def _hermite(P0, P1, m0, m1, t):
    t2, t3 = t * t, t * t * t
    return (2 * t3 - 3 * t2 + 1) * P0 + (t3 - 2 * t2 + t) * m0 + (-2 * t3 + 3 * t2) * P1 + (t3 - t2) * m1


def rail_nodes(cage, r):
    """Kilavuz dugumleri (yerel), dugum tegetleri (Hermite) ve A / B cerceveleri. Tension (kayitta, 0..1):
    tegetler (1 - tension) ile kisalir; 0 yumusak egri, buyudukce noktalar arasi daha duz / gergin."""
    fa = _ring_frame(cage.V[r["a"]], r["ta"])
    N = [fa[0]] + [np.asarray(p, float) for p in r["pts"]]
    fb = None
    if r.get("b") is not None:
        fb = _ring_frame(cage.V[r["b"]], r["tb"])
        N.append(fb[0])
    N = np.array(N)
    M = np.zeros_like(N)
    for i in range(1, len(N) - 1):
        M[i] = (N[i + 1] - N[i - 1]) / 2
    M[0] = fa[3] * np.linalg.norm(N[1] - N[0])             # A'ya dik cikis
    M[-1] = fb[3] * np.linalg.norm(N[-1] - N[-2]) if fb is not None else N[-1] - N[-2]
    M *= max(1.0 - min(max(float(r.get("tension", 0.0)), 0.0), 1.0), 0.02)   # yon ayni kalir (halka cerceveleri)
    return N, M, fa, fb


def rail_seg(N, M, i):
    """i. parcanin Hermite tegetleri. Ara dugumun tegeti iki komsusundan (Catmull-Rom); aralar esit degilse
    parcanin kendi boyuna gore olceklenir: yoksa uzun araliktan sonraki kisa parcada egri tasip gobek yapardi."""
    l = np.linalg.norm(np.diff(N, axis=0), axis=1)
    m0, m1 = M[i], M[i + 1]
    if 0 < i:
        m0 = m0 * (2 * l[i] / max(l[i - 1] + l[i], 1e-12))
    if i + 1 < len(N) - 1:
        m1 = m1 * (2 * l[i] / max(l[i] + l[i + 1], 1e-12))
    return m0, m1


def rail_curve(N, M, per=16):
    """Egri cizgisi (yerel noktalar), parca basina per adim."""
    pts = [N[0]]
    for i in range(len(N) - 1):
        m0, m1 = rail_seg(N, M, i)
        for k in range(1, per + 1):
            pts.append(_hermite(N[i], N[i + 1], m0, m1, k / per))
    return np.array(pts)


def rail_place(cage, r):
    """Halkalari (ve tek uclu ucun kapagini) kilavuza gore yerlestir: her halka kendi dugumunde egriye dik."""
    N, M, fa, fb = rail_nodes(cage, r)
    cA, u, v, t = fa
    qa = (cage.V[r["a"]] - cA) @ np.array([u, v, t]).T     # A'nin sekli kendi cercevesinde
    frames = [(u, v, t)]
    for i in range(1, len(N)):                             # burulmasiz tasima: tegetten tegete en kisa donme
        ti = M[i] / max(np.linalg.norm(M[i]), 1e-12)
        ui = _rotate_to(frames[-1][2], ti) @ frames[-1][0]
        ui = ui - ui.dot(ti) * ti
        ui /= max(np.linalg.norm(ui), 1e-12)
        frames.append((ui, np.cross(ti, ui), ti))
    L = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(N, axis=0), axis=1))])
    s = L / L[-1] if L[-1] > 1e-12 else np.linspace(0, 1, len(N))
    qb = None
    if fb is not None:                                     # B'nin sekli, oraya tasinan cercevede
        qb = (cage.V[r["b"]] - fb[0]) @ np.array(frames[-1]).T
    for i, ring in enumerate(r["rings"], 1):
        q = qa if qb is None else (1 - s[i]) * qa + s[i] * qb
        cage.V[ring] = N[i] + r["sc"][i - 1] * (q @ np.array(frames[i]))
    if r.get("b") is None and r.get("cap"):                # tek uclu: uc kapagi kati tasinir
        cage.V[r["cap"]] = N[-1] + r["sc"][-1] * (np.array(r["capq"]) @ np.array(frames[-1]))


def rail_cable(N, h, d, fixed_b, elastic, ta, tb=None, iters=30, friction=12.0, bend=6.0):
    """Kilavuzu kablo gibi cek: N dugumler (0: A, hep sabit; son: B, fixed_b ise sabit), h cekilen dugum, d
    tasima; ta / tb: A'dan cikis, B'ye giris yonu (uclarda hayali dugum: kablo dik cikar). Serbest dugumler
    bukulme enerjisi + boy degisimi cezasiyla en aza indirilir (Levenberg-Marquardt): aradaki kivrimlar acilir /
    gevser, gerilince dumduz olup esit uzar. friction: tek uclu rail'de cekilen noktanin ardindaki serbest kuyruk
    masadaki kablo gibi eski yerinde kalmaya direnir: ortadan cekilince orada bombe olur, serbest uc az oynar, hat
    gerilince surukleyerek gelir (ucu itip cekerken ve iki uclu rail'de kuyruk yok). elastic: boy sertligi, 0 zincir (boylar korunur) .. 1 lastik (kolay uzar / kisalir). Yeni dugumler."""
    N = np.asarray(N, float)
    n = len(N)
    L = np.linalg.norm(np.diff(N, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(L)])
    ta = np.asarray(ta, float)
    t = N[h] + d
    # baslangic: tasima yay uzunluguna gore sabit uca dogru azalarak; boyundan kisa kalan bolum yana bombe
    # yapar (dumduz sikisan kablo dengede kalir, bukulmezdi): yon surukleme ya da egrinin eski seklinden
    w = np.clip(s / max(s[h], 1e-12), 0.0, 1.0)
    w[h:] = (s[-1] - s[h:]) / max(s[-1] - s[h], 1e-12)    # tek uclu kuyruk da serbest uca dogru azalir
    w[h] = 1.0
    X = N + np.outer(w, d)
    for a, b in ((0, h), (h, n - 1)) if fixed_b else ((0, h),):
        c = X[b] - X[a]
        D, Ls = np.linalg.norm(c), s[b] - s[a]
        if b - a < 2 or D >= Ls or D < 1e-12:
            continue
        k = c / D
        for v in [d] + list(N[a:b + 1] - N[a]) + list(np.eye(3)):
            v = v - v.dot(k) * k
            if np.linalg.norm(v) > 1e-6 * max(Ls, 1.0):
                X[a:b + 1] += np.outer(np.sin(np.pi * (s[a:b + 1] - s[a]) / Ls),
                                       v / np.linalg.norm(v) * 0.5 * math.sqrt(Ls * Ls - D * D))
                break
    X[0], X[h] = N[0], t
    if fixed_b:
        X[-1] = N[-1]
    Y = np.vstack([[N[0] - ta * L[0]], X] + ([[N[-1] + np.asarray(tb, float) * L[-1]]] if fixed_b else []))
    m = len(Y)
    free = np.ones(m, bool)
    free[[0, 1, h + 1]] = False
    if fixed_b:
        free[[m - 2, m - 1]] = False
    cols = np.repeat(free, 3)
    Y0 = Y.copy()
    Y0[1:n + 1] = N                                        # surtunme: kuyruk dugumlerinin eski yeri
    fw = np.zeros(m)
    if not fixed_b:
        fw[h + 2:n + 1] = friction
    fw = np.repeat(fw, 3)[cols]
    Jf = np.diag(fw)
    sk = 10.0 ** (1.5 - 3.0 * min(max(float(elastic), 0.0), 1.0))   # boy cezasi: 30 (zincir) .. 0.03 (lastik)
    nl = len(L)
    Jb = np.zeros((3 * (m - 2), 3 * m))                    # bukulme satirlari sabit (dogrusal)
    for j in range(1, m - 1):
        for q, c in ((j - 1, 1.0), (j, -2.0), (j + 1, 1.0)):
            Jb[3 * (j - 1):3 * j, 3 * q:3 * q + 3] += bend * c * np.eye(3)
    Jb = Jb[:, cols]
    ri = np.repeat(np.arange(nl), 3)
    ca = 3 * (np.arange(nl) + 1)[:, None] + np.arange(3)   # parca i: Y[i + 1] -> Y[i + 2]

    def resid(Y):
        v = Y[2:] - Y[1:-1] if not fixed_b else Y[2:-1] - Y[1:-2]   # gercek parcalar (hayali uclar haric)
        l = np.linalg.norm(v, axis=1)
        return np.concatenate([sk * (l - L), bend * (Y[:-2] - 2 * Y[1:-1] + Y[2:]).ravel(),
                               fw * (Y[free] - Y0[free]).ravel()]), v, l

    def jac(v, l):
        Jl = np.zeros((nl, 3 * m))
        u = (sk * v / np.maximum(l, 1e-12)[:, None]).ravel()
        Jl[ri, ca.ravel()] = -u
        Jl[ri, (ca + 3).ravel()] = u
        return np.vstack([Jl[:, cols], Jb, Jf])
    lam = 1e-2
    r, v, l = resid(Y)
    cost = r.dot(r)
    for _ in range(iters):
        J = jac(v, l)
        A = J.T @ J
        g = J.T @ r
        while True:
            step = np.linalg.solve(A + lam * np.diag(np.diag(A) + 1e-9), -g)
            Y2 = Y.copy()
            Y2[free] += step.reshape(-1, 3)
            r2, v2, l2 = resid(Y2)
            c2 = r2.dot(r2)
            if c2 < cost or lam > 1e8:
                break
            lam *= 4
        if c2 >= cost:
            break
        Y, r, v, l, cost, lam = Y2, r2, v2, l2, c2, max(lam / 3, 1e-6)
        if np.abs(step).max() < 1e-6:
            break
    return Y[1:n + 1]


def _rail_cols(cage, X, Y):
    """X ve Y halkalari arasindaki dortgenler (her birinde iki X, iki Y noktasi)."""
    xs, ys = set(X), set(Y)
    return [fi for fi, f in enumerate(cage.F)
            if len(f) == 4 and sum(v in xs for v in f) == 2 and sum(v in ys for v in f) == 2]


def rail_insert(cage, r, i):
    """Uctan uca listede i ile i + 1 halkalari arasina yeni halka: aradaki dortgenler ikiye bolunur, yeni
    halkanin dugumu o parcanin ortasinda. Basarisizsa False."""
    full = rail_full(r)
    if not 0 <= i < len(full) - 1:
        return False
    X, Y = full[i], full[i + 1]
    cols = _rail_cols(cage, X, Y)
    n = len(X)
    if len(cols) != (n if r["closed"] else n - 1):
        return False
    N, M, _, _ = rail_nodes(cage, r)
    mid = _hermite(N[i], N[i + 1], *rail_seg(N, M, i), 0.5)
    base = len(cage.V)
    cage.V = np.vstack([cage.V, (cage.V[X] + cage.V[Y]) / 2])
    mx = {X[j]: base + j for j in range(n)}
    my = {Y[j]: base + j for j in range(n)}
    for fi in cols:
        f = cage.F[fi]
        cage.F[fi] = [my.get(v, v) for v in f]             # X tarafi: Y yerine yeni halka
        cage.F.append([mx.get(v, v) for v in f])           # Y tarafi: X yerine yeni halka
    for j in range(n):                                     # boyuna keskin kenarlar ikiye bolunur
        e = _edge(X[j], Y[j])
        if e in cage.C:
            sh = cage.C.pop(e)
            cage.C[_edge(X[j], base + j)] = sh
            cage.C[_edge(base + j, Y[j])] = sh
    sc = (r["sc"][i - 1] if i > 0 else 1.0) + (r["sc"][i] if i < len(r["sc"]) else 1.0)
    r["rings"].insert(i, list(range(base, base + n)))
    r["pts"].insert(i, [float(x) for x in mid])
    r["sc"].insert(i, sc / 2)
    return True


def rail_remove(cage, r, i):
    """Uctan uca listede i. ara halkayi kaldir: iki yanindaki dortgenler birlesir, noktalari silinir.
    Tek uclu kullanimda serbest uc kaldirilamaz. Basarisizsa False."""
    full = rail_full(r)
    if not 0 < i < len(full) - 1:
        return False
    X, Y, Z = full[i - 1], full[i], full[i + 1]
    before, after = _rail_cols(cage, X, Y), _rail_cols(cage, Y, Z)
    if len(before) != len(after) or not before:
        return False
    yx = dict(zip(Y, X))
    for fi in after:
        cage.F[fi] = [yx.get(v, v) for v in cage.F[fi]]
    for j in range(len(X)):
        sa, sb = cage.C.pop(_edge(X[j], Y[j]), None), cage.C.pop(_edge(Y[j], Z[j]), None)
        if sa is not None and sb is not None:
            cage.C[_edge(X[j], Z[j])] = max(sa, sb)
    drop = set(before)
    cage.F = [f for fi, f in enumerate(cage.F) if fi not in drop]
    del r["rings"][i - 1], r["pts"][i - 1], r["sc"][i - 1]
    cage.delete(set())                                     # bosta kalan halka noktalari (indeksler kayar)
    return True


def rail_resample(cage, r, k):
    """Ara halka sayisini k yap (tek uclu: serbest uc haric): en uzun parcaya ekleyerek / en kisa araliktaki
    halkayi kaldirarak; sonra ara dugumler eski egri boyunca esit araliga dizilir (egrinin sekli korunur)."""
    one = r.get("b") is None
    N, M, _, _ = rail_nodes(cage, r)
    old = rail_curve(N, M, 32)
    sc_old = [1.0] + list(r["sc"]) + ([] if one else [1.0])
    count = lambda: len(r["rings"]) - (1 if one else 0)
    while count() < k:
        N, _, _, _ = rail_nodes(cage, r)
        span = int(np.argmax(np.linalg.norm(np.diff(N, axis=0), axis=1)))
        if not rail_insert(cage, r, span):
            return False
    while count() > k:
        N, _, _, _ = rail_nodes(cage, r)
        gaps = np.linalg.norm(np.diff(N, axis=0), axis=1)
        i = min(range(1, len(N) - 1), key=lambda j: gaps[j - 1] + gaps[j])   # A, B / serbest uc kalir
        if not rail_remove(cage, r, i):
            return False
    seg = np.linalg.norm(np.diff(old, axis=0), axis=1)
    acc = np.concatenate([[0.0], np.cumsum(seg)])
    total = acc[-1]
    n_in = count()
    stops = n_in + 1                                       # A ... (k ara) ... B / serbest uc
    for j in range(1, n_in + 1):
        d = total * j / stops
        m = min(int(np.searchsorted(acc, d)), len(acc) - 1)
        m0 = max(m - 1, 0)
        w = (d - acc[m0]) / max(acc[m] - acc[m0], 1e-12)
        r["pts"][j - 1] = [float(x) for x in old[m0] + w * (old[m] - old[m0])]
        f = j / stops * (len(sc_old) - 1)                  # olcek eski dugum olcekleri arasinda
        a = min(int(f), len(sc_old) - 2)
        r["sc"][j - 1] = float(sc_old[a] + (f - a) * (sc_old[a + 1] - sc_old[a]))
    return True


def rail_bridge2(cage, A, B, closed, dA, dB, flip, region, k):
    """Iki uclu Rail bridge: A ile B arasina k ara halkali serit (Bridge gibi; bolge yuzleri silinir) ve
    kilavuz. Baslangic egrisi A'ya ve B'ye dik giren yumusak egri. Kilavuz kaydini dondurur."""
    ta, tb = dA.sum(0), -dB.sum(0)                         # A'dan cikis, B'ye giris yonu
    ta = ta / max(np.linalg.norm(ta), 1e-12)
    tb = tb / max(np.linalg.norm(tb), 1e-12)
    n, base = len(A), len(cage.V)
    cage.bridge(A, B, closed, k, 0.0, 0.0, None, None, flip)
    r = dict(a=list(A), b=list(B), rings=[list(range(base + j * n, base + (j + 1) * n)) for j in range(k)],
             pts=[], sc=[1.0] * k, ta=[float(x) for x in ta], tb=[float(x) for x in tb], closed=bool(closed),
             cap=None, capq=None)
    cA, cB = cage.V[A].mean(0), cage.V[B].mean(0)
    L = np.linalg.norm(cB - cA)
    curve = np.array([_hermite(cA, cB, ta * L, tb * L, t) for t in np.linspace(0, 1, 129)])
    acc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(curve, axis=0), axis=1))])
    for j in range(1, k + 1):                              # ara dugumler egri boyunca esit aralikli
        m = min(int(np.searchsorted(acc, acc[-1] * j / (k + 1))), len(acc) - 1)
        r["pts"].append([float(x) for x in curve[m]])
    cage.rails.append(r)
    if region:
        cage.delete(region)                                # kilavuz indeksleri delete ile yeni numaralara gecer
    r = cage.rails[-1]
    rail_place(cage, r)
    return r


def rail_bridge1(cage, region, A, k, length=None):
    """Tek uclu Rail bridge: secili bolge (siniri A halkasi) k ara halka ve serbest uc halkasiyla uzar; bolgenin
    yuzleri ucta kapak olur (extrude gibi). Baslangicta uc, A'nin normali yonunde length kadar ileride."""
    P = cage.V[A]
    n0 = sum(cage.normal(fi) for fi in region)
    cA, u, v, t = _ring_frame(P, n0 if np.linalg.norm(n0) > 1e-12 else [0, 0, 1.0])
    if length is None:
        length = 4.0 * float(np.linalg.norm(P - cA, axis=1).mean())
    layers = {}
    verts = sorted({x for fi in region for x in cage.F[fi]})
    capq = (cage.V[verts] - cA) @ np.array([u, v, t]).T
    cage.extrude(set(region), k + 1, layers)
    r = dict(a=list(A), b=None, rings=[[layers[x][j] for x in A] for j in range(1, k + 2)],
             pts=[[float(c) for c in cA + t * length * j / (k + 1)] for j in range(1, k + 2)],
             sc=[1.0] * (k + 1), ta=[float(c) for c in t], tb=None, closed=True,
             cap=[layers[x][-1] for x in verts], capq=capq.tolist())
    cage.rails.append(r)
    cage.delete(set())                                     # ic noktalarin ara katlari bosta kalir
    r = cage.rails[-1]
    rail_place(cage, r)
    return r

# ---------- Belgedeki SubD nesnesi ----------
# Kafes belgede Mesh::FeaturePython nesnesinde saklanir: CagePoints (noktalar) + CageFaces (indeksler, -1
# ayirici). execute() kafesten SubPatch yuzeyini uretir; nesne secilip makro calistirilinca yeniden duzenlenir.

def subd_mesh(cage, levels, rules=RULES):
    S, Q = subdivide(cage.V, cage.F, levels, cage.C, rules)
    T = np.concatenate([Q[:, [0, 1, 2]], Q[:, [0, 2, 3]]])
    return Mesh.Mesh(S[T].tolist())


# ---------- SubPatch yuzeyinden B-rep kati (STEP) ----------
# Bir Catmull-Clark adimindan sonra kafes tamamen dortgendir; her dortgen bir "yama" olur. Yama k adim daha
# bolunur, cocuk dortgenlerin yama icindeki (u, v) konumlari izlenerek (2^k + 1)^2 noktalik bir izgara kurulur
# ve bu izgaraya bikubik B-spline yuzey oturtulur. Komsu yamalar ayni sinir noktalarini paylastigi icin
# kenarlari cakisir; yuzler dikilip (sew) kapali kabuktan kati yapilir.

def patch_grids(cage, k=3, rules=RULES):
    V, F, C = catmull_clark(cage.V, cage.F, cage.C, rules)
    M = len(F)
    P = np.tile(np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float), (M, 1, 1))   # her dortgenin kose (u, v)'si
    for _ in range(k):
        V, F, C = catmull_clark(V, F, C, rules)
        P = _child_uv(P)
    n = 2 ** k
    grids = np.zeros((M, n + 1, n + 1, 3))
    parent = np.arange(len(F)) // (4 ** k)                  # cocuklar ebeveyn sirasinda bitisik bloklar
    ij = np.rint(P * n).astype(int)
    grids[np.repeat(parent, 4), ij[..., 0].ravel(), ij[..., 1].ravel()] = V[F].reshape(-1, 3)
    return grids


def _child_uv(P):
    """Dortgenlerin kose (u, v)'lerinden cocuklarinkiler. j. kosenin cocugu: [kose j, (j, j+1) kenar ortasi,
    merkez, (j-1, j) kenar ortasi] (catmull_clark cocuk sirasi)."""
    child = np.stack([P, (P + np.roll(P, -1, 1)) / 2, np.repeat(P.mean(1, keepdims=True), 4, 1),
                      (P + np.roll(P, 1, 1)) / 2], 2)
    return child.reshape(-1, 4, 2)


# ---------- G1 / G2: duzgun hucreler birebir ----------
# Catmull-Clark yuzeyi duzgun hucrelerde (4 kosesi de 4 yuzlu ic nokta, koselere degen keskin kenar yok) duz
# bikubik B-spline'dir: hucrenin 4 x 4 komsu kafes noktasi kutuplaridir. Boyle hucreler yaklasiksiz Bezier yamasina
# cevrilir; komsu duzgun yamalar egrilikte de (C2) sureklidir. Duzgun olmayan hucreler (3 / 5+ yuzlu nokta, acik ya
# da keskin kenar komsulugu) G2'de birkac seviye daha bolunur: cocuklarinin cogu yine duzgundur. Geriye kalan
# (giderek kuculen) hucrelere limit yuzey noktalarindan B-spline oturtulur. Gercek SubD yuzeyi 3 / 5+ yuzlu
# noktalarda zaten yalniz G1'dir; tam G2 orada matematiksel olarak yoktur.

_BSP2BEZ = np.array([[1, 4, 1, 0], [0, 4, 2, 0], [0, 2, 4, 0], [0, 1, 4, 1]], float) / 6.0   # kubik B-spline -> Bezier
CONTINUITY = ("G0", "G1", "G2")
PRECISION = {"Draft": (2, 1, 1e-2), "Normal": (3, 2, 1e-3), "Fine": (4, 3, 1e-4)}   # Samples, Refine, Tolerance
_NET = (((0, 0), (1, 0), (0, -1), (1, -1), (-1, -1), (-1, 0)),   # hucre kenari p -> q: karsi yuzde p'nin ve q'nun
        ((1, 0), (1, 1), (2, 0), (2, 1), (2, -1), (1, -1)),      # dis komsusu, p'nin kosegeni, p'nin oteki dis
        ((1, 1), (0, 1), (1, 2), (0, 2), (2, 2), (2, 1)),        # komsusu (izgara konumlari; a = (0, 0),
        ((0, 1), (0, 0), (-1, 1), (-1, 0), (-1, 2), (0, 2)))     # b = (1, 0), c = (1, 1), d = (0, 1))


def _rot(f, a):
    i = f.index(a)
    return f[i:] + f[:i]


def _regular_net(F, he, valence, sharp, fi):
    """fi dortgen hucresinin 4 x 4 kontrol kafesi (N[i][j]: i a -> b, j a -> d yonunde) ya da duzgun degilse None.
    F: yuz listeleri, he: yonlu kenar -> yuz, valence: noktadaki yuz sayisi, sharp: keskin kenarli noktalar."""
    f = F[fi]
    if any(valence[v] != 4 or v in sharp for v in f):
        return None
    G = {(0, 0): f[0], (1, 0): f[1], (1, 1): f[2], (0, 1): f[3]}
    for p, q, xp, yp, _, _ in _NET:
        g = he.get((G[q], G[p]))                           # kenarin karsi yuzu: [q, p, x, y]
        if g is None:
            return None
        g = _rot(F[g], G[q])
        G[xp], G[yp] = g[2], g[3]
    for p, _, xp, _, dp, wp in _NET:
        h = he.get((G[xp], G[p]))                          # p'nin kosegen yuzu: [x, p, w, kosegen]
        if h is None:
            return None
        h = _rot(F[h], G[xp])
        if h[2] != G[wp]:
            return None
        G[dp] = h[3]
    return [[G[(i - 1, j - 1)] for j in range(4)] for i in range(4)]


def _crop(V, F, C, keep):
    """keep yuzlerinin 2 halkalik komsulugu: (V, F, C, keep'in yeni indeksleri). Bir Catmull-Clark adimi keep
    yuzlerinin cocuklarini ve onlarin 1 halka komsularini burada da dogru hesaplar (kesik sinir yalniz dis halkayi
    bozar); boylece yalniz gereken bolge bolunur."""
    fs = np.zeros(len(F), bool)
    fs[keep] = True
    for _ in range(2):
        vs = np.zeros(len(V), bool)
        vs[F[fs].ravel()] = True
        fs = vs[F].any(1)
    faces = np.nonzero(fs)[0]
    used = np.unique(F[faces])
    remap = -np.ones(len(V), int)
    remap[used] = np.arange(len(used))
    fmap = -np.ones(len(F), int)
    fmap[faces] = np.arange(len(faces))
    C2 = {}
    for e, s in C.items():
        a, b = tuple(e)
        if remap[a] >= 0 and remap[b] >= 0:
            C2[_edge(int(remap[a]), int(remap[b]))] = s
    return V[used], remap[F[faces]], C2, fmap[keep]


def _limit(V, F, C, rules=RULES):
    """Dortgen mesh noktalarinin Catmull-Clark limit yuzeyindeki yerleri. Ic nokta (n yuzlu):
    (n^2 v + 4 kenar komsulari + kosegenler) / (n (n + 5)); acik / keskin kenar egrisi uzerinde (onceki + 4 v +
    sonraki) / 6; 3+ keskin kenar ya da sivri acik kenar kosesi yerinde kalir."""
    N = len(V)
    fv = F.ravel()
    nxt, prv, opp = (np.roll(F, -s, 1).ravel() for s in (1, -1, 2))
    n = np.bincount(fv, minlength=N).astype(float)
    es, fs = np.zeros((N, 3)), np.zeros((N, 3))
    np.add.at(es, fv, V[nxt] + V[prv])                     # her kenar komsusu iki yuzde: 2 kat
    np.add.at(fs, fv, V[opp])
    nn = np.maximum(n, 1)[:, None]
    L = (nn * nn * V + 2 * es + fs) / (nn * (nn + 5))
    keys = np.minimum(fv, nxt) * N + np.maximum(fv, nxt)
    uniq, cnt = np.unique(keys, return_counts=True)
    curve = {(int(k // N), int(k % N)) for k in uniq[cnt == 1]}
    curve |= {tuple(sorted(e)) for e, s in C.items() if s > 0}
    cc, cs = np.zeros(N), np.zeros((N, 3))
    for a, b in curve:
        cc[a] += 1
        cc[b] += 1
        cs[a] += V[b]
        cs[b] += V[a]
    two = cc == 2
    L[two] = (cs[two] + 4 * V[two]) / 6
    pin = cc >= 3
    if rules[0] == BOUNDARY_RULES[0]:
        pin |= two & (n == 1)
    L[pin] = V[pin]
    return L


def _fit_grids(V, F, C, cells, k, rules):
    """cells hucrelerinin (2^k + 1)^2 limit yuzeyi noktasi izgarasi (patch_grids gibi, yalniz bu hucreler)."""
    M = len(cells)
    P = np.tile(np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float), (M, 1, 1))
    act = np.asarray(cells, int)
    for _ in range(k):
        V, F, C, act = _crop(V, F, C, act)
        V, F, C = catmull_clark(V, F, C, rules)
        act = (4 * act[:, None] + np.arange(4)).ravel()   # dortgenin j. kose cocugu 4 f + j
        P = _child_uv(P)
    n = 2 ** k
    grids = np.zeros((M, n + 1, n + 1, 3))
    ij = np.rint(P * n).astype(int)
    L = _limit(V, F, C, rules)
    grids[np.repeat(np.arange(len(act)) // 4 ** k, 4), ij[..., 0].ravel(), ij[..., 1].ravel()] = \
        L[F[act]].reshape(-1, 3)
    return grids


def solid_patches(cage, k=3, rules=RULES, continuity="G0", refine=2):
    """B-rep yamalari: (birebir Bezier kutuplari [m, 4, 4, 3], B-spline oturtulacak izgaralar [m, n, n, 3]).
    G0: her yama alt bolme noktalarindan oturtulur (eski yontem). G1: duzgun yamalar birebir, otekiler limit
    noktalarindan. G2: duzgun olmayanlar refine seviye daha bolunur, cocuklarin duzgun olanlari da birebir."""
    if continuity not in CONTINUITY[1:]:
        return np.zeros((0, 4, 4, 3)), patch_grids(cage, k, rules)
    V, F, C = catmull_clark(cage.V, cage.F, cage.C, rules)
    act = np.arange(len(F))
    nets = []
    depth = max(int(refine), 0) if continuity == "G2" else 0
    for level in range(depth + 1):
        Fl = F.tolist()
        he = {(f[i], f[(i + 1) % 4]): fi for fi, f in enumerate(Fl) for i in range(4)}
        valence = np.bincount(F.ravel(), minlength=len(V))
        sharp = {v for e, s in C.items() if s > 0 for v in e}
        rest = []
        for fi in act.tolist():
            net = _regular_net(Fl, he, valence, sharp, fi)
            if net is None:
                rest.append(fi)
            else:
                nets.append(V[np.array(net)])
        act = np.array(rest, int)
        if level == depth or not len(act):
            break
        V, F, C, act = _crop(V, F, C, act)
        V, F, C = catmull_clark(V, F, C, rules)
        act = (4 * act[:, None] + np.arange(4)).ravel()
    n = 2 ** k + 1
    grids = _fit_grids(V, F, C, act, k, rules) if len(act) else np.zeros((0, n, n, 3))
    bez = np.einsum("ia,mabc,jb->mijc", _BSP2BEZ, np.array(nets).reshape(-1, 4, 4, 3), _BSP2BEZ)
    return bez, grids


def _cubic_interp(m):
    """m noktalik esit aralikli kubik interpolasyon: (A^-1, dugumler, katlar). Dugumler ortalama yontemiyle;
    uclar sabit (clamped) oldugundan ilk / son kutup ilk / son noktadir ve A ters cevrilince sinir egrisi
    yalniz sinir noktalarina bagli kalir. Esit aralik simetrik: ters yonde ayni egriyi verir."""
    t = np.linspace(0.0, 1.0, m)
    inner = [t[j:j + 3].mean() for j in range(1, m - 3)]
    U = [0.0] * 4 + inner + [1.0] * 4

    def N(i, p, x):
        if p == 0:
            return 1.0 if (U[i] <= x < U[i + 1]) or (x == 1.0 and U[i] < U[i + 1] == 1.0) else 0.0
        a = (x - U[i]) / (U[i + p] - U[i]) * N(i, p - 1, x) if U[i + p] > U[i] else 0.0
        b = (U[i + p + 1] - x) / (U[i + p + 1] - U[i + 1]) * N(i + 1, p - 1, x) if U[i + p + 1] > U[i + 1] else 0.0
        return a + b
    A = np.array([[N(j, 3, x) for j in range(m)] for x in t])
    return np.linalg.inv(A), [0.0] + inner + [1.0], [4] + [1] * len(inner) + [4]


def subd_solid(cage, k=3, tol=1e-3, rules=RULES, continuity="G0", refine=2):
    """SubPatch yuzeyinin B-spline yamalarindan kati (kafes kapali degilse kabuk) dondurur."""
    import Part
    bez, grids = solid_patches(cage, k, rules, continuity, refine)
    faces = []

    def face(P, mults, knots):
        s = Part.BSplineSurface()
        s.buildFromPolesMultsKnots([[App.Vector(*p) for p in row] for row in P.tolist()],
                                   mults, mults, knots, knots, False, False, 3, 3)
        faces.append(s.toShape())
    for P in bez:                                          # birebir yamalar: tek aralikli kubik Bezier
        face(P, [4, 4], [0.0, 1.0])
    if len(grids):
        Ainv, knots, mults = _cubic_interp(grids.shape[1])
        for P in np.einsum("ia,mabc,jb->mijc", Ainv, grids, Ainv):   # her yama icin A^-1 G A^-T
            face(P, mults, knots)
    shell = Part.Shell(faces)
    shell.sewShape(max(tol * 10, 1e-4))
    shells = shell.Shells
    if len(shells) == 1 and shells[0].isClosed():
        solid = Part.Solid(shells[0])
        if solid.Volume < 0:
            solid.reverse()
        return solid
    return shell


def _solid_props(obj, continuity=None):
    """Kati nesnesinin Continuity / Refine ozellikleri (eski katilarda yok: onlar G0 olarak kalir)."""
    if "Continuity" not in obj.PropertiesList:
        obj.addProperty("App::PropertyEnumeration", "Continuity", "SubD",
                        "G0 = fast fit, G1 = regular patches exact, G2 = refined exact patches (curvature continuous)")
        obj.Continuity = list(CONTINUITY)
        obj.Continuity = continuity or "G0"
    if "Refine" not in obj.PropertiesList:
        obj.addProperty("App::PropertyInteger", "Refine", "SubD",
                        "G2: extra subdivision levels around irregular points (more exact patches)").Refine = 2


def solid_settings(obj, precision=None, continuity=None):
    """Hassasiyet on ayarini (Draft / Normal / Fine) ve surekliligi kati nesnesine yaz."""
    _solid_props(obj, continuity)
    if precision in PRECISION:
        obj.Samples, obj.Refine, obj.Tolerance = PRECISION[precision]
    if continuity in CONTINUITY and obj.Continuity != continuity:
        obj.Continuity = continuity


def solid_precision(obj):
    """Katinin degerlerine uyan on ayar adi ya da "Custom"."""
    vals = (obj.Samples, getattr(obj, "Refine", 2), obj.Tolerance)
    return next((n for n, p in PRECISION.items() if p[0] == vals[0] and p[1] == vals[1]
                 and abs(p[2] - vals[2]) < 1e-12), "Custom")


class SubDSolid:
    """Belgedeki SubD nesnesinden uretilen parametrik B-rep kati: kafes degisince yeniden hesaplanir,
    STEP / IGES olarak disa aktarilabilir, Part / PartDesign araclariyla kullanilabilir."""

    def __init__(self, obj, source):
        obj.addProperty("App::PropertyLink", "Source", "SubD", "Source SubD object").Source = source
        obj.addProperty("App::PropertyInteger", "Samples", "SubD",
                        "Subdivision steps per patch (2^n + 1 points); higher = closer to the SubD surface").Samples = 3
        obj.addProperty("App::PropertyFloat", "Tolerance", "SubD", "B-spline approximation tolerance (mm)").Tolerance = 0.001
        solid_settings(obj, _param().GetString("Precision", "Normal"), _param().GetString("Continuity", "G2"))
        obj.Proxy = self

    def execute(self, obj):
        if is_subd(obj.Source):
            # kaynak yalniz tasindiysa (kafes ayni) yamalari yeniden kurma: hem hizli hem geri alma bozulmaz
            cont, refine = str(getattr(obj, "Continuity", "G0")), int(getattr(obj, "Refine", 2))
            key = (cage_key(obj.Source), obj.Samples, obj.Tolerance, cont, refine)
            if getattr(self, "_key", None) != key or obj.Shape.isNull():
                cage = cage_from_obj(obj.Source)
                if not cage.F:
                    return
                obj.Shape = subd_solid(cage, max(1, min(obj.Samples, 5)), max(obj.Tolerance, 1e-6),
                                       obj_rules(obj.Source), cont, max(0, min(refine, 5)))
                self._key = key
            # kafes nesnenin yerel koordinatinda: kati da kaynakla ayni yerde dursun
            if not obj.Placement.isSame(obj.Source.Placement, 1e-12):
                obj.Placement = obj.Source.Placement

    def dumps(self):
        return None

    def loads(self, state):
        return None


def repair_proxies(doc):
    """Python baglantisi (Proxy) kopmus SubD nesnelerini ve katilarini onarir, yeniden hesaplar. Belge acilirken
    FreeCAD modulu ice aktaramazsa (or. FreeCAD 26.3 'module not permitted': Mod klasoru disindaki eklenti)
    nesneler Proxy'siz kalir; o zaman kafes kaydedilir ama yuzey / kati guncellenmez (eski sekil gorunur)."""
    fixed = []
    for o in (doc.Objects if doc is not None else []):
        props = getattr(o, "PropertiesList", [])
        cls = SubDObject if "CageFaces" in props else SubDSolid if "Samples" in props and "Source" in props else None
        if cls is None or type(getattr(o, "Proxy", None)).__name__ == cls.__name__:
            continue
        o.Proxy = cls.__new__(cls)                         # __init__ ozellik ekler; onlar zaten var
        o.touch()
        fixed.append(o)
    for o in fixed:
        o.recompute()
    if fixed:
        App.Console.PrintMessage(T("SubD: repaired %d objects (%s)\n") % (len(fixed), ", ".join(o.Name for o in fixed)))
    return fixed


def make_solid(subd_obj):
    doc = subd_obj.Document
    o = doc.addObject("Part::FeaturePython", subd_obj.Name + "_Solid")
    SubDSolid(o, subd_obj)
    if App.GuiUp:
        o.ViewObject.Proxy = 0
    doc.recompute()
    return o


def _surface_props(obj):
    """Yuzey kurallari ozellikleri (eski nesnelerde yok; yoksa varsayilan kurallar gecerli)."""
    if "Boundary" not in obj.PropertiesList:
        obj.addProperty("App::PropertyEnumeration", "Boundary", "SubD",
                        "Open edge rule: Edge and corner = one-face corners stay sharp, Edge only = they round too")
        obj.Boundary = list(BOUNDARY_RULES)
        obj.Boundary = RULES[0]
    if "SmoothTriangles" not in obj.PropertiesList:
        obj.addProperty("App::PropertyBool", "SmoothTriangles", "SubD",
                        "Loop-like weights next to triangles (smoother triangle fans)").SmoothTriangles = RULES[1]


def obj_rules(obj):
    """SubD nesnesinin alt bolme kurallari: (acik kenar kurali, ucgen yumusatma)."""
    return (str(getattr(obj, "Boundary", RULES[0])), bool(getattr(obj, "SmoothTriangles", RULES[1])))


def set_surface(obj, levels, rules):
    """Seviye ve kurallari nesneye yaz (yeni / kopya nesneler kaynaklarininkini alir)."""
    _surface_props(obj)
    obj.Levels = levels
    if obj.Boundary != rules[0]:
        obj.Boundary = rules[0]
    if obj.SmoothTriangles != rules[1]:
        obj.SmoothTriangles = rules[1]


def _crease_props(obj):
    if "CageCreases" not in obj.PropertiesList:      # eski nesnelerde yok
        obj.addProperty("App::PropertyIntegerList", "CageCreases", "SubD", "Creased edges: point index pairs")
        obj.addProperty("App::PropertyFloatList", "CageSharpness", "SubD", "Sharpness of the creased edges")


def cage_to_obj(cage, obj):
    _crease_props(obj)
    obj.CagePoints = [App.Vector(*p) for p in cage.V.tolist()]
    obj.CageFaces = _poly_index(cage.F)
    edges = list(cage.C.items())
    obj.CageCreases = [int(v) for e, _ in edges for v in sorted(e)]
    obj.CageSharpness = [float(s) for _, s in edges]
    wires = sorted(sorted(e) for e in cage.free_edges())
    if wires or "CageWires" in obj.PropertiesList:
        if "CageWires" not in obj.PropertiesList:
            obj.addProperty("App::PropertyIntegerList", "CageWires", "SubD", "Free edges without faces: point index pairs")
        flat = [int(v) for e in wires for v in e]
        if list(obj.CageWires) != flat:
            obj.CageWires = flat
    if cage.rails or "Rails" in obj.PropertiesList:
        if "Rails" not in obj.PropertiesList:
            obj.addProperty("App::PropertyString", "Rails", "SubD", "Rail bridge guide curves (JSON)")
            obj.setEditorMode("Rails", 2)
        txt = json.dumps(cage.rails) if cage.rails else ""
        if obj.Rails != txt:
            obj.Rails = txt
    if cage.dims or "Dimensions" in obj.PropertiesList:
        if "Dimensions" not in obj.PropertiesList:
            obj.addProperty("App::PropertyString", "Dimensions", "SubD", "Distance measurements (JSON)")
            obj.setEditorMode("Dimensions", 2)
        txt = json.dumps(cage.dims) if cage.dims else ""
        if obj.Dimensions != txt:
            obj.Dimensions = txt
    if cage.links or "SketchLinks" in obj.PropertiesList:
        _link_props(obj)
        txt = json.dumps(cage.links) if cage.links else ""
        if obj.SketchLinks != txt:
            obj.SketchLinks = txt
        doc = obj.Document
        srcs = []
        for l in cage.links:
            for name, _ in l["src"]:
                o = doc.getObject(name)
                if o is not None and o not in srcs:
                    srcs.append(o)
        if list(obj.SketchSources) != srcs:               # bagimlilik: kaynak degisince nesne yeniden hesaplanir
            obj.SketchSources = srcs
    if cage.groups or "ParamGroups" in obj.PropertiesList:
        _group_props(obj)
        txt = json.dumps(cage.groups) if cage.groups else ""
        if obj.ParamGroups != txt:
            obj.ParamGroups = txt
        drop_handle_links(obj)


def _link_props(obj):
    if "SketchLinks" not in obj.PropertiesList:
        obj.addProperty("App::PropertyString", "SketchLinks", "SubD", "Sketch links of From edge polygons (JSON)")
        obj.setEditorMode("SketchLinks", 2)
        obj.addProperty("App::PropertyLinkListGlobal", "SketchSources", "SubD", "Objects the linked edges come from")
        obj.setEditorMode("SketchSources", 1)
        obj.addProperty("App::PropertyBool", "SketchLink", "SubD",
                        "Linked points follow the source sketch / edges when they change; off = unlinked").SketchLink = True


def cage_from_obj(obj):
    V = np.array([tuple(p) for p in obj.CagePoints], float).reshape(-1, 3)
    F, f = [], []
    for i in obj.CageFaces:
        if i < 0:
            F.append(f)
            f = []
        else:
            f.append(i)
    C = {}
    if "CageCreases" in obj.PropertiesList:
        pairs = obj.CageCreases
        for k, s in enumerate(obj.CageSharpness):
            C[_edge(pairs[2 * k], pairs[2 * k + 1])] = s
    W = set()
    if "CageWires" in obj.PropertiesList:
        w = obj.CageWires
        W = {_edge(w[2 * k], w[2 * k + 1]) for k in range(len(w) // 2) if w[2 * k] != w[2 * k + 1]}
    cage = Cage(V, F, C, W)
    if "Rails" in obj.PropertiesList and obj.Rails:
        try:
            cage.rails = json.loads(obj.Rails)
        except ValueError:
            cage.rails = []
    if "Dimensions" in obj.PropertiesList and obj.Dimensions:
        try:
            cage.dims = json.loads(obj.Dimensions)
        except ValueError:
            cage.dims = []
    if "SketchLinks" in obj.PropertiesList and obj.SketchLinks:
        try:
            cage.links = json.loads(obj.SketchLinks)
        except ValueError:
            cage.links = []
    if "ParamGroups" in obj.PropertiesList and obj.ParamGroups:
        try:
            cage.groups = json.loads(obj.ParamGroups)
        except ValueError:
            cage.groups = []
    return cage


# ---------- Part / PartDesign yuzeylerinden kafes (From face) ----------
# Secili B-rep yuzeyleri tek kafes olur. Her B-rep kenari bir kez orneklenir (duz kenar az, egri kenar aciya gore:
# 90 derece basina arc_segs; max_len'den uzun parca yok) ve her yuzde AYNI noktalar kullanilir: ortak kenarlar
# kaynaklanir. Dort kenarli yuz (duz ya da kavisli: silindir, B-spline...): yuzeyin UV parametresinde Coons
# izgarasi; karsi kenarlar ayni sayida orneklenir (bu kosul bitisik yuzlere yayilir). Duz yuz: tek sinir -> iceri
# halka + ortada cokgen, dis + ic sinir (or. kutu agzi) -> dortgen serit. Yuzlerin keskin birlesimi (or. agiz ile
# duvar) ve acik sinirdaki keskin koseler crease ile keskin kalir, teget birlesim yumusak kalir. exact: kafes
# tekrar tekrar duzeltilir ki SubD limit yuzeyi butun ornek noktalardan gecsin (sinir ve ic).

def _edge_turn(e):
    """Kenar boyunca teget donusu (radyan); duz kenarda ~0."""
    a, b = e.FirstParameter, e.LastParameter
    ts = [e.tangentAt(a + (b - a) * k / 16.0) for k in range(17)]
    ang = 0.0
    for t0, t1 in zip(ts, ts[1:]):
        if t0.Length > 1e-12 and t1.Length > 1e-12:
            ang += t0.getAngle(t1)
    return ang


def _edge_count(e, arc_segs, max_len):
    turn = _edge_turn(e)
    n = int(math.ceil(turn / (math.pi / 2) * arc_segs - 1e-9)) if turn > 1e-3 else 1
    if max_len > 0:
        n = max(n, int(math.ceil(e.Length / max_len - 1e-9)))
    return max(n, 1)


def _edge_points(e, n):
    """Kenarin n parcali ornegi, PARAMETRE yonunde. Her yuzde ayni cagrilir: ortak kenarin noktalari ayni."""
    return np.array([tuple(v) for v in e.discretize(Number=n + 1)], float)


def _chain(edges, counts, rev=None):
    """Sirali kapali kenar zinciri: her kenarin ornegi ZINCIR yonunde ((n+1, 3) dizi, uclar dahil). Kenarin
    parametre yonu zincirle ters olabilir (komsuluktan bulunur; Vertexes sirasina guvenilmez). Kapali kenarda
    (uclari ayni) komsuluk yonu vermez: rev[k] (bu yuzun telinde ters mi) ya da kenarin kendi Orientation'i."""
    ends = [(np.array(tuple(e.valueAt(e.FirstParameter))), np.array(tuple(e.valueAt(e.LastParameter))))
            for e in edges]
    nxt = ends[1 % len(edges)]                             # ilk kenarin sonu ikinci kenara degmeli
    fwd0 = min(np.linalg.norm(ends[0][1] - q) for q in nxt) <= min(np.linalg.norm(ends[0][0] - q) for q in nxt)
    prev = ends[0][0] if fwd0 else ends[0][1]
    out = []
    for k, (e, n, (p0, p1)) in enumerate(zip(edges, counts, ends)):
        q = _edge_points(e, n)
        if np.linalg.norm(p1 - p0) < 1e-7:                 # kapali kenar (silindir agzindaki cember)
            if (rev[k] if rev is not None else e.Orientation == "Reversed"):
                q = q[::-1]
        elif np.linalg.norm(p1 - prev) < np.linalg.norm(p0 - prev):
            q = q[::-1]
        out.append(q)
        prev = q[-1]
    return out


def auto_max_len(edges, arc_segs):
    """Otomatik en uzun parca: egri kenarlardaki adimin 1.5 kati (duz kenarlar da benzer aralikla orneklenmeli,
    yoksa noktalardan gecmeye zorlanan SubD uzun duz parcada dalgalanir); toplam ~400 parcayi asmasin. Egri kenar
    yoksa 0 (duz kenar tek parca)."""
    steps = [e.Length / _edge_count(e, arc_segs, 0) for e in edges if _edge_turn(e) > 1e-3]
    if not steps:
        return 0.0
    return max(1.5 * min(steps), sum(e.Length for e in edges) / 400.0)


def _face_normal(face, p):
    u, v = face.Surface.parameter(App.Vector(*p))
    return np.array(tuple(face.normalAt(u, v)))


def _pair_order(eo, ei):
    """Halka yuzun dis ve ic kenarlarini esle (ayni sayida): ic listenin donus + yonu, orta noktalara gore."""
    mid = lambda e: np.array(tuple(e.valueAt((e.FirstParameter + e.LastParameter) / 2)))
    mo, mi = [mid(e) for e in eo], [mid(e) for e in ei]
    best = None
    for rev in (False, True):
        lst = list(range(len(ei)))[::-1] if rev else list(range(len(ei)))
        for k in range(len(ei)):
            order = lst[k:] + lst[:k]
            d = sum(np.linalg.norm(mo[j] - mi[order[j]]) for j in range(len(eo)))
            if best is None or d < best[0]:
                best = (d, order)
    return best[1]


def _grid_patch(face, edges, counts, rev=None):
    """Dort kenarli yuz: kenar ornekleri arasinda UV'de Coons izgarasi, ic noktalar yuzeyin uzerinde. (V, F)"""
    segs = _chain(edges, counts, rev)
    n0, n1 = len(segs[0]) - 1, len(segs[1]) - 1
    if len(segs[2]) - 1 != n0 or len(segs[3]) - 1 != n1:
        return None
    B, R, Tp, L = segs[0], segs[1], segs[2][::-1], segs[3][::-1]   # alt, sag, ust (soldan saga), sol (alttan uste)
    S = face.Surface
    uv = lambda P: np.array([tuple(S.parameter(App.Vector(*p))) for p in P], float)
    UB, UR, UT, UL = uv(B), uv(R), uv(Tp), uv(L)
    for k, per in enumerate((S.UPeriod() if S.isUPeriodic() else 0, S.VPeriod() if S.isVPeriodic() else 0)):
        if per:                                            # periyodik yuzey (silindir): her kenar boyunca surekli
            for A, ref in ((UB, None), (UL, UB[0]), (UR, UB[-1]), (UT, UL[-1])):   # (tam cember 0 -> 2pi),
                A[:, k] = np.unwrap(A[:, k], period=per)   # baslangici komsu kosesiyle ayni turda
                if ref is not None:
                    A[:, k] += per * np.round((ref[k] - A[0, k]) / per)
    idx = lambda i, j: j * (n0 + 1) + i
    V = np.zeros(((n0 + 1) * (n1 + 1), 3))
    for j in range(n1 + 1):
        t = j / float(n1)
        for i in range(n0 + 1):
            s = i / float(n0)
            if j == 0:
                p = B[i]
            elif j == n1:
                p = Tp[i]
            elif i == 0:
                p = L[j]
            elif i == n0:
                p = R[j]
            else:
                w = ((1 - t) * UB[i] + t * UT[i] + (1 - s) * UL[j] + s * UR[j]
                     - ((1 - s) * (1 - t) * UB[0] + s * (1 - t) * UB[-1] + (1 - s) * t * UT[0] + s * t * UT[-1]))
                p = np.array(tuple(S.value(float(w[0]), float(w[1]))))
            V[idx(i, j)] = p
    F = [[idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)] for j in range(n1) for i in range(n0)]
    q = V[F[len(F) // 2]]                                  # yon: yuzun (disari) normaliyle ayni
    if np.cross(q[1] - q[0], q[3] - q[0]).dot(_face_normal(face, q.mean(0))) < 0:
        F = [f[::-1] for f in F]
    return V, F


def _planar_patch(face, wires):
    """Duz yuz. wires: [(kenarlar, sayilar)], ilki dis sinir. Tek sinir: iceri halka + ortada cokgen; dis + ic
    sinir: dortgen serit. (V, F)"""
    c = np.array(tuple(face.CenterOfMass))
    n = _face_normal(face, c)
    loop = lambda ec: np.vstack([q[:-1] for q in _chain(*ec)])
    ccw = lambda P: np.cross(P - c, np.roll(P - c, -1, 0)).sum(0).dot(n) > 0
    Po = loop(wires[0])
    if not ccw(Po):
        Po = Po[::-1]
    m = len(Po)
    if len(wires) == 1:
        R = c + 0.6 * (Po - c)
        F = [[i, (i + 1) % m, m + (i + 1) % m, m + i] for i in range(m)] + [list(range(m, 2 * m))]
        return np.vstack([Po, R]), F
    Pi = loop(wires[1])
    if len(Pi) != m:                                       # farkli yapida sinirlar: yay uzunluguna gore esit sayi
        m = max(m, len(Pi))
        Po, Pi = _resample(Po, m), _resample(Pi, m)
    if not ccw(Pi):
        Pi = Pi[::-1]
    Pi = np.roll(Pi, -int(np.argmin(np.linalg.norm(Pi - Po[0], axis=1))), 0)
    return np.vstack([Po, Pi]), [[i, (i + 1) % m, m + (i + 1) % m, m + i] for i in range(m)]


def _area2(P):
    """2B kapali cokgenin isaretli alani (saat yonu tersi > 0)."""
    return 0.5 * float(np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1]))


def _cross2(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def _seg_cross(p, q, a, b):
    """pq ile ab kesisiyor mu (ortak uc noktasi sayilmaz)."""
    d1, d2 = _cross2(a, b, p), _cross2(a, b, q)
    d3, d4 = _cross2(p, q, a), _cross2(p, q, b)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)) and min(abs(d1), abs(d2), abs(d3), abs(d4)) > 1e-12


def _offset2(P, d):
    """Kapali 2B halkayi soluna (bolge tarafina) d kadar kaydir: kose normali komsu kenarlarin ortalamasi."""
    t = np.roll(P, -1, 0) - P
    t /= np.maximum(np.linalg.norm(t, axis=1), 1e-12)[:, None]
    nl = np.stack([-t[:, 1], t[:, 0]], 1)                  # kenarin sol normali
    nv = nl + np.roll(nl, 1, 0)
    nv /= np.maximum(np.linalg.norm(nv, axis=1), 1e-12)[:, None]
    cosh = np.maximum((nv * nl).sum(1), 0.5)               # sivri kosede asiri uzamasin
    return P + nv * (d / cosh)[:, None]


def _earcut(X, poly):
    """Basit cokgenin (X noktalari, poly indeks listesi, saat yonu tersi; kopru noktalari iki kez gecebilir) kulak
    kirpma ucgenleri."""
    poly, tris = list(poly), []
    while len(poly) > 3:
        n, cut = len(poly), None
        for k in range(n):
            a, b, c = poly[k - 1], poly[k], poly[(k + 1) % n]
            if _cross2(X[a], X[b], X[c]) <= 1e-12:
                continue
            if any(p not in (a, b, c) and _cross2(X[a], X[b], X[p]) >= 0 and _cross2(X[b], X[c], X[p]) >= 0
                   and _cross2(X[c], X[a], X[p]) >= 0 for p in poly):
                continue
            cut = k
            break
        if cut is None:                                    # sayisal sorun: en disbukey koseyi kirp
            cut = max(range(n), key=lambda k: _cross2(X[poly[k - 1]], X[poly[k]], X[poly[(k + 1) % n]]))
        tris.append([poly[cut - 1], poly[cut], poly[(cut + 1) % n]])
        del poly[cut]
    tris.append(poly)
    return tris


def _flip_delaunay(X, tris):
    """Iki ucgenin ortak kenarini, karsi acilar toplami 180 dereceyi asiyorsa ve dortgen disbukeyse cevir."""
    def ang(o, a, b):
        u, v = X[a] - X[o], X[b] - X[o]
        return math.atan2(abs(u[0] * v[1] - u[1] * v[0]), u.dot(v))
    for _ in range(50):
        ef = {}
        for ti, t in enumerate(tris):
            for k in range(3):
                ef.setdefault(frozenset((t[k], t[k - 2])), []).append(ti)
        flipped, used = False, set()
        for e, ts in ef.items():
            if len(ts) != 2 or ts[0] in used or ts[1] in used:
                continue
            t1, t2 = tris[ts[0]], tris[ts[1]]
            c = next(v for v in t1 if v not in e)
            d = next(v for v in t2 if v not in e)
            k = t1.index(c)
            a, b = t1[(k + 1) % 3], t1[(k + 2) % 3]        # t1 = (c, a, b), t2 = (b, a, d) saat yonu tersi
            if ang(c, a, b) + ang(d, b, a) <= math.pi + 1e-9:
                continue
            if _cross2(X[c], X[a], X[d]) <= 1e-12 or _cross2(X[d], X[b], X[c]) <= 1e-12:
                continue
            tris[ts[0]], tris[ts[1]] = [c, a, d], [d, b, c]
            used.update(ts)
            flipped = True
        if not flipped:
            break
    return tris


def _merge_quads(X, tris):
    """Komsu ucgen ciftlerini, olusan dortgen disbukey ve acilari 90 dereceye yakinsa dortgene birlestir."""
    ef = {}
    for ti, t in enumerate(tris):
        for k in range(3):
            ef.setdefault(frozenset((t[k], t[k - 2])), []).append(ti)
    cand = []
    for e, ts in ef.items():
        if len(ts) != 2:
            continue
        t1, t2 = tris[ts[0]], tris[ts[1]]
        c = next(v for v in t1 if v not in e)
        d = next(v for v in t2 if v not in e)
        k = t1.index(c)
        a, b = t1[(k + 1) % 3], t1[(k + 2) % 3]
        q = [a, d, b, c]
        dev = 0.0
        for i in range(4):
            o, p, r = X[q[i]], X[q[i - 1]], X[q[(i + 1) % 4]]
            if _cross2(p, o, r) <= 1e-12:
                dev = None
                break
            u, v = p - o, r - o
            dev = max(dev, abs(math.degrees(math.atan2(abs(u[0] * v[1] - u[1] * v[0]), u.dot(v))) - 90))
        if dev is not None and dev < 55:
            cand.append((dev, ts[0], ts[1], q))
    used, faces = set(), []
    for dev, i, j, q in sorted(cand, key=lambda c: c[0]):
        if i not in used and j not in used:
            used.update((i, j))
            faces.append(q)
    return faces + [t for ti, t in enumerate(tris) if ti not in used]


def _holes_patch(face, wires):
    """Birden fazla delikli duz yuz: her deligin cevresine dortgen halka, dis sinir ile halkalar arasi kulak
    kirpma + Delaunay cevirme ile ucgenlenip dortgenlere birlestirilir. (V, F) ya da None"""
    c = np.array(tuple(face.CenterOfMass))
    n = _face_normal(face, c)
    ax = np.eye(3)[int(np.argmin(np.abs(n)))]
    u = np.cross(n, ax)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    to2 = lambda P: np.stack([(P - c).dot(u), (P - c).dot(v)], 1)
    loops3 = [np.vstack([q[:-1] for q in _chain(*w)]) for w in wires]
    loops2 = [to2(P) for P in loops3]
    for k in range(len(loops3)):                           # dis sinir saat yonu tersi, delikler saat yonu
        if (_area2(loops2[k]) > 0) != (k == 0):
            loops3[k], loops2[k] = loops3[k][::-1], loops2[k][::-1]
    # delik halkalari: deligin boyutu, komsu sinirlara aciklik ve ornek araligina gore genislik
    X2, X3, idx = [], [], []
    for P2, P3 in zip(loops2, loops3):
        idx.append(list(range(len(X2), len(X2) + len(P2))))
        X2 += list(P2)
        X3 += list(P3)
    region = [idx[0]]
    F = []
    for k in range(1, len(loops2)):
        P = loops2[k]
        others = np.vstack([loops2[j] for j in range(len(loops2)) if j != k])
        gap = min(np.linalg.norm(others - p, axis=1).min() for p in P)
        r = np.linalg.norm(P - P.mean(0), axis=1).mean()
        step = np.linalg.norm(np.roll(P, -1, 0) - P, axis=1).mean()
        d = min(0.35 * r, 0.3 * gap, step)
        if d <= 1e-9:
            return None
        O = _offset2(P, d)
        oi = list(range(len(X2), len(X2) + len(O)))
        X2 += list(O)
        X3 += [c + p[0] * u + p[1] * v for p in O]
        h, m = idx[k], len(O)
        F += [[h[i], h[(i + 1) % m], oi[(i + 1) % m], oi[i]] for i in range(m)]
        region.append(oi)
    X = np.array(X2)
    # delikleri kopruyle dis sinira bagla: en yakin, hicbir sinir kenarini kesmeyen nokta cifti
    poly, holes = list(region[0]), region[1:]
    while holes:
        segs = [(poly[i], poly[(i + 1) % len(poly)]) for i in range(len(poly))]
        for hl in holes:
            segs += [(hl[i], hl[(i + 1) % len(hl)]) for i in range(len(hl))]
        best = None
        for hi, hl in enumerate(holes):
            for b in hl:
                for ai, a in enumerate(poly):
                    dd = np.linalg.norm(X[a] - X[b])
                    if best is not None and dd >= best[0]:
                        continue
                    if any(_seg_cross(X[a], X[b], X[s], X[t]) for s, t in segs):
                        continue
                    best = (dd, hi, b, ai)
        if best is None:
            return None
        _, hi, b, ai = best
        hl = holes.pop(hi)
        k = hl.index(b)
        ring = hl[k:] + hl[:k]
        poly = poly[:ai + 1] + ring + [b, poly[ai]] + poly[ai + 1:]
    tris = _flip_delaunay(X, _earcut(X, poly))
    F += _merge_quads(X, tris)
    return np.array(X3), F


def _resample(P, m):
    """Kapali cokgeni yay uzunluguna gore m noktaya yeniden orneklendir."""
    seg = np.linalg.norm(np.roll(P, -1, 0) - P, axis=1)
    L = np.concatenate([[0], np.cumsum(seg)])
    t = np.linspace(0, L[-1], m, endpoint=False)
    return np.stack([np.interp(t, L, np.concatenate([P[:, d], P[:1, d]])) for d in range(3)], 1)


def cage_from_faces(faces, arc_segs=4, max_len=0.0, exact=True):
    """Secili yuzlerden (dunya koordinatinda) tek kafes ya da hata metni. max_len <= 0: auto_max_len."""
    edges, ekeys = [], {}

    def eid(e):
        # ayni kenar geometriyle taninir (uclar + orta nokta): Part.getShape her yuzu ayri kopyalar, isSame tutmaz
        r = lambda p: tuple(np.round(np.array(tuple(p)) / 1e-4).astype(np.int64))
        k = (frozenset((r(e.valueAt(e.FirstParameter)), r(e.valueAt(e.LastParameter)))),
             r(e.valueAt((e.FirstParameter + e.LastParameter) / 2)))
        if k not in ekeys:
            ekeys[k] = len(edges)
            edges.append(e)
        return ekeys[k]
    def reversed_in(o, e):
        """Yuzun kendi kenari o, ortak tablodaki e'ye gore bu yuzun telinde ters mi (kapali kenarlar icin)."""
        same = np.dot(tuple(o.tangentAt(o.FirstParameter)), tuple(e.tangentAt(e.FirstParameter))) > 0
        return (o.Orientation == "Reversed") == same
    info, revs = [], []
    for f in faces:
        planar = f.isPlanar() if hasattr(f, "isPlanar") else f.Surface.TypeId == "Part::GeomPlane"
        ws = [f.OuterWire] + [w for w in f.Wires if not w.isSame(f.OuterWire)]
        wl = [[eid(e) for e in w.OrderedEdges] for w in ws]
        revs.append([{i: reversed_in(o, edges[i]) for i, o in zip(w, ow.OrderedEdges)} for w, ow in zip(wl, ws)])
        if len(ws) == 1 and len(wl[0]) == 4:
            kind = "grid"
        elif planar and len(ws) <= 2:
            kind = "planar"
        elif planar:
            kind = "holes"
        else:
            return T("From face: curved faces need exactly four edges")
        info.append((f, kind, wl))
    if max_len <= 0:
        max_len = auto_max_len(edges, arc_segs)
    # parca sayilari: dort kenarli yuzde karsi kenarlar, halka yuzde eslenen dis / ic kenarlar ayni sayida
    parent = list(range(len(edges)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for f, kind, wl in info:
        if kind == "grid":
            parent[find(wl[0][0])] = find(wl[0][2])
            parent[find(wl[0][1])] = find(wl[0][3])
        elif len(wl) == 2 and len(wl[0]) == len(wl[1]):
            order = _pair_order([edges[i] for i in wl[0]], [edges[i] for i in wl[1]])
            wl[1][:] = [wl[1][j] for j in order]
            for a, b in zip(wl[0], wl[1]):
                parent[find(a)] = find(b)
    best = {}
    for i, e in enumerate(edges):
        best[find(i)] = max(best.get(find(i), 1), _edge_count(e, arc_segs, max_len))
    cnt = [best[find(i)] for i in range(len(edges))]
    # yamalar ve kaynak (ortak kenar noktalari ayni oldugundan konumla birlesir)
    Vs, Fs, off = [], [], 0
    for (f, kind, wl), rv in zip(info, revs):
        ws = [([edges[i] for i in w], [cnt[i] for i in w], [r[i] for i in w]) for w, r in zip(wl, rv)]
        res = (_grid_patch(f, *ws[0]) if kind == "grid" else
               _holes_patch(f, ws) if kind == "holes" else _planar_patch(f, ws))
        if res is None:
            return T("From face: this face could not be meshed")
        Vs.append(res[0])
        Fs += [[v + off for v in q] for q in res[1]]
        off += len(res[0])
    V = np.vstack(Vs)
    tol = 1e-6 * max(np.ptp(V, 0).max(), 1.0)
    key = lambda p: tuple(np.round(np.asarray(p, float) / tol).astype(np.int64))
    index, uniq, remap = {}, [], []
    for p in V:
        k = key(p)
        if k not in index:
            index[k] = len(uniq)
            uniq.append(p)
        remap.append(index[k])
    cage = Cage(np.array(uniq), [[remap[v] for v in q] for q in Fs])
    # keskinlik: iki yuzun keskin birlesimi (normaller > 10 derece) ve acik sinirdaki keskin koseler
    use = {}
    for fi, (f, kind, wl) in enumerate(info):
        for w in wl:
            for i in w:
                use.setdefault(i, []).append(fi)
    C, sharp = {}, math.cos(math.radians(10))
    for i, fl in use.items():
        if len(fl) != 2:
            continue
        e = edges[i]
        m = tuple(e.valueAt((e.FirstParameter + e.LastParameter) / 2))
        if _face_normal(info[fl[0]][0], m).dot(_face_normal(info[fl[1]][0], m)) < sharp:
            ids = [index.get(key(p)) for p in _edge_points(e, cnt[i])]
            for a, b in zip(ids, ids[1:]):
                if a is not None and b is not None and a != b:
                    C[_edge(a, b)] = CREASE_HARD
    ends = {}                                              # acik sinir noktasi -> oradan cikan teget yonleri
    for i, fl in use.items():
        if len(fl) == 1:
            e = edges[i]
            for prm, sgn in ((e.FirstParameter, 1.0), (e.LastParameter, -1.0)):
                k = index.get(key(tuple(e.valueAt(prm))))
                t = np.array(tuple(e.tangentAt(prm))) * sgn
                if k is not None and np.linalg.norm(t) > 1e-12:
                    ends.setdefault(k, []).append(t / np.linalg.norm(t))
    ef = cage.edge_faces()
    for k, ts in ends.items():
        if len(ts) == 2 and ts[0].dot(ts[1]) > -sharp:     # yonler zit degil: kose (tek yuzlu kose zaten sabit)
            at = [e for e in ef if k in e]
            for e in at:
                if len(ef[e]) == 1:
                    C[e] = CREASE_HARD
            inner = [e for e in at if len(ef[e]) == 2]
            if inner:
                C[inner[0]] = CREASE_HARD                  # 3 keskin kenar: kose noktasi yerinde kalir
    cage.C = C
    if exact:                                              # SubD limit yuzeyi ornek noktalardan gecsin
        # duz yuzun ic noktalari duzeltilmez: limit yuzeyi zaten duzlemde; zorlanirsa duzensiz agda yuzler katlanir
        fit = np.zeros(len(cage.V), bool)
        for i, e in enumerate(edges):
            for p in _edge_points(e, cnt[i]):
                k = index.get(key(p))
                if k is not None:
                    fit[k] = True
        off = 0
        for (f, kind, wl), nv in zip(info, [len(x) for x in Vs]):
            if kind == "grid":
                fit[remap[off:off + nv]] = True
            off += nv
        T0, P = cage.V.copy(), cage.V.copy()
        lv = 4 if len(P) < 600 else 3
        for _ in range(20):
            Sx, _ = subdivide(P, cage.F, lv, cage.C)
            P = P + np.where(fit[:, None], T0 - Sx[:len(P)], 0.0)
        cage.V = P
    return cage


# ---------- Kenarlardan / noktalardan cokgen (From face > Edge / Point) ----------
# Edge: uclari (CHAIN_TOL icinde) cakisan B-rep kenarlari tek zincir olur (farkli objelerden iki yarim daire de);
# her zincir tek cokgen: kapali zincir oldugu gibi, acik zincirin uclari duz kenarla kapanir (Max Length'e gore
# bolunur). Ayri yerlerdeki zincirler ayri cokgenler. Duzlemsel olmasi gerekmez. Point: secim sirasiyla tek cokgen.
CHAIN_TOL = 0.01         # mm: bu kadar yakin kenar uclari ayni nokta sayilir


def _edge_chains(edges, tol=CHAIN_TOL):
    """Kenarlari uclarindan zincirle: [([(kenar, ters mi)], kapali mi)] ya da hata metni (dallanan secim)."""
    nodes = []

    def node(p):
        for k, q in enumerate(nodes):
            if np.linalg.norm(p - q) <= tol:
                return k
        nodes.append(p)
        return len(nodes) - 1
    ends = [(node(np.array(tuple(e.valueAt(e.FirstParameter)))), node(np.array(tuple(e.valueAt(e.LastParameter)))))
            for e in edges]
    adj = {}
    for i, (a, b) in enumerate(ends):
        if a != b:
            adj.setdefault(a, []).append(i)
            adj.setdefault(b, []).append(i)
    if any(len(l) > 2 for l in adj.values()):
        return T("From face: the selected edges branch (more than two meet at a point)")
    chains, used = [], set()
    for i, (a, b) in enumerate(ends):
        if a == b:                                         # kendi basina kapali kenar (cember)
            chains.append(([(edges[i], False)], True))
            used.add(i)
    starts = [n for n, l in adj.items() if len(l) == 1] + list(adj)   # once acik zincirler, sonra halkalar
    for n0 in starts:
        for i in adj[n0]:
            if i in used:
                continue
            chain, n = [], n0
            while i is not None:
                used.add(i)
                a, b = ends[i]
                rev = a != n
                chain.append((edges[i], rev))
                n = a if rev else b
                i = next((k for k in adj[n] if k not in used), None)
            chains.append((chain, n == n0))
    return chains


def _on_line(P):
    """Uc noktadan az ya da hepsi bir dogru uzerinde (capraz giden cokgenin alani 0 olabilir, o sayilmaz)."""
    if len(P) < 3:
        return True
    s = np.linalg.svd(P - P.mean(0), compute_uv=False)
    return s[1] <= 1e-9 * max(s[0], 1.0)


def cage_from_edges(edges, arc_segs=4, max_len=0.0, counts=None):
    """Secili kenarlardan (dunya) zincir basina bir cokgen: Cage ya da hata metni. Duz kenar (ve acik zinciri
    kapatan kenar) tek parca; egri kenar Arc Segments'e gore (max_len <= 0: auto_max_len, yalniz egrilerde).
    counts (kenar basina parca) verilirse o kullanilir (Sketch bagi: ayni noktalar ayni sirada); kullanilan sayilar
    donen kafesin .counts'una yazilir."""
    chains = _edge_chains(edges)
    if isinstance(chains, str):
        return chains
    if max_len <= 0:
        max_len = auto_max_len(edges, arc_segs)
    ecount = list(counts) if counts else \
        [_edge_count(e, arc_segs, max_len) if _edge_turn(e) > 1e-3 else 1 for e in edges]
    pos = {id(e): k for k, e in enumerate(edges)}
    V, F = [], []
    for chain, closed in chains:
        P = []
        for e, rev in chain:
            q = _edge_points(e, ecount[pos[id(e)]])
            q = list(q[::-1] if rev else q)
            P += q[1:] if P else q
        if closed:
            P = P[:-1]                                     # acik zincir: uclari tek duz kenarla kapanir
        P = np.array(P)
        if _on_line(P):
            continue                                       # tek duz kenar: alani yok
        F.append(list(range(len(V), len(V) + len(P))))
        V += list(P)
    if not F:
        return T("From face: no polygon can be made from the selected edges")
    c = Cage(np.array(V), F)
    c.counts = ecount
    return c


def cage_from_points(points):
    """Noktalardan (dunya, secim sirasinda) tek cokgen: Cage ya da hata metni. Ust uste gelen ardisik noktalar
    tek sayilir."""
    P = []
    for p in points:
        if not P or np.linalg.norm(p - P[-1]) > CHAIN_TOL:
            P.append(np.asarray(p, float))
    if len(P) > 1 and np.linalg.norm(P[0] - P[-1]) <= CHAIN_TOL:
        P.pop()
    P = np.array(P)
    if _on_line(P):
        return T("From face: select 3 or more points that are not on one line")
    return Cage(P, [list(range(len(P)))])


# ---------- Sketch bagi (From face > Edge) ----------
# From edge ile kurulan cokgen kaynak kenarlarini hatirlar (Cage.links, belgede SketchLinks); kaynak (or. Sketcher
# cizimi) degisince:
# - kafes hala yalniz o duz cokgense bastan kurulur (nokta sayisi da yeni egriye gore);
# - kafes duzenlendiyse (extrude, bevel...) bagli noktalar egrinin ayni orneklerine gider ve cokgen duzleminde
#   tam ustlerinde / altlarinda duran noktalar (extrude kopyalari) ayni kadar kayar. Diger noktalar (bevel / inset
#   halkalari) yerinde kalir; hesap yok, yalniz tasima. Bagli nokta, kaydedilen yerinde duran kafes noktasidir:
#   elle tasinan nokta bagdan cikar, hicbiri yerinde degilse bag biter (kaydi silinir).
# Kayit (yerel koordinat): src [[nesne, alt eleman]], segs, max, counts (kenar basina parca), pts (bagli
# noktalarin son yeri).

LINK_ENDED = "ended"


def link_new(src, vals, cage, verts, counts):
    """From edge kafesinden bag kaydi: verts = cokgenin noktalari (yerel, kafesteki indeksleri, sirayla); counts:
    cage_from_edges'in kenar basina parca sayilari."""
    return dict(src=[list(s) for s in src], segs=int(vals["Arc Segments"]), max=float(vals["Max Length"]),
                counts=list(counts), pts=cage.V[verts].tolist())


def link_edges(link, doc):
    """Bagin kaynak kenarlari (dunya) ya da hata metni (nesne silindi / kenar yok)."""
    import Part
    out = []
    for name, sub in link["src"]:
        o = doc.getObject(name)
        if o is None:
            return T("Sketch link: %s is gone") % name
        try:
            sh = Part.getShape(o, sub, needSubElement=True)
        except Exception:
            sh = None
        if sh is None or sh.isNull() or sh.ShapeType != "Edge":
            return T("Sketch link: %s.%s is gone") % (o.Label, sub)
        out.append(sh)
    return out


def link_update(cage, link, doc, to_local, alone=True):
    """Bagi kaynagin simdiki haline getirir (kafes ve kayit degisir). Dondurur: None (degisiklik yok), True
    (degisti), LINK_ENDED (bagli noktalarin hicbiri yerinde degil: bag bitti) ya da hata metni."""
    old = np.array(link["pts"], float).reshape(-1, 3)
    if not len(old) or not len(cage.V):
        return LINK_ENDED
    tol = 1e-6 * max(np.ptp(np.vstack([cage.V, old]), 0).max(), 1.0)
    keep = [k for k, p in enumerate(old) if (np.linalg.norm(cage.V - p, axis=1) < tol).any()]
    if not keep:
        return LINK_ENDED
    edges = link_edges(link, doc)
    if isinstance(edges, str):
        return edges
    if alone and len(cage.V) == len(old) and len(keep) == len(old) and not cage.C:
        # duzenlenmemis duz cokgen: bastan kur (egri uzayip kisalinca nokta sayisi da ayarlanir)
        fresh = cage_from_edges(edges, link["segs"], link["max"])
        if isinstance(fresh, str):
            return fresh
        V = to_local(fresh.V)
        if len(V) == len(old) and np.allclose(V, old, atol=1e-9, rtol=0):
            return None
        cage.V, cage.F, cage.C = V, _match_faces(cage, fresh.F, V), {}
        link.update(counts=list(fresh.counts), pts=V.tolist())
        return True
    counts = link.get("counts")
    if isinstance(counts, dict):                           # onceki surumun kaydi
        counts = counts.get("edges")
    c = cage_from_edges(edges, link["segs"], link["max"], counts)
    if isinstance(c, str):
        return c
    new = to_local(c.V)
    if len(new) != len(old):
        return T("Sketch link: the curve no longer gives the same points")
    if np.allclose(new[keep], old[keep], atol=1e-9, rtol=0):
        return None
    n = np.linalg.svd(old - old.mean(0), full_matrices=False)[2][-1]   # cokgen normali: kopyalar bu yonde
    flat = lambda X: X - np.outer(X.dot(n), n)
    Pf, Of = flat(cage.V), flat(old)
    P = cage.V.copy()
    for k in keep:
        P[np.linalg.norm(Pf - Of[k], axis=1) < tol] += new[k] - old[k]
    cage.V = P
    pts = old.copy()
    pts[keep] = new[keep]
    link["pts"] = pts.tolist()
    return True


def _match_faces(cage, F, V):
    """Yeni cokgenler eskilerle ayni yone baksin (From edge kameraya cevirmisti): her biri merkezine en yakin
    eski cokgenin yonunu alir."""
    out = []
    for f in F:
        P = V[f]
        n = np.cross(P - P.mean(0), np.roll(P - P.mean(0), -1, 0)).sum(0)
        best = min(cage.F, key=lambda g: np.linalg.norm(cage.V[g].mean(0) - P.mean(0)))
        Q = cage.V[best]
        m = np.cross(Q - Q.mean(0), np.roll(Q - Q.mean(0), -1, 0)).sum(0)
        out.append(f[::-1] if n.dot(m) < 0 else f)
    return out


def links_update(cage, doc, to_local):
    """Kafesin butun baglarini kaynaklarina getirir, biten baglari siler: (degisti mi, [hata metinleri], bitti mi)."""
    changed, errs, ended, keep = False, [], False, []
    for link in cage.links:
        r = link_update(cage, link, doc, to_local, alone=len(cage.links) == 1)
        if r == LINK_ENDED:
            ended = True
            continue
        keep.append(link)
        if isinstance(r, str):
            errs.append(r)
        elif r:
            changed = True
    cage.links = keep
    return changed, errs, ended


def obj_links_update(obj):
    """Nesne yeniden hesaplanirken: bagli kaynak degistiyse bagli noktalar onu izler; hicbiri yerinde degilse bag
    silinir. Modelleyici bu nesnenin kafesini duzenliyorsa o yapar (SubDSession.links_refresh)."""
    if "SketchLinks" not in obj.PropertiesList or not obj.SketchLinks or not obj.SketchLink:
        return
    s = _sessions.get(obj.Document.Name)
    if s is not None and s.mode != "object" and s.obj is obj:
        _later(lambda: _quiet(s.links_refresh))
        return
    cage = cage_from_obj(obj)
    inv = obj.Placement.inverse()
    to_local = lambda X: np.array([tuple(inv.multVec(App.Vector(*p))) for p in X], float).reshape(-1, 3)
    changed, errs, ended = links_update(cage, obj.Document, to_local)
    for e in errs:
        App.Console.PrintWarning("%s: %s\n" % (obj.Label, e))
    if changed or ended:
        cage_to_obj(cage, obj)


# ---------- Parametrik grup (datum tutamagi) ----------
# Secili noktalar bir grup olur; belgede ayni konumda bir datum plane (SubD_PG_0001) olusur. Noktalar tutamagin
# YEREL koordinatinda saklanir (loc); tutamak FreeCAD'in Placement / Attachment araclariyla tasininca (dunya
# matrisi H degisince) noktalar farki izler: yeni = nokta + agirlik * (H_simdi * loc - H_son * loc). Tutamak SubD'ye
# bagli degil ve SubD de tutamaga FreeCAD bagimliligi ile bagli degil (LCS / yuz gibi SubD'ye bagli referanslarda
# dongu olurdu): tutamagin Placement'i degisince belge gozlemcisi SubD'yi gunceller. Uye nokta, kaydedilen yerinde duran
# kafes noktasidir (Sketch bagi gibi): elle tasinan nokta gruptan cikar, hicbiri yerinde degilse grup biter.
# Kayit: name (datum), pts (uyelerin son yeri, nesne yerel), loc (tutamak yerel), w (agirlik 0..1), H (son 4x4).

GROUP_PREFIX = "SubD_PG_"


def _group_props(obj):
    if "ParamGroups" not in obj.PropertiesList:
        obj.addProperty("App::PropertyString", "ParamGroups", "SubD", "Parametric groups of points (JSON)")
        obj.setEditorMode("ParamGroups", 2)
    drop_handle_links(obj)


def drop_handle_links(obj):
    """Eski surumler SubD'yi tutamaklara baglardi (GroupHandles): tutamak SubD'nin ya da ona bagli bir seyin (LCS ...)
    uzerinden referans alinca bagimlilik dongusu olusuyordu. Bag kalktı; tutamak degisince SubD'yi gozlemci izler."""
    if "GroupHandles" in obj.PropertiesList:
        obj.removeProperty("GroupHandles")


def _mat4(pl):
    return np.array(pl.toMatrix().A, float).reshape(4, 4)


def group_color(n):
    """n. grubun rengi: altin oranla donen ton, birbirini izleyen gruplar rastgele gorunur ama hep ayirt edilir."""
    return [round(c, 3) for c in colorsys.hsv_to_rgb((0.13 + 0.61803398875 * n) % 1.0, 0.55, 0.95)]


def group_rgb(g):
    """Grup kaydinin rengi (eski kayitlarda yok: adindaki siraya gore)."""
    if g.get("color"):
        return tuple(g["color"])
    tail = str(g.get("name", ""))[len(GROUP_PREFIX):]
    return tuple(group_color(int(tail) if tail.isdigit() else 0))


# Agacin Description sutunu (Label2) grubun renginde bir kare gosterir; Label'dan ayri oldugu icin ad degisince kalir.
GROUP_SQUARES = ((0.0, "\U0001F7E5"), (0.08, "\U0001F7E7"), (0.14, "\U0001F7E8"),
                 (0.33, "\U0001F7E9"), (0.60, "\U0001F7E6"), (0.76, "\U0001F7EA"))


def group_square(rgb):
    """Renge en yakin tonlu renkli kare (kirmizi, turuncu, sari, yesil, mavi, mor)."""
    hue = colorsys.rgb_to_hsv(*rgb)[0]
    return min(GROUP_SQUARES, key=lambda hs: min(abs(hue - hs[0]), 1.0 - abs(hue - hs[0])))[1]


def group_axes(V, faces, inside):
    """Grubun yerel ekseni (nesne yereli, sutunlar x y z; z = dista bakan normal) ya da None (belirsizse: nesnenin
    eksenleri). V: grubun noktalari, faces: secili cokgenlerin kose koordinatlari (nokta / kenar seciminde bos:
    en iyi uyan duzlem), inside: kafesin ortasi (normal buna ters yone bakar). x = yukari x normal (duvar boyunca
    yatay), y = normal x x (duvar boyunca yukari)."""
    V = np.asarray(V, float).reshape(-1, 3)
    n = np.zeros(3)
    tot = 0.0
    for P in faces:
        P = np.asarray(P, float)
        a = np.cross(P - P.mean(0), np.roll(P - P.mean(0), -1, 0)).sum(0)   # alan agirlikli (Newell) normal
        n += a
        tot += np.linalg.norm(a)
    if tot > 1e-12:
        if np.linalg.norm(n) < 0.3 * tot:                  # normaller birbirini goturuyor (silindir halkasi gibi)
            return None
    elif len(V) >= 3:
        u, sv, vt = np.linalg.svd(V - V.mean(0), full_matrices=False)
        if sv[1] < 1e-9 * max(sv[0], 1e-12):                # tek dogru uzerinde
            return None
        n = vt[2]
        if np.dot(n, V.mean(0) - np.asarray(inside, float)) < 0:
            n = -n
    else:
        return None
    ln = np.linalg.norm(n)
    if ln < 1e-12:
        return None
    n = n / ln
    up = np.array([0.0, 0.0, 1.0])
    if abs(n @ up) > 0.95:                                  # tavan / taban: x nesnenin X'i
        x = np.array([1.0, 0.0, 0.0]) - n * n[0]
    else:
        x = np.cross(up, n)
    x = x / np.linalg.norm(x)
    return np.c_[x, np.cross(n, x), n]


def group_new(name, pts, loc, H, color=None):
    return dict(name=name, pts=np.asarray(pts, float).tolist(), loc=np.asarray(loc, float).tolist(),
                w=[1.0] * len(pts), H=np.asarray(H, float).flatten().tolist(), color=list(color) if color else None)


def group_index(cage, g):
    """Grubun uyelerinin kafes indeksleri (pts ile ayni sirada; bulunamayan None). Kilitliyken uye, kayitli konumunda
    duran kafes noktasidir; kilit aciksa (elle duzenleniyor) kayittaki idx kullanilir (nokta sayisi degismediyse)."""
    V = cage.V
    if g.get("locked", True) is False:
        idx = g.get("idx") or []
        ok = g.get("nv") == len(V)
        return [i if ok and 0 <= i < len(V) else None for i in idx]
    old = np.array(g["pts"], float).reshape(-1, 3)
    if not len(old) or not len(V):
        return [None] * len(old)
    tol = 1e-6 * max(np.ptp(np.vstack([V, old]), 0).max(), 1.0)
    out = []
    for p in old:
        d = np.linalg.norm(V - p, axis=1)
        j = int(d.argmin())
        out.append(j if d[j] < tol else None)
    return out


def group_index_all(cage, g):
    """group_index'in esi, ama her uye icin ayni yerdeki TUM kafes noktalari (Knife / ayirma sonrasi birlesmemis
    ciftler): [[j, ...], ...]; kilit aciksa idx tek tek. Ciftlerin hepsi birlikte tasinir, karisma olmaz."""
    if g.get("locked", True) is False:
        return [[j] if j is not None else [] for j in group_index(cage, g)]
    V = cage.V
    old = np.array(g["pts"], float).reshape(-1, 3)
    if not len(old) or not len(V):
        return [[] for _ in old]
    tol = 1e-6 * max(np.ptp(np.vstack([V, old]), 0).max(), 1.0)
    return [[int(j) for j in np.where(np.linalg.norm(V - p, axis=1) < tol)[0]] for p in old]


def group_member_sets(V, groups):
    """{grup adi: uye indeks kumesi}; kafes noktalari bir kez sozluge konur (cok grup / cok nokta icin hizli)."""
    out, at = {}, None
    for g in groups:
        if g.get("locked", True) is False:
            idx = g.get("idx") or []
            out[g["name"]] = set(idx) if g.get("nv") == len(V) else set()
            continue
        if at is None:
            at = {}
            for i, r in enumerate(np.round(V, 5) + 0.0):
                at.setdefault(r.tobytes(), []).append(i)
        P = np.round(np.array(g["pts"], float).reshape(-1, 3), 5) + 0.0
        out[g["name"]] = {i for r in P for i in at.get(r.tobytes(), ())}
    return out


def group_update(cage, g, H_now, R_inv):
    """Grubu tutamagin simdiki haline getirir. Dondurur: None (degisiklik yok), True (degisti), LINK_ENDED.
    Kilit aciksa grup tutamagi izlemez (noktalar serbest duzenleniyor)."""
    if g.get("locked", True) is False:
        return None
    old = np.array(g["pts"], float).reshape(-1, 3)
    idx = group_index_all(cage, g)
    keep = [k for k, j in enumerate(idx) if j]
    if not keep:
        return LINK_ENDED
    loc, w = np.array(g["loc"], float).reshape(-1, 3), np.array(g["w"], float)
    H_last = np.array(g["H"], float).reshape(4, 4)
    lost = len(keep) != len(old)
    if lost:                                               # elle tasinan noktalar gruptan cikti
        g["pts"], g["loc"], g["w"] = old[keep].tolist(), loc[keep].tolist(), w[keep].tolist()
        old, loc, w = old[keep], loc[keep], w[keep]
        idx = [idx[k] for k in keep]
    if np.allclose(H_now, H_last, atol=1e-9, rtol=0):
        return True if lost else None
    at = lambda H: loc @ H[:3, :3].T + H[:3, 3]
    delta = (at(H_now) - at(H_last)) @ R_inv.T             # dunya farki -> nesne yereli
    V = cage.V.copy()
    done = set()
    for k, js in enumerate(idx):
        for j in js:
            if j not in done:                              # ayni yerdeki kayitlar bir noktayi iki kez surmesin
                done.add(j)
                V[j] = V[j] + w[k] * delta[k]
    cage.V = V
    g["pts"] = [V[js[0]].tolist() for js in idx]
    g["H"] = H_now.flatten().tolist()
    return True


def groups_update(cage, doc, placement):
    """Kafesin butun gruplarini tutamaklarina getirir: (degisti mi, [hata metinleri], bitti mi)."""
    changed, errs, ended, keep = False, [], False, []
    R_inv = _mat4(placement)[:3, :3].T
    for g in cage.groups:
        h = doc.getObject(g["name"])
        if h is None:
            errs.append(T("Group handle %s is gone: the group is released") % g["name"])
            ended = True
            continue
        try:
            H = _mat4(h.getGlobalPlacement())
        except Exception:
            H = _mat4(h.Placement)
        r = group_update(cage, g, H, R_inv)
        if r == LINK_ENDED:
            ended = True
            continue
        keep.append(g)
        if r:
            changed = True
    cage.groups = keep
    return changed, errs, ended


def group_handle(doc, obj, name, placement, extent=10.0, color=None):
    """SubD grubunun datum plane tutamagini belgeye ekler (Part::DatumPlane, ad = etiket, SubDGroup = SubD nesnesi)."""
    h = doc.addObject("Part::DatumPlane", name)
    h.Label = name
    h.addProperty("App::PropertyString", "SubDGroup", "SubD", "SubD object this handle belongs to").SubDGroup = obj.Name
    h.setEditorMode("SubDGroup", 1)
    h.Placement = placement
    if App.GuiUp:
        _quiet(lambda: setattr(h.ViewObject, "Visibility", True))
        for attr in ("Length", "Width"):
            _quiet(lambda attr=attr: setattr(h.ViewObject, attr, max(float(extent), 10.0)))
        if color:                                           # tutamak da grubun renginde
            _quiet(lambda: setattr(h.ViewObject, "ShapeColor", tuple(color)))
    if color:
        h.Label2 = group_square(color)                      # agacta Description sutununda renkli kare
    return h


def groups_reconcile(cage, obj):
    """Geri al / ileri al kafesi eski bir haline dondurdukten sonra grup kayitlari ile belgedeki tutamaklari
    eslestirir: kayitta olmayan tutamak (grup olusturma geri alindi) silinir, tutamagi olmayan kayit (silinen
    olusturma ileri alindi) icin tutamak yeniden kurulur (Attachment ayari kaybolur, yer kalir). Tutamagin kayitli
    yeri, noktalarin yeri ile ayni sayilir: tutamak tasinmis olsa da geri alinan noktalar yerinde kalir, bundan
    sonraki tutamak hareketleri o yerden olculur. Dondurur: mesaj listesi."""
    doc, msgs = obj.Document, []
    names = {g["name"] for g in cage.groups}
    for h in list(doc.Objects):
        if "SubDGroup" in getattr(h, "PropertiesList", []) and h.SubDGroup == obj.Name and h.Name not in names:
            name = h.Name                                  # silinen nesnenin ozelligine erisilemez
            doc.removeObject(name)
            msgs.append(T("Group handle %s removed") % name)
    for g in cage.groups:
        h = doc.getObject(g["name"])
        if h is None:
            H = np.array(g["H"], float).reshape(4, 4)
            pts = np.array(g["pts"], float).reshape(-1, 3)
            pl = App.Placement(App.Matrix(*H.flatten().tolist()))
            group_handle(doc, obj, g["name"], pl, float(np.ptp(pts, 0).max()) if len(pts) > 1 else 10.0,
                         group_rgb(g))
            msgs.append(T("Group handle %s created again") % g["name"])
        else:
            Hr = np.array(g["H"], float).reshape(4, 4)
            try:
                Hn = _mat4(h.getGlobalPlacement())
            except Exception:
                Hn = _mat4(h.Placement)
            if g.get("locked", True) and not np.allclose(Hn, Hr, atol=1e-7, rtol=0):
                # kilitli grup: noktalarla birlikte tutamak da kayitli yerine doner
                set_handle_placement(h, App.Placement(App.Matrix(*Hr.flatten().tolist())))
                try:
                    Hn = _mat4(h.getGlobalPlacement())
                except Exception:
                    Hn = _mat4(h.Placement)
            g["H"] = Hn.flatten().tolist()
    return msgs


def obj_groups_update(obj):
    """Nesne yeniden hesaplanirken: tutamak tasindiysa grubun noktalari onu izler (obj_links_update'in esi)."""
    if "ParamGroups" not in obj.PropertiesList or not obj.ParamGroups:
        return
    s = _sessions.get(obj.Document.Name)
    if s is not None and s.mode != "object" and s.obj is obj:
        _later(lambda: _quiet(s.groups_refresh))
        return
    cage = cage_from_obj(obj)
    changed, errs, ended = groups_update(cage, obj.Document, obj.Placement)
    for e in errs:
        App.Console.PrintWarning("%s: %s\n" % (obj.Label, e))
    if changed or ended:
        cage_to_obj(cage, obj)


def cage_key(obj):
    """Kafesin degisip degismedigini anlamak icin ozet (Placement haric)."""
    return (hash(tuple((p.x, p.y, p.z) for p in obj.CagePoints)), tuple(obj.CageFaces),
            tuple(getattr(obj, "CageCreases", ())), tuple(getattr(obj, "CageSharpness", ())), obj.Levels,
            obj_rules(obj), tuple(getattr(obj, "CageWires", ())))


def is_subd(obj):
    return obj is not None and "CageFaces" in getattr(obj, "PropertiesList", [])


class SubDViewProvider:
    """SubD nesnesinin gorunum sinifi: grup tutamaklari agacta SubD'nin altinda toplanir."""

    def __init__(self, vobj):
        vobj.Proxy = self
        self.Object = vobj.Object

    def attach(self, vobj):
        self.Object = vobj.Object

    def claimChildren(self):
        obj = getattr(self, "Object", None)
        if obj is None or not _alive(obj):
            return []
        return [h for h in obj.Document.Objects
                if "SubDGroup" in getattr(h, "PropertiesList", []) and h.SubDGroup == obj.Name]

    def dumps(self):
        return None

    def loads(self, state):
        return None


def ensure_view_provider(obj):
    """Eski belgelerdeki SubD nesnesine de gorunum sinifini ver (agac duzeni)."""
    if App.GuiUp and obj.ViewObject is not None and not isinstance(obj.ViewObject.Proxy, SubDViewProvider):
        _quiet(lambda: SubDViewProvider(obj.ViewObject))


class SubDObject:
    def __init__(self, obj):
        obj.addProperty("App::PropertyVectorList", "CagePoints", "SubD", "Cage points (mm)")
        obj.addProperty("App::PropertyIntegerList", "CageFaces", "SubD", "Cage faces: point indices, -1 separated")
        obj.addProperty("App::PropertyInteger", "Levels", "SubD", "SubPatch subdivision level").Levels = 3
        _crease_props(obj)
        _surface_props(obj)
        obj.Proxy = self

    def execute(self, obj):
        try:
            obj_links_update(obj)
        except Exception as e:                             # bag hatasi yuzeyi engellemesin
            App.Console.PrintWarning("%s: sketch link: %s\n" % (obj.Label, e))
        try:
            obj_groups_update(obj)
        except Exception as e:                             # grup hatasi yuzeyi engellemesin
            App.Console.PrintWarning("%s: group: %s\n" % (obj.Label, e))
        # Nesne yalniz tasindiysa (kafes ayni) yuzeyi yeniden uretme: FreeCAD mesh'i Placement ile zaten tasir.
        # Yeniden atamak hem gereksiz hem de geri almayi bozuyordu (eski mesh yeni Placement'i geri getiriyordu).
        key = cage_key(obj)
        if getattr(self, "_key", None) == key and obj.Mesh.CountPoints:
            return
        cage = cage_from_obj(obj)
        if cage.F:
            m = subd_mesh(cage, max(1, min(obj.Levels, 5)), obj_rules(obj))
            m.Placement = obj.Placement                    # yeni mesh atanirken nesnenin Placement'i sifirlanmasin
            obj.Mesh = m
            self._key = key
        elif obj.Mesh.CountPoints:                         # yalniz serbest kenar / nokta: eski yuzey kalmasin
            m = Mesh.Mesh()
            m.Placement = obj.Placement
            obj.Mesh = m
            self._key = key

    def dumps(self):
        return None

    def loads(self, state):
        return None


# ---------- Etkilesimli oturum ----------

def _sep(*nodes):
    s = coin.SoSeparator()
    for n in nodes:
        s.addChild(n)
    return s


def _graphics_view(view):
    """3B gorunumun QGraphicsView'i. Belge kapatilip yenisi acilinca yeni gorunum eskisinin bellek adresine
    gelebiliyor; o zaman PySide onbellekteki silinmis eski sarmalayiciyi dondurur ("Internal C++ object already
    deleted"). Gecersizse sarmalayiciyi onbellekten dusur ve yeniden iste."""
    gv = view.graphicsView()
    try:
        import shiboken6 as shiboken
    except ImportError:
        try:
            import shiboken2 as shiboken
        except ImportError:
            return gv
    if not shiboken.isValid(gv):
        shiboken.invalidate(gv)
        gv = view.graphicsView()
    return gv


def _scene_parent(view):
    # Bu FreeCAD surumunde yalnizca SoFCUnifiedSelection > ObjectGroup altindaki
    # dugumler ciziliyor; sahne kokune eklenenler gorunmuyor.
    sg = view.getSceneGraph()
    for i in range(sg.getNumChildren()):
        sel = sg.getChild(i)
        if sel.getTypeId().getName().getString() == "SoFCUnifiedSelection":
            for j in range(sel.getNumChildren()):
                if sel.getChild(j).getName().getString() == "ObjectGroup":
                    return sel.getChild(j)
    return sg


def _set_coords(node, pts):
    node.point.setValues(0, len(pts), pts.tolist())
    node.point.setNum(len(pts))


def _set_index(node, idx):
    node.coordIndex.setValues(0, len(idx), idx)
    node.coordIndex.setNum(len(idx))


def _poly_index(faces):
    idx = []
    for f in faces:
        idx += list(f) + [-1]
    return [int(i) for i in idx]


# ---------- Tasima gizmosu (FcYBTraytracerWB rig_gizmo.py'nin Tasi kismi) ----------

GIZMO_PX = 90.0          # varsayilan ok boyu, ekranda piksel (yakinlasip uzaklasinca sabit kalir); Options'ta ayarlanir
PICK_PX = 8.0            # tutamagi yakalama uzakligi (px)
BEVEL_WEIGHT = 0.4       # bevel / inset icin ek direnc (DRAG_GAIN ile carpilir: 0.2, Shift ile 0.06)
DRAG_GAIN = 0.5         # direnc: tasima fare hareketinin bu oraninda (1 = imleci birebir izler)
FINE_GAIN = 0.15         # Shift basiliyken ince ayar orani
GZ_RGB = ((0.94, 0.27, 0.27), (0.13, 0.77, 0.37), (0.23, 0.51, 0.96))
HOT = (1.0, 1.0, 0.6)
AXES = (App.Vector(1, 0, 0), App.Vector(0, 1, 0), App.Vector(0, 0, 1))
PLANE_A, PLANE_B = 0.22, 0.42          # duzlem karesinin iki eksendeki baslangic / bitisi (ok boyunun orani)
FREE_R = 0.07                          # merkez kupun yari boyu
MIRROR = np.array([-1.0, 1.0, 1.0])    # X = 0 duzlemine gore ayna

FALLOFFS = ("Smooth", "Linear", "Peak", "Dome")


def _falloff(t, mode):
    """Azalma egrisi: t = uzaklik / yaricap (0..1) -> etki (merkezde 1, sinirda 0). mode: FALLOFFS sirasi."""
    t = np.clip(t, 0.0, 1.0)
    if mode == 1:
        return 1.0 - t
    if mode == 2:
        return (1.0 - t) ** 2
    if mode == 3:
        return np.sqrt(1.0 - t * t)
    return 0.5 * (1.0 + np.cos(math.pi * t))
_last_box = {"Width": 100.0, "Depth": 100.0, "Height": 100.0,           # Box paneli son degerleri (mm, dilim)
             "Segments X": 1.0, "Segments Y": 1.0, "Segments Z": 1.0}
def _param():
    """Modelleyici ayarlari (FreeCAD parametreleri; dil ayariyla ayni yer)."""
    return App.ParamGet(subd_i18n.PARAM)


def _chain_order(edges):
    """Kenarlar tek bir acik zincir ya da kapali halka ise noktalarini sirayla dondurur: (sira, kapali); degilse None."""
    adj = {}
    for e in edges:
        x, y = tuple(e)
        adj.setdefault(x, []).append(y)
        adj.setdefault(y, []).append(x)
    if not adj or any(len(n) > 2 for n in adj.values()):
        return None
    ends = [x for x, n in adj.items() if len(n) == 1]
    if len(ends) not in (0, 2):
        return None
    start = min(ends) if ends else min(adj)
    order, prev, cur = [start], None, start
    while True:
        nxt = [x for x in adj[cur] if x != prev]
        if not nxt or nxt[0] == start:
            break
        order.append(nxt[0])
        prev, cur = cur, nxt[0]
    if len(order) != len(adj):
        return None
    return order, not ends


TWIST_FALLOFFS = ("Linear", "Smooth", "Ease in", "Ease out", "Peak", "Dome")
AXIS_OPS = ("array", "twist")          # ekseni cizgiyle / dugmeyle verilen son islemler (Draw axis, Axis X / Y / Z)


def _twist_weight(t, mode):
    """Twist'te acinin eksen boyunca orani: t = konum (0 baslangic .. 1 bitis) -> 0..1. Egriler FALLOFFS gibi
    adlandirilir ama eksen boyuncadir. 0-3 rampa (bir uc sabit, oteki uc tam aci): Linear, Smooth, Ease in,
    Ease out. 4-5 can (iki uc sabit, orta tam): Peak ve Dome, Magnet'in ayni adli egrileri ortadan olculur."""
    t = np.clip(t, 0.0, 1.0)
    if mode == 1:
        return 0.5 * (1.0 - np.cos(math.pi * t))
    if mode == 2:
        return t * t
    if mode == 3:
        return 1.0 - (1.0 - t) ** 2
    if mode >= 4:
        return _falloff(np.abs(2.0 * t - 1.0), 2 if mode == 4 else 3)
    return t


# --- Atanabilir kisayollar: komutlarin varsayilan tusu yok; kullanici bir tus atayabilir ---
# subd_modeler.set_user_key("twist", "Shift+T")   (bos metin atamayi kaldirir). Atama FreeCAD parametresinde
# (Key_twist / Key_array) saklanir ve acik modelleyicide hemen gecerlidir; dugme yazisi sonraki acilista gorunur.
# Atanan tus, ayni tusu kullanan yerlesik kisayolun onune gecer.
USER_KEYS = ("array", "twist", "dim")
USER_KEY_DEFAULT = {"dim": "U"}       # Distance: U (set_user_key("dim", "") ile kaldirilir, baska tus atanabilir)


def user_key(name):
    """name komutuna atanmis tus (tasinabilir metin, orn. "Shift+T"); atanmamissa ""."""
    return _param().GetString("Key_" + name, USER_KEY_DEFAULT.get(name, ""))


def set_user_key(name, seq=""):
    """Komuta tus ata ("Shift+T", "Ctrl+Alt+T" ...); bos metin atamayi kaldirir. Yazilan tusu dondurur."""
    if name not in USER_KEYS:
        raise ValueError("name must be one of %s" % (USER_KEYS,))
    seq = QtGui.QKeySequence(seq).toString(QtGui.QKeySequence.PortableText) if seq else ""
    _param().SetString("Key_" + name, seq)
    return seq


def _key_label(name):
    """Dugmenin sagindaki yazi: "Shift+T" -> "⇧T", "Ctrl+T" -> "^T"; atanmamissa bos."""
    return user_key(name).replace("Ctrl+", "^").replace("Shift+", "⇧").replace("Alt+", "⌥").replace("Meta+", "◆")


def _user_command(ev):
    """Tus olayina atanmis komut adi ya da None."""
    try:
        seq = QtGui.QKeySequence(ev.keyCombination()).toString(QtGui.QKeySequence.PortableText)
    except Exception:
        return None
    if not seq:
        return None
    return next((n for n in USER_KEYS if user_key(n) == seq), None)


def _active_wb():
    """Etkin calisma tezgahinin adi (orn. "SubDModelerWorkbench") ya da None."""
    try:
        return Gui.activeWorkbench().name()
    except Exception:
        return None


# ---------- Distance (olcum): iki nokta arasi mesafe ----------
# Kalici olculer kafesle saklanir (Cage.dims, belgede Dimensions ozelligi). Kayit: dict(a=uc, b=uc); uc (anchor) turleri:
#   ["v", i] kafes noktasi | ["e", i, j, t] kenar uzerinde (t: i'den j'ye oran) | ["f", [i, j, ...]] cokgen merkezi
#   ["w", x, y, z] dunya noktasi (baska nesne, sabit) | ["s", x, y, z] yuzey noktasi (yalniz gecici olcu)
DIM_KEEP = ("v", "e", "f", "w")
DIM_TOL = 0.1       # hedefli olcunun varsayilan toleransi (mm)
# Kayda eklenebilen alanlar: name (ad), hide (ekranda gizli), target (hedef mm), tol (tolerans mm)


def _dim_remap(d, remap):
    """Cage._compact sonrasi nokta indekslerini yeni numaralara cevir; noktasi silinen olcu kopar (None)."""
    try:
        out = []
        for a in (d["a"], d["b"]):
            k = a[0]
            if k == "v":
                out.append(["v", remap[a[1]]])
            elif k == "e":
                out.append(["e", remap[a[1]], remap[a[2]], a[3]])
            elif k == "f":
                out.append(["f", [remap[v] for v in a[1]]])
            else:
                out.append(list(a))
    except KeyError:
        return None
    d["a"], d["b"] = out
    return d


def _dim_point(a, V, to_world):
    """Olcu ucunun dunya noktasi (numpy 3) ya da None (nokta artik yok). V: kafes noktalari (yerel)."""
    k = a[0]
    try:
        if k == "v":
            return to_world(V[a[1]])
        if k == "e":
            return to_world(V[a[1]] + a[3] * (V[a[2]] - V[a[1]]))
        if k == "f":
            return to_world(V[list(a[1])].mean(0))
        return np.array(a[1:4], float)
    except (IndexError, KeyError, ValueError):
        return None


_last_shell = {"Thickness": 2.0, "Detail filter": 6.0, "Mode": 0.0, "Fill": 1.0, "Sharp": 1.0, "Flat rim": 1.0}   # Shell paneli
_knife_opts = {"split": False, "objects": False}   # Knife (line) secenekleri (ikisi kapali: yalniz keser)
_mline_opts = {"flip": False}                       # Mirror (line): kaynak yan ters
_last_reduce = {"Degrees": 5.0, "Delete points": 1.0}   # Reduce points paneli
_last_merge = {"Fixed": 0.0, "Distance": 1.0}   # Merge points paneli
MERGE_AUTO = 1e-3       # Merge points, Automatic: bu kadar mm yakin noktalar "ayni yerde" sayilir
_last_cyl = {"Diameter": 100.0, "Height": 100.0, "Segments": 16.0, "Height Segments": 1.0}   # Cylinder paneli
_clipboard = None                       # kopyalanan cokgenler (Cage, dunya koordinati); bkz. copy_faces
PLANE_EDGE = 0.2        # duzlem normali ile bakis arasindaki |cos| bundan kucukse kare yandan gorunuyor sayilir


def _zebra_texture():
    """Zebra seritleri icin doku: yansima yonune gore (kure ortam haritasi) siyah / beyaz seritler."""
    g = coin.SoGroup()
    g.addChild(coin.SoTextureCoordinateEnvironment())
    tex = coin.SoTexture2()
    w = 64
    row = bytes(255 if (x // 4) % 2 else 20 for x in range(w))  # 8 siyah + 8 beyaz serit
    tex.image.setValue(coin.SbVec2s(w, 2), 1, row * 2)
    tex.model = coin.SoTexture2.MODULATE
    g.addChild(tex)
    return g


# Matcap: yuzey normaline gore boyanmis "isikli kure" dokusu (kure ortam haritasi). Isik, parlama ve kenar kararmasi
# (sahte AO) dokunun icinde: sahne isiklarindan bagimsiz ve cok hizli. Ayar: taban rengi, ortam, ana / dolgu isigi,
# parlama (guc, siddet), kenar kararmasi, kenar parlamasi; "env": taban yerine stüdyo yansimasi (metal).
MATCAPS = {
    "Clay": dict(base=(0.80, 0.74, 0.68), amb=0.30, key=0.72, fill=0.16, spec=(14.0, 0.06), ao=0.55, rim=0.04),
    "Plastic": dict(base=(0.62, 0.66, 0.72), amb=0.24, key=0.74, fill=0.18, spec=(70.0, 0.55), ao=0.45, rim=0.22),
    "Black plastic": dict(base=(0.14, 0.14, 0.15), amb=0.40, key=1.10, fill=0.50, spec=(40.0, 0.60), ao=0.40, rim=0.45),
    "Ceramic": dict(base=(0.93, 0.92, 0.89), amb=0.34, key=0.60, fill=0.14, spec=(180.0, 0.90), ao=0.35, rim=0.30),
    "Metal": dict(base=(0.80, 0.81, 0.84), amb=0.0, key=0.0, fill=0.0, spec=(120.0, 0.50), ao=0.30, rim=0.0, env=True),
}


def _matcap_image(name, size=256):
    """MATCAPS ayarindan size x size RGB doku (bytes; ilk satir kurenin alti, Coin dokusu gibi)."""
    c = MATCAPS[name]
    y, x = np.mgrid[0:size, 0:size]
    x, y = (x + 0.5) / size * 2 - 1, (y + 0.5) / size * 2 - 1
    r = np.maximum(np.hypot(x, y), 1.0)                    # daire disi kenar rengini alir (doku suzmesinde dikis yok)
    x, y = x / r, y / r
    n = np.dstack((x, y, np.sqrt(np.maximum(1 - x * x - y * y, 0.0))))
    light = lambda *v: np.array(v) / np.linalg.norm(v)
    key, fill = light(-0.45, 0.65, 0.65), light(0.75, 0.05, 0.45)
    half = (key + (0, 0, 1)) / np.linalg.norm(key + (0, 0, 1))
    base = np.array(c["base"])
    if c.get("env"):                                       # yansima yonune gore stüdyo: tavan, ufuk bandi, koyu zemin
        ry = 2 * n[..., 2] * n[..., 1]
        sky = 0.45 + 0.50 * np.clip(ry, 0, 1) ** 0.7
        ground = 0.10 + 0.25 * np.clip(-ry - 0.15, 0, 1)
        t = np.clip((ry + 0.04) / 0.08, 0, 1)
        env = ground + t * t * (3 - 2 * t) * (sky - ground) + 0.30 * np.exp(-(ry / 0.05) ** 2)
        rgb = env[..., None] * base
    else:
        dif = c["amb"] + c["key"] * np.maximum(n @ key, 0) + c["fill"] * np.maximum(n @ fill, 0)
        rgb = dif[..., None] * base
    p, k = c["spec"]
    rgb = rgb + (k * np.maximum(n @ half, 0) ** p)[..., None]
    rgb = rgb + (c["rim"] * (1 - n[..., 2]) ** 4 * np.clip(y * 0.5 + 0.5, 0, 1))[..., None]
    t = np.clip(n[..., 2] / 0.55, 0, 1)                    # kenara dogru karar: sahte AO
    rgb = rgb * (1 - c["ao"] * (1 - t * t * (3 - 2 * t)))[..., None]
    return (np.clip(rgb, 0, 1) * 255 + 0.5).astype(np.uint8).tobytes()


def _matcap_shader():
    """Dokuyu piksel basina okur: kure ortam haritasi koseler arasinda dogrusal gittiginden keskin dokuda ucgenler
    basamak basamak gorunuyordu. Yedek: shader derlenmezse ortam haritasi calismaya devam eder."""
    vs = coin.SoVertexShader()
    vs.sourceType = coin.SoShaderObject.GLSL_PROGRAM
    vs.sourceProgram = ("varying vec3 vN; void main(){ vN = gl_NormalMatrix * gl_Normal; "
                        "gl_ClipVertex = gl_ModelViewMatrix * gl_Vertex; gl_Position = ftransform(); }")
    fs = coin.SoFragmentShader()
    fs.sourceType = coin.SoShaderObject.GLSL_PROGRAM
    fs.sourceProgram = ("uniform sampler2D tex; varying vec3 vN; void main(){ vec3 n = normalize(vN); "
                        "if (!gl_FrontFacing) n = -n; "
                        "gl_FragColor = vec4(texture2D(tex, n.xy * 0.5 + 0.5).rgb, gl_FrontMaterial.diffuse.a); }")
    prog = coin.SoShaderProgram()
    prog.shaderObject.set1Value(0, vs)
    prog.shaderObject.set1Value(1, fs)
    return prog


# Ortam yansitma: equirectangular (2:1) ortam resmi piksel basina yansima yonunden okunur. Radiance .hdr ya da
# JPG / PNG; HDR 8 bit dokuya tonemap ile iner (Coin dokusu 8 bit). Ortam DUNYAYA bagli, Z yukari: model donunce
# yansima kayar.
ENV_FILTER = "Environment image (*.hdr *.jpg *.jpeg *.png)"
ENV_RADIUS = 1.0        # ortam kuresinin yaricapi / modelin kosegeni (kucuk: yansima konumla daha cok degisir)


def _read_hdr(path):
    """Radiance .hdr (RGBE, duz ya da RLE satirlar): H x W x 3 kayan nokta, ilk satir resmin ustu."""
    with open(path, "rb") as fh:
        data = fh.read()
    if not data.startswith((b"#?", b"RADIANCE")):
        raise ValueError("not a Radiance HDR file")
    end = data.find(b"\n\n")
    if end < 0:
        raise ValueError("broken HDR header")
    i = data.index(b"\n", end + 2)
    parts = data[end + 2:i].split()
    if len(parts) != 4 or parts[0] != b"-Y" or parts[2] != b"+X":
        raise ValueError("unsupported HDR orientation")
    H, W = int(parts[1]), int(parts[3])
    i += 1
    out = np.zeros((H, W, 4), np.uint8)
    if not 8 <= W <= 32767 or data[i] != 2 or data[i + 1] != 2 or data[i + 2] & 0x80:   # eski tip: duz piksel
        out[:] = np.frombuffer(data, np.uint8, H * W * 4, i).reshape(H, W, 4)
    else:
        for y in range(H):
            if data[i] != 2 or data[i + 1] != 2 or (data[i + 2] << 8 | data[i + 3]) != W:
                raise ValueError("broken HDR scanline")
            i += 4
            for c in range(4):
                x = 0
                while x < W:
                    n = data[i]
                    i += 1
                    if n > 128:                            # tekrar
                        n -= 128
                        out[y, x:x + n, c] = data[i]
                        i += 1
                    else:                                  # duz
                        out[y, x:x + n, c] = np.frombuffer(data, np.uint8, n, i)
                        i += n
                    x += n
    e = out[..., 3].astype(float)
    return out[..., :3] * np.where(e > 0, np.ldexp(1.0, (e - 136).astype(int)), 0.0)[..., None]


def _env_read(path, maxw=2048):
    """Ortam resmi -> (dogrusal RGB, H x W x 3, ilk satir resmin ustu; HDR mi). Hata: ValueError."""
    if path.lower().endswith(".hdr"):
        px, hdr = _read_hdr(path), True
    else:
        img = QtGui.QImage(path)
        if img.isNull():
            raise ValueError("cannot read the image")
        img = img.convertToFormat(QtGui.QImage.Format_RGB888)
        px = np.frombuffer(img.constBits(), np.uint8, img.sizeInBytes()).reshape(img.height(), img.bytesPerLine())
        px = (px[:, :img.width() * 3].reshape(img.height(), img.width(), 3) / 255.0) ** 2.2
        hdr = False
    k = max(int(math.ceil(px.shape[1] / maxw)), 1)         # buyuk resim: bellek ve doku siniri
    return px[::k, ::k], hdr


def _env_blur(px, sigma):
    """Panoramayi sigma radyan bulaniklastirir: FFT; yatayda sarmal ve kutuplara dogru genisleyen, dikeyde sabit."""
    if sigma <= 0:
        return px
    H, W = px.shape[:2]
    lat = (0.5 - (np.arange(H) + 0.5) / H) * math.pi
    sx = sigma * W / (2 * math.pi) / np.maximum(np.cos(lat), 0.02)
    f = np.fft.rfftfreq(W)
    out = np.fft.irfft(np.fft.rfft(px, axis=1) * np.exp(-2 * (math.pi * f[None, :] * sx[:, None]) ** 2)[..., None],
                       n=W, axis=1)
    pad = np.concatenate([out[::-1], out, out[::-1]])      # kutuplarda yansitma
    f = np.fft.rfftfreq(3 * H)
    out = np.fft.irfft(np.fft.rfft(pad, axis=0) * np.exp(-2 * (math.pi * f * sigma * H / math.pi) ** 2)[:, None, None],
                       n=3 * H, axis=0)[H:2 * H]
    return np.maximum(out, 0.0)


def _env_atlas(px, hdr, rough):
    """Tek dokuda iki katman (genislik, 2 x yukseklik, RGB bytes; ilk satir altta): ust yari yansima (rough kadar bulanik),
    alt yari yuzeyin rengini aydinlatan cok bulanik ortam."""
    if rough > 0.02 and px.shape[1] > 1024:                # bulanik yansimada yuksek cozunurluk bosa gider
        px = px[::2, ::2]
    H, W = px.shape[:2]
    k = max(W // 256, 1)
    low = _env_blur(px[::k, ::k], 0.55)
    dif = low[np.arange(H) * low.shape[0] // H][:, np.arange(W) * low.shape[1] // W]
    ref = _env_blur(px, rough * rough * 0.9)               # algi icin karesel: dusuk degerlerde ince ayar
    if hdr:                                                # ortalama parlaklik orta gri civari, sonra tonemap
        lum = px @ np.array([0.2126, 0.7152, 0.0722])
        g = 0.5 / max(float(np.exp(np.log(lum + 1e-4).mean())), 1e-4)
        ref, dif = (lambda a: g * a / (1 + g * a))(ref), (lambda a: g * a / (1 + g * a))(dif)
    data = (np.clip(np.concatenate([dif[::-1], ref[::-1]]), 0, 1) ** (1 / 2.2) * 255 + 0.5).astype(np.uint8)
    return W, 2 * H, data.tobytes()


# gl_ClipVertex: yoksa kesit duzlemleri bu golgelendiriciyle cizilen yuzeyi kesmez
ENV_VERT = ("varying vec3 vN; varying vec3 vP; void main(){ vN = gl_Normal; vP = gl_Vertex.xyz; "
            "gl_ClipVertex = gl_ModelViewMatrix * gl_Vertex; gl_Position = ftransform(); }")
ENV_FRAG = """uniform sampler2D tex; uniform vec3 px; uniform vec3 py; uniform vec3 pz; uniform vec3 ctr; uniform float rad;
uniform float rough; uniform float metal; uniform float texel; uniform vec3 base;
varying vec3 vN; varying vec3 vP;
vec2 pano(vec3 d){
  vec3 w = normalize(vec3(dot(px, d), dot(py, d), dot(pz, d)));
  return vec2(0.5 - atan(w.y, w.x) / 6.2831853, asin(clamp(w.z, -1.0, 1.0)) / 3.1415927 + 0.5);
}
void main(){
  mat3 M = mat3(gl_ModelViewMatrix[0].xyz, gl_ModelViewMatrix[1].xyz, gl_ModelViewMatrix[2].xyz);
  vec3 I = gl_ProjectionMatrix[3][3] > 0.5 ? vec3(0.0, 0.0, -1.0) : (gl_ModelViewMatrix * vec4(vP, 1.0)).xyz;
  I = normalize(I * M);
  vec3 n = normalize(vN);
  if (dot(n, I) > 0.0) n = -n;
  vec3 r = reflect(I, n);
  vec3 o = vP - ctr;                                     // resim modeli saran kureye sarili: isin kureye carptigi yon
  float b = dot(o, r);
  r = o + (sqrt(max(b * b - dot(o, o) + rad * rad, 0.0)) - b) * r;
  vec2 uv = pano(r), ud = pano(n);
  vec3 refl = texture2D(tex, vec2(uv.x, clamp(uv.y * 0.5 + 0.5, 0.5 + texel, 1.0 - texel))).rgb;
  vec3 dif = texture2D(tex, vec2(ud.x, clamp(ud.y * 0.5, texel, 0.5 - texel))).rgb;
  dif = mix(vec3(dot(dif, vec3(0.2126, 0.7152, 0.0722))), dif, 0.4);   // yuzey rengi ortamin rengine fazla boyanmasin
  float f0 = mix(0.08, 1.0, metal);
  float f = f0 + (1.0 - f0) * pow(1.0 - max(dot(-I, n), 0.0), 5.0) * (1.0 - 0.7 * rough);
  gl_FragColor = vec4((1.0 - metal) * (1.0 - f) * base * dif * 1.5 + f * refl, gl_FrontMaterial.diffuse.a);
}"""

def _env_nodes():
    """(grup, doku, [px, py, pz], ctr, rad, rough, metal, texel, base): shader nesne normalini kullanir, Placement
    donmesi parametreden; ctr, rad: resmin sarili oldugu kure (yerel); rough, metal: gorunum; texel: atlasta katman
    siniri; base: yuzey rengi, ortamin bulanik kopyasiyla aydinlanir; ustune Fresnel yansima gelir. Sonsuz uzak ortamda duz yuzeydeki paralel isinlar resmin tek noktasina bakar ve
    yuzey cizgi cizgi gorunurdu; sonlu kurede yuzeyin konumu da ise girer."""
    g = coin.SoGroup()
    q = coin.SoComplexity()
    q.textureQuality = 0.3                                 # mipmap yok: atan dikisinde cizgi cikmasin
    g.addChild(q)
    tex = coin.SoTexture2()
    tex.model = coin.SoTexture2.REPLACE
    g.addChild(tex)
    vs, fs = coin.SoVertexShader(), coin.SoFragmentShader()
    vs.sourceType = fs.sourceType = coin.SoShaderObject.GLSL_PROGRAM
    vs.sourceProgram, fs.sourceProgram = ENV_VERT, ENV_FRAG
    prog = coin.SoShaderProgram()
    prog.shaderObject.set1Value(0, vs)
    prog.shaderObject.set1Value(1, fs)
    rows = []
    for name in ("px", "py", "pz"):
        p = coin.SoShaderParameter3f()
        p.name = name
        rows.append(p)
    ctr, base = coin.SoShaderParameter3f(), coin.SoShaderParameter3f()
    ctr.name, base.name = "ctr", "base"
    rad, rough, metal, texel = (coin.SoShaderParameter1f() for _ in range(4))
    rad.name, rough.name, metal.name, texel.name = "rad", "rough", "metal", "texel"
    rad.value = 1000.0
    for i, p in enumerate(rows + [ctr, rad, rough, metal, texel, base]):
        fs.parameter.set1Value(i, p)
    g.addChild(prog)
    return g, tex, rows, ctr, rad, rough, metal, texel, base


def _zebra_button(on):
    """FreeCAD arac cubugundaki Zebra dugmesini (SubD_Zebra) oturumun durumuna getir."""
    cmd = Gui.Command.get("SubD_Zebra")
    acts = cmd.getAction() if cmd is not None else []
    if acts and acts[0].isChecked() != on:
        acts[0].blockSignals(True)
        acts[0].setChecked(on)
        acts[0].blockSignals(False)


def _later(fn):
    QtCore.QTimer.singleShot(0, fn)


def _gain():
    fine = QtWidgets.QApplication.keyboardModifiers() & QtCore.Qt.ShiftModifier
    return FINE_GAIN if fine else DRAG_GAIN


def _left_held():
    return bool(QtWidgets.QApplication.mouseButtons() & QtCore.Qt.LeftButton)


def _ctrl_held():
    return bool(QtWidgets.QApplication.keyboardModifiers() & QtCore.Qt.ControlModifier)


def _snap(x, step):
    return round(x / step) * step if step > 0 else x


def _vv(v):
    W, H = v.getSize()
    return v.getCameraNode().getViewVolume(float(W) / H), W, H


class _Handle(object):
    """Tek bir tutamak: renk/kalinlik dugumleri (vurgu icin) + cizim."""

    def __init__(self, rgb, width):
        self.rgb = rgb
        self.width = width
        self.sep = coin.SoSeparator()
        self.col = coin.SoBaseColor()
        self.col.rgb = rgb
        self.style = coin.SoDrawStyle()
        self.style.lineWidth = width
        self.sep.addChild(self.col)
        self.sep.addChild(self.style)

    def hot(self, on):
        self.col.rgb = HOT if on else self.rgb
        self.style.lineWidth = self.width + (1.5 if on else 0.0)

    def lines(self, pts, idx):
        co = coin.SoCoordinate3()
        co.point.setValues(0, len(pts), pts)
        ls = coin.SoIndexedLineSet()
        ls.coordIndex.setValues(0, len(idx), idx)
        self.sep.addChild(co)
        self.sep.addChild(ls)
        return co


class _FillHandle(_Handle):
    """Dolgulu tutamak (duzlem karesi, merkez kup): cizgi + yari saydam dolgu; soluk = tutulamaz."""

    def __init__(self, rgb, width, alpha):
        _Handle.__init__(self, rgb, width)
        self.alpha = alpha
        self.mat = coin.SoMaterial()
        self.mat.diffuseColor = rgb
        self.mat.emissiveColor = rgb
        self.mat.transparency = alpha
        self.fill = coin.SoSeparator()
        self.fill.addChild(self.mat)
        self.sep.addChild(self.fill)
        self._dim = False

    def hot(self, on):
        _Handle.hot(self, on and not self._dim)
        c = HOT if (on and not self._dim) else self.rgb
        self.mat.diffuseColor = c
        self.mat.emissiveColor = c

    def dim(self, off):
        if off == self._dim:
            return
        self._dim = off
        self.mat.transparency = 0.93 if off else self.alpha
        self.col.rgb = tuple(0.35 + 0.3 * c for c in self.rgb) if off else self.rgb


def _plane_corners(i):
    a, b = AXES[(i + 1) % 3], AXES[(i + 2) % 3]
    return [a * u + b * w for u, w in ((PLANE_A, PLANE_A), (PLANE_B, PLANE_A), (PLANE_B, PLANE_B), (PLANE_A, PLANE_B))]


def _plane_square(i):
    h = _FillHandle(GZ_RGB[i], 1.5, 0.6)
    pts = [(p.x, p.y, p.z) for p in _plane_corners(i)]
    h.lines(pts, [0, 1, 2, 3, 0, -1])
    co = coin.SoCoordinate3()
    co.point.setValues(0, 4, pts)
    fs = coin.SoFaceSet()
    fs.numVertices.setValue(4)
    sh = coin.SoShapeHints()
    sh.vertexOrdering = coin.SoShapeHints.UNKNOWN_ORDERING
    for n in (sh, co, fs):
        h.fill.addChild(n)
    return h.sep, h


def _free_cube():
    h = _FillHandle((0.92, 0.92, 0.95), 1.5, 0.15)
    r = FREE_R
    c = [(x, y, z) for x in (-r, r) for y in (-r, r) for z in (-r, r)]
    h.lines(c, [0, 1, 3, 2, 0, -1, 4, 5, 7, 6, 4, -1, 0, 4, -1, 1, 5, -1, 2, 6, -1, 3, 7, -1])
    cu = coin.SoCube()
    cu.width = cu.height = cu.depth = 2 * r
    h.fill.addChild(cu)
    return h.sep, h


def _arrow(i):
    h = _Handle(GZ_RGB[i], 3.0)
    s = coin.SoSeparator()
    r = coin.SoRotation()
    r.rotation.setValue(coin.SbRotation(coin.SbVec3f(1, 0, 0), coin.SbVec3f(*AXES[i])))
    s.addChild(r)
    s.addChild(h.sep)
    h.lines([(0, 0, 0), (1, 0, 0)], [0, 1, -1])
    t = coin.SoTransform()
    t.translation = (1.05, 0, 0)
    t.rotation.setValue(coin.SbVec3f(0, 0, 1), -math.pi / 2)
    cn = coin.SoCone()
    cn.bottomRadius = 0.06
    cn.height = 0.18
    h.sep.addChild(t)
    h.sep.addChild(cn)
    return s, h


def _stretch_handle(i):
    """Stretch tutamagi: eksen cizgisi + ucunda kucuk kup (LightWave Stretch)."""
    h = _Handle(GZ_RGB[i], 3.0)
    s = coin.SoSeparator()
    r = coin.SoRotation()
    r.rotation.setValue(coin.SbRotation(coin.SbVec3f(1, 0, 0), coin.SbVec3f(*AXES[i])))
    s.addChild(r)
    s.addChild(h.sep)
    h.lines([(0, 0, 0), (1, 0, 0)], [0, 1, -1])
    t = coin.SoTranslation()
    t.translation = (1.04, 0, 0)
    cu = coin.SoCube()
    cu.width = cu.height = cu.depth = 0.12
    h.sep.addChild(t)
    h.sep.addChild(cu)
    return s, h


def _size_cube():
    """Size tutamagi: merkezde orantili olcekleme kupu."""
    h = _FillHandle((0.95, 0.80, 0.25), 1.5, 0.35)
    r = 0.11
    c = [(x, y, z) for x in (-r, r) for y in (-r, r) for z in (-r, r)]
    h.lines(c, [0, 1, 3, 2, 0, -1, 4, 5, 7, 6, 4, -1, 0, 4, -1, 1, 5, -1, 2, 6, -1, 3, 7, -1])
    cu = coin.SoCube()
    cu.width = cu.height = cu.depth = 2 * r
    h.fill.addChild(cu)
    return h.sep, h


RING_N = 64


def _ring_pts(i):
    """i ekseni etrafindaki birim cemberin noktalari (eksene dik duzlemde)."""
    a = AXES[(i + 1) % 3]
    b = AXES[(i + 2) % 3]
    return [a * math.cos(2 * math.pi * k / RING_N) + b * math.sin(2 * math.pi * k / RING_N) for k in range(RING_N)]


def _ring(i):
    h = _Handle(GZ_RGB[i], 3.5)
    h.lines([(p.x, p.y, p.z) for p in _ring_pts(i)], list(range(RING_N)) + [0, -1])
    return h.sep, h


def _rotation_matrix(axis, ang):
    """Birim eksen etrafinda ang radyanlik donus matrisi (Rodrigues)."""
    x, y, z = axis
    c, s, t = math.cos(ang), math.sin(ang), 1 - math.cos(ang)
    return np.array([[t * x * x + c, t * x * y - s * z, t * x * z + s * y],
                     [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
                     [t * x * z - s * y, t * y * z + s * x, t * z * z + c]])


def _ray(v, pos):
    """Imlec pikselinden (sol alt kokenli) dunya isini: (baslangic, birim yon)."""
    p = v.getPoint((pos[0], pos[1]))
    cn = v.getCameraNode()
    if cn.isOfType(coin.SoPerspectiveCamera.getClassTypeId()):
        o = App.Vector(*cn.position.getValue().getValue())
        d = p - o
    else:
        d = App.Vector(*v.getViewDirection())
        o = p - d * 1e7
    d.normalize()
    return o, d


def _to_screen(vv, W, H, p):
    sx, sy, _ = vv.projectToScreen(coin.SbVec3f(p.x, p.y, p.z)).getValue()
    return sx * W, sy * H


def _seg_dist(px, py, a, b):
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = 0.0 if L2 < 1e-9 else max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / L2))
    return math.hypot(px - (ax + t * dx), py - (ay + t * dy))


def _view_dir(v, p):
    cn = v.getCameraNode()
    if cn.isOfType(coin.SoPerspectiveCamera.getClassTypeId()):
        d = p - App.Vector(*cn.position.getValue().getValue())
    else:
        d = App.Vector(*v.getViewDirection())
    if d.Length > 1e-12:
        d.normalize()
    return d


def _in_poly(x, y, poly):
    inside = False
    n = len(poly)
    for k in range(n):
        (x1, y1), (x2, y2) = poly[k], poly[(k + 1) % n]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside


def _line_param(o, d, p0, a):
    w = o - p0
    b = d.dot(a)
    den = 1.0 - b * b
    if abs(den) < 1e-6:
        return None
    return (w.dot(a) - b * w.dot(d)) / den


def _plane_hit(o, d, c, n):
    den = d.dot(n)
    if abs(den) < 1e-6:
        return None
    t = (c - o).dot(n) / den
    return o + d * t


class MoveGizmo(object):
    """Secimin merkezinde gizmo. Tasi: X/Y/Z oklari, XY/YZ/XZ duzlem kareleri, merkezde serbest tasima kupu.
    Dondur: X/Y/Z eksenleri etrafinda halkalar. Yalnizca cizimdir; fare olaylarini oturum yakalar ve
    pick() ile tutamagi sorar."""

    def __init__(self, view):
        self.view = view
        self.mode = "move"   # "move" | "rotate"
        self.center = None
        self.size = 1.0
        self.px = min(max(_param().GetFloat("GizmoSize", GIZMO_PX), 30.0), 400.0)   # ok boyu (px)
        self.hover = None
        self.handles = {}
        self.root = coin.SoSwitch()
        ann = coin.SoAnnotation()          # derinlik testi olmadan: nesnenin icinde kalsa da gorunur
        self.pos = coin.SoTranslation()
        self.scale = coin.SoScale()
        ann.addChild(self.pos)
        ann.addChild(self.scale)
        self.modes = coin.SoSwitch()
        ann.addChild(self.modes)
        g = coin.SoSeparator()
        for i in range(3):
            s, h = _arrow(i)
            g.addChild(s)
            self.handles[("move", i)] = h
        for i in range(3):
            s, h = _plane_square(i)
            g.addChild(s)
            self.handles[("plane", i)] = h
        s, h = _free_cube()
        g.addChild(s)
        self.handles[("free", 0)] = h
        self.modes.addChild(g)
        g = coin.SoSeparator()
        for i in range(3):
            s, h = _ring(i)
            g.addChild(s)
            self.handles[("rotate", i)] = h
        self.modes.addChild(g)
        g = coin.SoSeparator()                         # Stretch: uclari kuplu eksenler; Size: merkez kup
        for i in range(3):
            s, h = _stretch_handle(i)
            g.addChild(s)
            self.handles[("stretch", i)] = h
        s, h = _size_cube()
        g.addChild(s)
        self.handles[("size", 0)] = h
        self.modes.addChild(g)
        self.root.addChild(ann)
        self.root.whichChild = -1
        self.set_mode("move")

    def set_mode(self, mode):
        self.mode = mode
        self.modes.whichChild = {"move": 0, "rotate": 1, "scale": 2}[mode]

    def place(self, center):
        self.center = center
        if center is None:
            self.root.whichChild = -1
            return
        self.root.whichChild = 0
        self.pos.translation = (center.x, center.y, center.z)
        self.update_scale()

    def update_scale(self):
        if self.center is None:
            return
        vv, W, H = _vv(self.view)
        p = self.center
        size = vv.getWorldToScreenScale(coin.SbVec3f(p.x, p.y, p.z), self.px / float(H))
        if abs(size - self.size) > 1e-9 * max(size, 1.0):
            self.size = size
            self.scale.scaleFactor.setValue(size, size, size)
        vd = _view_dir(self.view, p)
        for i in range(3):
            self.handles[("plane", i)].dim(abs(vd.dot(AXES[i])) < PLANE_EDGE)

    def pick(self, pos):
        """Imlecin PICK_PX yakinindaki tutamak: (tur, eksen) ya da None."""
        if self.center is None:
            return None
        vv, W, H = _vv(self.view)
        x, y = pos
        if self.mode == "rotate":
            best, bd = None, PICK_PX
            for i in range(3):
                pts = [_to_screen(vv, W, H, self.center + q * self.size) for q in _ring_pts(i)[::2]]
                for k in range(len(pts)):
                    dd = _seg_dist(x, y, pts[k], pts[(k + 1) % len(pts)])
                    if dd < bd:
                        best, bd = ("rotate", i), dd
            return best
        c = _to_screen(vv, W, H, self.center)
        k = max(self.px / GIZMO_PX, 0.6)                   # merkez kupun yakalama yaricapi gizmoyla buyur
        if self.mode == "scale":
            if math.hypot(x - c[0], y - c[1]) < 12.0 * k:
                return ("size", 0)
            best, bd = None, PICK_PX
            for i in range(3):
                e = _to_screen(vv, W, H, self.center + AXES[i] * (self.size * 1.1))
                dd = _seg_dist(x, y, c, e)
                if dd < bd:
                    best, bd = ("stretch", i), dd
            return best
        if math.hypot(x - c[0], y - c[1]) < 10.0 * k:  # oncelik: merkez kup, sonra kareler, sonra oklar
            return ("free", 0)
        for i in range(3):
            if self.handles[("plane", i)]._dim:
                continue
            poly = [_to_screen(vv, W, H, self.center + q * self.size) for q in _plane_corners(i)]
            if _in_poly(x, y, poly):
                return ("plane", i)
        best, bd = None, PICK_PX
        for i in range(3):
            e = _to_screen(vv, W, H, self.center + AXES[i] * (self.size * 1.2))
            dd = _seg_dist(x, y, c, e)
            if dd < bd:
                best, bd = ("move", i), dd
        return best

    def set_hover(self, key):
        if key == self.hover:
            return
        if self.hover in self.handles:
            self.handles[self.hover].hot(False)
        if key in self.handles:
            self.handles[key].hot(True)
        self.hover = key


HELP = ("1 / 2 / 3: point / edge / polygon mode\n"
        "Y: X symmetry on/off\n"
        "R: rotate gizmo  |  H: stretch / size gizmo\n"
        "C: crease selected edges\n"
        "Click: select  |  Shift+click: add\n"
        "Drag: box select (Shift add, Ctrl remove)\n"
        "L: loop from the selection (or under the cursor)  |  X: ghost\n"
        "E: extrude  |  K: knife (loop cut)  |  J: knife (line / plane)  |  Shift+V: mirror (line)\n"
        "Knife / loop stay on (yellow) until K / J / Esc / right click\n"
        "B: bevel (right = inset, up = shift)  |  I: inset\n"
        "Edge mode B: edge bevel (right = width)\n"
        "Edge mode E: extrude open edges  |  extrude: wheel = segments, X / Y / Z = axis\n"
        "O: object / edit mode  |  D: duplicate object\n"
        "T: move gizmo  |  G: free move (mouse)  |  Tab / S: SubPatch\n"
        "A: select all  |  Del: delete  |  P: make polygon\n"
        "Ctrl+C / X / V: copy / cut / paste polygons\n"
        "Ctrl+Z: undo  |  Ctrl+Y: redo  |  F: save to document\n"
        "Ctrl + gizmo / T: stepped (mm, degrees)\n"
        "Esc / right click: cancel")


class SubDSession:
    def __init__(self, obj=None, object_mode=False):
        self.view = Gui.ActiveDocument.ActiveView
        self.doc_name = Gui.ActiveDocument.Document.Name   # bu belge kapaninca oturum kapanir (_DocWatcher)
        self.obj = obj       # duzenlenen (Object modunda: son etkin) SubD nesnesi; yeni kafes ilk kayitta olusur
        self.dirty = False
        self.info = ""       # surukleme sirasinda gosterilen deger (delta / aci)
        self.cage = cage_from_obj(obj) if obj is not None else \
            Cage.cube([_last_box[n] for n in self.BOX_SIZE])   # nesnesiz Object modu: gorunmeyen yer tutucu
        self.hidden = []     # duzenlerken gizlenen nesneler (SubD ve ondan uretilen katilar); kapatinca geri
        self._groups_backfill_color()
        self._snapshot()
        if not object_mode:
            self._hide_edited()
        self.mode = "object" if object_mode else "face"   # "object" | "point" | "edge" | "face"
        self.edit_mode = "face"                            # Object modundan Edit'e donulecek secim modu
        self.sel = set()     # yuz indeksleri
        self.sel_v = set()   # nokta indeksleri
        self.last_v = None   # en son tiklanan secili nokta (Weld hedefi)
        self.sel_e = set()   # frozenset({a, b}) kenarlar
        self.sel_order = {}  # mod -> tiklama sirasiyla secilenler (Align: ilk tiklanan sabit kalir)
        self.sel_order_key = None   # sel_order'in gecerli oldugu (nokta, yuz) sayisi: numaralar degisince silinir
        self.sym = False     # X=0 duzlemine gore simetri
        self.ghost = False   # hayalet modu: arkadakiler gorunur / secilebilir
        self.snap = False    # serbest tasimada Part / PartDesign kose ve kenarlarina yakalama
        self.zebra = False   # zebra seritleri (yuzey surekliligi analizi)
        self.matcap = _param().GetString("Matcap", "")     # yuzey gorunumu (MATCAPS adi; bos = kapali)
        self.env_on = _param().GetBool("EnvOn", False)     # ortam yansitma (acikken matcap'in onune gecer)
        self.env_path = _param().GetString("EnvPath", "")  # ortam resmi
        self.env_loaded = None                             # dokuya yuklu resim
        self.env_px = None                                 # okunan resim (dogrusal, HDR mi); Roughness degisince gerekir
        self.env_rough = _param().GetFloat("EnvRoughness", 0.25)   # 0 ayna, 1 tamamen dagitik
        self.env_metal = _param().GetFloat("EnvMetal", 0.0)        # 0 plastik (renk + hafif yansima), 1 metal
        self.cage_color = QtGui.QColor(_param().GetString("CageColor", "#8c99b3"))  # Subpatch kapaliyken kafes yuzleri
        if not self.cage_color.isValid():
            self.cage_color = QtGui.QColor("#8c99b3")
        self.env_color = QtGui.QColor(_param().GetString("EnvColor", "#b8bfcc"))   # plastik yuzey rengi
        if not self.env_color.isValid():
            self.env_color = QtGui.QColor("#b8bfcc")
        self._snapdata = None
        self.smooth = obj is not None
        self.levels = obj.Levels if obj is not None else 3       # son (belge / kati) seviyesi
        self.rules = obj_rules(obj) if obj is not None else RULES
        self.edit_levels = _param().GetInt("EditLevels", 3)      # duzenlerken gorunen seviye (hiz icin dusuk)
        self.undo = []
        self.redo = []
        self.histories = {}  # nesne adi -> (undo, redo, cage_key): Object moduna gidip donunce adimlar kaybolmasin
        self.splits = []     # Knife "Split into two objects": (geri alma derinligi, yeni nesnenin adi)
        self.grab = None     # klavyeyle tasima (G, extrude)
        self.drag = None     # gizmo tutamagiyla tasima
        self.armed = None    # kalici arac ("knife" / "line" / "loop" / "magnet"): kapatilana kadar her tiklamada uygulanir
        self.magnet_radius = None    # Magnet yaricapi (mm); ilk acilista kafesin boyuna gore
        self.magnet_mode = 0         # Magnet azalma egrisi (FALLOFFS sirasi)
        self.magnet_radial = False   # Magnet etki alani: False kure, True bakis yonunde boydan boya silindir
        self.magnet_size = None      # Magnet: sag tusla yaricap ayari surerken
        self._magnet_snap = None     # Magnet: (kamera + kafes anahtari, _snap_prepare hedefleri, kafes dahil hedefler)
        self.op = None       # surmekte olan islem (panelde canli degerler)
        self.ebevel_segs = 1  # kenar bevel segment sayisi (son kullanilan)
        self.extrude_segs = 1        # extrude dilim sayisi (son kullanilan)
        self.bevel_each = False      # Bevel "Individual": her yuz ayri (son kullanilan)
        self.extrude_each = False    # extrude: her yuz kendi normalinde (son kullanilan)
        self.rail = None             # Rail bridge duzenleniyor: dict(i=kilavuz, sel=ara nokta, span=parca)
        self.rail_segs = 4           # Rail bridge baslangic ara halka sayisi (son kullanilan)
        self.rail_falloff = False    # Rail bridge: secili nokta tasinirken kilavuz kablo gibi gelsin (rail_cable)
        self.rail_elastic = 0.1      # ... esneklik: 0 zincir (boylar korunur), 1 lastik; 0.1 kablo
        self._last_click = None      # cift tiklama (kilavuzu yeniden duzenleme)
        self.box_vals = dict(_last_box)  # Box paneli: son kullanilan degerler (modelleyici kapansa da)
        self.cyl_vals = dict(_last_cyl)  # Cylinder paneli (iki topoloji ortak)
        self.bridge_vals = {"Number of Cuts": 0.0, "Twist": 0.0, "Smoothness": 0.0, "Profile Factor": 0.0,
                            "Auto direction": 1.0}
        self.last_op = None  # biten son islem (panelden degerle yeniden uygulanabilir)
        self.box = None      # bos alana basildi: surukleme kutu secimi olur
        self.rubber = None   # secim kutusu cizimi (QRubberBand)
        self.line_overlay = None  # knife cizgisi cizimi
        self.dim_p0 = None        # Distance: ilk nokta (hedef); ikinci tiklamada olcu olur
        self.dim_tmp = None       # gecici olcu (dict a, b)
        self.dim_hot = None       # imlecin altindaki kalici olcu (Distance acikken tiklayinca silinir)
        self.dim_last = None      # panelde gosterilen son olcunun iki dunya noktasi
        self.dim_live_q = None    # ilk nokta varken canli ikinci uc (dunya)
        self.dim_ext = None       # baska nesnelerin yakalanacak noktalari (arac acilinca toplanir)
        self.dim_root = self.dim_group = None   # olculerin sahne dugumu (_dim_scene)
        self.dim_list = None      # Distances liste penceresi (_DimList)
        self.dim_list_was_open = False
        self._hd_key, self._hd_names = None, []   # grup tutamaklarinin olculeri (_handle_dim_items)
        self._grp_sets = {}                       # grup adi -> uye indeksleri (_group_faces)
        self.grp_active = self.grp_dismissed = None   # grup paneli: acik grup, kapatilan secim
        self.dim_sel = None       # listede secili olcu (kafes sirasi)
        self.dim_hits = []        # kalici olculerin x dugmeleri (tiklayinca silinir)
        self._dim_imgs = {}       # olcu yazisi goruntuleri onbellegi
        self.dim_sig = None
        self._dim_cd = None
        self._smooth = None       # son SubPatch yuzeyi (noktalar, dortgenler): Distance yuzey yakalamasi
        self._line_cursor = None  # knife cizgisi yakalarken: gorunumun onceki imleci
        self.overlay = None
        self.suspended = False   # baska bir calisma tezgahinda: arac paneli / cubuklar / yakalama kapali (bkz. suspend)
        self.dock_was_open = False
        wb = _active_wb()
        self.wb_home = {"SubDModelerWorkbench"} | ({wb} if wb else set())   # modelleyicinin gorundugu tezgahlar
        self._pending = False
        self._born = None         # az once eklenen Box / Cylinder nesnesinin adi (bkz. _born_undo)
        self.mouse = (0, 0)
        self.gizmo = MoveGizmo(self.view)
        self._build_scene()
        self.set_placement(obj.Placement if obj is not None else App.Placement())
        self._build_dock()
        # pivy geri cagrilari: gizmo tutamagi tutulunca olay "islendi" sayilir, secime/gezinmeye gitmez
        self._cb_btn, self._cb_mov = self._on_button, self._on_motion
        self.view.addEventCallbackPivy(coin.SoMouseButtonEvent.getClassTypeId(), self._cb_btn)
        self.view.addEventCallbackPivy(coin.SoLocation2Event.getClassTypeId(), self._cb_mov)
        self.timer = QtCore.QTimer()      # yakinlastirinca gizmo ekranda ayni boyda kalsin
        self.timer.timeout.connect(self.gizmo.update_scale)
        self.timer.timeout.connect(self._ghost_tick)
        self.timer.start(50)
        self.autosave = QtCore.QTimer()   # degisiklikler belgeye duzenli yazilsin (FreeCAD kapanirsa kaybolmasin)
        self.autosave.timeout.connect(self._autosave)
        self.autosave.start(20000)
        self.wb_timer = QtCore.QTimer()   # baska calisma tezgahina gecildi mi (bkz. suspend)
        self.wb_timer.timeout.connect(self._wb_tick)
        self.wb_timer.start(300)
        self._mdi = Gui.getMainWindow().findChild(QtWidgets.QMdiArea)   # belge sekmesi degisince hemen karar ver
        self._mdi_cb = lambda *a: _later(lambda: _quiet(self._wb_tick))
        if self._mdi is not None:
            _quiet(lambda: self._mdi.subWindowActivated.connect(self._mdi_cb))
        self.dim_timer = QtCore.QTimer()  # gorunum donunce olculer ekranda yerinde kalsin
        self.dim_timer.timeout.connect(self._dim_tick)
        self.dim_timer.start(150)
        # Klavye Qt duzeyinde yakalanir: Coin'e giden olaylarin bazilarini (Tab, A, E...)
        # FreeCAD tuketiyor, Ctrl+Z gibi tuslar da FreeCAD kisayollarina gidiyor.
        gv = self.gv = _graphics_view(self.view)
        self.key_targets = [gv, gv.viewport()]
        self.key_filter = _KeyFilter(self)
        for w in self.key_targets:
            w.installEventFilter(self.key_filter)
        self.overlay = _Overlay(self)
        self.wheel_guard = _WheelGuard(gv.viewport())
        QtWidgets.QApplication.instance().installEventFilter(self.wheel_guard)
        self.history_guard = _HistoryGuard(self)
        QtWidgets.QApplication.instance().installEventFilter(self.history_guard)
        self.geo_switch.whichChild = -1 if object_mode else 0
        if object_mode and not any(is_subd(o) for o in App.getDocument(self.doc_name).Objects):
            self.info = T("SubD Modeler is ready: Create > Box or Cylinder to start")   # bos acilis (issue #7)
        self.refresh()
        if not object_mode:
            self.view.viewIsometric()
            self.view.fitAll()

    # --- Object / Edit modu ---
    # Object modunda belgedeki butun SubD nesneleri kendi halleriyle (Placement'lariyla) gorunur; secim FreeCAD'in
    # secimidir (agactan da secilebilir). Gizmo secili nesneleri tasir / dondurur (Placement) ya da olcekler
    # (kafes noktalari). Tab: secili (etkin) nesnenin kafesine girer / Object moduna doner.

    def _hide_edited(self):
        """Duzenlenen nesneyi ve ondan uretilen katilari gizle (yerine kendi onizlememiz gorunur)."""
        _tree_mark(self.obj)
        self.hidden = []
        if self.obj is not None:
            for o in [self.obj] + [o for o in self.obj.InList if "Samples" in o.PropertiesList]:
                if o.ViewObject.Visibility:
                    o.ViewObject.Visibility = False
                    self.hidden.append(o)

    def _show_edited(self):
        _tree_mark(None)
        if self.obj is not None and not self.hidden:
            self.obj.ViewObject.Visibility = True
        for o in self.hidden:
            try:
                o.ViewObject.Visibility = True
            except Exception:
                pass
        self.hidden = []

    def subd_objects(self):
        doc = App.ActiveDocument
        return [o for o in doc.Objects if is_subd(o)] if doc else []

    def selected_objects(self):
        """Object modunda secili SubD nesneleri (FreeCAD secimi; secili kati -> kaynagi)."""
        out = []
        for o in Gui.Selection.getSelection():
            if not is_subd(o) and "Samples" in getattr(o, "PropertiesList", []) and is_subd(o.Source):
                o = o.Source
            if is_subd(o) and o not in out:
                out.append(o)
        return out

    def set_object_selection(self, objs):
        Gui.Selection.clearSelection()
        for o in objs:
            Gui.Selection.addSelection(o.Document.Name, o.Name)

    def objects_center(self, objs=None):
        objs = self.selected_objects() if objs is None else objs
        cs = []
        for o in objs:
            bb = o.Mesh.BoundBox
            if bb.isValid():
                cs.append((bb.Center.x, bb.Center.y, bb.Center.z))
        return np.mean(cs, 0) if cs else None

    def enter_object_mode(self):
        """Kafesi nesneye yaz, nesneyi goster, Object moduna gec (nesne secili kalir)."""
        if self.mode == "object" or self.drag or self.grab:
            return
        if len(self.cage.V):                               # yuzsuz (yalniz kenar / nokta) kafes de yazilir
            self.save(commit=True)
        else:
            self._drop_empty()
        self._show_edited()
        self.drop_tool()
        self.dim_sel = None
        self.edit_mode = self.mode
        self.mode = "object"
        self.geo_switch.whichChild = -1
        self.op = self.last_op = None
        self.rail = None
        self._keep_history()
        self.undo, self.redo, self.splits = [], [], []
        self.set_object_selection([self.obj] if self.obj is not None else [])
        self.refresh()

    def edit_object(self, o=None):
        """Object modundan: verilen ya da secili (son) nesnenin kafesine gir."""
        if self.drag or self.grab:
            return
        if o is None:
            objs = self.selected_objects()
            o = objs[-1] if objs else self.obj
        if o is None:
            return
        self._born = None
        Gui.Selection.clearSelection()
        self.obj = o
        self.cage = cage_from_obj(o)
        self.levels = o.Levels
        self.rules = obj_rules(o)
        self._sync_surface_widgets()
        h = self.histories.pop(o.Name, None)
        self.undo, self.redo = h[:2] if h is not None and h[2] == cage_key(o) else ([], [])
        self.op, self.last_op = None, None
        self.splits = []
        self.rail = None
        self.sel, self.sel_v, self.sel_e = set(), set(), set()
        self.dirty = False
        self.set_placement(o.Placement)
        self._snapshot()
        self._hide_edited()
        self.mode = self.edit_mode
        self.geo_switch.whichChild = 0
        self.dim_sel = None
        self._dim_end()
        self.refresh()

    def toggle_object_mode(self):
        if self.mode == "object":
            self.edit_object()
        else:
            self.enter_object_mode()

    def add_box(self):
        """Box dugmesi: her zaman YENI bir SubD nesnesi (Edit modunda once mevcut kafes kaydedilir; eskiden
        Edit modunda mevcut kafesin yerine kutu konuyordu ve model kayboluyordu)."""
        if self.drag or self.grab:
            return
        if self.mode != "object":
            self.enter_object_mode()
        self.object_new_box()

    BOX_SIZE = ("Width", "Depth", "Height")                # X, Y, Z (mm)
    BOX_SEGS = ("Segments X", "Segments Y", "Segments Z")

    def object_new_box(self):
        """Yeni bir SubD kutusu ekle, hemen kafesine gir ve Box panelini ac (LightWave Box numeric): en / derinlik
        / yukseklik ve X / Y / Z segment sayilari; deger yazip Enter = kutu yeniden kurulur. Son degerler
        hatirlanir. (Object modunda yalniz puruzsuz yuzey gorunur; tek kutu kafesinin yuzeyi kuremsi
        oldugundan 'kutu istedim, kure geldi' gibi gorunuyordu.)"""
        vals = dict(self.box_vals)
        vals.update(zip(self.BOX_SIZE, (b.value() for b in self.box_size)))   # Options paneldeki olculer
        self._new_cage_object("SubD box", self._box_cage(vals))
        self.last_op = self.new_op("box", vals)
        self.op = None
        self.info = T("Box: size / segments in the panel")
        self.refresh()

    def _new_cage_object(self, title, cage):
        """Kafesten yeni bir SubD nesnesi kur ve hemen kafesine gir (Polygon modu); panel islemi icin geri
        alma durumu da kaydedilir (new_op'un 'before'i)."""
        doc = App.ActiveDocument
        doc.openTransaction(title)
        o = doc.addObject("Mesh::FeaturePython", "SubD")
        SubDObject(o)
        if App.GuiUp:
            SubDViewProvider(o.ViewObject)
        cage_to_obj(cage, o)
        set_surface(o, self.levels, self.rules)
        o.recompute()
        doc.commitTransaction()
        self.edit_mode = "face"
        self.edit_object(o)
        self.push()
        self._born = o.Name

    CYL_FIELDS = ("Diameter", "Height", "Segments", "Height Segments")

    def add_cylinder(self, topo=False):
        """Create > Cylinder / Topo cylinder: her zaman YENI bir SubD nesnesi, Z ekseni boyunca. Panelde cap /
        yukseklik / segment / yukseklik segmenti; deger yazip Enter = silindir yeniden kurulur. Cylinder:
        kapaklar merkeze ucgen yelpaze; Topo cylinder: kapaklar dortgen izgara (segment 4'un katina yuvarlanir)."""
        if self.drag or self.grab:
            return
        if self.mode != "object":
            self.enter_object_mode()
        vals = dict(self.cyl_vals)
        self._new_cage_object("SubD cylinder", self._cyl_cage(vals, topo))
        self.last_op = self.new_op("cylinder", vals, topo=topo)
        self.op = None
        self.info = T("Cylinder: size / segments in the panel")
        self.refresh()

    def _cyl_cage(self, v, topo):
        """Cylinder degerlerini sinirla (olcu > 0, segment 3-128, topo: 4'un kati, yukseklik segmenti 1-64)."""
        for n in ("Diameter", "Height"):
            v[n] = max(float(v[n]), 0.001)
        segs = min(max(int(round(v["Segments"])), 3), 128)
        if topo:
            segs = max(4 * int(round(segs / 4.0)), 4)
        v["Segments"] = float(segs)
        v["Height Segments"] = float(min(max(int(round(v["Height Segments"])), 1), 64))
        self.cyl_vals = dict(v)
        _last_cyl.update(v)
        return Cage.cylinder(v["Diameter"], v["Height"], segs, int(v["Height Segments"]), topo)

    def _box_cage(self, v):
        """Box degerlerini sinirla (olcu > 0, segment 1-64) ve kafesi kur; degerler sonraki kutular icin hatirlanir."""
        for n in self.BOX_SIZE:
            v[n] = max(float(v[n]), 0.001)
        for n in self.BOX_SEGS:
            v[n] = float(min(max(int(round(v[n])), 1), 64))
        self.box_vals = dict(v)
        _last_box.update(v)
        for b, n in zip(self.box_size, self.BOX_SIZE):     # Options paneldeki kutu olculeri de ayni olsun
            b.blockSignals(True)
            b.setValue(v[n])
            b.blockSignals(False)
        return Cage.box([v[n] for n in self.BOX_SIZE], [int(v[n]) for n in self.BOX_SEGS])

    def object_duplicate(self):
        """Secili nesnelerin kopyalari (ayni yerde); kopyalar secili olur, gizmoyla tasinabilir."""
        objs = self.selected_objects()
        if not objs:
            return
        doc = objs[0].Document
        doc.openTransaction("SubD duplicate")
        copies = []
        for o in objs:
            n = doc.addObject("Mesh::FeaturePython", "SubD")
            SubDObject(n)
            if App.GuiUp:
                SubDViewProvider(n.ViewObject)
            cage_to_obj(cage_from_obj(o), n)
            set_surface(n, o.Levels, obj_rules(o))
            n.Placement = o.Placement
            n.Label = o.Label
            n.recompute()
            copies.append(n)
        doc.commitTransaction()
        self.set_object_selection(copies)
        self.refresh()

    def object_mirror(self, axis=0):
        """Secili nesnelerin dunya X / Y / Z = 0 duzlemine gore aynalanmis kopyalari (kafes dunyada aynalanir)."""
        objs = self.selected_objects()
        if not objs:
            return
        flip = np.ones(3)
        flip[axis] = -1.0
        doc = objs[0].Document
        doc.openTransaction("SubD mirror")
        copies = []
        for o in objs:
            c = cage_from_obj(o)
            pl = o.Placement
            V = np.array([tuple(pl.multVec(App.Vector(*p))) for p in c.V]) * flip
            m = Cage(V, [f[::-1] for f in c.F], dict(c.C))  # aynalama yonu cevirir: yuzleri ters sirala
            n = doc.addObject("Mesh::FeaturePython", "SubD")
            SubDObject(n)
            if App.GuiUp:
                SubDViewProvider(n.ViewObject)
            cage_to_obj(m, n)
            set_surface(n, o.Levels, obj_rules(o))
            n.Label = o.Label + "_Mirror" + "XYZ"[axis]
            n.recompute()
            copies.append(n)
        doc.commitTransaction()
        self.set_object_selection(copies)
        self.refresh()

    def object_solids(self):
        """Object modunda: secili her SubD nesnesi icin B-rep kati olustur (varsa guncelle)."""
        objs = self.selected_objects()
        if not objs:
            return
        doc = objs[0].Document
        doc.openTransaction("SubD solid")
        for o in objs:
            if not any("Samples" in x.PropertiesList for x in o.InList):
                make_solid(o)
        doc.recompute()
        doc.commitTransaction()
        self.refresh()

    def object_delete(self):
        objs = self.selected_objects()
        if not objs:
            return
        self._remove_objects(objs)
        Gui.Selection.clearSelection()
        self.refresh()

    def _remove_objects(self, objs):
        """SubD nesnelerini (bagli katilariyla) belgeden tek geri alma adimiyla siler."""
        doc = objs[0].Document
        doc.openTransaction("SubD delete")
        for o in objs:
            for s in [x for x in o.InList if "Samples" in x.PropertiesList]:
                doc.removeObject(s.Name)
            if o is self.obj:
                self.obj = None
            doc.removeObject(o.Name)
        doc.commitTransaction()

    def _drop_empty(self):
        """Edit oturumu biterken butun cokgenler silinmisse nesne de gider (kafes kaydedilmediginden eski haliyle
        agacta kaliyordu). Oturum icinde Ctrl+Z calissin diye silme aninda degil burada; FreeCAD Undo nesneyi
        oturum basindaki haliyle, gorunur olarak geri getirir."""
        if not _alive(self.obj):
            return
        self._show_edited()
        if self.edit_start is not None:
            self._restore_start()
        self._remove_objects([self.obj])

    def object_pick(self, pos):
        """Imlecin altindaki SubD nesnesi (kati tiklandiysa kaynagi)."""
        info = self.view.getObjectInfo((int(pos[0]), int(pos[1])))
        if not info:
            return None
        o = App.ActiveDocument.getObject(info.get("Object", ""))
        if o is not None and not is_subd(o) and "Samples" in o.PropertiesList and is_subd(o.Source):
            o = o.Source
        return o if is_subd(o) else None

    def object_click(self, pos, shift):
        o = self.object_pick(pos)
        objs = self.selected_objects()
        if o is None:
            objs = objs if shift else []
        elif shift:
            objs = [x for x in objs if x is not o] if o in objs else objs + [o]
        else:
            objs = [o]
        self.set_object_selection(objs)
        self.refresh()

    def object_box_select(self, p0, p1, add=False, sub=False):
        lo, hi = np.minimum(p0, p1), np.maximum(p0, p1)
        found = []
        for o in self.subd_objects():
            bb = o.Mesh.BoundBox
            if not bb.isValid():
                continue
            x, y = tuple(self.view.getPointOnViewport(bb.Center))[:2]
            if lo[0] <= x <= hi[0] and lo[1] <= y <= hi[1]:
                found.append(o)
        objs = self.selected_objects()
        if sub:
            objs = [o for o in objs if o not in found]
        elif add:
            objs = objs + [o for o in found if o not in objs]
        else:
            objs = found
        self.set_object_selection(objs)
        self.refresh()

    def objects_move(self, g, d):
        for o, pl0 in zip(g["objs"], g["pl0"]):
            o.Placement = App.Placement(pl0.Base + App.Vector(*d), pl0.Rotation)
        self.gizmo.place(g["c0"] + App.Vector(*d))
        self.refresh(fast=True)

    def objects_rotate(self, g, axis, deg):
        rot = App.Placement(App.Vector(), App.Rotation(axis, deg), g["c0"])   # c0 etrafinda
        for o, pl0 in zip(g["objs"], g["pl0"]):
            o.Placement = rot.multiply(pl0)
        self.refresh(fast=True)

    def objects_scale(self, g, kk):
        """Dunya eksenlerinde olcek: Placement olceklenemedigi icin kafes noktalari olceklenir."""
        c = np.array(tuple(g["c0"]))
        for o, pl0, c0 in zip(g["objs"], g["pl0"], g["cage0"]):
            Pw = np.array([tuple(pl0.multVec(App.Vector(*p))) for p in c0.V])
            Pw = c + (Pw - c) * kk
            inv = pl0.inverse()
            cage = c0.copy()
            cage.V = np.array([tuple(inv.multVec(App.Vector(*p))) for p in Pw])
            cage_to_obj(cage, o)
            o.recompute()
        self.refresh(fast=True)

    def object_panel_spec(self):
        """Object modunda panel: tek nesne -> konum, donus (Z / Y / X), boyut; birden fazla -> toplam boyut."""
        objs = self.selected_objects()
        if not objs:
            return None
        edit = self.drag is None
        bb = None
        for o in objs:
            b = o.Mesh.BoundBox
            if bb is None:
                bb = App.BoundBox(b)
            else:
                bb.add(b)
        size = [(T("Size ") + n, v, " mm", 3, lambda val, j=j: self.set_object_size(j, val))
                for j, (n, v) in enumerate(zip("XYZ", (bb.XLength, bb.YLength, bb.ZLength)))]
        if len(objs) > 1:
            return T("%d objects") % len(objs), size, edit
        o = objs[0]
        pl = o.Placement
        ypr = pl.Rotation.getYawPitchRoll()
        pos = [(T("Position ") + n, v, " mm", 3, lambda val, j=j: self.set_object_position(j, val))
               for j, (n, v) in enumerate(zip("XYZ", tuple(pl.Base)))]
        rot = [(T("Rotation ") + n, ypr[j], "°", 2, lambda val, j=j: self.set_object_rotation(j, val))
               for j, n in enumerate("ZYX")]
        return o.Label, pos + rot + size, edit

    def set_object_position(self, j, val):
        objs = self.selected_objects()
        if len(objs) != 1:
            return
        o = objs[0]
        o.Document.openTransaction("SubD placement")
        b = list(tuple(o.Placement.Base))
        b[j] = val
        o.Placement = App.Placement(App.Vector(*b), o.Placement.Rotation)
        o.Document.recompute()
        o.Document.commitTransaction()
        self.refresh()

    def set_object_rotation(self, j, val):
        """j: 0 = Z (yaw), 1 = Y (pitch), 2 = X (roll), derece."""
        objs = self.selected_objects()
        if len(objs) != 1:
            return
        o = objs[0]
        ypr = list(o.Placement.Rotation.getYawPitchRoll())
        ypr[j] = val
        r = App.Rotation()
        r.setYawPitchRoll(*ypr)
        o.Document.openTransaction("SubD placement")
        o.Placement = App.Placement(o.Placement.Base, r)
        o.Document.recompute()
        o.Document.commitTransaction()
        self.refresh()

    def set_object_size(self, j, val):
        """Secili nesne(ler)in dunya kutusunun j eksenindeki boyunu, kutunun ortasina gore olcekleyerek val yap."""
        objs = self.selected_objects()
        if not objs or val <= 0:
            return
        bb = None
        for o in objs:
            b = o.Mesh.BoundBox
            if bb is None:
                bb = App.BoundBox(b)
            else:
                bb.add(b)
        cur = (bb.XLength, bb.YLength, bb.ZLength)[j]
        if cur < 1e-9:
            return
        kk = np.ones(3)
        kk[j] = val / cur
        objs[0].Document.openTransaction("SubD size")
        g = dict(objs=objs, pl0=[App.Placement(o.Placement) for o in objs],
                 cage0=[cage_from_obj(o) for o in objs], c0=bb.Center)
        self.objects_scale(g, kk)
        objs[0].Document.recompute()
        objs[0].Document.commitTransaction()
        self.refresh()

    # --- sahne ---
    def _build_scene(self):
        self.root = coin.SoSeparator()
        ps = coin.SoPickStyle()
        ps.style = coin.SoPickStyle.UNPICKABLE  # secimi kendimiz yapiyoruz
        self.root.addChild(ps)

        def hints(crease):
            h = coin.SoShapeHints()
            h.vertexOrdering = coin.SoShapeHints.COUNTERCLOCKWISE
            h.creaseAngle = crease
            h.faceType = coin.SoShapeHints.UNKNOWN_FACE_TYPE   # icbukey kafes cokgenleri dogru ucgenlensin
            return h

        def material(diffuse, spec=(0.3, 0.3, 0.3), shine=0.4, transp=0.0):
            m = coin.SoMaterial()
            m.diffuseColor = diffuse
            m.specularColor = spec
            m.shininess = shine
            m.transparency = transp
            return m

        def offset():
            o = coin.SoPolygonOffset()
            o.factor, o.units = 1.0, 1.0
            return o

        self.c_coord = coin.SoCoordinate3()
        self.s_coord = coin.SoCoordinate3()
        self.s_faces = coin.SoIndexedFaceSet()
        self.c_faces = coin.SoIndexedFaceSet()
        self.sel_faces = coin.SoIndexedFaceSet()
        self.c_lines = coin.SoIndexedLineSet()
        self.sel_lines = coin.SoIndexedLineSet()

        # hayalet modunda bu iki malzeme yari saydam olur: arkadaki kenar / noktalar gorunur
        self.smooth_mat = material((0.72, 0.75, 0.80), (0.5, 0.5, 0.5), 0.6)
        self.cage_mat = material(tuple(self.cage_color.getRgbF()[:3]))
        # FreeCAD saydam yuzeyleri her karede CPU'da ucgen ucgen siralar (21 bin dortgende 16 -> 55 ms): goruntu temiz ama
        # dondururken agir. Yalniz bu yuzeylerde hareket surerken siralama kapanir (_ghost_tick), dururken geri gelir.
        self.ghost_blend = coin.SoTransparencyType()
        self.ghost_blend.value = coin.SoTransparencyType.SORTED_OBJECT_SORTED_TRIANGLE_BLEND
        self._ghost_cam, self._ghost_still = None, 0
        self.s_switch = coin.SoSwitch()
        self.zebra_switch = coin.SoSwitch()             # yansima seritleri: ortam (kure) dokusu
        self.zebra_switch.addChild(_zebra_texture())
        self.matcap_switch = coin.SoSwitch()               # matcap: isik dokuda, sahne isiklari yerine (REPLACE)
        mc = coin.SoGroup()
        mc.addChild(coin.SoTextureCoordinateEnvironment())
        self.matcap_tex = coin.SoTexture2()
        self.matcap_tex.model = coin.SoTexture2.REPLACE
        mc.addChild(self.matcap_tex)
        mc.addChild(_matcap_shader())
        self.matcap_switch.addChild(mc)
        self.s_normal = coin.SoNormal()                    # kose normalleri: surface_normals
        nb = coin.SoNormalBinding()
        nb.value = coin.SoNormalBinding.PER_VERTEX_INDEXED
        self.env_switch = coin.SoSwitch()
        eg, self.env_tex, self.env_rows, self.env_ctr, self.env_rad, self.env_rough_p, self.env_metal_p, \
            self.env_texel, self.env_base = _env_nodes()
        self.env_base.value = tuple(self.env_color.getRgbF()[:3])
        self.env_metal_p.value = self.env_metal
        self.env_rough_p.value = self.env_rough
        self.env_switch.addChild(eg)
        self.s_switch.addChild(_sep(self.ghost_blend, hints(0.5), self.smooth_mat, offset(), self.zebra_switch, self.matcap_switch,
                                    self.env_switch, self.s_coord, self.s_normal, nb, self.s_faces))
        if self.matcap not in MATCAPS:
            self.matcap = ""
        if self.env_on and self._env_load() is not None:   # hatirlanan resim tasinmis / silinmis
            self.env_on = False
        self._matcap_apply()
        self.c_switch = coin.SoSwitch()
        self.grp_faces = coin.SoIndexedFaceSet()           # grup cokgenleri: her grubun kendi renginde
        self.grp_mat = coin.SoMaterial()
        self.grp_mat.specularColor = (0.3, 0.3, 0.3)
        self.grp_mat.shininess = 0.4
        grp_bind = coin.SoMaterialBinding()
        grp_bind.value = coin.SoMaterialBinding.PER_FACE
        self._grp_key = None
        self.c_switch.addChild(_sep(
            _sep(self.ghost_blend, hints(0.0), self.cage_mat, offset(), self.c_coord, self.c_faces),
            _sep(hints(0.0), self.grp_mat, grp_bind, offset(), self.c_coord, self.grp_faces),
            _sep(hints(0.0), material((1.0, 0.55, 0.15)), self.c_coord, self.sel_faces)))

        def lines(color, width, node):
            lm = coin.SoLightModel()
            lm.model = coin.SoLightModel.BASE_COLOR
            bc = coin.SoBaseColor()
            bc.rgb = color
            ds = coin.SoDrawStyle()
            ds.lineWidth = width
            return _sep(lm, bc, ds, self.c_coord, node)

        def points(color, size, coord, node):
            lm = coin.SoLightModel()
            lm.model = coin.SoLightModel.BASE_COLOR
            bc = coin.SoBaseColor()
            bc.rgb = color
            ds = coin.SoDrawStyle()
            ds.pointSize = size
            return _sep(lm, bc, ds, coord, node)

        # kafes / yuzey nesnenin yerel koordinatinda: duzenlenen nesnenin Placement'i ile dunyaya tasinir.
        # Object modunda bu grup gizlenir (nesneler belgedeki halleriyle gorunur). Gizmo dunya koordinatinda.
        self.geo_switch = coin.SoSwitch()
        geo = coin.SoSeparator()
        self.geo_switch.addChild(geo)
        self.geo_switch.whichChild = 0
        self.xform = coin.SoTransform()
        geo.addChild(self.xform)
        self.p_switch = coin.SoSwitch()  # nokta modunda tum kafes noktalari
        self.p_switch.addChild(points((0.10, 0.15, 0.35), 6, self.c_coord, coin.SoPointSet()))
        self.selp_coord = coin.SoCoordinate3()
        self.selp_set = coin.SoPointSet()
        geo.addChild(self.p_switch)
        geo.addChild(points((1.0, 0.45, 0.0), 10, self.selp_coord, self.selp_set))

        self.crease_lines = coin.SoIndexedLineSet()
        for n in (self.s_switch, self.c_switch,
                  lines((0.10, 0.15, 0.35), 1.5, self.c_lines),
                  lines((0.85, 0.20, 0.85), 3.0, self.crease_lines),
                  lines((1.0, 0.45, 0.0), 3.0, self.sel_lines)):
            geo.addChild(n)
        # Rail bridge kilavuzlari: egri borunun icinden gecer, ustte cizilir (SoAnnotation). Soluk: duzenlenmeyen
        # kilavuzlar (cift tiklayinca duzenlenir); parlak: duzenlenen, ara noktalari tutamak, secili parca turuncu.
        self.rail_dcoord, self.rail_dlines = coin.SoCoordinate3(), coin.SoIndexedLineSet()
        self.rail_coord, self.rail_lines = coin.SoCoordinate3(), coin.SoIndexedLineSet()
        self.rail_span = coin.SoIndexedLineSet()
        self.rail_hcoord, self.rail_hset = coin.SoCoordinate3(), coin.SoPointSet()
        self.rail_scoord, self.rail_sset = coin.SoCoordinate3(), coin.SoPointSet()

        def draw(color, width, coord, node, point=False):
            lm = coin.SoLightModel()
            lm.model = coin.SoLightModel.BASE_COLOR
            bc = coin.SoBaseColor()
            bc.rgb = color
            ds = coin.SoDrawStyle()
            if point:
                ds.pointSize = width
            else:
                ds.lineWidth = width
            return _sep(lm, bc, ds, coord, node)
        ann = coin.SoAnnotation()
        for n in (draw((0.62, 0.55, 0.30), 1.5, self.rail_dcoord, self.rail_dlines),
                  draw((0.98, 0.80, 0.15), 2.5, self.rail_coord, self.rail_lines),
                  draw((1.0, 0.45, 0.0), 4.0, self.rail_coord, self.rail_span),
                  draw((0.98, 0.80, 0.15), 9, self.rail_hcoord, self.rail_hset, True),
                  draw((1.0, 0.30, 0.0), 13, self.rail_scoord, self.rail_sset, True)):
            ann.addChild(n)
        geo.addChild(ann)
        self.root.addChild(self.geo_switch)
        self.root.addChild(self.gizmo.root)
        self.parent = _scene_parent(self.view)
        self.parent.addChild(self.root)

    def _build_dock(self):
        self.dock = QtWidgets.QDockWidget(T("SubD Options"))
        self.dock.setObjectName("SubDModelerOptions_%s" % getattr(self, "doc_name", ""))   # saveState: adsiz dock uyarisi
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        self.status = QtWidgets.QLabel()
        # Araclar ve modlar gorunumdeki menulerde; bu Options paneli olculer, ayarlar ve kisayol listesi
        lay.addWidget(self._measure_box())
        lay.addWidget(self._align_group())
        lay.addWidget(self._bevel_group())
        lay.addWidget(self._box_group())
        lay.addWidget(self._surface_group())
        lay.addWidget(self._solid_group())
        lay.addWidget(self._gizmo_group())
        lay.addWidget(self._panel_group())
        keys = QtWidgets.QLabel(T(HELP))
        keys.setStyleSheet("color:#9a9a9a;")
        lay.addWidget(keys)
        lay.addWidget(self.status)
        lay.addStretch()
        # kaydirilabilir: icerik uzun olsa da dock FreeCAD penceresini ekrandan tasirmasin
        sa = QtWidgets.QScrollArea()
        sa.setWidgetResizable(True)
        sa.setFrameShape(QtWidgets.QFrame.NoFrame)
        sa.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        sa.setWidget(w)
        self.dock.setWidget(sa)
        Gui.getMainWindow().addDockWidget(QtCore.Qt.RightDockWidgetArea, self.dock)
        self.dock.hide()   # gunluk is viewport'taki islem panelinde; buyuk panel Options (N) ile acilir

    @staticmethod
    def _spin(lo, hi, value=0.0, suffix=" mm", decimals=3):
        b = QtWidgets.QDoubleSpinBox()
        b.setRange(lo, hi)
        b.setDecimals(decimals)
        b.setSuffix(suffix)
        b.setKeyboardTracking(False)
        b.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        b.setValue(value)
        return b

    def _measure_box(self):
        """Olcu paneli: secimin merkezi ve boyutu (yazip Enter = uygula), yakalama adimlari, canli deger."""
        g = QtWidgets.QGroupBox(T("Measure"))
        grid = QtWidgets.QGridLayout(g)
        for j, (t, col) in enumerate((("X", "#ef4444"), ("Y", "#22c55e"), ("Z", "#3b82f6"))):
            lb = QtWidgets.QLabel(t)
            lb.setStyleSheet("color:%s; font-weight:bold;" % col)
            lb.setAlignment(QtCore.Qt.AlignCenter)
            grid.addWidget(lb, 0, j + 1)
        self.center_boxes, self.size_boxes = [], []
        for row, (name, boxes, fn, lo) in enumerate(((T("Center"), self.center_boxes, self.set_center, -1e6),
                                                     (T("Size"), self.size_boxes, self.set_size, 0.0))):
            grid.addWidget(QtWidgets.QLabel(name), row + 1, 0)
            for j in range(3):
                b = self._spin(lo, 1e6)
                b.editingFinished.connect(lambda j=j, b=b, fn=fn: fn(j, b.value()))
                grid.addWidget(b, row + 1, j + 1)
                boxes.append(b)
        grid.addWidget(QtWidgets.QLabel("Crease (C)"), 6, 0)
        self.sharp_box = self._spin(0.0, CREASE_HARD, CREASE_HARD, "", 2)
        self.sharp_box.setToolTip(T("Sharpness: %g = fully sharp, 1-3 = semi-sharp (soft chamfer)") % CREASE_HARD)
        grid.addWidget(self.sharp_box, 6, 1)
        b = QtWidgets.QPushButton(T("Apply / remove"))
        b.clicked.connect(lambda checked=False: self.crease_selected())
        grid.addWidget(b, 6, 2, 1, 2)
        grid.addWidget(QtWidgets.QLabel(T("Ctrl step")), 3, 0)
        self.step_box = self._spin(0.0, 1e4, 1.0)
        self.angle_box = self._spin(0.0, 180.0, 15.0, " °", 2)
        grid.addWidget(self.step_box, 3, 1)
        grid.addWidget(self.angle_box, 3, 2)
        self.info_label = QtWidgets.QLabel()
        self.info_label.setStyleSheet("color:#f59e0b; font-weight:bold;")
        grid.addWidget(self.info_label, 4, 0, 1, 4)
        self.dims_label = QtWidgets.QLabel()
        grid.addWidget(self.dims_label, 5, 0, 1, 4)
        return g

    @staticmethod
    def _align_icon(vertical):
        """Hizalama dugmesinin simgesi: bir cizgi uzerine dizilmis uc nokta (yatay ya da dikey); onbellekte."""
        if vertical in _align_icons:
            return _align_icons[vertical]
        pm = QtGui.QPixmap(20, 20)
        pm.fill(QtCore.Qt.transparent)
        p = QtGui.QPainter(pm)
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        p.setPen(QtGui.QPen(QtGui.QColor("#9a9a9a"), 1.5, QtCore.Qt.DashLine))
        p.drawLine(QtCore.QPointF(10, 2), QtCore.QPointF(10, 18)) if vertical else \
            p.drawLine(QtCore.QPointF(2, 10), QtCore.QPointF(18, 10))
        p.setPen(QtCore.Qt.NoPen)
        p.setBrush(QtGui.QColor("#f5c518"))
        for k in (4.5, 10, 15.5):
            p.drawEllipse(QtCore.QPointF(10, k) if vertical else QtCore.QPointF(k, 10), 2.6, 2.6)
        p.end()
        _align_icons[vertical] = QtGui.QIcon(pm)
        return _align_icons[vertical]

    def _align_group(self):
        """Secili noktalari yatay / dikey hizala: kucuk kare simgeli dugmeler (align_points)."""
        g = QtWidgets.QGroupBox(T("Align"))
        row = QtWidgets.QHBoxLayout(g)
        for vertical, tip in ((False, "Align horizontally: the selected points get the height of the first clicked "
                                      "one on screen"),
                              (True, "Align vertically: the selected points line up with the first clicked one on "
                                     "screen")):
            b = QtWidgets.QToolButton()
            b.setIcon(self._align_icon(vertical))
            b.setIconSize(QtCore.QSize(20, 20))
            b.setFixedSize(28, 28)
            b.setToolTip(T(tip))
            b.clicked.connect(lambda checked=False, v=vertical: self.align_points(v))
            row.addWidget(b)
        row.addStretch()
        return g

    def align_points(self, vertical):
        """Secili (en az iki) noktayi ekranda yatay ya da dikey hizala: ekranin dikey (yatay hizada) ya da yatay
        (dikey hizada) yonune en yakin dunya ekseninde hepsi ilk tiklanan noktaya gelir (o eksende olcek 0).
        Tiklanmis nokta yoksa (yalniz kutu secimi) ortalamaya gelir."""
        if self.drag or self.grab:
            return
        cam = self.view.getCameraNode().orientation.getValue()
        up = np.array(cam.multVec(coin.SbVec3f(0, 1, 0)).getValue())
        right = np.array(cam.multVec(coin.SbVec3f(1, 0, 0)).getValue())
        axis = int(np.argmax(np.abs(right if vertical else up)))
        msg = (T("Align vertically: %s") if vertical else T("Align horizontally: %s")) % "XYZ"[axis]
        if self.rail is not None:                          # kilavuz duzenleniyor: secili kilavuz dugumleri
            self._rail_align(axis, msg)
            return
        if self.mode in ("face", "edge"):
            self._align_centers(axis, msg)
            return
        verts = self.selected_verts()
        if len(verts) < 2:
            self.info = T("Align: select two or more points")
            self.refresh()
            return
        P = self.l2w(self.cage.V[verts])
        k = np.ones(3)
        k[axis] = 0.0
        first = self._first_selected(self.sel_v)
        c = self.cage.V[first].copy() if first is not None else self.w2l(P.mean(0))
        self.push()
        self.last_op = None
        self.apply_scale(self.move_state(verts), c, k)
        self.info = msg
        self.refresh()

    def _align_centers(self, axis, msg):
        """Poligon / kenar modunda: secili poligonlarin (kenarlarin) merkezlerini hizala. Birbirine degen secili
        elemanlar tek parca sayilir (ortak noktalar yirtilmasin); her parca bozulmadan tasinir, merkezi (noktalarinin
        ortalamasi) ilk tiklanan elemanin parcasinin merkezine, tiklanmis eleman yoksa parcalarin ortalamasina gelir."""
        items = ([set(self.cage.F[fi]) for fi in self.sel] if self.mode == "face" else [set(e) for e in self.sel_e])
        groups = []
        for it in items:                                   # ortak noktasi olanlari birlestir
            hit = [g for g in groups if g & it]
            for g in hit:
                groups.remove(g)
                it |= g
            groups.append(it)
        if len(groups) < 2:
            self.info = T("Align: select two or more separate polygons or edges")
            self.refresh()
            return
        centers = [self.l2w(self.cage.V[sorted(g)]).mean(0)[axis] for g in groups]
        first = self._first_selected(self.sel if self.mode == "face" else self.sel_e)
        if first is None:
            val = float(np.mean(centers))
        else:
            first = set(self.cage.F[first]) if self.mode == "face" else set(first)
            val = float(next(c for g, c in zip(groups, centers) if g & first))
        self.push()
        self.last_op = None
        for g, c in zip(groups, centers):
            d = np.zeros(3)
            d[axis] = val - c
            self.apply_move(self.move_state(sorted(g)), self.w2l_dir(d))
        self.info = msg
        self.refresh()

    def _rail_align(self, axis, msg):
        """Secili kilavuz dugumlerini (Shift ile, A / B dahil) dunya ekseninde ayni degere getir: secimde A ya da
        B varsa (sabit uclar) onun degeri, yoksa ilk secilen dugumun degeri. Halkalar yeniden yerlesir."""
        r = self._rail()
        if r is None:
            return
        N, _, _, _ = rail_nodes(self.cage, r)
        pick = self._rail_pick(r, N)
        movable = [j for j in pick if self._rail_movable(r, N, j)]
        if len(pick) < 2 or not movable:
            self.info = T("Align: Shift-click two or more guide points (the end points too)")
            self.refresh()
            return
        W = self.l2w(N[pick])
        fixed = [k for k, j in enumerate(pick) if not self._rail_movable(r, N, j)]
        val = W[fixed, axis].mean() if fixed else W[0, axis]

        def set_(r):
            for j in movable:
                w = self.l2w(np.array(r["pts"][j - 1], float))
                w[axis] = val
                r["pts"][j - 1] = [float(x) for x in self.w2l(w)]
        if self._rail_change(set_):
            self.info = msg

    def _bevel_group(self):
        g = QtWidgets.QGroupBox(T("Bevel / Smooth Shift (B, I = inset only)"))
        row = QtWidgets.QHBoxLayout(g)
        self.inset_box = self._spin(-1e4, 1e4, 0.0)
        self.shift_box = self._spin(-1e4, 1e4, 0.0)
        for name, b in (("Inset", self.inset_box), ("Shift", self.shift_box)):
            row.addWidget(QtWidgets.QLabel(name))
            row.addWidget(b)
        b = QtWidgets.QPushButton(T("Apply"))
        b.setToolTip(T("Bevel the selected faces with these values"))
        b.clicked.connect(lambda checked=False: self.bevel(values=(self.inset_box.value(), self.shift_box.value())))
        row.addWidget(b)
        return g

    def _surface_group(self):
        """SubD yuzeyi: duzenleme / son seviye, acik kenar kurali, ucgen yumusatma. Edit modunda duzenlenen
        nesneye (kayitta), Object modunda secili nesnelere hemen uygulanir; yeni nesneler de bunlari alir."""
        g = QtWidgets.QGroupBox(T("SubD surface"))
        grid = QtWidgets.QGridLayout(g)
        grid.setVerticalSpacing(2)
        self.edit_lv_box, self.final_lv_box = QtWidgets.QSpinBox(), QtWidgets.QSpinBox()
        for col, (name, b, tip) in enumerate((
                ("Edit level", self.edit_lv_box, "Subdivision level shown while editing (lower = faster)"),
                ("Final level", self.final_lv_box, "Subdivision level of the SubD object in the document"))):
            b.setRange(1, 5)
            b.setToolTip(T(tip))
            grid.addWidget(QtWidgets.QLabel(T(name)), 0, 2 * col)
            grid.addWidget(b, 0, 2 * col + 1)
        self.boundary_box = QtWidgets.QComboBox()
        for r in BOUNDARY_RULES:
            self.boundary_box.addItem(T(r), r)
        self.boundary_box.setToolTip(T("Open edges: Edge and corner = one-face corners stay sharp, "
                                       "Edge only = they round too"))
        grid.addWidget(QtWidgets.QLabel(T("Open edges")), 1, 0)
        grid.addWidget(self.boundary_box, 1, 1, 1, 3)
        self.smooth_tri_box = QtWidgets.QCheckBox(T("Smooth triangles"))
        self.smooth_tri_box.setToolTip(T("Loop-like weights next to triangles: smoother triangle fans "
                                         "(e.g. cylinder caps)"))
        grid.addWidget(self.smooth_tri_box, 2, 0, 1, 2)
        self.zebra_box = QtWidgets.QCheckBox(T("Zebra stripes"))
        self.zebra_box.setToolTip(T("Reflection stripes on the surface: kinks show as broken stripes (G0), "
                                    "corners in the stripes as G1, smooth stripes as G2"))
        grid.addWidget(self.zebra_box, 2, 2, 1, 2)
        self.env_box = QtWidgets.QCheckBox(T("Environment reflection"))
        self.env_box.setToolTip(T("The smooth preview reflects an image of the surroundings: an equirectangular (2:1) "
                                  "panorama, best as .hdr; JPG and PNG also work"))
        grid.addWidget(self.env_box, 3, 0, 1, 2)
        self.env_btn = QtWidgets.QPushButton(T("Load image..."))
        self.env_btn.setToolTip(self.env_path or T("Environment image"))
        grid.addWidget(self.env_btn, 3, 2, 1, 2)
        for row, name, tip, val, fn in (
                (4, "Roughness", "0: mirror, 1: fully blurred reflections", self.env_rough,
                 lambda v: self.set_env_look(rough=v)),
                (5, "Metal", "0: plastic (surface color with a soft reflection), 1: metal (reflection only)",
                 self.env_metal, lambda v: self.set_env_look(metal=v))):
            sp = QtWidgets.QDoubleSpinBox()
            sp.setRange(0.0, 1.0)
            sp.setSingleStep(0.05)
            sp.setDecimals(2)
            sp.setKeyboardTracking(False)
            sp.setValue(val)
            sp.setToolTip(T(tip))
            sp.valueChanged.connect(fn)
            grid.addWidget(QtWidgets.QLabel(T(name)), row, 0, 1, 2)
            grid.addWidget(sp, row, 2, 1, 2)
        self.env_color_btn = QtWidgets.QPushButton()
        self.env_color_btn.setToolTip(T("Surface color for the plastic look (no effect at Metal 1)"))
        self.env_color_btn.clicked.connect(lambda checked=False: self.env_pick_color())
        self._env_swatch()
        grid.addWidget(QtWidgets.QLabel(T("Surface color")), 6, 0, 1, 2)
        grid.addWidget(self.env_color_btn, 6, 2, 1, 2)
        self.cage_color_btn = QtWidgets.QPushButton()
        self.cage_color_btn.setToolTip(T("Color of the cage faces when Subpatch is off"))
        self.cage_color_btn.setStyleSheet("background: %s;" % self.cage_color.name())
        self.cage_color_btn.clicked.connect(lambda checked=False: self.cage_pick_color())
        grid.addWidget(QtWidgets.QLabel(T("Cage color")), 7, 0, 1, 2)
        grid.addWidget(self.cage_color_btn, 7, 2, 1, 2)
        self.link_box = QtWidgets.QCheckBox(T("Follow sketch"))
        self.link_box.setToolTip(T("From edge polygons follow their source sketch / edges when it changes; off = "
                                   "unlinked (the cage keeps its shape)"))
        grid.addWidget(self.link_box, 8, 0, 1, 4)
        self.link_box.toggled.connect(self.set_sketch_link)
        self.env_box.toggled.connect(lambda on: self.set_env(bool(on)))
        self.env_btn.clicked.connect(lambda checked=False: self.env_pick())
        self._sync_surface_widgets()
        self.edit_lv_box.valueChanged.connect(self.set_edit_levels)
        self.final_lv_box.valueChanged.connect(lambda v: self.change_surface(levels=v))
        self.boundary_box.currentIndexChanged.connect(
            lambda i: self.change_surface(rules=(BOUNDARY_RULES[i], self.rules[1])))
        self.smooth_tri_box.toggled.connect(lambda on: self.change_surface(rules=(self.rules[0], bool(on))))
        self.zebra_box.toggled.connect(lambda on: self.toggle_zebra() if bool(on) != self.zebra else None)
        return g

    def _solid_group(self):
        """Kati (STEP) ayarlari: hassasiyet on ayari ve sureklilik. Yeni katilar bunlari alir; duzenlenen (ya da
        Object modunda secili) nesnenin katisi varsa o da guncellenir."""
        g = QtWidgets.QGroupBox(T("Solid (STEP)"))
        grid = QtWidgets.QGridLayout(g)
        self.precision_box, self.continuity_box = QtWidgets.QComboBox(), QtWidgets.QComboBox()
        for n in PRECISION:
            self.precision_box.addItem(T(n), n)
        self.precision_box.setToolTip(T("Draft / Normal / Fine: samples per patch, G2 refine levels and "
                                        "sewing tolerance"))
        for n, tip in zip(CONTINUITY, ("G0: fast, every patch fitted (tiny kinks between patches)",
                                       "G1: regular patches exact, the rest fitted to the limit surface",
                                       "G2: as G1, irregular patches refined into more exact ones (curvature "
                                       "continuous except at 3 / 5+ edge points, where SubD itself is G1)")):
            self.continuity_box.addItem(n, n)
            self.continuity_box.setItemData(self.continuity_box.count() - 1, T(tip), QtCore.Qt.ToolTipRole)
        prm = _param()
        self.precision_box.setCurrentIndex(max(0, list(PRECISION).index(prm.GetString("Precision", "Normal"))
                                               if prm.GetString("Precision", "Normal") in PRECISION else 1))
        self.continuity_box.setCurrentIndex(CONTINUITY.index(prm.GetString("Continuity", "G2"))
                                            if prm.GetString("Continuity", "G2") in CONTINUITY else 2)
        grid.addWidget(QtWidgets.QLabel(T("Precision")), 0, 0)
        grid.addWidget(self.precision_box, 0, 1)
        grid.addWidget(QtWidgets.QLabel(T("Continuity")), 0, 2)
        grid.addWidget(self.continuity_box, 0, 3)
        self.precision_box.currentIndexChanged.connect(lambda i: self.set_solid_settings())
        self.continuity_box.currentIndexChanged.connect(lambda i: self.set_solid_settings())
        return g

    def _gizmo_group(self):
        """Gizmo boyu (ekranda piksel): Move / Rotate / Stretch gizmolarinin hepsi; kaydedilir."""
        g = QtWidgets.QGroupBox(T("Gizmo"))
        row = QtWidgets.QHBoxLayout(g)
        row.addWidget(QtWidgets.QLabel(T("Size")))
        self.gizmo_box = QtWidgets.QSpinBox()
        self.gizmo_box.setRange(30, 400)
        self.gizmo_box.setSingleStep(10)
        self.gizmo_box.setSuffix(" px")
        self.gizmo_box.setValue(int(round(self.gizmo.px)))
        self.gizmo_box.setToolTip(T("Length of the gizmo arrows on screen (default %d px)") % GIZMO_PX)
        self.gizmo_box.valueChanged.connect(self.set_gizmo_size)
        row.addWidget(self.gizmo_box)
        b = QtWidgets.QPushButton(T("Default"))
        b.clicked.connect(lambda checked=False: self.gizmo_box.setValue(int(GIZMO_PX)))
        row.addWidget(b)
        row.addStretch()
        return g

    def set_gizmo_size(self, px):
        self.gizmo.px = float(px)
        _param().SetFloat("GizmoSize", float(px))
        self.gizmo.update_scale()

    def _panel_group(self):
        """Islem panelinin (secim / arac degerleri) yeri: gorunumun sol ya da sag alti; kaydedilir. Kucuk
        ekranlarda sol alttaki panel modelin onunu kapatabiliyor."""
        g = QtWidgets.QGroupBox(T("Tool panel"))
        row = QtWidgets.QHBoxLayout(g)
        row.addWidget(QtWidgets.QLabel(T("Position")))
        self.panel_side_box = QtWidgets.QComboBox()
        self.panel_side_box.addItem(T("Left"), "left")
        self.panel_side_box.addItem(T("Right"), "right")
        self.panel_side_box.setToolTip(T("Where the panel with the values of the running tool and of the selection "
                                         "sits in the 3D view: bottom left or bottom right"))
        self.panel_side_box.setCurrentIndex(1 if _param().GetString("PanelSide", "left") == "right" else 0)
        self.panel_side_box.currentIndexChanged.connect(self.set_panel_side)
        row.addWidget(self.panel_side_box)
        row.addStretch()
        return g

    def set_panel_side(self, index):
        _param().SetString("PanelSide", "right" if index == 1 else "left")
        if self.overlay is not None:
            self.overlay.relayout()

    def _sync_surface_widgets(self):
        """Options'taki yuzey ayarlarini oturumun degerlerine getir (baska nesneye gecince)."""
        if not hasattr(self, "edit_lv_box"):
            return
        for b, val in ((self.edit_lv_box, self.edit_levels), (self.final_lv_box, self.levels)):
            b.blockSignals(True)
            b.setValue(int(val))
            b.blockSignals(False)
        for b, fn in ((self.boundary_box, lambda: b.setCurrentIndex(
                           BOUNDARY_RULES.index(self.rules[0]) if self.rules[0] in BOUNDARY_RULES else 0)),
                      (self.smooth_tri_box, lambda: b.setChecked(self.rules[1])),
                      (self.zebra_box, lambda: b.setChecked(self.zebra)),
                      (self.env_box, lambda: b.setChecked(self.env_on)),
                      (self.link_box, lambda: b.setChecked(bool(getattr(self.obj, "SketchLink", True))))):
            b.blockSignals(True)
            fn()
            b.blockSignals(False)

    def set_sketch_link(self, on):
        """Options > Follow sketch: duzenlenen nesnenin bagi acik / kapali (nesnenin SketchLink ozelligi)."""
        if self.obj is None or not _alive(self.obj):
            return
        _link_props(self.obj)
        self.obj.SketchLink = bool(on)
        if on:
            self.links_refresh()

    def set_edit_levels(self, v):
        self.edit_levels = min(max(int(v), 1), 5)
        _param().SetInt("EditLevels", self.edit_levels)
        self.refresh()

    def change_surface(self, levels=None, rules=None):
        """Son seviye / alt bolme kurallari. Edit modunda kafes kaydedilince nesneye yazilir; Object modunda
        secili nesnelere hemen uygulanir."""
        if levels is not None:
            self.levels = min(max(int(levels), 1), 5)
        if rules is not None:
            self.rules = rules
        if self.mode == "object":
            objs = self.selected_objects()
            if objs:
                doc = objs[0].Document
                doc.openTransaction("SubD surface")
                for o in objs:
                    set_surface(o, self.levels, self.rules)
                    o.recompute()
                doc.commitTransaction()
        else:
            self.dirty = True
        self.refresh()

    def _solids_of(self, objs):
        return [x for o in objs for x in o.InList if "Samples" in x.PropertiesList and x.Source == o]

    def set_solid_settings(self):
        """Hassasiyet / sureklilik: sonraki katilar icin kaydedilir; duzenlenen ya da secili nesnelerin
        katilari hemen guncellenir."""
        prec, cont = self.precision_box.currentData(), self.continuity_box.currentData()
        _param().SetString("Precision", prec)
        _param().SetString("Continuity", cont)
        objs = self.selected_objects() if self.mode == "object" else [self.obj] if self.obj is not None else []
        solids = self._solids_of([o for o in objs if o is not None and o.Document is not None])
        if solids:
            doc = solids[0].Document
            doc.openTransaction("SubD solid settings")
            for x in solids:
                solid_settings(x, prec, cont)
                x.recompute()
            doc.commitTransaction()
            self.info = T("Solid: %s, %s") % (T(prec), cont)
        self.refresh()

    def _box_group(self):
        g = QtWidgets.QGroupBox(T("New box (width x depth x height)"))
        row = QtWidgets.QHBoxLayout(g)
        self.box_size = [self._spin(0.001, 1e6, self.box_vals[n]) for n in self.BOX_SIZE]
        for b in self.box_size:
            row.addWidget(b)
        b = QtWidgets.QPushButton(T("Create"))
        b.clicked.connect(lambda checked=False: self.add_box())   # her zaman yeni nesne (kafesi silmez)
        row.addWidget(b)
        return g

    def _update_measure(self, fast):
        verts = self.selected_verts()
        V = self.cage.V
        P = self.l2w(V[verts]) if verts else None          # dunya koordinatinda
        vals = [P.mean(0), P.max(0) - P.min(0)] if P is not None else [None, None]
        for boxes, val in zip((self.center_boxes, self.size_boxes), vals):
            for j, b in enumerate(boxes):
                b.setEnabled(val is not None)
                if val is not None and not b.hasFocus():
                    b.blockSignals(True)
                    b.setValue(float(val[j]))
                    b.blockSignals(False)
        self.info_label.setText(self.info)
        if fast or not self.cage.F:
            return
        # CAD icin onemli: SubPatch yuzeyi kafesin icinde kalir, gercek parca olcusu yuzeyinkidir
        S, _ = subdivide(V, self.cage.F, self.levels, self.cage.C, self.rules)
        cd, sd = np.ptp(V, 0), np.ptp(S, 0)
        self.dims_label.setText(T("Cage:  %.2f x %.2f x %.2f mm\nSurface:  %.2f x %.2f x %.2f mm") % (tuple(cd) + tuple(sd)))

    def _groups_backfill_color(self):
        """Renk eklenmeden once kurulan gruplara renk verir (kayda yazilir, tutamak da o renge boyanir)."""
        if not _alive(self.obj):
            return
        ensure_view_provider(self.obj)                     # grup tutamaklari agacta SubD'nin altinda
        if "GroupHandles" in self.obj.PropertiesList:      # eski surumun bagimlilik baglari (dongu yapardi)
            _quiet(lambda: drop_handle_links(self.obj))
        for g in self.cage.groups:
            if not g.get("color"):
                g["color"] = list(group_rgb(g))
                h = self.obj.Document.getObject(g["name"])
                if h is not None and App.GuiUp:
                    _quiet(lambda: setattr(h.ViewObject, "ShapeColor", tuple(g["color"])))
            h = self.obj.Document.getObject(g["name"])
            if h is not None and not h.Label2:              # renk onceden verilmis ama kare yok
                h.Label2 = group_square(group_rgb(g))
            if h is not None:
                _quiet(lambda: handle_dim_setup(h))         # baska bir gruba bagliysa mesafe ifadesi

    def _group_faces(self, fast):
        """Tum kosesi bir grupta olan cokgenleri grubun rengine boyar. Uyeler kayitli konumlarindan bulunur;
        surukleme sirasinda (fast) uyelik yeniden hesaplanmaz."""
        gs, F, V = self.cage.groups, self.cage.F, self.cage.V
        key = (len(V), len(F), tuple((g["name"], len(g["pts"]), g.get("locked", True)) for g in gs))
        if fast and key == self._grp_key:
            return
        self._grp_key = key
        faces, cols = [], []
        self._grp_sets = group_member_sets(V, gs) if gs else {}
        for g in gs:
            mem, col = self._grp_sets.get(g["name"], set()), group_rgb(g)
            for f in F:
                if mem and all(v in mem for v in f):
                    faces.append(f)
                    cols.append(col)
        _set_index(self.grp_faces, _poly_index(faces))
        if cols:
            self.grp_mat.diffuseColor.setValues(0, len(cols), cols)
            self.grp_mat.diffuseColor.setNum(len(cols))

    def refresh(self, fast=False):
        if self.mode == "object":
            if not (self.drag and fast):               # surukleme sirasinda gizmo nesneyle birlikte tasinir
                c = self.objects_center()
                self.gizmo.place(App.Vector(*c) if c is not None else None)
            self._dim_tick()
            if self.overlay is not None:
                self.overlay.sync()
            self.status.setText(T("Object mode: %d selected") % len(self.selected_objects()))
            return
        V, F = self.cage.V, self.cage.F
        _set_coords(self.c_coord, V if len(V) else np.zeros((1, 3)))
        _set_index(self.c_faces, _poly_index(F))
        _set_index(self.c_lines, _poly_index([list(e) for e in self.cage.edges()]))
        _set_index(self.crease_lines, _poly_index([list(e) for e in self.cage.C]))
        sel = [F[i] for i in sorted(self.sel) if i < len(F)] if self.mode == "face" else []
        _set_index(self.sel_faces, _poly_index(sel))
        self._group_faces(fast)
        if self.mode == "edge":
            _set_index(self.sel_lines, _poly_index([list(e) for e in self.sel_e]))
        else:
            _set_index(self.sel_lines, _poly_index([f + [f[0]] for f in sel]))
        pts = V[sorted(self.sel_v)] if self.mode == "point" and self.sel_v else np.zeros((0, 3))
        _set_coords(self.selp_coord, pts if len(pts) else np.zeros((1, 3)))
        self.selp_set.numPoints = len(pts)
        self.p_switch.whichChild = 0 if self.mode == "point" else -1
        if self.smooth and F:
            lv = self.edit_levels
            S, Q, cr = V, F, self.cage.C
            for _ in range(min(lv, 2) if fast else lv):
                S, Q, cr = catmull_clark(S, Q, cr, self.rules)
            _set_coords(self.s_coord, S)
            self._smooth = (S, Q)
            _set_index(self.s_faces, np.concatenate([Q, -np.ones((len(Q), 1), int)], 1).ravel().tolist())
            Nv = surface_normals(S, Q, cr)
            self.s_normal.vector.setValues(0, Nv.size // 3, Nv.reshape(-1, 3).tolist())
            self.s_normal.vector.setNum(Nv.size // 3)
            idx = np.arange(4 * len(Q)).reshape(-1, 4)
            idx = np.concatenate([idx, -np.ones((len(Q), 1), int)], 1).ravel().tolist()
            self.s_faces.normalIndex.setValues(0, len(idx), idx)
            self.s_faces.normalIndex.setNum(len(idx))
            lo, hi = S.min(0), S.max(0)                    # ortam kuresi modeli sarar
            self.env_ctr.value = tuple(((lo + hi) / 2).tolist())
            self.env_rad.value = float(max(np.linalg.norm(hi - lo), 1e-6)) * ENV_RADIUS
        self.s_switch.whichChild = 0 if self.smooth else -1
        self.c_switch.whichChild = -1 if self.smooth else 0
        verts = self.selected_verts() if not self.grab and self.armed != "magnet" and not self.suspended else []
        self.gizmo.place(App.Vector(*self.l2w(V[verts].mean(0))) if verts else None)
        self._rail_draw()
        self._update_measure(fast)
        self._dim_tick()
        if self.overlay is not None:
            self.overlay.sync()
        count = {"point": len(self.sel_v), "edge": len(self.sel_e), "face": len(self.sel)}[self.mode]
        self.status.setText(T("%d points, %d faces  |  %d selected") % (len(V), len(F), count)
                            + (T("  |  Symmetry X") if self.sym else "") + ("  |  SubPatch" if self.smooth else ""))

    # --- yerel (kafes) <-> dunya koordinati: duzenlenen nesnenin Placement'i ---
    # Kafes noktalari nesnenin yerel koordinatinda saklanir; ekran, isinlar ve gizmo dunya koordinatinda.
    def set_placement(self, pl):
        self._R = np.array(pl.Rotation.toMatrix().A).reshape(4, 4)[:3, :3]
        self._t = np.array(tuple(pl.Base))
        self.xform.translation = tuple(self._t)
        self.xform.rotation.setValue(*pl.Rotation.Q)
        for p, row in zip(self.env_rows, self._R):         # ortam yansitma: yerel yon -> dunya
            p.value = tuple(row.tolist())

    def l2w(self, P):
        return np.asarray(P) @ self._R.T + self._t

    def w2l(self, P):
        return (np.asarray(P) - self._t) @ self._R

    def l2w_dir(self, d):
        return np.asarray(d) @ self._R.T

    def w2l_dir(self, d):
        return np.asarray(d) @ self._R

    def view_dir(self):
        """Bakis yonu, kafesin yerel koordinatinda."""
        return self.w2l_dir(np.array(tuple(self.view.getViewDirection())))

    def world_matrix(self, kk):
        """Dunya eksenlerinde kk olcegi (Stretch / Size), yerel koordinatta: R^T diag(kk) R."""
        return self._R.T @ np.diag(kk) @ self._R

    # --- isin / secim ---
    def ray(self, pos):
        """Imlecten isin, kafesin yerel koordinatinda: (baslangic, birim yon)."""
        p = np.array(tuple(self.view.getPoint(*pos)))
        d = np.array(tuple(self.view.getViewDirection()))
        if self.view.getCameraType() == "Perspective":
            o = np.array(self.view.getCameraNode().position.getValue().getValue())
            d = p - o
        else:
            o = p - d * 1e5
        return self.w2l(o), self.w2l_dir(d / np.linalg.norm(d))

    def pick(self, pos):
        V, F = self.cage.V, self.cage.F
        tri, owner = [], []
        for fi, f in enumerate(F):
            for i in range(1, len(f) - 1):
                tri.append((f[0], f[i], f[i + 1]))
                owner.append(fi)
        if not tri:
            return None
        o, d = self.ray(pos)
        T = V[np.array(tri)]
        e1, e2 = T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]
        p = np.cross(d, e2)
        det = (e1 * p).sum(1)
        ok = np.abs(det) > 1e-12
        inv = 1.0 / np.where(ok, det, 1.0)
        s = o - T[:, 0]
        u = (s * p).sum(1) * inv
        q = np.cross(s, e1)
        v = (q * d).sum(1) * inv
        t = (e2 * q).sum(1) * inv
        hit = ok & (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 0)
        if not hit.any():
            return None
        k = int(np.argmin(np.where(hit, t, np.inf)))
        return owner[k], o + d * t[k]

    def ray_plane(self, pos, c, n):
        o, d = self.ray(pos)
        den = d.dot(n)
        if abs(den) < 1e-9:
            return None
        return o + d * ((c - o).dot(n) / den)

    def project(self, P):
        """Kafes (yerel) noktalarini ekran (viewport) koordinatlarina cevirir: olaylardaki Position ile ayni duzen."""
        W = self.l2w(np.asarray(P, float).reshape(-1, 3))
        return np.array([tuple(self.view.getPointOnViewport(App.Vector(*p)))[:2] for p in W], float)

    def depth(self, P):
        return np.asarray(P).dot(self.view_dir())

    def pick_point(self, pos, radius=12.0):
        V = self.cage.V
        if not len(V):
            return None
        d = np.linalg.norm(self.project(V) - np.array(pos, float), axis=1)
        near = np.where(d <= radius)[0]
        if not len(near):
            return None
        return int(near[np.argmin(self.depth(V[near]))])  # kameraya en yakin olan

    def pick_edge(self, pos, radius=8.0):
        edges = list(self.cage.edges())
        if not edges:
            return None
        V = self.cage.V
        S = self.project(V)
        E = np.array([sorted(e) for e in edges])
        a, b = S[E[:, 0]], S[E[:, 1]]
        p = np.array(pos, float)
        ab = b - a
        t = np.clip(((p - a) * ab).sum(1) / np.maximum((ab * ab).sum(1), 1e-12), 0, 1)
        d = np.linalg.norm(a + ab * t[:, None] - p, axis=1)
        near = np.where(d <= radius)[0]
        if not len(near):
            return None
        mids = (V[E[near, 0]] + V[E[near, 1]]) / 2
        return edges[near[np.argmin(self.depth(mids))]]

    def click(self, pos, shift):
        if self.mode == "object":
            self.object_click(pos, shift)
            return
        if self.mode == "point":
            hit, target = self.pick_point(pos), self.sel_v
        elif self.mode == "edge":
            hit, target = self.pick_edge(pos), self.sel_e
        else:
            hit = self.pick(pos)
            hit, target = (None if hit is None else hit[0]), self.sel
        if hit is None:
            if not shift:
                target.clear()
        elif shift:
            target ^= {hit}
        else:
            target.clear()
            target.add(hit)
        order = self._click_order()
        if not shift:
            del order[:]
        if hit in target:
            order.append(hit)
        elif hit in order:
            order.remove(hit)
        if self.mode == "point":
            self.last_v = hit if hit in self.sel_v else None
        self.last_op = None
        self.refresh()

    def _click_order(self):
        """Bu modda tiklama sirasiyla secilenler. Nokta / yuz sayisi degistiyse numaralar kaymis olabilir: silinir."""
        key = (len(self.cage.V), len(self.cage.F))
        if key != self.sel_order_key:
            self.sel_order, self.sel_order_key = {}, key
        return self.sel_order.setdefault(self.mode, [])

    def _first_selected(self, target):
        """Hala secili olan ilk tiklanan eleman; tiklanmis eleman yoksa None."""
        return next((x for x in self._click_order() if x in target), None)

    # --- islemler ---
    def push(self):
        self.undo.append((self.cage.copy(), set(self.sel), set(self.sel_v), set(self.sel_e)))
        cut = len(self.undo) - 100
        if cut > 0:                                        # en eski adimlar dusuyor: ayirma kayitlari da kayar
            del self.undo[:cut]
            self.splits = [(k - cut, n, c) for k, n, c in self.splits if k - cut > 0]
        self.redo = []
        self.dirty = True

    # --- son islem: paneldeki degerle yeniden uygulama (LightWave Numeric gibi) ---
    # Her islem bir sozlukle izlenir: kind, before (islemden onceki durum), args (yeniden yapmak icin
    # gerekenler), values (paneldeki degerler). Islem surerken self.op, bitince self.last_op olur.
    # Panelde bir deger degisince: once durumu 'before'a dondur, sonra islemi yeni degerle tekrar uygula.

    def new_op(self, kind, values, **args):
        self.op = dict(kind=kind, before=self.undo[-1], args=args, values=dict(values))
        return self.op

    def redo_last(self, name, value):
        if self.drag:
            return
        if self.grab:                                      # arac surerken sayi yazildi: araci bitir, degerle uygula
            self.end_grab(True)
        op = self.last_op
        if op is None:
            return
        if op["kind"] == "knife":                          # iki alan, tek deger: konum orani
            value = value / 100.0 if name == "Position" else value / max(op["args"]["L"], 1e-9)
            name = "t"
        op["values"][name] = value
        self._reapply(op)

    def _apply_op(self, op):
        # tasima / dondurme / olcek degerleri DUNYA eksenlerinde; kafes yerel: donusturulerek uygulanir
        k, a, v = op["kind"], op["args"], op["values"]
        if k == "box":
            self.cage = self._box_cage(v)
            self.sel, self.sel_v, self.sel_e = set(), set(), set()
        elif k == "cylinder":
            self.cage = self._cyl_cage(v, a["topo"])
            self.sel, self.sel_v, self.sel_e = set(), set(), set()
        elif k == "move":
            self.apply_move(self.move_state(a["verts"]), self.w2l_dir(np.array([v["X"], v["Y"], v["Z"]])))
        elif k == "rotate":
            self.apply_rotation(self.move_state(a["verts"]), a["c"], self.w2l_dir(np.array(a["axis"])),
                                math.radians(v["Angle"]))
        elif k == "scale":
            self.apply_scale(self.move_state(a["verts"]), a["c"], np.array([v["X"], v["Y"], v["Z"]]))
        elif k == "magnet":
            v["Radius"] = max(v["Radius"], 1e-3)
            self.magnet_radius, self.magnet_mode, self.magnet_radial = v["Radius"], int(v["Mode"]), bool(v["Radial"])
            st = self._magnet_state(a["c"], v["Radius"], int(v["Mode"]), a["mask"], a["axis"] if v["Radial"] else None)
            self._magnet_apply(st, self.w2l_dir(np.array([v["X"], v["Y"], v["Z"]])))
        elif k == "extrude":
            segs = min(max(int(round(v["Segments"])), 1), 64)
            v["Segments"] = float(segs)
            self.extrude_segs = segs
            st = self._extrude_build(op["before"], segs, a["axes"], bool(v.get("Each normal")))
            if st is not None:
                st["op"] = op
                self.extrude_set(st, v)
        elif k == "bevel":
            self.bevel(inset_only=a["inset_only"], values=(v["Inset"], v.get("Shift", 0.0)), push=False,
                       each=bool(v.get("Individual")))
        elif k == "ebevel":
            segs = min(max(int(round(v["Segments"])), 1), 32)
            v["Segments"] = float(segs)
            self.ebevel_segs = segs
            if v["Width"] > 0:
                self.cage.bevel_edges(a["edges"], v["Width"], segs)
            self.sel, self.sel_v, self.sel_e = set(), set(), set()
        elif k == "mirror":
            self._mirror(a["axis"], v["Center"], a["faces"])
        elif k == "bridge":
            err = self._bridge(v)
            if err:
                self.info = err
        elif k == "reduce":
            self._reduce(v, a["verts"])
        elif k == "merge":
            self._merge(v, a["verts"])
        elif k == "hole":
            err = self._hole(v, a.get("cap", False))
            if err:
                self.info = err
        elif k == "shell":
            err = self._shell(v, a["faces"])
            if err:
                self.info = err
        elif k in ("fromface", "fromedge", "frompoint"):
            c = self._from_cage(k, a, v)
            if isinstance(c, str):
                self.info = c
                return
            c.V = self.w2l(c.V)
            n0 = len(self.cage.V) if a["append"] else 0
            if a["append"]:
                new = self.cage.append(c)
            else:
                self.cage, new = c, list(range(len(c.F)))
            if k == "fromedge" and a.get("src"):           # kaynak degisince bagli noktalar onu izler
                self.cage.links.append(link_new(a["src"], v, self.cage, list(range(n0, n0 + len(c.V))), c.counts))
            self.sel, self.sel_v, self.sel_e = set(new) if self.mode == "face" else set(), set(), set()
            self.info = T("From face: %d polygons, %d points") % (len(c.F), len(c.V))
        elif k == "plane":
            n = np.array([v["Nx"], v["Ny"], v["Nz"]], float)
            if np.linalg.norm(n) > 1e-9:
                a["side"] = self._plane_cut(n / np.linalg.norm(n), v["Offset"], a["verts"],
                                           _knife_opts["split"] or _knife_opts["objects"])
            self.sel, self.sel_v, self.sel_e = set(), set(), set()   # clear_selection son islemi de silerdi
        elif k == "mplane":
            n = np.array([v["Nx"], v["Ny"], v["Nz"]], float)
            if np.linalg.norm(n) > 1e-9:
                ln = np.linalg.norm(n)
                if not self._mirror_line(n / ln, v["Offset"] / ln, a["src"]):
                    self.info = T("Mirror: no polygons on the source side (Flip side)")
            self.sel, self.sel_v, self.sel_e = set(), set(), set()
        elif k == "array":
            err = self._array(v, a["faces"])
            if err:
                self.info = err
        elif k == "twist":
            err = self._twist(v, a["faces"])
            if err:
                self.info = err
        elif k == "circle":
            err = self._circle(v, a)
            if err:
                self.info = err
        elif k == "knife":
            items, maxt = self._knife_cut(*a["edge"])
            self.slide_set(items, min(max(v["t"], 0.0), maxt))
        elif k == "addpoints":
            n = min(max(int(round(v["Count"])), 1), 32)
            v["Count"] = float(n)
            self.add_count = n
            new = self.cage.split_edges(a["edges"], n)
            self.sel, self.sel_e = set(), set()
            self.sel_v = {p for _, ids in new.values() for p in ids}
            self.info = T("Add points: %d points") % len(self.sel_v)
        elif k == "addpoint":
            t = min(max(v["Position"], 0.1), 99.9) / 100.0
            v["Position"] = t * 100.0
            at = {a["edge"]: (a["a"], t)}
            if a["mirror"] is not None:
                at[a["medge"]] = (a["mirror"], t)
            new = self.cage.split_edges((), at=at)
            self.sel, self.sel_e = set(), set()
            self.sel_v = {p for _, ids in new.values() for p in ids}

    def panel_spec(self):
        """Viewport'taki arac panelinin icerigi: (baslik, [(ad, deger, birim, ondalik, geri cagri)], duzenlenebilir)."""
        if self.armed == "dim":
            return self._dim_panel()
        if self.mode == "object":
            return self.object_panel_spec()
        if self.rail is not None:
            spec = self._rail_panel()
            if spec is not None:
                return spec
        live = bool(self.grab or self.drag)
        op = self.op if live else self.last_op
        if self.armed == "line" and not live and (op is None or op["kind"] != "plane"):
            return T("Knife (line)"), [], True, self._knife_checks()   # cizmeden once secenekler
        if self.armed == "mline" and not live and (op is None or op["kind"] != "mplane"):
            return T("Mirror (line)"), [], True, self._mline_checks()
        if self.armed == "magnet" and not live and (op is None or op["kind"] != "magnet"):
            return T("Magnet"), [("Radius", self.magnet_radius, " mm", 3, self.magnet_set_radius)], True, \
                self._magnet_checks()                      # surtmeden once secenekler
        if op is not None:
            k, v, a = op["kind"], op["values"], op["args"]
            r = lambda name: (lambda val: self.redo_last(name, val))
            mm = lambda n, d=3: (n, v[n], " mm", d, r(n))
            if k == "move":
                return T("Move"), [mm("X"), mm("Y"), mm("Z")], self.drag is None
            if k == "rotate":
                return T("Rotate ") + a["letter"], [("Angle", v["Angle"], "°", 2, r("Angle"))], self.drag is None
            if k == "scale":
                return "Stretch / Size", [(n, v[n], " x", 4, r(n)) for n in "XYZ"], self.drag is None
            if k == "magnet":
                return T("Magnet"), [mm("Radius"), mm("X"), mm("Y"), mm("Z")], not live, self._magnet_checks()
            if k == "extrude":
                axes = a["axes"]
                fields = [mm(n) for n in "XYZ" if n in axes] if axes else [mm("Distance")]
                fields.append(("Segments", v["Segments"], "", 0, r("Segments")))
                checks = [(n, n in axes, lambda on, n=n: self.extrude_axis(n, on)) for n in "XYZ"]
                if a.get("face"):
                    checks.append(("Individual faces (own normals)", bool(v.get("Each normal")),
                                   self.extrude_each_normal, True))
                return "Extrude", fields, self.drag is None, checks, "Axis (none = normal)"
            if k == "bevel":
                return ("Inset" if a["inset_only"] else "Bevel"), \
                    [mm("Inset")] + ([] if a["inset_only"] else [mm("Shift")]), self.drag is None, \
                    [("Individual polygons", bool(v.get("Individual")), self.bevel_individual)]
            if k == "ebevel":
                return T("Edge bevel"), [mm("Width"), ("Segments", v["Segments"], "", 0, r("Segments"))], \
                    self.drag is None
            if k == "mirror":
                return T("Mirror ") + "XYZ"[a["axis"]], [mm("Center")], self.drag is None
            if k == "fromface":
                return T("From face"), [("Arc Segments", v["Arc Segments"], "", 0, r("Arc Segments")),
                                        mm("Max Length"), ("Exact Edge", v["Exact Edge"], "", 0, r("Exact Edge"))], \
                    self.drag is None
            if k == "fromedge":
                return T("From edges"), [("Arc Segments", v["Arc Segments"], "", 0, r("Arc Segments")),
                                         mm("Max Length")], self.drag is None
            if k == "frompoint":
                return T("From points"), [], self.drag is None
            if k == "shell":
                mode = int(v["Mode"])
                checks = [(n, mode == i, lambda on, n=n: self.shell_option(n, on)) for i, n in enumerate(self.SHELL_MODES)]
                checks += [("Fill", bool(v["Fill"]), lambda on: self.shell_option("Fill", on)),
                           ("Sharp", bool(v["Sharp"]), lambda on: self.shell_option("Sharp", on)),
                           ("Flat rim", bool(v.get("Flat rim", 1.0)), lambda on: self.shell_option("Flat rim", on))]
                return T("Shell"), [mm("Thickness"), mm("Detail filter")], self.drag is None, checks
            if k == "reduce":
                return T("Reduce points"), [("Degrees", v["Degrees"], "°", 2, r("Degrees"))], self.drag is None, \
                    [("Delete points", bool(v["Delete points"]), self.reduce_option)]
            if k == "merge":
                return T("Merge points"), [mm("Distance")] if v["Fixed"] else [], self.drag is None, \
                    [("Fixed distance", bool(v["Fixed"]), self.merge_option)]
            if k == "hole":
                return T("Circle on face" if a.get("cap") else "Hole"), [mm("Radius"), ("Segments", v["Segments"], "", 0, r("Segments")), mm("Support")], \
                    self.drag is None
            if k == "bridge":
                return T("Bridge"), [(n, v[n], "", d, r(n)) for n, d in zip(self.BRIDGE, (0, 0, 2, 2))], \
                    self.drag is None, [("Auto direction", bool(v.get("Auto direction", 1.0)), self.bridge_option)]
            if k == "box":
                return T("Box"), [mm(n) for n in self.BOX_SIZE] + \
                    [(n, v[n], "", 0, r(n)) for n in self.BOX_SEGS], self.drag is None
            if k == "cylinder":
                return T("Topo cylinder" if a["topo"] else "Cylinder"), [mm("Diameter"), mm("Height")] +                     [(n, v[n], "", 0, r(n)) for n in ("Segments", "Height Segments")], self.drag is None
            if k == "plane":
                return T("Knife (plane: N.p = Offset)"), [("Normal X", v["Nx"], "", 4, r("Nx")),
                                                       ("Normal Y", v["Ny"], "", 4, r("Ny")),
                                                       ("Normal Z", v["Nz"], "", 4, r("Nz")),
                                                       mm("Offset")], self.drag is None, self._knife_checks()
            if k == "mplane":
                return T("Mirror (plane: N.p = Offset)"), [("Normal X", v["Nx"], "", 4, r("Nx")),
                                                         ("Normal Y", v["Ny"], "", 4, r("Ny")),
                                                         ("Normal Z", v["Nz"], "", 4, r("Nz")),
                                                         mm("Offset")], self.drag is None, self._mline_checks()
            if k == "array":
                return self._array_panel(v)
            if k == "twist":
                return self._twist_panel(v)
            if k == "circle":
                return self._circle_panel(v)
            if k == "addpoints":
                return T("Add points"), [("Count", v["Count"], "", 0, r("Count"))], self.drag is None
            if k == "addpoint":
                return T("Add points"), [("Position", v["Position"], " %", 2, r("Position"))], self.drag is None
            if k == "knife":
                t, L = v["t"], a["L"]
                return "Knife", [("Position", t * 100, " %", 2, r("Position")),
                                 ("Distance", t * L, " mm", 3, r("Distance"))], self.drag is None
        spec = self._group_panel()
        if spec is not None:
            return spec
        verts = self.selected_verts()
        if verts:
            P = self.l2w(self.cage.V[verts])               # dunya koordinatinda
            c, s = P.mean(0), P.max(0) - P.min(0)
            return T("Selection"), [(T("Center ") + n, c[j], " mm", 3, lambda val, j=j: self.set_center(j, val))
                             for j, n in enumerate("XYZ")] + \
                            [(T("Size ") + n, s[j], " mm", 3, lambda val, j=j: self.set_size(j, val))
                             for j, n in enumerate("XYZ")], True
        return None

    def set_center(self, axis, value):
        """Secimi, merkezi (noktalarin ortalamasi, dunya koordinati) verilen degere gelecek sekilde
        dunya ekseni boyunca tasir."""
        verts = self.selected_verts()
        if not verts or self.drag or self.grab or self._locked_block(verts):
            return
        d = np.zeros(3)
        d[axis] = value - self.l2w(self.cage.V[verts]).mean(0)[axis]
        if abs(d[axis]) < 1e-9:
            return
        self.push()
        self.last_op = None
        self.apply_move(self.move_state(verts), self.w2l_dir(d))
        self.refresh()

    def set_size(self, axis, value):
        """Secimin o (dunya) eksendeki boyunu, kendi kutusunun ortasina gore olcekleyerek verilen degere getirir."""
        verts = self.selected_verts()
        if not verts or self.drag or self.grab or value <= 0 or self._locked_block(verts):
            return
        P = self.l2w(self.cage.V[verts])
        lo, hi = P.min(0), P.max(0)
        cur = hi[axis] - lo[axis]
        if cur < 1e-9 or abs(cur - value) < 1e-9:
            return
        k = np.ones(3)
        k[axis] = value / cur
        self.push()
        self.last_op = None
        self.apply_scale(self.move_state(verts), self.w2l((lo + hi) / 2), k)
        self.refresh()

    def selected_edges(self):
        """Kenar modunda secili kenarlar; yuz modunda secili yuzlerin kenarlari; nokta modunda iki ucu secili kenarlar."""
        if self.mode == "edge":
            return set(self.sel_e)
        if self.mode == "face":
            return {_edge(f[i], f[(i + 1) % len(f)]) for f in (self.cage.F[fi] for fi in self.sel)
                    for i in range(len(f))}
        return {e for e in self.cage.edges() if e <= self.sel_v}

    def reduce_points(self):
        """Modify > Reduce points (LightWave): secili noktalar / kenarlar / cokgenler (secim yoksa butun kafes)
        uzerinde, kenar boyunca neredeyse duz giden ara noktalari cikarir. Panel: Degrees (sapma siniri),
        Delete points (bosta kalan noktalar silinsin; kapaliysa kafeste tek nokta olarak kalir)."""
        if self.mode == "object" or self.drag or self.grab:
            return
        verts = set(self.selected_verts()) or set(range(len(self.cage.V)))
        if self.sym:
            M = self.mirror_map()
            verts |= {M[v] for v in list(verts) if M.get(v) is not None}
        self.drop_tool()
        self.push()
        op = self.new_op("reduce", {"Degrees": _last_reduce["Degrees"], "Delete points": _last_reduce["Delete points"]},
                         verts=sorted(verts))
        self._apply_op(op)
        self.last_op, self.op = op, None
        self.refresh()

    def _reduce(self, v, verts):
        v["Degrees"] = min(max(float(v["Degrees"]), 0.0), 90.0)
        _last_reduce.update(v)
        n = self.cage.reduce_points(verts, v["Degrees"])
        if n and v["Delete points"]:
            self.cage.delete(set())                        # bosta kalan noktalar gider (indeksler kayar)
        self.sel, self.sel_v, self.sel_e = (set(), set(), set()) if n else (self.sel, self.sel_v, self.sel_e)
        self.info = T("Reduce points: %d removed") % n
        self.dirty = True

    def reduce_option(self, on):
        op = self.last_op
        if op is None or op["kind"] != "reduce":
            return
        op["values"]["Delete points"] = 1.0 if on else 0.0
        self._reapply(op)

    def merge_points(self):
        """Modify > Merge points (LightWave, M): secili noktalar / kenarlar / cokgenler (secim yoksa butun kafes)
        uzerinde cakisan noktalari birlestirir. Panel: Fixed distance kapali = Automatic (ayni yerdekiler); acik =
        Distance icindekiler."""
        if self.mode == "object" or self.drag or self.grab:
            return
        verts = set(self.selected_verts()) or set(range(len(self.cage.V)))
        if self.sym:
            M = self.mirror_map()
            verts |= {M[v] for v in list(verts) if M.get(v) is not None}
        self.drop_tool()
        self.push()
        op = self.new_op("merge", dict(_last_merge), verts=sorted(verts))
        self._apply_op(op)
        self.last_op, self.op = op, None
        self.refresh()

    def _merge(self, v, verts):
        v["Distance"] = max(float(v["Distance"]), 0.0)
        _last_merge.update(v)
        n = self.cage.merge_points(verts, v["Distance"] if v["Fixed"] else MERGE_AUTO)
        if n:
            self.sel, self.sel_v, self.sel_e = set(), set(), set()   # noktalarin numarasi kaydi
        self.info = T("Merge points: %d merged") % n
        self.dirty = True

    def merge_option(self, on):
        op = self.last_op
        if op is None or op["kind"] != "merge":
            return
        op["values"]["Fixed"] = 1.0 if on else 0.0
        self._reapply(op)

    def merge_polygons(self):
        """Polygon > Merge polygons (LightWave, Shift+Z): secili bitisik cokgenler tek cokgen olur. Bitisik degillerse
        ya da tek bir dis sinir cikmiyorsa uyarir, kafese dokunmaz. Symmetry: aynadaki esler de birlesir."""
        if self.mode == "object" or self.drag or self.grab:
            return
        sel = self.selected_faces()
        regions = self.cage.face_regions(sel)
        apart = len(regions) > 1
        if self.sym:
            mir = self.mirror_faces(sel, self.mirror_map()) - sel
            if apart:                                      # iki yarisi aynayla birbirine bagliysa tek parca
                regions = self.cage.face_regions(sel | mir)
                apart = len(regions) > 1
            elif mir:                                      # tek parca: aynasi da ayri birlesir (duzlemde kaynamaz)
                regions = [sel, mir]
        loops =[self.cage.face_outline(r) for r in regions]
        if len(sel) < 2:
            self.info = T("Merge polygons: select at least two adjacent polygons")
        elif apart:
            self.info = T("Merge polygons: the selected polygons are not adjacent")
        elif any(L is None for L in loops):
            self.info = T("Merge polygons: no single outline (hole, closed shape or mixed directions: try Flip)")
        else:
            self.drop_tool()
            self.push()
            new = self.cage.merge_faces(regions, loops)
            self.op = self.last_op = None
            self.sel, self.sel_v, self.sel_e = set(new) if self.mode == "face" else set(), set(), set()
            self.info = T("Merge polygons: %d polygons merged") % sum(len(r) for r in regions)
        self.refresh()

    # --- Add points / Connect / Subdivide (LightWave Add Points, Split, Subdivide): kenar ekleme (#3, #6) ---
    # Add points: secili kenarlar (nokta modunda iki ucu secili kenarlar) esit bolunur; secim yoksa arac acilir,
    # tiklanan kenarin o yerine nokta duser. Connect: ayni cokgendeki iki noktayi (ya da iki kenarin ortasini)
    # kenarla baglar, cokgen ikiye bolunur. Subdivide: secili cokgenler dortgenlere bolunur. Komsu cokgenler yeni
    # kenar noktalarini alir, yuzeyde bosluk kalmaz.

    def _mirror_edges(self, edges):
        M = self.mirror_map()
        alive = self.cage.edges()
        out = set(edges)
        for e in edges:
            me = frozenset(M.get(v) for v in e)
            if None not in me and me in alive:
                out.add(me)
        return out

    def add_points(self):
        if self.mode == "object" or self.drag or self.grab:
            return
        edges = self.selected_edges() if self.mode in ("edge", "point") else set()
        if not edges:
            self.toggle_tool("addpoint")
            return
        if self.sym:
            edges = self._mirror_edges(edges)
        self.drop_tool()
        self.push()
        op = self.new_op("addpoints", {"Count": float(getattr(self, "add_count", 1))}, edges=edges)
        self._apply_op(op)
        self.last_op, self.op = op, None
        self.refresh()

    def add_point_at(self, pos):
        """Add points araci: imlecin altindaki kenara, tiklanan yere nokta (ortaya yakinsa ortaya)."""
        h = self.hover_edge(pos)
        if h is None:
            return
        fi, u, w = h
        su, sw = self.project(self.cage.V[[u, w]])
        d = sw - su
        t = float(np.clip((np.array(pos[:2], float) - su).dot(d) / max(d.dot(d), 1e-12), 0.0, 1.0))
        if abs(t - 0.5) < 0.08:
            t = 0.5
        self.push()
        mirror = None
        if self.sym:
            M = self.mirror_map()
            mu, mw = M.get(u), M.get(w)
            if mu is not None and mw is not None and _edge(mu, mw) != _edge(u, w) and _edge(mu, mw) in self.cage.edges():
                mirror = mu
        op = self.new_op("addpoint", {"Position": t * 100.0}, edge=_edge(u, w), a=u, mirror=mirror,
                         medge=None if mirror is None else _edge(M[u], M[w]))
        self._apply_op(op)
        self.last_op, self.op = op, None
        self.refresh()

    def _connect_pairs(self, pts):
        """Secili noktalardan baglanacak ciftler: bir cokgende tam iki secili nokta varsa ve komsu degillerse."""
        pairs = set()
        for f in self.cage.F:
            s = [i for i, v in enumerate(f) if v in pts]
            if len(s) == 2 and s[1] - s[0] not in (1, len(f) - 1):
                pairs.add(tuple(sorted((f[s[0]], f[s[1]]))))
        return sorted(pairs)

    def connect(self):
        """Polygon > Connect (LightWave Split, Ctrl+L): nokta modu: ayni cokgendeki iki secili nokta arasina kenar;
        kenar modu: ayni cokgenin iki secili kenarinin ortalari arasina kenar."""
        if self.mode == "object" or self.drag or self.grab:
            return
        cage = self.cage
        if self.mode == "point":
            pts = set(self.sel_v)
            if self.sym:
                M = self.mirror_map()
                pts |= {M[v] for v in list(pts) if M.get(v) is not None}
            pairs = self._connect_pairs(pts)
            if not pairs:
                self.info = T("Connect: select two points of one polygon that are not neighbours")
                self.refresh()
                return
            self.drop_tool()
            self.push()
            new = cage.connect(pairs)
        elif self.mode == "edge":
            E = self._mirror_edges(self.sel_e) if self.sym else set(self.sel_e)
            pairs = []
            for f in cage.F:
                s = [e for e in (_edge(f[i], f[(i + 1) % len(f)]) for i in range(len(f))) if e in E]
                if len(s) == 2:
                    pairs.append(tuple(s))
            if not pairs:
                self.info = T("Connect: select two edges of one polygon")
                self.refresh()
                return
            self.drop_tool()
            self.push()
            mid = {e: ids[0] for e, (a, ids) in cage.split_edges({e for p in pairs for e in p}, 1).items()}
            new = cage.connect([(mid[a], mid[b]) for a, b in pairs])
        else:
            self.info = T("Connect: select points or edges")
            self.refresh()
            return
        self.op = self.last_op = None
        self.sel, self.sel_v = set(), set()
        self.sel_e = set(new) if self.mode == "edge" else set()
        if self.mode == "point":
            self.sel_v = {v for e in new for v in e}
        self.info = T("Connect: %d edges added") % len(new)
        self.refresh()

    def subdivide(self):
        """Polygon > Subdivide (LightWave, Shift+D): secili cokgenler dortgenlere bolunur (n kenarli -> n dortgen)."""
        if self.mode == "object" or self.drag or self.grab:
            return
        sel = self.selected_faces()
        if not sel:
            self.info = T("Subdivide: select polygons")
            self.refresh()
            return
        if self.sym:
            sel = self.mirror_faces(sel, self.mirror_map())
        self.drop_tool()
        self.push()
        new = self.cage.subdivide_faces(sel)
        self.op = self.last_op = None
        self.sel, self.sel_v, self.sel_e = set(new) if self.mode == "face" else set(), set(), set()
        self.info = T("Subdivide: %d polygons -> %d") % (len(sel), len(new))
        self.refresh()

    def weld(self):
        """Modify > Weld (LightWave, Ctrl+W): secili noktalar en son tiklanan secili noktada birlesir, o yerinde
        kalir. Symmetry: aynadaki esleri de hedefin aynasinda birlesir. Kutu seciminde sira yok: hedefe tiklanmali."""
        if self.mode == "object" or self.drag or self.grab:
            return
        t = self.last_v
        if self.mode != "point" or len(self.sel_v) < 2:
            self.info = T("Weld: select at least two points")
        elif t not in self.sel_v:
            self.info = T("Weld: Shift+click the target point last")
        else:
            self.drop_tool()
            self.push()
            to = {v: t for v in self.sel_v if v != t}
            if self.sym:
                M = self.mirror_map()
                mt = M.get(t)
                if mt is not None and mt not in to:
                    for v in list(to):
                        mv = M.get(v)
                        if mv is not None and mv not in to and mv not in (t, mt):
                            to[mv] = mt
            n = self.cage.merge_map(to)
            self.op = self.last_op = None
            self.sel, self.sel_v, self.sel_e = set(), set(), set()
            self.last_v = None
            self.info = T("Weld: %d points welded") % n
        self.refresh()

    def crease_selected(self):
        """Secili kenarlara panel keskinligini verir; hepsi zaten o keskinlikteyse keskinligi kaldirir."""
        edges = self.selected_edges()
        if not edges or self.drag or self.grab:
            return
        if self.sym:
            M = self.mirror_map()
            alive = self.cage.edges()
            for e in list(edges):
                me = frozenset(M.get(v) for v in e)
                if None not in me and me in alive:
                    edges.add(me)
        s = self.sharp_box.value()
        self.push()
        if s <= 0 or all(self.cage.C.get(e) == s for e in edges):
            for e in edges:
                self.cage.C.pop(e, None)
        else:
            for e in edges:
                self.cage.C[e] = s
        self.refresh()

    def make_solid(self):
        """Kaydeder ve bu SubD nesnesinden parametrik B-rep kati olusturur (zaten varsa onu gunceller)."""
        if not self.cage.F:
            return
        obj = self.save(full=False)
        solid = next((o for o in obj.InList if "Samples" in o.PropertiesList and o.Source == obj), None)
        if solid is None:
            obj.Document.openTransaction("SubD kati")
            solid = make_solid(obj)
            obj.Document.commitTransaction()
        else:
            solid.recompute()
        sh = solid.Shape
        App.Console.PrintMessage(T("SubD solid: %s, %d faces, valid=%s, closed=%s, volume=%.3f mm3\n") % (
            solid.Name, len(sh.Faces), sh.isValid(), sh.ShapeType == "Solid", getattr(sh, "Volume", 0.0)))
        self.info = T("Solid: %s (%d faces%s)") % (solid.Label, len(sh.Faces),
                                                "" if sh.ShapeType == "Solid" else T(", open shell"))
        self.refresh()
        return solid

    # --- belgeye yazma ve FreeCAD geri almasi ---
    # Duzenlerken kafes belgeye sik yazilir (otomatik kayit, F, Ctrl+S, kati, ayirma); bunlar FreeCAD'in geri alma
    # listesine GIRMEZ (eskiden her biri ayri "SubD" adimiydi: Object modunda Ctrl+Z son islemi degil rastgele bir
    # ara kaydi geri aliyordu). Edit oturumu bitince (Object moduna donus / kapatma) butun oturum TEK adim olur:
    # nesne once sessizce oturum basindaki haline (snapshot) dondurulur, sonra islem icinde son hal yazilir;
    # boylece FreeCAD'in Undo'su tam oturum oncesine doner.

    def _snapshot(self):
        """Edit oturumunun baslangic hali (FreeCAD geri almasinin donecegi yer)."""
        o = self.obj
        self.edit_start = None if not _alive(o) else (cage_from_obj(o), o.Levels, obj_rules(o), o.Mesh.copy())
        self._saved_key = cage_key(o) if _alive(o) else None

    def save(self, full=True, commit=False):
        """Kafesi belgedeki SubD nesnesine yazar (yoksa olusturur); nesne SubPatch yuzeyini uretir.
        full: bagli katilar da yeniden hesaplanir (otomatik kayitta atlanir, yavas olabilir). commit: Edit oturumu
        bitiyor, oturumdaki degisiklikler FreeCAD'in geri alma listesine tek adim olarak girer; yoksa sessiz yazar."""
        if not _alive(self.obj):
            self.obj = None
        if commit and self.mode != "object":               # duzenleme bitiyor: yerinde nokta kalmadiysa bag biter
            try:
                self.links_refresh()
            except Exception as e:
                App.Console.PrintWarning("SubD sketch link: %s\n" % e)
        doc = (App.listDocuments().get(self.doc_name) or App.ActiveDocument) if self.obj is None else self.obj.Document
        if self.obj is None:                               # yeni nesne: olusturma her zaman geri alinabilir
            doc.openTransaction("SubD")
            self.obj = doc.addObject("Mesh::FeaturePython", "SubD")
            SubDObject(self.obj)
            if App.GuiUp:
                SubDViewProvider(self.obj.ViewObject)
            self.obj.ViewObject.Visibility = False
            self._write(full)
            doc.commitTransaction()
            self._snapshot()
        elif commit and self.edit_start is not None and self._same_as_start():
            if self.dirty:                                 # ara kayit yazilmis ama kafes basa dondu: sessizce geri yaz
                self._write(full)
        elif commit and self.edit_start is not None:
            self._restore_start()
            doc.openTransaction(T("SubD edit %s") % self.obj.Label)
            self._write(full)
            doc.commitTransaction()
            self._snapshot()
        else:
            self._write(full)
        self.dirty = False
        return self.obj

    def links_refresh(self):
        """Sketch baglari (obj_links_update'in modelleyicideki esi): kaynak degistiyse bagli noktalar onu izler
        (geri alinabilir adim); hicbiri yerinde degilse bag biter."""
        if self.mode == "object" or not self.cage.links or self.drag or self.grab:
            return
        if self.obj is not None and not getattr(self.obj, "SketchLink", True):
            return
        doc = self.obj.Document if self.obj is not None else App.ActiveDocument
        c = self.cage.copy()
        changed, errs, ended = links_update(c, doc, self.w2l)
        if errs:
            self.info = errs[0]
        if changed:
            self.push()
            self.op = self.last_op = None
            self.cage = c
            self.info = T("Sketch link: the linked points followed the sketch")
        elif ended:
            self.cage.links = c.links
            self.info = T("Sketch link ended: none of the linked points is where it was")
        if changed or ended:
            self.dirty = True
        if changed or ended or errs:
            self.refresh()

    def groups_refresh(self):
        """Parametrik gruplar (obj_groups_update'in modelleyicideki esi): tutamak tasindiysa grubun noktalari onu
        izler (geri alinabilir adim); tutamak silindiyse ya da hicbir uye yerinde degilse grup biter."""
        if self.mode == "object" or not self.cage.groups or self.drag or self.grab or not _alive(self.obj):
            return
        c = self.cage.copy()
        changed, errs, ended = groups_update(c, self.obj.Document, self.obj.Placement)
        if errs:
            self.info = errs[0]
        if changed:
            self.push()
            self.op = self.last_op = None
            self.cage = c
            self.info = self.info or T("Group: the points followed the handle")
        elif ended:
            self.cage.groups = c.groups
            self.info = self.info or T("Group released: none of its points is where it was")
        if changed or ended:
            self.dirty = True
        if changed or ended or errs:
            self.refresh()

    # --- grup paneli: kilit, gruba ekle, gruptan cikar ---
    # Grup varsayilan kilitli: uyeleri tek basina tasinamaz (yalniz tutamak surer). Bir grubun uyesi secilince
    # panel acilir (Done ile kapanir): Unlock / Lock, Add to group, Remove from group. Kilit aciksa uyeler idx ile
    # izlenir (nokta sayisi degismemeli) ve serbest duzenlenir; Lock yeni yerlerinde tutamaga yeniden baglar.

    def _grp(self, name=None):
        name = name or self.grp_active
        return next((g for g in self.cage.groups if g["name"] == name), None)

    def _locked_block(self, verts):
        """Verilen noktalardan biri kilitli bir grubun uyesiyse True (ve bilgi mesaji)."""
        if not self.cage.groups:
            return False
        sets = group_member_sets(self.cage.V, self.cage.groups)
        vs = set(verts)
        for g in self.cage.groups:
            if g.get("locked", True) and vs & sets.get(g["name"], set()):
                self.info = T("Group %s is locked: Unlock it in the group panel to edit its points") % g["name"]
                self.refresh()
                return True
        return False

    def _group_expand(self, verts):
        """Kilitli bir grubun uyesi tasinacaksa grubun tamami tasinir. Dondurur: (noktalar, {grup: [(k, j)]}) --
        k: kayittaki sira, j: kafes indeksi."""
        moves = {}
        if not self.cage.groups:
            return list(verts), moves
        sets = group_member_sets(self.cage.V, self.cage.groups)
        vs = set(verts)
        for g in self.cage.groups:
            if g.get("locked", True) and set(verts) & sets.get(g["name"], set()):
                pairs = [(k, j) for k, js in enumerate(group_index_all(self.cage, g)) for j in js]
                moves[g["name"]] = pairs
                vs |= {j for _, j in pairs}
        return sorted(vs), moves

    def _group_follow(self, moves):
        """Tasinan / donen kilitli grubun tutamagini noktalarin yeni yerine getirir (iki yonlu baglanti)."""
        if not moves or not _alive(self.obj):
            return
        doc = self.obj.Document
        P4 = _mat4(self.obj.Placement)
        P4i = np.linalg.inv(P4)
        done = False
        for name, pairs in moves.items():
            g, h = self._grp(name), doc.getObject(name)
            if g is None or h is None or not pairs:
                continue
            ks, js = [k for k, _ in pairs], [j for _, j in pairs]
            old = np.array(g["pts"], float).reshape(-1, 3)[ks]
            new = self.cage.V[js]
            R, t = rigid_fit(old, new)
            M = np.eye(4)
            M[:3, :3], M[:3, 3] = R, t
            Hold = np.array(g["H"], float).reshape(4, 4)
            Hn = P4 @ M @ P4i @ Hold
            if np.allclose(Hn, Hold, atol=1e-9, rtol=0):
                continue
            ok = set_handle_placement(h, App.Placement(App.Matrix(*Hn.flatten().tolist())))
            try:
                Hg = _mat4(h.getGlobalPlacement())
            except Exception:
                Hg = _mat4(h.Placement)
            allp = np.array(g["pts"], float).reshape(-1, 3)
            allp[ks] = new
            g["pts"], g["H"] = allp.tolist(), Hg.flatten().tolist()
            if not ok:
                self.info = T("Group %s: the handle could not follow") % name
            done = True
        if done:
            self.dirty = True

    def _group_panel(self):
        gs = self.cage.groups
        if not gs:
            self.grp_active = None
            return None
        verts = set(self.selected_verts())
        hit = [g["name"] for g in gs if verts & self._grp_sets.get(g["name"], set())]
        if hit and (self.grp_dismissed is None or frozenset(verts) != self.grp_dismissed):
            if self.grp_active not in hit:
                self.grp_active = hit[0]                   # uyesine tiklanan grup acilir
            self.grp_dismissed = None
        g = self._grp()
        if g is None:
            self.grp_active = None
            return None
        locked = g.get("locked", True)
        n = len(self._grp_sets.get(g["name"], ())) if locked else len(g.get("idx") or [])
        title = "%s %s %s (%d)" % ("\U0001F512" if locked else "\U0001F513", g["name"], group_square(group_rgb(g)), n)
        buttons = [("Unlock" if locked else "Lock", self.group_toggle_lock), ("Add to group", self.group_add),
                   ("Remove from group", self.group_remove), ("Done", self.group_done)]
        return title, [], True, [], None, buttons

    def group_done(self):
        self.grp_active = None
        self.grp_dismissed = frozenset(self.selected_verts())

    def group_toggle_lock(self):
        g = self._grp()
        if g is None or self.drag or self.grab:
            return
        V = self.cage.V
        self.push()
        if g.get("locked", True):                          # kilidi ac: uyeler indeksle izlenir
            pairs = [(k, j) for k, js in enumerate(group_index_all(self.cage, g)) for j in js]
            for key in ("pts", "loc", "w"):
                g[key] = [g[key][k] for k, _ in pairs]
            g["idx"], g["nv"], g["locked"] = [j for _, j in pairs], len(V), False
            self.info = T("Group %s unlocked: its points can be edited, Lock binds them again") % g["name"]
        else:                                              # kilitle: yeni yerlerinde tutamaga yeniden bagla
            h = self.obj.Document.getObject(g["name"]) if _alive(self.obj) else None
            idx = g.get("idx") or []
            if h is None or g.get("nv") != len(V) or not idx:
                self.cage.groups = [x for x in self.cage.groups if x is not g]
                self.grp_active = None
                self.info = T("Group %s released: its points or handle changed, it cannot be locked again") % g["name"]
            else:
                try:
                    H = _mat4(h.getGlobalPlacement())
                except Exception:
                    H = _mat4(h.Placement)
                P = V[idx]
                g.update(pts=P.tolist(), loc=((self.l2w(P) - H[:3, 3]) @ H[:3, :3]).tolist(), w=[1.0] * len(idx),
                         H=H.flatten().tolist(), locked=True)
                g.pop("idx", None)
                g.pop("nv", None)
                self.info = T("Group %s locked") % g["name"]
        self.dirty = True
        self.refresh()

    def group_add(self):
        """Secili noktalar (cokgen / kenar secimi de olur) aktif gruba katilir; baska bir gruptakiler katilmaz."""
        g = self._grp()
        verts = self.selected_verts()
        if g is None or not verts or self.drag or self.grab:
            return
        V = self.cage.V
        taken = set().union(*group_member_sets(V, self.cage.groups).values())
        new = [v for v in verts if v not in taken]
        if not new:
            self.info = T("Group: the selected points already belong to a group")
            self.refresh()
            return
        self.push()
        P = V[new]
        if g.get("locked", True):
            H = np.array(g["H"], float).reshape(4, 4)      # son eslenen tutamak: yeni uyeler grupla birlikte kayar
            loc = (self.l2w(P) - H[:3, 3]) @ H[:3, :3]
        else:
            loc = np.zeros_like(P)                         # kilitlenince hesaplanir
            g["idx"] = list(g.get("idx") or []) + [int(v) for v in new]
        g["pts"] = list(g["pts"]) + P.tolist()
        g["loc"] = list(g["loc"]) + loc.tolist()
        g["w"] = list(g["w"]) + [1.0] * len(new)
        self.info = T("Group %s: %d points added") % (g["name"], len(new))
        self.dirty = True
        self.refresh()

    def group_remove(self):
        """Secili noktalardan aktif grubun uyeleri gruptan cikar (kafeste yerinde kalirlar)."""
        g = self._grp()
        verts = set(self.selected_verts())
        if g is None or not verts or self.drag or self.grab:
            return
        idx = group_index_all(self.cage, g)
        drop = [k for k, js in enumerate(idx) if verts & set(js)]
        if not drop:
            self.info = T("Group %s: none of the selected points belongs to it") % g["name"]
            self.refresh()
            return
        self.push()
        keep = [k for k in range(len(idx)) if k not in set(drop)]
        for key in ("pts", "loc", "w") + (("idx",) if g.get("locked", True) is False else ()):
            g[key] = [g[key][k] for k in keep]
        if not keep:                                       # son uye de cikti: grup biter (tutamak belgede kalir)
            self.cage.groups = [x for x in self.cage.groups if x is not g]
            self.grp_active = None
            self.info = T("Group %s released: no points left") % g["name"]
        else:
            self.info = T("Group %s: %d points removed") % (g["name"], len(drop))
        self.dirty = True
        self.refresh()

    def make_group(self):
        """Secili noktalar (cokgen / kenar secimi de olur) bir parametrik grup olur: belgede ayni konumda bir datum
        plane (SubD_PG_0001) olusur; onu tasiyinca / Attachment ile baglayinca noktalar izler."""
        verts = self.selected_verts()
        if self.mode == "object" or not verts or self.drag or self.grab:
            return
        used = set().union(*group_member_sets(self.cage.V, self.cage.groups).values()) if self.cage.groups else set()
        if any(v in used for v in verts):
            self.info = T("Group: some of the points already belong to a group")
            self.refresh()
            return
        obj = self.obj if _alive(self.obj) else self.save(full=False)
        doc = obj.Document
        names = {o.Name for o in doc.Objects}
        n = 1
        while "%s%04d" % (GROUP_PREFIX, n) in names:
            n += 1
        name = "%s%04d" % (GROUP_PREFIX, n)
        color = group_color(n)
        P = self.cage.V[verts]                              # nesne yereli
        W = self.l2w(P)
        center = W.mean(0)
        doc.openTransaction(T("SubD group"))
        faces = [self.cage.V[self.cage.F[fi]] for fi in sorted(self.sel) if fi < len(self.cage.F)] if self.mode == "face" else []
        Rl = group_axes(P, faces, self.cage.V.mean(0))     # secilenin kendi ekseni: z = dista bakan normal
        rot = obj.Placement.Rotation
        if Rl is not None:
            rot = rot.multiply(App.Placement(App.Matrix(Rl[0, 0], Rl[0, 1], Rl[0, 2], 0, Rl[1, 0], Rl[1, 1], Rl[1, 2], 0,
                                                       Rl[2, 0], Rl[2, 1], Rl[2, 2], 0, 0, 0, 0, 1)).Rotation)
        h = group_handle(doc, obj, name, App.Placement(App.Vector(*center), rot),
                         float(np.ptp(W, 0).max()) if len(W) > 1 else 10.0, color)
        H = _mat4(h.getGlobalPlacement())
        loc = (W - H[:3, 3]) @ H[:3, :3]
        self.push()
        self.cage.groups.append(group_new(name, P, loc, H, color))
        self.save(full=False)
        doc.commitTransaction()
        doc.recompute()
        self.info = T("Group %s: %d points") % (name, len(verts))
        self.refresh()

    def _restore_start(self):
        """Nesneyi sessizce oturum basindaki haline dondurur (geri alma listesine girmez)."""
        c0, lv0, r0, m0 = self.edit_start
        cage_to_obj(c0, self.obj)
        set_surface(self.obj, lv0, r0)
        m0 = m0.copy()
        m0.Placement = self.obj.Placement
        self.obj.Mesh = m0
        if hasattr(self.obj.Proxy, "_key"):                # yuzey onbellegi son ara kaydi tutuyor: yeniden uretilsin
            self.obj.Proxy._key = None

    def _write(self, full):
        cage_to_obj(self.cage, self.obj)
        set_surface(self.obj, self.levels, self.rules)
        self.obj.recompute()
        if full:
            for o in self.obj.InList:
                o.recompute()
        self._saved_key = cage_key(self.obj)

    def _same_as_start(self):
        c0, lv0, r0, _ = self.edit_start
        c = self.cage
        return (lv0 == self.levels and r0 == self.rules and c0.F == c.F and c0.C == c.C
                and c0.free_edges() == c.free_edges() and c0.V.shape == c.V.shape and np.array_equal(c0.V, c.V))

    def doc_undone(self):
        """FreeCAD'in Undo / Redo'su (Edit menusu, arac cubugu, Object modunda Ctrl+Z) belgeyi degistirdi.
        Duzenlenen nesneye dokunduysa kafes belgeden yeniden okunur (yoksa modelleyici eski kafesi tutar ve
        sonraki otomatik kayit geri almayi bozardi)."""
        if self.mode == "object":
            if not _alive(self.obj):
                self.obj = None
            self.refresh()
            return
        o = self.obj
        if not _alive(o):                                  # nesnenin olusturulmasi geri alindi: kafes de gitsin
            self.redo = []
            self.obj, self.undo, self.splits, self.dirty = None, [], [], False   # (yoksa Object moduna gecerken
            self.cage = Cage()                                                    # yeni nesne olarak yazilirdi)
            self.op = self.last_op = None
            self.hidden = []
            self.enter_object_mode()
            return
        if cage_key(o) == self._saved_key:
            return                                         # baska bir nesne / islem geri alindi
        self.cage = cage_from_obj(o)
        self.levels, self.rules = o.Levels, obj_rules(o)
        self._sync_surface_widgets()
        self.undo, self.redo, self.splits, self.op, self.last_op = [], [], [], None, None
        self.rail = None
        self.sel, self.sel_v, self.sel_e = set(), set(), set()
        self.dirty = False
        self.set_placement(o.Placement)
        self._snapshot()
        self.info = T("Document undo: cage reloaded")
        self.refresh()

    def _autosave(self):
        if not self.suspended and self.mode != "object" and self.dirty and self.drag is None and self.grab is None and len(self.cage.V):
            self.save(full=False)

    def do_undo(self, record=True):
        """Son adimi geri alir. record=False: basarisiz islemin kendi kaydini atmak (Redo listesine girmez)."""
        if record and self._born_undo():
            return
        if not self.undo:
            self.info = T("Nothing to undo in this edit (O: object mode, Ctrl+Z undoes whole edits)")
            self.refresh()
            return
        if record:
            self.redo.append(self._history_state())
        self._history_restore(self.undo.pop())
        if self._undo_splits():                            # ayrilan nesne silindi: ileri alinacak kafes eksik kalirdi
            self.redo = []
        self.refresh()

    def _born_undo(self):
        """Yeni eklenen Box / Cylinder'da kafes olusturuldugu gibiyken (ya da geri alinacak adim kalmadiysa) Ctrl+Z
        nesnenin kendisini geri alir: nesne silinir, Object moduna gecilir. Aksi halde baslangic kafesi ayni haliyle
        geri yuklenir ve hicbir sey olmamis gibi gorunurdu."""
        o = self.obj
        if self._born is None or self.mode == "object" or not _alive(o) or o.Name != self._born:
            return False
        if self.undo:
            try:
                a, b = self.undo[0][0], self.cage
                same = len(self.undo) == 1 and len(a.V) == len(b.V) and a.F == b.F and np.allclose(a.V, b.V)
            except Exception:
                same = False
            if not same:
                return False
        self._drop_empty()                                 # nesneyi belgeden siler (tek geri alma adimi)
        self._show_edited()
        self.drop_tool()
        self.dim_sel = None
        self.mode = "object"
        self.geo_switch.whichChild = -1
        self.op = self.last_op = None
        self.rail = None
        self.undo, self.redo, self.splits = [], [], []
        self._born = None
        Gui.Selection.clearSelection()
        self.info = T("Undo: the new object was removed")
        self.refresh()
        return True

    def do_redo(self):
        if not self.redo:
            self.info = T("Nothing to redo")
            self.refresh()
            return
        self.undo.append(self._history_state())
        self._history_restore(self.redo.pop())
        self.refresh()

    def step_history(self, redo=False):
        """Undo / Redo dugmeleri: tuslarla ayni yol (Edit modunda adim adim, Object modunda belge)."""
        self.key("redo" if redo else "z", True)

    def _history_state(self):
        return (self.cage.copy(), set(self.sel), set(self.sel_v), set(self.sel_e))

    def _history_restore(self, state):
        self.cage, self.sel, self.sel_v, self.sel_e = state
        self.grab = None
        self.op = self.last_op = None
        self.dirty = True                                  # kafes kaydedilmis olabilir (or. Knife ayirma): yeniden yaz
        if _alive(self.obj) and (self.cage.groups or any(h.SubDGroup == self.obj.Name for h in self._handles())):
            msgs = groups_reconcile(self.cage, self.obj)
            if msgs:
                self.info = msgs[0]
        r = self._rail()
        if r is None:                                      # kilavuzun olusturulmasi geri alindi
            self.rail = None
        else:
            self.rail["pick"] = []
            if self.rail["sel"] is not None and self.rail["sel"] >= len(r["pts"]):
                self.rail["sel"] = None

    def _handles(self):
        """Belgedeki bu SubD'nin grup tutamaklari (bos olabilir: nesne henuz yok)."""
        if not _alive(self.obj):
            return []
        return [h for h in self.obj.Document.Objects if "SubDGroup" in getattr(h, "PropertiesList", [])]

    def _keep_history(self):
        """Object moduna gecerken adimlar saklanir; nesne degismediyse geri donunce Ctrl+Z / Redo devam eder.
        Ayirma kesiminin oncesi saklanmaz: ayrilan nesne artik belgenin kendi geri almasinda."""
        if _alive(self.obj) and self.undo + self.redo:
            k = self.splits[-1][0] if self.splits else 0
            self.histories[self.obj.Name] = (self.undo[k:], self.redo, cage_key(self.obj))

    def _undo_splits(self):
        """Geri alma, "Split into two objects" kesiminin oncesine indiyse o kesimin olusturdugu nesneyi de sil
        (yoksa ayrilan parca hem kafeste hem yeni nesnede kalirdi). Silinen olduysa True."""
        removed = False
        while self.splits and len(self.undo) < self.splits[-1][0]:
            _, name, count = self.splits.pop()
            doc = self.obj.Document if _alive(self.obj) else App.ActiveDocument
            o = doc.getObject(name) if doc is not None else None
            if o is None:
                continue
            if doc.UndoCount == count and doc.UndoNames and doc.UndoNames[0] == "SubD split":
                doc.undo()                                 # olusturma islemi hala en ustte: iz birakmadan geri al
            else:
                doc.openTransaction("SubD split undo")
                doc.removeObject(name)
                doc.commitTransaction()
            removed = True
        return removed

    def clear_selection(self):
        self.sel.clear()
        self.sel_v.clear()
        self.sel_e.clear()
        self.sel_order = {}
        self.last_op = None

    def toggle_smooth(self):
        self.smooth = not self.smooth
        self.refresh()

    def set_mode(self, mode):
        if mode == "object":
            self.enter_object_mode()
            return
        if self.mode == "object":                          # Object modundan secim moduna: secili nesneye gir
            self.edit_mode = mode
            self.edit_object()
            return
        self.mode = mode
        self.last_op = None
        self.rail = None
        self.refresh()

    def toggle_ghost(self):
        """Hayalet modu: yuzeyler yari saydam, arkadaki kenar ve noktalar gorunur ve secilebilir."""
        self.ghost = not self.ghost
        for m, t in ((self.smooth_mat, 0.7), (self.cage_mat, 0.75)):
            m.transparency = t if self.ghost else 0.0
        self.refresh()

    def _ghost_tick(self):
        """Zamanlayicidan: hayalet modunda kamera ya da surukleme hareket ederken hizli harmanlama, 200 ms durunca
        ucgen siralamasi (temiz goruntu)."""
        if not self.ghost or self.mode == "object":
            return
        cam = self.view.getCameraNode()
        if cam is None:
            return
        key = (cam.orientation.getValue().getValue(), cam.position.getValue().getValue(),
               cam.height.getValue() if hasattr(cam, "height") else cam.heightAngle.getValue())
        moving = key != self._ghost_cam or self.drag is not None or self.grab is not None
        self._ghost_cam = key
        self._ghost_still = 0 if moving else self._ghost_still + 1
        T_ = coin.SoTransparencyType
        want = T_.SORTED_OBJECT_BLEND if self._ghost_still < 4 else T_.SORTED_OBJECT_SORTED_TRIANGLE_BLEND
        if self.ghost_blend.value.getValue() != want:
            self.ghost_blend.value = want

    def toggle_zebra(self):
        """Zebra seritleri: puruzsuz yuzeye yansiyan siyah / beyaz seritler. Kirik serit = kirik yuzey (G0),
        seritte kose = G1, yumusak serit = G2. Subpatch kapaliysa acilir."""
        self.zebra = not self.zebra
        self.zebra_switch.whichChild = 0 if self.zebra else -1
        self._matcap_apply()
        if self.zebra and not self.smooth:
            self.toggle_smooth()
        self._sync_surface_widgets()
        _zebra_button(self.zebra)
        self.refresh()

    def set_matcap(self, name):
        """Yuzey gorunumu: MATCAPS adi ya da "" (kapali); secim hatirlanir. Zebra analiz icin matcap'in onune gecer
        (zebra acikken matcap gizli, kapaninca geri gelir); matcap secmek zebrayi kapatir. Subpatch kapaliysa acilir."""
        self.matcap = name
        _param().SetString("Matcap", name)
        if name and self.env_on:
            self.set_env(False)
        self._matcap_apply()
        if name and self.zebra:
            self.toggle_zebra()
        if name and not self.smooth:
            self.toggle_smooth()
        self.refresh()

    def _matcap_apply(self):
        if self.matcap:
            self.matcap_tex.image.setValue(coin.SbVec2s(256, 256), 3, _matcap_image(self.matcap))
        self.matcap_switch.whichChild = 0 if self.matcap and not self.zebra and not self.env_on else -1
        self.env_switch.whichChild = 0 if self.env_on and not self.zebra else -1

    def _env_load(self):
        """env_path resmini dokuya yukler; hata mesajini dondurur (yoksa None)."""
        if self.env_loaded == self.env_path:
            return None
        try:
            self.env_px = _env_read(self.env_path)
        except (OSError, ValueError, IndexError) as e:
            return str(e)
        self._env_upload()
        self.env_loaded = self.env_path
        return None

    def _env_upload(self):
        w, h, data = _env_atlas(*self.env_px, self.env_rough)
        self.env_tex.image.setValue(coin.SbVec2s(w, h), 3, data)
        self.env_texel.value = 1.0 / h

    def cage_pick_color(self):
        """Options > Cage color: Subpatch kapaliyken kafes yuzlerinin rengi; secim hatirlanir (secili yuzler turuncu kalir)."""
        c = QtWidgets.QColorDialog.getColor(self.cage_color, Gui.getMainWindow(), T("Cage color"))
        if not c.isValid():                                # vazgecildi
            return
        self.cage_color = c
        _param().SetString("CageColor", c.name())
        self.cage_mat.diffuseColor = tuple(c.getRgbF()[:3])
        self.cage_color_btn.setStyleSheet("background: %s;" % c.name())
        self.refresh()

    def env_pick_color(self):
        """Options > Surface color: secim hatirlanir; yalniz ortam yansitmada kullanilir (Metal 1'de etkisiz)."""
        c = QtWidgets.QColorDialog.getColor(self.env_color, Gui.getMainWindow(), T("Surface color"))
        if not c.isValid():                                # vazgecildi
            return
        self.env_color = c
        _param().SetString("EnvColor", c.name())
        self.env_base.value = tuple(c.getRgbF()[:3])
        self._env_swatch()
        self.refresh()

    def _env_swatch(self):
        self.env_color_btn.setStyleSheet("background: %s;" % self.env_color.name())

    def set_env_look(self, rough=None, metal=None):
        """Options: Roughness (bulaniklik; doku yeniden kurulur) ve Metal (plastik - metal karisimi)."""
        if rough is not None and rough != self.env_rough:
            self.env_rough = rough
            _param().SetFloat("EnvRoughness", rough)
            self.env_rough_p.value = rough
            if self.env_px is not None:
                self._env_upload()
        if metal is not None and metal != self.env_metal:
            self.env_metal = metal
            _param().SetFloat("EnvMetal", metal)
            self.env_metal_p.value = metal
        self.refresh()

    def set_env(self, on):
        """Options > Environment reflection: resim yoksa once sorar. Zebra analiz icin onune gecer (kapaninca geri
        gelir); ortam acikken matcap gizli kalir. Subpatch kapaliysa acilir. Secim hatirlanir."""
        if on and not self.env_path and not self.env_pick(False):
            on = False
        if on:
            err = self._env_load()
            if err:
                self.info = T("Environment image: %s") % err
                on = False
        self.env_on = on
        _param().SetBool("EnvOn", on)
        self._matcap_apply()
        if on and self.zebra:
            self.toggle_zebra()
        if on and not self.smooth:
            self.toggle_smooth()
        self._sync_surface_widgets()
        self.refresh()

    def env_pick(self, apply=True):
        """Ortam resmi sec; secilmediyse False. apply: secildiyse hemen ac."""
        start = os.path.dirname(self.env_path) if self.env_path else ""
        path = QtWidgets.QFileDialog.getOpenFileName(Gui.getMainWindow(), T("Environment image"), start,
                                                     T(ENV_FILTER))[0]
        if not path:
            return False
        self.env_path, self.env_loaded = path, None
        _param().SetString("EnvPath", path)
        self.env_btn.setToolTip(path)
        if apply:
            self.set_env(True)
        return True

    def switch_language(self):
        """EN <-> TR: dili kaydet, modelleyiciyi kapatip (kafes belgeye yazilir) ayni nesneyle yeniden ac."""
        subd_i18n.set_language("tr" if subd_i18n.language() == "en" else "en")
        object_mode = self.mode == "object"
        self.close()
        start(self.obj, object_mode=object_mode)

    def toggle_options(self):
        """Options paneli (sagdaki dock: olculer, keskinlik, bevel, kutu, kisayollar) ac-kapa."""
        self.dock.setVisible(not self.dock.isVisible())
        self.refresh()

    def toggle_gizmo(self, mode="rotate"):
        """Gizmo modu: R = Dondur, H = Stretch / Size; ayni tusa tekrar basmak Tasi'ya doner."""
        self.drop_tool()                                   # gizmo tutamaklari knife tiklamalarinin altinda kalmasin
        self.gizmo.set_hover(None)
        self.gizmo.set_mode(mode if self.gizmo.mode != mode else "move")
        self.refresh()

    def set_gizmo(self, mode):
        self.drop_tool()
        self.gizmo.set_hover(None)
        self.gizmo.set_mode(mode)
        self.refresh()

    def toggle_sym(self):
        self.sym = not self.sym
        self.refresh()

    def select_all(self):
        if self.mode == "point":
            allv = set(range(len(self.cage.V)))
            self.sel_v = set() if self.sel_v == allv else allv
        elif self.mode == "edge":
            alle = self.cage.edges()
            self.sel_e = set() if self.sel_e == alle else alle
        else:
            allf = set(range(len(self.cage.F)))
            self.sel = set() if self.sel == allf else allf
        self.last_op = None
        self.refresh()

    def selected_verts(self):
        if self.mode == "object":
            return []
        if self.mode == "point":
            return sorted(self.sel_v)
        if self.mode == "edge":
            return sorted({v for e in self.sel_e for v in e})
        return sorted({v for fi in self.sel for v in self.cage.F[fi]})

    # --- simetri (X = 0 duzlemi) ---
    def mirror_map(self):
        """Her nokta icin X'e gore aynadaki esini dondurur (duzlemdekiler kendisi). Es, aynalanmis konuma tol
        icindeki en yakin noktadir; komsu hucrelere de bakilir, yoksa hucre sinirina dusen kucuk kaymalar
        (orn. 0.006 mm) esi kacirir."""
        V = self.cage.V
        if not len(V):
            return {}
        tol = 1e-4 * max(np.ptp(V, 0).max(), 1.0)
        cells = np.floor(V / tol).astype(np.int64)
        grid = {}
        for i, c in enumerate(map(tuple, cells)):
            grid.setdefault(c, []).append(i)
        Q = V * MIRROR
        qcells = np.floor(Q / tol).astype(np.int64)
        offs = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
        out = {}
        for i, (cx, cy, cz) in enumerate(qcells.tolist()):
            cand = [j for a, b, c in offs for j in grid.get((cx + a, cy + b, cz + c), ())]
            if not cand:
                out[i] = None
                continue
            d = np.linalg.norm(V[cand] - Q[i], axis=1)
            k = int(np.argmin(d))
            out[i] = cand[k] if d[k] <= tol else None
        return out

    def mirror_faces(self, faces, M):
        lookup = {frozenset(f): i for i, f in enumerate(self.cage.F)}
        out = set(faces)
        for fi in faces:
            mf = [M.get(v) for v in self.cage.F[fi]]
            if None not in mf and frozenset(mf) in lookup:
                out.add(lookup[frozenset(mf)])
        return out

    # --- extrude (yuz ya da acik kenar) ---
    # Yonler: X / Y / Z eksenlerinden hicbiri secili degilse NORMAL yonunde (Distance): yuzlerde secili yuzlerin
    # ortalama normali, kenarlarda yuzeyin devami (yuz duzleminde disari, noktaya gore). Eksen seciliyse yalniz o
    # DUNYA eksen(ler)inde (X / Y / Z mm). Segments: yan duvar dilim sayisi (tekerlek). Her extrude normalle
    # (eksensiz) baslar, segment sayisi hatirlanir; paneldeki kutucuklar ve degerlerle islem yeniden yapilir.

    def extrude(self, lazy=False):
        """lazy: araç düğmeden başlatıldı; çekme, fare görünüme girip ilk hareket ettiğinde başlar."""
        if self.drag or self.grab or not (self.mode == "face" and self.sel or self.mode == "edge" and self.sel_e):
            return
        self.drop_tool()
        self.push()
        axes = set()                                       # her extrude normalle baslar; eksen panelden secilir
        each = self.extrude_each and self.mode == "face"   # son secim hatirlanir (eksen her seferinde bos baslar)
        op = self.new_op("extrude", dict(Distance=0.0, X=0.0, Y=0.0, Z=0.0, Segments=float(self.extrude_segs),
                                         **{"Each normal": float(each)}), axes=axes, dw=np.zeros(3),
                         face=self.mode == "face")
        st = self._extrude_build(op["before"], self.extrude_segs, axes, each)
        if st is None:
            self.op = None
            self.do_undo(False)
            self.info = T("Extrude: select open (boundary) edges")
            self.refresh()
            return
        st.update(op=op, lazy=lazy)
        self._extrude_anchor(st)
        self.grab = st
        self.gizmo.place(None)
        self.extrude_set(st, op["values"])

    def _extrude_build(self, before, segs, axes, each=False):
        """Durumu extrude oncesine dondurur ve extrude topolojisini kurar (yeni katlar tabanda); noktalari
        extrude_set yerlestirir. Acik kenar yoksa None. each: her yuz kendi normalinde (esit kalinlik)."""
        c, s0, sv0, se0 = before
        self.cage, self.sel, self.sel_v, self.sel_e = c.copy(), set(s0), set(sv0), set(se0)
        cage = self.cage
        M = self.mirror_map() if self.sym else {}
        layers = {}
        m_axis = None
        if self.mode == "face":
            orig = {v for fi in self.sel for v in cage.F[fi]}
            region = self.mirror_faces(self.sel, M) if self.sym else set(self.sel)
            n = sum(cage.normal(fi) for fi in self.sel)
            if np.linalg.norm(n) > 1e-9:
                m_axis = n / np.linalg.norm(n)
            if each:                                       # her yuz ayri extrude, kendi normalinde (Blender
                verts, rows, dirs = [], [], []             # Extrude Individual Faces): ortak nokta her yuzde ayri
                for fi in sorted(region):
                    n = cage.normal(fi)
                    n = n / max(np.linalg.norm(n), 1e-12)
                    L = {}
                    cage.extrude({fi}, segs, L)
                    for v in sorted(L):
                        verts.append(v)
                        rows.append(L[v][1:])
                        dirs.append(n)
                dirs = np.array(dirs).reshape(-1, 3)
            else:
                # normaller birbirini goturuyorsa (or. butun yuzler) nokta basina Smooth Shift yonleri
                nv = cage.bevel_dirs(region)[0] if m_axis is None else None
                cage.extrude(region, segs, layers)
                verts = sorted(layers)
                rows = [layers[v][1:] for v in verts]
                if m_axis is None:
                    dirs = nv
                else:
                    dirs = np.tile(m_axis, (len(verts), 1))
                    if self.sym:                           # aynadaki taraf aynalanmis normalle
                        dirs[[i for i, v in enumerate(verts) if v not in orig]] *= MIRROR
        else:
            edges = set(self.sel_e)
            if self.sym:
                alive = cage.edges()
                edges |= {me for me in (frozenset(M.get(v) for v in e) for e in self.sel_e)
                          if None not in me and me in alive}
            res = cage.extrude_edges(edges, segs, layers, self.view_dir())
            if res is None:
                return None
            verts, dirs = res
            orig = {v for e in self.sel_e for v in e}
            om = [i for i, v in enumerate(verts) if v in orig]
            m = dirs[om].mean(0) if om else np.zeros(3)
            if np.linalg.norm(m) > 0.3:                    # yonler kabaca ayni: fare bu eksende olculur
                m_axis = m / np.linalg.norm(m)
            self.sel_e = {_edge(layers[a][-1], layers[b][-1]) for a, b in map(tuple, se0)
                          if a in layers and b in layers and len(layers[a]) > 1}
            rows = [layers[v][1:] for v in verts]
        side = np.array([1 if v in orig else -1 for v in verts])
        plane = [i for i, v in enumerate(verts) if M.get(v) == v]
        if plane:
            dirs = dirs.copy()
            dirs[plane, 0] = 0.0
        base = cage.V[verts].copy()
        return dict(kind="extrude", before=before, segs=segs, axes=axes, verts=verts, base=base,
                    idx=np.array(rows).T, dirs=dirs, side=side, plane=plane,
                    m_axis=m_axis, c=base[side > 0].mean(0) if (side > 0).any() else base.mean(0))

    def _extrude_anchor(self, st):
        """Fare olcusunun baslangici: secili eksenlere gore cizgi / duzlem (yerel), olcunun basladigi degerler."""
        vd, c, axes = self.view_dir(), st["c"], st["axes"]

        def facing(ax):                                    # ax eksenini iceren, bakisa en dik duzlemin normali
            n = np.cross(ax, np.cross(vd, ax))
            return n / np.linalg.norm(n) if np.linalg.norm(n) > 1e-6 else vd
        if not axes:
            n = facing(st["m_axis"]) if st["m_axis"] is not None else None
        elif len(axes) == 1:
            n = facing(self.w2l_dir(np.eye(3)["XYZ".index(next(iter(axes)))]))
        elif len(axes) == 2:
            n = self.w2l_dir(np.eye(3)["XYZ".index(({"X", "Y", "Z"} - axes).pop())])
        else:
            n = vd
        p0 = self.ray_plane(self.mouse, c, n) if n is not None else None
        if n is not None and p0 is None:                   # duzlem yandan gorunuyor: bakis duzlemi
            n = vd
            p0 = self.ray_plane(self.mouse, c, n)
        vv, W, H = _vv(self.view)
        st.update(n=n, p0=p0, m0=tuple(self.mouse), v0=dict(st["op"]["values"]),
                  wpp=vv.getWorldToScreenScale(coin.SbVec3f(*self.l2w(c)), 1.0 / float(H)))

    def extrude_move(self):
        g = self.grab
        if g.get("lazy"):                                  # düğmeden başladı: ölçü buradan başlar
            g["lazy"] = False
            self._extrude_anchor(g)
            return
        v, v0, axes = g["op"]["values"], g["v0"], g["axes"]
        step = self.step_box.value() if _ctrl_held() else 0.0
        if g["p0"] is None:                                # ortak yon yok (or. delik kenari): saga = mesafe
            d, t = None, (self.mouse[0] - g["m0"][0]) * g["wpp"] * _gain()
        else:
            p = self.ray_plane(self.mouse, g["c"], g["n"])
            if p is None:
                return
            d = p - g["p0"]
            t = float(d.dot(g["m_axis"])) if not axes and g["m_axis"] is not None else 0.0
        if not axes:
            v["Distance"] = _snap(v0["Distance"] + t, step)
        else:
            dw = self.l2w_dir(d) if d is not None else np.zeros(3)
            for j, a in enumerate("XYZ"):
                v[a] = _snap(v0[a] + dw[j], step) if a in axes else 0.0
        self.extrude_set(g, v)

    def extrude_set(self, st, v):
        """Katlari yerlestir: normal yonde Distance ya da secili dunya eksenlerinde (X, Y, Z); kat k = k / segs."""
        axes = st["axes"]
        if not axes:
            D = st["dirs"] * v["Distance"]
        else:
            d = self.w2l_dir(np.array([v[a] if a in axes else 0.0 for a in "XYZ"]))
            D = np.where(st["side"][:, None] < 0, d * MIRROR, d)
            D[st["plane"], 0] = 0.0                        # simetri duzlemindeki noktalar duzlemde kalir
        V = self.cage.V
        s = len(st["idx"])
        for k, row in enumerate(st["idx"], 1):
            V[row] = st["base"] + D * (k / s)
        o = st["side"] > 0
        st["op"]["args"]["dw"] = self.l2w_dir(D[o].mean(0)) if o.any() else np.zeros(3)   # dunya, asil taraf
        st["op"]["args"]["m"] = self.l2w_dir(st["m_axis"]) if st["m_axis"] is not None else None
        if not axes:
            self.info = T("Extrude  %.3f mm   %d segments") % (v["Distance"], s)
        else:
            self.info = T("Extrude  %s mm   %d segments") % (
                "  ".join("%s %.3f" % (a, v[a]) for a in "XYZ" if a in axes), s)
        self.refresh(fast=True)

    def extrude_wheel(self, steps):
        """Extrude surerken fare tekerlegi: segment sayisi (kenar bevel gibi); mesafe ayni kalir."""
        g = self.grab
        segs = min(max(g["segs"] + steps, 1), 64)
        if segs == g["segs"]:
            return
        self.extrude_segs = segs
        g["op"]["values"]["Segments"] = float(segs)
        st = self._extrude_build(g["before"], segs, g["axes"], bool(g["op"]["values"].get("Each normal")))
        for k in ("op", "lazy", "n", "p0", "m0", "v0", "wpp"):
            st[k] = g[k]
        self.grab = st
        self.extrude_set(st, st["op"]["values"])

    def extrude_axis(self, a, on=None):
        """Paneldeki X / Y / Z kutucugu (on = isaretli mi): eksen kilidi degisir (surerken: olcu yeniden baslar;
        bitmisse: islem yeni eksenlerle yeniden yapilir, cekilen mesafe korunur). on=None: klavyeden X / Y / Z
        (Blender gibi): yalniz o eksen; ayni tusa tekrar basmak normal yonune dondurur."""
        g = self.grab if self.grab is not None and self.grab.get("kind") == "extrude" else None
        op = g["op"] if g is not None else self.last_op
        if op is None or op["kind"] != "extrude" or self.drag:
            return
        axes, v, dw = op["args"]["axes"], op["values"], op["args"]["dw"]
        was = set(axes)
        if on is None:
            axes.clear()
            if was != {a}:
                axes.add(a)
        else:
            (axes.add if on else axes.discard)(a)
        if axes == was:
            return
        if axes:                                           # eksen kilidi tek yon: "Each normal" kalkar
            v["Each normal"] = 0.0
        for j, n in enumerate("XYZ"):                      # o anki yer degistirmeden yeni degerler
            v[n] = float(dw[j]) if n in axes else 0.0
        if not axes:                                       # eksenden normale: mesafe = yer degistirmenin normal bileseni
            m = op["args"].get("m")
            v["Distance"] = float(dw.dot(m)) if m is not None else float(np.linalg.norm(dw))
        if g is not None:
            self._extrude_anchor(g)
            self.extrude_set(g, v)
        else:
            self._reapply(op)

    def extrude_each_normal(self, on):
        """Extrude panelindeki "Each normal": isaretliyse secili her yuz kendi normalinde kayar (kosede esit
        kalinlik), degilse bolge ortalama normalde. Eksen kilidini kaldirir; mesafe korunur, secim hatirlanir."""
        g = self.grab if self.grab is not None and self.grab.get("kind") == "extrude" else None
        op = g["op"] if g is not None else self.last_op
        if op is None or op["kind"] != "extrude" or self.drag or not op["args"].get("face"):
            return
        v, axes = op["values"], op["args"]["axes"]
        if bool(v.get("Each normal")) == on and not axes:
            return
        self.extrude_each = on
        v["Each normal"] = float(on)
        if axes:                                           # eksenden normale: mesafe = yer degistirmenin normal bileseni
            m, dw = op["args"].get("m"), op["args"]["dw"]
            v["Distance"] = float(dw.dot(m)) if m is not None else float(np.linalg.norm(dw))
            axes.clear()
            for n in "XYZ":
                v[n] = 0.0
        if g is not None:
            st = self._extrude_build(g["before"], g["segs"], axes, on)
            st.update(op=op, lazy=g.get("lazy"))
            self.grab = st
            self._extrude_anchor(st)                       # olcu buradan yeniden baslar (eksen degisimi gibi)
            self.extrude_set(st, v)
        else:
            self._reapply(op)

    def _reapply(self, op):
        c, a, b, e = op["before"]
        self.cage, self.sel, self.sel_v, self.sel_e = c.copy(), set(a), set(b), set(e)
        self._apply_op(op)
        self._group_follow((op.get("args") or {}).get("grp"))
        self.dirty = True
        self.refresh()

    def bevel(self, inset_only=False, values=None, lazy=False, push=True, each=None):
        """LightWave Bevel / Smooth Shift: secili yuzleri bolge olarak extrude eder; yeni noktalar yuz
        normalleri boyunca kaydirilir (Shift) ve sinirdan iceri alinir (Inset). values=(inset, shift) mm
        verilirse dogrudan uygular; verilmezse fareyle: saga = inset, yukari = shift (I: yalniz inset).
        each (Individual): her yuz kendi basina bevel'lenir (None = son kullanilan)."""
        if not self.sel or self.mode != "face" or self.drag or self.grab:
            return
        each = self.bevel_each if each is None else bool(each)
        if push:                                           # push=False: son islem yeniden uygulaniyor
            self.drop_tool()
            self.push()
            self.new_op("bevel", dict(Inset=0.0, Shift=0.0, Individual=float(each)), inset_only=inset_only)
        st = self._bevel_build(inset_only, each)
        if values is not None:
            self.bevel_set(st, *values)
            if push:
                self.last_op, self.op = self.op, None
            self.refresh()
            return
        vv, W, H = _vv(self.view)
        c = st["base"].mean(0)
        st["wpp"] = vv.getWorldToScreenScale(coin.SbVec3f(*self.l2w(c)), 1.0 / float(H))   # 1 piksel kac mm
        st["p0"] = tuple(self.mouse)
        st["lazy"] = lazy
        self.grab = st
        self.gizmo.place(None)
        self.bevel_set(st, 0.0, 0.0)

    def edge_bevel(self, values=None, lazy=False):
        """Kenar modunda B: secili kenarlari pahla (Rounder). Fare saga = genislik (canli); panelde mm."""
        if self.mode != "edge" or not self.sel_e or self.drag or self.grab:
            return
        self.drop_tool()
        edges = set(self.sel_e)
        if self.sym:
            M, alive = self.mirror_map(), self.cage.edges()
            for e in list(edges):
                me = frozenset(M.get(v) for v in e)
                if None not in me and me in alive:
                    edges.add(me)
        self.push()
        segs = self.ebevel_segs
        op = self.new_op("ebevel", dict(Width=0.0, Segments=float(segs)), edges=edges)
        st = dict(kind="ebevel", edges=edges, before=self.undo[-1], op=op, segs=segs, w=0.0)
        if values is not None:
            w, st["segs"] = values if isinstance(values, tuple) else (values, segs)
            op["values"]["Segments"] = float(st["segs"])
            self.ebevel_set(st, w)
            self.last_op, self.op = self.op, None
            self.refresh()
            return
        vv, W, H = _vv(self.view)
        c = self.cage.V[sorted({v for e in edges for v in e})].mean(0)
        st["wpp"] = vv.getWorldToScreenScale(coin.SbVec3f(*self.l2w(c)), 1.0 / float(H))
        st["p0"] = tuple(self.mouse)
        st["lazy"] = lazy
        self.grab = st
        self.gizmo.place(None)
        self.info = T("Edge bevel: drag right, wheel = segments (Ctrl: mm step, Shift: fine)")
        self.refresh()

    def ebevel_set(self, st, w):
        """Durumu bevel oncesine dondur ve w genislik, st['segs'] segmentle yeniden pahla (canli)."""
        st["w"] = w
        self.cage = st["before"][0].copy()
        self.sel, self.sel_v, self.sel_e = set(), set(), set()
        ok = w <= 0 or self.cage.bevel_edges(st["edges"], w, st["segs"])
        if self.op is not None and self.op["kind"] == "ebevel":
            self.op["values"].update(Width=w, Segments=float(st["segs"]))
        self.info = (T("Edge bevel  %.3f mm   %d segments") % (w, st["segs"]) if ok
                     else T("Edge bevel: points on open edges are not supported"))
        self.refresh(fast=True)

    def ebevel_wheel(self, steps):
        """Kenar bevel surerken fare tekerlegi: segment sayisi (Blender gibi)."""
        g = self.grab
        g["segs"] = min(max(g["segs"] + steps, 1), 32)
        self.ebevel_segs = g["segs"]                       # sonraki bevel'ler de bununla baslasin
        self.ebevel_set(g, g["w"])

    def ebevel_move(self):
        g = self.grab
        if g.get("lazy"):
            g["p0"], g["lazy"] = tuple(self.mouse), False
            return
        step = self.step_box.value() if _ctrl_held() else 0.0
        w = (self.mouse[0] - g["p0"][0]) * g["wpp"] * _gain() * BEVEL_WEIGHT
        self.ebevel_set(g, max(_snap(w, step), 0.0))

    def _bevel_build(self, inset_only, each):
        """Secili yuzleri extrude eder, bevel durumunu dondurur. each: her yuz ayri extrude edilir ve yonleri
        yalniz kendisinden hesaplanir (komsu secili yuzlerle ortak kenarda da kendi yan duvari olur)."""
        region = self.mirror_faces(self.sel, self.mirror_map()) if self.sym else set(self.sel)
        groups = [{fi} for fi in sorted(region)] if each else [region]
        new, N, O = [], [], []
        for g in groups:
            Nv, Ov = self.cage.bevel_dirs(g)
            new += self.cage.extrude(g)
            N.append(Nv)
            O.append(Ov)
        return dict(kind="bevel", new=new, base=self.cage.V[new].copy(), n=np.vstack(N), o=np.vstack(O),
                    inset_only=inset_only)

    def bevel_individual(self, on):
        """Bevel panelindeki "Individual": surerken ya da son bevel'de hemen uygulanir, sonrakiler icin hatirlanir."""
        self.bevel_each = bool(on)
        g = self.grab if self.grab is not None and self.grab.get("kind") == "bevel" else None
        op = self.op if g is not None else self.last_op
        if op is None or op["kind"] != "bevel" or self.drag:
            return
        v = op["values"]
        v["Individual"] = float(on)
        if g is None:
            self._reapply(op)
            return
        c, a, b, e = op["before"]                          # surerken: bevel oncesine don, yeniden kur
        self.cage, self.sel, self.sel_v, self.sel_e = c.copy(), set(a), set(b), set(e)
        st = self._bevel_build(g["inset_only"], self.bevel_each)
        for k in ("wpp", "p0", "lazy", "param", "vals0"):
            if k in g:
                st[k] = g[k]
        self.grab = st
        self.bevel_set(st, v["Inset"], v["Shift"])

    def bevel_set(self, st, inset, shift):
        self.cage.V[st["new"]] = st["base"] + st["n"] * shift + st["o"] * inset
        p = st.get("param")
        self.info = "%sInset %.3f mm   %sShift %.3f mm%s" % (
            "▸ " if p == "inset" else "", inset, "▸ " if p == "shift" else "", shift,
            "" if st.get("inset_only") or self.grab is not st else "   " + T("(B: switch)"))
        if self.op is not None and self.op["kind"] == "bevel":
            self.op["values"].update(Inset=inset, Shift=shift)
        for b, val in ((self.inset_box, inset), (self.shift_box, shift)):
            b.blockSignals(True)
            b.setValue(val)
            b.blockSignals(False)
        self.refresh(fast=True)

    def bevel_move(self):
        g = self.grab
        if g.get("lazy"):                                  # duğmeden başladı: ölçü buradan başlar
            g["p0"], g["lazy"] = tuple(self.mouse), False
            return
        # eksen kilidi: fare ilk hangi yone giderse yalniz o deger degisir (saga-sola Inset, yukari-asagi Shift);
        # B'ye tekrar basinca digerine gecilir (bevel_switch). Deger kilidin basladigi andaki degerden ilerler.
        dx, dy = self.mouse[0] - g["p0"][0], self.mouse[1] - g["p0"][1]
        if g.get("param") is None:
            if max(abs(dx), abs(dy)) < 6:
                return
            g["param"] = "shift" if abs(dy) > abs(dx) and not g["inset_only"] else "inset"
        step = self.step_box.value() if _ctrl_held() else 0.0
        k = g["wpp"] * _gain() * BEVEL_WEIGHT              # direnc: fare hareketinin kucuk bir orani (Shift = ince)
        inset, shift = g.get("vals0", (0.0, 0.0))
        if g["param"] == "inset":
            inset = _snap(inset + dx * k, step)
        else:
            shift = _snap(shift + dy * k, step)
        self.bevel_set(g, inset, shift)

    def bevel_switch(self):
        """Bevel surerken B / I: aktif degeri Inset <-> Shift degistirir; oteki son degerinde kalir."""
        g = self.grab
        if g["inset_only"]:
            return
        v = self.op["values"]
        g["vals0"] = (v["Inset"], v["Shift"])
        g["p0"] = tuple(self.mouse)
        g["param"] = "shift" if g.get("param") == "inset" else "inset"
        self.bevel_set(g, *g["vals0"])

    def _plane_cut(self, n, d, verts, split=False):
        """Duzlemle kes; simetride aynalanmis duzlemle de. verts: kesilecek yuzlerin noktalari (None = hepsi).
        split: kesimden sonra iki yan kesim cizgisinden ayrilir (Knife line "Split in two"); split_plane'in
        (bir yanin yuzleri, ayrildi mi) sonucunu dondurur (yoksa None)."""
        verts = None if verts is None else set(verts)
        new = self.cage.plane_cut(n, d, verts)
        if verts is not None:
            verts |= set(new)
        nm = n * MIRROR
        if self.sym and not np.allclose(nm, n, atol=1e-9):  # duzlem X=0'a gore kendi aynasi degilse
            if verts is not None:
                M = self.mirror_map()
                verts |= {M[v] for v in list(verts) if M.get(v) is not None}
            new = self.cage.plane_cut(nm, d, verts)
            if verts is not None:
                verts |= set(new)
        return self.cage.split_plane(n, d, verts) if split else None

    def _mirror_line(self, n, d, src=None):
        """Mirror (line): n.p = d duzlemine (yerel) gore ayna. Kafes once duzlemle kesilir, kaynak olmayan yandaki
        yuzler silinir, kalan yari duzleme gore aynalanir; duzlemdeki noktalar ortak kalir (iki yari tek parca).
        Kaynak yan: src (secimin merkezi) hangi yandaysa o, yoksa duzlemin negatif yani (Front'ta yukaridan asagi
        cizilen cizgide sol); Flip ters cevirir. Yeni yuz sayisini dondurur (0: kaynak yanda yuz yok)."""
        c = self.cage
        n = np.asarray(n, float)
        if src is not None and float(n.dot(src)) - d > 0:
            n, d = -n, -d
        if _mline_opts["flip"]:
            n, d = -n, -d
        size = max(np.ptp(c.V, 0).max(), 1.0)
        tol, eps = 1e-4 * size, 1e-7 * size
        s = c.V @ n - d
        near = np.abs(s) < tol
        c.V[near] -= np.outer(s[near], n)                  # duzleme cok yakin noktalar tam duzleme
        c.plane_cut(n, d)
        s = c.V @ n - d
        far = {fi for fi, f in enumerate(c.F) if (s[f] > eps).any()}
        if len(far) == len(c.F):
            return 0
        if far:
            c.delete(far)
            s = c.V @ n - d
        keep = np.abs(s) < 10 * eps
        return len(c.mirror(set(range(len(c.F))), c.V - 2.0 * np.outer(s, n), keep))

    def _mline_checks(self):
        return [("Flip side", _mline_opts["flip"], self.mline_option)]

    def mline_option(self, on):
        """Mirror (line) Flip: sonraki aynalarda gecerli; son islem bir ayna ise hemen ona da uygulanir."""
        _mline_opts["flip"] = bool(on)
        op = self.last_op
        if op is None or op["kind"] != "mplane" or self.box is not None:
            self.refresh()
            return
        self._reapply(op)

    def _knife_checks(self):
        """Knife (line) panelindeki secenekler: (ad, isaretli, fn(on))."""
        return [("Split in two", _knife_opts["split"], lambda on: self.knife_option("split", on)),
                ("Split into two objects", _knife_opts["objects"], lambda on: self.knife_option("objects", on))]

    def knife_option(self, name, on):
        """Knife (line) secenegi: sonraki kesimlerde gecerli. Son islem bir kesimse hemen ona da uygulanir."""
        _knife_opts[name] = bool(on)
        op = self.last_op
        if op is None or op["kind"] != "plane" or self.box is not None:
            self.refresh()
            return
        self._reapply(op)
        if _knife_opts["objects"]:
            self._split_objects(op["args"].get("side"))

    def _split_objects(self, side):
        """Kesimin bir yanini (kucuk olani) yeni bir SubD nesnesine tasir; duzenlenen nesne oteki yanla kalir.
        side: split_plane'in (yuzler, ayrildi mi) sonucu."""
        if not side:
            return
        pos, closed = side
        if not closed:
            self.info = T("Knife: the cut does not go all the way around, the model cannot be split into objects")
            self.refresh()
            return
        rest = len(self.cage.F) - len(pos)
        part = set(pos) if len(pos) <= rest else set(range(len(self.cage.F))) - set(pos)
        if not part or len(part) == len(self.cage.F):
            return
        sub = self.cage.subset(part)
        self.cage.delete(part)
        obj = self.save(full=False)
        doc = obj.Document
        doc.openTransaction("SubD split")
        o = doc.addObject("Mesh::FeaturePython", "SubD")
        SubDObject(o)
        if App.GuiUp:
            SubDViewProvider(o.ViewObject)
        cage_to_obj(sub, o)
        set_surface(o, self.levels, self.rules)
        o.Placement = obj.Placement                        # kafes yerel koordinatta: ayni yerde dursun
        o.Label = obj.Label + "_Split"
        o.recompute()
        doc.commitTransaction()
        self.splits.append((len(self.undo), o.Name, doc.UndoCount))   # bu kesim geri alinirsa nesne de silinir
        self.last_op = None                                # yeni nesne olustu: panelden yeniden kesilmez
        self.sel, self.sel_v, self.sel_e = set(), set(), set()
        _later(lambda: _tree_mark(self.obj))               # agac yeni nesne icin satirlarini yeniden kurdu
        self.info = T("Knife: split into %s (%d polygons)") % (o.Label, len(sub.F))
        self.refresh()

    def line_snap_data(self, cage):
        """Snap aciksa knife cizgisinin yakalayacagi noktalar: kafes noktalari, kenar ortalari, kenarlar (ekran
        konumlariyla); kapaliysa None. Duzlem bakis yonunu icerdiginden arkada kalan noktalar da gecerli hedef."""
        if not self.snap or not len(cage.V) or not (cage.F or cage.W):
            return None
        E = np.array([sorted(e) for e in cage.edges()])
        return dict(V=cage.V.copy(), E=E, S=self.project(cage.V), M=self.project(cage.V[E].mean(1)))

    def line_snap(self, pos, data):
        """Imlec altindaki hedef: (ekran konumu, dunya noktasi, "point" / "mid" / "edge") ya da None. Nokta ve kenar
        ortasi 12 piksel, kenar 8 piksel icinde; nokta ve kenar ortasi kenara gore onceliklidir."""
        if data is None:
            return None
        p = np.array(pos, float)
        V, E, S, M = data["V"], data["E"], data["S"], data["M"]
        ds, dm = np.linalg.norm(S - p, axis=1), np.linalg.norm(M - p, axis=1)
        i, j = int(np.argmin(ds)), int(np.argmin(dm))
        if min(ds[i], dm[j]) < 12:
            if ds[i] <= dm[j]:
                return tuple(S[i]), self.l2w(V[i][None])[0], "point"
            return tuple(M[j]), self.l2w(V[E[j]].mean(0)[None])[0], "mid"
        a, b = S[E[:, 0]], S[E[:, 1]]
        ab = b - a
        t = np.clip(((p - a) * ab).sum(1) / np.maximum((ab * ab).sum(1), 1e-12), 0, 1)
        d = np.linalg.norm(a + ab * t[:, None] - p, axis=1)
        k = int(np.argmin(d))
        if d[k] < 8:
            return tuple(a[k] + t[k] * ab[k]), self.l2w((V[E[k, 0]] + t[k] * (V[E[k, 1]] - V[E[k, 0]]))[None])[0], "edge"
        return None

    def line_plane(self, p0, p1, w0=None, w1=None):
        """Ekrandaki cizgiden kesim duzlemi (iki uc ve bakis yonunden gecen): (n, d, gosterilecek p1) ya da None.
        Ctrl: cizgi acisi 15 derece adim; ekranda gosterilen cizgi de bu yuvarlanmis aciyla cizilir. w0, w1: yakalanan
        ucun tam dunya noktasi (piksele yuvarlanmaz, duzlem tam oradan gecer)."""
        p0, p1 = np.array(p0, float), np.array(p1, float)
        if _ctrl_held():
            v = p1 - p0
            ang = math.radians(_snap(math.degrees(math.atan2(v[1], v[0])), 15.0))
            p1 = p0 + np.linalg.norm(v) * np.array([math.cos(ang), math.sin(ang)])
        P0 = (np.array(tuple(self.view.getPoint(int(round(p0[0])), int(round(p0[1])))))   # numpy sayisi kabul etmiyor
              if w0 is None else w0)
        P1 = np.array(tuple(self.view.getPoint(int(round(p1[0])), int(round(p1[1]))))) if w1 is None else w1
        if self.view.getCameraType() == "Perspective":
            eye = np.array(self.view.getCameraNode().position.getValue().getValue())
            n = np.cross(P0 - eye, P1 - eye)
        else:
            n = np.cross(P1 - P0, np.array(tuple(self.view.getViewDirection())))
        ln = np.linalg.norm(n)
        if ln < 1e-12:
            return None
        n /= ln
        for e in np.eye(3) if w1 is None else ():          # dunya eksenine cok yakinsa tam eksene oturt (yakalanan uc
            if abs(n.dot(e)) > math.cos(math.radians(2.0)):   # kacmasin diye o zaman degil)
                n = e * np.sign(n.dot(e))
        nl = self.w2l_dir(n)                             # kesim kafesin yerel koordinatinda yapilir
        return nl, float(nl.dot(self.w2l(P0))), p1

    def line_live(self, b, pos):
        """Knife cizgisi surulurken (LightWave gibi canli): durumu kesim oncesine dondur, yeni duzlemle kes."""
        if b.get("tool") == "array":
            self.array_line_live(b, pos)
            return
        self._pending = False
        if b.get("cancelled"):
            return
        mline = b.get("tool") == "mline"
        if "before" not in b:                              # ilk hareket: geri alma noktasi ve kesilecek yuzler
            verts = set(self.selected_verts()) if self.mode == "face" and self.sel else None
            sv = self.selected_verts()                     # Mirror (line): secim kaynak yani belirler
            b["src"] = self.cage.V[sv].mean(0) if mline and sv else None
            if verts is not None and self.sym and not mline:
                M = self.mirror_map()
                verts |= {M[v] for v in list(verts) if M.get(v) is not None}
            self.push()
            b["before"], b["verts"] = self.undo[-1], verts
        c = b["before"][0]
        self.cage = c.copy()
        self.sel, self.sel_v, self.sel_e = set(), set(), set()
        # hedefler kesim oncesi kafesten: canli kesimin yeni noktalari cizginin ucunu kendine cekmesin. Ctrl: aci adimi
        hit = None if _ctrl_held() else self.line_snap(pos, b.get("snap"))
        self._line_mark(hit)
        pl = self.line_plane(b["p0"], hit[0] if hit else pos, b.get("w0"), hit[1] if hit else None)
        b["plane"] = pl
        if pl is not None and mline:
            n, d, p1 = pl
            b["new"] = self._mirror_line(n, d, b["src"])
            self._line_show(b["p0"], p1)
            self.info = (T("Mirror  N (%.3f, %.3f, %.3f)   Offset %.3f mm") % (n[0], n[1], n[2], d) if b["new"]
                         else T("Mirror: no polygons on the source side (Flip side)"))
        elif pl is not None:
            n, d, p1 = pl
            b["side"] = self._plane_cut(n, d, b["verts"], _knife_opts["split"] or _knife_opts["objects"])
            self._line_show(b["p0"], p1)
            self.info = T("Knife  N (%.3f, %.3f, %.3f)   Offset %.3f mm") % (n[0], n[1], n[2], d)
        self.refresh(fast=True)

    def line_finish(self, b):
        """Birakildi: canli kesim kalir ve 'son islem' olur (panelde Normal / Konum); hicbir sey kesilmediyse geri al."""
        if b.get("tool") == "array":                       # eksen verildi: arac kapanir, panel eksen alanlarini gosterir
            if not b.get("cancelled") and "old" not in b:  # surukleme yok: tiklama = bakis yonunde eksen
                self.array_line_live(b, b["p0"])
            self._line_hide()
            if not b.get("cancelled"):
                self.armed = None
                self.dirty = True
                self.info = ""
                self.refresh()
            return
        self._line_hide()
        if b.get("cancelled"):
            return
        if "before" not in b:
            self.info = ""
            self.refresh()
            return
        c = b["before"][0]
        pl = b.get("plane")
        if b.get("tool") == "mline":
            if pl is None or not b.get("new"):
                self.cage, self.sel, self.sel_v, self.sel_e = self.undo.pop()
                self.info = T("Mirror: no polygons on the source side (Flip side)") if pl is not None else ""
                self.refresh()
                return
            n, d, _ = pl
            self.last_op = dict(kind="mplane", before=b["before"], args=dict(src=b["src"]),
                                values=dict(Nx=n[0], Ny=n[1], Nz=n[2], Offset=d))
            self.info = T("Mirror: %d polygons") % b["new"]
            self.refresh()
            return
        if pl is None or (len(self.cage.F) == len(c.F) and len(self.cage.V) == len(c.V)):
            self.cage, self.sel, self.sel_v, self.sel_e = self.undo.pop()
            self.info = T("Knife: the line does not cut the model")
            self.refresh()
            return
        n, d, _ = pl
        self.last_op = dict(kind="plane", before=b["before"], args=dict(verts=b["verts"], side=b.get("side")),
                            values=dict(Nx=n[0], Ny=n[1], Nz=n[2], Offset=d))
        self.info = ""
        if _knife_opts["objects"]:
            self._split_objects(b.get("side"))
            return
        self.refresh()

    def line_cancel(self):
        """Esc: canli kesimi geri al."""
        b, self.box = self.box, None
        self._line_hide()
        if b is not None and b.get("tool") == "array":     # Array ekseni: onceki eksene don, arac acik kalir
            b["cancelled"] = True
            op = self.last_op
            if op is not None and op["kind"] in AXIS_OPS and "old" in b:
                op["values"].update(b["old"])
                self._reapply(op)
            else:
                self.refresh()
            return
        if b is not None:
            b["cancelled"] = True                          # bekleyen canli adimlar bir sey yapmasin
            if "before" in b:
                self.cage, self.sel, self.sel_v, self.sel_e = self.undo.pop()
        self.info = ""
        self.refresh()

    def _line_show(self, p0, p1):
        self._line_widget().set_line(self._line_q(p0), self._line_q(p1))

    def _line_hide(self):
        if self.line_overlay is not None:
            self.line_overlay.set_line(None, None)
        self._line_mark(None)

    def _line_widget(self):
        if self.line_overlay is None:
            self.line_overlay = _LineOverlay(self.gv.viewport())
        return self.line_overlay

    def _line_q(self, p):
        """Coin olay konumu (sol alt kokenli, aygit pikseli) -> Qt (sol ust, mantiksal)."""
        dpr = self.gv.viewport().devicePixelRatioF()
        return QtCore.QPointF(p[0] / dpr, (self.view.getSize()[1] - p[1]) / dpr)

    def _line_mark(self, hit):
        """Yakalama gostergesi: hedef varken imlec arti olur ve hedefte isaret cizilir; yoksa imlec eski haline doner."""
        if hit is not None:
            vp = self.gv.viewport()
            if self._line_cursor is None:
                self._line_cursor = vp.cursor()
                vp.setCursor(QtCore.Qt.CrossCursor)
            self._line_widget().set_mark(self._line_q(hit[0]), hit[2])
        else:
            if self._line_cursor is not None:
                cur, self._line_cursor = self._line_cursor, None
                self.gv.viewport().setCursor(cur)
            if self.line_overlay is not None:
                self.line_overlay.set_mark(None, None)

    def _line_hover(self):
        """Knife (line) acik, cizgi cekilmiyor: imlecin altindaki hedefi goster."""
        if self.armed == "dim":
            self._dim_hover()
            return
        self._pending = False
        on = self.armed in ("line", "mline", "array") and self.box is None and self.mode != "object"
        self._line_mark(self.line_snap(self.mouse, self.line_snap_data(self._line_cage())) if on else None)

    # --- kalici araclar (LightWave gibi): dugme sari yanar, kapatilana kadar her tiklama araci uygular ---
    # "knife": kenara tikla = loop cut + kaydirma; "line": cizgi cek = duzlemle kes; "loop": kenara tikla = halka sec.
    # Kapatma: ayni dugme / tus (K, J), Esc, sag tik; baska bir arac (extrude, bevel, tasima...) da kapatir.

    def toggle_tool(self, name):
        if self.drag or self.grab or self.box is not None or self.mode == "object":
            return
        self.armed = None if self.armed == name else name
        self._line_mark(None)
        self._magnet_ring(None)
        self._magnet_snap = None
        if self.armed == "magnet" and self.magnet_radius is None:
            P = self.l2w(self.cage.V) if len(self.cage.V) else np.zeros((1, 3))
            self.magnet_radius = max(float(np.linalg.norm(np.ptp(P, 0))) * 0.25, 1.0)
        self.info = {None: "",
                     "knife": T("Knife: click edges to cut (K / Esc / right click: off)"),
                     "line": T("Knife: drag cut lines (Ctrl: 15° steps; J / Esc / right click: off)"),
                     "mline": T("Mirror: drag the mirror line, the other side is replaced by the mirror "
                                "(Ctrl: 15° steps; Shift+V / Esc / right click: off)"),
                     "loop": T("Loop: click edges to select loops (Esc / right click: off)"),
                     "addpoint": T("Add points: click an edge to add a point there, near the middle it snaps to "
                                   "the middle (Esc / right click: off)"),
                     "magnet": T("Magnet: drag to pull the surface softly, a selection limits it "
                                 "(right drag: radius; Esc / right click: off)")}[self.armed]
        self.refresh()

    def drop_tool(self):
        self.rail = None                                   # baska arac basladi: kilavuz duzenleme biter
        if self.armed:
            if self.armed == "dim":
                self._dim_end()
            self.armed = None
            self.info = ""
            self.magnet_size = None
            self._line_mark(None)
            self._magnet_ring(None)

    # --- Magnet (LightWave gibi): basilan yerin cevresindeki noktalar uzakliga gore yumusakca tasinir ---
    # Etki alani yaricapli bir kure ya da bakis yonunde boydan boya silindir (Radial); secim varsa yalniz secili
    # noktalar etkilenir. Noktalar tasinir, yeni cokgen olusmaz. Sag tusla surukleme yaricapi ayarlar.

    def _magnet_checks(self):
        shape = [(n, self.magnet_radial == bool(i), lambda on, i=i: self.magnet_option("Radial", on, i))
                 for i, n in enumerate(("Spherical", "Radial"))]
        return shape + [(n, self.magnet_mode == i, lambda on, i=i: self.magnet_option("Mode", on, i))
                        for i, n in enumerate(FALLOFFS)]

    def magnet_option(self, name, on, value):
        """Etki alani / azalma egrisi (ikisi de tek secim); son Magnet islemi varsa yeni secimle yeniden yapilir."""
        if on:
            if name == "Radial":
                self.magnet_radial = bool(value)
            else:
                self.magnet_mode = value
            op = self.last_op
            if op is not None and op["kind"] == "magnet":
                op["values"][name] = float(value)
                self._reapply(op)
                return
        self.refresh()

    def magnet_set_radius(self, value):
        self.magnet_radius = max(value, 1e-3)
        self.refresh()

    def _magnet_center(self, pos):
        """Imlecin altindaki kafes noktasi; kafes yoksa imlecin kafes ortasindan gecen ekran duzlemindeki yeri."""
        hit = self.pick(pos)
        if hit is not None:
            return hit[1]
        return self.ray_plane(pos, self.cage.V.mean(0), self.view_dir())

    def _magnet_state(self, c, radius, mode, mask, axis=None):
        """Etkilenen noktalar ve etkileri: move_state + w (noktalar), mw (aynadaki esler). axis (yerel birim yon)
        verilirse uzaklik o eksene dik olculur: silindir. Simetride merkez bir taraftadir: iki esi de alan icindeyse
        yakin olan tasinir, oteki onu aynalar (iki taraf ayni kalsin)."""
        R = self.cage.V - np.asarray(c, float)
        if axis is not None:
            R = R - np.outer(R @ axis, axis)
        dist = np.linalg.norm(self.l2w_dir(R), axis=1)
        idx = np.nonzero(dist < radius)[0] if mask is None else np.array(
            [v for v in mask if dist[v] < radius], dtype=int)
        M = self.mirror_map() if self.sym else {}
        chosen = {int(v) for v in idx}
        verts = [int(v) for v in idx if not (M.get(int(v)) in chosen and M[int(v)] != v and
                                             (dist[M[int(v)]], M[int(v)]) < (dist[v], v))]
        st = self.move_state(verts)
        w = _falloff(dist[verts] / radius, mode) if verts else np.zeros(0)
        by = dict(zip((M.get(v) for v in verts), w))
        st["w"] = w
        st["mw"] = np.array([by[m] for m in st["mirrored"]], float)
        return st

    def _magnet_apply(self, st, d):
        self.apply_transform(st, lambda P: P + st["w"][:, None] * d, lambda P: P + st["mw"][:, None] * (d * MIRROR))

    def _magnet_ring(self, c):
        """Etki alani halkasi (yerel merkez) ya da None: kaldir."""
        if c is None:
            if self.line_overlay is not None:
                self.line_overlay.set_ring(None, 0.0)
            return
        cam = self.view.getCameraNode().orientation.getValue()
        right = np.array(cam.multVec(coin.SbVec3f(1, 0, 0)).getValue())
        w = self.l2w(np.asarray(c, float))
        a, b = self.project(np.array([c, self.w2l(w + right * self.magnet_radius)]))
        dpr = self.gv.viewport().devicePixelRatioF()
        self._line_widget().set_ring(self._line_q(a), float(np.linalg.norm(b - a)) / dpr)

    def _magnet_hover(self):
        self._pending = False
        on = self.armed == "magnet" and self.box is None and not self.grab and self.mode != "object"
        if not on or not len(self.cage.V):
            self._magnet_ring(None)
            self._line_mark(None)
            return
        self._magnet_snap_ready()                          # basmadan once hedef gorunsun: basinca merkez oraya oturur
        hit = self._snap_target(self.mouse)
        self._magnet_mark(hit)
        self._magnet_ring(self.w2l(hit[0]) if hit is not None else self._magnet_center(self.mouse))

    def _magnet_mark(self, hit):
        self._line_mark(None if hit is None else (self.project(self.w2l(hit[0]))[0], None,
                                                  "point" if hit[1] == "vertex" else "edge"))

    def _magnet_snap_ready(self):
        """Snap acikken hedefleri hazirlar, kamera ya da kafes degismedikce yeniden toplamaz (fare her
        kipirdadiginda cagrilir). Donus: _snap_prepare hedefleri; Snap kapaliysa None."""
        if not self.snap:
            self._snapdata = self._magnet_snap = None
            return None
        ref = self.project(np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float))
        key = (ref.tobytes(), self.cage.V.tobytes(), len(self.cage.F))
        if self._magnet_snap is None or self._magnet_snap[0] != key:
            self._snap_prepare([0])
            base = self._snapdata
            self._snap_with_cage(base)
            self._magnet_snap = (key, base, self._snapdata)
        self._snapdata = self._magnet_snap[2]              # cekme bunu kendi hedefleriyle degistirmis olabilir
        return self._magnet_snap[1]

    def _snap_with_cage(self, base, st=None):
        """Yakalama hedefleri: base (_snap_prepare) + kafesin noktalari ve kenarlari; st (_magnet_state) varsa
        onun tasidigi noktalar haric."""
        V = self.cage.V
        skip = set(st["verts"]) | set(st["mirrored"]) if st is not None else set()
        W, S = self.l2w(V), self.project(V)
        keep = [i for i in range(len(V)) if i not in skip]
        E = np.array([sorted(e) for e in self.cage.edges() if not skip & set(e)], dtype=int).reshape(-1, 2)
        self._snapdata = dict(P=np.vstack([base["P"], W[keep]]), Ps=np.vstack([base["Ps"], S[keep]]),
                              A=np.vstack([base["A"], W[E[:, 0]]]), As=np.vstack([base["As"], S[E[:, 0]]]),
                              B=np.vstack([base["B"], W[E[:, 1]]]), Bs=np.vstack([base["Bs"], S[E[:, 1]]]))

    def magnet_start(self, pos):
        if self.drag or self.grab or self.mode == "object" or not len(self.cage.V):
            return
        self.mouse = pos
        base = self._magnet_snap_ready()                   # Part / PartDesign ve kafesin kendi hedefleri
        hit = self._snap_target(pos)
        c = self.w2l(hit[0]) if hit is not None else self._magnet_center(pos)   # ilk tiklama hedefe oturur
        if c is None:
            return
        mask = self.selected_verts() or None
        n = self.view_dir()
        axis = n / np.linalg.norm(n)
        st = self._magnet_state(c, self.magnet_radius, self.magnet_mode, mask, axis if self.magnet_radial else None)
        if not len(st["verts"]):
            self.info = T("Magnet: no points within the radius")
            self.refresh()
            return
        if base:                                           # cekerken tasinan noktalar hedef olmasin
            self._snap_with_cage(base, st)
        self.push()
        self.new_op("magnet", dict(Radius=self.magnet_radius, X=0.0, Y=0.0, Z=0.0, Mode=float(self.magnet_mode),
                                   Radial=float(self.magnet_radial)), c=np.array(c, float), mask=mask, axis=axis)
        self.grab = dict(kind="magnet", state=st, c=c, n=n, p0=self.ray_plane(pos, c, n), mask=mask, axis=axis,
                         V0=self.cage.V.copy(), d=np.zeros(3), snap=base)
        self.gizmo.place(None)
        self._magnet_ring(c)
        self.refresh()

    def magnet_move(self):
        g = self.grab
        if not _left_held():                               # tus gorunum disinda birakildi
            self.magnet_end()
            return
        p = self.ray_plane(self.mouse, g["c"], g["n"])
        if p is None or g["p0"] is None:
            return
        g["d"] = p - g["p0"]
        if _ctrl_held():                                   # Ctrl: ortho, yalniz ekranda yatay ya da dikey (buyuk olan)
            cam = self.view.getCameraNode().orientation.getValue()
            axes = [self.w2l_dir(np.array(cam.multVec(coin.SbVec3f(*a)).getValue())) for a in ((1, 0, 0), (0, 1, 0))]
            g["d"] = max((ax * g["d"].dot(ax) for ax in axes), key=np.linalg.norm)
        hit = self._snap_target(self.mouse)
        if hit is not None:                                # son nokta: merkez hedefe tam oturur
            g["d"] = self.w2l(hit[0]) - g["c"]
        self._magnet_mark(hit)
        dw = self.l2w_dir(g["d"])
        self.op["values"].update(X=float(dw[0]), Y=float(dw[1]), Z=float(dw[2]))
        self.info = "Δ  X %.3f   Y %.3f   Z %.3f mm" % tuple(dw) + (
            "   " + T("Snap: %s") % T(hit[1]) if hit is not None else "")
        self._magnet_apply(g["state"], g["d"])

    def magnet_wheel(self, steps):
        """Surerken fare tekerlegi: yaricap."""
        g = self.grab
        self.magnet_radius *= 1.1 ** steps
        self.op["values"]["Radius"] = self.magnet_radius
        self.cage.V[:] = g["V0"]
        g["state"] = self._magnet_state(g["c"], self.magnet_radius, self.magnet_mode, g["mask"],
                                        g["axis"] if self.magnet_radial else None)
        if g["snap"]:
            self._snap_with_cage(g["snap"], g["state"])
        self._magnet_ring(g["c"])
        self._magnet_apply(g["state"], g["d"])

    def magnet_size_start(self, pos):
        """Sag tus basildi: surukleme yaricapi ayarlar, kipirdamadan birakmak araci kapatir."""
        self.magnet_size = dict(p0=pos, r0=self.magnet_radius, c=self._magnet_center(pos), moved=False)

    def magnet_size_move(self):
        self._pending = False
        g = self.magnet_size
        if g is None:
            return
        dx = self.mouse[0] - g["p0"][0]                    # saga buyur, sola kuculur (150 piksel: e kat)
        if not g["moved"] and abs(dx) < 4:
            return
        g["moved"] = True
        self.magnet_radius = max(g["r0"] * math.exp(dx / 150.0), 1e-3)
        if self.last_op is not None and self.last_op["kind"] == "magnet":
            self.last_op = None                            # panel son cekmeyi degil yeni yaricapi gostersin
        self.info = T("Radius %.3f mm") % self.magnet_radius
        if g["c"] is not None:
            self._magnet_ring(g["c"])
        self.refresh()

    def magnet_size_end(self):
        g, self.magnet_size = self.magnet_size, None
        if g is None:
            return
        if g["moved"]:
            self.info = ""
        else:
            self.drop_tool()
        self.refresh()

    def magnet_end(self):
        g = self.grab
        if g is None or g.get("kind") != "magnet":
            return
        self._line_mark(None)
        self.end_grab(bool(np.any(g["d"])))                # kipirdamadiysa geri alma adimi kalmasin
        self.info = ""

    def hover_edge(self, pos):
        """Imlecin altindaki yuz ve o yuzun imlece en yakin kenari: (yuz, u, v) ya da None."""
        hit = self.pick(pos)
        if hit is None:
            return None
        fi, p = hit
        f = self.cage.F[fi]
        V = self.cage.V

        def dist(a, b):
            ab = V[b] - V[a]
            t = np.clip((p - V[a]).dot(ab) / ab.dot(ab), 0, 1)
            return np.linalg.norm(p - (V[a] + t * ab))
        i = min(range(len(f)), key=lambda i: dist(f[i], f[(i + 1) % len(f)]))
        return fi, f[i], f[(i + 1) % len(f)]

    def loop_select_from_selection(self):
        """Seciliden halka (LightWave Select Loop): kenar modu: her secili kenarin halkasi; nokta modu: iki ucu secili
        kenarlarin halkalari; yuz modu: secili dortgenlerin yonunde dortgen halkasi. Bir sey bulduysa True."""
        cage = self.cage
        if self.mode == "edge":
            seeds = set(self.sel_e)
        elif self.mode == "point":
            seeds = {e for e in cage.edges() if e <= self.sel_v}
        else:
            # her secili dortgenin iki olasi halkasindan oteki secili yuzleri en cok iceren(ler)i: secili yuzler
            # bitisik olmasa da ayni seritteyse (or. bir halkadan aralikli yuzler) halka o yonde secilir
            add = set()
            for fi in sorted(self.sel):
                f = cage.F[fi] if fi < len(cage.F) else []
                if len(f) != 4 or fi in add:
                    continue
                rings = [set(cage.face_ring(fi, f[i], f[i + 1])) for i in (0, 1)]
                hits = [len(r & self.sel) for r in rings]
                best = max(hits)
                if best > 1:
                    for r, h in zip(rings, hits):
                        if h == best:
                            add |= r
            self.sel |= add
            return bool(add)
        if not seeds:
            return False
        loop = set()
        for e in seeds:
            loop |= cage.edge_loop(e)
        if self.mode == "edge":
            self.sel_e |= loop
        else:
            self.sel_v |= {v for e in loop for v in e}
        return True

    def loop_select(self, pos=None):
        """L: secim varsa halkalari seciliden genisletir (loop_select_from_selection); yoksa imlecin altindaki
        kenardan halka secer (eklenir). Kenar modu: kenar halkasi; nokta modu: o halkanin noktalari; yuz modu:
        imlece en yakin kenari kesen dortgen halkasi (knife'in kesecegi halka). pos: Select loop aracinin tiklamasi."""
        if pos is None and self.loop_select_from_selection():
            self.last_op = None
            self.refresh()
            return
        pos = self.mouse if pos is None else pos
        if self.mode == "face":
            h = self.hover_edge(pos)
            if h is None:
                return
            self.sel |= set(self.cage.face_ring(*h)) or {h[0]}
        else:
            e = self.pick_edge(pos)
            if e is None:
                h = self.hover_edge(pos)
                if h is None:
                    return
                e = _edge(h[1], h[2])
            loop = self.cage.edge_loop(e)
            if self.mode == "edge":
                self.sel_e |= loop
            else:
                self.sel_v |= {v for e in loop for v in e}
        self.last_op = None
        self.refresh()

    def visible(self, P):
        """Noktalar kameradan gorunuyor mu (kafes yuzlerinin arkasinda kalmiyor mu)."""
        V, F = self.cage.V, self.cage.F
        tri = [(f[0], f[i], f[i + 1]) for f in F for i in range(1, len(f) - 1)]
        if not tri or not len(P):
            return np.ones(len(P), bool)
        T = V[np.array(tri)]
        e1, e2 = T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]
        persp = self.view.getCameraType() == "Perspective"
        vd = self.view_dir()
        eye = self.w2l(np.array(self.view.getCameraNode().position.getValue().getValue())) if persp else None
        eps = 1e-4 * max(np.ptp(V, 0).max(), 1.0)
        out = np.ones(len(P), bool)
        for k, p in enumerate(P):
            if persp:
                d = p - eye
                dist = np.linalg.norm(d)
                d, o = d / dist, eye
            else:
                d, dist = vd, 1e5
                o = p - d * dist
            q = np.cross(d, e2)
            det = (e1 * q).sum(1)
            ok = np.abs(det) > 1e-12
            inv = 1.0 / np.where(ok, det, 1.0)
            s = o - T[:, 0]
            u = (s * q).sum(1) * inv
            r = np.cross(s, e1)
            v = (r * d).sum(1) * inv
            t = (e2 * r).sum(1) * inv
            hit = ok & (u >= -1e-9) & (v >= -1e-9) & (u + v <= 1 + 1e-9) & (t > 0) & (t < dist - 50 * eps)
            out[k] = not hit.any()
        return out

    def box_select(self, p0, p1, add=False, sub=False):
        """Ekranda cizilen kutunun icindekileri secer. Hayalet modu kapaliysa yalniz gorunenler."""
        if self.mode == "object":
            self.object_box_select(p0, p1, add, sub)
            return
        V, F = self.cage.V, self.cage.F
        if not len(V):
            return
        lo, hi = np.minimum(p0, p1), np.maximum(p0, p1)
        inside = lambda S: np.all((S >= lo) & (S <= hi), axis=1)
        vd = self.view_dir()
        lift = -vd * 1e-3 * max(np.ptp(V, 0).max(), 1.0)   # yuzeydeki noktayi kameraya biraz yaklastir
        if self.mode == "point":
            idx = np.nonzero(inside(self.project(V)))[0]
            if not self.ghost:
                idx = idx[self.visible(V[idx] + lift)]
            found, target = {int(i) for i in idx}, self.sel_v
        elif self.mode == "edge":
            edges = list(self.cage.edges())
            E = np.array([sorted(e) for e in edges])
            S = self.project(V)
            idx = np.nonzero(inside(S[E[:, 0]]) & inside(S[E[:, 1]]))[0]
            if not self.ghost and len(idx):
                idx = idx[self.visible((V[E[idx, 0]] + V[E[idx, 1]]) / 2 + lift)]
            found, target = {edges[i] for i in idx}, self.sel_e
        else:
            C = np.array([V[f].mean(0) for f in F])
            idx = np.nonzero(inside(self.project(C)))[0]
            if not self.ghost and len(idx):
                # yalniz onu tamamen kapali olanlar elenir; yuzun yonune bakilmaz: ondeki yuzler silinince
                # delikten gorunen arka yuzlerin ic tarafi da secilebilmeli. Yuz, merkezi ya da koselerine
                # yakin noktalarindan biri gorunuyorsa gorunur sayilir (delikten bakinca merkez cogu zaman kapali).
                keep = []
                for i in idx:
                    P = V[F[i]]
                    samples = np.vstack([C[i], 0.75 * P + 0.25 * C[i]]) + lift
                    keep.append(self.visible(samples).any())
                idx = idx[np.array(keep, bool)]
            found, target = {int(i) for i in idx}, self.sel
        if sub:
            target -= found
        else:
            if not add:
                target.clear()
                del self._click_order()[:]
            target |= found
        self.last_op = None
        self.refresh()

    def _knife_cut(self, fi, u, v):
        """Loop cut + simetride aynadaki halka. Kaydirma icin ([(yeni nokta, a, b)], en buyuk t) dondurur."""
        M = self.mirror_map() if self.sym else {}
        items = self.cage.loop_cut(fi, u, v)
        if not items:
            return [], 1.0
        mu, mv = M.get(u), M.get(v)
        if mu is None or mv is None:
            return items, 1.0
        if {mu, mv} != {u, v}:
            # aynadaki kenar ilk kesimle bolunmediyse onu da kes (taraf: u'nun esi)
            faces = self.cage.edge_faces().get(_edge(mu, mv))
            if faces:
                items += self.cage.loop_cut(faces[0], mu, mv)
            return items, 1.0
        # kenar simetri duzlemini kesiyor: ayni halkaya aynalanmis ikinci kesim. Ikinci halkanin noktalari ayni
        # (a, b) uclari arasinda 1 - t'de durur (b'den bakinca t); ust uste binmesinler diye t < 0.5
        first = {m: (a, b) for m, a, b in items}
        m0, a0, b0 = items[0]
        faces = self.cage.edge_faces().get(_edge(m0, b0))
        if not faces:
            return items, 1.0
        for m2, x, y in self.cage.loop_cut(faces[0], m0, b0):
            m1 = x if x in first else y
            a, b = first[m1]
            items.append((m2, b, a))
        return items, 0.49

    def knife(self, pos=None):
        """Imlecin altindaki kenardan loop cut; sonra yeni halka fareyle kayar (Blender gibi), tikla = onay,
        Esc = ortada birak. Paneldeki Konum % / Mesafe mm ile tam yerine konur."""
        pos = self.mouse if pos is None else pos
        h = self.hover_edge(pos)
        if h is None or self.grab or self.drag:
            return
        fi, u, v = h
        V = self.cage.V
        L = float(np.linalg.norm(V[v] - V[u]))
        su, sv = self.project([V[u], V[v]])
        self.push()
        items, maxt = self._knife_cut(fi, u, v)
        if not items:
            self.undo.pop()
            return
        t0 = min(0.5, maxt)
        op = self.new_op("knife", dict(t=t0), edge=(fi, u, v), L=L, maxt=maxt)
        op["before"] = self.undo[-1]
        self.clear_selection()
        self.grab = dict(kind="slide", items=items, p0=tuple(pos), sdir=sv - su, L=L, maxt=maxt, t0=t0, op=op)
        self.gizmo.place(None)
        self.slide_set(items, t0)

    def slide_set(self, items, t):
        if not items:
            return
        m, a, b = (list(x) for x in zip(*items))
        V = self.cage.V
        V[m] = V[a] + t * (V[b] - V[a])
        if self.op is not None and self.op["kind"] == "knife":
            self.op["values"]["t"] = t
            self.info = "Knife  %.1f %%   %.3f mm" % (t * 100, t * self.op["args"]["L"])
        self.refresh(fast=True)

    def slide_move(self):
        g = self.grab
        s2 = float(g["sdir"].dot(g["sdir"]))
        if s2 < 1e-9:
            return
        t = g["t0"] + float((np.array(self.mouse, float) - g["p0"]).dot(g["sdir"])) / s2
        if _ctrl_held() and g["L"] > 1e-9:                 # Ctrl: mesafe mm adimina
            t = _snap(t * g["L"], self.step_box.value()) / g["L"]
        self.slide_set(g["items"], min(max(t, 0.0), g["maxt"]))

    def delete(self):
        """Del: yuz modunda secili cokgenler; kenar modunda secili kenarlar (yuz kenariysa o yuzler de), nokta
        modunda secili noktalar (onlari kullanan cokgen ve kenarlarla birlikte)."""
        if self.drag or self.grab:
            return
        if self.sel and self.mode == "face":
            self.drop_tool()
            self.push()
            self.cage.delete(self.mirror_faces(self.sel, self.mirror_map()) if self.sym else self.sel)
            self.clear_selection()
            self.refresh()
        elif self.sel_e and self.mode == "edge" or self.sel_v and self.mode == "point":
            self.drop_tool()
            self.push()
            edges, points = set(self.sel_e), set(self.sel_v)
            if self.sym:
                M = self.mirror_map()
                points |= {M[v] for v in points if M.get(v) is not None}
                edges |= {me for me in (frozenset(M.get(v) for v in e) for e in edges)
                          if None not in me and len(me) == 2}
            self.cage.delete_parts(edges if self.mode == "edge" else (), points if self.mode == "point" else ())
            self.clear_selection()
            self.refresh()

    # --- kopyala / kes / yapistir (cokgenler) ---
    # Pano modul duzeyinde: modelleyici kapanip acilsa da, baska bir SubD nesnesine gecilse de kalir. Noktalar
    # DUNYA koordinatinda saklanir; yapistirilan parca kopyalandigi yerde durur (Placement'i farkli nesnede de).

    def selected_faces(self):
        """Yuz modunda secili yuzler; nokta / kenar modunda butun noktalari secime giren yuzler."""
        if self.mode == "face":
            return set(self.sel)
        verts = set(self.selected_verts())
        return {fi for fi, f in enumerate(self.cage.F) if verts.issuperset(f)}

    def copy_faces(self, cut=False):
        """Ctrl+C / Ctrl+X: secili cokgenleri (keskin kenarlariyla) panoya alir; cut: sonra siler. Kenar / nokta
        modunda butun noktalari secili cokgenler cokgen, geri kalan secili kenarlar (nokta modunda iki ucu secili
        kenarlar) serbest kenar, hicbirine girmeyen secili noktalar serbest nokta olarak alinir."""
        global _clipboard
        if self.mode == "object" or self.drag or self.grab:
            return
        M = self.mirror_map() if self.sym else {}
        sel = self.selected_faces()
        if self.sym:
            sel = self.mirror_faces(sel, M)
        fe = {_edge(f[i], f[(i + 1) % len(f)]) for f in (self.cage.F[fi] for fi in sel) for i in range(len(f))}
        wires, points = set(), set()
        if self.mode == "edge":
            wires = set(self.sel_e)
        elif self.mode == "point":
            pts = set(self.sel_v)
            if self.sym:
                pts |= {M[v] for v in pts if M.get(v) is not None}
            wires = {e for e in self.cage.edges() if e <= pts}
            points = pts
        if self.sym and self.mode == "edge":
            wires |= {me for me in (frozenset(M.get(v) for v in e) for e in wires) if None not in me and len(me) == 2}
        wires -= fe
        points -= {v for fi in sel for v in self.cage.F[fi]} | {v for e in wires for v in e}
        if not sel and not wires and not points:
            self.info = T("Copy: select polygons, edges or points first")
            self.refresh()
            return
        sub = self.cage.subset(sel, wires, points)
        sub.V = self.l2w(sub.V)
        _clipboard = sub
        self.info = (T("Cut: %d polygons, %d edges, %d points") if cut else
                     T("Copied: %d polygons, %d edges, %d points")) % (len(sub.F), len(sub.W), len(points))
        if cut:
            self.push()
            free = wires & self.cage.free_edges()          # kopyalanan yuz kenarlari keserken yerinde kalir
            self.cage.delete_parts(free, points, sel)
            self.clear_selection()
        self.refresh()

    def paste_faces(self):
        """Ctrl+V: panodaki cokgenleri bu kafese ayri parca olarak ekler (kopyalandigi yerde) ve secer;
        T / gizmo ile hemen tasinabilir. Object modunda yeni bir SubD nesnesi olur."""
        if _clipboard is None or not len(_clipboard.V) or self.drag or self.grab:
            return
        if self.mode == "object":
            self.object_paste()
            return
        self.push()
        c = _clipboard.copy()
        c.V = self.w2l(c.V)
        n0 = len(self.cage.V)
        new = self.cage.append(c)
        self._select_new(new, set(range(n0, len(self.cage.V))))
        self.last_op = None
        self.info = T("Pasted: %d polygons, %d edges (T: move)") % (len(new), len(c.W))
        self.refresh()

    def _select_new(self, faces, verts):
        """Yapistirilan parcayi o anki modda sec (kenar modunda serbest kenarlar da, nokta modunda serbest
        noktalar da)."""
        self.sel, self.sel_v, self.sel_e = set(), set(), set()
        if self.mode == "face":
            self.sel = set(faces)
        elif self.mode == "point":
            self.sel_v = set(verts)
        else:
            self.sel_e = {e for e in self.cage.edges() if e <= verts}

    def paste_as_object(self):
        """Panodaki cokgenlerden yeni bir SubD nesnesi (or. sisenin agzindan kapak) ve hemen kafesine gir."""
        if _clipboard is None or not len(_clipboard.V) or self.drag or self.grab:
            return
        if self.mode != "object":
            self.enter_object_mode()
        self.object_paste()

    # --- Mirror X / Y / Z (LightWave Mirror) ---
    # Secili cokgenlerin (secim yoksa butun kafesin) DUNYA ekseni duzlemine (eksen = Center, varsayilan 0) gore
    # aynasi ayni kafese eklenir; duzlemdeki noktalar kaynaklanir, yarim model tek parca olur. Center panelden.

    def mirror_polygons(self, axis):
        if self.mode == "object":
            self.object_mirror(axis)
            return
        if self.drag or self.grab or not self.cage.F:
            return
        self.drop_tool()
        faces = self.selected_faces() or set(range(len(self.cage.F)))
        self.push()
        op = self.new_op("mirror", dict(Center=0.0), axis=axis, faces=faces)
        self._apply_op(op)
        self.last_op, self.op = op, None
        self.refresh()

    def _mirror(self, axis, center, faces):
        W = self.l2w(self.cage.V)
        d = W[:, axis] - center
        tol = 1e-4 * max(np.ptp(W, 0).max(), 1.0)
        keep = np.abs(d) < tol
        W[keep, axis] = center                             # duzleme cok yakin noktalar tam duzleme
        self.cage.V[keep] = self.w2l(W[keep])
        Wm = W.copy()
        Wm[:, axis] = 2 * center - W[:, axis]
        used = {v for fi in faces for v in self.cage.F[fi]}
        new = self.cage.mirror(faces, self.w2l(Wm), keep)
        self.sel, self.sel_v, self.sel_e = set(new) if self.mode == "face" else set(), set(), set()
        self.info = T("Mirror %s: %d polygons, %d points merged") % ("XYZ"[axis], len(new),
                                                                    sum(1 for v in used if keep[v]))

    # --- Array (LightWave Array): secili cokgenlerin kopyalari AYRI PARCA olarak eklenir; asil parca yerinde kalir ---
    # Radial: eksen (nokta + yon, dunya koordinati) etrafinda Count eleman, Start..End acisi arasinda; tam daire
    # (End - Start >= 360) son elemani atlar. Eksen boyunca Offset ile helezon. Rectangular: dunya X / Y / Z
    # boyunca Count x Offset. Offset turu (ikisinde de): Between (elemanlar arasi), Total (ilk ile son arasi),
    # Multiply (aralik her adimda Factor ile carpilir). Eksen panelden (nokta + yon, X / Y / Z dugmeleri) ya da
    # ekranda cizgi cekerek verilir (yakalama acik; kisa tiklama = bakis yonunde eksen).
    ARRAY_MAX = 1000                                       # en cok kopya sayisi (donmayi onler)
    ARRAY_AXIS = ("Px", "Py", "Pz", "Dx", "Dy", "Dz")

    def array(self, mode=0):
        """Polygon > Array: mode 0 = Radial, 1 = Rectangular. Secili cokgenler (yoksa uyari) kopyalanir."""
        if self.mode == "object":
            self.info = T("Array: only in edit mode")
            self.refresh()
            return
        if self.drag or self.grab:
            return
        faces = self.selected_faces()
        if not faces:
            self.info = T("Array: select polygons")
            self.refresh()
            return
        self.drop_tool()
        self.push()
        op = self.new_op("array", self._array_defaults(faces, mode), faces=faces)
        self._apply_op(op)
        self.last_op, self.op = op, None
        self.refresh()

    def _array_defaults(self, faces, mode):
        """Baslangic degerleri: eksen dunya orijininden Z yonunde, Rectangular aralik = secimin boyutu."""
        cage = self.cage
        used = sorted({v for fi in faces for v in cage.F[fi]})
        size = np.ptp(self.l2w(cage.V[used]), 0)
        big = max(float(np.ptp(self.l2w(cage.V), 0).max()), 1.0)
        off = [float(s) if s > 1e-6 else big * 0.1 for s in size]
        return dict(Mode=float(mode), Count=6.0, Start=0.0, End=360.0, Offset=0.0, Rotate=1.0,
                    NX=3.0, NY=1.0, NZ=1.0, OX=off[0], OY=off[1], OZ=off[2], OffMode=0.0, Factor=1.2,
                    Px=0.0, Py=0.0, Pz=0.0, Dx=0.0, Dy=0.0, Dz=1.0)

    @staticmethod
    def _array_steps(n, off, kind, factor):
        """n elemanin ilkten uzakligi: 0 Between, 1 Total, 2 Multiply (aralik her adimda factor ile carpilir)."""
        if kind == 1:
            return [off * i / (n - 1) if n > 1 else 0.0 for i in range(n)]
        if kind == 2:
            out, x = [0.0], 0.0
            for i in range(1, n):
                x += off * factor ** (i - 1)
                out.append(x)
            return out
        return [off * i for i in range(n)]

    def _array(self, v, faces):
        """Kopyalari ekler; hata iletisi ya da None dondurur. Kafes 'before' durumundan baslar (bkz. _reapply)."""
        cage = self.cage
        sub = cage.subset(faces)
        if not len(sub.V):
            return T("Array: select polygons")
        W = self.l2w(sub.V)
        kind, fac = int(v["OffMode"]), v["Factor"]
        copies = []                                        # her kopya: yeni dunya noktalari
        if int(v["Mode"]) == 0:
            n = min(max(int(round(v["Count"])), 1), self.ARRAY_MAX)
            v["Count"] = float(n)
            d = np.array([v["Dx"], v["Dy"], v["Dz"]], float)
            ln = float(np.linalg.norm(d))
            if ln < 1e-9:
                return T("Array: the axis direction is zero")
            d /= ln
            p = np.array([v["Px"], v["Py"], v["Pz"]], float)
            c = W.mean(0)
            a0, a1 = v["Start"], v["End"]
            full = abs(a1 - a0) >= 360.0 - 1e-6           # tam daire: son eleman ilkinin ustune binmesin
            step = 0.0 if n == 1 else (a1 - a0) / (n if full else n - 1)
            hs = self._array_steps(n, v["Offset"], kind, fac)
            for i in range(n):
                th = a0 + i * step
                if abs(math.remainder(th, 360.0)) < 1e-9 and abs(hs[i]) < 1e-9:
                    continue                               # asil parcanin kendisi
                R = _axis_rot(d, math.radians(th))
                if v["Rotate"] >= 0.5:
                    copies.append((W - p) @ R.T + p + hs[i] * d)
                else:                                      # yalniz yer degistirir, yonu ayni kalir
                    copies.append(W + (R @ (c - p) + p - c) + hs[i] * d)
        else:
            nn = [min(max(int(round(v[k])), 1), self.ARRAY_MAX) for k in ("NX", "NY", "NZ")]
            for k, x in zip(("NX", "NY", "NZ"), nn):
                v[k] = float(x)
            if nn[0] * nn[1] * nn[2] > self.ARRAY_MAX:
                return T("Array: too many copies (max %d)") % self.ARRAY_MAX
            st = [self._array_steps(x, v[k], kind, fac) for x, k in zip(nn, ("OX", "OY", "OZ"))]
            for ix in range(nn[0]):
                for iy in range(nn[1]):
                    for iz in range(nn[2]):
                        t = np.array([st[0][ix], st[1][iy], st[2][iz]])
                        if np.abs(t).max() >= 1e-9:
                            copies.append(W + t)
        if not copies:
            return T("Array: nothing to add (raise Count or set an Offset)")
        new = []
        for Wn in copies:
            part = sub.copy()
            part.V = self.w2l(Wn)
            new += cage.append(part)
        self.sel, self.sel_v, self.sel_e = set(new) if self.mode == "face" else set(), set(), set()
        self.info = T("Array: %d copies, %d polygons") % (len(copies), len(new))
        return None

    def _array_panel(self, v):
        rect = int(v["Mode"]) == 1
        kind = int(v["OffMode"])
        r = lambda name: (lambda val: self.redo_last(name, val))
        if rect:
            fields = [("Count X", v["NX"], "", 0, r("NX")), ("Count Y", v["NY"], "", 0, r("NY")),
                      ("Count Z", v["NZ"], "", 0, r("NZ")),
                      ("Offset X", v["OX"], " mm", 3, r("OX")), ("Offset Y", v["OY"], " mm", 3, r("OY")),
                      ("Offset Z", v["OZ"], " mm", 3, r("OZ"))]
        else:
            fields = [("Count", v["Count"], "", 0, r("Count")), ("Start angle", v["Start"], "°", 2, r("Start")),
                      ("End angle", v["End"], "°", 2, r("End")),
                      ("Offset along axis", v["Offset"], " mm", 3, r("Offset"))]
        if kind == 2:
            fields.append(("Factor", v["Factor"], " x", 3, r("Factor")))
        if not rect:
            fields += [("Axis point " + n, v["P" + n.lower()], " mm", 3, r("P" + n.lower())) for n in "XYZ"]
            fields += [("Axis direction " + n, v["D" + n.lower()], "", 4, r("D" + n.lower())) for n in "XYZ"]
        checks = [("Radial", not rect, lambda on: self.redo_last("Mode", 0.0)),
                  ("Rectangular", rect, lambda on: self.redo_last("Mode", 1.0)),
                  ("Between", kind == 0, lambda on: self.redo_last("OffMode", 0.0)),
                  ("Total", kind == 1, lambda on: self.redo_last("OffMode", 1.0)),
                  ("Multiply", kind == 2, lambda on: self.redo_last("OffMode", 2.0))]
        buttons, side = [], []
        if not rect:
            checks.append(("Rotate copies", v["Rotate"] >= 0.5, lambda on: self.redo_last("Rotate", float(bool(on)))))
            buttons = [("Draw axis", self.array_pick)]
            side = [("Axis X", lambda: self.array_axis_dir((1.0, 0.0, 0.0))),   # kutucuklarin saginda alt alta
                    ("Axis Y", lambda: self.array_axis_dir((0.0, 1.0, 0.0))),
                    ("Axis Z", lambda: self.array_axis_dir((0.0, 0.0, 1.0)))]
        return T("Array"), fields, self.drag is None, checks, None, buttons, side

    def array_axis_dir(self, d):
        """Axis X / Y / Z dugmeleri: eksen yonunu dunya eksenine al (nokta ayni kalir)."""
        op = self._axis_op()
        if op is None or self.drag or self.grab:
            return
        op["values"].update(Dx=float(d[0]), Dy=float(d[1]), Dz=float(d[2]))
        self._reapply(op)

    def array_pick(self):
        """Draw axis: ekranda cizgi cekerek ekseni ver (bir kez; sonra arac kapanir). Yakalama Mirror (line) ile aynidir."""
        op = self._axis_op()
        if op is None or self.drag or self.grab or self.box is not None:
            return
        self.armed = None if self.armed == "array" else "array"
        self._line_mark(None)
        self.info = T("Array: drag the axis line, a click gives an axis along the view "
                      "(Ctrl: 15° steps; Esc / right click: off)") if self.armed else ""
        self.refresh()

    def _line_cage(self):
        """Cizgi yakalama hedeflerinin kafesi: Array ekseni cizilirken kopyalar degil, array oncesi kafes."""
        op = self.last_op
        if self.armed == "array" and op is not None and op["kind"] in AXIS_OPS:
            return op["before"][0]
        return self.cage

    def array_axis_from(self, p0, p1, w0=None, w1=None):
        """Ekrandaki cizgiden Array ekseni: (dunya noktasi, birim yon, gosterilecek p1) ya da None. Kisa cizgi
        (tiklama) eksenin bakis yonunde olmasi demektir. Ctrl: 15 derece adim. w0, w1: yakalanan ucun tam dunya
        noktasi. Dunya eksenine 2 dereceden yakin yon tam eksene oturur (kisa cizgide hep + yonde)."""
        p0, p1 = np.array(p0, float), np.array(p1, float)
        moved = abs(p1[0] - p0[0]) + abs(p1[1] - p0[1]) > 4
        if moved and _ctrl_held():
            v = p1 - p0
            ang = math.radians(_snap(math.degrees(math.atan2(v[1], v[0])), 15.0))
            p1 = p0 + np.linalg.norm(v) * np.array([math.cos(ang), math.sin(ang)])
            w1 = None
        P0 = (np.array(tuple(self.view.getPoint(int(round(p0[0])), int(round(p0[1]))))) if w0 is None
              else np.array(w0, float))
        if moved:
            P1 = (np.array(tuple(self.view.getPoint(int(round(p1[0])), int(round(p1[1]))))) if w1 is None
                  else np.array(w1, float))
            d = P1 - P0
            if np.linalg.norm(d) < 1e-9:
                return None
        else:
            d, p1 = np.array(tuple(self.view.getViewDirection()), float), p0
        d = d / np.linalg.norm(d)
        if w1 is None or not moved:
            for e in np.eye(3):
                if abs(d.dot(e)) > math.cos(math.radians(2.0)):
                    d = e * (np.sign(d.dot(e)) if moved else 1.0)
        return P0, d, p1

    def array_line_live(self, b, pos):
        """Eksen cizgisi surulurken canli: array oncesi durumdan yeni eksenle yeniden kopyala."""
        self._pending = False
        op = self.last_op
        if b.get("cancelled") or op is None or op["kind"] not in AXIS_OPS:
            return
        v = op["values"]
        if "old" not in b:                                 # Esc icin onceki eksen
            b["old"] = {k: v[k] for k in self.ARRAY_AXIS}
        hit = None if _ctrl_held() else self.line_snap(pos, b.get("snap"))
        self._line_mark(hit)
        ax = self.array_axis_from(b["p0"], hit[0] if hit else pos, b.get("w0"), hit[1] if hit else None)
        if ax is None:
            return
        P, D, p1 = ax
        v.update(Px=float(P[0]), Py=float(P[1]), Pz=float(P[2]), Dx=float(D[0]), Dy=float(D[1]), Dz=float(D[2]))
        c, a, bb, e = op["before"]
        self.cage, self.sel, self.sel_v, self.sel_e = c.copy(), set(a), set(bb), set(e)
        self._apply_op(op)
        self._line_show(b["p0"], p1)
        self.refresh(fast=True)

    # --- Twist (LightWave Twist): secili cokgenler bir eksen etrafinda burulur ---
    # Aci eksen boyunca Start..End arasinda bir egriyle (TWIST_FALLOFFS) 0'dan Angle'a cikar: Linear / Smooth / Ease in /
    # Ease out'ta bir uc sabit kalir (Reverse: oteki uc), Peak / Dome'da iki uc sabit, orta tam doner. Auto range:
    # Start / End secimin eksen boyunca uclaridir. Cuts: Start..End arasi eksene dik esit aralikli duzlemlerle secili
    # cokgenler once bolunur (burulma yumusak gorunsun diye); komsu cokgenler yeni noktalari alir, kafes kapali kalir.
    # Eksen Array ile ayni sekilde verilir (Draw axis, Axis X / Y / Z, sayi). Fare saga = saat yonu, Ctrl = adim.
    # Simetri yok sayilir. Kisayol atanabilir (set_user_key).
    TWIST_MAX_CUTS = 64

    def twist(self):
        """Modify > Twist: fareyi surukle = aci (canli), tikla = onayla, Esc = vazgec; sonra panelden sayiyla."""
        if self.mode == "object":
            self.info = T("Twist: only in edit mode")
            self.refresh()
            return
        if self.drag or self.grab:
            return
        faces = self.selected_faces()
        if not faces:
            self.info = T("Twist: select polygons")
            self.refresh()
            return
        if self._locked_block(sorted({v for fi in faces for v in self.cage.F[fi]})):
            return
        self.drop_tool()
        self.push()
        op = self.new_op("twist", self._twist_defaults(faces), faces=faces)
        self._apply_op(op)
        self.grab = dict(kind="twist", lazy=True, p0=tuple(self.mouse))
        self.gizmo.place(None)
        self.info = T("Twist: move the mouse to twist (Ctrl: steps), click to accept, Esc to cancel")
        self.refresh()

    def _twist_defaults(self, faces):
        """Baslangic: eksen secimin merkezinden, secimin en uzun oldugu dunya ekseni yonunde; aci 0."""
        cage = self.cage
        used = sorted({v for fi in faces for v in cage.F[fi]})
        W = self.l2w(cage.V[used])
        c = W.mean(0)
        D = np.eye(3)[int(np.argmax(np.ptp(W, 0)))]
        s = (W - c) @ D
        return dict(Angle=0.0, Start=float(s.min()), End=float(s.max()), Falloff=0.0, Reverse=0.0, Auto=1.0,
                    Cuts=0.0, Px=float(c[0]), Py=float(c[1]), Pz=float(c[2]),
                    Dx=float(D[0]), Dy=float(D[1]), Dz=float(D[2]))

    def _twist(self, v, faces):
        """Kafesi 'before' durumundan burar; hata iletisi ya da None dondurur."""
        cage = self.cage
        D = np.array([v["Dx"], v["Dy"], v["Dz"]], float)
        ln = float(np.linalg.norm(D))
        if ln < 1e-9:
            return T("Twist: the axis direction is zero")
        D /= ln
        P = np.array([v["Px"], v["Py"], v["Pz"]], float)
        moved = {x for fi in faces for x in cage.F[fi]}
        if not moved:
            return T("Twist: select polygons")
        s = (self.l2w(cage.V[sorted(moved)]) - P) @ D
        if v["Auto"] >= 0.5:
            v["Start"], v["End"] = float(s.min()), float(s.max())
        a, b = v["Start"], v["End"]
        span = b - a
        if span < 1e-9:
            return T("Twist: the selection has no length along the axis (set Start and End)")
        cuts = min(max(int(round(v["Cuts"])), 0), self.TWIST_MAX_CUTS)
        v["Cuts"] = float(cuts)
        if cuts:                                           # eksene dik duzlemlerle bol; yeni noktalar da burulur
            vmap = {}                                      # nokta -> ilk secili yuzler (bolunce yuz numaralari kayar)
            for fi in faces:
                for x in cage.F[fi]:
                    vmap.setdefault(x, set()).add(fi)
            nl = self.w2l_dir(D)
            base = float(nl.dot(self.w2l(P)))
            added = set()
            for i in range(1, cuts + 1):
                new = set(cage.plane_cut(nl, base + a + span * i / (cuts + 1), moved))
                added |= new
                moved |= new
            sel = set()                                    # secim: secili yuzlerden tureyen yuzler (komsular degil)
            if self.mode == "face":
                for fi, f in enumerate(cage.F):
                    old = [x for x in f if x not in added]
                    if old and set.intersection(*(vmap.get(x, set()) for x in old)):
                        sel.add(fi)
            self.sel, self.sel_v, self.sel_e = sel, set(), set()
        idx = sorted(moved)
        W = self.l2w(cage.V[idx])
        t = ((W - P) @ D - a) / span
        if v["Reverse"] >= 0.5:
            t = 1.0 - t
        th = math.radians(v["Angle"]) * _twist_weight(t, int(v["Falloff"]))
        r = W - P
        par = np.outer(r @ D, D)
        perp = r - par
        cs, sn = np.cos(th)[:, None], np.sin(th)[:, None]
        cage.V[idx] = self.w2l(P + par + perp * cs + np.cross(D, perp) * sn)
        self.info = T("Twist %.2f°   %s") % (v["Angle"], T(TWIST_FALLOFFS[int(v["Falloff"])])) + (
            "   " + T("(symmetry is ignored)") if self.sym else "")
        return None

    def twist_move(self):
        """Twist surerken fare: yatay hareket aciyi belirler (saga = saat yonu, eksene gore); Ctrl: Step acisi."""
        g = self.grab
        op = self.op
        if g.get("lazy"):                                  # dugmeden basladi: olcu buradan baslar
            g["p0"], g["lazy"] = tuple(self.mouse), False
            return
        if op is None or op["kind"] != "twist":
            return
        v = op["values"]
        D = np.array([v["Dx"], v["Dy"], v["Dz"]], float)
        toward = float(D.dot(np.array(tuple(self.view.getViewDirection()), float))) < 0.0   # eksen bize bakiyorsa
        ang = (self.mouse[0] - g["p0"][0]) * 0.5 * (_gain() / DRAG_GAIN) * (-1.0 if toward else 1.0)
        v["Angle"] = float(_snap(ang, self.angle_box.value() if _ctrl_held() else 0.0))
        c, a, b, e = op["before"]
        self.cage, self.sel, self.sel_v, self.sel_e = c.copy(), set(a), set(b), set(e)
        self._apply_op(op)
        self.refresh(fast=True)

    def twist_set(self, name, value):
        """Panelden bir Twist degeri: Start / End elle verilince Auto range kapanir."""
        if self.grab and self.grab.get("kind") == "twist":
            self.end_grab(True)
        op = self.last_op
        if op is not None and op["kind"] == "twist" and name in ("Start", "End"):
            op["values"]["Auto"] = 0.0
        self.redo_last(name, value)

    def _twist_panel(self, v):
        r = lambda name: (lambda val: self.twist_set(name, val))
        fields = [("Angle", v["Angle"], "°", 2, r("Angle")), ("Start", v["Start"], " mm", 3, r("Start")),
                  ("End", v["End"], " mm", 3, r("End")), ("Cuts", v["Cuts"], "", 0, r("Cuts"))]
        fields += [("Axis point " + n, v["P" + n.lower()], " mm", 3, r("P" + n.lower())) for n in "XYZ"]
        fields += [("Axis direction " + n, v["D" + n.lower()], "", 4, r("D" + n.lower())) for n in "XYZ"]
        fo = int(v["Falloff"])
        checks = [(n, fo == i, lambda on, i=i: self.twist_set("Falloff", float(i)))
                  for i, n in enumerate(TWIST_FALLOFFS)]
        checks += [("Reverse", v["Reverse"] >= 0.5, lambda on: self.twist_set("Reverse", float(bool(on)))),
                   ("Auto range", v["Auto"] >= 0.5, lambda on: self.twist_set("Auto", float(bool(on))))]
        side = [("Axis X", lambda: self.array_axis_dir((1.0, 0.0, 0.0))),
                ("Axis Y", lambda: self.array_axis_dir((0.0, 1.0, 0.0))),
                ("Axis Z", lambda: self.array_axis_dir((0.0, 0.0, 1.0)))]
        return T("Twist"), fields, self.drag is None, checks, None, [("Draw axis", self.array_pick)], side

    # --- Circle (Blender LoopTools Circle): secili nokta / kenar halkasi / cokgen bolgesi daire olur ---
    # Cokgenlerde: bolgenin siniri (tek halka) daireye dizilir, ic noktalar acisi ve yaricapi orantili tasinir.
    # Kenarlarda: tek halka ya da zincir (zincir de tam daire olur). Noktalarda: merkez etrafindaki aciya gore
    # siralanir. Daire en uygun duzlemde, merkezi secimin agirlik merkezindedir. Panel: Radius (mm, elle verilince
    # Auto radius kapanir), Strength (%), Auto radius (ortalama yaricap), Even spacing (esit aralik). Simetri yok sayilir.

    def circle(self):
        """Modify > Circle: secimi daireye cevirir; yaricapi panelden sayiyla verebilirsin."""
        if self.mode == "object":
            self.info = T("Circle: only in edit mode")
            self.refresh()
            return
        if self.drag or self.grab:
            return
        cage, sort, inner = self.cage, False, []
        if self.mode == "face":
            faces = set(self.sel)
            cnt = {}
            for fi in faces:
                f = cage.F[fi]
                for i in range(len(f)):
                    e = _edge(f[i], f[(i + 1) % len(f)])
                    cnt[e] = cnt.get(e, 0) + 1
            res = _chain_order([e for e, c in cnt.items() if c == 1]) if faces else None
            if res is None or not res[1]:
                self.info = T("Circle: select polygons that form one region without holes, an edge loop or points")
                self.refresh()
                return
            order = res[0]
            inner = sorted({x for fi in faces for x in cage.F[fi]} - set(order))
        elif self.mode == "edge":
            res = _chain_order(self.sel_e) if self.sel_e else None
            if res is None:
                self.info = T("Circle: select one edge loop or chain, polygons or points")
                self.refresh()
                return
            order = res[0]
        else:
            pts = sorted(self.sel_v)
            res = _chain_order([e for e in cage.edges() if e <= self.sel_v]) if len(pts) >= 3 else None
            if res is not None and len(res[0]) == len(pts):
                order = res[0]
            else:
                order, sort = pts, True
        if len(order) < 3:
            self.info = T("Circle: select at least 3 points")
            self.refresh()
            return
        if self._locked_block(sorted(set(order) | set(inner))):
            return
        self.drop_tool()
        self.push()
        op = self.new_op("circle", dict(Radius=0.0, Strength=100.0, Auto=1.0, Even=1.0),
                         order=order, inner=inner, sort=sort)
        self._apply_op(op)
        self.last_op, self.op = op, None
        self.refresh()

    def _circle(self, v, a):
        """Kafesi 'before' durumundan daireye cevirir; hata iletisi ya da None dondurur."""
        cage = self.cage
        order, inner = list(a["order"]), list(a["inner"])
        W = self.l2w(cage.V[order])
        Wi = self.l2w(cage.V[inner]) if inner else np.zeros((0, 3))
        n_pts = len(order)
        c = W.mean(0)
        _, _, vt = np.linalg.svd(W - c)
        if a["sort"]:                                      # noktalar: merkez etrafindaki aciya gore sirala
            nrm = np.cross(vt[0], vt[1])
            perm = np.argsort(np.arctan2((W - c) @ vt[1], (W - c) @ vt[0]))
            order, W = [order[i] for i in perm], W[perm]
        else:
            nrm = vt[2]
            area = sum(float(nrm.dot(np.cross(W[i] - c, W[(i + 1) % n_pts] - c))) for i in range(n_pts))
            if area < 0:                                   # halka normale gore saat yonunun tersine gitsin
                nrm = -nrm
        d0 = W[0] - c
        d0 = d0 - nrm * float(d0 @ nrm)
        if np.linalg.norm(d0) < 1e-9:
            d0 = vt[0] - nrm * float(vt[0] @ nrm)
        u = d0 / np.linalg.norm(d0)
        w = np.cross(nrm, u)
        d = W - c
        x, y = d @ u, d @ w
        th = np.arctan2(y, x)
        ang = [0.0]
        for i in range(1, n_pts):
            ang.append(ang[-1] + float((th[i] - th[i - 1]) % (2 * math.pi)))
        ang = np.array(ang)
        if ang[-1] >= 2 * math.pi - 1e-6:
            return T("Circle: the points do not go around a centre")
        r = np.hypot(x, y)
        if v["Auto"] >= 0.5:
            v["Radius"] = float(r.mean())
        R = float(v["Radius"])
        if R < 1e-9:
            return T("Circle: the radius is zero")
        k = min(max(float(v["Strength"]) / 100.0, 0.0), 1.0)
        tau = 2 * math.pi
        phi = tau * np.arange(n_pts) / n_pts if v["Even"] >= 0.5 else ang
        pol = lambda ph, rho: c + rho[:, None] * (np.cos(ph)[:, None] * u + np.sin(ph)[:, None] * w)
        new_b = W + k * (pol(phi, np.full(n_pts, R)) - W)
        new_i = None
        if inner:                                          # ic noktalar: aci parcali dogrusal, yaricap sinira orantili
            di = Wi - c
            xi, yi, hi = di @ u, di @ w, di @ nrm
            thi = np.arctan2(yi, xi) % tau
            ang_e, phi_e = np.append(ang, tau), np.append(phi, phi[0] + tau)
            j = np.clip(np.searchsorted(ang_e, thi, side="right") - 1, 0, n_pts - 1)
            t = (thi - ang_e[j]) / np.maximum(ang_e[j + 1] - ang_e[j], 1e-12)
            phi_i = phi_e[j] + t * (phi_e[j + 1] - phi_e[j])
            rho = np.empty(len(inner))
            for q in range(len(inner)):                    # isin eski sinirla kesistigi yer: R0
                p0, p1 = np.array([x[j[q]], y[j[q]]]), np.array([x[(j[q] + 1) % n_pts], y[(j[q] + 1) % n_pts]])
                e, dd = np.array([math.cos(thi[q]), math.sin(thi[q])]), p1 - p0
                den = e[0] * dd[1] - e[1] * dd[0]
                lam = (p0[0] * dd[1] - p0[1] * dd[0]) / den if abs(den) > 1e-12 else 0.0
                if lam <= 1e-9:                            # kesisim yok: iki ucun yaricapini ara degerle
                    lam = float(r[j[q]] + t[q] * (r[(j[q] + 1) % n_pts] - r[j[q]]))
                rho[q] = math.hypot(xi[q], yi[q]) * R / max(lam, 1e-9)
            tgt = pol(phi_i, rho) + hi[:, None] * nrm
            new_i = Wi + k * (tgt - Wi)
        cage.V[order] = self.w2l(new_b)
        if new_i is not None:
            cage.V[inner] = self.w2l(new_i)
        self.info = T("Circle: radius %.3f mm") % R + ("   " + T("(symmetry is ignored)") if self.sym else "")
        return None

    def circle_set(self, name, value):
        """Panelden bir Circle degeri: Radius elle verilince Auto radius kapanir."""
        op = self.last_op
        if op is not None and op["kind"] == "circle" and name == "Radius":
            op["values"]["Auto"] = 0.0
        self.redo_last(name, value)

    def _circle_panel(self, v):
        r = lambda name: (lambda val: self.circle_set(name, val))
        fields = [("Radius", v["Radius"], " mm", 3, r("Radius")), ("Strength", v["Strength"], " %", 1, r("Strength"))]
        checks = [("Auto radius", v["Auto"] >= 0.5, lambda on: self.circle_set("Auto", float(bool(on)))),
                  ("Even spacing", v["Even"] >= 0.5, lambda on: self.circle_set("Even", float(bool(on))))]
        return T("Circle"), fields, self.drag is None, checks

    def _axis_op(self):
        """Eksen verilebilen son islem (Array / Twist) ya da None; Twist surerken once onaylanir."""
        if self.grab and self.grab.get("kind") == "twist":
            self.end_grab(True)
        op = self.last_op
        return op if op is not None and op["kind"] in AXIS_OPS else None

    def user_command(self, name):
        """Kullanicinin tus atadigi komut (set_user_key)."""
        if self.drag is not None or self.grab or self.box is not None or name not in USER_KEYS:
            return
        getattr(self, name)()

    # --- Distance (olcum): iki nokta arasi mesafe ---
    # Yakalama (Q'dan bagimsiz, hep acik): kafes noktasi, kenar ortasi, cokgen merkezi, kenar; baska nesnelerin kose /
    # kenarlari; yoksa yuzey noktasi (yalniz gecici olcu). Keep acikken olcu kafese baglanir (Cage.dims): noktalar
    # tasinince / dondurulunce olcu canli guncellenir, geri alma ve belgeye kayit kendiliginden olur, uzerine
    # tiklayinca silinir. Keep kapaliyken ve Object modunda olcu gecicidir: sonraki olcuye ya da aractan cikinca kalkar.
    DIM_PICK = 12.0

    def dim_keep(self):
        return _param().GetBool("DimKeep", False)

    def dim_set_keep(self, on):
        _param().SetBool("DimKeep", bool(on))
        self.refresh()

    def dim(self):
        """Distance aracini ac / kapat (Edit ve Object modunda)."""
        if self.drag or self.grab or self.box is not None:
            return
        if self.armed == "dim":
            self.drop_tool()
            self.refresh()
            return
        self.drop_tool()
        self.armed = "dim"
        self._dim_end()
        self.info = T("Distance: click two points (Esc / right click: off)")
        self.refresh()

    def _dim_end(self):
        self.dim_p0 = self.dim_tmp = self.dim_ext = self.dim_live_q = self.dim_hot = None
        self.dim_sig = None
        self._dim_cd = None
        _quiet(self._dim_tick)

    def dim_escape(self):
        """Esc / sag tik: bekleyen ilk nokta varsa vazgec, yoksa araci kapat."""
        if self.dim_p0 is not None:
            self.dim_p0 = self.dim_live_q = None
            self._line_hide()
            self.info = T("Distance: click two points (Esc / right click: off)")
        else:
            self.drop_tool()
        self._dim_tick()
        self.refresh()

    def _cam_sig(self):
        cam = self.view.getCameraNode()
        try:
            h = cam.height.getValue()
        except AttributeError:
            h = cam.heightAngle.getValue()
        return (tuple(cam.position.getValue().getValue()), tuple(self.view.getCameraOrientation().Q), h,
                tuple(self.view.getSize()))

    def _proj_np(self, P):
        """Dunya noktalari -> ekran konumu (Coin olay duzeni: sol alt kokenli aygit pikseli), tek seferde."""
        P = np.asarray(P, float).reshape(-1, 3)
        vv, W, H = _vv(self.view)
        M = np.array(vv.getMatrix().getValue()).reshape(4, 4)
        h = np.c_[P, np.ones(len(P))] @ M
        h = h[:, :3] / h[:, 3:4]
        return np.c_[(h[:, 0] + 1) / 2 * W, (h[:, 1] + 1) / 2 * H]

    def _dim_world(self, a):
        if a[0] in ("v", "e", "f") and self.mode == "object":
            return None                                    # Object modunda kafes yok: gecici olcu dunya noktalariyla
        return _dim_point(a, self.cage.V, lambda p: self.l2w(np.asarray(p)[None])[0])

    def _dim_text(self, P, Q):
        d = np.asarray(Q) - np.asarray(P)
        return T("Distance %.3f mm   dX %.3f   dY %.3f   dZ %.3f") % (float(np.linalg.norm(d)), d[0], d[1], d[2])

    # yakalama hedefleri

    def _dim_cage_data(self):
        c = self.cage
        key = (self._cam_sig(), len(c.V), len(c.F), hash(c.V.tobytes()))
        d = self._dim_cd
        if d is not None and d["key"] == key:
            return d
        V = c.V
        E = np.array([sorted(e) for e in c.edges()], int).reshape(-1, 2)
        Fc = np.array([V[f].mean(0) for f in c.F]).reshape(-1, 3)
        none = np.zeros((0, 2))
        self._dim_cd = d = dict(key=key, E=E, S=self._proj_np(self.l2w(V)),
                                M=self._proj_np(self.l2w(V[E].mean(1))) if len(E) else none,
                                C=self._proj_np(self.l2w(Fc)) if len(Fc) else none)
        return d

    def _dim_cage_snap(self, p):
        c = self.cage
        if not len(c.V):
            return None
        d = self._dim_cage_data()
        V, E, S = c.V, d["E"], d["S"]
        R = self.DIM_PICK
        ds = np.linalg.norm(S - p, axis=1)
        i = int(np.argmin(ds))
        dm = np.linalg.norm(d["M"] - p, axis=1) if len(E) else np.array([1e9])
        j = int(np.argmin(dm))
        if min(ds[i], dm[j]) < R:                          # nokta ve kenar ortasi, yakin olan
            if ds[i] <= dm[j]:
                return dict(s=tuple(S[i]), w=self.l2w(V[i][None])[0], kind="point", a=["v", i])
            return dict(s=tuple(d["M"][j]), w=self.l2w(V[E[j]].mean(0)[None])[0], kind="mid",
                        a=["e", int(E[j, 0]), int(E[j, 1]), 0.5])
        if len(d["C"]):
            dc = np.linalg.norm(d["C"] - p, axis=1)
            k = int(np.argmin(dc))
            if dc[k] < R:
                f = [int(x) for x in c.F[k]]
                return dict(s=tuple(d["C"][k]), w=self.l2w(V[f].mean(0)[None])[0], kind="face", a=["f", f])
        if len(E):
            a, b = S[E[:, 0]], S[E[:, 1]]
            ab = b - a
            t = np.clip(((p - a) * ab).sum(1) / np.maximum((ab * ab).sum(1), 1e-12), 0, 1)
            dd = np.linalg.norm(a + ab * t[:, None] - p, axis=1)
            k = int(np.argmin(dd))
            if dd[k] < 8:
                i0, i1 = int(E[k, 0]), int(E[k, 1])
                w = self.l2w((V[i0] + t[k] * (V[i1] - V[i0]))[None])[0]
                return dict(s=tuple(a[k] + t[k] * ab[k]), w=w, kind="edge", a=["e", i0, i1, float(t[k])])
        return None

    def _dim_ext_data(self):
        """Baska nesnelerin yakalanacak noktalari ve kenarlari (dunya); arac acilinca bir kez toplanir."""
        if self.dim_ext is not None:
            return self.dim_ext
        import Part
        P, A, B = [], [], []
        for o in App.ActiveDocument.Objects:
            try:
                if not o.ViewObject.Visibility and not (is_subd(o) and self.mode == "object"):
                    continue
            except Exception:
                continue
            if is_subd(o):
                if o is self.obj and self.mode != "object":    # duzenlenen kafes zaten yakalaniyor
                    continue
                try:
                    cg = cage_from_obj(o)
                    if not len(cg.V):
                        continue
                    pl = o.Placement
                    R = np.array(pl.Rotation.toMatrix().A).reshape(4, 4)[:3, :3]
                    Pw = cg.V @ R.T + np.array(tuple(pl.Base))
                    E = np.array([sorted(e) for e in cg.edges()], int).reshape(-1, 2)
                except Exception:
                    continue
                P += Pw.tolist()
                A += Pw[E[:, 0]].tolist()
                B += Pw[E[:, 1]].tolist()
                continue
            if "Samples" in o.PropertiesList or not hasattr(o, "Shape") or o.TypeId.startswith("App::"):
                continue                                   # Origin eksenleri, gruplar, SubD'den turetilenler
            try:
                sh = Part.getShape(o)
            except Exception:
                continue
            if sh.isNull():
                continue
            P += [tuple(x.Point) for x in sh.Vertexes]
            for e in sh.Edges:
                q = [tuple(x) for x in e.discretize(Number=2 if _edge_turn(e) < 1e-3 else 33)]
                A += q[:-1]
                B += q[1:]
        self.dim_ext = dict(P=np.array(P, float).reshape(-1, 3), A=np.array(A, float).reshape(-1, 3),
                            B=np.array(B, float).reshape(-1, 3), key=None)
        return self.dim_ext

    def _dim_ext_snap(self, p):
        d = self._dim_ext_data()
        if not len(d["P"]) and not len(d["A"]):
            return None
        sig = self._cam_sig()
        if d["key"] != sig:
            none = np.zeros((0, 2))
            d.update(Ps=self._proj_np(d["P"]) if len(d["P"]) else none, As=self._proj_np(d["A"]) if len(d["A"]) else none,
                     Bs=self._proj_np(d["B"]) if len(d["B"]) else none, key=sig)
        if len(d["Ps"]):
            dist = np.linalg.norm(d["Ps"] - p, axis=1)
            i = int(np.argmin(dist))
            if dist[i] < self.DIM_PICK:
                return dict(s=tuple(d["Ps"][i]), w=d["P"][i], kind="point", a=["w"] + d["P"][i].tolist())
        if len(d["As"]):
            ab = d["Bs"] - d["As"]
            t = np.clip(((p - d["As"]) * ab).sum(1) / np.maximum((ab * ab).sum(1), 1e-9), 0, 1)
            dist = np.linalg.norm(d["As"] + t[:, None] * ab - p, axis=1)
            i = int(np.argmin(dist))
            if dist[i] < 10:
                w = d["A"][i] + t[i] * (d["B"][i] - d["A"][i])
                return dict(s=tuple(d["As"][i] + t[i] * ab[i]), w=w, kind="edge", a=["w"] + w.tolist())
        return None

    def _dim_ray_smooth(self, pos):
        """Imlecin altindaki SubPatch yuzey noktasi (dunya) ya da None."""
        S, Q = self._smooth
        o, d = self.ray((int(pos[0]), int(pos[1])))
        Q = np.asarray(Q)
        tri = np.concatenate([Q[:, [0, 1, 2]], Q[:, [0, 2, 3]]])
        Tr = S[tri]
        e1, e2 = Tr[:, 1] - Tr[:, 0], Tr[:, 2] - Tr[:, 0]
        pv = np.cross(d, e2)
        det = (e1 * pv).sum(1)
        ok = np.abs(det) > 1e-12
        inv = 1.0 / np.where(ok, det, 1.0)
        sv = o - Tr[:, 0]
        u = (sv * pv).sum(1) * inv
        qv = np.cross(sv, e1)
        v = (qv * d).sum(1) * inv
        t = (e2 * qv).sum(1) * inv
        hit = ok & (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 0)
        if not hit.any():
            return None
        k = int(np.argmin(np.where(hit, t, np.inf)))
        return self.l2w((o + d * t[k])[None])[0]

    def _dim_surface(self, pos):
        w = None
        if self.mode != "object" and len(self.cage.F):
            if self.smooth and self._smooth is not None:
                w = self._dim_ray_smooth(pos)
            else:
                h = self.pick((int(pos[0]), int(pos[1])))
                if h is not None:
                    w = self.l2w(h[1][None])[0]
        if w is None:
            info = self.view.getObjectInfo((int(pos[0]), int(pos[1])))
            if info and "x" in info:
                w = np.array([info["x"], info["y"], info["z"]], float)
        if w is None:
            return None
        return dict(s=tuple(pos), w=w, kind="surf", a=["s"] + [float(x) for x in w])

    def _dim_snap(self, pos):
        p = np.array(pos, float)
        hit = self._dim_cage_snap(p) if self.mode != "object" else None
        if hit is None:
            hit = self._dim_ext_snap(p)
        if hit is None:
            hit = self._dim_surface(pos)
        return hit

    # olculer

    def _dim_resolve(self, d):
        P, Q = self._dim_world(d["a"]), self._dim_world(d["b"])
        return None if P is None or Q is None else (P, Q)

    def _dim_obj_all(self):
        """Object modunda: SubD nesnelerinin kayitli olculeri: (nesne, sira, kayit, P, Q); P / Q None olabilir."""
        out = []
        for o in self.subd_objects():
            if "Dimensions" not in o.PropertiesList or not o.Dimensions:
                continue
            try:
                ds = json.loads(o.Dimensions)
                V = np.array([tuple(q) for q in o.CagePoints], float).reshape(-1, 3)
            except ValueError:
                continue
            pl = o.Placement
            to_w = lambda q, pl=pl: np.array(tuple(pl.multVec(App.Vector(*q))))
            for i, d in enumerate(ds):
                out.append((o, i, d, _dim_point(d["a"], V, to_w), _dim_point(d["b"], V, to_w)))
        return out

    def _dim_items(self):
        """Ekranda gosterilecek olculer: dict(P, Q, kind, idx, cur); Q yoksa ikinci nokta bekleniyor."""
        out = []
        if self.mode != "object":
            for i, d in enumerate(self.cage.dims):
                pq = self._dim_resolve(d)
                if pq is not None and not d.get("hide"):
                    out.append(dict(P=pq[0], Q=pq[1], kind="keep", idx=i, d=d))
        else:
            for o, i, d, P, Q in self._dim_obj_all():
                if P is not None and Q is not None and not d.get("hide"):
                    out.append(dict(P=P, Q=Q, kind="keep", idx=None, d=d))
        out += self._handle_dim_items()
        if self.dim_tmp is not None:
            pq = self._dim_resolve(self.dim_tmp)
            if pq is not None:
                out.append(dict(P=pq[0], Q=pq[1], kind="tmp", idx=None))
        if self.dim_p0 is not None:
            P = self._dim_world(self.dim_p0["a"])
            if P is not None:
                out.append(dict(P=P, Q=self.dim_live_q, kind="live", idx=None))
        return out

    def _handle_dim_items(self):
        """Birbirine bagli grup tutamaklarinin arasindaki olculer (oklu; degeri tiklayinca degistirilir)."""
        doc = App.ActiveDocument
        if doc is None:
            return []
        key = (doc.Name, len(doc.Objects))
        if self._hd_key != key:                            # her adimda tum nesneleri gezme: yalniz sayi degisince
            self._hd_key = key
            self._hd_names = [o.Name for o in doc.Objects if _is_handle(o)]
        edit = self.mode != "object" and _alive(self.obj)
        out = []
        for name in self._hd_names:
            h = doc.getObject(name)
            if h is None or (edit and h.SubDGroup != self.obj.Name):
                continue
            g = handle_dim_geom(h)
            if g is not None:
                out.append(dict(P=np.array(tuple(g[0])), Q=np.array(tuple(g[1])), kind="handle", idx=None,
                                hd=name, val=float(g[2])))
        return out

    def _dim_hit_label(self, pos):
        """Grup olcusunun yazisi imlecin altindaysa tutamagin adi, yoksa None."""
        for h in reversed(self.dim_hits):
            if h.get("hd") is None:
                continue
            c = self._proj_np(h["mid"])[0]
            if abs(pos[0] - c[0]) <= h["w"] / 2 and abs(pos[1] - c[1]) <= h["h"] / 2:
                return h["hd"]
        return None

    def handle_dim_edit(self, name):
        """Grup olcusunun degerini yazdirir: tutamagin Distance ozelligi (ifadeyle ofsete bagli) degisir, tutamak ve
        grubun noktalari kayar."""
        doc = App.ActiveDocument
        h = doc.getObject(name) if doc is not None else None
        if h is None:
            return
        handle_dim_setup(h)
        if "Distance" not in h.PropertiesList:
            return
        ref = handle_ref(h)
        v, ok = QtWidgets.QInputDialog.getDouble(
            Gui.getMainWindow(), T("Group distance"),
            T("Distance %s - %s (mm)") % (h.Label, ref.Label if ref is not None else ""),
            float(h.Distance.Value), -1e7, 1e7, 3)
        if ok and abs(v - float(h.Distance.Value)) > 1e-9:
            doc.openTransaction(T("Group distance"))
            h.Distance = v
            doc.recompute()
            doc.commitTransaction()
        self.refresh()

    def _dim_scene(self):
        """Olculerin 3B sahnedeki dugumu: kamerayla birlikte cizilir, gorunumu dondururken titremez / geride kalmaz."""
        if self.dim_root is None:
            r = self.dim_root = coin.SoSeparator()
            ps = coin.SoPickStyle()
            ps.style = coin.SoPickStyle.UNPICKABLE         # tiklamalar olculerden gecer
            lm = coin.SoLightModel()
            lm.model = coin.SoLightModel.BASE_COLOR
            db = coin.SoDepthBuffer()                      # olcu modelin icinden de gorunsun
            db.test = False
            for n in (ps, lm, db):
                r.addChild(n)
            self.dim_group = coin.SoSeparator()
            r.addChild(self.dim_group)
        if self.parent.findChild(self.dim_root) < 0:
            self.parent.addChild(self.dim_root)
        return self.dim_group

    DIM_COL = {"keep": (163, 230, 53), "tmp": (245, 158, 11), "live": (34, 211, 238), "handle": (255, 140, 40),
               "sel": (96, 165, 250), "bad": (248, 113, 113)}   # sel: listede secili, bad: hedeften tolerans disi
    DIM_TEXT_SCALE = 1.5          # olcu yazisi ve x dugmesi: gizmo boyut carpanina ek carpan

    def _dim_label(self, text, kind, hot, with_x, k):
        """Olcu yazisi kutusu: [x | 123.456 mm] (x yalniz kalici olcuyu silmek icin). Sahne icin piksel goruntusu:
        (RGBA baytlari, genislik, yukseklik, x bolgesinin genisligi); aygit pikseli."""
        k = k * self.DIM_TEXT_SCALE
        key = (text, kind, hot, with_x, round(k, 2))
        img = self._dim_imgs.get(key)
        if img is not None:
            return img
        dpr = self.gv.viewport().devicePixelRatioF()
        f = QtGui.QFont()
        f.setPixelSize(max(int(round(13 * k * dpr)), 9))
        f.setBold(True)
        fm = QtGui.QFontMetrics(f)
        pad = int(round(4 * k * dpr))
        xw = fm.height() + pad if with_x else 0
        w, h = 2 * pad + xw + fm.horizontalAdvance(text), fm.height() + pad
        qi = QtGui.QImage(w, h, QtGui.QImage.Format_ARGB32)
        qi.fill(0)
        col = QtGui.QColor(*self.DIM_COL[kind])
        p = QtGui.QPainter(qi)
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        p.setFont(f)
        p.setBrush(QtGui.QColor(24, 24, 28, 225))
        p.setPen(QtGui.QPen(col, max(1.0, k * dpr)))
        rad = 4 * k * dpr
        p.drawRoundedRect(QtCore.QRectF(0.5, 0.5, w - 1, h - 1), rad, rad)
        if with_x:
            if hot:                                        # imlec x'in uzerinde: tiklayinca silinir
                p.setPen(QtCore.Qt.NoPen)
                p.setBrush(QtGui.QColor(220, 60, 60, 235))
                p.drawRoundedRect(QtCore.QRectF(1, 1, xw, h - 2), rad, rad)
            p.setPen(QtGui.QColor(255, 255, 255) if hot else QtGui.QColor(200, 200, 205))
            p.drawText(QtCore.QRect(0, 0, xw, h), QtCore.Qt.AlignCenter, "\u00d7")
            p.setPen(QtGui.QPen(col, 1))
            p.drawLine(xw, 3, xw, h - 3)
        p.setPen(col)
        p.drawText(QtCore.QRect(xw, 0, w - xw, h), QtCore.Qt.AlignCenter, text)
        p.end()
        rgba = qi.convertToFormat(QtGui.QImage.Format_RGBA8888).mirrored(False, True)   # SoImage: alt satir once
        if len(self._dim_imgs) > 80:
            self._dim_imgs.clear()
        self._dim_imgs[key] = img = (bytes(rgba.constBits())[:w * h * 4], w, h, xw)
        return img

    def _dim_hit_x(self, pos):
        """Kalici olcunun x dugmesi imlecin (Coin konumu) altindaysa olcunun sirasi, yoksa None."""
        if self.mode == "object" or not self.dim_hits:
            return None
        for h in reversed(self.dim_hits):
            c = self._proj_np(h["mid"])[0]
            x0, y0 = c[0] - h["w"] / 2, c[1] - h["h"] / 2
            if x0 <= pos[0] <= x0 + h["xw"] and y0 <= pos[1] <= y0 + h["h"]:
                return h["idx"]
        return None

    def _dim_tick(self):
        """Sahnedeki olculeri noktalar ya da gorunen metin degisince yeniden kur (kamera donerken is yok)."""
        if self.suspended:
            return
        if self.dim_list is not None and self.dim_list.isVisible():
            self.dim_list.refresh(self._dim_rows())
        items = self._dim_items()
        if not items and self.dim_root is None:
            return
        k = min(max(self.gizmo.px / GIZMO_PX, 0.7), 3.0)   # yazi boyu gizmo boyutuyla orantili
        names = self.dim_names()
        meta = []
        for it in items:
            d = it.get("d") or {}
            Q = it["Q"]
            val = it["val"] if it.get("val") is not None else None if Q is None else float(np.linalg.norm(Q - it["P"]))
            tgt = d.get("target")
            bad = val is not None and tgt is not None and abs(val - tgt) > float(d.get("tol", DIM_TOL))
            sel = it["kind"] == "keep" and it["idx"] is not None and it["idx"] == self.dim_sel and self.mode != "object"
            meta.append(("bad" if bad else "sel" if sel else it["kind"], d.get("name", "") if names else "", val))
        vd = np.array(tuple(self.view.getViewDirection()))
        sig = (k, self.dim_hot, tuple(np.round(vd, 1)), tuple((it["kind"], it["idx"], tuple(np.round(it["P"], 3)),
                                       None if it["Q"] is None else tuple(np.round(it["Q"], 3)), m)
                                      for it, m in zip(items, meta)))
        if sig == self.dim_sig:
            return
        self.dim_sig = sig
        g = self._dim_scene()
        g.removeAllChildren()
        self.dim_hits = []
        nh = 0
        for it, (ck, name, val) in zip(items, meta):
            kind, P, Q = it["kind"], it["P"], it["Q"]
            O = None
            if kind == "handle" and Q is not None:         # grup olcusu: olcu cizgisi merkezden yana kayik, uzanti cizgili
                u = Q - P
                ln = float(np.linalg.norm(u))
                if ln > 1e-9:
                    pp = np.cross(u / ln, vd)
                    if np.linalg.norm(pp) < 1e-9:
                        pp = np.cross(u / ln, np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0]))
                    pp = pp / np.linalg.norm(pp)
                    side = 1.0 if nh % 2 == 0 else -1.0    # birbirini izleyen olculer karsi yanlara, giderek uzaga
                    O = pp * side * max(ln, 40.0) * 0.22 * (nh // 2 + 1)
                    nh += 1
            grp = coin.SoSeparator()
            bc = coin.SoBaseColor()
            bc.rgb = tuple(c / 255.0 for c in self.DIM_COL[ck])
            ds = coin.SoDrawStyle()
            ds.lineWidth = 1.6 * k * (2.0 if ck == "sel" else 1.0)
            ds.pointSize = 7 * k
            if kind == "live":
                ds.linePattern = 0xF0F0
            pts = [tuple(float(x) for x in P)] + ([tuple(float(x) for x in Q)] if Q is not None else [])
            if O is not None:                              # P -> P+O -> Q+O -> Q
                pts = [tuple(float(x) for x in e) for e in (P, P + O, Q + O, Q)]
            co = coin.SoCoordinate3()
            co.point.setValues(0, len(pts), pts)
            for n in (bc, ds, co):
                grp.addChild(n)
            if Q is not None:
                ls = coin.SoLineSet()
                ls.numVertices = len(pts) if O is not None else 2
                grp.addChild(ls)
            ps = coin.SoPointSet()
            ps.numPoints = len(pts)
            grp.addChild(ps)
            if O is not None:
                P, Q = P + O, Q + O                        # ok ucu ve yazi olcu cizgisinde
            if Q is not None and kind == "handle":         # iki ucta ok ucu: ok yonu ile bakis yonune dik
                u = Q - P
                ln = float(np.linalg.norm(u))
                if ln > 1e-9:
                    u = u / ln
                    pp = np.cross(u, np.array(tuple(self.view.getViewDirection())))
                    if np.linalg.norm(pp) < 1e-9:
                        pp = np.cross(u, np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0]))
                    pp = pp / np.linalg.norm(pp) * ln * 0.032
                    L = u * ln * 0.08
                    ends = [P + L + pp, P, P + L - pp, Q - L + pp, Q, Q - L - pp]
                    co2 = coin.SoCoordinate3()
                    co2.point.setValues(0, 6, [tuple(float(x) for x in e) for e in ends])
                    al = coin.SoLineSet()
                    al.numVertices.setValues(0, 2, [3, 3])
                    grp.addChild(co2)
                    grp.addChild(al)
            if Q is not None:
                with_x = kind == "keep" and it["idx"] is not None and self.mode != "object"
                hot = with_x and it["idx"] == self.dim_hot
                text = ("%s  " % name if name else "") + "%.3f mm" % val
                data, w, h, xw = self._dim_label(text, ck, hot, with_x, k)
                mid = (P + Q) / 2
                tr = coin.SoTranslation()
                tr.translation = tuple(float(x) for x in mid)
                im = coin.SoImage()
                im.image.setValue(coin.SbVec2s(w, h), 4, data)
                im.horAlignment = coin.SoImage.CENTER
                im.vertAlignment = coin.SoImage.HALF
                grp.addChild(tr)
                grp.addChild(im)
                if with_x:
                    self.dim_hits.append(dict(idx=it["idx"], mid=mid, w=w, h=h, xw=xw))
                elif kind == "handle":
                    self.dim_hits.append(dict(idx=None, hd=it["hd"], mid=mid, w=w, h=h, xw=0))
            g.addChild(grp)

    def _dim_hot_check(self):
        """Distance kapaliyken: imlec kalici olcunun x dugmesinin uzerine gelince x kirmiziya doner."""
        self._pending = False
        h = self._dim_hit_x(self.mouse)
        if h != self.dim_hot:
            self.dim_hot = h
            self._dim_tick()

    def _dim_hover(self):
        """Distance acik: imlecin altindaki hedefi goster, ilk nokta varsa canli mesafeyi yaz."""
        self._pending = False
        pos = self.mouse
        self.dim_hot = self._dim_hit_x(pos) if self.dim_p0 is None else None
        hit = None if self.dim_hot is not None else self._dim_snap(pos)
        self._line_mark((hit["s"], hit["w"], hit["kind"]) if hit else None)
        if self.dim_p0 is not None:
            P = self._dim_world(self.dim_p0["a"])
            self.dim_live_q = hit["w"] if hit else None
            if hit and P is not None:
                self.dim_last = (P, hit["w"])
                self.info = self._dim_text(P, hit["w"])
        self._dim_tick()
        if self.dim_p0 is not None and self.overlay is not None:
            self.overlay.sync()                            # panelin sayilari canli

    def dim_click(self, pos):
        self._pending = False
        hit = self._dim_snap(pos)
        if hit is None:
            self.info = T("Distance: nothing to snap to here")
            self.refresh()
            return
        if self.dim_p0 is None:                            # ilk nokta (onceki gecici olcu kalkar)
            self.dim_p0, self.dim_tmp = hit, None
            self.dim_live_q = None
            self.info = T("Distance: click the second point (Esc: cancel)")
            self._dim_tick()
            self.refresh()
            return
        d = dict(a=self.dim_p0["a"], b=hit["a"])
        self.dim_p0 = self.dim_live_q = None
        P, Q = self._dim_world(d["a"]), self._dim_world(d["b"])
        if (self.dim_keep() and self.mode != "object" and d["a"][0] in DIM_KEEP and d["b"][0] in DIM_KEEP):
            self.push()
            self.cage.dims.append(d)
            self.last_op = None
            self.dirty = True
            self.dim_tmp = None
        else:                                              # gecici: yuzey noktasi var, Keep kapali ya da Object modu
            self.dim_tmp = d
        if P is not None and Q is not None:
            self.dim_last = (P, Q)
            self.info = self._dim_text(P, Q)
        self._line_hide()
        self._dim_tick()
        self.refresh()

    def _dim_remove(self, i):
        if 0 <= i < len(self.cage.dims):
            self.push()
            del self.cage.dims[i]
            self.last_op = None
            self.dirty = True
        self.dim_hot = None
        self._dim_tick()
        self.refresh()

    def dim_clear(self):
        """Panel: gecici ve kalici butun olculeri sil (kalicilar geri alinabilir)."""
        self.dim_tmp = self.dim_p0 = self.dim_live_q = None
        self._line_hide()
        if self.mode != "object" and self.cage.dims:
            self.push()
            self.cage.dims = []
            self.last_op = None
            self.dirty = True
        self.info = ""
        self._dim_tick()
        self.refresh()

    # --- Distances listesi: ad, hedef, goster / gizle, sil, CSV ---
    def dim_names(self):
        return _param().GetBool("DimNames", True)

    def dim_set_names(self, on):
        _param().SetBool("DimNames", bool(on))
        self.dim_sig = None
        self._dim_tick()

    def _dim_rows(self):
        """Liste penceresi icin satirlar. Edit modu: duzenlenen nesnenin olculeri (duzenlenebilir); Object modu:
        butun SubD nesnelerinin olculeri (yalniz gosterim)."""
        def row(obj, idx, d, P, Q, editable):
            ok = P is not None and Q is not None
            v = (Q - P) if ok else None
            return dict(obj=obj, idx=idx, name=d.get("name", ""), hide=bool(d.get("hide")), target=d.get("target"),
                        tol=float(d.get("tol", DIM_TOL)), dist=float(np.linalg.norm(v)) if ok else None, d3=v,
                        editable=editable)
        rows = []
        if self.mode != "object":
            label = self.obj.Label if self.obj is not None else ""
            for i, d in enumerate(self.cage.dims):
                pq = self._dim_resolve(d) or (None, None)
                rows.append(row(label, i, d, pq[0], pq[1], True))
        else:
            for o, i, d, P, Q in self._dim_obj_all():
                rows.append(row(o.Label, i, d, P, Q, False))
        doc = App.ActiveDocument                           # grup tutamaklari arasi olculer (Distance duzenlenebilir)
        for it in self._handle_dim_items():
            h = doc.getObject(it["hd"])
            ref = handle_ref(h)
            v = it["Q"] - it["P"]
            rows.append(dict(obj=h.SubDGroup, idx="h:" + it["hd"], name="%s > %s" % (h.Label, ref.Label if ref is not None else ""),
                             hide=False, target=None, tol=DIM_TOL, dist=float(it["val"]), d3=v,
                             editable=self.mode != "object", handle=True, hd=it["hd"]))
        return rows

    def handle_dim_set(self, name, v):
        """Listeden bir grup olcusunun mesafesi yazildi: tutamagin Distance ozelligi degisir (geri alinabilir)."""
        doc = App.ActiveDocument
        h = doc.getObject(name) if doc is not None else None
        if h is None:
            return
        handle_dim_setup(h)
        if "Distance" in h.PropertiesList and abs(v - float(h.Distance.Value)) > 1e-9:
            doc.openTransaction(T("Group distance"))
            h.Distance = v
            doc.recompute()
            doc.commitTransaction()
        self.refresh()

    def dim_edit(self, idx, field, value):
        """Listeden bir olcunun ad / gizli / hedef / tolerans degeri degisti (geri alinabilir)."""
        if self.mode == "object" or not 0 <= idx < len(self.cage.dims):
            return
        d = self.cage.dims[idx]
        empty = (value is None) or (field == "name" and value == "") or (field == "hide" and not value)
        if (empty and field not in d) or (not empty and d.get(field) == value):
            return
        self.push()
        if empty:
            d.pop(field, None)
        else:
            d[field] = value
        self.last_op = None
        self.dirty = True
        self.dim_sig = None
        self._dim_tick()

    def dim_select(self, idx):
        """Listede secilen olcu 3B gorunumde mavi ve kalin cizilir."""
        if idx != self.dim_sel:
            self.dim_sel = idx
            self.dim_sig = None
            self._dim_tick()

    def dim_show_all(self, on):
        if self.mode == "object" or not self.cage.dims:
            return
        if all(bool(d.get("hide")) != bool(on) for d in self.cage.dims):
            return
        self.push()
        for d in self.cage.dims:
            if on:
                d.pop("hide", None)
            else:
                d["hide"] = True
        self.last_op = None
        self.dirty = True
        self.dim_sig = None
        self._dim_tick()

    def dim_list_toggle(self):
        if self.dim_list is None:
            self.dim_list = _DimList(self)
        on = not self.dim_list.isVisible()
        self.dim_list.setVisible(on)
        if on:
            self.dim_list.raise_()
            self.dim_list.refresh(self._dim_rows(), force=True)
        if self.overlay is not None:
            self.overlay.sync()

    def dim_export(self):
        """Olculeri CSV olarak kaydet (Edit modu: acik nesne, Object modu: hepsi)."""
        import csv
        rows = self._dim_rows()
        if not rows:
            self.info = T("Distances: nothing to export")
            self.refresh()
            return
        doc = App.ActiveDocument
        base = os.path.dirname(doc.FileName) if doc is not None and doc.FileName else os.path.expanduser("~")
        path, _ = QtWidgets.QFileDialog.getSaveFileName(Gui.getMainWindow(), T("Export distances"),
                                                       os.path.join(base, "distances.csv"), "CSV (*.csv)")
        if not path:
            return
        try:
            with open(path, "w", newline="", encoding="utf-8-sig") as f:      # BOM: Excel Turkce karakterleri tanir
                w = csv.writer(f)
                w.writerow(["Object", "Name", "Distance (mm)", "dX (mm)", "dY (mm)", "dZ (mm)", "Target (mm)",
                            "Tolerance (mm)", "Deviation (mm)", "Status"])
                for r in rows:
                    ok = r["dist"] is not None
                    tgt = r["target"]
                    dev = r["dist"] - tgt if ok and tgt is not None else None
                    status = "" if dev is None else "OK" if abs(dev) <= r["tol"] else "OUT"
                    f4 = lambda x: "" if x is None else "%.4f" % x
                    w.writerow([r["obj"], r["name"] or T("Distance %d") % (r["idx"] + 1), f4(r["dist"]),
                                *(f4(r["d3"][i]) if ok else "" for i in range(3)), f4(tgt),
                                f4(r["tol"]) if tgt is not None else "", f4(dev), status])
            self.info = T("Distances exported: %s") % path
        except OSError as e:
            self.info = T("Distances: could not write the file (%s)") % e
        self.refresh()

    def _dim_panel(self):
        P, Q = self.dim_last if self.dim_last is not None else (None, None)
        d = np.zeros(3) if P is None or Q is None else np.asarray(Q) - np.asarray(P)
        noop = lambda v: None
        fields = [("Distance", float(np.linalg.norm(d)), " mm", 3, noop)] + \
                 [(n, float(d[i]), " mm", 3, noop) for i, n in enumerate(("dX", "dY", "dZ"))]
        checks = [("Keep (fixed)", self.dim_keep(), self.dim_set_keep)] if self.mode != "object" else []
        return T("Distance"), fields, False, checks, None, [("Clear all", self.dim_clear)]

    # --- Bridge (Blender Bridge Edge Loops) ---
    # Kenar modu: iki acik kenar halkasi / zinciri; yuz modu: iki yuz bolgesi (yuzler silinir, sinirlari baglanir).
    # Iki tarafin kenar sayisi esit olmali. Noktalar en yakin eslesmeye gore hizalanir (kapali halkada donus ve yon),
    # Twist eslesmeyi kaydirir. Panel: Number of Cuts, Twist, Smoothness (yuzeye teget gecis), Profile Factor.
    BRIDGE = ("Number of Cuts", "Twist", "Smoothness", "Profile Factor")

    def bridge(self):
        if self.drag or self.grab or self.mode not in ("edge", "face"):
            return
        self.drop_tool()
        self.push()
        op = self.new_op("bridge", dict(self.bridge_vals))
        err = self._bridge(op["values"])
        if err:
            self.op = None
            self.do_undo(False)
            self.info = err
        else:
            self.last_op, self.op = op, None
        self.refresh()

    def bridge_option(self, on):
        """Bridge panelindeki "Auto direction": ada (bagimsiz) yuz uclarinin yonu oteki uca cevrilsin mi."""
        op = self.last_op
        if op is None or op["kind"] != "bridge":
            return
        op["values"]["Auto direction"] = 1.0 if on else 0.0
        self._reapply(op)

    def flip_faces(self):
        """Modify > Flip (Shift+F, LightWave Flip): secili cokgenlerin (secim yoksa hepsinin) yonunu, yani
        normalini cevirir. Kenar keskinlikleri ayni kalir."""
        if self.mode == "object" or self.drag or self.grab or not self.cage.F:
            return
        faces = self.selected_faces() or set(range(len(self.cage.F)))
        if self.sym:
            faces = self.mirror_faces(faces, self.mirror_map())
        self.drop_tool()
        self.push()
        self.last_op = None
        for fi in faces:
            self.cage.F[fi] = self.cage.F[fi][::-1]
        self.info = T("Flip: %d polygons") % len(faces)
        self.refresh()

    def _bridge(self, v):
        """Bridge'i uygular; olmazsa hata metni dondurur (kafes degismez)."""
        cage = self.cage
        v["Number of Cuts"] = float(min(max(int(round(v["Number of Cuts"])), 0), 64))
        v["Twist"] = float(int(round(v["Twist"])))
        self.bridge_vals = dict(v)
        prep = self._bridge_prep(int(v["Twist"]), bool(v.get("Auto direction", 1.0)))
        if isinstance(prep, str):
            return prep
        A, B, ca, dA, dB, flip, region = prep
        new = cage.bridge(A, B, ca, int(v["Number of Cuts"]), max(v["Smoothness"], 0.0), v["Profile Factor"],
                          dA, dB, flip)
        if region:
            cage.delete(region)                            # yeni yuzler sonda: indeksleri silinenler kadar kayar
            new = [fi - len(region) for fi in new]
        self.sel, self.sel_v, self.sel_e = set(new) if self.mode == "face" else set(), set(), set()
        self.info = T("Bridge: %d polygons") % len(new)
        return None

    def _bridge_ends(self):
        """Secimin acik sinir zincirleri: yuz modunda secili bolgenin siniri, kenar modunda secili acik kenarlar.
        (zincirler, kenar -> yuzler, bolge) ya da hata metni."""
        cage = self.cage
        ef = cage.edge_faces()
        region = set(self.sel) if self.mode == "face" else set()
        if region:                                         # bolge siniri: bolgede tek yuzu olan kenarlar
            edges = {e for e, fs in ef.items() if sum(f in region for f in fs) == 1}
        else:
            edges = set(self.sel_e)
            if any(len(ef.get(e, ())) != 1 for e in edges):
                return T("Bridge: select open (boundary) edges")
        return cage.edge_chains(edges), ef, region, edges

    def _bridge_prep(self, twist=0, auto=True):
        """Bridge / Rail bridge icin iki ucun eslenmesi: (A, B, kapali, dA, dB, flip, bolge) ya da hata metni.
        dA / dB: uclarda yuzeyden disari birim yonler; flip: yeni seridin yuz sirasi."""
        cage = self.cage
        res = self._bridge_ends()
        if isinstance(res, str):
            return res
        chains, ef, region, edges = res
        err = T("Bridge: select two loops (edges or face regions) with the same number of edges")
        if not chains or len(chains) != 2:
            return err
        (A, ca), (B, cb) = chains
        if ca != cb or len(A) != len(B) or len(A) < 2:
            return err
        PA, PB = cage.V[A], cage.V[B]
        if ca:                                             # kapali: en iyi donus + yon (merkezlere gore)
            QA, QB = PA - PA.mean(0), PB - PB.mean(0)
            n = len(A)
            best = None
            for rev in (False, True):
                order = B[::-1] if rev else B
                Q = QB[::-1] if rev else QB
                for k in range(n):
                    d = float(((QA - np.roll(Q, -k, 0)) ** 2).sum())
                    if best is None or d < best[0]:
                        best = (d, order[k:] + order[:k])
            B = best[1]
            t = int(twist) % n
            B = B[t:] + B[:t]
        elif (np.linalg.norm(PA[0] - PB[0]) + np.linalg.norm(PA[-1] - PB[-1]) >
              np.linalg.norm(PA[0] - PB[-1]) + np.linalg.norm(PA[-1] - PB[0])):
            B = B[::-1]
        # yuzeyden disari yonler (Smoothness icin): yuz modunda silinen bolgenin normali, kenar modunda yuzeyin devami
        if region:
            acc = {}
            for fi in region:
                for x in cage.F[fi]:
                    acc[x] = acc.get(x, np.zeros(3)) + cage.normal(fi)
            out = {x: d / max(np.linalg.norm(d), 1e-12) for x, d in acc.items()}
        else:
            out = cage._outward(cage._open_directed(edges))
        dA = np.array([out.get(x, np.zeros(3)) for x in A])
        dB = np.array([out.get(x, np.zeros(3)) for x in B])
        # ada (hicbir yere bagli olmayan, or. yapistirilmis) yuz bolgesinin normali keyfi: Auto direction ile
        # oteki uca donuk alinir (yoksa egri arkasindan dolanir). Yuzeye gomulu uclarin yonune dokunulmaz.
        island = lambda C: bool(region) and all(len(ef.get(_edge(C[i], C[(i + 1) % len(C)]), ())) == 1
                                                for i in range(len(C) if ca else len(C) - 1))
        isA, isB = island(A), island(B)
        if auto:
            gap = PB.mean(0) - PA.mean(0)
            if isA and dA.sum(0).dot(gap) < 0:
                dA = -dA
            if isB and dB.sum(0).dot(gap) > 0:
                dB = -dB
        # yuz yonu gomulu uctaki komsu yuzden (ada yuzun yonu keyfi): kalan komsu yuz kenari ayni yonde
        # geciyorsa serit ters sirali. Serit A kenarini A[0] -> A[1], B kenarini B[1] -> B[0] yonunde gecer.
        C = B if isA and not isB else A
        e0 = _edge(C[0], C[1])
        fs = ef.get(e0, ())
        keep = [f for f in fs if f not in region]
        f = cage.F[(keep or fs)[0]]
        same = any(f[i] == C[0] and f[(i + 1) % len(f)] == C[1] for i in range(len(f)))
        flip = same if keep else not same                  # silinen yuzun yerini alan seridin yonu onunla ayni
        if C is B:
            flip = not flip
        return A, B, ca, dA, dB, flip, region

    # --- Rail bridge (kilavuzlu kopru): kilavuz egri, egriye dik halkalar; ara noktalar tasinir / eklenir ---
    # Kilavuz kafeste saklanir (Cage.rails): geri alma ve belgeye kayit kendiliginden. Duzenlerken tiklama ara
    # noktayi (tutamak) ya da egri parcasini secer; gizmo / panel secili noktayi tasir, I parcaya halka ekler, Del
    # noktayi siler, Enter biter. Sonradan kilavuza cift tiklayinca yeniden duzenlenir.

    def rail_bridge(self):
        """Polygon > Rail bridge: iki secili bolge / acik kenar halkasi arasinda ya da tek bolgeden uzayan,
        kilavuz egrili boru; egri secili poligona dik cikar."""
        if self.drag or self.grab or self.mode not in ("edge", "face"):
            return
        res = self._bridge_ends()
        if isinstance(res, str):
            self.info = res
            self.refresh()
            return
        chains, _, region, _ = res
        if len(chains) == 2:
            prep = self._bridge_prep(0, True)
            if isinstance(prep, str):
                self.info = prep
                self.refresh()
                return
            self.drop_tool()
            self.push()
            rail_bridge2(self.cage, *prep, self.rail_segs)
        elif len(chains) == 1 and chains[0][1] and region:
            self.drop_tool()
            self.push()
            rail_bridge1(self.cage, region, chains[0][0], self.rail_segs)
        else:
            self.info = T("Rail bridge: select two polygon regions or edge loops, or one region to extend")
            self.refresh()
            return
        self.sel, self.sel_v, self.sel_e = set(), set(), set()
        self.last_op = None
        self.rail = dict(i=len(self.cage.rails) - 1, sel=None, span=None)
        self.info = T("Rail bridge: click a point and move it, I: add segment, Del: delete point, Enter: done")
        self.refresh()

    def _rail(self):
        """Duzenlenen kilavuz kaydi (yoksa None)."""
        if self.rail is None or self.rail["i"] >= len(self.cage.rails):
            return None
        return self.cage.rails[self.rail["i"]]

    @staticmethod
    def _rail_movable(r, N, j):
        """N'deki j dugumu tasinabilir mi: ara noktalar ve tek uclu kullanimda serbest uc; A ve B sabit."""
        return 0 < j < (len(N) if r.get("b") is None else len(N) - 1)

    def _rail_pick(self, r, N):
        """Secili kilavuz dugumleri (N sirasi): Shift ile secilenler ve tasinan secili nokta."""
        pick = [j for j in self.rail.get("pick", []) if j < len(N)]
        h = self.rail["sel"]
        if h is not None and h + 1 < len(N) and h + 1 not in pick:
            pick.append(h + 1)
        return pick

    def _rail_draw(self):
        dim, act = [], None
        if self.mode != "object":
            for ri, r in enumerate(self.cage.rails):
                if not rail_ok(self.cage, r):
                    continue
                N, M, _, _ = rail_nodes(self.cage, r)
                C = rail_curve(N, M)
                if self.rail is not None and ri == self.rail["i"]:
                    act = (r, N, C)
                else:
                    dim.append(C)
        idx, o = [], 0
        for C in dim:
            idx += list(range(o, o + len(C))) + [-1]
            o += len(C)
        _set_coords(self.rail_dcoord, np.vstack(dim) if dim else np.zeros((1, 3)))
        _set_index(self.rail_dlines, idx)
        H, S, span, C = np.zeros((0, 3)), np.zeros((0, 3)), [], np.zeros((1, 3))
        if act is not None:
            r, N, C = act
            H = N                                          # butun dugumler; A / B (sabit) yalniz hizalamada secilir
            h, sp = self.rail["sel"], self.rail["span"]
            pick = self._rail_pick(r, N)
            S = N[pick] if pick else S
            if sp is not None and sp < len(N) - 1:
                span = list(range(sp * 16, sp * 16 + 17)) + [-1]
        _set_coords(self.rail_coord, C)
        _set_index(self.rail_lines, list(range(len(C))) + [-1] if act is not None else [])
        _set_index(self.rail_span, span)
        _set_coords(self.rail_hcoord, H if len(H) else np.zeros((1, 3)))
        self.rail_hset.numPoints = len(H)
        _set_coords(self.rail_scoord, S if len(S) else np.zeros((1, 3)))
        self.rail_sset.numPoints = len(S)
        if self.rail is not None and not self.grab:        # gizmo secili (tasinabilir) kilavuz noktasinda
            h = self.rail["sel"] if act is not None else None
            self.gizmo.place(App.Vector(*self.l2w(N[h + 1])) if h is not None and h + 1 < len(N) else None)

    def _rail_curve_dist(self, C, pos):
        """Ekranda pos'un egriye (yerel noktalar) en yakin uzakligi (piksel) ve en yakin noktanin sirasi."""
        Q = self.project(C)
        p = np.array(pos[:2], float)
        best = (1e18, 0)
        for k in range(len(Q) - 1):
            a, b = Q[k], Q[k + 1]
            ab = b - a
            t = min(max((p - a).dot(ab) / max(ab.dot(ab), 1e-12), 0.0), 1.0)
            d = float(np.linalg.norm(a + t * ab - p))
            if d < best[0]:
                best = (d, k)
        return best

    def rail_pick(self, pos, radius=8.0):
        """Imlecin altindaki kilavuz (cift tiklama): sira numarasi ya da None."""
        best = None
        for ri, r in enumerate(self.cage.rails):
            if not rail_ok(self.cage, r):
                continue
            N, M, _, _ = rail_nodes(self.cage, r)
            d, _ = self._rail_curve_dist(rail_curve(N, M), pos)
            if d < radius and (best is None or d < best[0]):
                best = (d, ri)
        return None if best is None else best[1]

    def rail_edit(self, ri):
        """Kaydedilmis kilavuzu yeniden duzenlemeye al (cift tiklama)."""
        r = self.cage.rails[ri]
        if not rail_ok(self.cage, r):
            self.info = T("Rail bridge: the tube was changed, its guide curve can no longer be edited")
            self.refresh()
            return
        self.drop_tool()
        self.sel, self.sel_v, self.sel_e = set(), set(), set()
        self.last_op = None
        self.rail = dict(i=ri, sel=None, span=None)
        self.info = T("Rail bridge: click a point and move it, I: add segment, Del: delete point, Enter: done")
        self.refresh()

    def rail_click(self, pos):
        """Kilavuz duzenlenirken tiklama: ara nokta (tutamak) ya da egri parcasi secilir; bosluk secimi kaldirir."""
        r = self._rail()
        if r is None:
            self.rail = None
            self.refresh()
            return
        N, M, _, _ = rail_nodes(self.cage, r)
        shift = bool(QtWidgets.QApplication.keyboardModifiers() & QtCore.Qt.ShiftModifier)
        pick = self._rail_pick(r, N) if shift else []
        self.rail["span"] = None
        d = np.linalg.norm(self.project(N) - np.array(pos[:2], float), axis=1)
        j = int(np.argmin(d))
        if d[j] < 12.0:                                    # dugum: Shift ile secime ekle / cikar (hizalama icin)
            if shift and j in pick:
                pick.remove(j)
            else:
                pick.append(j)
            movable = [k for k in pick if self._rail_movable(r, N, k)]
            self.rail["sel"] = movable[-1] - 1 if movable else None
            self.rail["pick"] = pick
            self.refresh()
            return
        if shift:
            self.refresh()
            return
        self.rail["sel"], self.rail["pick"] = None, []
        d, k = self._rail_curve_dist(rail_curve(N, M), pos)
        if d < 8.0:
            self.rail["span"] = min(k // 16, len(N) - 2)   # rail_curve: parca basina 16 adim
        self.refresh()

    def _rail_change(self, fn, msg=None):
        """Kilavuzda bir degisiklik: geri alma kaydi, uygula, halkalari yerlestir. fn False donerse geri alinir."""
        r = self._rail()
        if r is None or self.drag or self.grab:
            return False
        self.push()
        if fn(r) is False:
            self.undo.pop()
            if msg:
                self.info = msg
            self.refresh()
            return False
        rail_place(self.cage, r)
        self.refresh()
        return True

    def rail_add(self):
        """Secili parcaya (yoksa secili noktadan sonraki, o da yoksa en uzun parcaya) ara halka ekle (I)."""
        r = self._rail()
        if r is None:
            return
        N, _, _, _ = rail_nodes(self.cage, r)
        sp, h = self.rail["span"], self.rail["sel"]
        if sp is None and h is not None:
            sp = min(h + 1, len(N) - 2)                    # secili noktadan sonraki parca (uc noktada: onceki)
        if sp is None:
            sp = int(np.argmax(np.linalg.norm(np.diff(N, axis=0), axis=1)))
        if self._rail_change(lambda r: rail_insert(self.cage, r, sp),
                             T("Rail bridge: the tube was changed, its guide curve can no longer be edited")):
            self.rail["sel"], self.rail["span"], self.rail["pick"] = sp, None, []   # yeni nokta secili: hemen tasinir

    def rail_delete(self):
        """Secili ara noktayi ve halkasini sil (Del); tek uclu kullanimda serbest uc silinmez."""
        r = self._rail()
        if r is None or self.rail["sel"] is None:
            return
        h = self.rail["sel"]
        if self._rail_change(lambda r: rail_remove(self.cage, r, h + 1),
                             T("Rail bridge: this point cannot be deleted")):
            self.rail["sel"], self.rail["pick"] = None, []

    def rail_segments(self, val):
        """Panel: ara halka sayisi; noktalar egri boyunca esit araliga dizilir (egrinin sekli korunur)."""
        k = min(max(int(round(val)), 0), 64)
        self.rail_segs = max(k, 1)
        self._rail_change(lambda r: rail_resample(self.cage, r, k))
        if self.rail is not None:
            self.rail["sel"] = self.rail["span"] = None
            self.rail["pick"] = []

    def rail_point_set(self, axis, val):
        """Panel: secili ara noktanin dunya X / Y / Z konumu (mm)."""
        h = self.rail["sel"] if self.rail is not None else None
        if h is None:
            return

        def set_(r):
            p = np.array(r["pts"][h], float)
            w = self.l2w(p)
            w[axis] = val
            r["pts"] = self._rail_move(r, r["pts"], self.w2l(w) - p).tolist()
        self._rail_change(set_)

    def _rail_move(self, r, pts, d):
        """Secili ara nokta d (yerel) kadar tasininca butun ara noktalar: Move neighbors kapaliysa yalniz o nokta,
        aciksa kilavuz kablo gibi davranir (rail_cable, Elasticity). Uc halkalar (A, B) hic oynamaz."""
        P = np.array(pts, float)
        h = self.rail["sel"]
        if h is None or h >= len(P):
            return P
        if not self.rail_falloff:
            P[h] += d
            return P
        N, _, fa, fb = rail_nodes(self.cage, r)
        N = N.copy()
        N[1:1 + len(P)] = P                                # surukleme baslangicindaki dugumler
        return rail_cable(N, h + 1, d, fb is not None, self.rail_elastic, fa[3],
                          fb[3] if fb is not None else None)[1:1 + len(P)]

    def rail_falloff_option(self, on):
        self.rail_falloff = bool(on)
        self.refresh()

    def rail_tension_set(self, val):
        """Panel: egrinin gerginligi (0 yumusak .. 1 noktalar arasi duz), bu kilavuzda saklanir."""
        v = min(max(float(val), 0.0), 1.0)
        self._rail_change(lambda r: r.__setitem__("tension", v))

    def rail_elastic_set(self, val):
        self.rail_elastic = min(max(float(val), 0.0), 1.0)
        self.refresh()

    def rail_scale(self, val):
        """Panel: secili noktadaki halkanin olcegi (1 = uc halkalarin boyu)."""
        h = self.rail["sel"] if self.rail is not None else None
        if h is None:
            return
        self._rail_change(lambda r: r["sc"].__setitem__(h, max(float(val), 0.01)))

    def rail_finish(self):
        self.rail = None
        self.info = ""
        self.refresh()

    def _rail_panel(self):
        r = self._rail()
        if r is None:
            self.rail = None
            return None
        k = len(r["rings"]) - (1 if r.get("b") is None else 0)
        fields = [("Segments", float(k), "", 0, self.rail_segments),
                  ("Tension", float(r.get("tension", 0.0)), "", 2, self.rail_tension_set)]
        h = self.rail["sel"]
        if h is not None and h < len(r["pts"]):
            w = self.l2w(np.array(r["pts"][h], float))
            fields += [(a, float(w[j]), " mm", 3, (lambda val, j=j: self.rail_point_set(j, val))) for j, a in enumerate("XYZ")]
            fields.append(("Scale", float(r["sc"][h]), " x", 3, self.rail_scale))
        if self.rail_falloff:
            fields.append(("Elasticity", float(self.rail_elastic), "", 2, self.rail_elastic_set))
        checks = [("Move neighbors", self.rail_falloff, self.rail_falloff_option)]
        buttons = [("Add segment (I)", self.rail_add), ("Delete point (Del)", self.rail_delete),
                   ("Done (Enter)", self.rail_finish),
                   ("Align horizontally (Shift-click guide points, the end points too)",
                    lambda: self.align_points(False), self._align_icon(False)),
                   ("Align vertically (Shift-click guide points, the end points too)",
                    lambda: self.align_points(True), self._align_icon(True))]
        return T("Rail bridge"), fields, self.drag is None, checks, None, buttons

    # --- From face (alt menu: Face / Edge / Point): FreeCAD'de secili Part / PartDesign elemanlarindan kafes ---
    # Object modunda yeni SubD nesnesi olur (dunya koordinati), kafeste duzenlerken ayni kafese eklenir.
    # Face (cage_from_faces): yuzeyler; yuz secili degilse kapali duz kenar halkalari duz yuz olur. Panel: Arc
    # Segments (90 derece basina), Max Length (mm, 0 = otomatik), Exact Edge (1: SubD kenari yuzey sinirindan tam
    # gecer, 0: kafes noktalari sinirda). Edge (cage_from_edges): zincir basina bir cokgen, Arc Segments / Max
    # Length. Point (cage_from_points): secim sirasiyla tek cokgen. Edge / Point cokgenleri kameraya bakar.
    FACE = ("Arc Segments", "Max Length", "Exact Edge")

    def _selected_part(self, names=None):
        """FreeCAD secimindeki Part / PartDesign alt elemanlari, dunya koordinatinda (Body yerlesimi dahil), secim
        sirasiyla: {"Face": [...], "Edge": [...], "Vertex": [...]}. names (sozluk): ayni sirayla (nesne, alt eleman)
        adlari da doldurulur (Sketch bagi)."""
        import Part
        items = [(so.Object, sub) for so in Gui.Selection.getSelectionEx("", 0) if not is_subd(so.Object)
                 for sub in so.SubElementNames]
        order = {(o, s): k for k, (d, o, s) in enumerate(_pick_order.picks)}
        items.sort(key=lambda it: order.get((it[0].Name, it[1]), len(order)))
        out = {"Face": [], "Edge": [], "Vertex": []}
        for obj, sub in items:
            try:
                sh = Part.getShape(obj, sub, needSubElement=True)
            except Exception:
                continue
            if not sh.isNull() and sh.ShapeType in out:
                out[sh.ShapeType].append(sh)
                if names is not None:
                    names.setdefault(sh.ShapeType, []).append((obj.Name, sub))
        return out

    def _selected_part_faces(self):
        """Secili yuzeyler; yuz yoksa secili kapali kenar halkalarindan duz yuz (ic ice iki halka: halka yuz).
        Kopya alinmaz: ortak kenarlar (isSame) taninsin. Hata metni dondurebilir."""
        import Part
        sel = self._selected_part()
        faces, edges = sel["Face"], [] if sel["Face"] else sel["Edge"]
        if edges:
            groups = Part.sortEdges(edges)
            wires = [Part.Wire(g) for g in groups]
            if any(not w.isClosed() for w in wires):
                return T("From face: the selected edges must form closed loops")
            try:
                f = Part.makeFace(wires, "Part::FaceMakerBullseye")
                faces += list(f.Faces)
            except Exception:
                return T("From face: no flat face can be made from the selected edges")
        return faces

    def from_face(self):
        if self.drag or self.grab:
            return
        faces = self._selected_part_faces()
        if not faces or isinstance(faces, str):
            self.info = faces or T("From face: select faces or a closed edge loop of a Part / PartDesign object")
            self.refresh()
            return
        self._create_from("fromface", {"Arc Segments": 4.0, "Max Length": 0.0, "Exact Edge": 1.0}, faces=faces)

    def from_edges(self):
        if self.drag or self.grab:
            return
        names = {}
        edges = self._selected_part(names)["Edge"]
        if not edges:
            self.info = T("From face: select edges of a Part / PartDesign object")
            self.refresh()
            return
        self._create_from("fromedge", {"Arc Segments": 4.0, "Max Length": 0.0}, edges=edges,
                          look=self._look(), src=names.get("Edge", []))

    def from_points(self):
        if self.drag or self.grab:
            return
        pts = [np.array(tuple(v.Point)) for v in self._selected_part()["Vertex"]]
        if not pts:
            self.info = T("From face: select 3 or more points of a Part / PartDesign object")
            self.refresh()
            return
        self._create_from("frompoint", {}, points=pts, look=self._look())

    def _look(self):
        return np.array(tuple(self.view.getViewDirection()))   # dunya bakis yonu: islem tekrarinda da ayni kalsin

    def _create_from(self, kind, vals, **args):
        test = self._from_cage(kind, args, dict(vals))
        if isinstance(test, str):
            self.info = test
            self.refresh()
            return
        if self.mode == "object":                          # yeni nesne (Paste as object gibi), hemen kafesine gir
            doc = App.ActiveDocument
            doc.openTransaction("SubD from face")
            o = doc.addObject("Mesh::FeaturePython", "SubD")
            SubDObject(o)
            if App.GuiUp:
                SubDViewProvider(o.ViewObject)
            cage_to_obj(test, o)
            set_surface(o, self.levels, self.rules)
            o.recompute()
            doc.commitTransaction()
            self.edit_mode = "face"
            self.edit_object(o)
            self.push()
            op = self.new_op(kind, vals, append=False, **args)
        else:
            self.drop_tool()
            self.push()
            op = self.new_op(kind, vals, append=True, **args)
        self._apply_op(op)
        self.last_op, self.op = op, None
        self.refresh()

    def _from_cage(self, kind, a, v):
        """From face islemlerinin kafesi (dunya) ya da hata metni; Edge / Point cokgenleri a["look"]a bakar."""
        if kind == "fromface":
            return self._face_cage(a["faces"], v)
        if kind == "frompoint":
            c = cage_from_points(a["points"])
        else:
            segs = min(max(int(round(v["Arc Segments"])), 1), 32)
            v["Arc Segments"] = float(segs)
            if v["Max Length"] <= 0:                       # otomatik: panelde hesaplanan deger gorunsun
                v["Max Length"] = float(auto_max_len(a["edges"], segs))
            c = cage_from_edges(a["edges"], segs, v["Max Length"])
        if not isinstance(c, str):
            for i, f in enumerate(c.F):
                P = c.V[f]
                n = np.cross(P - P.mean(0), np.roll(P - P.mean(0), -1, 0)).sum(0)
                if n.dot(a["look"]) > 0:
                    c.F[i] = f[::-1]
        return c

    def _face_cage(self, faces, v):
        """Yuzeylerden (dunya) tek kafes ya da hata metni; degerleri sinirlar, otomatik Max Length'i yazar."""
        segs = min(max(int(round(v["Arc Segments"])), 1), 32)
        v["Arc Segments"] = float(segs)
        v["Exact Edge"] = 1.0 if v["Exact Edge"] >= 0.5 else 0.0
        if v["Max Length"] <= 0:                           # otomatik: panelde hesaplanan deger gorunsun
            v["Max Length"] = float(auto_max_len([e for f in faces for e in f.Edges], segs))
        return cage_from_faces(faces, segs, v["Max Length"], v["Exact Edge"] > 0)

    # --- Make polygon (LightWave Make Polygon / Blender F): secimden yeni yuz ---
    # Nokta modu: secili noktalar (aralarindaki kenarlarin sirasi, yoksa duzlemlerinde aciya gore). Kenar modu:
    # kapali halka -> tek cokgen; iki ayri kenar / acik zincir -> uclari birlestiren dortgen / cokgen (capraz
    # olmayan esleme). Yon komsu yuzlerle ayni (ortak kenari ters sirada gecer), komsu yoksa kameraya bakar.

    def make_polygon(self):
        if self.drag or self.grab or self.mode not in ("point", "edge"):
            if self.mode == "face":
                self.info = T("Make polygon: select points or edges")
                self.refresh()
            return
        res = self._polygon_from_selection()
        if isinstance(res, str):
            self.info = res
            self.refresh()
            return
        cage = self.cage
        ef = cage.edge_faces()
        keys = {frozenset(f) for f in cage.F}
        if frozenset(res) in keys:
            self.info = T("Make polygon: this polygon already exists")
            self.refresh()
            return
        edges = [_edge(res[i], res[(i + 1) % len(res)]) for i in range(len(res))]
        if any(len(ef.get(e, ())) >= 2 for e in edges):
            self.info = T("Make polygon: an edge already has two polygons")
            self.refresh()
            return
        # yon: komsu yuz ortak kenari a -> b geciyorsa yeni yuz b -> a gecmeli
        vote = 0
        for i, e in enumerate(edges):
            a, b = res[i], res[(i + 1) % len(res)]
            for fi in ef.get(e, ()):
                f = cage.F[fi]
                same = any(f[k] == a and f[(k + 1) % len(f)] == b for k in range(len(f)))
                vote += -1 if same else 1
        if vote == 0:                                      # komsu yok: kameraya baksin
            P = cage.V[res]
            n = np.cross(P - P.mean(0), np.roll(P - P.mean(0), -1, 0)).sum(0)
            vote = 1 if n.dot(-self.view_dir()) >= 0 else -1
        if vote < 0:
            res = res[::-1]
        self.drop_tool()
        self.push()
        cage.F.append(list(res))
        self.last_op = None
        self.info = T("Make polygon: %d sides") % len(res)
        self.refresh()

    def _polygon_from_selection(self):
        """Secimden sirali nokta listesi ya da hata metni."""
        cage = self.cage
        err = T("Make polygon: select 3 or more points, a closed edge loop, or two edges")
        if self.mode == "edge":
            chains = cage.edge_chains(self.sel_e) if self.sel_e else None
            if not chains:
                return err
            if len(chains) == 1:
                verts = chains[0][0]
            elif len(chains) == 2 and not chains[0][1] and not chains[1][1]:
                A, B = chains[0][0], chains[1][0]           # iki acik zincir: capraz olmayan baglanti
                V = cage.V
                fwd = np.linalg.norm(V[A[-1]] - V[B[-1]]) + np.linalg.norm(V[B[0]] - V[A[0]])
                rev = np.linalg.norm(V[A[-1]] - V[B[0]]) + np.linalg.norm(V[B[-1]] - V[A[0]])
                verts = A + (B[::-1] if fwd <= rev else B)
            else:
                return err
        else:
            verts = sorted(self.sel_v)
            if len(verts) < 3:
                return err
            # aralarindaki kenarlar tek zincir olusturuyorsa o sira; yoksa ortak duzlemde aciya gore
            inner = {e for e in cage.edges() if e <= self.sel_v}
            chains = cage.edge_chains(inner) if inner else None
            if chains and len(chains) == 1 and len(chains[0][0]) == len(verts):
                verts = chains[0][0]
            else:
                P = cage.V[verts]
                c = P.mean(0)
                _, _, vt = np.linalg.svd(P - c)
                u, w = vt[0], vt[1]
                ang = np.arctan2((P - c) @ w, (P - c) @ u)
                verts = [verts[i] for i in np.argsort(ang)]
        if len(verts) < 3 or len(set(verts)) != len(verts):
            return err
        return list(verts)

    # --- Hole: secili yuzlerin (1, 2, 4...) ortasina daire delik; panel: Radius, Segments, Support ---
    HOLE = ("Radius", "Segments", "Support")

    SHELL_MODES = ("Inside", "Outside", "Centered")

    def shell(self):
        """Modify > Shell: Object modunda secili nesnenin tamami, Edit modunda secili cokgenler (secim yoksa hepsi)
        kabuklanir. Object modunda nesnenin kafesine girilir ki panelden kalinlik / suzgec ayarlanabilsin."""
        if self.drag or self.grab:
            return
        if self.mode == "object":
            objs = self.selected_objects()
            if not objs:
                self.info = T("Shell: select a SubD object or polygons")
                self.refresh()
                return
            self.edit_mode = "face"
            self.edit_object(objs[-1])
            self.sel = set(range(len(self.cage.F)))
        elif self.mode != "face":
            self.set_mode("face")
        region = set(self.sel) if self.sel else set(range(len(self.cage.F)))
        self.drop_tool()
        self.push()
        op = self.new_op("shell", dict(_last_shell), faces=sorted(region))
        err = self._shell(op["values"], op["args"]["faces"])
        if err:
            self.op = None
            self.do_undo(False)
            self.info = err
        else:
            self.last_op, self.op = op, None
        self.refresh()

    def _shell(self, v, faces):
        v["Thickness"] = max(float(v["Thickness"]), 0.001)
        v["Detail filter"] = max(float(v["Detail filter"]), 0.0)
        res = self.cage.shell(faces, v["Thickness"], int(v["Mode"]), v["Detail filter"], bool(v["Fill"]),
                              bool(v["Sharp"]), self.rules, bool(v.get("Flat rim", 1.0)))
        if isinstance(res, str):
            return res
        _last_shell.update(v)                              # yalniz basarili degerler sonraki kabuga gecer
        inner, open_inner, thin = res
        self.sel, self.sel_v, self.sel_e = set(inner), set(), set(open_inner)
        self.info = T("Shell: %d polygons, %d open edges to weave") % (len(inner), len(open_inner)) + \
            (T("  |  %d inner points merged / shortened (fine details, crossings)") % thin if thin else "")
        return None

    def shell_option(self, name, on):
        """Shell panelindeki kutucuk: yon (Inside / Outside / Centered tek secim), Fill, Sharp; islem yeniden yapilir."""
        op = self.last_op
        if op is None or op["kind"] != "shell":
            return
        v = op["values"]
        if name in self.SHELL_MODES:
            v["Mode"] = float(self.SHELL_MODES.index(name)) if on else 0.0
        else:
            v[name] = 1.0 if on else 0.0
        self._reapply(op)

    def hole(self, cap=False):
        if self.drag or self.grab or self.mode != "face":
            return
        self.drop_tool()
        self.push()
        op = self.new_op("hole", dict(Radius=0.0, Segments=0.0, Support=-1.0), cap=cap)   # 0 / -1: otomatik
        err = self._hole(op["values"], cap)
        if err:
            self.op = None
            self.do_undo(False)
            self.info = err
        else:
            self.last_op, self.op = op, None
        self.refresh()

    def _hole(self, v, cap=False):
        res = self.cage.hole(set(self.sel), v["Radius"], int(round(v["Segments"])), v["Support"], cap)
        if isinstance(res, str):
            return res
        new, v["Radius"], segs, v["Support"] = res
        v["Segments"] = float(segs)
        self.sel, self.sel_v, self.sel_e = set(new), set(), set()
        self.info = T("Circle on face: radius %.2f mm, %d segments" if cap else "Hole: radius %.2f mm, %d segments") % (
            v["Radius"], segs)
        return None

    def object_paste(self):
        doc = App.ActiveDocument
        doc.openTransaction("SubD paste")
        o = doc.addObject("Mesh::FeaturePython", "SubD")
        SubDObject(o)
        if App.GuiUp:
            SubDViewProvider(o.ViewObject)
        cage_to_obj(_clipboard, o)                         # dunya koordinati: Placement birim
        set_surface(o, self.levels, self.rules)
        o.recompute()
        doc.commitTransaction()
        c = _clipboard                                     # yuzsuz parca kenar / nokta modunda acilir
        self.edit_mode = "face" if c.F else ("edge" if c.W else "point")
        self.edit_object(o)
        self._select_new(range(len(self.cage.F)), set(range(len(self.cage.V))))
        self.info = T("Pasted as a new object: %d polygons, %d edges") % (len(self.cage.F), len(self.cage.W))
        self.refresh()

    # --- tasima (grab) ---
    def start_grab(self, verts, axis=None, pushed=False, lazy=False, op=None):
        if not verts or self.drag or self.grab:
            return
        moves = {}
        if op is None:                                     # duz tasima: kilitli grubun uyesi tum grubu surer
            verts, moves = self._group_expand(verts)
        elif self._locked_block(verts):
            return
        self.drop_tool()
        if not pushed:
            self.push()
        if op is None:
            op = self.new_op("move", dict(X=0.0, Y=0.0, Z=0.0), verts=list(verts), grp=moves)
        V = self.cage.V
        c = V[verts].mean(0)
        vd = self.view_dir()
        n = vd
        if axis is not None:
            n = np.cross(axis, np.cross(vd, axis))
            n = n / np.linalg.norm(n) if np.linalg.norm(n) > 1e-6 else vd
        p0 = self.ray_plane(self.mouse, c, n)
        if p0 is None and not lazy:
            return
        self.grab = dict(state=self.move_state(verts), c=c, n=n, p0=p0, axis=axis, lazy=lazy,
                         anchor=self._snap_prepare(verts))
        self.gizmo.place(None)

    # --- Snap: serbest tasimada (G, gizmo merkez kupu) gorunur Part / PartDesign sekillerine yakalama ---
    # Hedef: imlece 12 piksel icindeki kose noktasi, yoksa 10 piksel icindeki kenar noktasi. Secimin imlece en yakin
    # noktasi (capa) hedefe tam oturur. Hedefler surukleme basinda toplanir (kamera surukleme boyunca sabit).

    def toggle_snap(self):
        self.snap = not self.snap
        if not self.snap:
            self._line_mark(None)
        self.refresh()

    def _snap_prepare(self, verts):
        """Snap aciksa hedefleri topla ve capayi (secimin imlece en yakin noktasi, yerel) dondur; kapaliysa None."""
        self._snapdata = None
        if not self.snap or not len(verts):
            return None
        import Part
        pts, A, B = [], [], []
        for o in App.ActiveDocument.Objects:
            if (is_subd(o) or "Samples" in o.PropertiesList or not hasattr(o, "Shape")
                    or o.TypeId.startswith("App::")):              # Origin eksen / duzlemleri, gruplar
                continue
            try:
                if not o.ViewObject.Visibility:
                    continue
                sh = Part.getShape(o)
            except Exception:
                continue
            if sh.isNull():
                continue
            pts += [tuple(x.Point) for x in sh.Vertexes]
            for e in sh.Edges:
                q = [tuple(x) for x in e.discretize(Number=2 if _edge_turn(e) < 1e-3 else 33)]
                A += q[:-1]
                B += q[1:]
        vv, W, H = _vv(self.view)
        scr = lambda P: np.array([_to_screen(vv, W, H, App.Vector(*p)) for p in P], float).reshape(-1, 2)
        P, A, B = (np.array(x, float).reshape(-1, 3) for x in (pts, A, B))
        self._snapdata = dict(P=P, Ps=scr(P), A=A, B=B, As=scr(A), Bs=scr(B))
        loc = self.cage.V[list(verts)]
        S = scr(self.l2w(loc))
        return loc[int(np.argmin(np.linalg.norm(S - np.array(self.mouse, float), axis=1)))].copy()

    def _snap_target(self, pos):
        """Imlec altindaki yakalama hedefi: (dunya noktasi, "vertex" / "edge") ya da None."""
        d = getattr(self, "_snapdata", None)
        if not d:
            return None
        p = np.array(pos, float)
        if len(d["Ps"]):
            dist = np.linalg.norm(d["Ps"] - p, axis=1)
            i = int(np.argmin(dist))
            if dist[i] < 12:
                return d["P"][i], "vertex"
        if len(d["As"]):
            ab = d["Bs"] - d["As"]
            t = np.clip(((p - d["As"]) * ab).sum(1) / np.maximum((ab * ab).sum(1), 1e-9), 0, 1)
            dist = np.linalg.norm(d["As"] + t[:, None] * ab - p, axis=1)
            i = int(np.argmin(dist))
            if dist[i] < 10:
                return d["A"][i] + t[i] * (d["B"][i] - d["A"][i]), "edge"
        return None

    def move_state(self, verts):
        """Tasinacak noktalar ve baslangic konumlari; simetride aynadaki esler ve duzlemdeki noktalar da."""
        V = self.cage.V
        verts = list(verts)
        on_plane, mirrored = [], []
        if self.sym:
            M, chosen = self.mirror_map(), set(verts)
            on_plane = [v for v in verts if M.get(v) == v]
            mirrored = sorted({M[v] for v in verts if M.get(v) not in (None, v)} - chosen)
        return dict(verts=verts, orig=V[verts].copy(), on_plane=on_plane, plane_x=V[on_plane, 0].copy(),
                    mirrored=mirrored, morig=V[mirrored].copy())

    def apply_move(self, st, d):
        self.apply_transform(st, lambda P: P + d, lambda P: P + d * MIRROR)

    def apply_rotation(self, st, c, axis, ang):
        """Secimi c merkezi etrafinda dondurur; aynadaki esler aynalanmis donusle (S R S) doner."""
        R = _rotation_matrix(axis, ang)
        Rm = R * np.outer(MIRROR, MIRROR)
        cm = c * MIRROR
        self.apply_transform(st, lambda P: c + (P - c) @ R.T, lambda P: cm + (P - cm) @ Rm.T)

    def apply_scale(self, st, c, kk):
        """Stretch / Size: c (yerel) etrafinda, DUNYA eksenlerinde kk olcegi; aynadaki esler aynalanmis matrisle."""
        A = self.world_matrix(kk)
        Am = A * np.outer(MIRROR, MIRROR)
        cm = c * MIRROR
        self.apply_transform(st, lambda P: c + (P - c) @ A.T, lambda P: cm + (P - cm) @ Am.T)

    def apply_transform(self, st, f, fm):
        V = self.cage.V
        V[st["verts"]] = f(st["orig"])
        V[st["on_plane"], 0] = st["plane_x"]              # simetri duzlemindeki noktalar duzlemde kalir
        V[st["mirrored"]] = fm(st["morig"])               # aynadaki es noktalar
        self.refresh(fast=True)

    def grab_move(self):
        g = self.grab
        if g is None:
            return
        if g.get("kind") == "magnet":
            self.magnet_move()
            return
        if g.get("kind") == "bevel":
            self.bevel_move()
            return
        if g.get("kind") == "slide":
            self.slide_move()
            return
        if g.get("kind") == "ebevel":
            self.ebevel_move()
            return
        if g.get("kind") == "extrude":
            self.extrude_move()
            return
        if g.get("kind") == "twist":
            self.twist_move()
            return
        if g.get("lazy"):                                  # düğmeden başladı: ölçü buradan başlar
            g["p0"] = self.ray_plane(self.mouse, g["c"], g["n"])
            g["lazy"] = g["p0"] is None
            return
        p = self.ray_plane(self.mouse, g["c"], g["n"])
        if p is None:
            return
        d = p - g["p0"]
        step = self.step_box.value() if _ctrl_held() else 0.0
        vals = self.op["values"] if self.op is not None else {}
        if g["axis"] is not None:
            t = _snap(d.dot(g["axis"]), step)
            d = g["axis"] * t
            self.info = T("Distance %.3f mm") % t
            vals["Distance"] = t
        else:
            d = np.array([_snap(x, step) for x in d])
        hit = self._snap_target(self.mouse) if g["axis"] is None and g.get("anchor") is not None else None
        if hit is not None:                                # capa noktasi hedefe tam otursun
            d = self.w2l(hit[0]) - g["anchor"]
        dw = self.l2w_dir(d)                               # panelde ve gostergede dunya eksenleri
        if g["axis"] is None:
            self.info = "Δ  X %.3f   Y %.3f   Z %.3f mm" % tuple(dw) + (
                "   " + T("Snap: %s") % T(hit[1]) if hit is not None else "")
        vals.update(X=float(dw[0]), Y=float(dw[1]), Z=float(dw[2]))
        self.apply_move(g["state"], d)

    def end_grab(self, ok):
        g, self.grab = self.grab, None
        if g is not None and g.get("kind") == "slide":
            if not ok:                                     # Esc / sag tik: kesim kalir, ortaya doner
                self.slide_set(g["items"], g["t0"])
            ok = True
        if ok:
            self.last_op, self.op = self.op, None
            self._group_follow(((self.last_op or {}).get("args") or {}).get("grp"))
            self.refresh()
        else:
            self.op = None
            self.do_undo(False)

    # --- gizmo ile tasima ---
    def drag_begin(self, key, pos):
        if self.mode == "object":                          # secili nesneleri tasi / dondur / olcekle
            objs = self.selected_objects()
            c = self.objects_center(objs)
            if not objs or c is None:
                return
            objs[0].Document.openTransaction("SubD objects")
            self.info = ""
            self.op = dict(kind="object", values=dict(X=0.0, Y=0.0, Z=0.0, Angle=0.0), args={}, before=None)
            c0 = App.Vector(*c)
            self.drag = dict(key=key, objs=objs, pl0=[App.Placement(o.Placement) for o in objs],
                             cage0=[cage_from_obj(o) for o in objs], c0=c0, c=c,
                             free_n=_view_dir(self.view, c0), acc=0.0, vec=App.Vector())
            self.drag["ref"] = self.drag["ref0"] = self.drag_measure(pos)
            return
        if self.rail is not None:                          # kilavuzun secili ara noktasi tasinir
            r = self._rail()
            if r is None or self.rail["sel"] is None or key[0] not in ("move", "plane", "free"):
                return
            self.push()
            p = np.array(r["pts"][self.rail["sel"]], float)
            c0 = App.Vector(*self.l2w(p))
            self.op = dict(kind="railpoint", values=dict(X=0.0, Y=0.0, Z=0.0), args={}, before=None)
            self.drag = dict(key=key, rail=np.array(r["pts"], float), c0=c0, c=p, free_n=_view_dir(self.view, c0),
                             acc=0.0, vec=App.Vector(), anchor=None)
            self.drag["ref"] = self.drag["ref0"] = self.drag_measure(pos)
            return
        verts = self.selected_verts()
        kind, i = key
        if not verts:
            return
        moves = {}
        if kind in ("stretch", "size"):                    # olcek kilitli grubu bozar
            if self._locked_block(verts):
                return
        else:                                              # tasima / donme: kilitli grubun uyesi tum grubu surer
            verts, moves = self._group_expand(verts)
        self.push()
        self.info = ""
        c = self.cage.V[verts].mean(0)                     # yerel merkez (islem bununla yeniden yapilir)
        c0 = App.Vector(*self.l2w(c))                      # gizmo hesaplari dunya koordinatinda
        if kind == "rotate":
            self.new_op("rotate", dict(Angle=0.0), verts=verts, c=c, axis=tuple(AXES[i]), letter="XYZ"[i], grp=moves)
        elif kind in ("stretch", "size"):
            self.new_op("scale", dict(X=1.0, Y=1.0, Z=1.0), verts=verts, c=c)
        else:
            self.new_op("move", dict(X=0.0, Y=0.0, Z=0.0), verts=verts, grp=moves)
        self.drag = dict(key=key, state=self.move_state(verts), c0=c0, c=c, free_n=_view_dir(self.view, c0),
                         acc=0.0, vec=App.Vector(),
                         anchor=self._snap_prepare(verts) if kind == "free" else None)
        self.drag["ref"] = self.drag["ref0"] = self.drag_measure(pos)

    def drag_measure(self, pos):
        """Ok -> eksen parametresi, duzlem karesi / merkez kup -> duzlemdeki nokta."""
        g = self.drag
        o, d = _ray(self.view, pos)
        kind, i = g["key"]
        if kind in ("move", "stretch"):
            return _line_param(o, d, g["c0"], AXES[i])
        if kind == "size":
            return tuple(pos)
        if kind == "plane":
            return _plane_hit(o, d, g["c0"], AXES[i])
        if kind == "rotate":
            n = AXES[i]
            hit = _plane_hit(o, d, g["c0"], n)
            if hit is None or abs(d.dot(n)) < 0.08:        # halka yandan gorunuyor: ekranda merkez etrafindaki aci
                vv, W, H = _vv(self.view)
                c = _to_screen(vv, W, H, g["c0"])
                return ("screen", math.atan2(pos[1] - c[1], pos[0] - c[0]), d.dot(n))
            return ("plane", hit - g["c0"])
        return _plane_hit(o, d, g["c0"], g["free_n"])

    def drag_apply(self):
        self._pending = False
        g = self.drag
        if g is None:
            return
        cur = self.drag_measure(self.mouse)
        if cur is None or g["ref"] is None:
            return
        kind, i = g["key"]
        if kind == "rotate":
            r = g["ref"]
            g["ref"] = cur
            if r[0] != cur[0]:                             # olcu turu degisti (halka yandan gorunur oldu): yeni referans
                return
            n = AXES[i]
            if r[0] == "plane":
                a, b = r[1], cur[1]
                ang = math.atan2(n.dot(a.cross(b)), a.dot(b))
            else:
                ang = (cur[1] - r[1] + math.pi) % (2 * math.pi) - math.pi
                if r[2] > 0:
                    ang = -ang
            g["acc"] += ang * _gain()
            deg = math.degrees(g["acc"])
            if _ctrl_held():                               # Ctrl: aci adimina yakala
                deg = _snap(deg, self.angle_box.value())
            self.info = T("Angle %s  %.2f°") % ("XYZ"[i], deg)
            self.op["values"]["Angle"] = deg
            if self.mode == "object":
                self.objects_rotate(g, n, deg)
            else:
                self.apply_rotation(g["state"], g["c"], self.w2l_dir(np.array(tuple(n))), math.radians(deg))
            return
        if kind in ("stretch", "size"):
            # Stretch: imlecin eksen uzerindeki uzakligi / tutuldugu andaki uzaklik; Size: saga 150 px = 2 kat
            if kind == "stretch":
                if abs(g["ref0"]) < 1e-9:
                    return
                f = cur / g["ref0"]
            else:
                f = 1.0 + (cur[0] - g["ref0"][0]) / 150.0
            if QtWidgets.QApplication.keyboardModifiers() & QtCore.Qt.ShiftModifier:
                f = 1.0 + (f - 1.0) * 0.25                 # Shift: ince ayar
            if _ctrl_held():
                f = _snap(f, 0.05)                         # Ctrl: %5 adim
            f = max(f, 0.001)
            kk = np.ones(3)
            if kind == "stretch":
                kk[i] = f
            else:
                kk[:] = f
            self.op["values"].update(X=kk[0], Y=kk[1], Z=kk[2])
            self.info = ("Stretch %s  x%.4f" % ("XYZ"[i], f)) if kind == "stretch" else "Size  x%.4f" % f
            if self.mode == "object":
                self.objects_scale(g, kk)
            else:
                self.apply_scale(g["state"], g["c"], kk)
            return
        # her adimin farki o anki oranla (direnc; Shift = ince) carpilip toplanir: Shift ortada basilinca ziplamaz.
        # Ctrl: toplam tasima mm adimina yakalanir (ham deger birikmeye devam eder, Ctrl birakilinca kaldigi yerden).
        step = self.step_box.value() if _ctrl_held() else 0.0
        if kind == "move":
            g["acc"] += (cur - g["ref"]) * _gain()
            d = np.array(tuple(AXES[i])) * _snap(g["acc"], step)
        else:
            g["vec"] += (cur - g["ref"]) * (1.0 if kind == "free" else _gain())
            d = np.array([_snap(x, step) for x in g["vec"]])
        g["ref"] = cur
        hit = self._snap_target(self.mouse) if kind == "free" and g.get("anchor") is not None else None
        if hit is not None:                                # capa noktasi hedefe tam otursun (dunya)
            d = hit[0] - self.l2w(g["anchor"])
        self.info = "Δ  X %.3f   Y %.3f   Z %.3f mm" % tuple(d) + (
            "   " + T("Snap: %s") % T(hit[1]) if hit is not None else "")
        self.op["values"].update(X=float(d[0]), Y=float(d[1]), Z=float(d[2]))
        if g.get("rail") is not None:
            r = self._rail()
            if r is not None:
                r["pts"] = self._rail_move(r, g["rail"], self.w2l_dir(np.asarray(d, float))).tolist()
                rail_place(self.cage, r)
                self.refresh(fast=True)
            return
        if self.mode == "object":
            self.objects_move(g, d)
        else:
            self.apply_move(g["state"], self.w2l_dir(d))       # d dunya ekseninde, kafes yerel

    def drag_end(self):
        if self.drag is None:
            return
        self.drag_apply()
        self.gizmo.set_hover(None)
        g, self.drag = self.drag, None
        if self.mode == "object":                          # nesne donusumu: tek FreeCAD geri alma adimi
            doc = g["objs"][0].Document
            doc.recompute()
            doc.commitTransaction()
            self.op = None
            self.refresh()
            return
        if g.get("rail") is not None:                      # kilavuz noktasi: panel kilavuzu gostermeye devam eder
            self.op = None
            self.refresh()
            return
        self.last_op, self.op = self.op, None
        self._group_follow(((self.last_op or {}).get("args") or {}).get("grp"))
        self.refresh()

    # --- olaylar (Coin olayi icinde sahneye dokunulmaz; isler Qt zamanlayicisiyla olaydan sonra yapilir) ---
    def _schedule(self, fn):
        if not self._pending:
            self._pending = True
            _later(fn)

    def _on_button(self, cb):
        ev = cb.getEvent()
        pos = ev.getPosition().getValue()
        self.mouse = pos
        down = ev.getState() == coin.SoButtonEvent.DOWN
        left = ev.getButton() == coin.SoMouseButtonEvent.BUTTON1
        if self.drag is not None:
            if left and not down:
                cb.setHandled()
                _later(self.drag_end)
            return
        if self.box is not None and left and not down:
            cb.setHandled()
            self.finish_box(pos)
            return
        if left and not down and self.grab and self.grab.get("kind") == "magnet":
            cb.setHandled()
            _later(self.magnet_end)
            return
        if not down:
            return
        if self.grab:
            cb.setHandled()
            _later(lambda ok=left: self.end_grab(ok))
            return
        if not left:
            if self.armed == "magnet":                     # sag tusu _KeyFilter aliyor (FreeCAD menusu acilmasin)
                return
            if self.armed == "dim" and ev.getButton() == coin.SoMouseButtonEvent.BUTTON2:
                _later(self.dim_escape)                    # sag tik: bekleyen nokta varsa vazgec, yoksa kapat
                return
            if self.armed and ev.getButton() == coin.SoMouseButtonEvent.BUTTON2:   # sag tik: kalici araci kapat
                self.drop_tool()
                _later(self.refresh)
            return
        if self.dim_hits and self.dim_p0 is None and self.mode != "object":
            hd = self._dim_hit_label(pos)
            if hd is not None:                             # grup olcusunun degeri: tiklayinca degistirilir
                cb.setHandled()
                _later(lambda: self.handle_dim_edit(hd))
                return
            i = self._dim_hit_x(pos)
            if i is not None:                              # Distance olcusunun x dugmesi: olcuyu sil
                cb.setHandled()
                _later(lambda: self._dim_remove(i))
                return
        # kalici arac (sari dugme) acik kaldikca her tiklama onu uygular; secim / kutu secimi yapilmaz
        if self.armed in ("line", "mline", "array"):       # knife / ayna / array ekseni cizgisi baslar
            cb.setHandled()
            self.box = dict(p0=pos, press=pos, moved=False, kind="line", tool=self.armed, shift=False, ctrl=False)
            self.box["snap"] = data = self.line_snap_data(self._line_cage())
            hit = self.line_snap(pos, data)
            if hit is not None:                            # cizgi yakalanan hedeften baslar
                self.box.update(p0=hit[0], w0=hit[1])
            return
        if self.armed == "dim":                            # Distance: her iki modda, ilk / ikinci nokta ya da silme
            cb.setHandled()
            _later(lambda: self.dim_click(pos))
            return
        if self.armed in ("knife", "loop", "addpoint"):
            cb.setHandled()
            fn = {"knife": self.knife, "loop": self.loop_select, "addpoint": self.add_point_at}[self.armed]
            _later(lambda: fn(pos))
            return
        if self.armed == "magnet":
            cb.setHandled()
            _later(lambda: self.magnet_start(pos))
            return
        key = self.gizmo.pick(pos)
        if key is not None:
            cb.setHandled()                                # secime ve gezinmeye gitmesin
            self.drag_begin(key, pos)
            return
        if self.rail is not None:                          # kilavuz duzenleniyor: tiklama nokta / parca secer
            cb.setHandled()
            _later(lambda: self.rail_click(pos))
            return
        if self.mode != "object":                          # kilavuza cift tiklama: yeniden duzenle
            now, last = time.monotonic(), self._last_click
            self._last_click = (now, pos)
            if last is not None and now - last[0] < 0.4 and abs(pos[0] - last[1][0]) + abs(pos[1] - last[1][1]) < 8:
                ri = self.rail_pick(pos)
                if ri is not None:
                    cb.setHandled()
                    self._last_click = None
                    _later(lambda: self.rail_edit(ri))
                    return
        if (self.mode == "object" and self.object_pick(pos) is None
                and self.view.getObjectInfo((int(pos[0]), int(pos[1])))):
            return     # SubD olmayan sekil (Part / PartDesign): FreeCAD kendisi secsin (yuzey -> From face)
        # bos alana basildi: birakinca tiklama, surukleyince kutu secimi (Shift ekler, Ctrl cikarir)
        cb.setHandled()
        mods = QtWidgets.QApplication.keyboardModifiers()
        self.box = dict(p0=pos, moved=False,
                        shift=bool(ev.wasShiftDown() or mods & QtCore.Qt.ShiftModifier),
                        ctrl=bool(ev.wasCtrlDown() or _ctrl_held()))

    def _on_motion(self, cb):
        pos = cb.getEvent().getPosition().getValue()
        self.mouse = pos
        if self.magnet_size is not None:
            cb.setHandled()
            if QtWidgets.QApplication.mouseButtons() & QtCore.Qt.RightButton:
                self._schedule(self.magnet_size_move)
            else:                                          # sag tus gorunum disinda birakildi
                _later(self.magnet_size_end)
        elif self.drag is not None:
            if not _left_held():
                _later(self.drag_end)                      # sol tus gorunum disinda birakildi
                return
            cb.setHandled()
            self._schedule(self.drag_apply)
        elif self.grab:
            self._schedule(self._grab_step)
        elif self.box is not None:
            if not _left_held():                           # tus gorunum disinda birakildi: normal birakma gibi bitir
                self.finish_box(pos)
                return
            cb.setHandled()
            b = self.box
            p0 = b.get("press", b["p0"])                   # yakalanan baslangic imlecten uzak olabilir
            if b["moved"] or abs(pos[0] - p0[0]) + abs(pos[1] - p0[1]) > 4:
                b["moved"] = True
                if b.get("kind") == "line":
                    b["pos"] = pos                         # canli kesim: hareketler birlestirilir
                    self._schedule(lambda: self.line_live(b, b["pos"]))
                else:
                    _later(lambda: self._rubber_show(b["p0"], pos))
        else:
            key = self.gizmo.pick(pos)
            _later(lambda: self.gizmo.set_hover(key))
            if self.armed == "magnet":                     # imlec / isaret Magnet'in: _line_hover silmesin
                self._schedule(self._magnet_hover)
            elif self.armed in ("line", "mline", "array", "dim") or self._line_cursor is not None:
                self._schedule(self._line_hover)
            elif self.dim_hits:                            # x dugmesinin uzerine gelince kirmizi
                self._schedule(self._dim_hot_check)

    def finish_box(self, pos):
        """Sol tus birakildi: surukleme olduysa kutu secimi, olmadiysa tiklama secimi."""
        b, self.box = self.box, None
        if b is None:
            return
        if b.get("kind") == "line":
            if b["moved"]:
                _later(lambda: (self.line_live(b, pos), self.line_finish(b)))
            else:
                _later(lambda: self.line_finish(b))
            return
        if b["moved"]:
            _later(lambda: (self._rubber_hide(), self.box_select(b["p0"], pos, b["shift"], b["ctrl"])))
        else:
            _later(lambda: self.click(b["p0"], b["shift"]))

    def _rubber_show(self, p0, p1):
        """Secim kutusunu ciz: Coin olay konumu (sol alt kokenli, aygit pikseli) -> Qt (sol ust, mantiksal)."""
        vp = self.gv.viewport()
        if self.rubber is None:
            self.rubber = QtWidgets.QRubberBand(QtWidgets.QRubberBand.Rectangle, vp)
        dpr = vp.devicePixelRatioF()
        H = self.view.getSize()[1]
        a = QtCore.QPoint(int(p0[0] / dpr), int((H - p0[1]) / dpr))
        b = QtCore.QPoint(int(p1[0] / dpr), int((H - p1[1]) / dpr))
        self.rubber.setGeometry(QtCore.QRect(a, b).normalized())
        self.rubber.show()

    def _rubber_hide(self):
        if self.rubber is not None:
            self.rubber.hide()

    def _grab_step(self):
        self._pending = False
        self.grab_move()

    def key(self, k, ctrl=False):
        if self.drag is not None:
            return
        if self.box is not None and self.box.get("kind") == "line":
            if k == "escape":                              # knife cizgisi cekilirken: canli kesimi iptal
                self.line_cancel()
            return
        if self.grab:
            if k == "escape":
                self.end_grab(False)
            elif k in ("return", "enter"):
                self.end_grab(True)
            elif k in ("b", "i") and self.grab.get("kind") == "bevel":
                self.bevel_switch()
            elif k in ("x", "y", "z") and not ctrl and self.grab.get("kind") == "extrude":
                self.extrude_axis(k.upper())               # yalniz bu eksen; tekrar = normal yon
                self.refresh()
            return
        if k == "escape" and self.armed == "dim":          # Distance: ilk noktadan vazgec / araci kapat
            self.dim_escape()
            return
        if k == "q":                                       # Snap ac / kapat: her modda, kilavuz duzenleme de surer
            self.toggle_snap()
            return
        if self.mode == "object":                          # Object modu: yalniz nesne komutlari
            if ctrl and k in ("z", "redo"):
                doc = App.ActiveDocument
                if doc.UndoCount if k == "z" else doc.RedoCount:
                    doc.undo() if k == "z" else doc.redo()
                    doc.recompute()                       # katilar geri alinan kafese gore yenilensin
                self.doc_undone()
            elif ctrl and k == "v":                       # panodaki cokgenler -> yeni SubD nesnesi
                self.paste_faces()
            elif ctrl and k in ("c", "x", "w", "l"):      # nesne modunda kopyalanacak / baglanacak nokta yok
                pass
            elif ctrl:
                self.set_mode({"g": "point", "h": "face"}.get(k, "face"))
            elif k in ("tab", "o"):                   # O (ya da Tab): secili nesnenin kafesine gir
                self.edit_object()
            elif k in ("1", "2", "3"):
                self.set_mode({"1": "point", "2": "edge", "3": "face"}[k])
            elif k == "d":
                self.object_duplicate()
            elif k in ("delete", "backspace"):
                self.object_delete()
            elif k == "a":
                objs = self.subd_objects()
                self.set_object_selection([] if len(self.selected_objects()) == len(objs) else objs)
                self.refresh()
            elif k == "escape":
                Gui.Selection.clearSelection()
                self.refresh()
            elif k == "t":
                self.set_gizmo("move")
            elif k == "r":
                self.toggle_gizmo("rotate")
            elif k == "h":
                self.toggle_gizmo("scale")
            elif k == "n":
                self.toggle_options()
            return
        if k == "o":                                       # Edit -> Object modu (Tab LightWave gibi SubPatch)
            self.enter_object_mode()
            return
        if self.rail is not None and not ctrl:             # kilavuz duzenlenirken
            if k == "i":
                self.rail_add()
                return
            if k in ("delete", "backspace"):
                self.rail_delete()
                return
            if k in ("return", "enter", "escape"):
                self.rail_finish()
                return
            if k not in ("t", "x", "s", "tab", "n", "y"):   # baska bir arac / mod: kilavuz duzenleme biter
                self.rail_finish()
        if ctrl:                                           # LightWave: ^G nokta, ^H cokgen
            {"z": self.do_undo, "redo": self.do_redo, "g": lambda: self.set_mode("point"), "h": lambda: self.set_mode("face"),
             "c": self.copy_faces, "x": lambda: self.copy_faces(cut=True), "v": self.paste_faces, "w": self.weld,
             "l": self.connect}[k]()
        elif k == "n":
            self.toggle_options()
        elif k == "escape" and self.armed:
            self.drop_tool()
            self.refresh()
        elif k in ("tab", "s"):
            self.toggle_smooth()
        elif k == "e":
            self.extrude()
        elif k == "m":
            self.merge_points()
        elif k == "Z":
            self.merge_polygons()
        elif k == "k":
            self.toggle_tool("knife")
        elif k == "t":                                     # LightWave Move: gizmo (menudeki Move dugmesi gibi)
            self.set_gizmo("move")
        elif k == "g":                                     # Blender G: fareyle serbest tasima (gizmosuz)
            self.start_grab(self.selected_verts())
        elif k == "a":
            self.select_all()
        elif k in ("delete", "backspace"):
            self.delete()
        elif k == "f":
            self.save()
        elif k == "F":                                     # Shift+F: LightWave Flip
            self.flip_faces()
        elif k == "D":                                     # Shift+D: LightWave Subdivide
            self.subdivide()
        elif k in ("1", "2", "3"):
            self.set_mode({"1": "point", "2": "edge", "3": "face"}[k])
        elif k == "y":
            self.toggle_sym()
        elif k == "r":
            self.toggle_gizmo("rotate")
        elif k == "h":
            self.toggle_gizmo("scale")
        elif k == "c":
            self.crease_selected()
        elif k == "l":
            self.loop_select()
        elif k == "j":
            self.toggle_tool("line")
        elif k == "V":                                     # Shift+V: Mirror (line), LightWave Mirror gibi
            self.toggle_tool("mline")
        elif k == "p":                                     # LightWave Make Polygon
            self.make_polygon()
        elif k == "x":
            self.toggle_ghost()
        elif k == "b":
            self.edge_bevel() if self.mode == "edge" else self.bevel()
        elif k == "i":
            self.bevel(inset_only=True)
        elif k == "escape":
            self.clear_selection()
            self.refresh()

    # --- Calisma tezgahi degisince askiya al ---
    # Modelleyici yalniz SubD tezgahinda (ve baslatildigi tezgahta) gorunur. Baska tezgaha gecince (sketch duzenlenirken
    # FreeCAD kendisi Sketcher'a gecer) arac paneli, cubuklar ve fare / tus / Undo yakalama kalkar ki o tezgahin
    # araclari bozulmasin; oturum ve kafes kalir (sketch baglari ve otomatik kayit surer). Geri donunce her sey yerine gelir.
    def _wb_tick(self):
        """Modelleyici yalniz kendi belgesinde ve SubD (ya da baslatildigi) tezgahinda etkin; baska belgeye ya da
        tezgaha gecilince askiya alinir, geri donunce kaldigi yerden surer (ikinci dosyada ilkinin menuleri /
        tus ve fare yakalamasi calismasin)."""
        if getattr(self, "closed", False) or _sessions.get(self.doc_name) is not self:
            return                                         # kapanmis oturum geri canlanmasin (bekleyen geri cagri)
        now = _active_wb()
        gd = Gui.ActiveDocument
        if now is None or gd is None:
            return                                         # gecis sirasinda gecici: karar verme
        if now in self.wb_home and gd.Document.Name == self.doc_name:
            if self.suspended:
                self.resume()
        elif not self.suspended:
            self.suspend()

    def focus_ok(self):
        """Odak / pencere etkinlestirme yalniz modelleyici etkinken ve kendi belgesi onundeyken serbest; yoksa kullanicinin
        actigi diger belgenin sekmesi geri alinirdi."""
        gd = Gui.ActiveDocument
        return (not self.suspended) and gd is not None and gd.Document.Name == self.doc_name

    def suspend(self):
        global _session
        if self.suspended or getattr(self, "closed", False):
            return
        self.suspended = True
        if _session is self:
            _session = None
        if self.grab:                                      # yarim kalan arac / tasima: vazgec
            self.end_grab(False)
        if self.drag is not None:
            _quiet(self.drag_end)
        self.box = None
        self.drop_tool()
        _quiet(self._rubber_hide)
        _quiet(self._line_hide)
        self.dock_was_open = self.dock.isVisible()
        _quiet(self.dock.hide)
        _quiet(lambda: self.view.removeEventCallbackPivy(coin.SoMouseButtonEvent.getClassTypeId(), self._cb_btn))
        _quiet(lambda: self.view.removeEventCallbackPivy(coin.SoLocation2Event.getClassTypeId(), self._cb_mov))
        for w in self.key_targets:
            _quiet(lambda w=w: w.removeEventFilter(self.key_filter))
        app = QtWidgets.QApplication.instance()
        _quiet(lambda: app.removeEventFilter(self.wheel_guard))
        _quiet(lambda: app.removeEventFilter(self.history_guard))
        _quiet(self.timer.stop)
        _quiet(self.dim_timer.stop)
        if self.dim_root is not None:
            _quiet(lambda: self.parent.removeChild(self.dim_root))   # olculer de baska tezgahta gorunmesin
        self.dim_list_was_open = self.dim_list is not None and self.dim_list.isVisible()
        if self.dim_list is not None:
            _quiet(self.dim_list.hide)
        _quiet(lambda: self.gizmo.place(None))
        if self.overlay is not None:
            _quiet(lambda: self.overlay.set_visible(False))
        _quiet(lambda: self.parent.removeChild(self.root))   # kafes / onizleme sahneden de kalksin

    def resume(self):
        global _session
        if not self.suspended or getattr(self, "closed", False):
            return
        self.suspended = False
        _session = self
        if self.parent.findChild(self.root) < 0:
            _quiet(lambda: self.parent.addChild(self.root))
        _quiet(lambda: self.view.addEventCallbackPivy(coin.SoMouseButtonEvent.getClassTypeId(), self._cb_btn))
        _quiet(lambda: self.view.addEventCallbackPivy(coin.SoLocation2Event.getClassTypeId(), self._cb_mov))
        for w in self.key_targets:
            _quiet(lambda w=w: w.installEventFilter(self.key_filter))
        app = QtWidgets.QApplication.instance()
        _quiet(lambda: app.installEventFilter(self.wheel_guard))
        _quiet(lambda: app.installEventFilter(self.history_guard))
        _quiet(lambda: self.timer.start(50))
        _quiet(lambda: self.dim_timer.start(150))
        self.dim_sig = None
        if self.dim_root is not None:
            _quiet(self._dim_scene)                          # sahneye geri tak
        if self.dim_list_was_open and self.dim_list is not None:
            _quiet(self.dim_list.show)
        _quiet(self._dim_tick)
        if self.overlay is not None:
            _quiet(lambda: self.overlay.set_visible(True))
        if self.dock_was_open:
            _quiet(self.dock.show)
        self.refresh()
        if self.overlay is not None:
            self.overlay.refocus()

    def close(self, doc_gone=False):
        """Kafesi kaydet ve modelleyiciyi kaldir. Her adim ayri korunur: belge / 3B gorunum once kapandiysa
        widget'lari ve sahnesi zaten silinmistir; bir adimin hatasi kalanlari (or. arac cubugu) engellemesin,
        oturum her durumda kapansin (yoksa modelleyici bir daha acilamiyordu). doc_gone: belge kapaniyor,
        kaydetme / gosterme yok."""
        global _session
        if _session is self:
            _session = None
        if _sessions.get(self.doc_name) is self:
            del _sessions[self.doc_name]
        self.closed = True
        _quiet(lambda: _zebra_button(False))
        _quiet(self.timer.stop)
        _quiet(self.autosave.stop)
        _quiet(lambda: self.wb_timer.stop())
        if getattr(self, "_mdi", None) is not None:
            _quiet(lambda: self._mdi.subWindowActivated.disconnect(self._mdi_cb))
        _quiet(self.dim_timer.stop)
        if not doc_gone and self.mode != "object":         # Object modunda yazilacak kafes yok
            try:
                if len(self.cage.V):
                    self.save(commit=True)
                else:
                    self._drop_empty()
                self._show_edited()
            except Exception as e:
                App.Console.PrintError(T("SubD could not be saved: %s\n") % e)
        if not doc_gone:                                   # askidayken kaldirilmis olsa da zararsiz: kapanan oturumdan iz kalmasin
            _quiet(lambda: self.view.removeEventCallbackPivy(coin.SoMouseButtonEvent.getClassTypeId(), self._cb_btn))
            _quiet(lambda: self.view.removeEventCallbackPivy(coin.SoLocation2Event.getClassTypeId(), self._cb_mov))
        if True:
            # belge kapanirken 3B gorunum (gv + viewport) silinirken gelir: ona dokunmak cokertir, filtre zaten gider
            for w in (self.key_targets[2:] if doc_gone else self.key_targets):
                _quiet(lambda w=w: w.removeEventFilter(self.key_filter))
            _quiet(lambda: QtWidgets.QApplication.instance().removeEventFilter(self.wheel_guard))
            _quiet(lambda: QtWidgets.QApplication.instance().removeEventFilter(self.history_guard))
        if not doc_gone:                                   # belge kapanirken sahne zaten yok edilir; pivy'ye dokunmak cokertir
            _quiet(lambda: self.parent.removeChild(self.root))
            if self.dim_root is not None:
                _quiet(lambda: self.parent.removeChild(self.dim_root))
        if self.overlay is not None:
            _quiet(lambda: self.overlay.close(doc_gone))       # arac paneli ana pencerede: her durumda kaldirilir
            self.overlay = None
        if not doc_gone:                                   # (viewport() / silinmis widget'a dokunmak cokertir)
            _quiet(lambda: self._line_mark(None))
            for w in (self.rubber, self.line_overlay):
                if w is not None:
                    _quiet(w.deleteLater)
        self.rubber = self.line_overlay = None
        _quiet(self.dock.close)
        _quiet(self.dock.deleteLater)
        if self.dim_list is not None:
            _quiet(self.dim_list.close)
            _quiet(self.dim_list.deleteLater)


_MARK_ROLE = int(QtCore.Qt.UserRole) + 1077   # isaretli agac satirinin gorunmez etiketi


def _tree_mark(obj):
    """Model agacinda duzenlenen nesneyi isaretle (kalin, sari zemin); None: isareti kaldir. FreeCAD agaci
    QTreeWidget'tir; satir nesnenin etiketiyle, belgenin altinda aranir (yalniz yazi tipi / zemin degisir).
    Satir referansi tutulmaz: agac nesne eklenince satirlarini yeniden kurar (eski referans silinmis olur ve
    isaret kalkmazdi); isaretli satirlar gorunmez etiketle bulunur."""
    def walk(item, under_doc):
        if item.data(0, _MARK_ROLE):                       # eski isareti kaldir
            f = item.font(0)
            f.setBold(False)
            item.setFont(0, f)
            item.setBackground(0, QtGui.QBrush())
            item.setData(0, _MARK_ROLE, None)
        if obj is not None:
            under_doc = under_doc or item.text(0) == obj.Document.Label
            if under_doc and item.text(0) == obj.Label:
                f = item.font(0)
                f.setBold(True)
                item.setFont(0, f)
                item.setBackground(0, QtGui.QBrush(QtGui.QColor(232, 185, 35, 120)))
                item.setData(0, _MARK_ROLE, True)
        for i in range(item.childCount()):
            walk(item.child(i), under_doc)
    if obj is not None and not _alive(obj):
        obj = None
    try:
        for tw in Gui.getMainWindow().findChildren(QtWidgets.QTreeWidget):
            for i in range(tw.topLevelItemCount()):
                walk(tw.topLevelItem(i), False)
    except Exception:
        pass


def redges_of(cage, region):
    """Bolge yuzlerinin kenarlari."""
    return {_edge(f[i], f[(i + 1) % len(f)]) for fi in region for f in [cage.F[fi]] for i in range(len(f))}


def _seg_tri_hits(V, E, T):
    """E kenarlarindan T ucgenlerini (ortak noktasi olmayanlari) kesenler: [((a, b), (t0, t1, t2)), ...]."""
    out = []
    if not len(E) or not len(T):
        return out
    A, B = V[E[:, 0]], V[E[:, 1]]
    P = V[T]
    tmin, tmax = P.min(1), P.max(1)
    smin, smax = np.minimum(A, B), np.maximum(A, B)
    e1, e2 = P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]
    for k in range(len(E)):
        cand = np.nonzero(((tmin <= smax[k]) & (tmax >= smin[k])).all(1))[0]
        if not len(cand):
            continue
        a, b = E[k]
        cand = cand[~((T[cand] == a) | (T[cand] == b)).any(1)]
        if not len(cand):
            continue
        d = B[k] - A[k]
        q = np.cross(d, e2[cand])
        det = (e1[cand] * q).sum(1)
        ok = np.abs(det) > 1e-12
        inv = 1.0 / np.where(ok, det, 1.0)
        sv = A[k] - P[cand, 0]
        u = (sv * q).sum(1) * inv
        r = np.cross(sv, e1[cand])
        v = (r * d).sum(1) * inv
        t = (e2[cand] * r).sum(1) * inv
        h = ok & (u > 1e-6) & (v > 1e-6) & (u + v < 1 - 1e-6) & (t > 1e-6) & (t < 1 - 1e-6)
        for j in cand[h]:
            out.append(((int(a), int(b)), tuple(int(x) for x in T[j])))
    return out


def _alive(o):
    """Belge nesnesi hala belgede mi (silinmis / olusturmasi geri alinmis nesneye dokunmak hata verir)."""
    try:
        return o is not None and o.Document is not None and o.Document.getObject(o.Name) is not None
    except Exception:
        return False


def _quiet(fn):
    """Kapatirken: silinmis Qt / Coin nesnesine dokunma hatalarini yut (RuntimeError: already deleted...)."""
    try:
        fn()
    except Exception:
        pass


# ---------- Tutamak Attachment'i: konum korunur ----------
# FreeCAD'in Attachment penceresi referansi secince ofseti 0 yapar: tutamak (ve grubun noktalari) referansin
# ustune atlar. Tutamak bagli degilken bir referansa baglanirsa onceki yer korunur: ofset, o yeri verecek sekilde
# yazilir (ofset = bagli_yer^-1 * onceki_yer). Pencere acikken onizleme icin hemen, kapaninca (OK ofseti pencerenin
# 0 degerleriyle ezer) bir kez daha uygulanir. Kullanici kendi ofsetini yazdiysa (0 degil) dokunulmaz.

_ATTACH_PROPS = ("MapMode", "AttachmentSupport", "Support")
_attach_before = {}      # (belge, tutamak) -> (once bagli miydi, onceki Placement)
_attach_pending = {}     # (belge, tutamak) -> (korunacak Placement, bitis zamani)


def _is_handle(obj):
    return "SubDGroup" in getattr(obj, "PropertiesList", []) and "MapMode" in obj.PropertiesList


def _is_attached(obj):
    sup = obj.AttachmentSupport if "AttachmentSupport" in obj.PropertiesList else getattr(obj, "Support", [])
    return str(obj.MapMode) != "Deactivated" and len(sup) > 0


def handle_attach_before(obj):
    if not _is_handle(obj):
        return
    key = (obj.Document.Name, obj.Name)
    if key not in _attach_before:
        _attach_before[key] = (_is_attached(obj), App.Placement(obj.Placement))


def handle_attach_after(obj):
    if not _is_handle(obj):
        return
    key = (obj.Document.Name, obj.Name)
    was = _attach_before.get(key)
    if was is None:
        return
    now = _is_attached(obj)
    if now and not was[0]:
        _attach_pending[key] = (was[1], time.time() + 600)
        _attach_poll_start()
    if now or was[0]:                                       # karar verildi (bagli ya da bagi cozuldu)
        _attach_before.pop(key, None)


def handle_ref(h):
    """Bagli tutamagin referansi (baska bir grup tutamagi, LCS, duzlem, nesne yuzu ...); bagli degilse None."""
    if not _is_handle(h) or not _is_attached(h):
        return None
    sup = h.AttachmentSupport if "AttachmentSupport" in h.PropertiesList else getattr(h, "Support", [])
    ref = sup[0][0] if sup else None
    return ref


def handle_dim_setup(h):
    """Bir grup tutamagi bir referansa (baska tutamak, LCS ...) baglandiysa aralarindaki mesafeyi FreeCAD ifadesine baglar:
    tutamakta Distance ozelligi (Data sekmesi) olusur, AttachmentOffset.Base.z = Distance. Distance degisince
    tutamak ve grubun noktalari kayar. Zaten kuruluysa dokunmaz. Dondurur: kuruldu mu."""
    if handle_ref(h) is None or "Distance" in h.PropertiesList:
        return False
    h.addProperty("App::PropertyDistance", "Distance", "SubD", "Distance to the group it is attached to, along its normal")
    h.Distance = h.AttachmentOffset.Base.z
    h.setExpression("AttachmentOffset.Base.z", "Distance")
    h.recompute()
    return True


def handle_dim_geom(h):
    """Iki grup tutamagi arasindaki olcunun dunya noktalari (P: referans duzlemindeki iz, Q: tutamak) ve mesafe
    (referans duzleminin normali boyunca = ofsetin z'si); referans bir grup tutamagi degilse None."""
    if handle_ref(h) is None:
        return None
    off = h.AttachmentOffset
    frame = h.Placement.multiply(off.inverse())             # ofsetsiz baglanti duzlemi
    axis = frame.Rotation.multVec(App.Vector(0, 0, 1))
    q = h.Placement.Base
    return q - axis * off.Base.z, q, off.Base.z


def rigid_fit(P, Q):
    """P noktalarini Q'ya en iyi oturtan donusum (Kabsch): Q ~ P @ R.T + t. 3'ten az ya da bir dogru uzerindeki
    noktalarda yalniz oteleme (donme belirsiz)."""
    P, Q = np.asarray(P, float).reshape(-1, 3), np.asarray(Q, float).reshape(-1, 3)
    cp, cq = P.mean(0), Q.mean(0)
    R = np.eye(3)
    if len(P) >= 3:
        U, S, Vt = np.linalg.svd((P - cp).T @ (Q - cq))
        if S[1] > 1e-9 * max(S[0], 1e-12):
            d = np.sign(np.linalg.det(Vt.T @ U.T)) or 1.0
            R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, cq - R @ cp


def set_handle_placement(h, pl):
    """Tutamagin kuresel yerini pl yapar (bagliysa ofseti ayarlanir; Distance ifadesi varsa o da). Dondurur: basarildi mi."""
    try:
        if _is_attached(h):
            h.recompute()
            frame = h.Placement.multiply(h.AttachmentOffset.inverse())
            off = frame.inverse().multiply(pl)
            if "Distance" in h.PropertiesList and h.ExpressionEngine:
                h.Distance = off.Base.z
            h.AttachmentOffset = off
        else:
            h.Placement = pl
        h.recompute()
        return h.Placement.isSame(pl, 1e-6)
    except Exception:
        return False


def keep_handle_place(h, p_old):
    """Bagli tutamagin ofsetini, Placement'i p_old olacak sekilde yazar. Dondurur: degisti mi."""
    if not _is_attached(h) or not h.AttachmentOffset.isSame(App.Placement(), 1e-9):
        return False
    if h.Placement.isSame(p_old, 1e-6):
        return False
    h.recompute()                                           # ofset 0: Placement = referansin kendisi
    A = App.Placement(h.Placement)
    if A.isSame(p_old, 1e-6):
        return False
    h.AttachmentOffset = A.inverse().multiply(p_old)
    h.recompute()
    return True


_attach_timer = []


def _attach_poll_start():
    if _attach_timer:
        return
    t = QtCore.QTimer()
    t.setInterval(300)
    t.timeout.connect(lambda: _quiet(_attach_poll))
    _attach_timer.append(t)
    t.start()


def _attach_poll():
    dialog_open = bool(Gui.Control.activeDialog())
    for key, (p_old, until) in list(_attach_pending.items()):
        doc = App.listDocuments().get(key[0])
        h = doc.getObject(key[1]) if doc is not None else None
        if h is None or not _is_attached(h) or time.time() > until:
            _attach_pending.pop(key, None)
            continue
        if dialog_open:
            keep_handle_place(h, p_old)                     # onizleme: tutamak atlamasin
            continue
        _attach_pending.pop(key, None)                      # pencere kapandi: son durum
        doc.openTransaction(T("SubD keep handle place"))
        if h.AttachmentOffset.isSame(App.Placement(), 1e-9) and not h.Placement.isSame(p_old, 1e-6):
            keep_handle_place(h, p_old)
        handle_dim_setup(h)                                 # grup tutamagina baglandiysa mesafe ifadesi
        doc.commitTransaction()
    if not _attach_pending:
        while _attach_timer:
            _attach_timer.pop().stop()


class _DocWatcher(object):
    """Duzenlenen belge kapaninca modelleyiciyi de kapat: 3B gorunum silinince alt cubuk onunla gidiyor ama ana
    penceredeki arac paneli ve oturum ortada kaliyordu."""

    def slotDeletedDocument(self, doc):
        s = _sessions.get(doc.Name)
        if s is not None and s.doc_name == doc.Name:
            s.close(doc_gone=True)

    def slotUndoDocument(self, doc):
        s = _sessions.get(doc.Name)
        if s is not None and s.doc_name == doc.Name:
            _later(lambda: _quiet(s.doc_undone))

    def slotRedoDocument(self, doc):
        self.slotUndoDocument(doc)

    def slotBeforeChangeObject(self, obj, prop):
        if prop in _ATTACH_PROPS:
            _quiet(lambda: handle_attach_before(obj))

    def slotChangedObject(self, obj, prop):
        if prop in _ATTACH_PROPS:
            _quiet(lambda: handle_attach_after(obj))
        if prop == "Placement":
            _quiet(lambda: _handle_moved(obj))
        # Sketch bagi: kaynagin sekli degisince bagli SubD'ler hemen izlesin. Sketcher'da Cancel sketch'i geri yukleyip
        # yalniz onu yeniden hesapliyor; SubD "Touched" kalip duzenlenmis haliyle duruyordu.
        if prop != "Shape" or not App.GuiUp:
            return
        try:
            subs = [o for o in obj.InList if "SketchLinks" in o.PropertiesList and o.SketchLinks]
        except Exception:
            return
        if subs:
            key = (obj.Document.Name, obj.Name)
            if key not in _link_pending:
                _link_pending.add(key)
                _later(lambda: _quiet(lambda: _links_follow(*key)))

    def slotStartSaveDocument(self, doc, filename):
        # Ctrl+S: kafesin son hali dosyaya girsin (yoksa en son otomatik kayit / F ile yazilan hal kaydedilirdi)
        s = _sessions.get(doc.Name)
        if (s is not None and s.doc_name == doc.Name and s.mode != "object" and len(s.cage.V)
                and s.drag is None and s.grab is None and (s.dirty or s.obj is None)):
            try:
                s.save(full=False)
            except Exception as e:
                App.Console.PrintError(T("SubD could not be saved: %s\n") % e)


_handle_pending = set()  # yeri degisen grup tutamaklarinin SubD'leri (belge, nesne): olay dongusunde bir kez islenir


def _handle_moved(h):
    """Grup tutamagi tasindi: SubD'nin grubu izlesin (SubD tutamaga bagli degil: bagimlilik dongusu yok)."""
    if "SubDGroup" not in getattr(h, "PropertiesList", []) or not h.SubDGroup:
        return
    key = (h.Document.Name, h.SubDGroup)
    if key not in _handle_pending:
        _handle_pending.add(key)
        _later(lambda: _quiet(lambda: _handles_follow(*key)))


def _handles_follow(doc_name, name):
    _handle_pending.discard((doc_name, name))
    doc = App.listDocuments().get(doc_name)
    obj = doc.getObject(name) if doc is not None else None
    if obj is None or "ParamGroups" not in obj.PropertiesList or not obj.ParamGroups:
        return
    s = _sessions.get(doc_name)
    if s is not None and s.mode != "object" and s.obj is obj:
        s.groups_refresh()
        return
    obj_groups_update(obj)
    obj.recompute()


_link_pending = set()    # sekli degisen bag kaynaklari (belge, nesne): olay dongusunde bir kez islenir


def _links_follow(doc_name, name):
    """Kaynak degisti: onu kullanan SubD'ler izler (modelleyicinin duzenledigi kafes oturumda, digerleri recompute)."""
    _link_pending.discard((doc_name, name))
    doc = App.listDocuments().get(doc_name)
    src = doc.getObject(name) if doc is not None else None
    if src is None:
        return
    s = _sessions.get(doc_name)
    for o in src.InList:
        if "SketchLinks" not in o.PropertiesList or not o.SketchLinks or not o.SketchLink:
            continue
        if s is not None and s.mode != "object" and s.obj is o:
            s.links_refresh()
        else:
            o.recompute()                                  # kaynak ayniysa execute bir sey yapmaz


App.addDocumentObserver(_DocWatcher())


class _PickOrder(object):
    """FreeCAD secim sirasi (From face > Point icin): getSelectionEx nesneye gore grupluyor, farkli nesnelerden
    karisik secilen noktalarin sirasi kayboluyordu."""

    def __init__(self):
        self.picks = []      # (belge, nesne, alt eleman), secildigi sirayla

    def addSelection(self, doc, obj, sub, pnt):
        self.picks.append((doc, obj, sub))

    def removeSelection(self, doc, obj, sub):
        if (doc, obj, sub) in self.picks:
            self.picks.remove((doc, obj, sub))

    def clearSelection(self, doc):
        self.picks = []


_pick_order = _PickOrder()
if App.GuiUp:
    _quiet(lambda: Gui.Selection.removeObserver(getattr(Gui, "_subd_pick_order", None)))   # modul yeniden yuklendi
    Gui.Selection.addObserver(_pick_order)
    Gui._subd_pick_order = _pick_order


_KEYS = {
    QtCore.Qt.Key_Tab: "tab", QtCore.Qt.Key_S: "s", QtCore.Qt.Key_E: "e", QtCore.Qt.Key_K: "k",
    QtCore.Qt.Key_T: "t", QtCore.Qt.Key_G: "g", QtCore.Qt.Key_A: "a", QtCore.Qt.Key_F: "f",
    QtCore.Qt.Key_Z: "z", QtCore.Qt.Key_Delete: "delete", QtCore.Qt.Key_Backspace: "backspace",
    QtCore.Qt.Key_Escape: "escape", QtCore.Qt.Key_Return: "return", QtCore.Qt.Key_Enter: "enter",
    QtCore.Qt.Key_1: "1", QtCore.Qt.Key_2: "2", QtCore.Qt.Key_3: "3", QtCore.Qt.Key_Y: "y", QtCore.Qt.Key_R: "r",
    QtCore.Qt.Key_C: "c", QtCore.Qt.Key_B: "b", QtCore.Qt.Key_I: "i",
    QtCore.Qt.Key_H: "h", QtCore.Qt.Key_N: "n", QtCore.Qt.Key_L: "l", QtCore.Qt.Key_X: "x",
    QtCore.Qt.Key_J: "j", QtCore.Qt.Key_D: "d", QtCore.Qt.Key_O: "o", QtCore.Qt.Key_V: "v", QtCore.Qt.Key_P: "p",
    QtCore.Qt.Key_M: "m", QtCore.Qt.Key_W: "w", QtCore.Qt.Key_Q: "q",
}
# Ctrl+Z geri al, Ctrl+G nokta modu, Ctrl+H cokgen modu (LightWave), Ctrl+C / X / V cokgen kopyala / kes / yapistir,
# Ctrl+W weld, Ctrl+Y / Ctrl+Shift+Z ileri al
_CTRL_KEYS = ("z", "redo", "g", "h", "c", "x", "v", "w", "l")
_CTRL_ONLY = ("z", "v", "w")   # tek basina bizim degil (FreeCAD'e gider)


class _KeyFilter(QtCore.QObject):
    def __init__(self, session):
        super().__init__()
        self.session = session

    def eventFilter(self, obj, ev):
        if getattr(self.session, "closed", False) or getattr(self.session, "suspended", False):
            return False                                   # kapanmis / askidaki oturum olay yakalamasin
        t = ev.type()
        g = self.session.grab
        if t == QtCore.QEvent.Wheel and g is not None and g.get("kind") in ("ebevel", "extrude"):
            # kenar bevel / extrude surerken tekerlek yakinlastirmaz, segment sayisini degistirir
            d = ev.angleDelta().y()
            if d:
                (self.session.ebevel_wheel if g["kind"] == "ebevel" else self.session.extrude_wheel)(1 if d > 0 else -1)
            return True
        if t == QtCore.QEvent.Wheel and g is not None and g.get("kind") == "magnet":
            if ev.angleDelta().y():                        # Magnet surerken tekerlek: yaricap
                self.session.magnet_wheel(1 if ev.angleDelta().y() > 0 else -1)
            return True
        if t == QtCore.QEvent.MouseButtonRelease and ev.button() == QtCore.Qt.LeftButton and g is not None \
                and g.get("kind") == "magnet":             # birakma sahneye ulasmayabilir (kutu secimi gibi)
            _later(self.session.magnet_end)
            return False
        # Magnet'in sag tusu Qt duzeyinde: FreeCAD menusunu gezinme kodu sag tus birakilinca Coin geri cagrisindan once
        # aciyor (setHandled yetmez); basis ve birakmayi gormezse menu de acilmaz.
        if t in (QtCore.QEvent.MouseButtonPress, QtCore.QEvent.MouseButtonDblClick) \
                and ev.button() == QtCore.Qt.RightButton and self.session.armed == "magnet" \
                and g is None and self.session.drag is None:
            s = self.session
            vp = s.gv.viewport()
            dpr = vp.devicePixelRatioF()
            p = vp.mapFrom(obj, ev.position().toPoint()) if obj is not vp else ev.position().toPoint()
            s.mouse = (int(round(p.x() * dpr)), int(round(s.view.getSize()[1] - p.y() * dpr)))
            if s.magnet_size is None:
                s.magnet_size_start(s.mouse)
            return True
        if t == QtCore.QEvent.MouseButtonRelease and ev.button() == QtCore.Qt.RightButton \
                and self.session.magnet_size is not None:
            _later(self.session.magnet_size_end)
            return True
        if t == QtCore.QEvent.ContextMenu and self.session.armed == "magnet":
            return True                                    # sag tus Magnet'in: FreeCAD menusu acilmasin
        if t == QtCore.QEvent.Resize and self.session.overlay is not None:
            _later(self.session.overlay.relayout)          # gorunum boyutlaninca menuler yerinde kalsin
            return False
        if (t == QtCore.QEvent.MouseButtonRelease and ev.button() == QtCore.Qt.LeftButton
                and self.session.box is not None):
            # basisi biz yakalayinca FreeCAD gezinmesi birakmayi sahneye iletmiyor: Qt'de yakala
            s = self.session
            vp = s.gv.viewport()
            dpr = vp.devicePixelRatioF()
            p = vp.mapFrom(obj, ev.position().toPoint()) if obj is not vp else ev.position().toPoint()
            s.finish_box((p.x() * dpr, s.view.getSize()[1] - p.y() * dpr))
            return False
        if t not in (QtCore.QEvent.ShortcutOverride, QtCore.QEvent.KeyPress):
            return False
        cmd = _user_command(ev)                            # kullanicinin atadigi tus (varsayilan yok), yerlesiklerden once
        if cmd is not None:
            if t == QtCore.QEvent.ShortcutOverride:
                ev.accept()
                return True
            if not ev.isAutoRepeat():
                self.session.user_command(cmd)
            return True
        k = _KEYS.get(ev.key())
        ctrl = bool(ev.modifiers() & QtCore.Qt.ControlModifier)
        if k in ("f", "d") and ev.modifiers() & QtCore.Qt.ShiftModifier:
            k = k.upper()                                  # Shift+F: Flip (F: kaydet), Shift+D: Subdivide
        if k == "v" and not ctrl and ev.modifiers() & QtCore.Qt.ShiftModifier:
            k = "V"                                        # Shift+V: Mirror (line) (Ctrl+V FreeCAD'in)
        if ctrl and (k == "y" or (k == "z" and ev.modifiers() & QtCore.Qt.ShiftModifier)):
            k = "redo"
        extr = g is not None and g.get("kind") == "extrude"   # extrude surerken Z = eksen (Ctrl+Z degil)
        if k == "z" and not ctrl and not extr and ev.modifiers() & QtCore.Qt.ShiftModifier:
            k = "Z"                                        # Shift+Z: Merge polygons
        if k is None or (ctrl and k not in _CTRL_KEYS) or (not ctrl and k in _CTRL_ONLY and not (extr and k == "z")):
            return False
        if t == QtCore.QEvent.ShortcutOverride:
            ev.accept()  # tusu FreeCAD kisayollarindan once biz alalim
            return True
        if not ev.isAutoRepeat():
            self.session.key(k, ctrl)
        return True


def _std_history(action):
    """FreeCAD'in Undo / Redo eylemi mi: "z" / "redo" / None. Arac cubugundaki eylemin adi bos, komut adi ipucunda."""
    if action is None:
        return None
    for name, k in (("Std_Undo", "z"), ("Std_Redo", "redo")):
        if action.objectName() == name or name in action.toolTip():
            return k
    return None


class _HistoryGuard(QtCore.QObject):
    """Edit modunda FreeCAD'in Undo / Redo'su (Edit menusu, arac cubugu, odak gorunumde degilken Ctrl+Z / Ctrl+Y)
    belgeyi degil kafesi bir adim geri alir: oturum belgeye tek adim olarak yaziliyor, belge geri alinsa butun
    duzenleme ve adimlar giderdi. Araç cubugu dugmesinin oku (liste) belgenin kendi adimlarini gostermeye devam eder."""
    INPUT = (QtWidgets.QLineEdit, QtWidgets.QTextEdit, QtWidgets.QPlainTextEdit, QtWidgets.QAbstractSpinBox)

    def __init__(self, session):
        super().__init__()
        self.session = session

    def _active(self):
        s = self.session
        doc = App.ActiveDocument
        return (s.mode != "object" and doc is not None and doc.Name == s.doc_name
                and QtWidgets.QApplication.activeModalWidget() is None)

    def _key(self, obj, ev):
        mods = ev.modifiers()
        if not mods & QtCore.Qt.ControlModifier or mods & QtCore.Qt.AltModifier:
            return None
        if not isinstance(obj, QtWidgets.QWidget) or isinstance(obj, self.INPUT) or obj in self.session.key_targets:
            return None
        if ev.key() == QtCore.Qt.Key_Y:
            return "redo"
        if ev.key() == QtCore.Qt.Key_Z:
            return "redo" if mods & QtCore.Qt.ShiftModifier else "z"
        return None

    def eventFilter(self, obj, ev):
        if getattr(self.session, "closed", False) or getattr(self.session, "suspended", False):
            return False
        t = ev.type()
        if t in (QtCore.QEvent.ShortcutOverride, QtCore.QEvent.KeyPress):
            k = self._key(obj, ev)
            if k is None or not self._active():
                return False
            if t == QtCore.QEvent.ShortcutOverride:
                ev.accept()
            elif not ev.isAutoRepeat():
                self.session.key(k, True)
            return True
        if t not in (QtCore.QEvent.MouseButtonPress, QtCore.QEvent.MouseButtonRelease):
            return False
        if ev.button() != QtCore.Qt.LeftButton:
            return False
        pos = ev.position().toPoint()
        if isinstance(obj, QtWidgets.QToolButton):
            k = None if pos.x() >= obj.width() - 16 else _std_history(obj.defaultAction())   # sagdaki ok: belge listesi
        elif isinstance(obj, QtWidgets.QMenu):
            k = _std_history(obj.actionAt(pos))
        else:
            return False
        if k is None or not self._active():
            return False
        if t == QtCore.QEvent.MouseButtonRelease:
            if isinstance(obj, QtWidgets.QMenu):
                obj.close()
            _later(lambda: _quiet(lambda: self.session.key(k, True)))
        return True


# ---------- Viewport'a gomulu LightWave tarzi menuler ----------
# Sol: baslikli araç sutunu (kisayol sagda). Alt: durum cubugu (secim, modlar, simetri, SubPatch, Options).
# Ikisi de 3B gorunumun viewport'unun cocuk widget'i: gorunumun ustunde durur, tiklamalar gorunume gitmez.
# Dugmeler odak almaz (NoFocus), boylece klavye kisayollari gorunumde calismaya devam eder.

OVERLAY_CSS = """
QFrame#lwPanel, QFrame#lwBar { background: rgba(52, 52, 52, 238); border: 1px solid #1e1e1e; }
QLabel#lwHead { color: #e6e6e6; background: #2a2a2a; padding: 2px; font-weight: bold; }
QPushButton { background: #5a5a5a; color: #ececec; border: 1px solid #2b2b2b; border-radius: 2px;
              padding: 0px; min-height: 19px; }
QPushButton:hover { background: #6b6b6b; }
QPushButton:pressed { background: #4a4a4a; }
QPushButton:checked { background: #3f6fb0; border-color: #24426b; }
QPushButton[tool="true"]:checked { background: #e8b923; border-color: #8a6a00; }
QPushButton:disabled { color: #8a8a8a; background: #4c4c4c; }
QLabel:disabled { color: #7a7a7a; }
QLabel { color: #e0e0e0; background: transparent; }
QLabel#lwInfo { color: #f59e0b; font-weight: bold; background: rgba(40, 40, 40, 200); border-radius: 3px;
                padding: 3px 8px; }
QToolTip { color: #f0f0f0; background: #2a2a2a; border: 1px solid #1e1e1e; padding: 2px; }
"""


MENU_CSS = """
QMenu { background: #3a3a3a; color: #ececec; border: 1px solid #1e1e1e; padding: 2px; }
QMenu::item { padding: 3px 18px 3px 10px; }
QMenu::item:selected { background: #3f6fb0; }
"""


# Arac sutunu ve alt cubuk dugmelerinin ipuclari; acilir menu ogeleri "Dugme>Oge". Kisayol dugmenin uzerinde yaziyor.
TOOL_TIPS = {
    "Save": "Write the cage to the document now (it is also saved automatically every 20 s)",
    "Solid (STEP)": "Make a parametric B-rep solid from the SubD surface, ready for STEP export and PartDesign",
    "Close": "Write the cage to the document and close the modeler",
    "Box": "Add a box; set its size and segments in the panel",
    "Cylinder": "Add a new cylinder object along Z",
    "Cylinder>Cylinder": "Cylinder with triangle-fan caps",
    "Cylinder>Topo cylinder": "Cylinder with all-quad caps (segments are a multiple of 4)",
    "From face": "Make a cage from selected Part / PartDesign elements",
    "From face>Face": "Turn the selected faces into one cage",
    "From face>Edge": "Each chain of selected edges becomes one polygon",
    "From face>Point": "The selected points become one polygon, in the order they were picked",
    "Make polygon": "Make a polygon from 3 or more selected points, a closed edge loop, or two edges",
    "Duplicate": "Duplicate the selected objects",
    "Select all": "Select everything, or deselect if all is already selected",
    "Deselect": "Clear the selection; also cancels the running tool",
    "Select loop": "Stays on: click to select the loop under the cursor",
    "Extrude": "Pull out the selected faces or open edges; mouse wheel = segments, X / Y / Z = axis",
    "Bevel": "Polygon mode: inset and shift the selected faces. Edge mode: round the selected edges",
    "Inset": "Inset the selected faces into their own plane",
    "Knife (loop)": "Stays on: click an edge to cut a loop around the model, then slide it",
    "Knife (line)": "Stays on: draw a line to cut the model live with the plane through it",
    "Mirror (line)": "Stays on: draw a line, the model is mirrored across the plane through it; the far side is "
                     "replaced (source: the selection's side, else left of a downward line; Flip side in the panel)",
    "Add points": "Split the selected edges into equal parts (Count in the panel); with nothing selected it stays on: "
                  "click an edge to add a point there",
    "Connect": "Add an edge between two points of one polygon, or between the middles of two of its edges",
    "Subdivide": "Split the selected polygons into quads (a polygon with n sides becomes n quads)",
    "Array": "Copy the selected polygons as separate parts: around an axis (Radial) or along X / Y / Z (Rectangular); "
             "switch in the panel",
    "Circle": "Arrange the selected points, an edge loop / chain or the boundary of a polygon region on a circle "
              "(inner points follow). Type the Radius in the panel (Auto radius = average of the current one)",
    "Twist": "Twist the selected polygons around an axis: the angle grows along the axis by the Falloff curve. Move the "
             "mouse for the angle (Ctrl: steps), click to accept, then fine-tune in the panel",
    "Distance": "Measure the distance between two points: cage points, edge middles / edges, polygon centers, other objects (corners / edges) or the surface. Keep: the measurement stays, follows the points and is deleted with the × next to its text; otherwise it is temporary. Works in Object mode too (temporary)",
    "Clear all": "Delete all distances (undo brings the fixed ones back)",
    "Distances": "List of the distances: give each a name and a target, show / hide, delete, export as CSV. The list stays with the object between sessions",
    "Draw axis": "Drag a line in the view to set the axis (snaps like Mirror (line)); a click gives an axis along the view",
    "Axis X": "Axis direction along world X (the axis point stays)",
    "Axis Y": "Axis direction along world Y (the axis point stays)",
    "Axis Z": "Axis direction along world Z (the axis point stays)",
    "Bridge": "Connect two open edge loops or two face regions",
    "Bridge>Bridge": "Connect two open edge loops / chains or two face regions with a strip",
    "Bridge>Rail bridge": "Pipe with an editable guide curve, between two regions or growing out of one",
    "Hole": "Cut a round hole in the middle of the selected faces and re-mesh around it",
    "Circle on face": "Like Hole, but the circle stays filled and selected, ready for Extrude or Bevel",
    "Merge polygons": "Merge the selected adjacent polygons into one",
    "Parametric group": "Make the selected points a group with a datum plane handle (SubD_PG_0001) in the FreeCAD tree: "
                        "move the datum or attach it (Attachment, offset, expressions) and the points follow",
    "Delete": "Delete the selected elements (object mode: the selected objects)",
    "Copy": "Copy the selected polygons with their creases; edges and points are copied as free edges / points",
    "Cut": "Copy the selected polygons, edges or points, then remove them",
    "Paste": "Paste polygons, edges and points in place into this cage",
    "Paste as object": "Paste the copied polygons, edges and points as a new SubD object",
    "Move": "Move gizmo: arrows along X / Y / Z, squares in a plane, cube in the screen plane",
    "Rotate": "Rotate gizmo",
    "Stretch / Size": "Stretch the selection along an axis or size it uniformly",
    "Magnet": "Stays on: drag on the cage to pull it softly inside the radius; a selection limits it",
    "Crease": "Make the selected edges sharp, or soft again; exact dimensions need sharp edges",
    "Reduce points": "Remove in-between points along nearly straight edges (selection or whole cage)",
    "Merge points": "Merge points at the same position (selection or whole cage)",
    "Weld": "Weld the selected points into the last one clicked",
    "Flip": "Flip the normals of the selected polygons (no selection: all)",
    "Shell": "Give the object or selected polygons a wall thickness",
    "Mirror": "Mirror the selected polygons (no selection: whole cage) across a world plane and weld; "
              "object mode: mirrored copies",
    "Mirror>Mirror X": "Mirror across the plane perpendicular to X",
    "Mirror>Mirror Y": "Mirror across the plane perpendicular to Y",
    "Mirror>Mirror Z": "Mirror across the plane perpendicular to Z",
    "Undo": "Undo one step (object mode: one whole edit)",
    "Redo": "Redo the step you undid (object mode: one whole edit)",
    "Object": "Object mode: move, copy and select whole SubD objects",
    "Edit": "Edit the cage of the selected SubD object",
    "Point": "Point selection mode",
    "Edge": "Edge selection mode",
    "Polygon": "Polygon selection mode",
    "Symmetry": "X symmetry: tools also act on the mirrored side",
    "Subpatch": "Show the smooth SubD surface or only the cage",
    "Ghost": "See through the model and select hidden elements",
    "Snap": "Free moves snap to the vertices and edges of visible Part / PartDesign shapes",
    "Options": "Show or hide the Options panel",
}


def _tool_tip(name):
    tip = TOOL_TIPS.get(name)
    return T(tip) if tip else ""


def _min_width(b, px):
    """Dugmenin en az genisligi (kenarlik dahil). Stille: FreeCAD 27 temasinin QPushButton { min-width: 80px }
    kurali setMinimumWidth / setFixedSize degerini eziyor, butun dugmeler 82 piksel oluyordu."""
    b.setStyleSheet("QPushButton { min-width: %dpx; }" % max(px - 2, 0))


def _lw_button(text, key="", width=None):
    """LightWave dugmesi: solda ad, sagda soluk kisayol."""
    b = QtWidgets.QPushButton()
    b.setFocusPolicy(QtCore.Qt.NoFocus)
    lay = QtWidgets.QHBoxLayout(b)
    lay.setContentsMargins(5, 0, 5, 0)
    lay.setSpacing(6)
    for t, css in ((text, ""), (key, "color:#a8a8a8;")):
        lb = QtWidgets.QLabel(t)
        lb.setStyleSheet(css)
        lb.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents)
        lay.addWidget(lb, 1 if t is text else 0)
    fm = b.fontMetrics()
    full = max(width or 0, fm.horizontalAdvance(text + key) + 24)                # yazi hicbir zaman kesilmesin
    _min_width(b, full)
    # alt cubuk sigmayinca kisayol yazisi gizlenir (_Overlay._fit_bar): (kisayol etiketi, tam / kisa genislik)
    b.lw_key = (b.findChildren(QtWidgets.QLabel)[1], full, max(width or 0, fm.horizontalAdvance(text) + 20)) \
        if key else None
    return b


class _LineOverlay(QtWidgets.QWidget):
    """Knife cizgisi: viewport'un ustunde, fare olaylarini gecirir."""

    def __init__(self, parent):
        super().__init__(parent)
        self.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents)
        self.a = self.b = self.mark = self.kind = self.ring = None
        self.ring_r = 0.0
        self.resize(parent.size())
        self.show()

    def set_line(self, a, b):
        self.a, self.b = a, b
        self._changed()

    def set_mark(self, pt, kind):
        """Yakalama hedefi: "point" kare, "mid" ucgen, "face" daire, "surf" elmas, "edge" capraz."""
        self.mark, self.kind = pt, kind
        self._changed()

    def set_ring(self, pt, r):
        """Magnet etki alani: merkez ve yaricap (mantiksal piksel)."""
        self.ring, self.ring_r = pt, r
        self._changed()

    def _changed(self):
        self.resize(self.parent().size())
        self.lower()                                       # panellerin ve alt cubugun altinda kalsin
        self.update()

    def paintEvent(self, ev):
        if self.a is None and self.mark is None and self.ring is None:
            return
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        if self.ring is not None:
            p.setPen(QtGui.QPen(QtGui.QColor(245, 158, 11), 1.5, QtCore.Qt.DashLine))
            p.setBrush(QtCore.Qt.NoBrush)
            p.drawEllipse(self.ring, self.ring_r, self.ring_r)
        if self.a is not None:
            pen = QtGui.QPen(QtGui.QColor(245, 158, 11), 2, QtCore.Qt.DashLine)
            p.setPen(pen)
            p.drawLine(self.a, self.b)
            p.setBrush(QtGui.QColor(245, 158, 11))
            for pt in (self.a, self.b):
                p.drawEllipse(pt, 3, 3)
        if self.mark is not None:
            p.setPen(QtGui.QPen(QtGui.QColor(34, 211, 238), 2))
            p.setBrush(QtCore.Qt.NoBrush)
            x, y, r = self.mark.x(), self.mark.y(), 6.0
            if self.kind == "point":
                p.drawRect(QtCore.QRectF(x - r, y - r, 2 * r, 2 * r))
            elif self.kind == "mid":
                p.drawPolygon(QtGui.QPolygonF([QtCore.QPointF(x, y - r), QtCore.QPointF(x + r, y + r),
                                               QtCore.QPointF(x - r, y + r)]))
            elif self.kind == "face":
                p.drawEllipse(QtCore.QPointF(x, y), r, r)
            elif self.kind == "surf":
                p.drawPolygon(QtGui.QPolygonF([QtCore.QPointF(x, y - r), QtCore.QPointF(x + r, y),
                                               QtCore.QPointF(x, y + r), QtCore.QPointF(x - r, y)]))
            else:
                p.drawLine(QtCore.QPointF(x - r, y - r), QtCore.QPointF(x + r, y + r))
                p.drawLine(QtCore.QPointF(x - r, y + r), QtCore.QPointF(x + r, y - r))
        p.end()


class _DimList(QtWidgets.QDockWidget):
    """Distances penceresi: olculerin adi, mesafesi, hedefi; goster / gizle, sil, CSV. Edit modunda duzenlenir
    (geri alinabilir), Object modunda butun SubD nesnelerinin olculeri salt okunur gosterilir."""
    HEAD = ("", "Name", "Distance", "dX", "dY", "dZ", "Target", "Tol.", "Dev.", "")

    def __init__(self, session):
        super().__init__(T("Distances"), Gui.getMainWindow())
        self.s = session
        self.setObjectName("SubDDistances")
        self.key = None
        self.map = []                                      # tablo satiri -> kafes sirasi
        self._busy = False
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        top = QtWidgets.QHBoxLayout()
        self.names = QtWidgets.QCheckBox(T("Show names"))
        self.names.setToolTip(T("Show the name in front of the distance in the 3D view"))
        self.names.setChecked(session.dim_names())
        self.names.toggled.connect(session.dim_set_names)
        top.addWidget(self.names)
        top.addStretch()
        for text, fn in ((T("Show all"), lambda: session.dim_show_all(True)),
                         (T("Hide all"), lambda: session.dim_show_all(False)),
                         (T("Export CSV..."), session.dim_export)):
            b = QtWidgets.QPushButton(text)
            b.setFocusPolicy(QtCore.Qt.NoFocus)
            b.clicked.connect(lambda checked=False, fn=fn: _later(fn))
            top.addWidget(b)
        lay.addLayout(top)
        self.table = QtWidgets.QTableWidget(0, len(self.HEAD))
        self.table.setHorizontalHeaderLabels([T(h) if h else "" for h in self.HEAD])
        self.table.verticalHeader().hide()
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.table.setAlternatingRowColors(True)
        hh = self.table.horizontalHeader()
        hh.setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
        self.table.itemChanged.connect(self._changed)
        self.table.itemSelectionChanged.connect(self._selected)
        lay.addWidget(self.table)
        self.hint = QtWidgets.QLabel()
        self.hint.setStyleSheet("color:#9a9a9a;")
        lay.addWidget(self.hint)
        self.setWidget(w)
        Gui.getMainWindow().addDockWidget(QtCore.Qt.BottomDockWidgetArea, self)
        self.setFloating(True)
        self.resize(720, 280)

    def _ph(self, idx):
        return T("Distance %d") % (idx + 1)

    def refresh(self, rows, force=False):
        if self.table.state() == QtWidgets.QAbstractItemView.EditingState:
            return                                         # yazarken liste yeniden kurulmasin
        edit = self.s.mode != "object"
        key = (edit, tuple((r["obj"], r["idx"], r["name"], r["hide"], r["target"], r["tol"]) for r in rows))
        if force or key != self.key:
            self.key = key
            self._build(rows, edit)
        self._values(rows)

    def _build(self, rows, edit):
        self._busy = True
        t = self.table
        t.clearContents()
        t.setRowCount(len(rows))
        self.map = [r["idx"] for r in rows]
        ro = QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable
        rw = ro | QtCore.Qt.ItemIsEditable
        for i, r in enumerate(rows):
            on = QtWidgets.QTableWidgetItem()
            on.setFlags(ro | (QtCore.Qt.ItemIsUserCheckable if edit else QtCore.Qt.NoItemFlags))
            on.setCheckState(QtCore.Qt.Unchecked if r["hide"] else QtCore.Qt.Checked)
            on.setToolTip(T("Show / hide this distance in the 3D view"))
            t.setItem(i, 0, on)
            hd = r.get("handle")
            name = QtWidgets.QTableWidgetItem(r["name"] or (self._ph(r["idx"]) if not hd else ""))
            if not r["name"] and not hd:
                name.setForeground(QtGui.QBrush(QtGui.QColor(150, 150, 150)))
            if not edit:
                name.setText("%s > %s" % (r["obj"], r["name"] or self._ph(r["idx"])))
            name.setFlags(ro if hd or not edit else rw)
            t.setItem(i, 1, name)
            if hd:
                on.setFlags(ro)
                name.setToolTip(T("Group distance: the Distance value can be edited here"))
            for c in (2, 3, 4, 5, 8):
                it = QtWidgets.QTableWidgetItem()
                it.setFlags(rw if hd and edit and c == 2 else ro)
                it.setTextAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
                t.setItem(i, c, it)
            for c, v in ((6, r["target"]), (7, r["tol"] if r["target"] is not None else None)):
                it = QtWidgets.QTableWidgetItem("" if v is None else "%.3f" % v)
                it.setFlags(rw if edit and not hd else ro)
                it.setTextAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
                t.setItem(i, c, it)
            if edit and not hd:
                b = QtWidgets.QPushButton("\u00d7")
                b.setFixedWidth(26)
                b.setToolTip(T("Delete this distance"))
                b.clicked.connect(lambda checked=False, k=r["idx"]: _later(lambda: self.s._dim_remove(k)))
                t.setCellWidget(i, 9, b)
        self.hint.setText(T("Name / Target / Tol. are editable (double click); the selected row is drawn blue in the view")
                          if edit else T("Object mode: all distances of all SubD objects (edit them in Edit mode)"))
        sel = self.s.dim_sel
        if sel in self.map:
            t.selectRow(self.map.index(sel))
        self._busy = False

    def _values(self, rows):
        self._busy = True
        t = self.table
        red, green, plain = QtGui.QColor(248, 113, 113), QtGui.QColor(110, 200, 90), QtGui.QBrush()
        for i, r in enumerate(rows[:t.rowCount()]):
            ok = r["dist"] is not None
            vals = {2: "%.3f" % r["dist"] if ok else "-"}
            for c, k in ((3, 0), (4, 1), (5, 2)):
                vals[c] = "%.3f" % r["d3"][k] if ok else "-"
            dev = r["dist"] - r["target"] if ok and r["target"] is not None else None
            vals[8] = "" if dev is None else "%+.3f" % dev
            for c, text in vals.items():
                it = t.item(i, c)
                if it is not None and it.text() != text:
                    it.setText(text)
            if dev is not None:
                col = green if abs(dev) <= r["tol"] else red
                for c in (2, 8):
                    t.item(i, c).setForeground(QtGui.QBrush(col))
            else:
                for c in (2, 8):
                    t.item(i, c).setForeground(plain)
        self._busy = False

    def _changed(self, it):
        if self._busy or not it:
            return
        row, col = it.row(), it.column()
        if not 0 <= row < len(self.map):
            return
        idx, s = self.map[row], self.s
        if isinstance(idx, str):                           # grup olcusu: yalniz Distance yazilabilir
            if col == 2:
                try:
                    num = float(it.text().strip().replace(",", "."))
                except ValueError:
                    _later(lambda: self.refresh(s._dim_rows(), force=True))
                    return
                _later(lambda: s.handle_dim_set(idx[2:], num))
            return
        if col == 0:
            val = ("hide", it.checkState() != QtCore.Qt.Checked)
        elif col == 1:
            txt = it.text().strip()
            val = ("name", "" if txt == self._ph(idx) else txt)
        elif col in (6, 7):
            txt = it.text().strip().replace(",", ".")
            try:
                num = float(txt) if txt else None
            except ValueError:
                _later(lambda: self.refresh(s._dim_rows(), force=True))   # gecersiz sayi: eski deger geri gelsin
                return
            if col == 6:
                val = ("target", num)
            else:
                val = ("tol", DIM_TOL if num is None else abs(num))
        else:
            return
        _later(lambda: s.dim_edit(idx, *val))

    def _selected(self):
        if self._busy:
            return
        rows = self.table.selectionModel().selectedRows()
        r = rows[0].row() if rows else None
        k = self.map[r] if r is not None and r < len(self.map) and self.s.mode != "object" else None
        self.s.dim_select(k if isinstance(k, int) else None)


class _WheelGuard(QtCore.QObject):
    """Gorunumdeki panelin sayi kutularinda fare tekerlegi (uygulama duzeyinde, FreeCAD'in filtrelerinden ve 3B
    gorunumden once): imlec kutunun ustundeyse her tik degeri bir adim degistirir ve olay tuketilir. Kutuya
    birakinca FreeCAD odaksiz kutunun tekerlegini yutuyor, olay gorunume gidip zoom yapiyordu (bazen ikisi
    birden). Odak gorunumde kalir (kisayollar calisir)."""

    def __init__(self, viewport):
        super().__init__()
        self.vp = viewport
        self.last = None

    def eventFilter(self, obj, ev):
        if ev.type() != QtCore.QEvent.Wheel:
            return False
        w = QtWidgets.QApplication.widgetAt(QtGui.QCursor.pos())
        while w is not None and not isinstance(w, QtWidgets.QAbstractSpinBox):
            w = w.parentWidget()
        if w is None or not self.vp.isAncestorOf(w):
            return False
        key = (ev.timestamp(), ev.angleDelta().y())
        if key == self.last:                               # ayni tik ikinci kez (ust widget'a iletilen): yalniz tuket
            return True
        self.last = key
        if w.isEnabled():
            acc = (w.property("wheel_acc") or 0) + ev.angleDelta().y()   # hassas tekerlek: 120 = bir adim
            steps = int(acc / 120)
            w.setProperty("wheel_acc", acc - steps * 120)
            if steps:
                w.stepBy(steps)
        return True


class _ClickFocus(QtCore.QObject):
    """Viewport icindeki kutuya tiklaninca klavye odagi kutuya gecmiyor (3B gorunum geri aliyor): yazilan
    rakamlar gorunume gidip 1/2/3 mod degistiriyordu. Tiklamada odagi acikca kutuya ver, degeri sec."""

    def eventFilter(self, obj, ev):
        if ev.type() == QtCore.QEvent.MouseButtonPress and isinstance(obj, QtWidgets.QLineEdit):
            box = obj.parent()
            def grab():
                box.setFocus(QtCore.Qt.MouseFocusReason)
                obj.selectAll()
            _later(grab)
        return False


class _ToolPanel(QtWidgets.QFrame):
    """Viewport icindeki kucuk islem paneli (RaytracerWB deger paneli gibi): yapilan isleme gore alanlar.
    Islem surerken degerler canli akar; bitince alanlar duzenlenebilir: yazip Enter = islem bu degerle
    yeniden uygulanir. Islem yokken secimin merkezi / boyutu."""

    def __init__(self, parent, session):
        super().__init__(parent)
        self.s = session
        self.setObjectName("lwPanel")
        self.setStyleSheet(OVERLAY_CSS + """
            QDoubleSpinBox { background: #2c2c2c; color: #facc15; border: 1px solid #1e1e1e; padding: 1px 3px; }
            QDoubleSpinBox:focus { border-color: #f59e0b; }
            QDoubleSpinBox:disabled { color: #f59e0b; background: #333; }
            QDoubleSpinBox::up-button, QDoubleSpinBox::down-button { width: 0px; height: 0px; border: none; margin: 0px; padding: 0px; }
            QDoubleSpinBox::up-arrow, QDoubleSpinBox::down-arrow { image: none; width: 0px; height: 0px; }
            QCheckBox { color: #e0e0e0; spacing: 3px; }""")
        # Oklar (yukari / asagi dugmeleri) setButtonSymbols(NoButtons) ile degil, yukaridaki stil kuraliyla gizlenir:
        # Qt 6.11'de (FreeCAD 27.x) stil dosyali kutuda NoButtons duzenleme alanini 1 piksele daraltiyor, yazi gorunmuyordu.
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(5, 4, 5, 6)
        lay.setSpacing(3)
        self.title = QtWidgets.QLabel()
        self.title.setObjectName("lwHead")
        self.title.setAlignment(QtCore.Qt.AlignCenter)
        lay.addWidget(self.title)
        # istege bagli kutucuklar: spec'in 4. ogesi [(ad, isaretli, fn(on))], 5. ogesi satir basligi. Baslik varsa
        # kutucuklar baslikla tek satirda (extrude: Eksen X Y Z), yoksa alt alta (Knife line secenekleri)
        self.check_bar = QtWidgets.QWidget()
        self.check_lay = QtWidgets.QGridLayout(self.check_bar)
        self.check_lay.setContentsMargins(0, 0, 0, 0)
        self.check_lay.setHorizontalSpacing(8)
        self.check_lay.setVerticalSpacing(2)
        self.checks, self.check_fns, self.check_key = [], [], None
        self.side_fns = []                                 # kutucuklarin sagindaki dugmeler (spec'in 7. ogesi)
        lay.addWidget(self.check_bar)
        self.check_bar.hide()
        self.grid = QtWidgets.QGridLayout()
        self.grid.setHorizontalSpacing(6)
        self.grid.setVerticalSpacing(2)
        lay.addLayout(self.grid)
        # istege bagli dugmeler: spec'in 6. ogesi [(ad, fn())], alanlarin altinda (Rail bridge: segment ekle...)
        self.btn_bar = QtWidgets.QWidget()
        self.btn_lay = QtWidgets.QGridLayout(self.btn_bar)
        self.btn_lay.setContentsMargins(0, 2, 0, 0)
        self.btn_lay.setSpacing(3)
        self.btn_key, self.btn_fns = (), []
        lay.addWidget(self.btn_bar)
        self.btn_bar.hide()
        self.key = None
        self.rows = []
        self._focus_filter = _ClickFocus(self)
        self.setFixedWidth(220)
        self.hide()

    def show_spec(self, spec):
        """Paneli spec'e gore gunceller; boyutu degistiyse True (yeniden yerlestir)."""
        if spec is None:
            changed = self.isVisible()
            self.hide()
            self.key = None
            return changed
        title, fields, editable = spec[:3]
        checks = spec[3] if len(spec) > 3 else []
        label = spec[4] if len(spec) > 4 else None
        buttons = spec[5] if len(spec) > 5 else []
        side = spec[6] if len(spec) > 6 else []
        key = (title, tuple(f[0] for f in fields), tuple(c[0] for c in checks), label, tuple(b[0] for b in buttons))
        changed = key != self.key
        self.btn_fns = [b[1] for b in buttons]
        if tuple(b[0] for b in buttons) != self.btn_key:
            self._build_buttons([b[0] for b in buttons], [b[2] if len(b) > 2 else None for b in buttons])
        self.check_fns = [c[2] for c in checks]
        self.side_fns = [b[1] for b in side]
        if (label, tuple(c[0] for c in checks), tuple(b[0] for b in side)) != self.check_key:
            self._build_checks(label, [c[0] for c in checks], [len(c) > 3 and c[3] for c in checks],
                               [b[0] for b in side])
        for cb, c in zip(self.checks, checks):
            if cb.isChecked() != bool(c[1]):
                cb.blockSignals(True)
                cb.setChecked(bool(c[1]))
                cb.blockSignals(False)
        if self.check_bar.isVisible() != bool(checks):
            self.check_bar.setVisible(bool(checks))
        if changed:
            for lb, box in self.rows:                      # eski satirlari yerlesimden hemen cikar
                for wdg in (lb, box):
                    self.grid.removeWidget(wdg)
                    wdg.hide()
                    wdg.setParent(None)
                    wdg.deleteLater()
            self.rows = []
            for r, (name, val, suffix, dec, fn) in enumerate(fields):
                lb = QtWidgets.QLabel(T(name))
                box = QtWidgets.QDoubleSpinBox()
                box.setRange(-1e7, 1e7)
                box.setDecimals(dec)
                box.setSuffix(suffix)
                box.setKeyboardTracking(False)
                # NoButtons burada kullanilmaz: Qt 6.11'de stil dosyasiyla birlikte metin alanini 1 piksele daraltir
                # (oklar __init__ icindeki stil kuraliyla sifir genislige indirildi)
                box.setAlignment(QtCore.Qt.AlignRight)
                # fare tekerlegi alanin ustunde: deger bir adim degisir ve islem hemen yeniden yapilir (segment,
                # kesim sayisi, mm); adim: tam sayi / mm / derece 1, carpan (4 ondalik) 0.05, oteki 0.1
                box.setSingleStep(1.0 if dec == 0 or suffix.strip() in ("mm", "°") else 0.05 if dec >= 4 else 0.1)
                box.editingFinished.connect(lambda b=box, i=r: self._commit(i, b))
                box.valueChanged.connect(lambda val, b=box, i=r: self._commit(i, b))
                box.lineEdit().installEventFilter(self._focus_filter)
                box.installEventFilter(self._focus_filter)          # tekerlek: _WheelGuard (zoom yok)
                self.grid.addWidget(lb, r, 0)
                self.grid.addWidget(box, r, 1)
                lb.show()                                  # gorunur panele eklenen cocuk kendiliginden gorunmez;
                box.show()                                 # gorunmezse yukseklik hesabina girmez
                self.rows.append((lb, box))
            self.key = key
            self.title.setText(title)
        self.fns = [f[4] for f in fields]
        self.vals = [float(f[1]) for f in fields]
        for (lb, box), f in zip(self.rows, fields):
            box.setEnabled(editable)
            if not box.hasFocus():
                box.blockSignals(True)
                box.setValue(float(f[1]))
                box.blockSignals(False)
        if not self.isVisible():
            changed = True
            self.show()
            self.raise_()
        if changed:
            self.layout().activate()
            self.setFixedHeight(self.layout().sizeHint().height())
        return changed

    def _build_buttons(self, names, icons=None):
        """Yazili dugmeler iki sutunda; simgeli olanlar (icons[i] QIcon) en altta yan yana kucuk kareler, ad
        ipucu olur."""
        while self.btn_lay.count():
            it = self.btn_lay.takeAt(0)
            w = it.widget()
            if w is None and it.layout() is not None:
                while it.layout().count():
                    w2 = it.layout().takeAt(0).widget()
                    if w2 is not None:
                        w2.hide()
                        w2.setParent(None)
                        w2.deleteLater()
            if w is not None:
                w.hide()
                w.setParent(None)
                w.deleteLater()
        icons = icons or [None] * len(names)
        k, row = 0, None
        for i, (n, ic) in enumerate(zip(names, icons)):
            if ic is not None:
                b = QtWidgets.QToolButton()
                b.setIcon(ic)
                b.setIconSize(QtCore.QSize(20, 20))
                b.setFixedSize(28, 28)
                b.setToolTip(T(n))
                if row is None:
                    row = QtWidgets.QHBoxLayout()
                    row.setSpacing(3)
            else:
                b = QtWidgets.QPushButton(T(n))
                b.setToolTip(_tool_tip(n))
                if sum(1 for x in icons if x is None) == 1:    # tek yazili dugme: tam genislik
                    self.btn_lay.addWidget(b, 0, 0, 1, 2)
                else:
                    self.btn_lay.addWidget(b, k // 2, k % 2)
                k += 1
            b.setFocusPolicy(QtCore.Qt.NoFocus)            # kisayollar gorunumde kalsin
            b.clicked.connect(lambda checked=False, i=i: self._button(i))
            if ic is not None:
                row.addWidget(b)
            b.show()
        if row is not None:
            row.addStretch()
            self.btn_lay.addLayout(row, (k + 1) // 2, 0, 1, 2)
        self.btn_key = tuple(names)
        self.btn_bar.setVisible(bool(names))

    def _button(self, i):
        if i < len(self.btn_fns):
            fn = self.btn_fns[i]
            _later(lambda: (fn(), self.s.refresh()))

    def _build_checks(self, label, names, own=(), side=()):
        """own[i]: kutucuk etiketli satira degil, altinda kendi satirina (uzun adlar sigmiyordu). side: kutucuklarin
        sagina alt alta dizilen dugme adlari."""
        while self.check_lay.count():
            w = self.check_lay.takeAt(0).widget()
            if w is not None:
                w.hide()
                w.setParent(None)
                w.deleteLater()
        self.checks = []
        if label:
            lb = QtWidgets.QLabel(T(label))
            self.check_lay.addWidget(lb, 0, 0)
            lb.show()
        self.check_lay.setColumnStretch(0, 1 if label else 0)
        cols = 1 + sum(not (i < len(own) and own[i]) for i in range(len(names)))
        row = 1
        for i, n in enumerate(names):
            cb = QtWidgets.QCheckBox(T(n))
            cb.setFocusPolicy(QtCore.Qt.NoFocus)           # kisayollar gorunumde kalsin
            cb.toggled.connect(lambda on, i=i: self._check(i, on))
            if label and i < len(own) and own[i]:
                self.check_lay.addWidget(cb, row, 0, 1, cols)
                row += 1
            elif label:
                self.check_lay.addWidget(cb, 0, i + 1)
            else:
                self.check_lay.addWidget(cb, i, 0)
            cb.show()                                      # gorunur panele eklenen cocuk kendiliginden gorunmez;
            self.checks.append(cb)                         # gorunmezse panel yuksekligine girmez (panel ezilirdi)
        if side:
            box = QtWidgets.QWidget()
            col = QtWidgets.QVBoxLayout(box)
            col.setContentsMargins(0, 0, 0, 0)
            col.setSpacing(2)
            for i, n in enumerate(side):
                b = QtWidgets.QPushButton(T(n))
                b.setToolTip(_tool_tip(n))
                b.setFocusPolicy(QtCore.Qt.NoFocus)        # kisayollar gorunumde kalsin
                b.clicked.connect(lambda checked=False, i=i: self._side_button(i))
                col.addWidget(b)
            self.check_lay.addWidget(box, 0, 1, max(len(names), 1), 1, QtCore.Qt.AlignTop | QtCore.Qt.AlignRight)
            box.show()
        self.check_lay.setColumnStretch(1, 1 if side else 0)
        self.check_key = (label, tuple(names), tuple(side))

    def _side_button(self, i):
        if i < len(self.side_fns):
            fn = self.side_fns[i]
            _later(lambda: (fn(), self.s.refresh()))

    def _check(self, i, on):
        if i < len(self.check_fns):
            fn = self.check_fns[i]
            _later(lambda: (fn(on), self.s.refresh()))

    def _commit(self, i, box):
        val = box.value()
        if i >= len(self.fns) or abs(val - self.vals[i]) < 1e-12:   # odak kaybi: deger degismediyse bir sey yapma
            return
        self.vals[i] = val
        fn = self.fns[i]
        _later(lambda: self.s.focus_ok() and self.s.gv.setFocus())     # kisayollar gorunumde calismaya devam etsin
        _later(lambda: fn(val))


_align_icons = {}         # Align dugmelerinin simgeleri (yatay / dikey), bir kez cizilir
_mode_icons = {}          # Point / Edge / Polygon dugmelerinin simgeleri, bir kez cizilir


def _mode_icon(kind):
    """Secim modu simgesi (Blender'daki gibi): tel kup; "point" bir kose, "edge" bir kenar, "face" on yuz sari."""
    if kind in _mode_icons:
        return _mode_icons[kind]
    pm = QtGui.QPixmap(40, 40)                             # 2x cizilir: olcekli ekranda da keskin
    pm.setDevicePixelRatio(2.0)
    pm.fill(QtCore.Qt.transparent)
    p = QtGui.QPainter(pm)
    p.setRenderHint(QtGui.QPainter.Antialiasing)
    P = QtCore.QPointF
    front = [P(2.5, 7.5), P(12.5, 7.5), P(12.5, 17.5), P(2.5, 17.5)]
    back = [P(q.x() + 5, q.y() - 5) for q in front]
    hi = QtGui.QColor("#f5c518")
    if kind == "face":
        p.setPen(QtCore.Qt.NoPen)
        p.setBrush(QtGui.QColor(245, 197, 24, 170))
        p.drawPolygon(QtGui.QPolygonF(front))
    p.setPen(QtGui.QPen(QtGui.QColor("#a8a8a8"), 1.0))
    p.setBrush(QtCore.Qt.NoBrush)
    for sq in (back, front):
        p.drawPolygon(QtGui.QPolygonF(sq))
    for a, b in zip(front, back):
        p.drawLine(a, b)
    if kind == "edge":
        p.setPen(QtGui.QPen(hi, 2.4))
        p.drawLine(front[0], front[1])
    elif kind == "point":
        p.setPen(QtCore.Qt.NoPen)
        p.setBrush(hi)
        p.drawEllipse(front[1], 2.6, 2.6)
    p.end()
    _mode_icons[kind] = QtGui.QIcon(pm)
    return _mode_icons[kind]
_last_tab = ["create"]   # son secilen sekme; modelleyici yeniden acilinca (or. dil degisince) ayni sekme

TAB_CSS = """
QToolButton { background: #5a5a5a; color: #ececec; border: 1px solid #2b2b2b; border-radius: 2px;
              padding: 3px 10px; min-width: 46px; }
QToolButton:hover { background: #6b6b6b; }
QToolButton:checked { background: #3f6fb0; border-color: #24426b; }
"""


class _TabBar(object):
    """Arac panelinin ustunde Create / Modify sekme dugmeleri (LightWave'in ust sekmeleri gibi): panelde hangi
    menunun gorunecegini secer."""
    TABS = (("create", "Create"), ("modify", "Modify"))

    def __init__(self, overlay):
        self.w = QtWidgets.QWidget()
        self.w.setStyleSheet(TAB_CSS)
        lay = QtWidgets.QHBoxLayout(self.w)
        lay.setContentsMargins(4, 4, 4, 3)
        lay.setSpacing(3)
        self.buttons = {}
        for key, text in self.TABS:
            b = QtWidgets.QToolButton()
            b.setText(T(text))
            b.setCheckable(True)
            b.setFocusPolicy(QtCore.Qt.NoFocus)            # kisayollar gorunumde kalsin
            b.clicked.connect(lambda checked=False, k=key: (overlay.set_tab(k), overlay.refocus()))
            lay.addWidget(b)
            self.buttons[key] = b
        lay.addStretch(1)

    def sync(self, tab):
        for k, b in self.buttons.items():
            b.setChecked(k == tab)


class _ToolDock(QtWidgets.QDockWidget):
    """Arac sutununu tasiyan FreeCAD paneli: kenara yerlesir, sekme olur ya da yuzer; kapatilamaz (araclar
    kaybolmasin, modelleyiciyle gelir gider). Yeri ve yuzme durumu tercihlerde saklanir. Paneli surukleyince pencere
    odagi ona gecer: kisayollar calismaya devam etsin diye _KeyFilter buraya da takilir, odak gorunume geri verilir."""

    def __init__(self, overlay, content):
        super().__init__(T("SubD tools"), Gui.getMainWindow())
        self.setObjectName("SubDModelerTools")
        self.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable | QtWidgets.QDockWidget.DockWidgetFloatable)
        self.setWidget(content)
        self.topLevelChanged.connect(lambda floating: overlay.refocus())
        mw = self.parent()
        owned = [getattr(getattr(v, "overlay", None), "tool_dock", None) for v in _sessions.values()]
        for old in mw.findChildren(QtWidgets.QDockWidget, "SubDModelerTools"):   # duzgun kapanmamis oturumdan
            if old is not self and not any(old is o for o in owned):             # baska belgenin oturumununki kalir
                _quiet(lambda old=old: mw.removeDockWidget(old))
                _quiet(old.deleteLater)
        p = _param()
        area = QtCore.Qt.RightDockWidgetArea if p.GetInt("ToolDockArea", 1) == 2 else QtCore.Qt.LeftDockWidgetArea
        mw.addDockWidget(area, self)
        if p.GetBool("ToolDockFloating", False):
            self.setFloating(True)
            x, y, w, h = (p.GetInt("ToolDock" + n, -1) for n in "XYWH")
            if w > 0 and h > 0 and QtWidgets.QApplication.screenAt(QtCore.QPoint(x + 20, y + 20)) is not None:
                self.setGeometry(x, y, w, h)               # ekran sokulduysa eski yere acilip kaybolmasin
        else:
            near = [d for d in mw.findChildren(QtWidgets.QDockWidget) if d is not self and d.isVisible()
                    and not d.isFloating() and mw.dockWidgetArea(d) == area]
            if near:                                       # Model / Tasks yaninda sekme: tam yukseklik
                mw.tabifyDockWidget(near[0], self)
        self.show()
        self.raise_()

    def remember(self):
        p = _param()
        p.SetBool("ToolDockFloating", self.isFloating())
        if self.isFloating():
            g = self.geometry()
            for n, v in zip("XYWH", (g.x(), g.y(), g.width(), g.height())):
                p.SetInt("ToolDock" + n, v)
        else:
            p.SetInt("ToolDockArea", 2 if self.parent().dockWidgetArea(self) == QtCore.Qt.RightDockWidgetArea else 1)


class _Overlay(object):
    def __init__(self, session):
        self.s = s = session
        self.vp = s.gv.viewport()
        self.toggles = {}    # ad -> (dugme, durum fonksiyonu)
        self.gated = []      # (dugme, secim gerekir mi, yalniz hangi modda: None / "edit" / "object")
        self.tool_buttons = []  # kalici arac dugmeleri (Knife, Select loop): aciksa sari, yazi koyu

        # Arac sutunu bolumlerden olusur (baslik + dugmeler). Bolumler LightWave gibi sekmelere ayrilir: panelin
        # ustundeki Create / Modify dugmeleri hangisinin gorunecegini secer. Secili sekme panele sigmazsa relayout
        # bolumlerini yan yana sutunlara akitir. Alt cubuk, islem paneli ve mesajlar gorunumde kalir.
        self.left = QtWidgets.QFrame()
        self.left.setObjectName("lwPanel")
        self.left.setStyleSheet(OVERLAY_CSS)
        self.cols_lay = QtWidgets.QHBoxLayout(self.left)
        self.cols_lay.setContentsMargins(4, 4, 4, 4)
        self.cols_lay.setSpacing(4)
        self.cols = []       # sutun widget'lari (gerektikce eklenir)
        self.sections = []   # (sekme, bolum widget'i), sirayla
        self.groups = None   # son yerlesim: (sekme, her sutundaki bolum sayisi)
        self.tab = _last_tab[0]
        cur = {}

        def head(t, tab="create"):
            sec = QtWidgets.QWidget(self.left)             # baslangicta da panelin CSS'i gecerli olsun (yoksa ilk
            sec.hide()                                     # olcumde bolumler uzun cikar, sutun gereksiz bolunurdu)
            lay = cur["lay"] = QtWidgets.QVBoxLayout(sec)
            lay.setContentsMargins(0, 0, 0, 0)
            lay.setSpacing(2)
            lb = QtWidgets.QLabel(T(t))
            lb.setObjectName("lwHead")
            lb.setAlignment(QtCore.Qt.AlignCenter)
            lay.addSpacing(4)
            lay.addWidget(lb)
            self.sections.append((tab, sec))

        def btn(text, key, fn, toggle=None, sel=False, mode=None, tool=None):
            b = _lw_button(T(text), key, 112)
            b.setToolTip(_tool_tip(text))
            b.clicked.connect(lambda checked=False: (fn(), self.sync(), self.refocus()))
            if tool:                                     # kalici arac: acikken sari (mavi = acik modlar)
                b.setProperty("tool", True)
                toggle = lambda: s.armed == tool
                self.tool_buttons.append(b)
            if toggle:
                b.setCheckable(True)
                self.toggles[text] = (b, toggle)
            if sel or mode:
                self.gated.append((b, sel, mode))
            cur["lay"].addWidget(b)

        def menu_btn(text, items):
            """Acilir dugme (LightWave alt menusu gibi): tiklayinca secenekler dugmenin saginda acilir."""
            b = _lw_button(T(text), "▸", 112)
            b.setToolTip(_tool_tip(text))
            m = QtWidgets.QMenu(b)
            m.setStyleSheet(MENU_CSS)
            m.setToolTipsVisible(True)
            for label, fn in items:
                a = m.addAction(T(label))
                a.setToolTip(_tool_tip(text + ">" + label))
                a.triggered.connect(lambda checked=False, fn=fn: (fn(), self.sync(), self.refocus()))
            b.clicked.connect(lambda checked=False: m.popup(b.mapToGlobal(QtCore.QPoint(b.width() + 2, 0))))
            cur["lay"].addWidget(b)

        obj = lambda: s.mode == "object"
        head("File")
        btn("Save", "F", s.save, mode="edit")
        btn("Solid (STEP)", "", lambda: s.object_solids() if obj() else s.make_solid())
        btn("Close", "", lambda: _later(s.close))
        head("Create")
        btn("Box", "", s.add_box)
        menu_btn("Cylinder", [("Cylinder", lambda: s.add_cylinder()), ("Topo cylinder", lambda: s.add_cylinder(True))])
        menu_btn("From face", [("Face", s.from_face), ("Edge", s.from_edges),   # FreeCAD'de secili Part /
                               ("Point", s.from_points)])                     # PartDesign elemanlarindan
        btn("Make polygon", "P", s.make_polygon, sel=True, mode="edit")
        head("Object")
        btn("Duplicate", "D", s.object_duplicate, sel=True, mode="object")
        head("Select")
        btn("Select all", "A", lambda: s.key("a"))
        btn("Deselect", "Esc", lambda: s.key("escape"))
        btn("Select loop", "L", lambda: s.toggle_tool("loop"), mode="edit", tool="loop")
        head("Polygon")
        btn("Extrude", "E", lambda: s.extrude(lazy=True), sel=True, mode="edit")
        btn("Bevel", "B", lambda: s.edge_bevel(lazy=True) if s.mode == "edge" else s.bevel(lazy=True),
            sel=True, mode="edit")
        btn("Inset", "I", lambda: s.bevel(inset_only=True, lazy=True), sel=True, mode="edit")
        btn("Knife (loop)", "K", lambda: s.toggle_tool("knife"), mode="edit", tool="knife")
        btn("Knife (line)", "J", lambda: s.toggle_tool("line"), mode="edit", tool="line")
        btn("Add points", "", s.add_points, mode="edit", tool="addpoint")
        btn("Connect", "^L", s.connect, sel=True, mode="edit")
        btn("Subdivide", "⇧D", s.subdivide, sel=True, mode="edit")
        btn("Array", _key_label("array"), s.array, sel=True, mode="edit")   # Radial / Rectangular panelden secilir
        menu_btn("Bridge", [("Bridge", s.bridge), ("Rail bridge", s.rail_bridge)])
        btn("Hole", "", s.hole, sel=True, mode="edit")
        btn("Circle on face", "", lambda: s.hole(cap=True), sel=True, mode="edit")
        btn("Merge polygons", "⇧Z", s.merge_polygons, sel=True, mode="edit")
        btn("Parametric group", "", s.make_group, sel=True, mode="edit")
        btn("Delete", "Del", lambda: s.object_delete() if obj() else s.delete(), sel=True)
        head("Clipboard")
        btn("Copy", "^C", s.copy_faces, sel=True, mode="edit")
        btn("Cut", "^X", lambda: s.copy_faces(cut=True), sel=True, mode="edit")
        btn("Paste", "^V", s.paste_faces, mode="edit")
        btn("Paste as object", "", s.paste_as_object)
        head("Modify", tab="modify")
        btn("Move", "T", lambda: s.set_gizmo("move"), toggle=lambda: s.gizmo.mode == "move")
        btn("Rotate", "R", lambda: s.set_gizmo("rotate"), toggle=lambda: s.gizmo.mode == "rotate")
        btn("Stretch / Size", "H", lambda: s.set_gizmo("scale"), toggle=lambda: s.gizmo.mode == "scale")
        btn("Magnet", "", lambda: s.toggle_tool("magnet"), mode="edit", tool="magnet")
        btn("Twist", _key_label("twist"), s.twist, sel=True, mode="edit")    # tus atanabilir (set_user_key)
        btn("Circle", "", s.circle, sel=True, mode="edit")
        btn("Crease", "C", s.crease_selected, sel=True, mode="edit")
        btn("Reduce points", "", s.reduce_points, mode="edit")   # secim yoksa butun kafes
        btn("Merge points", "M", s.merge_points, mode="edit")    # secim yoksa butun kafes
        btn("Weld", "^W", s.weld, sel=True, mode="edit")         # son tiklanan noktada birlesir
        btn("Flip", "⇧F", s.flip_faces, mode="edit")             # secim yoksa butun cokgenler
        btn("Shell", "", s.shell)                     # Object: secili nesne; Edit: secili cokgenler (yoksa hepsi)
        # kafeste: secili cokgenler (secim yoksa hepsi) ayni kafese aynalanir; Object modunda: secili nesnelerin kopyasi
        menu_btn("Mirror", [("Mirror X", lambda: s.mirror_polygons(0)), ("Mirror Y", lambda: s.mirror_polygons(1)),
                            ("Mirror Z", lambda: s.mirror_polygons(2))])
        btn("Mirror (line)", "⇧V", lambda: s.toggle_tool("mline"), mode="edit", tool="mline")

        self.bar = QtWidgets.QFrame(self.vp)
        self.bar.setObjectName("lwBar")
        self.bar.setStyleSheet(OVERLAY_CSS)
        # iki grup: gecmis ve gorunum anahtarlari; sigmayinca ikinci grup alt satira gecer (_fit_bar)
        self.bar_lay = QtWidgets.QBoxLayout(QtWidgets.QBoxLayout.LeftToRight, self.bar)
        self.bar_lay.setContentsMargins(6, 3, 6, 3)
        self.bar_lay.setSpacing(3)
        groups = []
        for stretch in (0, 1):
            g = QtWidgets.QWidget()
            g.setObjectName("lwGroup")
            g.setStyleSheet("QWidget#lwGroup { background: transparent; }")   # tema QWidget zeminini boyamasin
            row = QtWidgets.QHBoxLayout(g)
            row.setContentsMargins(0, 0, 0, 0)
            row.setSpacing(3)
            self.bar_lay.addWidget(g, stretch)
            groups.append(row)
        row = groups[0]
        self.history_buttons = []
        for name, key, redo in (("Undo", "^Z", False), ("Redo", "^Y", True)):
            b = _lw_button(T(name), key)
            b.setToolTip(_tool_tip(name))
            b.clicked.connect(lambda checked=False, redo=redo: (s.step_history(redo), self.sync(), self.refocus()))
            row.addWidget(b)
            self.history_buttons.append((b, redo))
        self.bar_gap = QtWidgets.QSpacerItem(10, 0)        # gruplar arasi (yan yanayken)
        self.bar_lay.insertItem(1, self.bar_gap)

        def tog(text, key, fn, state, b=None):
            b = b or _lw_button(text, key)
            b.setToolTip(_tool_tip(text))
            b.setCheckable(True)
            b.clicked.connect(lambda checked=False: (fn(), self.sync(), self.refocus()))
            self.toggles[text] = (b, state)
            row.addWidget(b)
            return b

        # sol ust (Blender gibi): Object / Edit ve secim modlari, yaninda secim sayisi; alt cubuk tek satira sigsin
        self.modes = QtWidgets.QFrame(self.vp)
        self.modes.setObjectName("lwBar")
        self.modes.setStyleSheet(OVERLAY_CSS)
        row = QtWidgets.QHBoxLayout(self.modes)
        row.setContentsMargins(4, 3, 8, 3)
        row.setSpacing(0)                                  # iki parcali dugmeler bitisik
        tog("Object", "", lambda: s.mode == "object" or s.enter_object_mode(), lambda: s.mode == "object",
            _lw_button(T("Object"))).setToolTip(_tool_tip("Object") + " (O)")
        tog("Edit", "", lambda: s.mode != "object" or s.edit_object(), lambda: s.mode != "object",
            _lw_button(T("Edit"))).setToolTip(_tool_tip("Edit") + " (O)")
        row.addSpacing(8)
        for text, key, mode, tip in (("Point", "^G", "point", "Ctrl+G"), ("Edge", "2", "edge", "2"),
                                     ("Polygon", "^H", "face", "Ctrl+H")):
            b = QtWidgets.QPushButton()
            b.setFocusPolicy(QtCore.Qt.NoFocus)
            b.setIcon(_mode_icon(mode))
            b.setIconSize(QtCore.QSize(20, 20))
            b.setFixedSize(28, 24)
            _min_width(b, 28)
            tog(text, key, lambda m=mode: s.set_mode(m), lambda m=mode: s.mode == m, b)
            b.setToolTip("%s (%s)" % (_tool_tip(text), tip))   # Object modunda: secili nesneye o modda girer
        row.addSpacing(10)
        db = _lw_button(T("Distance"), _key_label("dim"))   # ölçü: her zaman gözün önünde, iki modda da
        db.setProperty("tool", True)
        self.tool_buttons.append(db)
        tog("Distance", _key_label("dim"), s.dim, lambda: s.armed == "dim", db).setToolTip(_tool_tip("Distance"))
        lb = QtWidgets.QPushButton("\u2261")                # Distance listesi (adlar, hedefler, goster / gizle, CSV)
        lb.setFocusPolicy(QtCore.Qt.NoFocus)
        lb.setFixedSize(28, 24)
        _min_width(lb, 28)
        tog("Distances", "", s.dim_list_toggle, lambda: s.dim_list is not None and s.dim_list.isVisible(),
            lb).setToolTip(_tool_tip("Distances"))
        row.addSpacing(10)
        self.sel_label = QtWidgets.QLabel()
        self.sel_label.setMinimumWidth(52)
        row.addWidget(self.sel_label)
        self.modes.show()

        row = groups[1]
        tog("Symmetry", "Y", s.toggle_sym, lambda: s.sym)
        tog("Subpatch", "Tab", s.toggle_smooth, lambda: s.smooth)
        tog("Ghost", "X", s.toggle_ghost, lambda: s.ghost)
        tog("Snap", "Q", s.toggle_snap, lambda: s.snap)
        mc = _lw_button(T("Matcap"), "▴")
        mc.setCheckable(True)
        mc.setToolTip(T("Surface look of the smooth preview (Edit mode)"))
        self.toggles["Matcap"] = (mc, lambda: bool(s.matcap) and not s.zebra and not s.env_on)
        mm = QtWidgets.QMenu(mc)
        mm.setStyleSheet(MENU_CSS)
        acts = []
        for key in ("",) + tuple(MATCAPS):
            a = mm.addAction(T(key or "Off"))
            a.setCheckable(True)
            a.triggered.connect(lambda checked=False, key=key: (s.set_matcap(key), self.sync(), self.refocus()))
            acts.append((a, key))
        mm.aboutToShow.connect(lambda: [a.setChecked(k == s.matcap) for a, k in acts])
        mm.aboutToHide.connect(self.sync)                  # iptal: tiklamayla degisen dugme durumu geri gelsin
        mc.clicked.connect(lambda checked=False: mm.popup(mc.mapToGlobal(QtCore.QPoint(0, -mm.sizeHint().height() - 2))))
        row.addWidget(mc)
        for vertical, tip in ((False, "Align horizontally: selected points (polygons, edges: their centers) get the "
                                      "height of the first clicked one on screen"),
                              (True, "Align vertically: selected points (polygons, edges: their centers) line up "
                                     "with the first clicked one on screen")):
            b = QtWidgets.QPushButton()                    # kucuk kare simgeli dugme (Options'taki ile ayni is)
            b.setFocusPolicy(QtCore.Qt.NoFocus)
            b.setIcon(s._align_icon(vertical))
            b.setIconSize(QtCore.QSize(16, 16))
            b.setFixedSize(22, 22)
            _min_width(b, 22)
            b.setToolTip(T(tip))
            b.clicked.connect(lambda checked=False, v=vertical: (s.align_points(v), self.sync()))
            row.addWidget(b)
        self.cc_label = QtWidgets.QLabel("Catmull-Clark")
        self.cc_label.setStyleSheet("padding-left: 10px;")
        row.addWidget(self.cc_label)
        row.addSpacing(14)
        self.step_label = QtWidgets.QLabel()
        self.step_label.setMinimumWidth(
            self.step_label.fontMetrics().horizontalAdvance(T("Step: %g mm / %g°") % (10.5, 15.5)) + 8)
        row.addWidget(self.step_label)
        row.addStretch(1)
        # bilgi / hata mesaji alt cubukta sigmiyordu: gorunumun icinde, sol menunun saginda, alt cubugun ustunde
        self.info_label = QtWidgets.QLabel(self.vp)
        self.info_label.setObjectName("lwInfo")
        self.info_label.setStyleSheet(OVERLAY_CSS)
        self.info_label.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents)
        tog("Options", "n", s.toggle_options, lambda: s.dock.isVisible())
        lang = _lw_button(subd_i18n.language().upper(), "")    # arayuz dili: EN <-> TR
        lang.setToolTip(T("Switch the interface language (reopens the modeler)"))
        lang.clicked.connect(lambda checked=False: _later(s.switch_language))
        row.addWidget(lang)

        self.panel = _ToolPanel(self.vp, s)
        self.tabs = _TabBar(self)
        self.tabs.sync(self.tab)

        self.scroll = QtWidgets.QScrollArea()              # panel kisalinca sutun panel boyunu zorlamasin
        self.scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scroll.setFocusPolicy(QtCore.Qt.NoFocus)
        self.scroll.setWidget(self.left)
        content = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(content)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(self.tabs.w)
        lay.addWidget(self.scroll, 1)
        self.tool_dock = _ToolDock(self, content)
        s.key_targets.append(self.tool_dock)
        self.tool_dock.installEventFilter(s.key_filter)    # boyut degisince sutunlar da yeniden akar

        self.left.show()
        self.bar.show()
        self.bar.raise_()
        self.relayout()
        self.refocus()
        self._ov_timer = QtCore.QTimer(self.vp)            # Overlay panelleri acilip kapaninca / boyut degisince kay
        self._ov_timer.timeout.connect(lambda: _quiet(self._ov_tick))
        self._ov_timer.start(400)

    def _ov_tick(self):
        if not self.bar.isVisible():
            return
        ov = self._overlay_rects()
        if sorted((k, r.x(), r.y(), r.width(), r.height()) for k, r in ov.items()) != getattr(self, "_ov_sig", None):
            self.relayout()

    def refocus(self):
        """Panele tiklandi: klavye yine gorunume gitsin (panel yuzuyorsa pencere odagi da ana pencereye donsun)."""
        def back():
            if not self.s.focus_ok():                      # askidayken ya da baska belge aktifken odagi calma
                return
            self.s.gv.window().activateWindow()
            self.s.gv.setFocus()
        _later(lambda: _quiet(back))

    def _fit_bar(self, w):
        """Alt cubuk gorunume sigmazsa (125 % olcekli 1920 px ekranda sigmiyordu, yazilar kesiliyordu; modlar
        sol uste tasindi) sirayla: Catmull-Clark yazisi, dugmelerdeki kisayol yazilari (ipucunda kalir), Step
        yazisi gizlenir; yine sigmazsa gorunum anahtarlari Undo / Redo'nun altindaki satira gecer."""
        if getattr(self, "_bar_fit_w", None) == w:
            return
        self._bar_fit_w = w
        buttons = [b for b in self.bar.findChildren(QtWidgets.QPushButton) if getattr(b, "lw_key", None)
                   and b.lw_key[0].text() != "▴"]
        lay = self.bar_lay

        def apply(level):
            self.cc_label.setHidden(level >= 1)
            for b in buttons:
                lab, full, short = b.lw_key
                if lab.isHidden() != (level >= 2) or not hasattr(b, "lw_tip"):
                    tip = getattr(b, "lw_tip", None)
                    if tip is None:
                        tip = b.lw_tip = b.toolTip()
                    lab.setHidden(level >= 2)
                    _min_width(b, full if level < 2 else short)
                    key = lab.text()
                    b.setToolTip(tip if level < 2 else ("%s (%s)" % (tip, key) if tip else key))
            self.step_label.setHidden(level >= 3)
            two = level >= 4
            lay.setDirection(QtWidgets.QBoxLayout.TopToBottom if two else QtWidgets.QBoxLayout.LeftToRight)
            self.bar_gap.changeSize(0 if two else 10, 0)
            lay.invalidate()
            lay.activate()
            return lay.minimumSize().width()

        for level in range(5):
            if apply(level) <= w:
                break

    OVERLAY_NAMES = ("OverlayLeft", "OverlayRight", "OverlayTop", "OverlayBottom")

    def _overlay_rects(self):
        """FreeCAD'in Overlay modundaki paneller (agac / ozellikler gorunumun ustunde yuzer): gorunumle kesisen acik
        olanlarin dikdortgenleri, gorunum koordinatinda. Cubuklarimiz onlarin altinda kalmasin diye kenara kayar."""
        out = {}
        try:
            mw = Gui.getMainWindow()
            vr = QtCore.QRect(self.vp.mapToGlobal(QtCore.QPoint(0, 0)), self.vp.size())
            for name in self.OVERLAY_NAMES:
                w = mw.findChild(QtWidgets.QTabWidget, name)
                if w is None or not w.isVisible():
                    continue
                r = QtCore.QRect(w.mapToGlobal(QtCore.QPoint(0, 0)), w.size()).intersected(vr)
                if not r.isEmpty():
                    out[name] = r.translated(-vr.topLeft())
        except Exception:
            return {}
        return out

    def relayout(self):
        w, h = self.vp.width(), self.vp.height()
        ov = self._overlay_rects()
        self._ov_sig = sorted((k, r.x(), r.y(), r.width(), r.height()) for k, r in ov.items())
        left = ov["OverlayLeft"].right() + 1 if "OverlayLeft" in ov else 0
        right = w - ov["OverlayRight"].left() if "OverlayRight" in ov else 0
        bottom = h - ov["OverlayBottom"].top() if "OverlayBottom" in ov else 0
        left, right = max(0, left), max(0, right)
        self._ins = (left, right, bottom)
        bw = max(240, w - left - right)
        self._fit_bar(bw)
        self.bar.setFixedWidth(bw)
        self.bar.adjustSize()
        self.bar.move(left, h - self.bar.height() - max(0, bottom))
        self.modes.adjustSize()                            # Object / Edit ve secim modlari: sol ust kose
        my = 8
        if "OverlayTop" in ov:                             # ustteki panel modlarin yerini kapliyorsa altina
            mr = QtCore.QRect(left + 8, 8, self.modes.width(), self.modes.height())
            if mr.intersects(ov["OverlayTop"]):
                my = ov["OverlayTop"].bottom() + 9
        self.modes.move(left + 8, my)
        self.modes.raise_()
        avail = self.scroll.viewport().height() - 8        # cerceve kenarlari
        self._flow_columns(avail if avail > 100 else 1 << 30)   # sekmesi arkada / henuz acilmamis panel: tek sutun
        bar = self.scroll.verticalScrollBar().sizeHint().width() if self.left.height() > avail > 92 else 0
        self.scroll.setMinimumWidth(self.left.width() + bar)
        # islem paneli: gorunumun solu ya da sagi (Options > Tool panel), alt cubugun hemen ustu
        py = h - self.bar.height() - max(0, bottom) - self.panel.height() - 8
        if self.panel_on_right():
            self.panel.move(max(left + 8, w - right - self.panel.width() - 8), py)
        else:
            self.panel.move(8 + left, py)
        self._place_info()

    def panel_on_right(self):
        return _param().GetString("PanelSide", "left") == "right"

    def _place_info(self):
        """Mesaj etiketi: gorunumun solunda, alt cubugun hemen ustunde; panel solda aciksa panelin ustunde, sagdaysa
        etiket altta kalir (panele yer birakir); bos ise gizli."""
        lab = self.info_label
        if not lab.text():
            lab.hide()
            return
        ins = getattr(self, "_ins", (0, 0, 0))
        x = 8 + ins[0]
        maxw = max(200, self.vp.width() - x - ins[1] - 12)
        panel_right = self.panel.isVisible() and self.panel_on_right()
        if panel_right:                                    # etiket sagdaki panelin altina girmesin
            maxw = max(200, self.panel.x() - x - 8)
        lab.setMaximumWidth(maxw)
        lab.setWordWrap(False)                             # sigiyorsa tek satir; sigmiyorsa satira bol
        lab.adjustSize()
        if lab.sizeHint().width() > maxw:
            lab.setWordWrap(True)
            lab.setFixedWidth(maxw)
            lab.adjustSize()
            lab.setMinimumWidth(0)
            lab.setMaximumWidth(maxw)
        bottom = self.vp.height() - self.bar.height() - max(0, ins[2]) - 8
        if self.panel.isVisible() and not panel_right:
            bottom = self.panel.y() - 6
        lab.move(x, bottom - lab.height())
        lab.show()
        lab.raise_()

    def set_tab(self, tab):
        """Create / Modify sekmesi: sol sutunda yalniz o sekmenin bolumleri (LightWave sekmeleri gibi)."""
        self.tab = _last_tab[0] = tab
        self.relayout()
        if self.tabs is not None:
            self.tabs.sync(tab)

    def _flow_columns(self, avail):
        """Secili sekmenin bolumlerini sirayla sutunlara yerlestir: sutun avail yuksekligi asacaksa sonraki
        bolum yeni sutuna. Oteki sekmelerin bolumleri gizlenir."""
        groups, used = [], None
        for tab, sec in self.sections:
            if tab != self.tab:
                continue
            sec.ensurePolished()
            hs = sec.sizeHint().height()
            if used is None or (used + hs > avail and groups[-1]):
                groups.append([])
                used = 0
            groups[-1].append(sec)
            used += hs
        key = (self.tab, [len(g) for g in groups])
        if key == self.groups:
            return
        self.groups = key
        for _, sec in self.sections:                       # eski sutunlardan cikar
            p = sec.parentWidget()
            if p is not None and p.layout() is not None:
                p.layout().removeWidget(sec)
            sec.hide()
        while len(self.cols) < len(groups):
            col = QtWidgets.QWidget()
            cl = QtWidgets.QVBoxLayout(col)
            cl.setContentsMargins(0, 0, 0, 0)
            cl.setSpacing(0)
            cl.addStretch(1)
            self.cols_lay.addWidget(col, 0, QtCore.Qt.AlignTop)
            self.cols.append(col)
        for col, g in zip(self.cols, groups):
            cl = col.layout()
            for sec in g:
                cl.insertWidget(cl.count() - 1, sec)
                sec.show()
            col.show()
        for col in self.cols[len(groups):]:
            col.hide()
        self.left.layout().activate()
        self.left.adjustSize()

    def sync(self):
        s = self.s
        if s.suspended:                                    # baska tezgahta: panel / mesaj geri gelmesin
            return
        for b, state in self.toggles.values():
            b.setChecked(bool(state()))
        for b in self.tool_buttons:                        # sari zeminde koyu yazi
            lit = b.isChecked()
            if b.property("lit") != lit:
                b.setProperty("lit", lit)
                name, key = b.findChildren(QtWidgets.QLabel)
                name.setStyleSheet("color:#1a1a1a;" if lit else "")
                key.setStyleSheet("color:#4a3a00;" if lit else "color:#a8a8a8;")
        is_obj = s.mode == "object"
        doc = App.ActiveDocument
        for b, redo in self.history_buttons:
            if is_obj:
                b.setEnabled(doc is not None and bool(doc.RedoCount if redo else doc.UndoCount))
            else:
                b.setEnabled(bool(s.redo if redo else s.undo))
        nobj = len(s.selected_objects()) if is_obj else 0
        has = nobj > 0 if is_obj else bool(s.selected_verts())
        for b, sel, mode in self.gated:
            ok = mode is None or (mode == "object") == is_obj
            b.setEnabled(ok and (has or not sel))
        count = nobj if is_obj else {"point": len(s.sel_v), "edge": len(s.sel_e), "face": len(s.sel)}[s.mode]
        self.sel_label.setText("Sel: %d" % count)
        self.step_label.setText(T("Step: %g mm / %g°") % (s.step_box.value(), s.angle_box.value()))
        info_changed = self.info_label.text() != s.info
        self.info_label.setText(s.info)
        if self.panel.show_spec(s.panel_spec()):
            self.relayout()
        elif info_changed:
            self._place_info()

    def set_visible(self, on):
        """Calisma tezgahi degisince arac paneli ve gorunumdeki cubuklar gizlenir / geri gelir (oturum kalir)."""
        if on:
            self.tool_dock.show()
            self.bar.show()
            self.modes.show()
            self.relayout()
            self.sync()
        else:
            for w in (self.tool_dock, self.bar, self.modes, self.panel, self.info_label):
                _quiet(w.hide)

    def close(self, doc_gone=False):
        _quiet(self._ov_timer.stop)
        d = self.tool_dock                                 # ana pencerede: gorunum silinse de kalir, once kaldir
        _quiet(d.remember)
        _quiet(lambda: d.parent().removeDockWidget(d))
        _quiet(d.deleteLater)
        if doc_gone:                                       # bar / panel ... 3B gorunumun cocuklari: gorunumle gitti, dokunma
            return
        for w in (self.bar, self.modes, self.panel, self.info_label):   # gorunumun cocuklari: silinmis olabilir
            _quiet(w.hide)
            _quiet(w.deleteLater)


_session = None        # su an etkin (gorunen) oturum
_sessions = {}         # belge adi -> oturum: her belgenin kendi oturumu var; baska belgeye gecince askiya alinir


def start(obj=None, object_mode=None):
    """Modelleyiciyi acar. Tek SubD nesnesi secili (ya da verilmis) ise onun kafesini duzenler; yoksa Object
    modunda acar. Bos belgede de bos acar (LightWave gibi; eskiden kutuyla baslardi ve her seferinde silinmesi
    gerekiyordu, issue #7): Create > Box / Cylinder yeni nesneyi kurar."""
    global _session
    if App.ActiveDocument is None or Gui.ActiveDocument is None:
        App.newDocument("SubD")
    doc_name = Gui.ActiveDocument.Document.Name
    old = _sessions.get(doc_name)
    if old is not None:                                    # bu belgenin oturumu zaten var: kapat, yeniden ac
        old.close()
    for o in list(_sessions.values()):                     # baska belgelerin oturumlari kapanmaz, askiya alinir
        if o.doc_name != doc_name:
            _quiet(o.suspend)
    repair_proxies(App.ActiveDocument)
    if object_mode is None:
        if obj is None:
            sel = [o for o in Gui.Selection.getSelection() if is_subd(o)]
            obj = sel[0] if len(sel) == 1 else None
            object_mode = obj is None
        else:
            object_mode = False
    if not object_mode and obj is not None:
        Gui.Selection.clearSelection()
    _session = SubDSession(obj, object_mode=object_mode)
    _sessions[_session.doc_name] = _session
    return _session
