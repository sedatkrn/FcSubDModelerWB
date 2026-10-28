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
#   K: loop cut + slide   J: LightWave knife (draw a line, live plane cut)   C: crease   Del: delete
#   Knife (K / J) and Select loop stay on (yellow button) until K / J / the button / Esc / right click
#   T: move gizmo   G: free move (mouse)   R: rotate gizmo   H: stretch / size gizmo   Y: X symmetry   Tab / S: SubPatch
#   Ctrl+C / Ctrl+X / Ctrl+V: copy / cut / paste polygons (object mode Ctrl+V: as a new SubD object)
#   Ctrl+Z: undo   F: save to the document   N: Options panel   Esc / right click: cancel
#   After every operation the small panel in the viewport shows its values; type a value + Enter and the
#   operation is redone with it. Ctrl while dragging = steps (mm / degrees).
#
# Code comments are mostly in Turkish (the author's language); the user interface is English by default
# with a Turkish option (subd_i18n.py).
import math

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


# ---------- Kafes (poligon mesh) ve duzenleme islemleri ----------

def _edge(a, b):
    return frozenset((a, b))


class Cage:
    def __init__(self, V=None, F=None, C=None):
        self.V = np.zeros((0, 3)) if V is None else V
        self.F = [] if F is None else F
        self.C = {} if C is None else C   # keskin kenarlar: frozenset((a, b)) -> keskinlik

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
        return Cage(self.V.copy(), [list(f) for f in self.F], dict(self.C))

    def edges(self):
        return {_edge(f[i], f[(i + 1) % len(f)]) for f in self.F for i in range(len(f))}

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

    def extrude_edges(self, edges, segs=1, layers=None):
        """Acik (tek yuzlu) kenarlari yuzeyin devami gibi disari uzatir: her kenardan segs dortgen. Ic kenarlar
        atlanir (SubD'de yuzgec / manifold olmayan kenar olurdu). (sirali taban noktalari, birim disari yonleri)
        ya da hic acik kenar yoksa None dondurur; layers extrude'daki gibi doldurulur."""
        directed = self._open_directed(edges)
        if not directed:
            return None
        verts = sorted({v for a, b, _ in directed for v in (a, b)})
        out = self._outward(directed)
        dirs = np.array([out[v] for v in verts])
        L = self._layers(verts, segs)
        for a, b, fi in directed:                          # yeni yuz b -> a sirasiyla: komsu yuzle ayni yonde
            self.F += [[L[b][k], L[a][k], L[a][k + 1], L[b][k + 1]] for k in range(segs)]
            e = _edge(a, b)
            if e in self.C:                                # keskinlik yeni acik kenara gecer
                self.C[_edge(L[a][-1], L[b][-1])] = self.C.pop(e)
        if layers is not None:
            layers.update(L)
        return verts, dirs

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

    def hole(self, region, radius=0.0, segs=0, support=-1.0):
        """Secili yuz bolgesinin ortasina daire delik acar (bolgenin yuzleri yeniden orulur). Bolgenin dis siniri
        (N nokta) delik cevresindeki M noktaya aci sirasiyla baglanir: M = N ise hepsi dortgen; M > N besgen,
        M < N ucgen olur. support > 0: delik kenarina destek halkasi (keskin, SubD'de daire kalir). radius / segs /
        support <= 0 ise otomatik. Daire bolge yuzeyine normal boyunca izdusurulur (hafif egri yuzeyde de olur).
        Donus: (yeni yuz indeksleri, kullanilan radius, segs, support) ya da hata metni (kafes degismez)."""
        if not region:
            return T("Hole: select faces")
        ef = self.edge_faces()
        edges = {e for e, fs in ef.items() if sum(f in region for f in fs) == 1}
        chains = self.edge_chains(edges)
        if not chains or len(chains) != 1 or not chains[0][1]:
            return T("Hole: the selected faces must form one patch without holes")
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
            return T("Hole: the selected faces are not flat")
        n = nsum / np.linalg.norm(nsum)
        c = csum / asum
        if any(self.normal(fi).dot(n) < 0.8 for fi in region):     # ortalamadan ~37 dereceden fazla sapan yuz
            return T("Hole: the selected faces are not flat")
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
        support = support if support >= 0 else 0.2 * radius
        if radius + support > 0.95 * dmin:
            return T("Hole: radius too large for the selected faces (max %.2f mm)") % (0.95 * dmin - support)
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
        self.delete(region)                                # yeni yuzler sonda: indeksleri silinenler kadar kayar
        return list(range(nf - len(region), len(self.F))), radius, segs, support

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

    def face_ring(self, fi, u, v):
        """(u, v) kenarindan gecen dortgen halkasi: {yuz: dondurulmus dortgen (q0q1 giris, q2q3 cikis kenari)}."""
        ef = self.edge_faces()

        def rotated(f, e):
            for i in range(len(f)):
                if {f[i], f[(i + 1) % len(f)]} == e:
                    return f[i:] + f[:i]

        ring = {}

        def walk(f, e):
            while f is not None and f not in ring and len(self.F[f]) == 4:
                q = rotated(self.F[f], e)
                ring[f] = q
                e = {q[2], q[3]}
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
        """Halkadaki her dortgeni giris / cikis kenarlarinin ortasindan ikiye boler."""
        # taraf: dortgende q0-q3 ve q1-q2 ayni tarafta; halka sirasiyla seed'in tarafi yayilir
        side = {seed}
        for q in ring.values():
            if q[1] in side or q[2] in side:
                side |= {q[1], q[2]}
            else:
                side |= {q[0], q[3]}
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
                q = ring[i]
                m, n = mid(q[0], q[1]), mid(q[2], q[3])
                split += [[q[0], m, n, q[3]], [m, q[1], q[2], n]]
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
        used = sorted({x for f in self.F for x in f})
        remap = {x: i for i, x in enumerate(used)}
        self.V = self.V[used]
        self.F = [[remap[x] for x in f] for f in self.F]
        alive = self.edges()
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
                    q = self.V[v] + off[k]
                    sl = -((q - ctr) @ pn) / (d @ pn)      # ic duvar boyunca duzleme
                    if abs(sl) <= 2 * t:
                        off[k] = off[k] + sl * d
        base = self.V[R].copy()
        if mode == 2:                                      # ortada: asil yuzey de yarim disari
            self.V[R] = base + t * k_even[:, None] * Nv
        # kopyada kisa kenarlari birlestir (akilli suzgec): destek halkasi seritleri, kalinliktan dar basamaklar ic
        # tarafta kopyalanmaz. Yoksa keskin katlanmadaki kisa kenarlar kalinlik kadar kayinca kafesler kesisir.
        # Acik kenar (delik agzi, acik alt kenar) ve secim siniri noktalari birlesmez: Fill ve orme bozulmasin.
        pos0, off0 = base + delta, off
        bverts = {x for e in border for x in e}
        fixed = np.array([v in pinned or v in bverts for v in R])
        par = list(range(len(R)))

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        thr = max(0.5 * radius, thickness) if radius > 1e-9 else 0.0
        if thr > 0:
            P0 = pos0 + off0
            lo, hi = P0.copy(), P0.copy()
            cand = sorted((np.linalg.norm(P0[ri[x]] - P0[ri[y]]), ri[x], ri[y]) for x, y in (tuple(e) for e in redges_of(self, region)))
            for ln, x, y in cand:
                if ln >= thr:
                    break
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
        roots = sorted({find(k) for k in range(len(R))})
        rid = {r: j for j, r in enumerate(roots)}
        grp = np.array([rid[find(k)] for k in range(len(R))])
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

    def subset(self, sel):
        """Secili yuzlerden yeni (bagimsiz) kafes: kullanilan noktalar yeniden numaralanir, iki ucu da kalan
        keskin kenarlar keskin kalir. Kopyala / yapistir icin."""
        faces = [self.F[fi] for fi in sorted(sel)]
        used = sorted({v for f in faces for v in f})
        remap = {v: i for i, v in enumerate(used)}
        F = [[remap[v] for v in f] for f in faces]
        sub = Cage(self.V[used].copy() if used else np.zeros((0, 3)), F)
        alive = sub.edges()
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
        return list(range(nf, len(self.F)))

    def delete(self, sel):
        self.F = [f for i, f in enumerate(self.F) if i not in sel]
        used = sorted({v for f in self.F for v in f})
        remap = {v: i for i, v in enumerate(used)}
        self.V = self.V[used] if used else np.zeros((0, 3))
        self.F = [[remap[v] for v in f] for f in self.F]
        alive = self.edges()
        self.C = {_edge(*(remap[v] for v in e)): s for e, s in self.C.items()
                  if all(v in remap for v in e) and _edge(*(remap[v] for v in e)) in alive}


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
    return Cage(V, F, C)


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


def _chain(edges, counts):
    """Sirali kapali kenar zinciri: her kenarin ornegi ZINCIR yonunde ((n+1, 3) dizi, uclar dahil). Kenarin
    parametre yonu zincirle ters olabilir (komsuluktan bulunur; Vertexes sirasina guvenilmez)."""
    ends = [(np.array(tuple(e.valueAt(e.FirstParameter))), np.array(tuple(e.valueAt(e.LastParameter))))
            for e in edges]
    nxt = ends[1 % len(edges)]                             # ilk kenarin sonu ikinci kenara degmeli
    fwd0 = min(np.linalg.norm(ends[0][1] - q) for q in nxt) <= min(np.linalg.norm(ends[0][0] - q) for q in nxt)
    prev = ends[0][0] if fwd0 else ends[0][1]
    out = []
    for e, n, (p0, p1) in zip(edges, counts, ends):
        q = _edge_points(e, n)
        if np.linalg.norm(p1 - prev) < np.linalg.norm(p0 - prev):
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


def _grid_patch(face, edges, counts):
    """Dort kenarli yuz: kenar ornekleri arasinda UV'de Coons izgarasi, ic noktalar yuzeyin uzerinde. (V, F)"""
    segs = _chain(edges, counts)
    n0, n1 = len(segs[0]) - 1, len(segs[1]) - 1
    if len(segs[2]) - 1 != n0 or len(segs[3]) - 1 != n1:
        return None
    B, R, Tp, L = segs[0], segs[1], segs[2][::-1], segs[3][::-1]   # alt, sag, ust (soldan saga), sol (alttan uste)
    S = face.Surface
    uv = lambda P: np.array([tuple(S.parameter(App.Vector(*p))) for p in P], float)
    UB, UR, UT, UL = uv(B), uv(R), uv(Tp), uv(L)
    for k, per in enumerate((S.UPeriod() if S.isUPeriodic() else 0, S.VPeriod() if S.isVPeriodic() else 0)):
        if per:                                            # periyodik yuzey (silindir): dikisin ayni tarafina al
            for A in (UB, UR, UT, UL):
                A[:, k] = UB[0, k] + (A[:, k] - UB[0, k] + per / 2) % per - per / 2
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
    info = []
    for f in faces:
        planar = f.isPlanar() if hasattr(f, "isPlanar") else f.Surface.TypeId == "Part::GeomPlane"
        ws = [f.OuterWire] + [w for w in f.Wires if not w.isSame(f.OuterWire)]
        wl = [[eid(e) for e in w.OrderedEdges] for w in ws]
        if len(ws) == 1 and len(wl[0]) == 4:
            kind = "grid"
        elif planar and len(ws) <= 2:
            kind = "planar"
        elif planar:
            return T("From face: faces with more than one hole are not supported yet")
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
    for f, kind, wl in info:
        res = (_grid_patch(f, [edges[i] for i in wl[0]], [cnt[i] for i in wl[0]]) if kind == "grid" else
               _planar_patch(f, [([edges[i] for i in w], [cnt[i] for i in w]) for w in wl]))
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
        T0, P = cage.V.copy(), cage.V.copy()
        lv = 4 if len(P) < 600 else 3
        for _ in range(20):
            Sx, _ = subdivide(P, cage.F, lv, cage.C)
            P = P + (T0 - Sx[:len(P)])
        cage.V = P
    return cage


def cage_key(obj):
    """Kafesin degisip degismedigini anlamak icin ozet (Placement haric)."""
    return (hash(tuple((p.x, p.y, p.z) for p in obj.CagePoints)), tuple(obj.CageFaces),
            tuple(getattr(obj, "CageCreases", ())), tuple(getattr(obj, "CageSharpness", ())), obj.Levels,
            obj_rules(obj))


def is_subd(obj):
    return obj is not None and "CageFaces" in getattr(obj, "PropertiesList", [])


class SubDObject:
    def __init__(self, obj):
        obj.addProperty("App::PropertyVectorList", "CagePoints", "SubD", "Cage points (mm)")
        obj.addProperty("App::PropertyIntegerList", "CageFaces", "SubD", "Cage faces: point indices, -1 separated")
        obj.addProperty("App::PropertyInteger", "Levels", "SubD", "SubPatch subdivision level").Levels = 3
        _crease_props(obj)
        _surface_props(obj)
        obj.Proxy = self

    def execute(self, obj):
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

GIZMO_PX = 90.0          # ok boyu, ekranda piksel (yakinlasip uzaklasinca sabit kalir)
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
_last_box = {"Width": 100.0, "Depth": 100.0, "Height": 100.0,           # Box paneli son degerleri (mm, dilim)
             "Segments X": 1.0, "Segments Y": 1.0, "Segments Z": 1.0}
def _param():
    """Modelleyici ayarlari (FreeCAD parametreleri; dil ayariyla ayni yer)."""
    return App.ParamGet(subd_i18n.PARAM)


_last_shell = {"Thickness": 2.0, "Detail filter": 6.0, "Mode": 0.0, "Fill": 1.0, "Sharp": 1.0, "Flat rim": 1.0}   # Shell paneli
_knife_opts = {"split": False, "objects": False}   # Knife (line) secenekleri (ikisi kapali: yalniz keser)
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
        size = vv.getWorldToScreenScale(coin.SbVec3f(p.x, p.y, p.z), GIZMO_PX / float(H))
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
        if self.mode == "scale":
            if math.hypot(x - c[0], y - c[1]) < 12.0:
                return ("size", 0)
            best, bd = None, PICK_PX
            for i in range(3):
                e = _to_screen(vv, W, H, self.center + AXES[i] * (self.size * 1.1))
                dd = _seg_dist(x, y, c, e)
                if dd < bd:
                    best, bd = ("stretch", i), dd
            return best
        if math.hypot(x - c[0], y - c[1]) < 10.0:      # oncelik: merkez kup, sonra kareler, sonra oklar
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
        "E: extrude  |  K: knife (loop cut)  |  J: knife (line / plane)\n"
        "Knife / loop stay on (yellow) until K / J / Esc / right click\n"
        "B: bevel (right = inset, up = shift)  |  I: inset\n"
        "Edge mode B: edge bevel (right = width)\n"
        "Edge mode E: extrude open edges  |  extrude: wheel = segments, X / Y / Z = axis\n"
        "O: object / edit mode  |  D: duplicate object\n"
        "T: move gizmo  |  G: free move (mouse)  |  Tab / S: SubPatch\n"
        "A: select all  |  Del: delete  |  P: make polygon\n"
        "Ctrl+C / X / V: copy / cut / paste polygons\n"
        "Ctrl+Z: undo  |  F: save to document\n"
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
            Cage.cube([_last_box[n] for n in self.BOX_SIZE])   # bos belge: Box olculerinde kutu
        self.hidden = []     # duzenlerken gizlenen nesneler (SubD ve ondan uretilen katilar); kapatinca geri
        self._snapshot()
        if not object_mode:
            self._hide_edited()
        self.mode = "object" if object_mode else "face"   # "object" | "point" | "edge" | "face"
        self.edit_mode = "face"                            # Object modundan Edit'e donulecek secim modu
        self.sel = set()     # yuz indeksleri
        self.sel_v = set()   # nokta indeksleri
        self.sel_e = set()   # frozenset({a, b}) kenarlar
        self.sym = False     # X=0 duzlemine gore simetri
        self.ghost = False   # hayalet modu: arkadakiler gorunur / secilebilir
        self.snap = False    # serbest tasimada Part / PartDesign kose ve kenarlarina yakalama
        self.zebra = False   # zebra seritleri (yuzey surekliligi analizi)
        self._snapdata = None
        self.smooth = obj is not None
        self.levels = obj.Levels if obj is not None else 3       # son (belge / kati) seviyesi
        self.rules = obj_rules(obj) if obj is not None else RULES
        self.edit_levels = _param().GetInt("EditLevels", 3)      # duzenlerken gorunen seviye (hiz icin dusuk)
        self.undo = []
        self.splits = []     # Knife "Split into two objects": (geri alma derinligi, yeni nesnenin adi)
        self.grab = None     # klavyeyle tasima (G, extrude)
        self.drag = None     # gizmo tutamagiyla tasima
        self.armed = None    # kalici arac ("knife" / "line" / "loop"): kapatilana kadar her tiklamada uygulanir
        self.op = None       # surmekte olan islem (panelde canli degerler)
        self.ebevel_segs = 1  # kenar bevel segment sayisi (son kullanilan)
        self.extrude_segs = 1        # extrude dilim sayisi (son kullanilan)
        self.box_vals = dict(_last_box)  # Box paneli: son kullanilan degerler (modelleyici kapansa da)
        self.cyl_vals = dict(_last_cyl)  # Cylinder paneli (iki topoloji ortak)
        self.bridge_vals = {"Number of Cuts": 0.0, "Twist": 0.0, "Smoothness": 0.0, "Profile Factor": 0.0}
        self.last_op = None  # biten son islem (panelden degerle yeniden uygulanabilir)
        self.box = None      # bos alana basildi: surukleme kutu secimi olur
        self.rubber = None   # secim kutusu cizimi (QRubberBand)
        self.line_overlay = None  # knife cizgisi cizimi
        self.overlay = None
        self._pending = False
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
        self.timer.start(50)
        self.autosave = QtCore.QTimer()   # degisiklikler belgeye duzenli yazilsin (FreeCAD kapanirsa kaybolmasin)
        self.autosave.timeout.connect(self._autosave)
        self.autosave.start(20000)
        # Klavye Qt duzeyinde yakalanir: Coin'e giden olaylarin bazilarini (Tab, A, E...)
        # FreeCAD tuketiyor, Ctrl+Z gibi tuslar da FreeCAD kisayollarina gidiyor.
        gv = self.gv = _graphics_view(self.view)
        self.key_targets = [gv, gv.viewport()]
        self.key_filter = _KeyFilter(self)
        for w in self.key_targets:
            w.installEventFilter(self.key_filter)
        self.overlay = _Overlay(self)
        self.geo_switch.whichChild = -1 if object_mode else 0
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
        if self.cage.F:
            self.save(commit=True)
        self._show_edited()
        self.drop_tool()
        self.edit_mode = self.mode
        self.mode = "object"
        self.geo_switch.whichChild = -1
        self.op = self.last_op = None
        self.undo, self.splits = [], []
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
        Gui.Selection.clearSelection()
        self.obj = o
        self.cage = cage_from_obj(o)
        self.levels = o.Levels
        self.rules = obj_rules(o)
        self._sync_surface_widgets()
        self.undo, self.op, self.last_op = [], None, None
        self.splits = []
        self.sel, self.sel_v, self.sel_e = set(), set(), set()
        self.dirty = False
        self.set_placement(o.Placement)
        self._snapshot()
        self._hide_edited()
        self.mode = self.edit_mode
        self.geo_switch.whichChild = 0
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
            o.ViewObject.Proxy = 0
        cage_to_obj(cage, o)
        set_surface(o, self.levels, self.rules)
        o.recompute()
        doc.commitTransaction()
        self.edit_mode = "face"
        self.edit_object(o)
        self.push()

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
                n.ViewObject.Proxy = 0
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
                n.ViewObject.Proxy = 0
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
        doc = objs[0].Document
        doc.openTransaction("SubD delete")
        for o in objs:
            for s in [x for x in o.InList if "Samples" in x.PropertiesList]:
                doc.removeObject(s.Name)
            if o is self.obj:
                self.obj = None
            doc.removeObject(o.Name)
        doc.commitTransaction()
        Gui.Selection.clearSelection()
        self.refresh()

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
        self.cage_mat = material((0.55, 0.60, 0.70))
        self.s_switch = coin.SoSwitch()
        self.zebra_switch = coin.SoSwitch()                # yansima seritleri: ortam (kure) dokusu
        self.zebra_switch.addChild(_zebra_texture())
        self.s_switch.addChild(_sep(hints(0.5), self.smooth_mat, offset(), self.zebra_switch, self.s_coord,
                                    self.s_faces))
        self.c_switch = coin.SoSwitch()
        self.c_switch.addChild(_sep(
            _sep(hints(0.0), self.cage_mat, offset(), self.c_coord, self.c_faces),
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
        self.root.addChild(self.geo_switch)
        self.root.addChild(self.gizmo.root)
        self.parent = _scene_parent(self.view)
        self.parent.addChild(self.root)

    def _build_dock(self):
        self.dock = QtWidgets.QDockWidget(T("SubD Options"))
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        self.status = QtWidgets.QLabel()
        # Araclar ve modlar gorunumdeki menulerde; bu Options paneli olculer, ayarlar ve kisayol listesi
        lay.addWidget(self._measure_box())
        lay.addWidget(self._bevel_group())
        lay.addWidget(self._box_group())
        lay.addWidget(self._surface_group())
        lay.addWidget(self._solid_group())
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
                      (self.zebra_box, lambda: b.setChecked(self.zebra))):
            b.blockSignals(True)
            fn()
            b.blockSignals(False)

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

    def refresh(self, fast=False):
        if self.mode == "object":
            if not (self.drag and fast):               # surukleme sirasinda gizmo nesneyle birlikte tasinir
                c = self.objects_center()
                self.gizmo.place(App.Vector(*c) if c is not None else None)
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
            S, Q = subdivide(V, F, min(lv, 2) if fast else lv, self.cage.C, self.rules)
            _set_coords(self.s_coord, S)
            _set_index(self.s_faces, np.concatenate([Q, -np.ones((len(Q), 1), int)], 1).ravel().tolist())
        self.s_switch.whichChild = 0 if self.smooth else -1
        self.c_switch.whichChild = -1 if self.smooth else 0
        verts = self.selected_verts() if not self.grab else []
        self.gizmo.place(App.Vector(*self.l2w(V[verts].mean(0))) if verts else None)
        self._update_measure(fast)
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
        self.last_op = None
        self.refresh()

    # --- islemler ---
    def push(self):
        self.undo.append((self.cage.copy(), set(self.sel), set(self.sel_v), set(self.sel_e)))
        cut = len(self.undo) - 100
        if cut > 0:                                        # en eski adimlar dusuyor: ayirma kayitlari da kayar
            del self.undo[:cut]
            self.splits = [(k - cut, n, c) for k, n, c in self.splits if k - cut > 0]
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
        elif k == "extrude":
            segs = min(max(int(round(v["Segments"])), 1), 64)
            v["Segments"] = float(segs)
            self.extrude_segs = segs
            st = self._extrude_build(op["before"], segs, a["axes"])
            if st is not None:
                st["op"] = op
                self.extrude_set(st, v)
        elif k == "bevel":
            self.bevel(inset_only=a["inset_only"], values=(v["Inset"], v.get("Shift", 0.0)), push=False)
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
        elif k == "hole":
            err = self._hole(v)
            if err:
                self.info = err
        elif k == "shell":
            err = self._shell(v, a["faces"])
            if err:
                self.info = err
        elif k == "fromface":
            c = self._face_cage(a["faces"], v)
            if isinstance(c, str):
                self.info = c
                return
            c.V = self.w2l(c.V)
            if a["append"]:
                new = self.cage.append(c)
            else:
                self.cage, new = c, list(range(len(c.F)))
            self.sel, self.sel_v, self.sel_e = set(new) if self.mode == "face" else set(), set(), set()
            self.info = T("From face: %d polygons, %d points") % (len(c.F), len(c.V))
        elif k == "plane":
            n = np.array([v["Nx"], v["Ny"], v["Nz"]], float)
            if np.linalg.norm(n) > 1e-9:
                a["side"] = self._plane_cut(n / np.linalg.norm(n), v["Offset"], a["verts"],
                                           _knife_opts["split"] or _knife_opts["objects"])
            self.sel, self.sel_v, self.sel_e = set(), set(), set()   # clear_selection son islemi de silerdi
        elif k == "knife":
            items, maxt = self._knife_cut(*a["edge"])
            self.slide_set(items, min(max(v["t"], 0.0), maxt))

    def panel_spec(self):
        """Viewport'taki arac panelinin icerigi: (baslik, [(ad, deger, birim, ondalik, geri cagri)], duzenlenebilir)."""
        if self.mode == "object":
            return self.object_panel_spec()
        live = bool(self.grab or self.drag)
        op = self.op if live else self.last_op
        if self.armed == "line" and not live and (op is None or op["kind"] != "plane"):
            return T("Knife (line)"), [], True, self._knife_checks()   # cizmeden once secenekler
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
            if k == "extrude":
                axes = a["axes"]
                fields = [mm(n) for n in "XYZ" if n in axes] if axes else [mm("Distance")]
                fields.append(("Segments", v["Segments"], "", 0, r("Segments")))
                checks = [(n, n in axes, lambda on, n=n: self.extrude_axis(n, on)) for n in "XYZ"]
                return "Extrude", fields, self.drag is None, checks, "Axis (none = normal)"
            if k == "bevel":
                return ("Inset" if a["inset_only"] else "Bevel"), \
                    [mm("Inset")] + ([] if a["inset_only"] else [mm("Shift")]), self.drag is None
            if k == "ebevel":
                return T("Edge bevel"), [mm("Width"), ("Segments", v["Segments"], "", 0, r("Segments"))], \
                    self.drag is None
            if k == "mirror":
                return T("Mirror ") + "XYZ"[a["axis"]], [mm("Center")], self.drag is None
            if k == "fromface":
                return T("From face"), [("Arc Segments", v["Arc Segments"], "", 0, r("Arc Segments")),
                                        mm("Max Length"), ("Exact Edge", v["Exact Edge"], "", 0, r("Exact Edge"))], \
                    self.drag is None
            if k == "shell":
                mode = int(v["Mode"])
                checks = [(n, mode == i, lambda on, n=n: self.shell_option(n, on)) for i, n in enumerate(self.SHELL_MODES)]
                checks += [("Fill", bool(v["Fill"]), lambda on: self.shell_option("Fill", on)),
                           ("Sharp", bool(v["Sharp"]), lambda on: self.shell_option("Sharp", on)),
                           ("Flat rim", bool(v.get("Flat rim", 1.0)), lambda on: self.shell_option("Flat rim", on))]
                return T("Shell"), [mm("Thickness"), mm("Detail filter")], self.drag is None, checks
            if k == "hole":
                return T("Hole"), [mm("Radius"), ("Segments", v["Segments"], "", 0, r("Segments")), mm("Support")], \
                    self.drag is None
            if k == "bridge":
                return T("Bridge"), [(n, v[n], "", d, r(n)) for n, d in zip(self.BRIDGE, (0, 0, 2, 2))], \
                    self.drag is None
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
            if k == "knife":
                t, L = v["t"], a["L"]
                return "Knife", [("Position", t * 100, " %", 2, r("Position")),
                                 ("Distance", t * L, " mm", 3, r("Distance"))], self.drag is None
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
        if not verts or self.drag or self.grab:
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
        if not verts or self.drag or self.grab or value <= 0:
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
        doc = App.ActiveDocument if self.obj is None else self.obj.Document
        if self.obj is None:                               # yeni nesne: olusturma her zaman geri alinabilir
            doc.openTransaction("SubD")
            self.obj = doc.addObject("Mesh::FeaturePython", "SubD")
            SubDObject(self.obj)
            if App.GuiUp:
                self.obj.ViewObject.Proxy = 0
            self.obj.ViewObject.Visibility = False
            self._write(full)
            doc.commitTransaction()
            self._snapshot()
        elif commit and self.edit_start is not None and self._same_as_start():
            if self.dirty:                                 # ara kayit yazilmis ama kafes basa dondu: sessizce geri yaz
                self._write(full)
        elif commit and self.edit_start is not None:
            c0, lv0, r0, m0 = self.edit_start
            cage_to_obj(c0, self.obj)                      # sessizce oturum basina
            set_surface(self.obj, lv0, r0)
            m0 = m0.copy()
            m0.Placement = self.obj.Placement
            self.obj.Mesh = m0
            if hasattr(self.obj.Proxy, "_key"):            # yuzey onbellegi son ara kaydi tutuyor: yeniden uretilsin
                self.obj.Proxy._key = None
            doc.openTransaction(T("SubD edit %s") % self.obj.Label)
            self._write(full)
            doc.commitTransaction()
            self._snapshot()
        else:
            self._write(full)
        self.dirty = False
        return self.obj

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
                and c0.V.shape == c.V.shape and np.array_equal(c0.V, c.V))

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
        self.undo, self.splits, self.op, self.last_op = [], [], None, None
        self.sel, self.sel_v, self.sel_e = set(), set(), set()
        self.dirty = False
        self.set_placement(o.Placement)
        self._snapshot()
        self.info = T("Document undo: cage reloaded")
        self.refresh()

    def _autosave(self):
        if self.mode != "object" and self.dirty and self.drag is None and self.grab is None and self.cage.F:
            self.save(full=False)

    def do_undo(self):
        if not self.undo:
            self.info = T("Nothing to undo in this edit (O: object mode, Ctrl+Z undoes whole edits)")
            self.refresh()
            return
        if self.undo:
            self.cage, self.sel, self.sel_v, self.sel_e = self.undo.pop()
            self.grab = None
            self.op = self.last_op = None
            self.dirty = True                              # kafes kaydedilmis olabilir (or. Knife ayirma): yeniden yaz
            self._undo_splits()
            self.refresh()

    def _undo_splits(self):
        """Geri alma, "Split into two objects" kesiminin oncesine indiyse o kesimin olusturdugu nesneyi de sil
        (yoksa ayrilan parca hem kafeste hem yeni nesnede kalirdi)."""
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

    def clear_selection(self):
        self.sel.clear()
        self.sel_v.clear()
        self.sel_e.clear()
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
        self.refresh()

    def toggle_ghost(self):
        """Hayalet modu: yuzeyler yari saydam, arkadaki kenar ve noktalar gorunur ve secilebilir."""
        self.ghost = not self.ghost
        for m, t in ((self.smooth_mat, 0.7), (self.cage_mat, 0.75)):
            m.transparency = t if self.ghost else 0.0
        self.refresh()

    def toggle_zebra(self):
        """Zebra seritleri: puruzsuz yuzeye yansiyan siyah / beyaz seritler. Kirik serit = kirik yuzey (G0),
        seritte kose = G1, yumusak serit = G2. Subpatch kapaliysa acilir."""
        self.zebra = not self.zebra
        self.zebra_switch.whichChild = 0 if self.zebra else -1
        if self.zebra and not self.smooth:
            self.toggle_smooth()
        self._sync_surface_widgets()
        _zebra_button(self.zebra)
        self.refresh()

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
        """Her nokta icin X'e gore aynadaki esini dondurur (duzlemdekiler kendisi)."""
        V = self.cage.V
        if not len(V):
            return {}
        tol = 1e-4 * max(np.ptp(V, 0).max(), 1.0)
        key = lambda p: tuple(np.round(p / tol).astype(np.int64))
        index = {key(p): i for i, p in enumerate(V)}
        return {i: index.get(key(p * [-1, 1, 1])) for i, p in enumerate(V)}

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
        op = self.new_op("extrude", dict(Distance=0.0, X=0.0, Y=0.0, Z=0.0, Segments=float(self.extrude_segs)),
                         axes=axes, dw=np.zeros(3))
        st = self._extrude_build(op["before"], self.extrude_segs, axes)
        if st is None:
            self.op = None
            self.do_undo()
            self.info = T("Extrude: select open (boundary) edges")
            self.refresh()
            return
        st.update(op=op, lazy=lazy)
        self._extrude_anchor(st)
        self.grab = st
        self.gizmo.place(None)
        self.extrude_set(st, op["values"])

    def _extrude_build(self, before, segs, axes):
        """Durumu extrude oncesine dondurur ve extrude topolojisini kurar (yeni katlar tabanda); noktalari
        extrude_set yerlestirir. Acik kenar yoksa None."""
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
            # normaller birbirini goturuyorsa (or. butun yuzler) nokta basina Smooth Shift yonleri
            nv = cage.bevel_dirs(region)[0] if m_axis is None else None
            cage.extrude(region, segs, layers)
            verts = sorted(layers)
            if m_axis is None:
                dirs = nv
            else:
                dirs = np.tile(m_axis, (len(verts), 1))
                if self.sym:                               # aynadaki taraf aynalanmis normalle
                    dirs[[i for i, v in enumerate(verts) if v not in orig]] *= MIRROR
        else:
            edges = set(self.sel_e)
            if self.sym:
                alive = cage.edges()
                edges |= {me for me in (frozenset(M.get(v) for v in e) for e in self.sel_e)
                          if None not in me and me in alive}
            res = cage.extrude_edges(edges, segs, layers)
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
        side = np.array([1 if v in orig else -1 for v in verts])
        plane = [i for i, v in enumerate(verts) if M.get(v) == v]
        if plane:
            dirs = dirs.copy()
            dirs[plane, 0] = 0.0
        base = cage.V[verts].copy()
        return dict(kind="extrude", before=before, segs=segs, axes=axes, verts=verts, base=base,
                    idx=np.array([layers[v][1:] for v in verts]).T, dirs=dirs, side=side, plane=plane,
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
        st = self._extrude_build(g["before"], segs, g["axes"])
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

    def _reapply(self, op):
        c, a, b, e = op["before"]
        self.cage, self.sel, self.sel_v, self.sel_e = c.copy(), set(a), set(b), set(e)
        self._apply_op(op)
        self.dirty = True
        self.refresh()

    def bevel(self, inset_only=False, values=None, lazy=False, push=True):
        """LightWave Bevel / Smooth Shift: secili yuzleri bolge olarak extrude eder; yeni noktalar yuz
        normalleri boyunca kaydirilir (Shift) ve sinirdan iceri alinir (Inset). values=(inset, shift) mm
        verilirse dogrudan uygular; verilmezse fareyle: saga = inset, yukari = shift (I: yalniz inset)."""
        if not self.sel or self.mode != "face" or self.drag or self.grab:
            return
        if push:                                           # push=False: son islem yeniden uygulaniyor
            self.drop_tool()
            self.push()
            self.new_op("bevel", dict(Inset=0.0, Shift=0.0), inset_only=inset_only)
        region = self.mirror_faces(self.sel, self.mirror_map()) if self.sym else set(self.sel)
        Nv, Ov = self.cage.bevel_dirs(region)
        new = self.cage.extrude(region)
        st = dict(kind="bevel", new=new, base=self.cage.V[new].copy(), n=Nv, o=Ov, inset_only=inset_only)
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
            o.ViewObject.Proxy = 0
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

    def line_plane(self, p0, p1):
        """Ekrandaki cizgiden kesim duzlemi (iki uc ve bakis yonunden gecen): (n, d, gosterilecek p1) ya da None.
        Ctrl: cizgi acisi 15 derece adim; ekranda gosterilen cizgi de bu yuvarlanmis aciyla cizilir."""
        p0, p1 = np.array(p0, float), np.array(p1, float)
        if _ctrl_held():
            v = p1 - p0
            ang = math.radians(_snap(math.degrees(math.atan2(v[1], v[0])), 15.0))
            p1 = p0 + np.linalg.norm(v) * np.array([math.cos(ang), math.sin(ang)])
        P0 = np.array(tuple(self.view.getPoint(int(round(p0[0])), int(round(p0[1])))))   # numpy sayisi kabul etmiyor
        P1 = np.array(tuple(self.view.getPoint(int(round(p1[0])), int(round(p1[1])))))
        if self.view.getCameraType() == "Perspective":
            eye = np.array(self.view.getCameraNode().position.getValue().getValue())
            n = np.cross(P0 - eye, P1 - eye)
        else:
            n = np.cross(P1 - P0, np.array(tuple(self.view.getViewDirection())))
        ln = np.linalg.norm(n)
        if ln < 1e-12:
            return None
        n /= ln
        for e in np.eye(3):                                # dunya eksenine cok yakinsa tam eksene oturt
            if abs(n.dot(e)) > math.cos(math.radians(2.0)):
                n = e * np.sign(n.dot(e))
        nl = self.w2l_dir(n)                               # kesim kafesin yerel koordinatinda yapilir
        return nl, float(nl.dot(self.w2l(P0))), p1

    def line_live(self, b, pos):
        """Knife cizgisi surulurken (LightWave gibi canli): durumu kesim oncesine dondur, yeni duzlemle kes."""
        self._pending = False
        if b.get("cancelled"):
            return
        if "before" not in b:                              # ilk hareket: geri alma noktasi ve kesilecek yuzler
            verts = set(self.selected_verts()) if self.mode == "face" and self.sel else None
            if verts is not None and self.sym:
                M = self.mirror_map()
                verts |= {M[v] for v in list(verts) if M.get(v) is not None}
            self.push()
            b["before"], b["verts"] = self.undo[-1], verts
        c = b["before"][0]
        self.cage = c.copy()
        self.sel, self.sel_v, self.sel_e = set(), set(), set()
        pl = self.line_plane(b["p0"], pos)
        b["plane"] = pl
        if pl is not None:
            n, d, p1 = pl
            b["side"] = self._plane_cut(n, d, b["verts"], _knife_opts["split"] or _knife_opts["objects"])
            self._line_show(b["p0"], p1)
            self.info = T("Knife  N (%.3f, %.3f, %.3f)   Offset %.3f mm") % (n[0], n[1], n[2], d)
        self.refresh(fast=True)

    def line_finish(self, b):
        """Birakildi: canli kesim kalir ve 'son islem' olur (panelde Normal / Konum); hicbir sey kesilmediyse geri al."""
        self._line_hide()
        if b.get("cancelled"):
            return
        if "before" not in b:
            self.info = ""
            self.refresh()
            return
        c = b["before"][0]
        pl = b.get("plane")
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
        if b is not None:
            b["cancelled"] = True                          # bekleyen canli adimlar bir sey yapmasin
            if "before" in b:
                self.cage, self.sel, self.sel_v, self.sel_e = self.undo.pop()
        self.info = ""
        self.refresh()

    def _line_show(self, p0, p1):
        vp = self.gv.viewport()
        if self.line_overlay is None:
            self.line_overlay = _LineOverlay(vp)
        dpr = vp.devicePixelRatioF()
        H = self.view.getSize()[1]
        q = lambda p: QtCore.QPointF(p[0] / dpr, (H - p[1]) / dpr)
        self.line_overlay.set_line(q(p0), q(p1))

    def _line_hide(self):
        if self.line_overlay is not None:
            self.line_overlay.set_line(None, None)

    # --- kalici araclar (LightWave gibi): dugme sari yanar, kapatilana kadar her tiklama araci uygular ---
    # "knife": kenara tikla = loop cut + kaydirma; "line": cizgi cek = duzlemle kes; "loop": kenara tikla = halka sec.
    # Kapatma: ayni dugme / tus (K, J), Esc, sag tik; baska bir arac (extrude, bevel, tasima...) da kapatir.

    def toggle_tool(self, name):
        if self.drag or self.grab or self.box is not None or self.mode == "object":
            return
        self.armed = None if self.armed == name else name
        self.info = {None: "",
                     "knife": T("Knife: click edges to cut (K / Esc / right click: off)"),
                     "line": T("Knife: drag cut lines (Ctrl: 15° steps; J / Esc / right click: off)"),
                     "loop": T("Loop: click edges to select loops (Esc / right click: off)")}[self.armed]
        self.refresh()

    def drop_tool(self):
        if self.armed:
            self.armed = None
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
        if self.sel and self.mode == "face":
            self.drop_tool()
            self.push()
            self.cage.delete(self.mirror_faces(self.sel, self.mirror_map()) if self.sym else self.sel)
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
        """Ctrl+C / Ctrl+X: secili cokgenleri (keskin kenarlariyla) panoya alir; cut: sonra siler."""
        global _clipboard
        if self.mode == "object" or self.drag or self.grab:
            return
        sel = self.selected_faces()
        if self.sym:
            sel = self.mirror_faces(sel, self.mirror_map())
        if not sel:
            self.info = T("Copy: select polygons first")
            self.refresh()
            return
        sub = self.cage.subset(sel)
        sub.V = self.l2w(sub.V)
        _clipboard = sub
        self.info = (T("Cut %d polygons") if cut else T("Copied %d polygons")) % len(sub.F)
        if cut:
            self.push()
            self.cage.delete(sel)
            self.clear_selection()
        self.refresh()

    def paste_faces(self):
        """Ctrl+V: panodaki cokgenleri bu kafese ayri parca olarak ekler (kopyalandigi yerde) ve secer;
        T / gizmo ile hemen tasinabilir. Object modunda yeni bir SubD nesnesi olur."""
        if _clipboard is None or not _clipboard.F or self.drag or self.grab:
            return
        if self.mode == "object":
            self.object_paste()
            return
        self.push()
        c = _clipboard.copy()
        c.V = self.w2l(c.V)
        new = self.cage.append(c)
        self.sel, self.sel_v, self.sel_e = set(), set(), set()
        if self.mode == "face":
            self.sel = set(new)
        elif self.mode == "point":
            self.sel_v = {v for fi in new for v in self.cage.F[fi]}
        else:
            self.sel_e = {_edge(f[i], f[(i + 1) % len(f)]) for f in (self.cage.F[fi] for fi in new)
                          for i in range(len(f))}
        self.last_op = None
        self.info = T("Pasted %d polygons (T: move)") % len(new)
        self.refresh()

    def paste_as_object(self):
        """Panodaki cokgenlerden yeni bir SubD nesnesi (or. sisenin agzindan kapak) ve hemen kafesine gir."""
        if _clipboard is None or not _clipboard.F or self.drag or self.grab:
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
            self.do_undo()
            self.info = err
        else:
            self.last_op, self.op = op, None
        self.refresh()

    def _bridge(self, v):
        """Bridge'i uygular; olmazsa hata metni dondurur (kafes degismez)."""
        cage = self.cage
        v["Number of Cuts"] = float(min(max(int(round(v["Number of Cuts"])), 0), 64))
        v["Twist"] = float(int(round(v["Twist"])))
        self.bridge_vals = dict(v)
        ef = cage.edge_faces()
        region = set(self.sel) if self.mode == "face" else set()
        if region:                                         # bolge siniri: bolgede tek yuzu olan kenarlar
            edges = {e for e, fs in ef.items() if sum(f in region for f in fs) == 1}
        else:
            edges = set(self.sel_e)
            if any(len(ef.get(e, ())) != 1 for e in edges):
                return T("Bridge: select open (boundary) edges")
        chains = cage.edge_chains(edges)
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
            t = int(v["Twist"]) % n
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
        # yuz yonu: A[0] -> A[1] kenarinda kalan komsu yuz ayni yonde gidiyorsa yeni yuzler ters sirali
        e0 = _edge(A[0], A[1])
        fs = ef.get(e0, ())
        keep = [f for f in fs if f not in region]
        f = cage.F[(keep or fs)[0]]
        same = any(f[i] == A[0] and f[(i + 1) % len(f)] == A[1] for i in range(len(f)))
        flip = same if keep else not same                  # silinen yuzun yerini alan seridin yonu onunla ayni
        new = cage.bridge(A, B, ca, int(v["Number of Cuts"]), max(v["Smoothness"], 0.0), v["Profile Factor"],
                          dA, dB, flip)
        if region:
            cage.delete(region)                            # yeni yuzler sonda: indeksleri silinenler kadar kayar
            new = [fi - len(region) for fi in new]
        self.sel, self.sel_v, self.sel_e = set(new) if self.mode == "face" else set(), set(), set()
        self.info = T("Bridge: %d polygons") % len(new)
        return None

    # --- From face: FreeCAD'de secili Part / PartDesign yuzeylerinden / kenar halkasindan kafes (cage_from_faces) ---
    # Object modunda yeni SubD nesnesi olur (dunya koordinati), kafeste duzenlerken ayni kafese eklenir. Panel:
    # Arc Segments (90 derece basina), Max Length (mm, 0 = otomatik), Exact Edge (1: SubD kenari yuzey sinirindan
    # tam gecer, 0: kafes noktalari sinirda).
    FACE = ("Arc Segments", "Max Length", "Exact Edge")

    def _selected_part_faces(self):
        """FreeCAD secimindeki Part / PartDesign yuzeyleri, dunya koordinatinda (Body yerlesimi dahil). Secili
        kenarlar kapali halka(lar) olusturuyorsa onlardan duz yuz yapilir (ic ice iki halka: halka yuz). Kopya
        alinmaz: ortak kenarlar (isSame) taninsin. Hata metni dondurebilir."""
        import Part
        faces, edges = [], []
        for so in Gui.Selection.getSelectionEx("", 0):
            if is_subd(so.Object):
                continue
            for sub in so.SubElementNames:
                try:
                    sh = Part.getShape(so.Object, sub, needSubElement=True)
                except Exception:
                    continue
                if not sh.isNull() and sh.ShapeType == "Face":
                    faces.append(sh)
                elif not sh.isNull() and sh.ShapeType == "Edge":
                    edges.append(sh)
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
        vals = {"Arc Segments": 4.0, "Max Length": 0.0, "Exact Edge": 1.0}
        test = self._face_cage(faces, dict(vals))
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
                o.ViewObject.Proxy = 0
            cage_to_obj(test, o)
            set_surface(o, self.levels, self.rules)
            o.recompute()
            doc.commitTransaction()
            self.edit_mode = "face"
            self.edit_object(o)
            self.push()
            op = self.new_op("fromface", vals, faces=faces, append=False)
        else:
            self.drop_tool()
            self.push()
            op = self.new_op("fromface", vals, faces=faces, append=True)
        self._apply_op(op)
        self.last_op, self.op = op, None
        self.refresh()

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
            self.do_undo()
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

    def hole(self):
        if self.drag or self.grab or self.mode != "face":
            return
        self.drop_tool()
        self.push()
        op = self.new_op("hole", dict(Radius=0.0, Segments=0.0, Support=-1.0))   # 0 / -1: otomatik
        err = self._hole(op["values"])
        if err:
            self.op = None
            self.do_undo()
            self.info = err
        else:
            self.last_op, self.op = op, None
        self.refresh()

    def _hole(self, v):
        res = self.cage.hole(set(self.sel), v["Radius"], int(round(v["Segments"])), v["Support"])
        if isinstance(res, str):
            return res
        new, v["Radius"], segs, v["Support"] = res
        v["Segments"] = float(segs)
        self.sel, self.sel_v, self.sel_e = set(new), set(), set()
        self.info = T("Hole: radius %.2f mm, %d segments") % (v["Radius"], segs)
        return None

    def object_paste(self):
        doc = App.ActiveDocument
        doc.openTransaction("SubD paste")
        o = doc.addObject("Mesh::FeaturePython", "SubD")
        SubDObject(o)
        if App.GuiUp:
            o.ViewObject.Proxy = 0
        cage_to_obj(_clipboard, o)                         # dunya koordinati: Placement birim
        set_surface(o, self.levels, self.rules)
        o.recompute()
        doc.commitTransaction()
        self.edit_mode = "face"
        self.edit_object(o)
        self.sel = set(range(len(self.cage.F)))
        self.info = T("Pasted %d polygons as a new object") % len(self.cage.F)
        self.refresh()

    # --- tasima (grab) ---
    def start_grab(self, verts, axis=None, pushed=False, lazy=False, op=None):
        if not verts or self.drag or self.grab:
            return
        self.drop_tool()
        if not pushed:
            self.push()
        if op is None:
            op = self.new_op("move", dict(X=0.0, Y=0.0, Z=0.0), verts=list(verts))
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
            self.refresh()
        else:
            self.op = None
            self.do_undo()

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
        verts = self.selected_verts()
        if not verts:
            return
        self.push()
        self.info = ""
        c = self.cage.V[verts].mean(0)                     # yerel merkez (islem bununla yeniden yapilir)
        c0 = App.Vector(*self.l2w(c))                      # gizmo hesaplari dunya koordinatinda
        kind, i = key
        if kind == "rotate":
            self.new_op("rotate", dict(Angle=0.0), verts=verts, c=c, axis=tuple(AXES[i]), letter="XYZ"[i])
        elif kind in ("stretch", "size"):
            self.new_op("scale", dict(X=1.0, Y=1.0, Z=1.0), verts=verts, c=c)
        else:
            self.new_op("move", dict(X=0.0, Y=0.0, Z=0.0), verts=verts)
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
        self.last_op, self.op = self.op, None
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
        if not down:
            return
        if self.grab:
            cb.setHandled()
            _later(lambda ok=left: self.end_grab(ok))
            return
        if not left:
            if self.armed and ev.getButton() == coin.SoMouseButtonEvent.BUTTON2:   # sag tik: kalici araci kapat
                self.drop_tool()
                _later(self.refresh)
            return
        # kalici arac (sari dugme) acik kaldikca her tiklama onu uygular; secim / kutu secimi yapilmaz
        if self.armed == "line":                           # knife cizgisi baslar; birakinca kesilir
            cb.setHandled()
            self.box = dict(p0=pos, moved=False, kind="line", shift=False, ctrl=False)
            return
        if self.armed in ("knife", "loop"):
            cb.setHandled()
            fn = self.knife if self.armed == "knife" else self.loop_select
            _later(lambda: fn(pos))
            return
        key = self.gizmo.pick(pos)
        if key is not None:
            cb.setHandled()                                # secime ve gezinmeye gitmesin
            self.drag_begin(key, pos)
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
        if self.drag is not None:
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
            if b["moved"] or abs(pos[0] - b["p0"][0]) + abs(pos[1] - b["p0"][1]) > 4:
                b["moved"] = True
                if b.get("kind") == "line":
                    b["pos"] = pos                         # canli kesim: hareketler birlestirilir
                    self._schedule(lambda: self.line_live(b, b["pos"]))
                else:
                    _later(lambda: self._rubber_show(b["p0"], pos))
        else:
            key = self.gizmo.pick(pos)
            _later(lambda: self.gizmo.set_hover(key))

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
            elif k in ("x", "y", "z") and self.grab.get("kind") == "extrude":
                self.extrude_axis(k.upper())               # yalniz bu eksen; tekrar = normal yon
                self.refresh()
            return
        if self.mode == "object":                          # Object modu: yalniz nesne komutlari
            if ctrl and k == "z":
                doc = App.ActiveDocument
                if doc.UndoCount:
                    doc.undo()
                    doc.recompute()                       # katilar geri alinan kafese gore yenilensin
                self.doc_undone()
            elif ctrl and k == "v":                       # panodaki cokgenler -> yeni SubD nesnesi
                self.paste_faces()
            elif ctrl and k in ("c", "x"):                # nesne modunda kopyalanacak cokgen yok
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
        if ctrl:                                           # LightWave: ^G nokta, ^H cokgen
            {"z": self.do_undo, "g": lambda: self.set_mode("point"), "h": lambda: self.set_mode("face"),
             "c": self.copy_faces, "x": lambda: self.copy_faces(cut=True), "v": self.paste_faces}[k]()
        elif k == "n":
            self.toggle_options()
        elif k == "escape" and self.armed:
            self.drop_tool()
            self.refresh()
        elif k in ("tab", "s"):
            self.toggle_smooth()
        elif k == "e":
            self.extrude()
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

    def close(self, doc_gone=False):
        """Kafesi kaydet ve modelleyiciyi kaldir. Her adim ayri korunur: belge / 3B gorunum once kapandiysa
        widget'lari ve sahnesi zaten silinmistir; bir adimin hatasi kalanlari (or. arac cubugu) engellemesin,
        oturum her durumda kapansin (yoksa modelleyici bir daha acilamiyordu). doc_gone: belge kapaniyor,
        kaydetme / gosterme yok."""
        global _session
        if _session is self:
            _session = None
        _quiet(lambda: _zebra_button(False))
        _quiet(self.timer.stop)
        _quiet(self.autosave.stop)
        if not doc_gone and self.mode != "object":         # Object modunda yazilacak kafes yok
            try:
                if self.cage.F:
                    self.save(commit=True)
                self._show_edited()
            except Exception as e:
                App.Console.PrintError(T("SubD could not be saved: %s\n") % e)
        _quiet(lambda: self.view.removeEventCallbackPivy(coin.SoMouseButtonEvent.getClassTypeId(), self._cb_btn))
        _quiet(lambda: self.view.removeEventCallbackPivy(coin.SoLocation2Event.getClassTypeId(), self._cb_mov))
        for w in self.key_targets:
            _quiet(lambda w=w: w.removeEventFilter(self.key_filter))
        _quiet(lambda: self.parent.removeChild(self.root))
        if self.overlay is not None:
            _quiet(self.overlay.close)
            self.overlay = None
        for w in (self.rubber, self.line_overlay):
            if w is not None:
                _quiet(w.deleteLater)
        self.rubber = self.line_overlay = None
        _quiet(self.dock.close)
        _quiet(self.dock.deleteLater)


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


class _DocWatcher(object):
    """Duzenlenen belge kapaninca modelleyiciyi de kapat: 3B gorunum silinince sol menu / alt cubuk onunla
    gidiyor ama ana penceredeki Create / Modify arac cubugu ve oturum ortada kaliyordu."""

    def slotDeletedDocument(self, doc):
        s = _session
        if s is not None and s.doc_name == doc.Name:
            s.close(doc_gone=True)

    def slotUndoDocument(self, doc):
        s = _session
        if s is not None and s.doc_name == doc.Name:
            _later(lambda: _quiet(s.doc_undone))

    def slotRedoDocument(self, doc):
        self.slotUndoDocument(doc)

    def slotStartSaveDocument(self, doc, filename):
        # Ctrl+S: kafesin son hali dosyaya girsin (yoksa en son otomatik kayit / F ile yazilan hal kaydedilirdi)
        s = _session
        if (s is not None and s.doc_name == doc.Name and s.mode != "object" and s.cage.F
                and s.drag is None and s.grab is None and (s.dirty or s.obj is None)):
            try:
                s.save(full=False)
            except Exception as e:
                App.Console.PrintError(T("SubD could not be saved: %s\n") % e)


App.addDocumentObserver(_DocWatcher())


_KEYS = {
    QtCore.Qt.Key_Tab: "tab", QtCore.Qt.Key_S: "s", QtCore.Qt.Key_E: "e", QtCore.Qt.Key_K: "k",
    QtCore.Qt.Key_T: "t", QtCore.Qt.Key_G: "g", QtCore.Qt.Key_A: "a", QtCore.Qt.Key_F: "f",
    QtCore.Qt.Key_Z: "z", QtCore.Qt.Key_Delete: "delete", QtCore.Qt.Key_Backspace: "backspace",
    QtCore.Qt.Key_Escape: "escape", QtCore.Qt.Key_Return: "return", QtCore.Qt.Key_Enter: "enter",
    QtCore.Qt.Key_1: "1", QtCore.Qt.Key_2: "2", QtCore.Qt.Key_3: "3", QtCore.Qt.Key_Y: "y", QtCore.Qt.Key_R: "r",
    QtCore.Qt.Key_C: "c", QtCore.Qt.Key_B: "b", QtCore.Qt.Key_I: "i",
    QtCore.Qt.Key_H: "h", QtCore.Qt.Key_N: "n", QtCore.Qt.Key_L: "l", QtCore.Qt.Key_X: "x",
    QtCore.Qt.Key_J: "j", QtCore.Qt.Key_D: "d", QtCore.Qt.Key_O: "o", QtCore.Qt.Key_V: "v", QtCore.Qt.Key_P: "p",
}
# Ctrl+Z geri al, Ctrl+G nokta modu, Ctrl+H cokgen modu (LightWave), Ctrl+C / X / V cokgen kopyala / kes / yapistir
_CTRL_KEYS = ("z", "g", "h", "c", "x", "v")
_CTRL_ONLY = ("z", "v")        # tek basina bizim degil (FreeCAD'e gider)


class _KeyFilter(QtCore.QObject):
    def __init__(self, session):
        super().__init__()
        self.session = session

    def eventFilter(self, obj, ev):
        t = ev.type()
        g = self.session.grab
        if t == QtCore.QEvent.Wheel and g is not None and g.get("kind") in ("ebevel", "extrude"):
            # kenar bevel / extrude surerken tekerlek yakinlastirmaz, segment sayisini degistirir
            d = ev.angleDelta().y()
            if d:
                (self.session.ebevel_wheel if g["kind"] == "ebevel" else self.session.extrude_wheel)(1 if d > 0 else -1)
            return True
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
        k = _KEYS.get(ev.key())
        ctrl = bool(ev.modifiers() & QtCore.Qt.ControlModifier)
        extr = g is not None and g.get("kind") == "extrude"   # extrude surerken Z = eksen (Ctrl+Z degil)
        if k is None or (ctrl and k not in _CTRL_KEYS) or (not ctrl and k in _CTRL_ONLY and not (extr and k == "z")):
            return False
        if t == QtCore.QEvent.ShortcutOverride:
            ev.accept()  # tusu FreeCAD kisayollarindan once biz alalim
            return True
        if not ev.isAutoRepeat():
            self.session.key(k, ctrl)
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
QLabel#lwInfo { color: #f59e0b; font-weight: bold; }
"""


MENU_CSS = """
QMenu { background: #3a3a3a; color: #ececec; border: 1px solid #1e1e1e; padding: 2px; }
QMenu::item { padding: 3px 18px 3px 10px; }
QMenu::item:selected { background: #3f6fb0; }
"""


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
    b.setMinimumWidth(max(width or 0, fm.horizontalAdvance(text + key) + 24))   # yazi hicbir zaman kesilmesin
    return b


class _LineOverlay(QtWidgets.QWidget):
    """Knife cizgisi: viewport'un ustunde, fare olaylarini gecirir."""

    def __init__(self, parent):
        super().__init__(parent)
        self.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents)
        self.a = self.b = None
        self.resize(parent.size())
        self.show()

    def set_line(self, a, b):
        self.a, self.b = a, b
        self.resize(self.parent().size())
        self.raise_()
        self.update()

    def paintEvent(self, ev):
        if self.a is None:
            return
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        pen = QtGui.QPen(QtGui.QColor(245, 158, 11), 2, QtCore.Qt.DashLine)
        p.setPen(pen)
        p.drawLine(self.a, self.b)
        p.setBrush(QtGui.QColor(245, 158, 11))
        for pt in (self.a, self.b):
            p.drawEllipse(pt, 3, 3)
        p.end()


class _ClickFocus(QtCore.QObject):
    """Viewport icindeki kutuya tiklaninca klavye odagi kutuya gecmiyor (3B gorunum geri aliyor): yazilan
    rakamlar gorunume gidip 1/2/3 mod degistiriyordu. Tiklamada odagi acikca kutuya ver, degeri sec."""

    def eventFilter(self, obj, ev):
        if ev.type() == QtCore.QEvent.Wheel and isinstance(obj, QtWidgets.QAbstractSpinBox):
            if not obj.hasFocus():                         # odaksiz kutu tekerlegi yok sayiyor: once odakla
                obj.setFocus(QtCore.Qt.MouseFocusReason)
            return False
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
            QCheckBox { color: #e0e0e0; spacing: 3px; }""")
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
        lay.addWidget(self.check_bar)
        self.check_bar.hide()
        self.grid = QtWidgets.QGridLayout()
        self.grid.setHorizontalSpacing(6)
        self.grid.setVerticalSpacing(2)
        lay.addLayout(self.grid)
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
        key = (title, tuple(f[0] for f in fields), tuple(c[0] for c in checks), label)
        changed = key != self.key
        self.check_fns = [c[2] for c in checks]
        if (label, tuple(c[0] for c in checks)) != self.check_key:
            self._build_checks(label, [c[0] for c in checks])
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
                box.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
                box.setAlignment(QtCore.Qt.AlignRight)
                # fare tekerlegi alanin ustunde: deger bir adim degisir ve islem hemen yeniden yapilir (segment,
                # kesim sayisi, mm); adim: tam sayi / mm / derece 1, carpan (4 ondalik) 0.05, oteki 0.1
                box.setSingleStep(1.0 if dec == 0 or suffix.strip() in ("mm", "°") else 0.05 if dec >= 4 else 0.1)
                box.editingFinished.connect(lambda b=box, i=r: self._commit(i, b))
                box.valueChanged.connect(lambda val, b=box, i=r: self._commit(i, b))
                box.lineEdit().installEventFilter(self._focus_filter)
                box.installEventFilter(self._focus_filter)          # tekerlek: odaklanip adim degissin
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

    def _build_checks(self, label, names):
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
        for i, n in enumerate(names):
            cb = QtWidgets.QCheckBox(T(n))
            cb.setFocusPolicy(QtCore.Qt.NoFocus)           # kisayollar gorunumde kalsin
            cb.toggled.connect(lambda on, i=i: self._check(i, on))
            if label:
                self.check_lay.addWidget(cb, 0, i + 1)
            else:
                self.check_lay.addWidget(cb, i, 0)
            cb.show()                                      # gorunur panele eklenen cocuk kendiliginden gorunmez;
            self.checks.append(cb)                         # gorunmezse panel yuksekligine girmez (panel ezilirdi)
        self.check_key = (label, tuple(names))

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
        _later(lambda: self.s.gv.setFocus())     # kisayollar gorunumde calismaya devam etsin
        _later(lambda: fn(val))


_last_tab = ["create"]   # son secilen sekme; modelleyici yeniden acilinca (or. dil degisince) ayni sekme

TAB_CSS = """
QToolBar { spacing: 3px; }
QToolButton { background: #5a5a5a; color: #ececec; border: 1px solid #2b2b2b; border-radius: 2px;
              padding: 3px 18px; min-width: 47px; }
QToolButton:hover { background: #6b6b6b; }
QToolButton:checked { background: #3f6fb0; border-color: #24426b; }
"""


class _TabBar(object):
    """FreeCAD arac cubugunda Create / Modify sekme dugmeleri (LightWave'in ust sekmeleri gibi): sol arac
    sutununda hangi menunun gorunecegini secer. Modelleyici kapaninca arac cubugu kaldirilir."""
    TABS = (("create", "Create"), ("modify", "Modify"))

    def __init__(self, overlay):
        self.mw = Gui.getMainWindow()
        for old in self.mw.findChildren(QtWidgets.QToolBar, "SubDModelerTabs"):   # duzgun kapanmamis oturumdan
            _quiet(lambda old=old: self.mw.removeToolBar(old))
            _quiet(old.deleteLater)
        self.tb = QtWidgets.QToolBar(T("SubD Modeler tabs"), self.mw)
        self.tb.setObjectName("SubDModelerTabs")
        self.tb.setStyleSheet(TAB_CSS)
        self.buttons = {}
        for key, text in self.TABS:
            b = QtWidgets.QToolButton()
            b.setText(T(text))
            b.setCheckable(True)
            b.setFocusPolicy(QtCore.Qt.NoFocus)            # kisayollar gorunumde kalsin
            b.clicked.connect(lambda checked=False, k=key: (overlay.set_tab(k), _later(overlay.s.gv.setFocus)))
            self.tb.addWidget(b)
            self.buttons[key] = b
        self.mw.addToolBar(self.tb)
        self.tb.show()

    def sync(self, tab):
        for k, b in self.buttons.items():
            b.setChecked(k == tab)

    def close(self):
        _quiet(lambda: self.mw.removeToolBar(self.tb))
        _quiet(self.tb.deleteLater)


class _Overlay(object):
    def __init__(self, session):
        self.s = s = session
        self.vp = s.gv.viewport()
        self.toggles = {}    # ad -> (dugme, durum fonksiyonu)
        self.gated = []      # (dugme, secim gerekir mi, yalniz hangi modda: None / "edit" / "object")
        self.tool_buttons = []  # kalici arac dugmeleri (Knife, Select loop): aciksa sari, yazi koyu

        # Arac sutunu bolumlerden olusur (baslik + dugmeler). Bolumler LightWave gibi sekmelere ayrilir: FreeCAD
        # arac cubugundaki Create / Modify dugmeleri hangisinin gorunecegini secer. Secili sekme gorunume
        # sigmazsa relayout bolumlerini yan yana sutunlara akitir.
        self.left = QtWidgets.QFrame(self.vp)
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
            b.clicked.connect(lambda checked=False: (fn(), self.sync()))
            if tool:                                       # kalici arac: acikken sari (mavi = acik modlar)
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
            m = QtWidgets.QMenu(b)
            m.setStyleSheet(MENU_CSS)
            for label, fn in items:
                m.addAction(T(label)).triggered.connect(
                    lambda checked=False, fn=fn: (fn(), self.sync(), _later(s.gv.setFocus)))
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
        btn("From face", "", s.from_face)              # FreeCAD'de secili Part / PartDesign yuzeyinden
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
        btn("Bridge", "", s.bridge, sel=True, mode="edit")
        btn("Hole", "", s.hole, sel=True, mode="edit")
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
        btn("Crease", "C", s.crease_selected, sel=True, mode="edit")
        btn("Shell", "", s.shell)                     # Object: secili nesne; Edit: secili cokgenler (yoksa hepsi)
        # kafeste: secili cokgenler (secim yoksa hepsi) ayni kafese aynalanir; Object modunda: secili nesnelerin kopyasi
        menu_btn("Mirror", [("Mirror X", lambda: s.mirror_polygons(0)), ("Mirror Y", lambda: s.mirror_polygons(1)),
                            ("Mirror Z", lambda: s.mirror_polygons(2))])

        self.bar = QtWidgets.QFrame(self.vp)
        self.bar.setObjectName("lwBar")
        self.bar.setStyleSheet(OVERLAY_CSS)
        row = QtWidgets.QHBoxLayout(self.bar)
        row.setContentsMargins(6, 3, 6, 3)
        row.setSpacing(3)
        self.sel_label = QtWidgets.QLabel()
        self.sel_label.setMinimumWidth(60)
        row.addWidget(self.sel_label)

        def tog(text, key, fn, state):
            b = _lw_button(text, key)
            b.setCheckable(True)
            b.clicked.connect(lambda checked=False: (fn(), self.sync()))
            self.toggles[text] = (b, state)
            row.addWidget(b)

        tog("Object", "O", s.toggle_object_mode, lambda: s.mode == "object")
        row.addSpacing(6)
        tog("Point", "^G", lambda: s.set_mode("point"), lambda: s.mode == "point")
        tog("Edge", "2", lambda: s.set_mode("edge"), lambda: s.mode == "edge")
        tog("Polygon", "^H", lambda: s.set_mode("face"), lambda: s.mode == "face")
        row.addSpacing(10)
        tog("Symmetry", "Y", s.toggle_sym, lambda: s.sym)
        tog("Subpatch", "Tab", s.toggle_smooth, lambda: s.smooth)
        tog("Ghost", "X", s.toggle_ghost, lambda: s.ghost)
        tog("Snap", "", s.toggle_snap, lambda: s.snap)
        cc = QtWidgets.QLabel("Catmull-Clark")
        cc.setStyleSheet("padding-left: 10px;")
        row.addWidget(cc)
        row.addSpacing(14)
        self.step_label = QtWidgets.QLabel()
        self.step_label.setMinimumWidth(
            self.step_label.fontMetrics().horizontalAdvance(T("Step: %g mm / %g°") % (10.5, 15.5)) + 8)
        row.addWidget(self.step_label)
        row.addSpacing(10)
        self.info_label = QtWidgets.QLabel()
        self.info_label.setObjectName("lwInfo")
        row.addWidget(self.info_label, 1)
        tog("Options", "n", s.toggle_options, lambda: s.dock.isVisible())
        lang = _lw_button(subd_i18n.language().upper(), "")    # arayuz dili: EN <-> TR
        lang.setToolTip(T("Switch the interface language (reopens the modeler)"))
        lang.clicked.connect(lambda checked=False: _later(s.switch_language))
        row.addWidget(lang)

        self.panel = _ToolPanel(self.vp, s)
        self.tabs = _TabBar(self)
        self.tabs.sync(self.tab)

        for w in (self.left, self.bar):
            w.show()
            w.raise_()
        self.relayout()

    def relayout(self):
        w, h = self.vp.width(), self.vp.height()
        self.bar.setFixedWidth(w)
        self.bar.adjustSize()
        self.bar.move(0, h - self.bar.height())
        self._flow_columns(h - self.bar.height() - 6 - 6 - 8)   # ust bosluk, alt bosluk, cerceve kenarlari
        self.left.move(6, 6)
        # islem paneli: sol menunun sagi, alt cubugun hemen ustu
        self.panel.move(self.left.x() + self.left.width() + 8, h - self.bar.height() - self.panel.height() - 8)

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
        nobj = len(s.selected_objects()) if is_obj else 0
        has = nobj > 0 if is_obj else bool(s.selected_verts())
        for b, sel, mode in self.gated:
            ok = mode is None or (mode == "object") == is_obj
            b.setEnabled(ok and (has or not sel))
        count = nobj if is_obj else {"point": len(s.sel_v), "edge": len(s.sel_e), "face": len(s.sel)}[s.mode]
        self.sel_label.setText("Sel: %d" % count)
        self.step_label.setText(T("Step: %g mm / %g°") % (s.step_box.value(), s.angle_box.value()))
        self.info_label.setText(s.info)
        if self.panel.show_spec(s.panel_spec()):
            self.relayout()

    def close(self):
        self.tabs.close()                                  # ana pencerede: gorunum silinse de kalir, once kaldir
        for w in (self.left, self.bar, self.panel):        # gorunumun cocuklari: gorunumle silinmis olabilir
            _quiet(w.hide)
            _quiet(w.deleteLater)


_session = None


def start(obj=None, object_mode=None):
    """Modelleyiciyi acar. Tek SubD nesnesi secili (ya da verilmis) ise onun kafesini duzenler; birden fazla
    secili ya da hic secili degilken belgede SubD nesneleri varsa Object modunda acar; hic yoksa yeni kutuyla
    Edit modunda baslar."""
    global _session
    if _session:
        _session.close()
    if App.ActiveDocument is None or Gui.ActiveDocument is None:
        App.newDocument("SubD")
    repair_proxies(App.ActiveDocument)
    if object_mode is None:
        if obj is None:
            sel = [o for o in Gui.Selection.getSelection() if is_subd(o)]
            obj = sel[0] if len(sel) == 1 else None
            object_mode = obj is None and any(is_subd(o) for o in App.ActiveDocument.Objects)
        else:
            object_mode = False
    if not object_mode and obj is not None:
        Gui.Selection.clearSelection()
    _session = SubDSession(obj, object_mode=object_mode)
    return _session
