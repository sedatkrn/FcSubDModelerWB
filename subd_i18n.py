# -*- coding: utf-8 -*-
"""SubD Modeler user interface language: English (default) or Turkish.

UI strings are written in English in the code and wrapped with T(); when Turkish is selected T() looks them
up in the TR table below (format placeholders such as %s / %.3f must stay the same).

The language is stored in FreeCAD's parameters:
    User parameter:BaseApp/Preferences/Mod/SubDModeler  ->  Language = "auto" | "en" | "tr"
"auto" (the default) follows FreeCAD's own language: Turkish if FreeCAD runs in Turkish, English otherwise.
"""
import FreeCAD as App

PARAM = "User parameter:BaseApp/Preferences/Mod/SubDModeler"
LANGUAGES = ("en", "tr")

_lang = None


def setting():
    return App.ParamGet(PARAM).GetString("Language", "auto")


def language():
    """Effective language code ("en" or "tr")."""
    global _lang
    if _lang is None:
        s = setting()
        if s in LANGUAGES:
            _lang = s
        else:
            gl = App.ParamGet("User parameter:BaseApp/Preferences/General").GetString("Language", "")
            _lang = "tr" if gl.lower().startswith(("turk", "tr")) else "en"
    return _lang


def set_language(lang):
    """Store the language ("auto", "en" or "tr"); the UI picks it up the next time the modeler opens."""
    global _lang
    App.ParamGet(PARAM).SetString("Language", lang)
    _lang = None


def T(text):
    return TR.get(text, text) if language() == "tr" else text


TR = {
    # --- viewport tool column ---
    "Create": "Oluştur",
    "Box": "Kutu",
    "Select": "Seçim",
    "Select all": "Tümünü seç",
    "Deselect": "Seçimi bırak",
    "Select loop": "Loop seç",
    "Polygon": "Çokgen",
    "Knife (loop)": "Knife (halka)",
    "Knife (line)": "Knife (çizgi)",
    "Delete": "Sil",
    "Modify": "Değiştir",
    "SubD Modeler tabs": "SubD Modeler sekmeleri",
    "Move": "Taşı",
    "Rotate": "Döndür",
    "File": "Dosya",
    "Save": "Kaydet",
    "Solid (STEP)": "Katı (STEP)",
    "Close": "Kapat",
    "Clipboard": "Pano",
    "Copy": "Kopyala",
    "Cut": "Kes",
    "Paste": "Yapıştır",
    "Paste as object": "Yeni nesneye",
    "Copy: select polygons first": "Kopyala: önce çokgen seçin",
    "Copied %d polygons": "%d çokgen kopyalandı",
    "Cut %d polygons": "%d çokgen kesildi",
    "Pasted %d polygons (T: move)": "%d çokgen yapıştırıldı (T: taşı)",
    "Pasted %d polygons as a new object": "%d çokgen yeni nesneye yapıştırıldı",
    # --- object mode ---
    "Object": "Nesne",
    "Duplicate": "Çoğalt",
    "Hole": "Delik",
    "Radius": "Yarıçap",
    "Support": "Destek",
    "Hole: select faces": "Delik: yüz seçin",
    "Hole: the selected faces must form one patch without holes":
        "Delik: seçili yüzler deliksiz tek bir parça olmalı",
    "Hole: the selected faces are not flat": "Delik: seçili yüzler düz değil",
    "Hole: radius too large for the selected faces (max %.2f mm)":
        "Delik: yarıçap seçili yüzlere sığmıyor (en çok %.2f mm)",
    "Hole: radius %.2f mm, %d segments": "Delik: yarıçap %.2f mm, %d segment",
    "Bridge": "Köprü",
    "Number of Cuts": "Kesim sayısı",
    "Twist": "Burgu",
    "Smoothness": "Yumuşaklık",
    "Profile Factor": "Profil çarpanı",
    "Bridge: %d polygons": "Köprü: %d çokgen",
    "Bridge: select open (boundary) edges": "Köprü: açık (sınır) kenarları seçin",
    "Bridge: select two loops (edges or face regions) with the same number of edges":
        "Köprü: kenar sayısı eşit iki halka seçin (kenar ya da yüz bölgesi)",
    "Mirror": "Ayna",
    "Mirror X": "Ayna X",
    "Mirror Y": "Ayna Y",
    "Mirror Z": "Ayna Z",
    "Mirror ": "Ayna ",
    "Mirror %s: %d polygons, %d points merged": "Ayna %s: %d çokgen, %d nokta birleşti",
    "Object mode: %d selected": "Nesne modu: %d seçili",
    "%d objects": "%d nesne",
    "Position ": "Konum ",
    "Rotation ": "Dönüş ",
    # --- bottom bar ---
    "Step: %g mm / %g°": "Adım: %g mm / %g°",
    "Switch the interface language (reopens the modeler)": "Arayüz dilini değiştir (modelleyici yeniden açılır)",
    # --- tool panel titles and fields ---
    "Selection": "Seçim",
    "Center ": "Merkez ",
    "Size ": "Boyut ",
    "Rotate ": "Döndür ",
    "Edge bevel": "Kenar bevel",
    "Knife (plane: N.p = Offset)": "Knife (düzlem: N.p = Konum)",
    "Distance": "Mesafe",
    "Angle": "Açı",
    "Width": "Genişlik",
    "Depth": "Derinlik",
    "Height": "Yükseklik",
    "Segments X": "Segment X",
    "Segments Y": "Segment Y",
    "Segments Z": "Segment Z",
    "Box: size / segments in the panel": "Kutu: ölçü / segment panelde",
    "Segments": "Segment",
    "Position": "Konum",
    "Offset": "Konum",
    # --- Options dock ---
    "SubD Options": "SubD Seçenekler",
    "Measure": "Ölçü",
    "Center": "Merkez",
    "Size": "Boyut",
    "Ctrl step": "Ctrl adım",
    "Sharpness: %g = fully sharp, 1-3 = semi-sharp (soft chamfer)":
        "Keskinlik: %g = tam keskin, 1-3 = yarı keskin (yumuşak pah)",
    "Apply / remove": "Uygula / kaldır",
    "Bevel / Smooth Shift (B, I = inset only)": "Bevel / Smooth Shift (B, I = yalnız inset)",
    "Apply": "Uygula",
    "Bevel the selected faces with these values": "Seçili yüzlere bu değerlerle bevel uygular",
    "New box (width x depth x height)": "Yeni kutu (en x derinlik x yükseklik)",
    "Cage:  %.2f x %.2f x %.2f mm\nSurface:  %.2f x %.2f x %.2f mm":
        "Kafes:  %.2f x %.2f x %.2f mm\nYüzey:  %.2f x %.2f x %.2f mm",
    "%d points, %d faces  |  %d selected": "%d nokta, %d yüz  |  %d seçili",
    "  |  Symmetry X": "  |  Simetri X",
    # --- status / live values ---
    "Solid: %s (%d faces%s)": "Katı: %s (%d yüz%s)",
    ", open shell": ", açık kabuk",
    "SubD solid: %s, %d faces, valid=%s, closed=%s, volume=%.3f mm3\n":
        "SubD katı: %s, %d yüz, geçerli=%s, kapalı=%s, hacim=%.3f mm3\n",
    "SubD could not be saved: %s\n": "SubD kaydedilemedi: %s\n",
    "SubD: repaired %d objects (%s)\n": "SubD: %d nesne onarıldı (%s)\n",
    "Edge bevel: drag right, wheel = segments (Ctrl: mm step, Shift: fine)":
        "Kenar bevel: fareyi sağa çek, tekerlek = segment (Ctrl: mm adım, Shift: ince)",
    "Edge bevel  %.3f mm   %d segments": "Kenar bevel  %.3f mm   %d segment",
    "Edge bevel: points on open edges are not supported": "Kenar bevel: açık kenarlı noktalar desteklenmiyor",
    "(B: switch)": "(B: değiştir)",
    "Knife: drag cut lines (Ctrl: 15° steps; J / Esc / right click: off)":
        "Knife: kesim çizgilerini sürükleyerek çek (Ctrl: 15° adım; J / Esc / sağ tık: kapat)",
    "Knife  N (%.3f, %.3f, %.3f)   Offset %.3f mm": "Knife  N (%.3f, %.3f, %.3f)   Konum %.3f mm",
    "Knife: the line does not cut the model": "Knife: çizgi modeli kesmiyor",
    "Loop: click edges to select loops (Esc / right click: off)":
        "Loop: halka seçmek için kenarlara tıkla (Esc / sağ tık: kapat)",
    "Knife: click edges to cut (K / Esc / right click: off)":
        "Knife: kesmek için kenarlara tıkla (K / Esc / sağ tık: kapat)",
    "Distance %.3f mm": "Mesafe %.3f mm",
    "Axis (none = normal)": "Eksen (boş = normal)",
    "Extrude: select open (boundary) edges": "Extrude: açık (sınır) kenarları seçin",
    "Extrude  %.3f mm   %d segments": "Extrude  %.3f mm   %d segment",
    "Extrude  %s mm   %d segments": "Extrude  %s mm   %d segment",
    "Angle %s  %.2f°": "Açı %s  %.2f°",
    # --- shortcut help ---
    "1 / 2 / 3: point / edge / polygon mode\n"
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
    "A: select all  |  Del: delete\n"
    "Ctrl+C / X / V: copy / cut / paste polygons\n"
    "Ctrl+Z: undo  |  F: save to document\n"
    "Ctrl + gizmo / T: stepped (mm, degrees)\n"
    "Esc / right click: cancel":
        "1 / 2 / 3: nokta / kenar / yüz modu\n"
        "Y: X simetrisi aç/kapa\n"
        "R: gizmo Döndür  |  H: gizmo Stretch / Size\n"
        "C: seçili kenarları keskinleştir (crease)\n"
        "Tık: seç  |  Shift+Tık: ekle\n"
        "Sürükle: kutu seçimi (Shift ekle, Ctrl çıkar)\n"
        "L: seçimden halka (ya da imlecin altındaki)  |  X: hayalet\n"
        "E: extrude  |  K: knife (loop cut)  |  J: knife (çizgi / düzlem)\n"
        "Knife / loop açık kalır (sarı): K / J / Esc / sağ tık kapatır\n"
        "B: bevel (sağa inset, yukarı shift)  |  I: inset\n"
        "Kenar modunda B: kenar bevel (sağa = genişlik)\n"
        "Kenar modunda E: açık kenarları extrude  |  extrude'da tekerlek = segment, X / Y / Z = eksen\n"
        "O: nesne / düzenleme modu  |  D: nesneyi çoğalt\n"
        "T: taşıma gizmosu  |  G: serbest taşı (fare)  |  Tab / S: SubPatch\n"
        "A: tümünü seç  |  Del: sil\n"
        "Ctrl+C / X / V: çokgen kopyala / kes / yapıştır\n"
        "Ctrl+Z: geri al  |  F: belgeye kaydet\n"
        "Ctrl + gizmo / T: adımlı (mm, derece)\n"
        "Esc / sağ tık: işlemi iptal",
    # --- workbench commands ---
    "SubD Modeler": "SubD Modeler",
    "Open the SubD modeler (edits the selected SubD object, or starts with a box)":
        "SubD modelleyiciyi aç (seçili SubD nesnesini düzenler ya da kutuyla başlar)",
    "SubD to solid": "SubD'den katı",
    "Make a parametric B-rep solid (STEP) from the selected SubD object":
        "Seçili SubD nesnesinden parametrik B-rep katı (STEP) oluştur",
    "Select a SubD object first.\n": "Önce bir SubD nesnesi seçin.\n",
    "Interface language": "Arayüz dili",
}
