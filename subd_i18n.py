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
    "Move": "Taşı",
    "Rotate": "Döndür",
    "File": "Dosya",
    "Save": "Kaydet",
    "Solid (STEP)": "Katı (STEP)",
    "Close": "Kapat",
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
    "Segments": "Segment",
    "Position": "Konum",
    "Offset": "Konum",
    # --- numeric dock ---
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
    "Edge bevel: drag right, wheel = segments (Ctrl: mm step, Shift: fine)":
        "Kenar bevel: fareyi sağa çek, tekerlek = segment (Ctrl: mm adım, Shift: ince)",
    "Edge bevel  %.3f mm   %d segments": "Kenar bevel  %.3f mm   %d segment",
    "Edge bevel: points on open edges are not supported": "Kenar bevel: açık kenarlı noktalar desteklenmiyor",
    "(B: switch)": "(B: değiştir)",
    "Knife: drag to draw the cut line (Ctrl: 15° steps, Esc: cancel)":
        "Knife: kesim çizgisini sürükleyerek çek (Ctrl: 15° adım, Esc: vazgeç)",
    "Knife  N (%.3f, %.3f, %.3f)   Offset %.3f mm": "Knife  N (%.3f, %.3f, %.3f)   Konum %.3f mm",
    "Knife: the line does not cut the model": "Knife: çizgi modeli kesmiyor",
    "Loop: click an edge of the loop": "Loop: halkadaki bir kenara tıkla",
    "Knife: click the edge to cut (Esc: cancel)": "Knife: kesilecek kenara tıkla (Esc: vazgeç)",
    "Distance %.3f mm": "Mesafe %.3f mm",
    "Angle %s  %.2f°": "Açı %s  %.2f°",
    # --- shortcut help ---
    "1 / 2 / 3: point / edge / polygon mode\n"
    "Y: X symmetry on/off\n"
    "R: rotate gizmo  |  H: stretch / size gizmo\n"
    "C: crease selected edges\n"
    "Click: select  |  Shift+click: add\n"
    "Drag: box select (Shift add, Ctrl remove)\n"
    "L: select loop under cursor  |  X: ghost\n"
    "E: extrude  |  K: knife (loop cut)  |  J: knife (line / plane)\n"
    "B: bevel (right = inset, up = shift)  |  I: inset\n"
    "Edge mode B: edge bevel (right = width)\n"
    "T / G: move  |  Tab / S: SubPatch\n"
    "A: select all  |  Del: delete\n"
    "Ctrl+Z: undo  |  F: save to document\n"
    "Ctrl + gizmo / T: stepped (mm, degrees)\n"
    "Esc / right click: cancel":
        "1 / 2 / 3: nokta / kenar / yüz modu\n"
        "Y: X simetrisi aç/kapa\n"
        "R: gizmo Döndür  |  H: gizmo Stretch / Size\n"
        "C: seçili kenarları keskinleştir (crease)\n"
        "Tık: seç  |  Shift+Tık: ekle\n"
        "Sürükle: kutu seçimi (Shift ekle, Ctrl çıkar)\n"
        "L: imlecin altındaki halkayı seç  |  X: hayalet\n"
        "E: extrude  |  K: knife (loop cut)  |  J: knife (çizgi / düzlem)\n"
        "B: bevel (sağa inset, yukarı shift)  |  I: inset\n"
        "Kenar modunda B: kenar bevel (sağa = genişlik)\n"
        "T / G: taşı  |  Tab / S: SubPatch\n"
        "A: tümünü seç  |  Del: sil\n"
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
