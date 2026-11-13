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
    "Cylinder": "Silindir",
    "Topo cylinder": "Topo silindir",
    "Cylinder: size / segments in the panel": "Silindir: ölçü / segment panelde",
    "Select": "Seçim",
    "Select all": "Tümünü seç",
    "Deselect": "Seçimi bırak",
    "Select loop": "Loop seç",
    "Polygon": "Çokgen",
    "Knife (loop)": "Knife (halka)",
    "Knife (line)": "Knife (çizgi)",
    "Add points": "Nokta ekle",
    "Connect": "Bağla",
    "Subdivide": "Böl",
    "Split the selected edges into equal parts (Count in the panel); with nothing selected it stays on: "
    "click an edge to add a point there":
        "Seçili kenarları eşit parçalara böler (panelde Count); seçim yoksa açık kalır: tıkladığın kenarın o yerine "
        "nokta ekler",
    "Add an edge between two points of one polygon, or between the middles of two of its edges":
        "Bir çokgenin iki noktası ya da iki kenarının ortası arasına kenar ekler",
    "Split the selected polygons into quads (a polygon with n sides becomes n quads)":
        "Seçili çokgenleri dörtgenlere böler (n kenarlı çokgen n dörtgen olur)",
    "Add points: click an edge to add a point there, near the middle it snaps to the middle (Esc / right click: off)":
        "Nokta ekle: kenara tıkla, oraya nokta eklenir; ortaya yakınsa ortaya oturur (Esc / sağ tık: kapat)",
    "Add points: %d points": "Nokta ekle: %d nokta",
    "Connect: select two points of one polygon that are not neighbours":
        "Bağla: aynı çokgenin komşu olmayan iki noktasını seç",
    "Connect: select two edges of one polygon": "Bağla: aynı çokgenin iki kenarını seç",
    "Connect: select points or edges": "Bağla: nokta ya da kenar seç",
    "Connect: %d edges added": "Bağla: %d kenar eklendi",
    "Subdivide: select polygons": "Böl: çokgen seç",
    "Subdivide: %d polygons -> %d": "Böl: %d çokgen -> %d",
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
    "Copy: select polygons, edges or points first": "Kopyala: önce çokgen, kenar ya da nokta seçin",
    "Copied: %d polygons, %d edges, %d points": "Kopyalandı: %d çokgen, %d kenar, %d nokta",
    "Cut: %d polygons, %d edges, %d points": "Kesildi: %d çokgen, %d kenar, %d nokta",
    "Pasted: %d polygons, %d edges (T: move)": "Yapıştırıldı: %d çokgen, %d kenar (T: taşı)",
    "Pasted as a new object: %d polygons, %d edges": "Yeni nesneye yapıştırıldı: %d çokgen, %d kenar",
    # --- object mode ---
    "Object": "Nesne",
    "Duplicate": "Çoğalt",
    "Make polygon": "Poligon ör",
    "Make polygon: select points or edges": "Poligon ör: nokta ya da kenar seçin",
    "Make polygon: this polygon already exists": "Poligon ör: bu poligon zaten var",
    "Make polygon: an edge already has two polygons": "Poligon ör: bir kenarın zaten iki poligonu var",
    "Make polygon: %d sides": "Poligon ör: %d kenarlı",
    "Make polygon: select 3 or more points, a closed edge loop, or two edges":
        "Poligon ör: 3 ya da daha fazla nokta, kapalı bir kenar halkası ya da iki kenar seçin",
    "From face": "Yüzeyden",
    "Arc Segments": "Yay segmenti",
    "Max Length": "En uzun",
    "Exact Edge": "Tam kenar",
    "From face: curved faces need exactly four edges": "Yüzeyden: kavisli yüzeylerin tam dört kenarı olmalı",
    "From face: this face could not be meshed": "Yüzeyden: bu yüzeyden kafes yapılamadı",
    "From face: the selected edges must form closed loops": "Yüzeyden: seçili kenarlar kapalı halka oluşturmalı",
    "From face: no flat face can be made from the selected edges":
        "Yüzeyden: seçili kenarlardan düz yüzey yapılamıyor",
    "From face: select faces or a closed edge loop of a Part / PartDesign object":
        "Yüzeyden: bir Part / PartDesign nesnesinin yüzeylerini ya da kapalı kenar halkasını seçin",
    "From face: faces with more than one hole are not supported yet":
        "Yüzeyden: birden fazla delikli yüzey henüz desteklenmiyor",
    "From face: select a flat face of a Part / PartDesign object (in object mode)":
        "Yüzeyden: bir Part / PartDesign nesnesinin düz yüzeyini seçin (nesne modunda)",
    "From face: %d polygons, %d points": "Yüzeyden: %d çokgen, %d nokta",
    "Face": "Yüzey",
    "Edge": "Kenar",
    "Point": "Nokta",
    "From edges": "Kenarlardan",
    "From points": "Noktalardan",
    "From face: select edges of a Part / PartDesign object":
        "Yüzeyden: bir Part / PartDesign nesnesinin kenarlarını seçin",
    "From face: select 3 or more points of a Part / PartDesign object":
        "Yüzeyden: bir Part / PartDesign nesnesinin 3 ya da daha fazla noktasını seçin",
    "From face: select 3 or more points that are not on one line":
        "Yüzeyden: aynı doğru üzerinde olmayan 3 ya da daha fazla nokta seçin",
    "From face: the selected edges branch (more than two meet at a point)":
        "Yüzeyden: seçili kenarlar dallanıyor (bir noktada ikiden fazla kenar buluşuyor)",
    "From face: no polygon can be made from the selected edges":
        "Yüzeyden: seçili kenarlardan çokgen yapılamıyor",
    "Snap: %s": "Yakala: %s",
    "vertex": "köşe",
    "edge": "kenar",
    "Hole": "Delik",
    "Radius": "Yarıçap",
    "Support": "Destek",
    "Hole: select faces": "Delik: yüz seçin",
    "Circle on face": "Yüzeye daire",
    "Circle on face: select faces": "Yüzeye daire: yüz seçin",
    "Circle on face: each selected patch must be closed, without holes inside":
        "Yüzeye daire: seçili her parça kapalı ve içi deliksiz olmalı",
    "Circle on face: the selected faces are not flat": "Yüzeye daire: seçili yüzler düz değil",
    "Circle on face: radius too large for the selected faces (max %.2f mm)":
        "Yüzeye daire: yarıçap seçili yüzler için çok büyük (en çok %.2f mm)",
    "Circle on face: radius %.2f mm, %d segments": "Yüzeye daire: yarıçap %.2f mm, %d segment",
    "Hole: each selected patch must be closed, without holes inside":
        "Delik: seçili her parça kapalı ve içi deliksiz olmalı",
    "Hole: the selected faces are not flat": "Delik: seçili yüzler düz değil",
    "Hole: radius too large for the selected faces (max %.2f mm)":
        "Delik: yarıçap seçili yüzlere sığmıyor (en çok %.2f mm)",
    "Hole: radius %.2f mm, %d segments": "Delik: yarıçap %.2f mm, %d segment",
    "Reduce points": "Nokta azalt",
    "Merge points": "Noktaları birleştir",
    "Merge points: %d merged": "Noktaları birleştir: %d nokta birleşti",
    "Fixed distance": "Sabit mesafe",
    "SubD tools": "SubD araçları",
    "Matcap": "Matcap",
    "Environment reflection": "Ortam yansıtma",
    "Environment image": "Ortam resmi",
    "Environment image: %s": "Ortam resmi: %s",
    "Environment image (*.hdr *.jpg *.jpeg *.png)": "Ortam resmi (*.hdr *.jpg *.jpeg *.png)",
    "Load image...": "Resim yükle...",
    "Roughness": "Pürüzlülük",
    "Surface color": "Yüzey rengi",
    "Cage color": "Kafes rengi",
    "Follow sketch": "Sketch'i izle",
    "From edge polygons follow their source sketch / edges when it changes; off = unlinked (the cage keeps its shape)":
        "From edge çokgenleri kaynak sketch / kenarlar değişince onları izler; kapalı = bağ yok (kafes şeklini korur)",
    "Sketch link: %s is gone": "Sketch bağı: %s yok",
    "Sketch link: %s.%s is gone": "Sketch bağı: %s.%s yok",
    "Sketch link: the curve no longer gives the same points": "Sketch bağı: eğri artık aynı noktaları vermiyor",
    "Sketch link: the linked points followed the sketch": "Sketch bağı: bağlı noktalar sketch'i izledi",
    "Sketch link ended: none of the linked points is where it was":
        "Sketch bağı bitti: bağlı noktaların hiçbiri yerinde değil",
    "Color of the cage faces when Subpatch is off": "Subpatch kapalıyken kafes yüzlerinin rengi",
    "Surface color for the plastic look (no effect at Metal 1)": "Plastik görünümün yüzey rengi (Metal 1'de etkisiz)",
    "0: mirror, 1: fully blurred reflections": "0: ayna, 1: tamamen bulanık yansıma",
    "0: plastic (surface color with a soft reflection), 1: metal (reflection only)":
        "0: plastik (yüzey rengi ve hafif yansıma), 1: metal (yalnız yansıma)",
    "The smooth preview reflects an image of the surroundings: an equirectangular (2:1) panorama, best as .hdr; JPG and PNG also work":
        "Pürüzsüz önizleme, ortam resmini yansıtır: eşdikdörtgen (2:1) panorama; en iyisi .hdr, JPG ve PNG de olur",
    "Off": "Kapalı",
    "Clay": "Kil",
    "Plastic": "Plastik",
    "Black plastic": "Siyah plastik",
    "Ceramic": "Seramik",
    "Metal": "Metal",
    "Surface look of the smooth preview (Edit mode)": "Pürüzsüz önizlemenin yüzey görünümü (Edit modu)",
    "Merge polygons": "Çokgenleri birleştir",
    "Parametric group": "Parametrik grup",
    "Make the selected points a group with a datum plane handle (SubD_PG_0001) in the FreeCAD tree: move the datum or attach it (Attachment, offset, expressions) and the points follow":
        "Seçili noktaları FreeCAD ağacında bir datum plane tutamağı (SubD_PG_0001) olan gruba çevirir: datum'u taşıyınca ya da bağlayınca (Attachment, ofset, ifadeler) noktalar onu izler",
    "Group handle %s is gone: the group is released": "Grup tutamağı %s yok: grup serbest bırakıldı",
    "Group: the points followed the handle": "Grup: noktalar tutamağı izledi",
    "Group released: none of its points is where it was": "Grup serbest kaldı: noktalarından hiçbiri yerinde değil",
    "Group: some of the points already belong to a group": "Grup: bazı noktalar zaten bir grupta",
    "SubD group": "SubD grubu",
    "Group distance": "Grup mesafesi",
    "Group distance: the Distance value can be edited here": "Grup ölçüsü: Distance değeri burada değiştirilebilir",
    "Unlock": "Kilidi aç",
    "Lock": "Kilitle",
    "Add to group": "Gruba ekle",
    "Remove from group": "Gruptan çıkar",
    "Done": "Tamam",
    "Group %s is locked: Unlock it in the group panel to edit its points": "%s grubu kilitli: noktalarını düzenlemek için grup panelinden kilidi açın",
    "Group %s unlocked: its points can be edited, Lock binds them again": "%s grubunun kilidi açıldı: noktaları düzenlenebilir, Kilitle yeniden bağlar",
    "Group %s released: its points or handle changed, it cannot be locked again": "%s grubu serbest kaldı: noktaları ya da tutamağı değişti, yeniden kilitlenemez",
    "Group %s locked": "%s grubu kilitlendi",
    "Group: the selected points already belong to a group": "Seçili noktalar zaten bir gruba ait",
    "Group %s: %d points added": "%s grubu: %d nokta eklendi",
    "Group %s: none of the selected points belongs to it": "%s grubu: seçili noktaların hiçbiri bu gruba ait değil",
    "Group %s released: no points left": "%s grubu serbest kaldı: nokta kalmadı",
    "Group %s: %d points removed": "%s grubu: %d nokta çıkarıldı",
    "Group %s: the handle could not follow": "%s grubu: tutamak noktaları izleyemedi",
    "Distance %s - %s (mm)": "Mesafe %s - %s (mm)",
    "SubD keep handle place": "SubD tutamağı yerini koru",
    "Group handle %s removed": "Grup tutamağı %s silindi",
    "Group handle %s created again": "Grup tutamağı %s yeniden oluşturuldu",
    "Group %s: %d points": "Grup %s: %d nokta",
    "Merge polygons: select at least two adjacent polygons": "Çokgenleri birleştir: en az iki bitişik çokgen seçin",
    "Merge polygons: the selected polygons are not adjacent": "Çokgenleri birleştir: seçili çokgenler bitişik değil",
    "Merge polygons: no single outline (hole, closed shape or mixed directions: try Flip)":
        "Çokgenleri birleştir: tek dış sınır çıkmıyor (delik, kapalı şekil ya da karışık yön: Çevir'i deneyin)",
    "Merge polygons: %d polygons merged": "Çokgenleri birleştir: %d çokgen birleşti",
    "Weld": "Kaynak",
    "Weld: select at least two points": "Kaynak: en az iki nokta seçin",
    "Weld: Shift+click the target point last": "Kaynak: hedef noktaya en son Shift+tıklayın",
    "Weld: %d points welded": "Kaynak: %d nokta birleşti",
    "Rail bridge": "Kılavuzlu köprü",
    "Rail bridge: select two polygon regions or edge loops, or one region to extend":
        "Kılavuzlu köprü: iki çokgen bölgesi / kenar halkası ya da uzatmak için tek bölge seçin",
    "Rail bridge: click a point and move it, I: add segment, Del: delete point, Enter: done":
        "Kılavuzlu köprü: noktaya tıklayıp taşıyın, I: segment ekle, Del: nokta sil, Enter: bitti",
    "Rail bridge: the tube was changed, its guide curve can no longer be edited":
        "Kılavuzlu köprü: boru değiştirilmiş, kılavuz eğrisi artık düzenlenemiyor",
    "Rail bridge: this point cannot be deleted": "Kılavuzlu köprü: bu nokta silinemez",
    "Add segment (I)": "Segment ekle (I)",
    "Delete point (Del)": "Nokta sil (Del)",
    "Done (Enter)": "Bitti (Enter)",
    "Move neighbors": "Komşular da gelsin",
    "Elasticity": "Esneklik",
    "Tension": "Gerginlik",
    "Rail Pipe": "Raylı boru",
    "Rail Groove": "Raylı kanal",
    "PartDesign: add a pipe along a guide curve from a flat face or closed sketch to a second one or a free end. Edit GuidePoints, Tension and EndScale in the property view":
        "PartDesign: düz bir yüzden ya da kapalı sketch'ten ikincisine ya da serbest uca, kılavuz eğri boyunca boru ekler. GuidePoints, Tension ve EndScale özellik panelinden düzenlenir",
    "PartDesign: cut a groove along a guide curve from a flat face or closed sketch (and to a second one). Edit GuidePoints, Tension and EndScale in the property view":
        "PartDesign: düz bir yüzden ya da kapalı sketch'ten (ve ikincisine) kılavuz eğri boyunca kanal oyar. GuidePoints, Tension ve EndScale özellik panelinden düzenlenir",
    "Rail Pipe: select one or two flat faces of a PartDesign Body.\n":
        "Raylı boru: bir PartDesign Body'nin bir ya da iki düz yüzünü seçin.\n",
    "Rail Pipe: the faces must belong to a PartDesign Body.\n": "Raylı boru: yüzler bir PartDesign Body'ye ait olmalı.\n",
    "Rail Pipe: the start face must be flat.\n": "Raylı boru: başlangıç yüzü düz olmalı.\n",
    "Rail Pipe: the first sweep failed; adjust GuidePoints in the property view.\n":
        "Raylı boru: ilk süpürme başarısız; GuidePoints'i özellik panelinden düzenleyin.\n",
    "Rail Pipe: select one or two flat faces or closed sketches of a PartDesign Body.\n":
        "Raylı boru: bir PartDesign Body'nin bir ya da iki düz yüzünü veya kapalı sketch'ini seçin.\n",
    "Edit Rail Pipe": "Raylı boruyu düzenle",
    "Click a guide point in the 3D view and move it with the gizmo (arrows: along an axis, squares: in a plane, cube: free) or drag the point itself. The white end points belong to the profiles and stay put. Click the curve to pick a span: Add point puts the new point where you clicked.":
        "3B görünümde bir kılavuz noktasına tıklayın ve gizmo ile taşıyın (oklar: eksen boyunca, kareler: düzlemde, küp: serbest) ya da noktanın kendisini sürükleyin. Beyaz uç noktalar profillere aittir, yerinde kalır. Bir eğri parçası seçmek için eğriye tıklayın: Nokta ekle yeni noktayı tıkladığınız yere koyar.",
    "Selected point": "Seçili nokta",
    "Add point": "Nokta ekle",
    "Delete point": "Nokta sil",
    "Guide": "Kılavuz",
    "The guide follows the dragged point like a cable": "Kılavuz, sürüklenen noktayı kablo gibi izler",
    "End scale": "Uç ölçeği",
    "Reversed": "Ters",
    "Rail Pipe: the sketch has no closed profile": "Raylı boru: sketch'te kapalı profil yok",
    "Rail Pipe: the start face is missing": "Raylı boru: başlangıç yüzü bulunamadı",
    "Rail Pipe: the start face must be flat": "Raylı boru: başlangıç yüzü düz olmalı",
    "Rail Pipe: a free end needs at least one guide point": "Raylı boru: serbest uç için en az bir kılavuz noktası gerekir",
    "Rail Pipe: two guide points are on top of each other": "Raylı boru: iki kılavuz noktası üst üste",
    "Rail Pipe: the sweep failed (%s)": "Raylı boru: süpürme başarısız (%s)",
    "Rail Pipe: the sweep gave no solid": "Raylı boru: süpürmeden katı çıkmadı",
    "Rail Pipe: the sweep gave an invalid solid": "Raylı boru: süpürme geçersiz bir katı verdi (eğri çok dar kıvrılıyor olabilir)",
    "Magnet": "Mıknatıs",
    "Magnet: drag to pull the surface softly, a selection limits it (right drag: radius; Esc / right click: off)":
        "Mıknatıs: yüzeyi yumuşakça çekmek için sürükleyin, seçim varsa yalnız onu etkiler (sağ tuşla sürükleme: yarıçap; Esc / sağ tık: kapat)",
    "Radius %.3f mm": "Yarıçap %.3f mm",
    "Spherical": "Küresel",
    "Radial": "Radyal",
    "Magnet: no points within the radius": "Mıknatıs: yarıçap içinde nokta yok",
    "Smooth": "Yumuşak",
    "Linear": "Doğrusal",
    "Peak": "Sivri",
    "Dome": "Kubbe",
    "Align": "Hizala",
    "Align horizontally: the selected points get the height of the first clicked one on screen":
        "Yatay hizala: seçili noktalar ekranda ilk tıklanan noktanın yüksekliğine gelir",
    "Align vertically: the selected points line up with the first clicked one on screen":
        "Dikey hizala: seçili noktalar ekranda ilk tıklanan noktanın altına / üstüne dizilir",
    "Align: select two or more points": "Hizala: iki ya da daha fazla nokta seçin",
    "Align: select two or more separate polygons or edges":
        "Hizala: birbirinden ayrı iki ya da daha fazla poligon veya kenar seçin",
    "Align horizontally: selected points (polygons, edges: their centers) get the height of the first clicked one on screen":
        "Yatay hizala: seçili noktalar (poligon ve kenarlarda merkezleri) ekranda ilk tıklananın yüksekliğine gelir",
    "Align vertically: selected points (polygons, edges: their centers) line up with the first clicked one on screen":
        "Dikey hizala: seçili noktalar (poligon ve kenarlarda merkezleri) ekranda ilk tıklananın altına / üstüne dizilir",
    "Align: Shift-click two or more guide points (the end points too)":
        "Hizala: Shift ile iki ya da daha fazla kılavuz noktası seçin (uç noktalar da olur)",
    "Align horizontally (Shift-click guide points, the end points too)":
        "Yatay hizala (kılavuz noktalarını Shift ile seçin, uç noktalar da olur)",
    "Align vertically (Shift-click guide points, the end points too)":
        "Dikey hizala (kılavuz noktalarını Shift ile seçin, uç noktalar da olur)",
    "Align horizontally: %s":"Yatay hizalandı: %s",
    "Align vertically: %s": "Dikey hizalandı: %s",
    "Scale": "Ölçek",
    "Flip": "Çevir",
    "Flip: %d polygons": "Çevir: %d çokgen",
    "Auto direction": "Otomatik yön",
    "Reduce points: %d removed": "Nokta azalt: %d nokta çıkarıldı",
    "Degrees": "Derece",
    "Delete points": "Noktaları sil",
    "Bridge": "Köprü",
    "Number of Cuts": "Kesim sayısı",
    "Twist": "Burgu",
    "Circle": "Daire",
    "Strength": "Güç",
    "Auto radius": "Otomatik yarıçap",
    "Even spacing": "Eşit aralık",
    "Circle: only in edit mode": "Daire: yalnız düzenleme modunda",
    "Circle: select polygons that form one region without holes, an edge loop or points": "Daire: delik içermeyen tek bölge oluşturan çokgenler, bir kenar halkası ya da noktalar seç",
    "Circle: select one edge loop or chain, polygons or points": "Daire: tek bir kenar halkası ya da zinciri, çokgenler veya noktalar seç",
    "Circle: select at least 3 points": "Daire: en az 3 nokta seç",
    "Circle: the points do not go around a centre": "Daire: noktalar bir merkezin etrafında dönmüyor",
    "Circle: the radius is zero": "Daire: yarıçap sıfır",
    "Circle: radius %.3f mm": "Daire: yarıçap %.3f mm",
    "Arrange the selected points, an edge loop / chain or the boundary of a polygon region on a circle (inner points follow). Type the Radius in the panel (Auto radius = average of the current one)": "Seçili noktaları, bir kenar halkasını / zincirini ya da çokgen bölgenin sınırını daireye dizer (iç noktalar eşlik eder). Yarıçapı panelden sayıyla yaz (Otomatik yarıçap = mevcut yarıçapların ortalaması)",
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
    "SubD mode": "SubD modu",
    "Edit the selected SubD object in the SubD modeler": "Seçili SubD nesnesini SubD modelleyicide düzenle",
    "Angle": "Açı",
    "Width": "Genişlik",
    "Depth": "Derinlik",
    "Height": "Yükseklik",
    "Segments X": "Segment X",
    "Segments Y": "Segment Y",
    "Segments Z": "Segment Z",
    "Box: size / segments in the panel": "Kutu: ölçü / segment panelde",
    "SubD Modeler is ready: Create > Box or Cylinder to start":
        "SubD Modeler hazır: başlamak için Oluştur > Kutu ya da Silindir",
    "Segments": "Segment",
    "SubD surface": "SubD yüzeyi",
    "Edit level": "Düzenleme seviyesi",
    "Final level": "Son seviye",
    "Subdivision level shown while editing (lower = faster)":
        "Düzenlerken gösterilen bölme seviyesi (düşük = hızlı)",
    "Subdivision level of the SubD object in the document": "Belgedeki SubD nesnesinin bölme seviyesi",
    "Open edges": "Açık kenarlar",
    "Edge and corner": "Kenar ve köşe",
    "Edge only": "Yalnız kenar",
    "Open edges: Edge and corner = one-face corners stay sharp, Edge only = they round too":
        "Açık kenarlar: Kenar ve köşe = tek yüzlü köşeler sivri kalır, Yalnız kenar = onlar da yuvarlanır",
    "Smooth triangles": "Üçgen yumuşatma",
    "Loop-like weights next to triangles: smoother triangle fans (e.g. cylinder caps)":
        "Üçgenlerin yanında Loop benzeri ağırlıklar: üçgen yelpazeler (ör. silindir kapakları) daha düzgün",
    "Zebra stripes": "Zebra şeritleri",
    "Reflection stripes on the surface: kinks show as broken stripes (G0), corners in the stripes as G1, "
    "smooth stripes as G2":
        "Yüzeye yansıyan şeritler: kırık şerit = kırık yüzey (G0), şeritte köşe = G1, yumuşak şerit = G2",
    "Solid (STEP)": "Katı (STEP)",
    "Gizmo": "Gizmo",
    "Default": "Varsayılan",
    "Length of the gizmo arrows on screen (default %d px)": "Gizmo oklarının ekrandaki boyu (varsayılan %d px)",
    "Tool panel": "Araç paneli",
    "Left": "Sol",
    "Right": "Sağ",
    "Where the panel with the values of the running tool and of the selection sits in the 3D view: bottom left or "
    "bottom right":
        "Çalışan aracın ve seçimin değerlerini gösteren panelin 3B görünümdeki yeri: sol alt ya da sağ alt",
    "Precision": "Hassasiyet",
    "Continuity": "Süreklilik",
    "Draft": "Taslak",
    "Normal": "Normal",
    "Fine": "Hassas",
    "Draft / Normal / Fine: samples per patch, G2 refine levels and sewing tolerance":
        "Taslak / Normal / Hassas: yama başına örnek, G2 inceltme seviyesi ve dikiş toleransı",
    "G0: fast, every patch fitted (tiny kinks between patches)":
        "G0: hızlı, her yama yaklaşık oturtulur (yamalar arasında çok küçük kırıklar)",
    "G1: regular patches exact, the rest fitted to the limit surface":
        "G1: düzgün yamalar birebir, diğerleri limit yüzeye oturtulur",
    "G2: as G1, irregular patches refined into more exact ones (curvature continuous except at 3 / 5+ edge "
    "points, where SubD itself is G1)":
        "G2: G1 gibi, düzensiz yamalar inceltilip daha çok birebir yamaya bölünür (eğrilik sürekli; yalnız "
        "3 / 5+ kenarlı noktalarda SubD'nin kendisi de G1'dir)",
    "Solid: %s, %s": "Katı: %s, %s",
    "SubD edit %s": "SubD düzenleme %s",
    "Document undo: cage reloaded": "Belge geri alındı: kafes yeniden yüklendi",
    "Nothing to undo in this edit (O: object mode, Ctrl+Z undoes whole edits)":
        "Bu düzenlemede geri alınacak bir şey yok (O: Object modu, orada Ctrl+Z bütün düzenlemeyi geri alır)",
    "Knife: the cut does not go all the way around, the model cannot be split into objects":
        "Knife: kesim modelin çevresini tam dolaşmıyor, iki nesneye ayrılamaz",
    "Shell": "Kabuk",
    "Thickness": "Kalınlık",
    "Detail filter": "Detay filtresi",
    "Inside": "İçeri",
    "Outside": "Dışarı",
    "Centered": "Ortada",
    "Fill": "Doldur",
    "Sharp": "Keskin kenar",
    "Flat rim": "Düz kenar",
    "Shell: nothing to shell": "Kabuk: kabuklanacak yüz yok",
    "Shell: thickness must be greater than 0": "Kabuk: kalınlık 0'dan büyük olmalı",
    "Shell: thickness %.2f mm is too large for this part (max %.2f mm)":
        "Kabuk: %.2f mm kalınlık bu parça için fazla (en çok %.2f mm)",
    "Shell: select a SubD object or polygons": "Kabuk: bir SubD nesnesi ya da poligon seçin",
    "Shell: %d polygons, %d open edges to weave": "Kabuk: %d poligon, örülecek %d açık kenar",
    "  |  wall kept at half thickness at %d points": "  |  %d noktada duvar yarı kalınlıkta tutuldu",
    "Split in two": "İkiye ayır",
    "Split into two objects": "İki objeye ayır",
    "Knife: split into %s (%d polygons)": "Knife: %s nesnesine ayrıldı (%d poligon)",
    "Diameter": "Çap",
    "Height Segments": "Yükseklik segment",
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
    "Mirror: drag the mirror line, the other side is replaced by the mirror "
    "(Ctrl: 15° steps; Shift+V / Esc / right click: off)":
        "Ayna: ayna çizgisini sürükleyerek çek, öbür taraf aynayla değiştirilir "
        "(Ctrl: 15° adım; Shift+V / Esc / sağ tık: kapat)",
    "Mirror  N (%.3f, %.3f, %.3f)   Offset %.3f mm": "Ayna  N (%.3f, %.3f, %.3f)   Konum %.3f mm",
    "Mirror: no polygons on the source side (Flip side)": "Ayna: kaynak tarafta çokgen yok (Tarafı çevir)",
    "Mirror: %d polygons": "Ayna: %d çokgen",
    "Mirror (line)": "Ayna (çizgi)",
    "Mirror (plane: N.p = Offset)": "Ayna (düzlem: N.p = Konum)",
    "Flip side": "Tarafı çevir",
    "Individual polygons": "Her çokgen ayrı",
    "Stays on: draw a line, the model is mirrored across the plane through it; the far side is replaced "
    "(source: the selection's side, else left of a downward line; Flip side in the panel)":
        "Açık kalır: bir çizgi çiz, model içinden geçen düzleme göre aynalanır; öbür taraf değiştirilir "
        "(kaynak: seçimin olduğu taraf, yoksa aşağı çizilen çizginin solu; paneldeki Tarafı çevir)",
    "Loop: click edges to select loops (Esc / right click: off)":
        "Loop: halka seçmek için kenarlara tıkla (Esc / sağ tık: kapat)",
    "Knife: click edges to cut (K / Esc / right click: off)":
        "Knife: kesmek için kenarlara tıkla (K / Esc / sağ tık: kapat)",
    "Distance %.3f mm": "Mesafe %.3f mm",
    "Axis (none = normal)": "Eksen (boş = normal)",
    "Individual faces (own normals)": "Yüzler ayrı (kendi normalinde)",
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
    "Esc / right click: cancel":
        "1 / 2 / 3: nokta / kenar / yüz modu\n"
        "Y: X simetrisi aç/kapa\n"
        "R: gizmo Döndür  |  H: gizmo Stretch / Size\n"
        "C: seçili kenarları keskinleştir (crease)\n"
        "Tık: seç  |  Shift+Tık: ekle\n"
        "Sürükle: kutu seçimi (Shift ekle, Ctrl çıkar)\n"
        "L: seçimden halka (ya da imlecin altındaki)  |  X: hayalet\n"
        "E: extrude  |  K: knife (loop cut)  |  J: knife (çizgi / düzlem)  |  Shift+V: ayna (çizgi)\n"
        "Knife / loop açık kalır (sarı): K / J / Esc / sağ tık kapatır\n"
        "B: bevel (sağa inset, yukarı shift)  |  I: inset\n"
        "Kenar modunda B: kenar bevel (sağa = genişlik)\n"
        "Kenar modunda E: açık kenarları extrude  |  extrude'da tekerlek = segment, X / Y / Z = eksen\n"
        "O: nesne / düzenleme modu  |  D: nesneyi çoğalt\n"
        "T: taşıma gizmosu  |  G: serbest taşı (fare)  |  Tab / S: SubPatch\n"
        "A: tümünü seç  |  Del: sil  |  P: poligon ör\n"
        "Ctrl+C / X / V: çokgen kopyala / kes / yapıştır\n"
        "Ctrl+Z: geri al  |  Ctrl+Y: ileri al  |  F: belgeye kaydet\n"
        "Ctrl + gizmo / T: adımlı (mm, derece)\n"
        "Esc / sağ tık: işlemi iptal",
    # --- workbench commands ---
    "SubD Modeler": "SubD Modeler",
    "Open the SubD modeler (edits the selected SubD object, otherwise opens in object mode)":
        "SubD modelleyiciyi aç (seçili SubD nesnesini düzenler, yoksa Object modunda açılır)",
    "SubD to solid": "SubD'den katı",
    "Make a parametric B-rep solid (STEP) from the selected SubD object":
        "Seçili SubD nesnesinden parametrik B-rep katı (STEP) oluştur",
    "Select a SubD object first.\n": "Önce bir SubD nesnesi seçin.\n",
    "Interface language": "Arayüz dili",
    # --- tool tips ---
    "Nothing to redo":
        "İleri alınacak bir şey yok",
    "Undo":
        "Geri al",
    "Redo":
        "İleri al",
    "Edit": "Düzenle",
    "Object mode: move, copy and select whole SubD objects":
        "Nesne modu: SubD nesnelerini bütün olarak taşı, kopyala, seç",
    "Edit the cage of the selected SubD object":
        "Seçili SubD nesnesinin kafesini düzenle",
    "Write the cage to the document now (it is also saved automatically every 20 s)":
        "Kafesi şimdi belgeye yaz (20 sn'de bir kendiliğinden de kaydedilir)",
    "Make a parametric B-rep solid from the SubD surface, ready for STEP export and PartDesign":
        "SubD yüzeyinden parametrik B-rep katı oluştur (STEP ve PartDesign için hazır)",
    "Write the cage to the document and close the modeler":
        "Kafesi belgeye yaz ve modelleyiciyi kapat",
    "Add a box; set its size and segments in the panel":
        "Kutu ekle; boyut ve segmentler panelde",
    "Add a new cylinder object along Z":
        "Z ekseni boyunca yeni silindir nesnesi ekle",
    "Cylinder with triangle-fan caps":
        "Kapakları üçgen yelpaze olan silindir",
    "Cylinder with all-quad caps (segments are a multiple of 4)":
        "Kapakları tamamen dörtgen silindir (segment 4'ün katı)",
    "Make a cage from selected Part / PartDesign elements":
        "Seçili Part / PartDesign elemanlarından kafes oluştur",
    "Turn the selected faces into one cage":
        "Seçili yüzleri tek kafese çevir",
    "Each chain of selected edges becomes one polygon":
        "Seçili kenarların her zinciri bir çokgen olur",
    "The selected points become one polygon, in the order they were picked":
        "Seçili noktalar, seçim sırasıyla bir çokgen olur",
    "Make a polygon from 3 or more selected points, a closed edge loop, or two edges":
        "3 ya da daha çok seçili noktadan, kapalı kenar halkasından ya da iki kenardan çokgen oluştur",
    "Duplicate the selected objects":
        "Seçili nesneleri çoğalt",
    "Select everything, or deselect if all is already selected":
        "Tümünü seç; zaten hepsi seçiliyse seçimi kaldır",
    "Clear the selection; also cancels the running tool":
        "Seçimi temizle; çalışan aracı da iptal eder",
    "Stays on: click to select the loop under the cursor":
        "Açık kalır: imlecin altındaki halkayı seçmek için tıkla",
    "Pull out the selected faces or open edges; mouse wheel = segments, X / Y / Z = axis":
        "Seçili yüzleri ya da açık kenarları dışarı çek; tekerlek = segment, X / Y / Z = eksen",
    "Polygon mode: inset and shift the selected faces. Edge mode: round the selected edges":
        "Poligon modu: seçili yüzleri içeri al ve kaydır. Kenar modu: seçili kenarları yuvarla",
    "Inset the selected faces into their own plane":
        "Seçili yüzleri kendi düzlemlerinde içeri al",
    "Stays on: click an edge to cut a loop around the model, then slide it":
        "Açık kalır: modelin çevresinde halka kesmek için kenara tıkla, sonra kaydır",
    "Stays on: draw a line to cut the model live with the plane through it":
        "Açık kalır: çizgi çiz, model bu çizgiden geçen düzlemle canlı kesilir",
    "Connect two open edge loops or two face regions":
        "İki açık kenar halkasını ya da iki yüz bölgesini bağla",
    "Connect two open edge loops / chains or two face regions with a strip":
        "İki açık kenar halkasını / zincirini ya da iki yüz bölgesini şeritle bağla",
    "Pipe with an editable guide curve, between two regions or growing out of one":
        "Düzenlenebilir kılavuz eğrili boru: iki bölge arasında ya da tek bölgeden uzayan",
    "Cut a round hole in the middle of the selected faces and re-mesh around it":
        "Seçili yüzlerin ortasına yuvarlak delik aç, çevresini yeniden ör",
    "Like Hole, but the circle stays filled and selected, ready for Extrude or Bevel":
        "Hole gibi, ama daire dolu ve seçili kalır; Extrude ya da Bevel'e hazır",
    "Merge the selected adjacent polygons into one":
        "Seçili komşu çokgenleri tek çokgende birleştir",
    "Delete the selected elements (object mode: the selected objects)":
        "Seçili elemanları sil (Object modu: seçili nesneleri)",
    "Copy the selected polygons with their creases; edges and points are copied as free edges / points":
        "Seçili çokgenleri keskinlikleriyle kopyala; kenarlar ve noktalar serbest kenar / nokta olarak kopyalanır",
    "Copy the selected polygons, edges or points, then remove them":
        "Seçili çokgen, kenar ya da noktaları kopyala, sonra kaldır",
    "Paste polygons, edges and points in place into this cage":
        "Çokgen, kenar ve noktaları bu kafese aynı yerine yapıştır",
    "Paste the copied polygons, edges and points as a new SubD object":
        "Kopyalanan çokgen, kenar ve noktaları yeni SubD nesnesi olarak yapıştır",
    "Move gizmo: arrows along X / Y / Z, squares in a plane, cube in the screen plane":
        "Taşıma gizmosu: oklar X / Y / Z boyunca, kareler bir düzlemde, küp ekran düzleminde",
    "Rotate gizmo":
        "Döndürme gizmosu",
    "Stretch the selection along an axis or size it uniformly":
        "Seçimi bir eksende uzat ya da orantılı boyutlandır",
    "Stays on: drag on the cage to pull it softly inside the radius; a selection limits it":
        "Açık kalır: yarıçap içinde yumuşakça çekmek için kafesi sürükle; seçim etkiyi sınırlar",
    "Make the selected edges sharp, or soft again; exact dimensions need sharp edges":
        "Seçili kenarları keskinleştir ya da yeniden yumuşat; tam ölçü için keskin kenar gerekir",
    "Remove in-between points along nearly straight edges (selection or whole cage)":
        "Neredeyse düz kenarlardaki ara noktaları çıkar (seçim ya da bütün kafes)",
    "Merge points at the same position (selection or whole cage)":
        "Aynı konumdaki noktaları birleştir (seçim ya da bütün kafes)",
    "Weld the selected points into the last one clicked":
        "Seçili noktaları son tıklanan noktada birleştir",
    "Flip the normals of the selected polygons (no selection: all)":
        "Seçili çokgenlerin normallerini çevir (seçim yoksa hepsi)",
    "Give the object or selected polygons a wall thickness":
        "Nesneye ya da seçili çokgenlere et kalınlığı ver",
    "Mirror the selected polygons (no selection: whole cage) across a world plane and weld; object mode: mirrored copies":
        "Seçili çokgenleri (seçim yoksa bütün kafes) dünya düzlemine göre aynala ve birleştir; Object modu: aynalı kopya",
    "Mirror across the plane perpendicular to X":
        "X'e dik düzleme göre aynala",
    "Mirror across the plane perpendicular to Y":
        "Y'ye dik düzleme göre aynala",
    "Mirror across the plane perpendicular to Z":
        "Z'ye dik düzleme göre aynala",
    "Undo one step (object mode: one whole edit)":
        "Bir adım geri al (Object modu: bütün bir düzenleme)",
    "Redo the step you undid (object mode: one whole edit)":
        "Geri alınan adımı ileri al (Object modu: bütün bir düzenleme)",
    "Point selection mode":
        "Nokta seçim modu",
    "Edge selection mode":
        "Kenar seçim modu",
    "Polygon selection mode":
        "Çokgen seçim modu",
    "X symmetry: tools also act on the mirrored side":
        "X simetrisi: araçlar aynalı tarafa da uygulanır",
    "Show the smooth SubD surface or only the cage":
        "Yumuşak SubD yüzeyini ya da yalnız kafesi göster",
    "See through the model and select hidden elements":
        "Modelin içini gör, gizli elemanları seç",
    "Free moves snap to the vertices and edges of visible Part / PartDesign shapes":
        "Serbest taşımalar görünür Part / PartDesign şekillerinin köşe ve kenarlarına yapışır",
    "Show or hide the Options panel":
        "Options panelini göster ya da gizle",
    "Array": "Dizi",
    "Rectangular": "Dikdörtgen",
    "Between": "Arası",
    "Total": "Toplam",
    "Multiply": "Çarpan",
    "Factor": "Çarpan katsayısı",
    "Rotate copies": "Kopyaları döndür",
    "Count X": "Adet X",
    "Count Y": "Adet Y",
    "Count Z": "Adet Z",
    "Offset X": "Aralık X",
    "Offset Y": "Aralık Y",
    "Offset Z": "Aralık Z",
    "Count": "Adet",
    "Start angle": "Başlangıç açısı",
    "End angle": "Bitiş açısı",
    "Offset along axis": "Eksen boyunca aralık",
    "Axis point X": "Eksen noktası X",
    "Axis point Y": "Eksen noktası Y",
    "Axis point Z": "Eksen noktası Z",
    "Axis direction X": "Eksen yönü X",
    "Axis direction Y": "Eksen yönü Y",
    "Axis direction Z": "Eksen yönü Z",
    "Axis X": "Eksen X",
    "Axis Y": "Eksen Y",
    "Axis Z": "Eksen Z",
    "Array: only in edit mode": "Dizi: yalnız düzenleme modunda",
    "Array: select polygons": "Dizi: çokgen seç",
    "Array: the axis direction is zero": "Dizi: eksen yönü sıfır",
    "Array: too many copies (max %d)": "Dizi: çok fazla kopya (en çok %d)",
    "Array: nothing to add (raise Count or set an Offset)": "Dizi: eklenecek bir şey yok (Adet'i artır ya da aralık ver)",
    "Array: %d copies, %d polygons": "Dizi: %d kopya, %d çokgen",
    "Array: drag the axis line, a click gives an axis along the view (Ctrl: 15° steps; Esc / right click: off)": "Dizi: eksen çizgisini sürükle, tıklamak bakış yönünde eksen verir (Ctrl: 15° adım; Esc / sağ tık: kapat)",
    "Draw axis": "Eksen çiz",
    "Copy the selected polygons as separate parts: around an axis (Radial) or along X / Y / Z (Rectangular); switch in the panel": "Seçili çokgenleri ayrı parça olarak kopyalar: eksen etrafında (Radyal) ya da X / Y / Z boyunca (Dikdörtgen); mod panelden seçilir",
    "Drag a line in the view to set the axis (snaps like Mirror (line)); a click gives an axis along the view": "Eksen için görünümde çizgi sürükle (Ayna (çizgi) gibi yakalar); tıklamak bakış yönünde eksen verir",
    "Axis direction along world X (the axis point stays)": "Eksen yönü dünya X boyunca (eksen noktası kalır)",
    "Axis direction along world Y (the axis point stays)": "Eksen yönü dünya Y boyunca (eksen noktası kalır)",
    "Axis direction along world Z (the axis point stays)": "Eksen yönü dünya Z boyunca (eksen noktası kalır)",
    "Ease in": "Yavaş başla",
    "Ease out": "Yavaş bitir",
    "Reverse": "Ters çevir",
    "Auto range": "Otomatik aralık",
    "Start": "Başlangıç",
    "End": "Bitiş",
    "Cuts": "Kesim",
    "Twist: only in edit mode": "Burma: yalnız düzenleme modunda",
    "Twist: select polygons": "Burma: çokgen seç",
    "Twist: the axis direction is zero": "Burma: eksen yönü sıfır",
    "Twist: the selection has no length along the axis (set Start and End)": "Burma: seçimin eksen boyunca uzunluğu yok (Başlangıç ve Bitiş'i ayarla)",
    "Twist: move the mouse to twist (Ctrl: steps), click to accept, Esc to cancel": "Burma: fareyi oynatarak burk (Ctrl: adım), tıkla onayla, Esc vazgeç",
    "Twist %.2f°   %s": "Burma %.2f°   %s",
    "(symmetry is ignored)": "(simetri yok sayılır)",
    "Twist the selected polygons around an axis: the angle grows along the axis by the Falloff curve. Move the mouse for the angle (Ctrl: steps), click to accept, then fine-tune in the panel": "Seçili çokgenleri bir eksen etrafında burar: açı eksen boyunca düşüş eğrisine göre artar. Açı için fareyi oynat (Ctrl: adım), tıklayınca onaylanır, sonra panelden ince ayar yap",
    "Keep (fixed)": "Kalsın (sabit)",
    "Clear all": "Hepsini sil",
    "Distance: click two points (Esc / right click: off)": "Mesafe: iki noktaya tıkla (Esc / sağ tık: kapat)",
    "Distance: click the second point (Esc: cancel)": "Mesafe: ikinci noktaya tıkla (Esc: vazgeç)",
    "Distance: nothing to snap to here": "Mesafe: burada yakalanacak bir şey yok",
    "Distance %.3f mm   dX %.3f   dY %.3f   dZ %.3f": "Mesafe %.3f mm   dX %.3f   dY %.3f   dZ %.3f",
    "Measure the distance between two points: cage points, edge middles / edges, polygon centers, other objects (corners / edges) or the surface. Keep: the measurement stays, follows the points and is deleted with the × next to its text; otherwise it is temporary. Works in Object mode too (temporary)": "İki nokta arasındaki mesafeyi ölçer: kafes noktaları, kenar ortaları / kenarlar, çokgen merkezleri, diğer nesnelerin köşe / kenarları ya da yüzey. Kalsın açıkken ölçü ekranda kalır, noktalarla birlikte değişir ve yazısının yanındaki × ile silinir; kapalıyken geçicidir. Object modunda da çalışır (geçici)",
    "Delete all distances (undo brings the fixed ones back)": "Bütün ölçüleri sil (sabit olanlar geri alınabilir)",
    "Distances": "Mesafeler",
    "Name": "Ad",
    "Target": "Hedef",
    "Tol.": "Tol.",
    "Dev.": "Sapma",
    "Show names": "Adları göster",
    "Show all": "Hepsini göster",
    "Hide all": "Hepsini gizle",
    "Export CSV...": "CSV olarak kaydet...",
    "Distance %d": "Mesafe %d",
    "Show the name in front of the distance in the 3D view": "Adı 3B görünümde mesafenin önünde göster",
    "Show / hide this distance in the 3D view": "Bu mesafeyi 3B görünümde göster / gizle",
    "Delete this distance": "Bu mesafeyi sil",
    "Name / Target / Tol. are editable (double click); the selected row is drawn blue in the view": "Ad / Hedef / Tol. çift tıklayarak düzenlenir; seçili satır görünümde mavi çizilir",
    "Object mode: all distances of all SubD objects (edit them in Edit mode)": "Object modu: bütün SubD nesnelerinin mesafeleri (düzenlemek için Edit moduna geç)",
    "Distances: nothing to export": "Mesafeler: kaydedilecek bir şey yok",
    "Export distances": "Mesafeleri kaydet",
    "Distances exported: %s": "Mesafeler kaydedildi: %s",
    "Distances: could not write the file (%s)": "Mesafeler: dosya yazılamadı (%s)",
    "List of the distances: give each a name and a target, show / hide, delete, export as CSV. The list stays with the object between sessions": "Mesafelerin listesi: her birine ad ve hedef ver, göster / gizle, sil, CSV olarak kaydet. Liste sonraki oturumlarda nesneyle birlikte kalır",
}
