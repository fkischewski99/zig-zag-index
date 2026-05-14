# Zigzag-Index: ZO-1 Zell-Zell-Kontakt-Analyse

Misst den **Zickzack-Index** von ZO-1-markierten Tight Junctions in Fluoreszenzbildern.
Der Zickzack-Index ist das Verhaeltnis von tatsaechlicher Membranlaenge zu einer Referenzstrecke zwischen zwei Tricellular Junctions (1.0 = gerade, >1.0 = gewellt/gezackt).

## Dateien

| Datei | Beschreibung |
|---|---|
| `Zigzag_Index.ijm` | ImageJ-Makro -- oeffnet LIF-Dateien, User setzt Punkte, Python misst Pfade |
| `Zigzag_Index_Auswahl.ijm` | Wie oben, aber mit Dialog zur Auswahl einzelner Serien |
| `zigzag_pathfind.py` | Python-Skript -- Dijkstra-Pfadmessung entlang der Membran |

## Voraussetzungen

### ImageJ / Fiji
- [Fiji](https://fiji.sc/) (ImageJ mit Plugins)
- **Bio-Formats**-Plugin (in Fiji standardmaessig enthalten)

### Python (>= 3.8)

```bash
pip install numpy scikit-image
```

## Installation

1. `zigzag_pathfind.py` in den ImageJ-Macros-Ordner kopieren (oder beliebig -- das Makro fragt nach dem Pfad, falls es die Datei nicht findet)
2. `Zigzag_Index.ijm` ueber **Plugins > Macros > Run...** ausfuehren, oder ebenfalls in den Macros-Ordner legen

Falls Python nicht als `python` erreichbar ist, den Pfad in Zeile 25 des Makros anpassen:
```
pythonCmd = "python3";   // oder "/usr/bin/python3" etc.
```

## Ablauf

```
LIF-Datei oeffnen
    |
    v
Kanal waehlen (ZO-1)
    |
    v
Referenzmodus waehlen: Gerade (Euklidisch) oder Kreisbogen
    |
    v
Pro Serie:
    |
    +-- 1. Z-Projektion (Max Intensity)
    |
    +-- 2. User setzt Punkte paarweise auf trizellulare Junctions
    |       Punkt 1+2 = erster Kontakt, Punkt 3+4 = zweiter, ...
    |
    +-- 3. Pfadmessung (Python, Dijkstra)
    |       - Kostenkarte aus Bildintensitaet (hell = billig)
    |       - Optimaler Pfad entlang der Membran
    |       - Zickzack-Index = Pfadlaenge / Referenzstrecke
    |
    +-- 4. Overlay-Visualisierung
    |       gruen = Membranpfad, rot = gerade Strecke,
    |       cyan = Kreisbogen (optional)
    |
    v
CSV-Export aller Ergebnisse
```

## Mathematik

### Euklidischer Zickzack-Index (Standard)

Der klassische Zickzack-Index vergleicht die Membranpfadlaenge mit der geraden Strecke:

```
ZI_euklid = L_membran / d_euklid
```

wobei:

```
d_euklid = sqrt((x2 - x1)^2 + (y2 - y1)^2)
```

**Problem:** Bei gebogenen Membranen (z.B. entlang einer gekruemmten Gewebegrenze) ist die gerade Strecke deutlich kuerzer als der Bogen. Der Index wird kuenstlich erhoeht, obwohl die Membran gar nicht gezackt ist.

### Kreisbogen-Zickzack-Index (optional)

Um die Makrokruemmung herauszurechnen, kann ein **Kreisbogen** als Referenz verwendet werden:

```
ZI_bogen = L_membran / L_bogen
```

#### Schritt 1: Least-Squares Circle Fit (Kasa-Methode)

Aus den N Dijkstra-Pfadpunkten (x_i, y_i) wird ein Kreis (a, b, r) gefittet.

Ansatz: (x - a)^2 + (y - b)^2 = r^2 wird linearisiert zu:

```
x_i^2 + y_i^2 = 2a * x_i + 2b * y_i + c
```

mit c = r^2 - a^2 - b^2. Das ergibt das ueberbestimmte Gleichungssystem:

```
| x_1  y_1  1 |   | 2a |     | x_1^2 + y_1^2 |
| x_2  y_2  1 | * | 2b |  =  | x_2^2 + y_2^2 |
| ...  ...  . |   | c  |     | ...            |
| x_N  y_N  1 |              | x_N^2 + y_N^2 |
```

Loesung via Least Squares (numpy.linalg.lstsq):
- Mittelpunkt: a = result[0] / 2, b = result[1] / 2
- Radius: r = sqrt(c + a^2 + b^2)

#### Schritt 2: Bogenlaenge

Die Bogenlaenge zwischen den zwei Endpunkten berechnet sich aus Radius und Sehnenlaenge:

```
L_bogen = r * 2 * arcsin(d_euklid / (2r))
```

Falls d_euklid >= 2r (Entartungsfall), wird auf d_euklid zurueckgefallen.

#### Interpretation

- Bei gerader Membran: r -> unendlich, L_bogen -> d_euklid, also ZI_bogen = ZI_euklid
- Bei gebogener, aber glatter Membran: L_bogen > d_euklid, also ZI_bogen < ZI_euklid (naeher an 1.0)
- Bei gebogener UND gezackter Membran: ZI_bogen > 1.0, aber kleiner als ZI_euklid

Der Kreisbogen-Index misst somit nur die tatsaechliche Zackigkeit, unabhaengig von der Makrokruemmung.

## Parameter

Im Makro (`Zigzag_Index.ijm`, Zeile 25):

| Parameter | Default | Beschreibung |
|---|---|---|
| `pythonCmd` | `"python"` | Python-Befehl |

## Python-Skript: CLI

```bash
# Standard (nur Euklidisch):
python zigzag_pathfind.py <bild.tif> <koordinaten.csv> <ergebnisse.csv> <pfade.csv> <px_breite> <px_hoehe>

# Mit Kreisbogen-Referenz:
python zigzag_pathfind.py <bild.tif> <koordinaten.csv> <ergebnisse.csv> <pfade.csv> <px_breite> <px_hoehe> circle <bogen_pfade.csv>
```

Berechnet den Zickzack-Index fuer jedes Koordinatenpaar via Dijkstra-Pfadfindung. Mit dem optionalen `circle`-Argument wird zusaetzlich ein Kreisbogen-Fit berechnet.

## CSV-Ausgabe

### Ohne Kreisbogen (Standard)

```
Name,Kanal,Zickzack,Gerade,Index
"Serie1",C1,15.67,12.34,1.2698
"Serie1",C1,10.23,8.91,1.1481

,,Mittelwert,,1.2090
```

### Mit Kreisbogen

```
Name,Kanal,Zickzack,Gerade,Index,Bogen,BogenIndex
"Serie1",C1,15.67,12.34,1.2698,14.50,1.0807
"Serie1",C1,10.23,8.91,1.1481,9.80,1.0439

,,Mittelwert,,1.2090,,1.0623
```

| Spalte | Beschreibung |
|---|---|
| Zickzack | Tatsaechliche Membranlaenge (Dijkstra-Pfad) |
| Gerade | Euklidische Distanz (Sehnenlaenge) |
| Index | Zickzack / Gerade |
| Bogen | Kreisbogenlaenge (Least-Squares Fit) |
| BogenIndex | Zickzack / Bogen |

## Overlay-Farben

| Element | Farbe |
|---|---|
| Membranpfad (Dijkstra) | Gruen |
| Gerade Referenzlinie | Rot |
| Kreisbogen-Referenz | Cyan (nur bei Kreisbogen-Modus) |
| Endpunkte | Magenta |
| Kontaktnummern | Gelb |
