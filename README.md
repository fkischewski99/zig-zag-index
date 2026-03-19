# Zigzag-Index: ZO-1 Zell-Zell-Kontakt-Analyse

Misst den **Zickzack-Index** von ZO-1-markierten Tight Junctions in Fluoreszenzbildern.
Der Zickzack-Index ist das Verhaeltnis von tatsaechlicher Membranlaenge zu gerader Strecke zwischen zwei Tricellular Junctions (1.0 = gerade, >1.0 = gewellt/gezackt).

## Dateien

| Datei | Beschreibung |
|---|---|
| `Zigzag_Index.ijm` | ImageJ-Makro -- oeffnet LIF-Dateien, User setzt Punkte, Python misst Pfade |
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
    |       - Zickzack-Index = Pfadlaenge / Euklidische Distanz
    |
    +-- 4. Overlay-Visualisierung (gruen = Pfad, rot = gerade Strecke)
    |
    v
CSV-Export aller Ergebnisse
```

## Parameter

Im Makro (`Zigzag_Index.ijm`, Zeile 25):

| Parameter | Default | Beschreibung |
|---|---|---|
| `pythonCmd` | `"python"` | Python-Befehl |

## Python-Skript: CLI

```bash
python zigzag_pathfind.py <bild.tif> <koordinaten.csv> <ergebnisse.csv> <pfade.csv> <px_breite> <px_hoehe>
```
Berechnet den Zickzack-Index fuer jedes Koordinatenpaar via Dijkstra-Pfadfindung.

## CSV-Ausgabe

```
Bild,Kanal,Nr,Gerade (microns),Membran (microns),Zickzack-Index
"Serie1",C1,1,12.34,15.67,1.2698
"Serie1",C1,2,8.91,10.23,1.1481
"Serie1",C1,Mittelwert,,,1.2090
"Serie1",C1,Std.Abw.,,,0.0861
Gesamt,,Mittelwert,,,1.2090
Gesamt,,Std.Abw.,,,0.0861
```
