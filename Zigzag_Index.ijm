// =============================================================================
// Zigzag_Index.ijm
// Zickzack-Index von ZO-1-Zell-Zell-Kontakten
//
// Verwendet Python + scikit-image fuer praezise Membranpfad-Erkennung.
// Dijkstra-Algorithmus auf Ridge-gefiltertem Bild findet den optimalen
// Pfad entlang der hellsten Membranstrukturen.
//
// Zickzack-Index = Gerade Strecke / Tatsaechliche Membranlaenge
//   1.0 = gerade Membran, <1.0 = Zickzackmuster
//
// Ablauf:
//   1. Bild oeffnen, Z-Projection, Kanal waehlen
//   2. User klickt Punktpaare (Start+Ende jedes Kontakts)
//   3. Python findet optimale Membranpfade (Dijkstra auf Ridge-Filter)
//   4. Zigzag-Index berechnen + Overlays
//
// Voraussetzungen:
//   pip install numpy scikit-image tifffile
// =============================================================================

// ---- Parameter ----
pythonCmd = "python";     // "python" oder "python3" oder voller Pfad
// ---- Ende Parameter ----

requires("1.53");

// Pfad zum Python-Skript (im gleichen Ordner wie dieses Makro)
macroDir = getDirectory("macros");
// Falls das Makro nicht aus dem Macros-Ordner laeuft, Skript im selben Verzeichnis suchen
scriptPath = macroDir + "zigzag_pathfind.py";
if (!File.exists(scriptPath)) {
    // Skript neben der .ijm Datei suchen
    scriptPath = getDirectory("imagej") + "zigzag_pathfind.py";
}
if (!File.exists(scriptPath)) {
    // Dialog: User soll Pfad angeben
    scriptPath = File.openDialog("zigzag_pathfind.py auswaehlen");
    if (scriptPath == "") exit("Kein Python-Skript ausgewaehlt.");
}

// ---- LIF-Datei oeffnen (Bio-Formats) ----
if (nImages == 0) {
    lifPath = File.openDialog("LIF-Datei auswaehlen");
    if (lifPath == "") exit("Keine Datei ausgewaehlt.");
    run("Bio-Formats Importer", "open=[" + lifPath + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
    if (nImages == 0)
        exit("Datei konnte nicht geoeffnet werden.");
}

origID = getImageID();
origTitle = getTitle();
chosenChannel = "";

// ---- Z-Projection (Maximum Intensity) ----
selectImage(origID);
if (nSlices > 1) {
    run("Z Project...", "projection=[Max Intensity]");
    zprojID = getImageID();
} else {
    run("Duplicate...", "title=[ZZ_zprojected]");
    zprojID = getImageID();
}

// ---- Channel Split + Auswahl ----
selectImage(zprojID);
getDimensions(w, h, channels, slices, frames);
if (channels > 1) {
    channelNames = newArray(channels);
    for (c = 0; c < channels; c++)
        channelNames[c] = "C" + (c + 1);
    Dialog.create("ZO-1 Kanal auswaehlen");
    Dialog.addChoice("ZO-1 Kanal:", channelNames, channelNames[0]);
    Dialog.show();
    chosenChannel = Dialog.getChoice();
    chosenIdx = parseInt(substring(chosenChannel, 1));

    selectImage(zprojID);
    run("Split Channels");
    list = getList("image.titles");
    zo1ID = 0;
    for (i = 0; i < list.length; i++) {
        if (startsWith(list[i], "C" + chosenIdx + "-")) {
            selectWindow(list[i]);
            zo1ID = getImageID();
        } else {
            for (c = 0; c < channels; c++) {
                if (startsWith(list[i], "C" + (c + 1) + "-") && (c + 1) != chosenIdx) {
                    selectWindow(list[i]);
                    close();
                }
            }
        }
    }
} else {
    zo1ID = zprojID;
}

// Pixelkalibrierung auslesen
selectImage(zo1ID);
getPixelSize(unit, pixelWidth, pixelHeight);
if (unit == "pixels" || unit == "pixel") {
    showMessage("Hinweis", "Bild ist nicht kalibriert. Ergebnisse werden in Pixeln angegeben.");
    unit = "px";
}

// ---- 1. User klickt Punktpaare ----
selectImage(zo1ID);
run("Select None");
setTool("multipoint");

waitForUser("Kontakte markieren (Punktpaare)",
    "Klicken Sie mit dem Multi-Point Tool auf die\n" +
    "Start- und Endpunkte der Zell-Zell-Kontakte:\n\n" +
    "  Punkt 1 + 2 = erster Kontakt\n" +
    "  Punkt 3 + 4 = zweiter Kontakt\n" +
    "  Punkt 5 + 6 = dritter Kontakt\n" +
    "  ... usw.\n\n" +
    "Setzen Sie die Punkte auf die trizellularen\n" +
    "Junctions an beiden Enden eines Kontakts.\n\n" +
    "Tipps:\n" +
    "- Zoomen mit +/- oder Mausrad\n" +
    "- Punkt entfernen: Alt+Klick\n" +
    "- Gerade Anzahl an Punkten setzen!\n\n" +
    "Wenn fertig: OK druecken.");

if (selectionType() != 10) {
    exit("Keine Punkte markiert. Bitte Makro erneut starten.");
}

getSelectionCoordinates(clickX, clickY);
nClicks = clickX.length;
run("Select None");

if (nClicks < 2) {
    exit("Mindestens 2 Punkte (1 Paar) noetig.");
}

if (nClicks % 2 != 0) {
    showMessage("Hinweis", "Ungerade Anzahl Punkte (" + nClicks + ").\nDer letzte Punkt wird ignoriert.");
    nClicks = nClicks - 1;
}

nPairs = nClicks / 2;
print("[INFO] " + nPairs + " Kontakt(e) markiert.");

// ---- 2. Bild + Koordinaten fuer Python speichern ----
tmpDir = getDirectory("temp");
tmpImage = tmpDir + "zz_input.tif";
tmpCoords = tmpDir + "zz_coords.csv";
tmpResults = tmpDir + "zz_results.csv";
tmpPaths = tmpDir + "zz_paths.csv";

// Bild als 8-bit Graustufen speichern
selectImage(zo1ID);
run("Duplicate...", "title=[ZZ_export]");
exportID = getImageID();
if (bitDepth() != 8) run("8-bit");
run("Grays");
saveAs("Tiff", tmpImage);
close();

// Koordinaten speichern
coordStr = "";
for (i = 0; i < nClicks; i++) {
    coordStr += d2s(clickX[i], 2) + "," + d2s(clickY[i], 2) + "\n";
}
File.saveString(coordStr, tmpCoords);

// ---- 3. Python-Skript ausfuehren ----
print("[INFO] Starte Python-Analyse...");

// Anfuehrungszeichen fuer Pfade mit Leerzeichen
pyCall = pythonCmd + " \"" + scriptPath + "\" \"" + tmpImage + "\" \"" +
         tmpCoords + "\" \"" + tmpResults + "\" \"" + tmpPaths + "\" " +
         d2s(pixelWidth, 6) + " " + d2s(pixelHeight, 6);

pyOutput = exec(pyCall);
print("[Python] " + pyOutput);

// Pruefen ob Ergebnisse existieren
if (!File.exists(tmpResults)) {
    exit("Python-Skript fehlgeschlagen.\n\n" +
         "Pruefen Sie:\n" +
         "1. Python installiert? ('" + pythonCmd + "' im PATH?)\n" +
         "2. scikit-image installiert? (pip install scikit-image)\n" +
         "3. tifffile installiert? (pip install tifffile)\n\n" +
         "Kommando: " + pyCall + "\n" +
         "Ausgabe: " + pyOutput);
}

// ---- 4. Ergebnisse einlesen ----
resultsStr = File.openAsString(tmpResults);
resultsLines = split(resultsStr, "\n");

pairNr = newArray(nPairs);
pairEuclid = newArray(nPairs);
pairLength = newArray(nPairs);
pairZigzag = newArray(nPairs);
nMeasured = 0;

for (i = 0; i < resultsLines.length; i++) {
    if (resultsLines[i] == "") continue;
    cols = split(resultsLines[i], ",");
    if (cols.length < 4) continue;
    idx = parseInt(cols[0]) - 1;
    if (idx < 0 || idx >= nPairs) continue;
    pairNr[idx] = parseInt(cols[0]);
    pairEuclid[idx] = parseFloat(cols[1]);
    pairLength[idx] = parseFloat(cols[2]);
    pairZigzag[idx] = parseFloat(cols[3]);
    if (pairLength[idx] > 0) nMeasured++;
}

if (nMeasured == 0) {
    exit("Kein Kontakt konnte gemessen werden.\n" +
         "Python-Ausgabe: " + pyOutput);
}

// ---- 5. Pfade einlesen und Overlays zeichnen ----
selectImage(zo1ID);
run("Remove Overlay");

pathsStr = File.openAsString(tmpPaths);
pathsLines = split(pathsStr, "\n");

// Pfade (gruen) zeichnen
for (i = 0; i < pathsLines.length - 1; i++) {
    if (pathsLines[i] == "" || pathsLines[i + 1] == "") continue;
    cols1 = split(pathsLines[i], ",");
    cols2 = split(pathsLines[i + 1], ",");
    if (cols1.length < 3 || cols2.length < 3) continue;

    // Nur aufeinanderfolgende Punkte des gleichen Paars verbinden
    nr1 = parseInt(cols1[0]);
    nr2 = parseInt(cols2[0]);
    if (nr1 != nr2) continue;

    px1 = parseFloat(cols1[1]);
    py1 = parseFloat(cols1[2]);
    px2 = parseFloat(cols2[1]);
    py2 = parseFloat(cols2[2]);

    setColor(0, 255, 0);
    Overlay.drawLine(px1, py1, px2, py2);
    Overlay.setPosition(0);
}

// Gerade Strecken (rot) + Endpunkte (magenta) + Nummern (gelb)
for (p = 0; p < nPairs; p++) {
    if (pairLength[p] <= 0) continue;

    sx = clickX[p * 2];
    sy = clickY[p * 2];
    ex = clickX[p * 2 + 1];
    ey = clickY[p * 2 + 1];

    // Gerade Strecke (rot)
    setColor(255, 0, 0);
    Overlay.drawLine(sx, sy, ex, ey);
    Overlay.setPosition(0);

    // Endpunkte (magenta)
    r = 4;
    setColor(255, 0, 255);
    Overlay.drawEllipse(sx - r, sy - r, 2*r, 2*r);
    Overlay.setPosition(0);
    Overlay.drawEllipse(ex - r, ey - r, 2*r, 2*r);
    Overlay.setPosition(0);

    // Nummer (gelb)
    setColor(255, 255, 0);
    Overlay.drawString("" + (p + 1), (sx + ex) / 2 + 5, (sy + ey) / 2 - 5);
    Overlay.setPosition(0);
}

Overlay.show();

// ---- Results-Tabelle ----
tableName = "Zigzag Index Results";
Table.create(tableName);

channelInfo = "alle";
if (chosenChannel != "") channelInfo = chosenChannel;

sumZZ = 0;
sumZZ2 = 0;
row = 0;

for (p = 0; p < nPairs; p++) {
    if (pairLength[p] <= 0) continue;
    Table.set("Bild", row, origTitle);
    Table.set("Kanal", row, channelInfo);
    Table.set("Nr", row, p + 1);
    Table.set("Gerade (" + unit + ")", row, pairEuclid[p]);
    Table.set("Membran (" + unit + ")", row, pairLength[p]);
    Table.set("Zickzack-Index", row, pairZigzag[p]);
    sumZZ += pairZigzag[p];
    sumZZ2 += pairZigzag[p] * pairZigzag[p];
    row++;
}

Table.update();

// ---- Zusammenfassung ----
meanZZ = sumZZ / nMeasured;
if (nMeasured > 1) {
    sdZZ = sqrt((sumZZ2 - nMeasured * meanZZ * meanZZ) / (nMeasured - 1));
} else {
    sdZZ = 0;
}

r0 = row;
Table.set("Bild", r0, "");
Table.set("Kanal", r0, "");
Table.set("Nr", r0, "");
Table.set("Gerade (" + unit + ")", r0, "");
Table.set("Membran (" + unit + ")", r0, "");
Table.set("Zickzack-Index", r0, "");

Table.set("Bild", r0 + 1, "");
Table.set("Kanal", r0 + 1, "");
Table.set("Nr", r0 + 1, "Mittelwert");
Table.set("Gerade (" + unit + ")", r0 + 1, "");
Table.set("Membran (" + unit + ")", r0 + 1, "");
Table.set("Zickzack-Index", r0 + 1, meanZZ);

Table.set("Bild", r0 + 2, "");
Table.set("Kanal", r0 + 2, "");
Table.set("Nr", r0 + 2, "Std.Abw.");
Table.set("Gerade (" + unit + ")", r0 + 2, "");
Table.set("Membran (" + unit + ")", r0 + 2, "");
Table.set("Zickzack-Index", r0 + 2, sdZZ);

Table.update();

// ---- Aufraeumen ----
// Temp-Dateien loeschen
File.delete(tmpImage);
File.delete(tmpCoords);
File.delete(tmpResults);
File.delete(tmpPaths);

selectImage(zo1ID);

print("================================================");
print("Zigzag-Index Analyse abgeschlossen");
print("Bild: " + origTitle);
print("Kanal: " + channelInfo);
print("Gemessene Kontakte: " + nMeasured);
print("Mittelwert Zigzag-Index: " + d2s(meanZZ, 4));
print("Standardabweichung: " + d2s(sdZZ, 4));
print("================================================");
