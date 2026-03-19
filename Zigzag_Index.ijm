// =============================================================================
// Zigzag_Index.ijm
// Zickzack-Index von ZO-1-Zell-Zell-Kontakten — Batch-Version
//
// Verwendet Python + scikit-image fuer praezise Membranpfad-Erkennung.
// Dijkstra-Algorithmus auf Ridge-gefiltertem Bild findet den optimalen
// Pfad entlang der hellsten Membranstrukturen.
//
// Zickzack-Index = Tatsaechliche Membranlaenge / Gerade Strecke
//   1.0 = gerade Membran, >1.0 = Zickzackmuster
//
// Ablauf:
//   1. LIF-Datei oeffnen, alle Serien durchiterieren
//   2. Kanal waehlen (einmal fuer alle Serien)
//   3. Pro Serie: User setzt Punkte paarweise auf trizellulare Junctions
//   5. Python findet optimale Membranpfade (Dijkstra auf Ridge-Filter)
//   6. CSV-Export aller Ergebnisse in eine Datei
//
// Voraussetzungen:
//   pip install numpy scikit-image
// =============================================================================

// ---- Parameter ----
pythonCmd = "python";     // "python" oder "python3" oder voller Pfad

// ---- Ende Parameter ----

requires("1.53");

// ---- Python-Skript finden ----
macroDir = getDirectory("macros");
scriptPath = macroDir + "zigzag_pathfind.py";
if (!File.exists(scriptPath)) {
    scriptPath = getDirectory("imagej") + "zigzag_pathfind.py";
}
if (!File.exists(scriptPath)) {
    scriptPath = File.openDialog("zigzag_pathfind.py auswaehlen");
    if (scriptPath == "") exit("Kein Python-Skript ausgewaehlt.");
}

// ---- LIF-Datei auswaehlen ----
lifPath = File.openDialog("LIF-Datei auswaehlen");
if (lifPath == "") exit("Keine Datei ausgewaehlt.");

// ---- Bio-Formats: Serien zaehlen ----
run("Bio-Formats Macro Extensions");
Ext.setId(lifPath);
Ext.getSeriesCount(seriesCount);

if (seriesCount == 0) exit("Keine Serien in der LIF-Datei gefunden.");
print("[INFO] LIF-Datei: " + lifPath);
print("[INFO] " + seriesCount + " Serie(n) gefunden.");

// ---- Kanal waehlen (einmal fuer alle Serien) ----
// Erste Serie oeffnen um Kanalinformationen zu lesen
Ext.setSeries(0);
Ext.getSizeC(nChannels);

chosenChannel = "";
chosenIdx = 1;
if (nChannels > 1) {
    channelNames = newArray(nChannels);
    for (c = 0; c < nChannels; c++)
        channelNames[c] = "C" + (c + 1);
    Dialog.create("ZO-1 Kanal auswaehlen");
    Dialog.addChoice("ZO-1 Kanal (gilt fuer alle Serien):", channelNames, channelNames[0]);
    Dialog.show();
    chosenChannel = Dialog.getChoice();
    chosenIdx = parseInt(substring(chosenChannel, 1));
} else {
    chosenChannel = "C1";
}

// ---- CSV-Speicherort waehlen ----
Dialog.create("CSV-Ergebnis speichern");
Dialog.addDirectory("Speicherort:", getDirectory("home"));
Dialog.addString("Dateiname:", "Zigzag_Results.csv", 30);
Dialog.show();
csvDir = Dialog.getString();
csvName = Dialog.getString();
csvPath = csvDir + csvName;

// ---- Vorbereitung ----
csvContent = "Name,Kanal,Zickzack,Gerade,Index\n";
tmpDir = getDirectory("temp");

// Gesamt-Statistik
totalSumZZ = 0;
totalSumZZ2 = 0;
totalMeasured = 0;

// ---- Hauptschleife: Alle Serien durchgehen ----
for (s = 0; s < seriesCount; s++) {
    // Serienname holen
    Ext.setSeries(s);
    Ext.getSeriesName(seriesName);
    print("================================================");
    print("[INFO] Serie " + (s + 1) + "/" + seriesCount + ": " + seriesName);

    // Serie oeffnen
    seriesStr = "";
    for (si = 0; si < seriesCount; si++) {
        if (si == s)
            seriesStr += " series_" + (si + 1);
    }
    run("Bio-Formats Importer", "open=[" + lifPath + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT" + seriesStr);

    if (nImages == 0) {
        print("[WARNUNG] Serie " + seriesName + " konnte nicht geoeffnet werden. Ueberspringe.");
        continue;
    }

    origID = getImageID();

    // ---- Z-Projection (Maximum Intensity) ----
    selectImage(origID);
    if (nSlices > 1) {
        run("Z Project...", "projection=[Max Intensity]");
        zprojID = getImageID();
        selectImage(origID);
        close();
    } else {
        zprojID = origID;
    }

    // ---- Channel Split + Auswahl ----
    selectImage(zprojID);
    getDimensions(w, h, channels, slices, frames);
    zo1ID = zprojID;

    if (channels > 1) {
        selectImage(zprojID);
        run("Split Channels");
        list = getList("image.titles");
        zo1ID = 0;
        for (i = 0; i < list.length; i++) {
            if (startsWith(list[i], "C" + chosenIdx + "-")) {
                selectWindow(list[i]);
                zo1ID = getImageID();
            } else {
                selectWindow(list[i]);
                close();
            }
        }
        if (zo1ID == 0) {
            print("[WARNUNG] Kanal " + chosenChannel + " nicht gefunden in Serie " + seriesName + ". Ueberspringe.");
            continue;
        }
    }

    // Pixelkalibrierung auslesen
    selectImage(zo1ID);
    getPixelSize(unit, pixelWidth, pixelHeight);
    if (unit == "pixels" || unit == "pixel") {
        unit = "px";
    }

    // ---- Manuelle Punktauswahl ----
    selectImage(zo1ID);
    setTool("multipoint");

    waitForUser("Serie: " + seriesName + " — Kontakte markieren",
        "Punkte paarweise auf trizellulare Junctions setzen:\n" +
        "  Punkt 1 + 2 = erster Kontakt\n" +
        "  Punkt 3 + 4 = zweiter Kontakt\n" +
        "  ... usw.\n\n" +
        "Tipps:\n" +
        "- Zoomen mit +/- oder Mausrad\n" +
        "- Punkt verschieben: Ziehen\n" +
        "- Punkt entfernen: Alt+Klick\n" +
        "- Gerade Anzahl an Punkten!\n\n" +
        "Wenn fertig: OK druecken.");

    // ---- Punkte auslesen ----
    if (selectionType() != 10) {
        print("[WARNUNG] Keine Punkte in Serie " + seriesName + ". Ueberspringe.");
        // Aufraeumen
        run("Remove Overlay");
        selectImage(zo1ID);
        close();
        continue;
    }

    getSelectionCoordinates(clickX, clickY);
    nClicks = clickX.length;
    run("Select None");
    run("Remove Overlay");

    if (nClicks < 2) {
        print("[WARNUNG] Weniger als 2 Punkte in Serie " + seriesName + ". Ueberspringe.");
        selectImage(zo1ID);
        close();
        continue;
    }

    if (nClicks % 2 != 0) {
        showMessage("Hinweis", "Ungerade Anzahl Punkte (" + nClicks + ") in Serie " + seriesName + ".\nDer letzte Punkt wird ignoriert.");
        nClicks = nClicks - 1;
    }

    nPairs = nClicks / 2;
    print("[INFO] " + nPairs + " Kontakt(e) markiert.");

    // ---- Bild + Koordinaten fuer Python speichern ----
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

    // ---- Python-Skript ausfuehren ----
    print("[INFO] Starte Python-Analyse fuer " + seriesName + "...");
    pyCall = pythonCmd + " \"" + scriptPath + "\" \"" + tmpImage + "\" \"" +
             tmpCoords + "\" \"" + tmpResults + "\" \"" + tmpPaths + "\" " +
             d2s(pixelWidth, 6) + " " + d2s(pixelHeight, 6);

    pyOutput = exec(pyCall);
    print("[Python] " + pyOutput);

    // Pruefen ob Ergebnisse existieren
    if (!File.exists(tmpResults)) {
        print("[WARNUNG] Python fehlgeschlagen fuer Serie " + seriesName + ".");
        print("[Python-Kommando] " + pyCall);
        print("[Python-Ausgabe] " + pyOutput);
        // Aufraeumen und weiter
        selectImage(zo1ID);
        close();
        continue;
    }

    // ---- Ergebnisse einlesen ----
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
        print("[WARNUNG] Kein Kontakt konnte gemessen werden in Serie " + seriesName + ".");
        selectImage(zo1ID);
        close();
        // Temp-Dateien loeschen
        if (File.exists(tmpImage)) File.delete(tmpImage);
        if (File.exists(tmpCoords)) File.delete(tmpCoords);
        if (File.exists(tmpResults)) File.delete(tmpResults);
        if (File.exists(tmpPaths)) File.delete(tmpPaths);
        continue;
    }

    // ---- Pfade einlesen und Overlays zeichnen ----
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

    // ---- Ergebnisse an CSV anhaengen ----
    seriesSumZZ = 0;
    seriesSumZZ2 = 0;

    for (p = 0; p < nPairs; p++) {
        if (pairLength[p] <= 0) continue;
        csvContent += "\"" + seriesName + "\"," + chosenChannel + "," +
                      d2s(pairLength[p], 2) + "," + d2s(pairEuclid[p], 2) + "," +
                      d2s(pairZigzag[p], 4) + "\n";
        seriesSumZZ += pairZigzag[p];
        seriesSumZZ2 += pairZigzag[p] * pairZigzag[p];
    }

    // Serien-Zusammenfassung
    seriesMean = seriesSumZZ / nMeasured;
    if (nMeasured > 1) {
        seriesSD = sqrt((seriesSumZZ2 - nMeasured * seriesMean * seriesMean) / (nMeasured - 1));
    } else {
        seriesSD = 0;
    }

    csvContent += "\n";
    csvContent += ",,Mittelwert,," + d2s(seriesMean, 4) + "\n";
    csvContent += "\n\n";

    // Gesamt-Statistik aktualisieren
    totalSumZZ += seriesSumZZ;
    totalSumZZ2 += seriesSumZZ2;
    totalMeasured += nMeasured;

    print("[INFO] Serie " + seriesName + ": " + nMeasured + " Kontakte, Mittelwert=" + d2s(seriesMean, 4) + ", SD=" + d2s(seriesSD, 4));

    // ---- User Ergebnis zeigen, dann naechste Serie ----
    if (s < seriesCount - 1) {
        waitForUser("Serie " + seriesName + " fertig",
            "Ergebnisse fuer " + seriesName + ":\n\n" +
            "Gemessene Kontakte: " + nMeasured + "\n" +
            "Mittelwert Zigzag-Index: " + d2s(seriesMean, 4) + "\n" +
            "Standardabweichung: " + d2s(seriesSD, 4) + "\n\n" +
            "OK druecken fuer naechste Serie.");
    }

    // ---- Aufraeumen: Bilder dieser Serie schliessen ----
    selectImage(zo1ID);
    close();

    // Temp-Dateien loeschen
    if (File.exists(tmpImage)) File.delete(tmpImage);
    if (File.exists(tmpCoords)) File.delete(tmpCoords);
    if (File.exists(tmpResults)) File.delete(tmpResults);
    if (File.exists(tmpPaths)) File.delete(tmpPaths);
}

// ---- Gesamt-Zusammenfassung an CSV anhaengen ----
if (totalMeasured > 0) {
    totalMean = totalSumZZ / totalMeasured;
    if (totalMeasured > 1) {
        totalSD = sqrt((totalSumZZ2 - totalMeasured * totalMean * totalMean) / (totalMeasured - 1));
    } else {
        totalSD = 0;
    }

    // ---- CSV schreiben ----
    File.saveString(csvContent, csvPath);
    print("================================================");
    print("CSV gespeichert: " + csvPath);
}

// ---- Zusammenfassung im Log ----
print("================================================");
print("Zigzag-Index Analyse abgeschlossen");
print("Verarbeitete Serien: " + seriesCount);
print("Gesamt gemessene Kontakte: " + totalMeasured);
if (totalMeasured > 0) {
    print("Gesamt Mittelwert Zigzag-Index: " + d2s(totalMean, 4));
    print("Gesamt Standardabweichung: " + d2s(totalSD, 4));
}
print("CSV: " + csvPath);
print("================================================");

summaryMsg = "Alle " + seriesCount + " Serien verarbeitet.\n\n" +
    "Gesamt gemessene Kontakte: " + totalMeasured + "\n";
if (totalMeasured > 0) {
    summaryMsg += "Mittelwert Zigzag-Index: " + d2s(totalMean, 4) + "\n" +
        "Standardabweichung: " + d2s(totalSD, 4) + "\n";
}
summaryMsg += "\nCSV gespeichert unter:\n" + csvPath;
showMessage("Analyse abgeschlossen", summaryMsg);
