// =============================================================================
// Zigzag_Index.ijm
// Zickzack-Index von ZO-1-Zell-Zell-Kontakten — Batch-Version
//
// Verwendet Python + scikit-image fuer praezise Membranpfad-Erkennung.
// Dijkstra-Algorithmus auf Ridge-gefiltertem Bild findet den optimalen
// Pfad entlang der hellsten Membranstrukturen.
//
// Zickzack-Index = Gerade Strecke / Tatsaechliche Membranlaenge
//   1.0 = gerade Membran, <1.0 = Zickzackmuster
//
// Ablauf:
//   1. LIF-Datei oeffnen, alle Serien durchiterieren
//   2. Kanal waehlen (einmal fuer alle Serien)
//   3. Pro Serie: Auto-Vorauswahl von Junction-Punkten (Find Maxima)
//   4. User kann Punkte bearbeiten (verschieben/loeschen/hinzufuegen)
//   5. Python findet optimale Membranpfade (Dijkstra auf Ridge-Filter)
//   6. CSV-Export aller Ergebnisse in eine Datei
//
// Voraussetzungen:
//   pip install numpy scikit-image tifffile
// =============================================================================

// ---- Parameter ----
pythonCmd = "python";     // "python" oder "python3" oder voller Pfad

// Find Maxima Parameter
gaussSigma = 2;           // Glaettung vor Find Maxima
maxNoise = 20;            // Prominenz-Schwelle fuer Find Maxima
edgeMargin = 10;          // Punkte naeher am Rand werden ignoriert
minPairDist = 20;         // Minimaler Abstand fuer Paarbildung (px)
maxPairDist = 100;        // Maximaler Abstand fuer Paarbildung (px)
targetPairs = 10;         // Anzahl gewuenschter Paare (= 20 Punkte)
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
csvPath = File.saveDialog("CSV-Ergebnis speichern", "Zigzag_Results.csv");
if (csvPath == "") exit("Kein Speicherort fuer CSV gewaehlt.");

// ---- Vorbereitung ----
csvContent = "Bild,Kanal,Nr,Gerade (microns),Membran (microns),Zickzack-Index\n";
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

    // ---- Auto-Vorauswahl: Find Maxima + Paarbildung ----
    selectImage(zo1ID);
    imgW = getWidth();
    imgH = getHeight();

    // Bild duplizieren und glaetten fuer Find Maxima
    run("Duplicate...", "title=[ZZ_maxima_temp]");
    maxTempID = getImageID();
    run("Gaussian Blur...", "sigma=" + gaussSigma);

    // Find Maxima → Liste
    run("Find Maxima...", "prominence=" + maxNoise + " output=[List]");

    // Maxima aus der Results-Tabelle lesen
    nMaxima = nResults;
    maxX = newArray(nMaxima);
    maxY = newArray(nMaxima);
    maxIntensity = newArray(nMaxima);

    // Intensitaeten aus dem geglaetteten Bild lesen
    selectImage(maxTempID);
    for (i = 0; i < nMaxima; i++) {
        maxX[i] = getResult("X", i);
        maxY[i] = getResult("Y", i);
        maxIntensity[i] = getPixel(round(maxX[i]), round(maxY[i]));
    }
    run("Clear Results");

    // Maxima-Temp schliessen
    selectImage(maxTempID);
    close();

    // Maxima am Bildrand ausschliessen
    filteredX = newArray(nMaxima);
    filteredY = newArray(nMaxima);
    filteredInt = newArray(nMaxima);
    nFiltered = 0;
    for (i = 0; i < nMaxima; i++) {
        if (maxX[i] >= edgeMargin && maxX[i] < imgW - edgeMargin &&
            maxY[i] >= edgeMargin && maxY[i] < imgH - edgeMargin) {
            filteredX[nFiltered] = maxX[i];
            filteredY[nFiltered] = maxY[i];
            filteredInt[nFiltered] = maxIntensity[i];
            nFiltered++;
        }
    }

    // Nach Intensitaet sortieren (absteigend) — Bubble Sort
    for (i = 0; i < nFiltered - 1; i++) {
        for (j = 0; j < nFiltered - 1 - i; j++) {
            if (filteredInt[j] < filteredInt[j + 1]) {
                // Swap
                tmpVal = filteredInt[j]; filteredInt[j] = filteredInt[j + 1]; filteredInt[j + 1] = tmpVal;
                tmpVal = filteredX[j]; filteredX[j] = filteredX[j + 1]; filteredX[j + 1] = tmpVal;
                tmpVal = filteredY[j]; filteredY[j] = filteredY[j + 1]; filteredY[j + 1] = tmpVal;
            }
        }
    }

    // ---- Paarbildung: Naechste-Nachbarn ----
    paired = newArray(nFiltered);
    for (i = 0; i < nFiltered; i++) paired[i] = 0;

    pairX = newArray(targetPairs * 2);
    pairY = newArray(targetPairs * 2);
    nFoundPairs = 0;

    for (i = 0; i < nFiltered && nFoundPairs < targetPairs; i++) {
        if (paired[i]) continue;

        // Naechsten unverpaarten Nachbarn suchen
        bestJ = -1;
        bestDist = 1e30;
        for (j = i + 1; j < nFiltered; j++) {
            if (paired[j]) continue;
            dx = filteredX[i] - filteredX[j];
            dy = filteredY[i] - filteredY[j];
            dist = sqrt(dx * dx + dy * dy);
            if (dist >= minPairDist && dist <= maxPairDist && dist < bestDist) {
                bestDist = dist;
                bestJ = j;
            }
        }

        if (bestJ >= 0) {
            pairX[nFoundPairs * 2] = filteredX[i];
            pairY[nFoundPairs * 2] = filteredY[i];
            pairX[nFoundPairs * 2 + 1] = filteredX[bestJ];
            pairY[nFoundPairs * 2 + 1] = filteredY[bestJ];
            paired[i] = 1;
            paired[bestJ] = 1;
            nFoundPairs++;
        }
    }

    nAutoPoints = nFoundPairs * 2;
    print("[INFO] " + nFoundPairs + " Paare automatisch gefunden.");

    // ---- Vorschlag visualisieren ----
    selectImage(zo1ID);
    run("Remove Overlay");

    if (nAutoPoints > 0) {
        // Overlay-Linien zwischen Paaren als visuelle Hilfe
        for (p = 0; p < nFoundPairs; p++) {
            sx = pairX[p * 2];
            sy = pairY[p * 2];
            ex = pairX[p * 2 + 1];
            ey = pairY[p * 2 + 1];

            // Verbindungslinie (gelb, gestrichelt)
            setColor(255, 255, 0);
            Overlay.drawLine(sx, sy, ex, ey);
            Overlay.setPosition(0);

            // Paar-Nummer
            setColor(255, 255, 0);
            Overlay.drawString("" + (p + 1), (sx + ex) / 2 + 5, (sy + ey) / 2 - 5);
            Overlay.setPosition(0);
        }
        Overlay.show();

        // Multi-Point Selection erstellen
        // Arrays auf exakte Groesse trimmen
        selX = newArray(nAutoPoints);
        selY = newArray(nAutoPoints);
        for (i = 0; i < nAutoPoints; i++) {
            selX[i] = pairX[i];
            selY[i] = pairY[i];
        }
        makeSelection("point", selX, selY);
    }

    setTool("multipoint");

    waitForUser("Serie: " + seriesName + " — Kontakte pruefen",
        "Automatisch " + nFoundPairs + " Paar(e) vorgeschlagen.\n\n" +
        "Punkte pruefen und ggf. bearbeiten:\n" +
        "  Punkt 1 + 2 = erster Kontakt\n" +
        "  Punkt 3 + 4 = zweiter Kontakt\n" +
        "  ... usw.\n\n" +
        "Punkte auf trizellulare Junctions setzen.\n\n" +
        "Tipps:\n" +
        "- Zoomen mit +/- oder Mausrad\n" +
        "- Punkt verschieben: Ziehen\n" +
        "- Punkt entfernen: Alt+Klick\n" +
        "- Neuen Punkt: Klick\n" +
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
        csvContent += "\"" + seriesName + "\"," + chosenChannel + "," + (p + 1) + "," +
                      d2s(pairEuclid[p], 2) + "," + d2s(pairLength[p], 2) + "," +
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

    csvContent += "\"" + seriesName + "\"," + chosenChannel + ",Mittelwert,,," + d2s(seriesMean, 4) + "\n";
    csvContent += "\"" + seriesName + "\"," + chosenChannel + ",Std.Abw.,,," + d2s(seriesSD, 4) + "\n";

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

    csvContent += "Gesamt,,Mittelwert,,," + d2s(totalMean, 4) + "\n";
    csvContent += "Gesamt,,Std.Abw.,,," + d2s(totalSD, 4) + "\n";

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

showMessage("Analyse abgeschlossen",
    "Alle " + seriesCount + " Serien verarbeitet.\n\n" +
    "Gesamt gemessene Kontakte: " + totalMeasured + "\n" +
    (totalMeasured > 0 ? "Mittelwert Zigzag-Index: " + d2s(totalMean, 4) + "\n" +
    "Standardabweichung: " + d2s(totalSD, 4) + "\n\n" : "\n") +
    "CSV gespeichert unter:\n" + csvPath);
