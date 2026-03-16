// =============================================================================
// Zigzag_Index_Auto.ijm
// Automatischer Zickzack-Index von ZO-1-Zell-Zell-Kontakten
//
// Wie Zigzag_Index.ijm, aber Tight Junctions werden automatisch erkannt.
// Keine manuelle Punktmarkierung noetig.
//
// Erkennung:
//   - AnalyzeSkeleton findet Branch Points (trizellulare Junctions)
//   - Branches zwischen zwei Branch Points = Zell-Zell-Kontakte
//   - Sackgassen (Junction-zu-Endpoint) und kurze Aeste werden gefiltert
//
// Zickzack-Index = Gerade Strecke / Tatsaechliche Membranlaenge
//   1.0 = gerade Membran, <1.0 = Zickzackmuster
//
// Ablauf:
//   1. Bild oeffnen, Z-Projection, Kanal automatisch waehlen
//   2. Preprocessing + Skeletonisierung (identisch zum manuellen Skript)
//   3. AnalyzeSkeleton: Branch Points + Branches automatisch erkennen
//   4. Junction-zu-Junction Branches filtern (= Zell-Zell-Kontakte)
//   5. Skeleton-Pfad tracen, Zigzag-Index berechnen, Overlays zeichnen
// =============================================================================

// ---- Parameter (anpassbar) ----
gaussSigma = 1.0;
medianRadius = 1;
rollingBallRadius = 50;
thresholdMethod = "Otsu";
minContactLength = 10;    // Min. Kontaktlaenge (kalibrierte Einheiten, bzw. Pixel falls unkalibriert)
zo1Channel = 0;           // 0 = erster/einziger Kanal, 1-N = bestimmter Kanal
snapRadius = 5;           // Snap-Radius (Pixel) fuer Endpunkte auf Skeleton
// ---- Ende Parameter ----

requires("1.53");

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

// ---- Z-Projection (Maximum Intensity) ----
selectImage(origID);
if (nSlices > 1) {
    run("Z Project...", "projection=[Max Intensity]");
    zprojID = getImageID();
} else {
    run("Duplicate...", "title=[ZZ_zprojected]");
    zprojID = getImageID();
}

// ---- Channel Auswahl (automatisch per Parameter, kein Dialog) ----
selectImage(zprojID);
getDimensions(w, h, channels, slices, frames);
channelInfo = "alle";
if (channels > 1) {
    if (zo1Channel == 0) zo1Channel = 1;
    channelInfo = "C" + zo1Channel;
    selectImage(zprojID);
    run("Split Channels");
    list = getList("image.titles");
    zo1ID = 0;
    for (i = 0; i < list.length; i++) {
        if (startsWith(list[i], "C" + zo1Channel + "-")) {
            selectWindow(list[i]);
            zo1ID = getImageID();
        } else {
            for (c = 0; c < channels; c++) {
                if (startsWith(list[i], "C" + (c + 1) + "-") && (c + 1) != zo1Channel) {
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
    unit = "px";
}

// ---- 1. Preprocessing ----
selectImage(zo1ID);
run("Select None");
run("Duplicate...", "title=[ZZ_work]");
workID = getImageID();

if (bitDepth() != 8) run("8-bit");
run("Grays");
run("Median...", "radius=" + medianRadius);
run("Gaussian Blur...", "sigma=" + gaussSigma);
run("Subtract Background...", "rolling=" + rollingBallRadius);
run("Enhance Contrast...", "saturated=0.35 normalize");

// ---- 2. Segmentierung ----
setAutoThreshold(thresholdMethod + " dark");
run("Convert to Mask");
run("Close-");
run("Close-");
run("Open");
run("Open");

// ---- 3. Skeletonisierung + AnalyzeSkeleton ----
run("Skeletonize");
skelID = getImageID();

if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
if (isOpen("Branch information")) { selectWindow("Branch information"); run("Close"); }
selectImage(skelID);
run("Analyze Skeleton (2D/3D)", "prune=[shortest branch] prune cycle method=[shortest branch]");

// Tagged skeleton finden (Pixelwerte: 30=Endpoint, 70=Junction, 127=Slab)
taggedID = 0;
list = getList("image.titles");
for (i = 0; i < list.length; i++) {
    if (indexOf(list[i], "Tagged skeleton") >= 0) {
        selectWindow(list[i]);
        taggedID = getImageID();
    }
}

// ---- 4. Branch-Information auswerten ----
if (!isOpen("Branch information")) {
    if (taggedID > 0) { selectImage(taggedID); close(); }
    closeImagesStartingWith("ZZ_");
    closeWindowIfOpen("Results");
    closeWindowIfOpen("Longest shortest paths");
    exit("Keine Branch-Information gefunden.\nKein zusammenhaengendes Skeleton vorhanden.");
}

selectWindow("Branch information");
nBranches = Table.size;

if (nBranches == 0) {
    closeWindowIfOpen("Branch information");
    if (taggedID > 0) { selectImage(taggedID); close(); }
    closeImagesStartingWith("ZZ_");
    closeWindowIfOpen("Results");
    closeWindowIfOpen("Longest shortest paths");
    exit("Keine Branches im Skeleton gefunden.");
}

// Branch-Daten einlesen
brV1x = newArray(nBranches);
brV1y = newArray(nBranches);
brV2x = newArray(nBranches);
brV2y = newArray(nBranches);
brLength = newArray(nBranches);
brEuclid = newArray(nBranches);

selectWindow("Branch information");
for (i = 0; i < nBranches; i++) {
    brV1x[i] = Table.get("V1 x", i);
    brV1y[i] = Table.get("V1 y", i);
    brV2x[i] = Table.get("V2 x", i);
    brV2y[i] = Table.get("V2 y", i);
    brLength[i] = Table.get("Branch length", i);
    brEuclid[i] = Table.get("Euclidean distance", i);
}

// AnalyzeSkeleton-Fenster schliessen
closeWindowIfOpen("Branch information");
closeWindowIfOpen("Results");
closeWindowIfOpen("Longest shortest paths");

// ---- 5. Junction-zu-Junction Branches filtern ----
// Ein Branch ist ein Zell-Zell-Kontakt wenn beide Endpunkte (V1, V2)
// auf Junction-Pixeln liegen (Wert 70 im Tagged Skeleton)
isContact = newArray(nBranches);
nContacts = 0;

if (taggedID > 0) {
    selectImage(taggedID);
    tagW = getWidth();
    tagH = getHeight();

    for (i = 0; i < nBranches; i++) {
        v1x = round(brV1x[i]);
        v1y = round(brV1y[i]);
        v2x = round(brV2x[i]);
        v2y = round(brV2y[i]);

        // V1: Junction? (Wert 70, mit 1px Toleranz fuer Cluster-Raender)
        v1Junction = false;
        selectImage(taggedID);
        for (dy = -1; dy <= 1 && !v1Junction; dy++) {
            for (dx = -1; dx <= 1 && !v1Junction; dx++) {
                nx = v1x + dx;
                ny = v1y + dy;
                if (nx >= 0 && nx < tagW && ny >= 0 && ny < tagH) {
                    if (getPixel(nx, ny) == 70) v1Junction = true;
                }
            }
        }

        // V2: Junction?
        v2Junction = false;
        selectImage(taggedID);
        for (dy = -1; dy <= 1 && !v2Junction; dy++) {
            for (dx = -1; dx <= 1 && !v2Junction; dx++) {
                nx = v2x + dx;
                ny = v2y + dy;
                if (nx >= 0 && nx < tagW && ny >= 0 && ny < tagH) {
                    if (getPixel(nx, ny) == 70) v2Junction = true;
                }
            }
        }

        // Nur Junction-zu-Junction UND ueber Mindestlaenge
        if (v1Junction && v2Junction && brLength[i] >= minContactLength) {
            isContact[i] = 1;
            nContacts++;
        }
    }

    selectImage(taggedID);
    close();
} else {
    // Fallback ohne Tagged Skeleton: nur Laengenfilter
    print("[WARNUNG] Kein Tagged Skeleton. Alle Branches ueber Mindestlaenge werden verwendet.");
    for (i = 0; i < nBranches; i++) {
        if (brLength[i] >= minContactLength) {
            isContact[i] = 1;
            nContacts++;
        }
    }
}

if (nContacts == 0) {
    closeImagesStartingWith("ZZ_");
    exit("Keine Zell-Zell-Kontakte gefunden.\n" +
         nBranches + " Branches erkannt, keiner erfuellt Junction-zu-Junction + Mindestlaenge (" + minContactLength + ").\n" +
         "Parameter minContactLength anpassen oder Preprocessing ueberpruefen.");
}

print("[INFO] " + nContacts + " Zell-Zell-Kontakte automatisch erkannt (von " + nBranches + " Branches).");

// ---- 6. Pfade tracen, Zigzag-Index berechnen, Overlays zeichnen ----
selectImage(zo1ID);
run("Remove Overlay");

pairLengths = newArray(nContacts);
pairEuclid  = newArray(nContacts);
pairZigzag  = newArray(nContacts);
pairStartX  = newArray(nContacts);
pairStartY  = newArray(nContacts);
pairEndX    = newArray(nContacts);
pairEndY    = newArray(nContacts);
nMeasured = 0;
ci = 0;

for (i = 0; i < nBranches; i++) {
    if (isContact[i] == 0) continue;

    sx = round(brV1x[i]);
    sy = round(brV1y[i]);
    ex = round(brV2x[i]);
    ey = round(brV2y[i]);

    // Snap auf naechsten Skeleton-Pixel
    selectImage(skelID);
    bestSx = sx; bestSy = sy; bestSD = 999999;
    for (dy = -snapRadius; dy <= snapRadius; dy++) {
        for (dx = -snapRadius; dx <= snapRadius; dx++) {
            nx = sx + dx; ny = sy + dy;
            if (nx >= 0 && nx < getWidth() && ny >= 0 && ny < getHeight()) {
                if (getPixel(nx, ny) > 0) {
                    d = sqrt(dx*dx + dy*dy);
                    if (d < bestSD) { bestSD = d; bestSx = nx; bestSy = ny; }
                }
            }
        }
    }

    bestEx = ex; bestEy = ey; bestED = 999999;
    for (dy = -snapRadius; dy <= snapRadius; dy++) {
        for (dx = -snapRadius; dx <= snapRadius; dx++) {
            nx = ex + dx; ny = ey + dy;
            if (nx >= 0 && nx < getWidth() && ny >= 0 && ny < getHeight()) {
                if (getPixel(nx, ny) > 0) {
                    d = sqrt(dx*dx + dy*dy);
                    if (d < bestED) { bestED = d; bestEx = nx; bestEy = ny; }
                }
            }
        }
    }

    if (bestSD > snapRadius || bestED > snapRadius) {
        print("[WARNUNG] Kontakt " + (ci+1) + ": Kein Skeleton in der Naehe. Uebersprungen.");
        pairLengths[ci] = -1;
        ci++;
        continue;
    }

    // Skeleton-Pfad tracen (BFS)
    traceResult = traceSkeleton(skelID, bestSx, bestSy, bestEx, bestEy);

    if (traceResult == "") {
        print("[WARNUNG] Kontakt " + (ci+1) + ": Kein zusammenhaengender Pfad. Uebersprungen.");
        pairLengths[ci] = -1;
        ci++;
        continue;
    }

    // Pfadlaenge berechnen (kalibriert)
    points = split(traceResult, ";");
    pathLength = 0;
    for (k = 0; k < points.length - 1; k++) {
        c1 = split(points[k], ",");
        c2 = split(points[k + 1], ",");
        pdx = (parseFloat(c2[0]) - parseFloat(c1[0])) * pixelWidth;
        pdy = (parseFloat(c2[1]) - parseFloat(c1[1])) * pixelHeight;
        pathLength += sqrt(pdx*pdx + pdy*pdy);
    }

    // Euklidische Distanz (kalibriert)
    edx = (bestEx - bestSx) * pixelWidth;
    edy = (bestEy - bestSy) * pixelHeight;
    euclidDist = sqrt(edx*edx + edy*edy);

    pairLengths[ci] = pathLength;
    pairEuclid[ci]  = euclidDist;
    pairStartX[ci]  = bestSx;
    pairStartY[ci]  = bestSy;
    pairEndX[ci]    = bestEx;
    pairEndY[ci]    = bestEy;

    if (pathLength > 0)
        pairZigzag[ci] = euclidDist / pathLength;
    else
        pairZigzag[ci] = 0;

    nMeasured++;

    // ---- Overlay zeichnen ----
    selectImage(zo1ID);

    // Membranpfad (gruen)
    for (k = 0; k < points.length - 1; k++) {
        c1 = split(points[k], ",");
        c2 = split(points[k + 1], ",");
        setColor(0, 255, 0);
        Overlay.drawLine(parseFloat(c1[0]), parseFloat(c1[1]),
                         parseFloat(c2[0]), parseFloat(c2[1]));
        Overlay.setPosition(0);
    }

    // Gerade Strecke (rot)
    setColor(255, 0, 0);
    Overlay.drawLine(bestSx, bestSy, bestEx, bestEy);
    Overlay.setPosition(0);

    // Endpunkte (magenta)
    r = 4;
    setColor(255, 0, 255);
    Overlay.drawEllipse(bestSx - r, bestSy - r, 2*r, 2*r);
    Overlay.setPosition(0);
    Overlay.drawEllipse(bestEx - r, bestEy - r, 2*r, 2*r);
    Overlay.setPosition(0);

    // Nummer (gelb)
    setColor(255, 255, 0);
    midX = (bestSx + bestEx) / 2;
    midY = (bestSy + bestEy) / 2;
    Overlay.drawString("" + (ci+1), midX + 5, midY - 5);
    Overlay.setPosition(0);

    ci++;
}

Overlay.show();

if (nMeasured == 0) {
    closeImagesStartingWith("ZZ_");
    exit("Kein Kontakt konnte gemessen werden.\n" +
         "Skeleton-Pfade konnten nicht getracet werden.\n" +
         "Preprocessing-Parameter ueberpruefen.");
}

// ---- Results-Tabelle ----
tableName = "Zigzag Index Results (Auto)";
Table.create(tableName);

sumZZ  = 0;
sumZZ2 = 0;
row = 0;

for (p = 0; p < nContacts; p++) {
    if (pairLengths[p] < 0) continue;
    Table.set("Bild", row, origTitle);
    Table.set("Kanal", row, channelInfo);
    Table.set("Nr", row, p + 1);
    Table.set("Gerade (" + unit + ")", row, pairEuclid[p]);
    Table.set("Membran (" + unit + ")", row, pairLengths[p]);
    Table.set("Zickzack-Index", row, pairZigzag[p]);
    sumZZ  += pairZigzag[p];
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
closeImagesStartingWith("ZZ_work");
closeImagesStartingWith("Tagged skeleton");
closeImagesStartingWith("Mask of");

selectImage(zo1ID);

print("================================================");
print("Zigzag-Index Analyse (AUTOMATISCH) abgeschlossen");
print("Bild: " + origTitle);
print("Kanal: " + channelInfo);
print("Erkannte Kontakte: " + nContacts);
print("Davon gemessen: " + nMeasured);
print("Mittelwert Zigzag-Index: " + d2s(meanZZ, 4));
print("Standardabweichung: " + d2s(sdZZ, 4));
print("================================================");


// =============================================================================
// Hilfsfunktionen
// =============================================================================

function traceSkeleton(imgID, startX, startY, endX, endY) {
    selectImage(imgID);
    w = getWidth();
    h = getHeight();

    maxSteps = w * h;
    if (maxSteps > 50000) maxSteps = 50000;

    visited = newArray(w * h);
    queueX  = newArray(maxSteps);
    queueY  = newArray(maxSteps);
    parent  = newArray(maxSteps);
    qHead = 0;
    qTail = 0;

    queueX[qTail] = startX;
    queueY[qTail] = startY;
    parent[qTail] = -1;
    visited[startY * w + startX] = 1;
    qTail++;

    found = false;
    foundIdx = -1;

    while (qHead < qTail && !found) {
        cx = queueX[qHead];
        cy = queueY[qHead];
        ci2 = qHead;
        qHead++;

        if (cx == endX && cy == endY) {
            found = true;
            foundIdx = ci2;
            break;
        }

        for (dy = -1; dy <= 1; dy++) {
            for (dx = -1; dx <= 1; dx++) {
                if (dx == 0 && dy == 0) continue;
                nx = cx + dx;
                ny = cy + dy;
                if (nx >= 0 && nx < w && ny >= 0 && ny < h) {
                    nIdx = ny * w + nx;
                    if (visited[nIdx] == 0 && getPixel(nx, ny) > 0) {
                        visited[nIdx] = 1;
                        if (qTail < maxSteps) {
                            queueX[qTail] = nx;
                            queueY[qTail] = ny;
                            parent[qTail] = ci2;
                            qTail++;
                        }
                    }
                }
            }
        }
    }

    if (!found) return "";

    pathStr = "";
    idx = foundIdx;
    revX = newArray(qTail);
    revY = newArray(qTail);
    revCount = 0;
    while (idx >= 0) {
        revX[revCount] = queueX[idx];
        revY[revCount] = queueY[idx];
        revCount++;
        idx = parent[idx];
    }

    for (k = revCount - 1; k >= 0; k--) {
        if (k < revCount - 1) pathStr += ";";
        pathStr += "" + revX[k] + "," + revY[k];
    }

    return pathStr;
}

function closeWindowIfOpen(name) {
    if (isOpen(name)) {
        selectWindow(name);
        run("Close");
    }
}

function closeImagesStartingWith(prefix) {
    list = getList("image.titles");
    for (i = 0; i < list.length; i++) {
        if (startsWith(list[i], prefix)) {
            selectWindow(list[i]);
            close();
        }
    }
}
