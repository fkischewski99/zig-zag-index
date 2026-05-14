"""
zigzag_pathfind.py
Findet den Membranpfad zwischen Punktpaaren via Dijkstra.

Aufruf:
  python zigzag_pathfind.py <image> <coords> <results> <paths> <px_w> <px_h> [circle [<arc_paths>]]

  circle      - Optional: Kreisbogen-Referenz berechnen (Least-Squares Circle Fit)
  arc_paths   - Optional: CSV-Datei fuer Kreisbogen-Overlaypunkte

Abhaengigkeiten: pip install numpy scikit-image
"""

import sys
import csv
import numpy as np
from skimage import io
from skimage.graph import route_through_array
from scipy.ndimage import gaussian_filter


def fit_circle(points_cal):
    """Least-Squares Circle Fit (Kasa-Methode).
    points_cal: Nx2 Array kalibrierter Koordinaten.
    Gibt (a, b, r) zurueck: Mittelpunkt und Radius."""
    x = points_cal[:, 0]
    y = points_cal[:, 1]
    A = np.column_stack([x, y, np.ones(len(x))])
    b = x**2 + y**2
    result, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    a_c = result[0] / 2.0
    b_c = result[1] / 2.0
    r = np.sqrt(result[2] + a_c**2 + b_c**2)
    return a_c, b_c, r


def arc_length_from_circle(r, euclid):
    """Bogenlaenge aus Radius und Sehnenlaenge.
    L = r * 2 * arcsin(d / (2r)). Fallback auf euclid wenn d >= 2r."""
    ratio = euclid / (2.0 * r)
    if ratio >= 1.0:
        return euclid
    return r * 2.0 * np.arcsin(ratio)


def generate_arc_points(a, b, r, start_cal, end_cal, path_cal, n_points=100):
    """Erzeugt Punkte entlang des Kreisbogens zwischen start und end.
    Waehlt die Bogenseite, auf der die Pfadpunkte liegen."""
    angle_start = np.arctan2(start_cal[1] - b, start_cal[0] - a)
    angle_end = np.arctan2(end_cal[1] - b, end_cal[0] - a)

    path_mid = path_cal[len(path_cal) // 2]
    angle_mid = np.arctan2(path_mid[1] - b, path_mid[0] - a)

    diff1 = (angle_end - angle_start) % (2 * np.pi)
    diff2 = (angle_start - angle_end) % (2 * np.pi)

    mid_in_forward = (angle_mid - angle_start) % (2 * np.pi) <= diff1

    if mid_in_forward:
        angles = np.linspace(angle_start, angle_start + diff1, n_points)
    else:
        angles = np.linspace(angle_start, angle_start - diff2, n_points)

    xs = a + r * np.cos(angles)
    ys = b + r * np.sin(angles)
    return np.column_stack([xs, ys])


def distance_to_line_segment(shape, start, end):
    """Abstand jedes Pixels zur Strecke start-end. (row, col)."""
    rows, cols = np.mgrid[0:shape[0], 0:shape[1]]
    r1, c1 = float(start[0]), float(start[1])
    r2, c2 = float(end[0]), float(end[1])
    dr, dc = r2 - r1, c2 - c1
    seg_len_sq = dr * dr + dc * dc
    if seg_len_sq < 1e-10:
        return np.sqrt((rows - r1)**2 + (cols - c1)**2)
    t = np.clip(((rows - r1) * dr + (cols - c1) * dc) / seg_len_sq, 0.0, 1.0)
    closest_r = r1 + t * dr
    closest_c = c1 + t * dc
    return np.sqrt((rows - closest_r)**2 + (cols - closest_c)**2)


def main():
    image_path = sys.argv[1]
    coords_path = sys.argv[2]
    results_path = sys.argv[3]
    paths_path = sys.argv[4]
    pixel_width = float(sys.argv[5])
    pixel_height = float(sys.argv[6])

    compute_arc = len(sys.argv) > 7 and sys.argv[7] == "circle"
    arc_paths_path = sys.argv[8] if len(sys.argv) > 8 else None

    # Bild laden
    img = io.imread(image_path)
    if img.ndim == 3:
        img = np.mean(img, axis=2)
    img = img.astype(np.float64)
    print(f"Bild: {img.shape[1]}x{img.shape[0]}", flush=True)

    # Leichte Glaettung gegen Pixelrauschen
    img_smooth = gaussian_filter(img, sigma=1.0)

    # Normalisieren auf [0, 1]
    img_min = img_smooth.min()
    img_max = img_smooth.max()
    if img_max - img_min > 0:
        img_norm = (img_smooth - img_min) / (img_max - img_min)
    else:
        img_norm = np.zeros_like(img_smooth)

    # Basis-Kostenkarte: exponentiell
    # Hell (Membran, norm=1): cost = e^0 = 1
    # Dunkel (Hintergrund, norm=0): cost = e^15 = 3.3 Mio
    # -> Pfad verlässt NIEMALS die hellen Pixel
    base_cost = np.exp(15.0 * (1.0 - img_norm))
    print("Kostenkarte berechnet", flush=True)

    # Koordinaten lesen
    coords = []
    with open(coords_path, 'r') as f:
        for row in csv.reader(f):
            if len(row) >= 2:
                coords.append((float(row[0]), float(row[1])))

    n_pairs = len(coords) // 2
    results = []
    all_paths = []
    all_arc_paths = []

    for p in range(n_pairs):
        x1, y1 = coords[p * 2]
        x2, y2 = coords[p * 2 + 1]

        start = (
            max(0, min(int(round(y1)), img.shape[0] - 1)),
            max(0, min(int(round(x1)), img.shape[1] - 1))
        )
        end = (
            max(0, min(int(round(y2)), img.shape[0] - 1)),
            max(0, min(int(round(x2)), img.shape[1] - 1))
        )

        print(f"Paar {p+1}/{n_pairs}: ({start[1]},{start[0]}) -> ({end[1]},{end[0]})", flush=True)

        # Proximity-Kosten: Pixel weit von der Verbindungslinie = teurer
        # Verhindert dass der Pfad ueber eine ANDERE Membran laeuft
        line_dist = distance_to_line_segment(img.shape, start, end)
        proximity = 1.0 + (line_dist ** 2) * 0.005

        cost = base_cost * proximity

        try:
            path_indices, _ = route_through_array(
                cost, start, end,
                fully_connected=True,
                geometric=True
            )
            path = np.array(path_indices)

            # Pfadlaenge kalibriert
            diffs = np.diff(path, axis=0).astype(np.float64)
            diffs[:, 0] *= pixel_height
            diffs[:, 1] *= pixel_width
            path_length = float(np.sum(np.sqrt(np.sum(diffs**2, axis=1))))

            # Euclidean Distance kalibriert
            dx = (end[1] - start[1]) * pixel_width
            dy = (end[0] - start[0]) * pixel_height
            euclid = float(np.sqrt(dx**2 + dy**2))

            zigzag = path_length / euclid if euclid > 0 else 0.0

            arc_len = 0.0
            arc_zigzag = 0.0
            if compute_arc and len(path) >= 4:
                path_cal = path.astype(np.float64).copy()
                path_cal[:, 0] *= pixel_height
                path_cal[:, 1] *= pixel_width
                try:
                    a_c, b_c, r_c = fit_circle(path_cal)
                    arc_len = float(arc_length_from_circle(r_c, euclid))
                    arc_zigzag = path_length / arc_len if arc_len > 0 else 0.0

                    start_cal = np.array([start[0] * pixel_height, start[1] * pixel_width])
                    end_cal = np.array([end[0] * pixel_height, end[1] * pixel_width])
                    arc_pts_cal = generate_arc_points(a_c, b_c, r_c, start_cal, end_cal, path_cal)
                    arc_pts_px = arc_pts_cal.copy()
                    arc_pts_px[:, 0] /= pixel_height
                    arc_pts_px[:, 1] /= pixel_width
                    for k in range(len(arc_pts_px)):
                        all_arc_paths.append((p + 1, int(round(arc_pts_px[k, 1])), int(round(arc_pts_px[k, 0]))))
                except Exception as e:
                    print(f"  [WARNUNG] Kreisbogen-Fit fehlgeschlagen: {e}", flush=True)
                    arc_len = euclid
                    arc_zigzag = zigzag
            elif compute_arc:
                arc_len = euclid
                arc_zigzag = zigzag

            if compute_arc:
                results.append((p + 1, euclid, path_length, zigzag, arc_len, arc_zigzag))
                print(f"  Laenge={path_length:.2f}, Euklid={euclid:.2f}, ZZ={zigzag:.4f}"
                      f", Bogen={arc_len:.2f}, BogenZZ={arc_zigzag:.4f}", flush=True)
            else:
                results.append((p + 1, euclid, path_length, zigzag))
                print(f"  Laenge={path_length:.2f}, Euklid={euclid:.2f}, ZZ={zigzag:.4f}", flush=True)

            # Pfadpunkte speichern
            for k in range(len(path)):
                all_paths.append((p + 1, int(path[k, 1]), int(path[k, 0])))

        except Exception as e:
            print(f"[FEHLER] Paar {p+1}: {e}", file=sys.stderr, flush=True)
            if compute_arc:
                results.append((p + 1, 0, 0, 0, 0, 0))
            else:
                results.append((p + 1, 0, 0, 0))

    with open(results_path, 'w', newline='') as f:
        for row in results:
            csv.writer(f).writerow(row)

    with open(paths_path, 'w', newline='') as f:
        for row in all_paths:
            csv.writer(f).writerow(row)

    if compute_arc and arc_paths_path:
        with open(arc_paths_path, 'w', newline='') as f:
            for row in all_arc_paths:
                csv.writer(f).writerow(row)

    print(f"OK: {n_pairs} Paare", flush=True)


if __name__ == '__main__':
    main()
