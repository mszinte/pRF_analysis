import os
import cortex
import tempfile
import numpy as np
from lxml import etree


from skimage.measure import find_contours
from skimage.morphology import binary_closing, disk

from cortex import db
from cortex.svgoverlay import get_overlay, _make_layer, _find_layer, parser
from cortex.dataset import Vertex


class ROIpack(object):
    def __init__(self, subject, roifile):
        self.subject = subject
        self.roifile = roifile
        self.rois = {}
        self.load_roifile()

    def load_roifile(self):
        if not os.path.exists(self.roifile):
            raise IOError(f"ROI file {self.roifile} not found")

        empty = Vertex(None, self.subject)

        if self.roifile.endswith("npz"):
            roidata = np.load(self.roifile)
            for roi in roidata.keys():
                self.rois[roi] = empty.copy(roidata[roi])
            roidata.close()
        else:
            raise ValueError("Only NPZ supported for this implementation")

    def to_svg(self, res=2048, filename=None):
        if filename is None:
            filename = tempfile.mktemp(
                suffix=".svg", prefix=self.subject + "-rois-"
            )

        # base flatmap image (same as quickflat)
        mpts, mpolys = db.get_surf(
            self.subject, "flat", merge=True, nudge=True
        )

        svgmpts = mpts[:, :2].copy()
        svgmpts -= svgmpts.min(0)
        svgmpts *= 1024 / svgmpts.max(0)[1]
        svgmpts[:, 1] = 1024 - svgmpts[:, 1]

        # number of vertices per hemisphere
        lh_pts, _ = db.get_surf(self.subject, "flat")[0]
        n_lh = len(lh_pts)

        coords_lh = svgmpts[:n_lh]
        coords_rh = svgmpts[n_lh:]

        svgroipack = get_overlay(
            self.subject, filename, mpts, mpolys
        )
        svg = etree.parse(svgroipack.svgfile, parser=parser)
        rois_layer = _find_layer(
            _find_layer(svg, "rois"), "shapes"
        )

        
        roi_names = list(self.rois.keys())
        masks_lh = [self.rois[r].left for r in roi_names]
        masks_rh = [self.rois[r].right for r in roi_names]

        label_lh, xgrid_lh, ygrid_lh = label_masks(masks_lh, coords_lh, res=res)
        label_rh, xgrid_rh, ygrid_rh = label_masks(masks_rh, coords_rh, res=res)

        
        for idx, roi_name in enumerate(roi_names):
            print(f"Adding {roi_name}")
            roi_layer = _make_layer(rois_layer, roi_name)

            for label_img, xgrid, ygrid in (
                (label_lh, xgrid_lh, ygrid_lh),
                (label_rh, xgrid_rh, ygrid_rh),
            ):
                img = (label_img == (idx + 1)).astype(np.uint8)

                if img.sum() == 0:
                    continue

                img = binary_closing(img, disk(2))

                contours = find_contours(img, 0.5)
                if not contours:
                    continue

                contour = max(contours, key=len)

                pts = np.column_stack([
                    xgrid[contour[:, 1].astype(int)],
                    ygrid[contour[:, 0].astype(int)],
                ])

                svgpath = etree.SubElement(roi_layer, "path")
                svgpath.attrib["style"] = (
                    "fill:none;"
                    "stroke:#000000;"
                    "stroke-width:1px;"
                    "stroke-linecap:round;"
                    "stroke-linejoin:round"
                )
                svgpath.attrib["d"] = svg_path_from_points(
                    pts,
                    close=np.linalg.norm(pts[0] - pts[-1]) < 2
                )

        with open(svgroipack.svgfile, "wb") as f:
            f.write(etree.tostring(svg, pretty_print=True))


def label_masks(mask_list, coords, res=2048, radius=2):
    """
    Convertit une liste de masques ROI (1D vertex) en image labelisée,
    avec une grille commune par hémisphère.
    """
    xmin, ymin = coords.min(0)
    xmax, ymax = coords.max(0)

    xgrid = np.linspace(xmin, xmax, res)
    ygrid = np.linspace(ymin, ymax, res)

    # mapping vertex -> pixel
    x = np.clip(((coords[:, 0] - xmin) / (xmax - xmin) * (res - 1)).astype(int), 0, res - 1)
    y = np.clip(((coords[:, 1] - ymin) / (ymax - ymin) * (res - 1)).astype(int), 0, res - 1)

    label_img = np.zeros((res, res), dtype=np.int32)

    for idx, mask in enumerate(mask_list):
        roi_idx = np.where(mask > 0)[0]
        for i in roi_idx:
            cx, cy = x[i], y[i]
            label_img[
                max(0, cy - radius):min(res, cy + radius + 1),
                max(0, cx - radius):min(res, cx + radius + 1),
            ] = idx + 1

    return label_img, xgrid, ygrid


def svg_path_from_points(pts, close=False):
    pts = np.nan_to_num(pts)
    d = f"M {pts[0,0]:.2f} {pts[0,1]:.2f} "
    d += " ".join(f"L {x:.2f} {y:.2f}" for x, y in pts[1:])
    if close:
        d += " Z"
    return d