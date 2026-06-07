"""
topology_render.py  —  Drawing helpers for plate-boundary topology overlays

Used by render_videos_cartopy.py and render_videos_pygmt.py to draw:
  - non-subduction plate boundaries as thin grey lines
  - subduction zones as red lines with triangular teeth pointing toward
    the overriding plate side (Left / Right of the line as it's walked
    in vertex order)

The teeth follow the standard GPlately / pyGPlates plotting convention:
the triangle's *base* sits on the trench line and the *apex* points
into the overriding plate.

This module is dependency-light on purpose — the caller passes in
`ccrs` (cartopy) or the active `pygmt.Figure`, so this file doesn't
need to import either at module scope.  That keeps the import graph
clean for non-rendering callers of plate_model_utils.
"""
from __future__ import annotations
import numpy as np


# ----------------------------------------------------------------------
# cartopy
# ----------------------------------------------------------------------

def draw_topologies_cartopy(ax, topology, ccrs,
                            other_pen=None, sz_pen=None,
                            polygon_pen=None,
                            tooth_size_deg=2.8,
                            tooth_half_base_deg=1.1,
                            tooth_spacing_pts=4):
    """Draw plate-boundary topology on a cartopy axis.

    Parameters
    ----------
    ax : cartopy GeoAxes
    topology : dict
        Result of plate_model_utils.topology_lines(t).
    ccrs : cartopy.crs module
        Passed in so this module doesn't hard-import cartopy.
    other_pen, sz_pen, polygon_pen : dict
        Matplotlib line kwargs for the non-SZ, SZ, and closed-plate-polygon
        backstop lines respectively.
    tooth_size_deg : float
        Perpendicular distance from the trench to the tooth apex.
    tooth_half_base_deg : float
        Half the along-trench length of the tooth base.
    tooth_spacing_pts : int
        Place a tooth every Nth segment.  Topology polylines are
        typically densely sampled at ~0.5–1° spacing, so 3 gives a
        tooth roughly every 2°.
    """
    if other_pen is None:
        # All cartopy boundaries are red, matching the pyGMT convention.
        # SZs are distinguished from ridges/transforms by the triangular
        # teeth, not by line colour.
        other_pen = dict(color="red", linewidth=0.5)
    if sz_pen is None:
        sz_pen = dict(color="red", linewidth=0.7)
    if polygon_pen is None:
        # Thin grey backstop — drawn underneath the red overlay so each
        # plate is visibly closed even where individual sub-segments are
        # unlabelled or dropped by feature-type filtering.  Pattern
        # mirrors GPlately's ``plot_all_topological_sections`` +
        # ``plot_plate_polygon_by_id`` combination.
        polygon_pen = dict(color="0.45", linewidth=0.35)

    pc = ccrs.PlateCarree()

    # Closed plate-polygon rings — drawn FIRST as a thin grey backstop
    # underneath the colour-coded SZ/other overlays AND under the COB
    # outline (zorder=3).  See plate_model_utils.plate_polygons for why
    # this is needed.
    for pts in topology.get("polygons", []):
        for chunk in _split_polyline_at_dateline(pts):
            ax.plot(chunk[:, 1], chunk[:, 0], transform=pc, zorder=2,
                    **polygon_pen)

    # Non-subduction plate boundaries — single colour, thin lines.
    # Split any polyline that crosses the dateline so matplotlib doesn't
    # draw a stray horizontal line across the entire map (cartopy with
    # the PlateCarree transform connects (-179°, +179°) literally rather
    # than wrapping the short way).
    for pts in topology.get("other", []):
        for chunk in _split_polyline_at_dateline(pts):
            ax.plot(chunk[:, 1], chunk[:, 0], transform=pc, zorder=4,
                    **other_pen)

    # Subduction zones — line + teeth (line also dateline-split).
    for pts, polarity in topology.get("subduction", []):
        for chunk in _split_polyline_at_dateline(pts):
            ax.plot(chunk[:, 1], chunk[:, 0], transform=pc, zorder=4,
                    **sz_pen)
        _draw_sz_teeth_cartopy(
            ax, pts, polarity, pc,
            tooth_size_deg=tooth_size_deg,
            tooth_half_base_deg=tooth_half_base_deg,
            tooth_spacing_pts=tooth_spacing_pts,
            color=sz_pen.get("color", "red"),
        )


def draw_cob_cartopy(ax, cob_lines, ccrs, pen=None):
    """Draw COB (continent-ocean boundary) land polygons as black outlines
    on a cartopy axis.

    Parameters
    ----------
    ax : cartopy GeoAxes
    cob_lines : list of (N, 2) [lat, lon] arrays
        Result of ``plate_model_utils.cob_polylines(t)``.
    ccrs : cartopy.crs module
        Passed in so this module doesn't hard-import cartopy.
    pen : dict, optional
        Matplotlib line kwargs.  Default: thin solid black line.

    Polygons that cross the antimeridian are split by
    ``_split_polyline_at_dateline`` so matplotlib doesn't draw a
    spurious horizontal segment across the whole map.  The COB lines
    sit at ``zorder=3`` — below the red topology overlay (zorder=4)
    and its subduction-zone teeth (zorder=5) — so plate-boundary
    detail remains visible where the two coincide along continental
    margins.
    """
    if pen is None:
        pen = dict(color="black", linewidth=0.6)
    pc = ccrs.PlateCarree()
    for pts in cob_lines:
        for chunk in _split_polyline_at_dateline(pts):
            ax.plot(chunk[:, 1], chunk[:, 0], transform=pc, zorder=3,
                    **pen)


def _split_polyline_at_dateline(pts, max_lon_jump_deg=180.0):
    """Yield contiguous chunks of `pts` (an (N, 2) [lat, lon] array)
    such that no consecutive longitude jump within a chunk exceeds
    `max_lon_jump_deg`.

    Topology polylines that genuinely cross the antimeridian appear in
    the resolved geometry as two adjacent vertices at e.g. (+179°, -179°);
    matplotlib + the PlateCarree transform draws a literal segment
    between them — a horizontal stripe across the whole map.  Splitting
    the polyline at such jumps suppresses the spurious connecting line.
    """
    if pts is None or len(pts) < 2:
        return
    lon = pts[:, 1]
    chunk_start = 0
    for i in range(1, len(pts)):
        if abs(lon[i] - lon[i - 1]) > max_lon_jump_deg:
            if i - chunk_start >= 2:
                yield pts[chunk_start:i]
            chunk_start = i
    if len(pts) - chunk_start >= 2:
        yield pts[chunk_start:]


def _draw_sz_teeth_cartopy(ax, pts, polarity, pc,
                           tooth_size_deg=1.5,
                           tooth_half_base_deg=0.5,
                           tooth_spacing_pts=3,
                           color="red"):
    """Filled triangle teeth along an SZ polyline.  Tip points into the
    overriding plate (Left or Right of the line in walking order)."""
    if len(pts) < 2:
        return
    # Default to "Right" (most common convention) when polarity is missing.
    side = polarity if polarity in ("Left", "Right") else "Right"
    sign = +1.0 if side == "Left" else -1.0

    for i in range(0, len(pts) - 1, tooth_spacing_pts):
        lat0, lon0 = pts[i]
        lat1, lon1 = pts[i + 1]
        # Skip vertex pairs that straddle the dateline — drawing a tooth
        # at their midpoint would put it at lon ≈ 0 on the wrong side
        # of the globe.
        if abs(lon1 - lon0) > 180.0:
            continue
        mlat = 0.5 * (lat0 + lat1)
        mlon = 0.5 * (lon0 + lon1)
        coslat = max(0.05, np.cos(np.radians(mlat)))

        # Tangent in METRIC space (longitude scaled by cos(lat) so the
        # perpendicular comes out perpendicular on the actual sphere
        # rather than on the lat/lon graticule).
        tx_lat = lat1 - lat0
        ty_mlon = (lon1 - lon0) * coslat
        norm = float(np.hypot(tx_lat, ty_mlon))
        if norm <= 1e-9:
            continue
        tx_lat /= norm
        ty_mlon /= norm

        # Left-perpendicular in metric space, converted back to
        # (Δlat, Δlon-degrees) for cartopy's PlateCarree transform.
        px_lat = -ty_mlon
        py_lon_deg = tx_lat / coslat
        tip_lat = mlat + sign * tooth_size_deg * px_lat
        tip_lon = mlon + sign * tooth_size_deg * py_lon_deg

        # Base sits on the line — two points either side of midpoint
        # along the tangent direction.
        b1_lat = mlat - tooth_half_base_deg * tx_lat
        b1_lon = mlon - tooth_half_base_deg * ty_mlon / coslat
        b2_lat = mlat + tooth_half_base_deg * tx_lat
        b2_lon = mlon + tooth_half_base_deg * ty_mlon / coslat

        ax.fill([b1_lon, b2_lon, tip_lon],
                [b1_lat, b2_lat, tip_lat],
                facecolor=color, edgecolor=color, linewidth=0.3,
                transform=pc, zorder=5)


# ----------------------------------------------------------------------
# pyGMT
# ----------------------------------------------------------------------

def draw_topologies_pygmt(fig, topology, projection, region,
                          other_pen="0.4p,red",
                          sz_pen="0.7p,red",
                          sz_fill="red",
                          polygon_pen="0.25p,gray55",
                          tooth_spacing="0.525c",
                          tooth_size="4.5p"):
    """Draw plate-boundary topology on a pygmt Figure.

    Non-SZ boundaries are thin red lines.  Subduction zones use GMT's
    'front' line style (``-Sf``): a line with teeth at a fixed gap on
    one side (``+l`` left, ``+r`` right) and filled triangles (``+t``).

    A thin grey closed-plate-polygon backstop (``polygon_pen``) is
    drawn FIRST so that plates visually close even when individual
    sub-segments are unlabelled or filtered out by feature-type checks.
    See plate_model_utils.plate_polygons for the rationale.
    """
    # Closed plate-polygon rings — backstop drawn first so the
    # colour-coded SZ/other overlays render on top.
    for pts in topology.get("polygons", []):
        if len(pts) < 2:
            continue
        fig.plot(x=pts[:, 1], y=pts[:, 0], pen=polygon_pen,
                 projection=projection, region=region)

    # Non-subduction plate boundaries — single thin-line plot per segment.
    for pts in topology.get("other", []):
        if len(pts) < 2:
            continue
        fig.plot(x=pts[:, 1], y=pts[:, 0], pen=other_pen,
                 projection=projection, region=region)

    # Subduction zones — partition by polarity so we can plot each
    # group in one go with the correct ``+l`` / ``+r`` front style.
    sz_left, sz_right = [], []
    for pts, polarity in topology.get("subduction", []):
        if len(pts) < 2:
            continue
        if polarity == "Left":
            sz_left.append(pts)
        else:  # "Right" or "Unknown" — both default to right-teeth.
            sz_right.append(pts)

    style_l = f"f{tooth_spacing}/{tooth_size}+l+t"
    style_r = f"f{tooth_spacing}/{tooth_size}+r+t"

    for pts in sz_left:
        fig.plot(x=pts[:, 1], y=pts[:, 0], pen=sz_pen,
                 style=style_l, fill=sz_fill,
                 projection=projection, region=region)
    for pts in sz_right:
        fig.plot(x=pts[:, 1], y=pts[:, 0], pen=sz_pen,
                 style=style_r, fill=sz_fill,
                 projection=projection, region=region)
