#!/usr/bin/env python3
"""
Standalone GDML geometry viewer.

This tool parses a GDML file, converts supported solids into triangle meshes,
and writes a self-contained WebGL HTML viewer. It is intentionally independent
of ROOT, Streamlit, Plotly, VTK, and Geant4 so it can run in a plain Python
environment.

Run:
    python3 gdml_viewer.py /path/to/geometry.gdml --serve

Notes on accuracy:
    The parser preserves the GDML volume hierarchy, placements, rotations,
    units, materials, and names. It supports box, tube, cutTube, cone, and union.
    General boolean subtraction/intersection is reported and approximated by the
    first operand in the pure-Python backend.

    For more accurate CSG tessellation, install pyg4ometry and run with:

        python3 gdml_viewer.py geometry.gdml --backend pyg4ometry --serve
"""

from __future__ import annotations

import argparse
import html
import json
import math
import os
import re
import struct
import sys
import tempfile
import webbrowser
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

Vec3 = Tuple[float, float, float]
Mat4 = List[List[float]]
Triangle = Tuple[Vec3, Vec3, Vec3]

LEN_UNITS = {
    None: 1.0, "": 1.0, "mm": 1.0, "millimeter": 1.0,
    "cm": 10.0, "m": 1000.0, "um": 0.001, "nm": 1e-6,
}
ANGLE_UNITS = {
    None: math.pi / 180.0, "": math.pi / 180.0,
    "deg": math.pi / 180.0, "degree": math.pi / 180.0,
    "rad": 1.0, "radian": 1.0,
}
QUALITY_SEGMENTS = {"fast": 12, "medium": 24, "high": 48}

MATERIAL_COLORS = {
    "Material_ECAL": (55, 130, 245, 185),
    "Material_HCAL": (235, 150, 35, 190),
    "G4_Si": (80, 195, 120, 220),
    "G4_Fe": (160, 165, 175, 225),
    "Galactic": (220, 230, 255, 18),
}


def tag_name(el: ET.Element | str) -> str:
    tag = el.tag if isinstance(el, ET.Element) else el
    return tag.rsplit("}", 1)[-1]


def child(el: ET.Element, name: str) -> Optional[ET.Element]:
    for item in list(el):
        if tag_name(item) == name:
            return item
    return None


def children(el: ET.Element, name: str) -> List[ET.Element]:
    return [item for item in list(el) if tag_name(item) == name]


def f(attrs: Dict[str, str], key: str, default: float = 0.0) -> float:
    try:
        return float(attrs.get(key, default))
    except (TypeError, ValueError):
        return default


def length(attrs: Dict[str, str], key: str, default: float = 0.0) -> float:
    return f(attrs, key, default) * LEN_UNITS.get(attrs.get("lunit") or attrs.get("unit"), 1.0)


def angle(attrs: Dict[str, str], key: str, default_deg: float = 0.0) -> float:
    return f(attrs, key, default_deg) * ANGLE_UNITS.get(attrs.get("aunit") or attrs.get("unit"), math.pi / 180.0)


def identity() -> Mat4:
    return [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]


def matmul(a: Mat4, b: Mat4) -> Mat4:
    out = [[0.0] * 4 for _ in range(4)]
    for i in range(4):
        for k in range(4):
            aik = a[i][k]
            if aik:
                for j in range(4):
                    out[i][j] += aik * b[k][j]
    return out


def translate(x: float, y: float, z: float) -> Mat4:
    m = identity()
    m[0][3], m[1][3], m[2][3] = x, y, z
    return m


def rot_x(a: float) -> Mat4:
    c, s = math.cos(a), math.sin(a)
    return [[1, 0, 0, 0], [0, c, -s, 0], [0, s, c, 0], [0, 0, 0, 1]]


def rot_y(a: float) -> Mat4:
    c, s = math.cos(a), math.sin(a)
    return [[c, 0, s, 0], [0, 1, 0, 0], [-s, 0, c, 0], [0, 0, 0, 1]]


def rot_z(a: float) -> Mat4:
    c, s = math.cos(a), math.sin(a)
    return [[c, -s, 0, 0], [s, c, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]


def euler_xyz(rx: float, ry: float, rz: float) -> Mat4:
    return matmul(rot_z(rz), matmul(rot_y(ry), rot_x(rx)))


def transform_point(m: Mat4, p: Vec3) -> Vec3:
    x, y, z = p
    return (
        m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3],
        m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3],
        m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3],
    )


def parse_vec(el: Optional[ET.Element], unit_map: Dict[Optional[str], float]) -> Vec3:
    if el is None:
        return (0.0, 0.0, 0.0)
    unit = unit_map.get(el.attrib.get("unit"), 1.0)
    return (f(el.attrib, "x") * unit, f(el.attrib, "y") * unit, f(el.attrib, "z") * unit)


def placement_matrix(physvol: ET.Element) -> Mat4:
    pos = parse_vec(child(physvol, "position"), LEN_UNITS)
    rot = parse_vec(child(physvol, "rotation"), ANGLE_UNITS)
    return matmul(translate(*pos), euler_xyz(*rot))


def add_tri(tris: List[Triangle], a: Vec3, b: Vec3, c: Vec3) -> None:
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    if nx * nx + ny * ny + nz * nz > 1e-16:
        tris.append((a, b, c))


def add_quad(tris: List[Triangle], a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> None:
    add_tri(tris, a, b, c)
    add_tri(tris, a, c, d)


def mesh_tube_like(
    rmin0: float, rmax0: float, rmin1: float, rmax1: float, zlen: float,
    start: float, delta: float, segments: int,
    low_plane: Optional[Vec3] = None, high_plane: Optional[Vec3] = None,
) -> List[Triangle]:
    if zlen <= 0 or max(rmax0, rmax1) <= 0:
        return []
    frac = min(1.0, abs(delta) / (2 * math.pi))
    n = max(3, int(math.ceil(segments * frac)))
    full = abs(delta) >= 2 * math.pi * 0.999
    phis = [start + delta * i / n for i in range(n + 1)]
    z0, z1 = -zlen / 2.0, zlen / 2.0

    def z_at(base: float, plane: Optional[Vec3], x: float, y: float) -> float:
        if plane is None or abs(plane[2]) < 1e-12:
            return base
        return base - (plane[0] * x + plane[1] * y) / plane[2]

    lo_o: List[Vec3] = []
    hi_o: List[Vec3] = []
    lo_i: List[Vec3] = []
    hi_i: List[Vec3] = []
    for phi in phis:
        co, si = math.cos(phi), math.sin(phi)
        x, y = rmax0 * co, rmax0 * si
        lo_o.append((x, y, z_at(z0, low_plane, x, y)))
        x, y = rmax1 * co, rmax1 * si
        hi_o.append((x, y, z_at(z1, high_plane, x, y)))
        x, y = rmin0 * co, rmin0 * si
        lo_i.append((x, y, z_at(z0, low_plane, x, y)))
        x, y = rmin1 * co, rmin1 * si
        hi_i.append((x, y, z_at(z1, high_plane, x, y)))

    tris: List[Triangle] = []
    has_hole = max(rmin0, rmin1) > 1e-9
    for i in range(n):
        j = i + 1
        add_quad(tris, lo_o[i], lo_o[j], hi_o[j], hi_o[i])
        if has_hole:
            add_quad(tris, lo_i[j], lo_i[i], hi_i[i], hi_i[j])
            add_quad(tris, lo_i[i], lo_i[j], lo_o[j], lo_o[i])
            add_quad(tris, hi_i[j], hi_i[i], hi_o[i], hi_o[j])
        else:
            add_tri(tris, (0, 0, z0), lo_o[j], lo_o[i])
            add_tri(tris, (0, 0, z1), hi_o[i], hi_o[j])
    if not full and has_hole:
        add_quad(tris, lo_i[0], lo_o[0], hi_o[0], hi_i[0])
        add_quad(tris, lo_o[-1], lo_i[-1], hi_i[-1], hi_o[-1])
    return tris


def mesh_box(attrs: Dict[str, str]) -> List[Triangle]:
    x, y, z = length(attrs, "x") / 2, length(attrs, "y") / 2, length(attrs, "z") / 2
    v = [(-x, -y, -z), (x, -y, -z), (x, y, -z), (-x, y, -z),
         (-x, -y, z), (x, -y, z), (x, y, z), (-x, y, z)]
    faces = [(0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1),
             (1, 5, 6, 2), (2, 6, 7, 3), (3, 7, 4, 0)]
    tris: List[Triangle] = []
    for a, b, c, d in faces:
        add_quad(tris, v[a], v[b], v[c], v[d])
    return tris


def mesh_solid_primitive(el: ET.Element, segments: int) -> List[Triangle]:
    attrs = el.attrib
    kind = tag_name(el)
    if kind == "box":
        return mesh_box(attrs)
    if kind == "tube":
        return mesh_tube_like(
            length(attrs, "rmin"), length(attrs, "rmax"),
            length(attrs, "rmin"), length(attrs, "rmax"),
            length(attrs, "z"), angle(attrs, "startphi"), angle(attrs, "deltaphi", 360), segments,
        )
    if kind == "cutTube":
        low = (f(attrs, "lowX"), f(attrs, "lowY"), f(attrs, "lowZ", -1.0))
        high = (f(attrs, "highX"), f(attrs, "highY"), f(attrs, "highZ", 1.0))
        return mesh_tube_like(
            length(attrs, "rmin"), length(attrs, "rmax"),
            length(attrs, "rmin"), length(attrs, "rmax"),
            length(attrs, "z"), angle(attrs, "startphi"), angle(attrs, "deltaphi", 360),
            segments, low, high,
        )
    if kind == "cone":
        return mesh_tube_like(
            length(attrs, "rmin1"), length(attrs, "rmax1"),
            length(attrs, "rmin2"), length(attrs, "rmax2"),
            length(attrs, "z"), angle(attrs, "startphi"), angle(attrs, "deltaphi", 360), segments,
        )
    return []


def color_for_material(name: str) -> Tuple[int, int, int, int]:
    for key, rgba in MATERIAL_COLORS.items():
        if name.startswith(key):
            return rgba
    h = 2166136261
    for ch in name:
        h = ((h ^ ord(ch)) * 16777619) & 0xFFFFFFFF
    return (80 + (h & 127), 80 + ((h >> 8) & 127), 80 + ((h >> 16) & 127), 195)


def semantic_name(name: str) -> str:
    return re.sub(r"0x[0-9a-fA-F]+$", "", name)


@dataclass
class VolumeDef:
    name: str
    material: str
    solid: str
    element: ET.Element
    seq: int


@dataclass
class Scene:
    positions: List[float] = field(default_factory=list)
    colors: List[int] = field(default_factory=list)
    ids: List[int] = field(default_factory=list)
    objects: List[Dict[str, object]] = field(default_factory=list)
    warnings: Dict[str, int] = field(default_factory=dict)
    bbox_min: List[float] = field(default_factory=lambda: [float("inf")] * 3)
    bbox_max: List[float] = field(default_factory=lambda: [float("-inf")] * 3)

    def warn(self, message: str) -> None:
        self.warnings[message] = self.warnings.get(message, 0) + 1


class GDMLSceneBuilder:
    def __init__(
        self,
        gdml: Path,
        quality: str,
        include_containers: bool,
        include_world: bool,
        max_objects: int,
        max_vertices: int,
    ) -> None:
        self.gdml = gdml
        self.quality = quality
        self.segments = QUALITY_SEGMENTS[quality]
        self.include_containers = include_containers
        self.include_world = include_world
        self.max_objects = max_objects
        self.max_vertices = max_vertices
        self.root = ET.parse(gdml).getroot()
        self.solids: Dict[str, ET.Element] = {}
        self.volumes: Dict[str, VolumeDef] = {}
        self.world = ""
        self.mesh_cache: Dict[str, Tuple[List[Triangle], str]] = {}
        self.scene = Scene()
        self.backend_name = "python"

    def parse(self) -> None:
        solids_el = child(self.root, "solids")
        structure = child(self.root, "structure")
        if solids_el is None or structure is None:
            raise ValueError("GDML must contain <solids> and <structure> sections.")
        for el in list(solids_el):
            if "name" in el.attrib:
                self.solids[el.attrib["name"]] = el
        counts: Dict[str, int] = {}
        for el in children(structure, "volume"):
            name = el.attrib.get("name", f"volume_{len(counts)}")
            seq = counts.get(name, 0)
            counts[name] = seq + 1
            mat = child(el, "materialref")
            solid = child(el, "solidref")
            self.volumes[name] = VolumeDef(
                name=name,
                material=mat.attrib.get("ref", "") if mat is not None else "",
                solid=solid.attrib.get("ref", "") if solid is not None else "",
                element=el,
                seq=seq,
            )
        setup = child(self.root, "setup")
        world = child(setup, "world") if setup is not None else None
        self.world = world.attrib.get("ref", "") if world is not None else next(reversed(self.volumes))

    def solid_mesh(self, ref: str) -> Tuple[List[Triangle], str]:
        if ref in self.mesh_cache:
            return self.mesh_cache[ref]
        el = self.solids.get(ref)
        if el is None:
            self.scene.warn(f"Missing solid reference: {ref}")
            self.mesh_cache[ref] = ([], "missing")
            return self.mesh_cache[ref]
        kind = tag_name(el)
        if kind in {"box", "tube", "cutTube", "cone"}:
            tris = mesh_solid_primitive(el, self.segments)
        elif kind in {"subtraction", "intersection"}:
            first = child(el, "first")
            first_ref = first.attrib.get("ref", "") if first is not None else ""
            tris, subkind = self.solid_mesh(first_ref)
            self.scene.warn(f"{kind} solid {ref} approximated by first operand {first_ref}")
            kind = f"{kind}/{subkind}"
        elif kind == "union":
            first = child(el, "first")
            second = child(el, "second")
            tris, subkind = self.solid_mesh(first.attrib.get("ref", "") if first is not None else "")
            tris2, _ = self.solid_mesh(second.attrib.get("ref", "") if second is not None else "")
            m = matmul(translate(*parse_vec(child(el, "position"), LEN_UNITS)), euler_xyz(*parse_vec(child(el, "rotation"), ANGLE_UNITS)))
            tris = tris + [(transform_point(m, a), transform_point(m, b), transform_point(m, c)) for a, b, c in tris2]
            kind = f"union/{subkind}"
        else:
            tris = []
            self.scene.warn(f"Unsupported solid type {kind}: {ref}")
        self.mesh_cache[ref] = (tris, kind)
        return self.mesh_cache[ref]

    def add_volume(self, volume: VolumeDef, matrix: Mat4, path: str) -> None:
        tris, solid_type = self.solid_mesh(volume.solid)
        if not tris:
            return
        obj_id = len(self.scene.objects)
        rgba = color_for_material(volume.material)
        first = len(self.scene.positions) // 3
        lo = [float("inf")] * 3
        hi = [float("-inf")] * 3
        for tri in tris:
            for point in tri:
                tp = transform_point(matrix, point)
                self.scene.positions.extend(tp)
                self.scene.colors.extend(rgba)
                self.scene.ids.append(obj_id)
                for i in range(3):
                    lo[i] = min(lo[i], tp[i])
                    hi[i] = max(hi[i], tp[i])
                    self.scene.bbox_min[i] = min(self.scene.bbox_min[i], tp[i])
                    self.scene.bbox_max[i] = max(self.scene.bbox_max[i], tp[i])
        count = len(self.scene.positions) // 3 - first
        self.scene.objects.append({
            "id": obj_id,
            "name": path.rsplit("/", 1)[-1],
            "semanticName": semantic_name(path.rsplit("/", 1)[-1]),
            "volume": volume.name,
            "semanticVolume": semantic_name(volume.name),
            "material": volume.material or "(none)",
            "solid": volume.solid,
            "solidType": solid_type,
            "firstVertex": first,
            "vertexCount": count,
            "bbox": [*lo, *hi],
            "path": path,
        })

    def traverse(self, volume: VolumeDef, matrix: Mat4, path: str) -> None:
        if len(self.scene.objects) >= self.max_objects or len(self.scene.positions) // 3 >= self.max_vertices:
            return
        phys = children(volume.element, "physvol")
        container = bool(phys) or volume.material.startswith("Galactic")
        render = (self.include_containers or not container) and (self.include_world or volume.name != self.world)
        if render:
            self.add_volume(volume, matrix, path)
        for pv in phys:
            ref_el = child(pv, "volumeref")
            if ref_el is None:
                continue
            ref = ref_el.attrib.get("ref", "")
            child_vol = self.volumes.get(ref)
            if child_vol is None:
                self.scene.warn(f"Missing volume reference: {ref}")
                continue
            pv_name = pv.attrib.get("name", ref)
            self.traverse(child_vol, matmul(matrix, placement_matrix(pv)), f"{path}/{pv_name}")
            if len(self.scene.objects) >= self.max_objects or len(self.scene.positions) // 3 >= self.max_vertices:
                self.scene.warn("Stopped at configured object/vertex limit")
                return

    def build(self) -> Dict[str, object]:
        self.parse()
        world = self.volumes.get(self.world)
        if world is None:
            raise ValueError(f"World volume {self.world!r} was not found.")
        self.traverse(world, identity(), self.world)
        if self.scene.bbox_min[0] == float("inf"):
            self.scene.bbox_min = [-1, -1, -1]
            self.scene.bbox_max = [1, 1, 1]
        return self.to_dict()

    def to_dict(self) -> Dict[str, object]:
        materials: Dict[str, int] = {}
        material_colors: Dict[str, List[int]] = {}
        solid_types: Dict[str, int] = {}
        bases: Dict[str, int] = {}
        for obj in self.scene.objects:
            materials[str(obj["material"])] = materials.get(str(obj["material"]), 0) + 1
            material_colors[str(obj["material"])] = list(color_for_material(str(obj["material"])))
            solid_types[str(obj["solidType"])] = solid_types.get(str(obj["solidType"]), 0) + 1
            bases[str(obj["semanticVolume"])] = bases.get(str(obj["semanticVolume"]), 0) + 1
        return {
            "meta": {
                "source": str(self.gdml),
                "sourceFile": self.gdml.name,
                "backend": self.backend_name,
                "quality": self.quality,
                "world": self.world,
                "objects": len(self.scene.objects),
                "vertices": len(self.scene.positions) // 3,
                "triangles": len(self.scene.positions) // 9,
                "bbox": [*self.scene.bbox_min, *self.scene.bbox_max],
                "materials": materials,
                "materialColors": material_colors,
                "solidTypes": solid_types,
                "semanticVolumeCounts": dict(sorted(bases.items())),
                "warnings": [f"{msg} (x{n})" if n > 1 else msg for msg, n in sorted(self.scene.warnings.items())],
            },
            "positions": self.scene.positions,
            "colors": self.scene.colors,
            "ids": self.scene.ids,
            "objects": self.scene.objects,
        }


def maybe_call(obj: object, names: Sequence[str]) -> object | None:
    for name in names:
        if hasattr(obj, name):
            value = getattr(obj, name)
            if callable(value):
                try:
                    return value()
                except TypeError:
                    continue
            return value
    return None


def numeric_triplet(value: object) -> Optional[Vec3]:
    seen: set[int] = set()
    return _numeric_triplet(value, seen)


def _numeric_triplet(value: object, seen: set[int]) -> Optional[Vec3]:
    if value is None:
        return None
    obj_id = id(value)
    if obj_id in seen:
        return None
    seen.add(obj_id)

    if hasattr(value, "eval"):
        try:
            value = value.eval()  # type: ignore[assignment]
        except TypeError:
            pass

    if hasattr(value, "pos"):
        point = _numeric_triplet(getattr(value, "pos"), seen)
        if point is not None:
            return point

    if all(hasattr(value, axis) for axis in ("x", "y", "z")):
        try:
            coords = []
            for axis in ("x", "y", "z"):
                attr = getattr(value, axis)
                coords.append(float(attr() if callable(attr) else attr))
            return (coords[0], coords[1], coords[2])
        except (TypeError, ValueError):
            pass

    if hasattr(value, "tolist"):
        try:
            value = value.tolist()  # type: ignore[assignment]
        except TypeError:
            pass

    if isinstance(value, (list, tuple)) and len(value) >= 3:
        try:
            return (float(value[0]), float(value[1]), float(value[2]))
        except (TypeError, ValueError):
            return None

    if not isinstance(value, (str, bytes)):
        try:
            seq = list(value)  # type: ignore[arg-type]
        except TypeError:
            seq = []
        if len(seq) >= 3:
            try:
                return (float(seq[0]), float(seq[1]), float(seq[2]))
            except (TypeError, ValueError):
                pass
    return None


def triangulate_vertex_face_mesh(vertices: object, faces: object) -> List[Triangle]:
    verts_raw = list(vertices)  # type: ignore[arg-type]
    verts: List[Vec3] = []
    for item in verts_raw:
        point = numeric_triplet(item)
        if point is not None:
            verts.append(point)
    if not verts:
        return []

    tris: List[Triangle] = []
    if hasattr(faces, "tolist"):
        faces = faces.tolist()  # type: ignore[assignment]
    if isinstance(faces, tuple):
        faces = list(faces)

    face_list: List[object]
    if isinstance(faces, list) and faces and all(isinstance(x, int) for x in faces):
        face_list = []
        i = 0
        while i < len(faces):
            n = int(faces[i])
            face_list.append(faces[i + 1 : i + 1 + n])
            i += n + 1
    else:
        face_list = list(faces)  # type: ignore[arg-type]

    for face in face_list:
        if hasattr(face, "tolist"):
            face = face.tolist()
        if not isinstance(face, (list, tuple)) or len(face) < 3:
            continue
        try:
            indices = [int(x) for x in face]
        except (TypeError, ValueError):
            continue
        for i in range(1, len(indices) - 1):
            a, b, c = indices[0], indices[i], indices[i + 1]
            if 0 <= a < len(verts) and 0 <= b < len(verts) and 0 <= c < len(verts):
                add_tri(tris, verts[a], verts[b], verts[c])
    return tris


def parse_off_mesh(path: Path) -> List[Triangle]:
    """Parse a simple OFF mesh written by pyg4ometry/pycgal."""

    def next_data_line(handle) -> str:
        while True:
            line = handle.readline()
            if not line:
                return ""
            line = line.strip()
            if line and not line.startswith("#"):
                return line

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header = next_data_line(handle)
        if header == "OFF":
            counts_line = next_data_line(handle)
        elif header.startswith("OFF"):
            counts_line = header[3:].strip()
        else:
            raise ValueError(f"Not an OFF mesh: {header!r}")

        counts = counts_line.split()
        if len(counts) < 2:
            raise ValueError(f"Malformed OFF counts line: {counts_line!r}")
        n_vertices = int(counts[0])
        n_faces = int(counts[1])

        vertices: List[Vec3] = []
        for _ in range(n_vertices):
            parts = next_data_line(handle).split()
            if len(parts) < 3:
                raise ValueError("Malformed OFF vertex line")
            vertices.append((float(parts[0]), float(parts[1]), float(parts[2])))

        tris: List[Triangle] = []
        for _ in range(n_faces):
            parts = next_data_line(handle).split()
            if len(parts) < 4:
                continue
            n = int(parts[0])
            indices = [int(x) for x in parts[1 : 1 + n]]
            if len(indices) < 3:
                continue
            for i in range(1, len(indices) - 1):
                a, b, c = indices[0], indices[i], indices[i + 1]
                if 0 <= a < len(vertices) and 0 <= b < len(vertices) and 0 <= c < len(vertices):
                    add_tri(tris, vertices[a], vertices[b], vertices[c])
    return tris


def extract_triangles_via_off(mesh: object) -> List[Triangle]:
    write_off = getattr(mesh, "writeOff", None)
    if not callable(write_off):
        return []
    fd, name = tempfile.mkstemp(suffix=".off")
    os.close(fd)
    path = Path(name)
    try:
        try:
            write_off(str(path))
        except TypeError:
            write_off(path)
        if not path.exists() or path.stat().st_size == 0:
            return []
        return parse_off_mesh(path)
    finally:
        try:
            path.unlink()
        except OSError:
            pass


def extract_triangles_from_pyg4ometry_mesh(mesh: object) -> List[Triangle]:
    """Best-effort adapter for pyg4ometry/pycsg mesh objects.

    pyg4ometry has changed mesh internals across releases. This adapter accepts
    the common forms: vertex/face arrays, polygon lists, and objects exposing a
    conversion method such as toVerticesAndPolygons().
    """

    for method in ("toVerticesAndPolygons", "toVerticesAndFaces", "toVerticesAndTriangles"):
        converted = maybe_call(mesh, [method])
        if isinstance(converted, tuple) and len(converted) >= 2:
            tris = triangulate_vertex_face_mesh(converted[0], converted[1])
            if tris:
                return tris
        if isinstance(converted, dict):
            vertices = converted.get("vertices") or converted.get("verts") or converted.get("points")
            faces = converted.get("polygons") or converted.get("faces") or converted.get("triangles")
            if vertices is not None and faces is not None:
                tris = triangulate_vertex_face_mesh(vertices, faces)
                if tris:
                    return tris

    vertices = maybe_call(mesh, ["vertices", "verts", "vertexList", "points"])
    faces = maybe_call(mesh, ["faces", "facets", "triangles", "polygons"])
    if vertices is not None and faces is not None:
        tris = triangulate_vertex_face_mesh(vertices, faces)
        if tris:
            return tris

    polygons = maybe_call(mesh, ["polygons"])
    if polygons is None and isinstance(mesh, (list, tuple)):
        polygons = mesh
    tris: List[Triangle] = []
    if polygons is not None:
        for polygon in list(polygons):  # type: ignore[arg-type]
            raw_vertices = maybe_call(polygon, ["vertices", "verts", "points"])
            if raw_vertices is None and isinstance(polygon, (list, tuple)):
                raw_vertices = polygon
            if raw_vertices is None:
                continue
            verts = [p for p in (numeric_triplet(v) for v in list(raw_vertices)) if p is not None]
            for i in range(1, len(verts) - 1):
                add_tri(tris, verts[0], verts[i], verts[i + 1])
    if tris:
        return tris

    return extract_triangles_via_off(mesh)


def short_repr(value: object, limit: int = 240) -> str:
    text = repr(value)
    if len(text) > limit:
        return text[: limit - 3] + "..."
    return text


def public_attrs(value: object, limit: int = 80) -> List[str]:
    return [name for name in dir(value) if not name.startswith("_")][:limit]


def describe_object(label: str, value: object, depth: int = 0) -> List[str]:
    indent = "  " * depth
    lines = [
        f"{indent}{label}: type={type(value).__module__}.{type(value).__name__}",
        f"{indent}{label}: repr={short_repr(value)}",
        f"{indent}{label}: attrs={', '.join(public_attrs(value))}",
    ]
    return lines


def debug_pyg4ometry_mesh(gdml_path: Path, solid_name: Optional[str] = None) -> int:
    try:
        from pyg4ometry import gdml as pyg4_gdml  # type: ignore
    except ModuleNotFoundError:
        print("pyg4ometry is not installed in this Python environment.", file=sys.stderr)
        return 3

    try:
        reader = pyg4_gdml.Reader(str(gdml_path), reduceNISTMaterialsToPredefined=True)
    except TypeError:
        reader = pyg4_gdml.Reader(str(gdml_path))
    registry = reader.getRegistry()
    solids = getattr(registry, "solidDict", {})
    if not solids:
        print("No pyg4ometry solids found in registry.", file=sys.stderr)
        return 4

    if solid_name:
        solid = solids.get(solid_name)
        if solid is None:
            print(f"Solid not found: {solid_name}", file=sys.stderr)
            print("First available solids:")
            for name in list(solids)[:20]:
                print(f"  {name}")
            return 4
    else:
        solid_name, solid = next(iter(solids.items()))

    print(f"Inspecting solid: {solid_name}")
    for line in describe_object("solid", solid):
        print(line)
    mesh_func = getattr(solid, "mesh", None)
    if not callable(mesh_func):
        print("solid.mesh is not callable.")
        return 5

    mesh = mesh_func()
    for line in describe_object("mesh", mesh):
        print(line)
    converted = maybe_call(mesh, ["toVerticesAndPolygons"])
    if converted is not None:
        for line in describe_object("mesh.toVerticesAndPolygons()", converted):
            print(line)
        if isinstance(converted, tuple):
            print(f"toVerticesAndPolygons tuple length: {len(converted)}")
            for i, item in enumerate(converted[:2]):
                try:
                    item_len = len(item)  # type: ignore[arg-type]
                except TypeError:
                    item_len = "unknown"
                print(f"  item {i}: type={type(item).__module__}.{type(item).__name__}, len={item_len}, repr={short_repr(item)}")
    tris = extract_triangles_from_pyg4ometry_mesh(mesh)
    print(f"extract_triangles_from_pyg4ometry_mesh -> {len(tris)} triangles")
    if tris:
        print(f"first triangle: {tris[0]}")
        return 0

    try:
        off_tris = extract_triangles_via_off(mesh)
        print(f"extract_triangles_via_off -> {len(off_tris)} triangles")
        if off_tris:
            print(f"first OFF triangle: {off_tris[0]}")
            return 0
    except Exception as exc:
        print(f"extract_triangles_via_off failed: {exc}")

    polygons = maybe_call(mesh, ["polygons"])
    if polygons is not None:
        try:
            polygon_list = list(polygons)  # type: ignore[arg-type]
        except TypeError:
            polygon_list = []
        print(f"polygons: {len(polygon_list)}")
        if polygon_list:
            poly = polygon_list[0]
            for line in describe_object("polygon[0]", poly):
                print(line)
            raw_vertices = maybe_call(poly, ["vertices", "verts", "points"])
            print(f"polygon[0] vertices attr: {short_repr(raw_vertices)}")
            if raw_vertices is not None:
                try:
                    verts = list(raw_vertices)  # type: ignore[arg-type]
                except TypeError:
                    verts = []
                print(f"polygon[0] vertex count: {len(verts)}")
                if verts:
                    vertex = verts[0]
                    for line in describe_object("vertex[0]", vertex):
                        print(line)
                    print(f"numeric_triplet(vertex[0]) -> {numeric_triplet(vertex)}")
                    if hasattr(vertex, "pos"):
                        pos = getattr(vertex, "pos")
                        for line in describe_object("vertex[0].pos", pos):
                            print(line)
                        print(f"numeric_triplet(vertex[0].pos) -> {numeric_triplet(pos)}")

    vertices = maybe_call(mesh, ["vertices", "verts", "vertexList", "points"])
    faces = maybe_call(mesh, ["faces", "facets", "triangles", "polygons"])
    print(f"mesh vertices attr: {short_repr(vertices)}")
    print(f"mesh faces attr: {short_repr(faces)}")
    return 6


class Pyg4ometrySceneBuilder(GDMLSceneBuilder):
    """Scene builder that asks pyg4ometry to tessellate each GDML solid.

    The XML traversal remains local because it gives stable names/paths and lets
    the same WebGL frontend work without depending on pyg4ometry in the browser.
    """

    def parse(self) -> None:
        self.backend_name = "pyg4ometry"
        super().parse()
        try:
            from pyg4ometry import gdml as pyg4_gdml  # type: ignore
        except ModuleNotFoundError as exc:
            raise RuntimeError(
                "pyg4ometry is not installed. Install it in this environment or use --backend python."
            ) from exc

        try:
            reader = pyg4_gdml.Reader(str(self.gdml), reduceNISTMaterialsToPredefined=True)
        except TypeError:
            reader = pyg4_gdml.Reader(str(self.gdml))
        registry = reader.getRegistry()
        solid_dict = getattr(registry, "solidDict", None)
        if solid_dict is None:
            raise RuntimeError("Could not find pyg4ometry registry.solidDict after reading GDML.")
        self.pyg4ometry_solids = solid_dict

    def solid_mesh(self, ref: str) -> Tuple[List[Triangle], str]:
        key = f"pyg4ometry:{ref}"
        if key in self.mesh_cache:
            return self.mesh_cache[key]
        solid = getattr(self, "pyg4ometry_solids", {}).get(ref)
        if solid is None:
            self.scene.warn(f"pyg4ometry missing solid reference: {ref}")
            return super().solid_mesh(ref)
        mesh_func = getattr(solid, "mesh", None)
        if not callable(mesh_func):
            self.scene.warn(f"pyg4ometry solid has no mesh() method: {ref}")
            return super().solid_mesh(ref)
        try:
            mesh = mesh_func()
            tris = extract_triangles_from_pyg4ometry_mesh(mesh)
        except Exception as exc:
            self.scene.warn(f"pyg4ometry mesh failed for {ref}: {exc}")
            return super().solid_mesh(ref)
        if not tris:
            self.scene.warn(f"pyg4ometry returned no extractable triangles for {ref}")
            return super().solid_mesh(ref)
        kind = f"pyg4ometry/{type(solid).__name__}"
        self.mesh_cache[key] = (tris, kind)
        return self.mesh_cache[key]


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GDML Viewer - {title}</title>
<style>
html, body {{ margin: 0; height: 100%; font-family: system-ui, sans-serif; }}
body {{
  --bg: #f7f8fb; --panel: #ffffff; --text: #1d2430; --muted: #475569;
  --border: #d9dde7; --soft: #f8fafc; --pill: #eef2ff; --canvas: #fbfcff;
  --warn-bg: #fff7ed; --warn-border: #fed7aa; --hud: rgba(255,255,255,.86);
  background: var(--bg); color: var(--text);
}}
body.dark {{
  --bg: #0f172a; --panel: #111827; --text: #e5e7eb; --muted: #a3b1c6;
  --border: #334155; --soft: #1f2937; --pill: #243044; --canvas: #08111f;
  --warn-bg: #2b2115; --warn-border: #7c4a16; --hud: rgba(15,23,42,.86);
}}
#app {{ display: grid; grid-template-columns: minmax(280px, 360px) 1fr; height: 100%; }}
aside {{ border-right: 1px solid var(--border); padding: 14px; overflow: auto; background: var(--panel); }}
main {{ position: relative; min-width: 0; }}
canvas {{ display: block; width: 100%; height: 100%; background: var(--canvas); }}
h1 {{ font-size: 18px; margin: 0 0 10px; }}
h2 {{ font-size: 13px; margin: 18px 0 8px; color: var(--muted); }}
.stat {{ display: grid; grid-template-columns: 1fr auto; gap: 8px; font-size: 13px; margin: 3px 0; }}
input, select, button {{ width: 100%; box-sizing: border-box; margin: 5px 0; color: var(--text); background: var(--soft); border: 1px solid var(--border); border-radius: 6px; padding: 6px; }}
button {{ cursor: pointer; }}
input[type=checkbox] {{ width: auto; }}
.row {{ display: flex; gap: 8px; align-items: center; }}
.row > * {{ flex: 1; }}
.pill {{ display: inline-block; padding: 2px 7px; border-radius: 99px; background: var(--pill); margin: 2px; font-size: 12px; }}
.warn {{ background: var(--warn-bg); border: 1px solid var(--warn-border); padding: 8px; border-radius: 6px; font-size: 12px; max-height: 120px; overflow: auto; }}
#details {{ font-size: 12px; white-space: pre-wrap; background: var(--soft); padding: 8px; border-radius: 6px; max-height: 180px; overflow: auto; }}
#hud {{ position: absolute; left: 12px; bottom: 12px; background: var(--hud); padding: 8px 10px; border-radius: 6px; font-size: 12px; border: 1px solid var(--border); }}
.axis-label {{ position: absolute; display: none; pointer-events: none; font-size: 13px; font-weight: 800; transform: translate(-50%, -50%); text-shadow: 0 0 3px var(--panel), 0 0 3px var(--panel); }}
.material-list {{ display: grid; gap: 4px; max-height: 230px; overflow: auto; padding-right: 4px; }}
.material-item {{ display: grid; grid-template-columns: auto 18px 1fr auto; gap: 7px; align-items: center; font-size: 12px; padding: 3px 0; }}
.swatch {{ width: 14px; height: 14px; border-radius: 3px; border: 1px solid var(--border); }}
.small-actions {{ display: grid; grid-template-columns: 1fr 1fr; gap: 6px; }}
.camera-actions {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 5px; }}
.camera-actions button {{ font-size: 11px; padding: 5px 2px; }}
.axis-key {{ display: flex; gap: 10px; margin: 5px 0; font-size: 12px; font-weight: 700; }}
.axis-x {{ color: #ef4444; }} .axis-y {{ color: #22c55e; }} .axis-z {{ color: #3b82f6; }}
.range-value {{ flex: 0 0 48px; text-align: right; font-variant-numeric: tabular-nums; }}
</style>
</head>
<body>
<div id="app">
<aside>
  <h1>GDML Viewer</h1>
  <div class="stat"><span>Source</span><b id="sourceName"></b></div>
  <div class="stat"><span>Objects</span><b id="objCount"></b></div>
  <div class="stat"><span>Triangles</span><b id="triCount"></b></div>
  <div class="stat"><span>Vertices</span><b id="vertexCount"></b></div>
  <div class="stat"><span>Backend</span><b id="backendName"></b></div>
  <div class="stat"><span>Quality</span><b id="qualityName"></b></div>
  <div class="stat"><span>World</span><b id="worldName"></b></div>
  <h2>View</h2>
  <button id="reset">Reset Camera</button>
  <div class="camera-actions">
    <button id="viewXYPos" title="Look toward the origin along -z">XY +Z</button>
    <button id="viewXYNeg" title="Look toward the origin along +z">XY -Z</button>
    <button id="viewXZPos" title="Look toward the origin along -y">XZ +Y</button>
    <button id="viewXZNeg" title="Look toward the origin along +y">XZ -Y</button>
    <button id="viewYZPos" title="Look toward the origin along -x">YZ +X</button>
    <button id="viewYZNeg" title="Look toward the origin along +x">YZ -X</button>
  </div>
  <label class="row"><span>Opacity</span><input id="opacity" type="range" min="5" max="100" value="78"></label>
  <label><input id="wireframe" type="checkbox"> Wireframe overlay</label>
  <label><input id="coordinateAxes" type="checkbox" checked> Coordinate axes</label>
  <div class="axis-key"><span class="axis-x">X</span><span class="axis-y">Y</span><span class="axis-z">Z</span></div>
  <label><input id="darkMode" type="checkbox"> Dark mode</label>
  <h2>Phi cutout</h2>
  <label><input id="phiCutEnabled" type="checkbox"> Remove azimuthal wedge</label>
  <label class="row"><span>Start phi</span><input id="phiCutStart" type="range" min="-180" max="180" value="-45" step="1"><span id="phiCutStartValue" class="range-value">-45°</span></label>
  <label class="row"><span>Width</span><input id="phiCutWidth" type="range" min="0" max="359" value="90" step="1"><span id="phiCutWidthValue" class="range-value">90°</span></label>
  <div style="font-size:11px;color:var(--muted)">The wedge is measured counter-clockwise about global +z from its start angle.</div>
  <h2>Filter</h2>
  <input id="search" placeholder="Search name, material, solid...">
  <select id="solidType"></select>
  <h2>Materials</h2>
  <div class="small-actions"><button id="allMaterials">All</button><button id="noMaterials">None</button></div>
  <div id="materials"></div>
  <h2>Picked Object</h2>
  <div id="details">Click a point in the detector view.</div>
  <h2>Accuracy Notes</h2>
  <div class="warn" id="warnings"></div>
</aside>
<main>
  <canvas id="gl"></canvas>
  <span id="axisLabelX" class="axis-label axis-x">X</span>
  <span id="axisLabelY" class="axis-label axis-y">Y</span>
  <span id="axisLabelZ" class="axis-label axis-z">Z</span>
  <div id="hud">Drag rotate. Wheel zoom. Shift-drag pan. Click inspect.</div>
</main>
</div>
<script>
const SCENE = {scene_json};
const canvas = document.getElementById('gl');
const gl = canvas.getContext('webgl', {{antialias: true, alpha: false}});
if (!gl) alert('WebGL is not available in this browser.');
document.getElementById('sourceName').textContent = SCENE.meta.sourceFile || (SCENE.meta.source || '').split('/').pop();
document.getElementById('objCount').textContent = SCENE.meta.objects.toLocaleString();
document.getElementById('triCount').textContent = SCENE.meta.triangles.toLocaleString();
document.getElementById('vertexCount').textContent = SCENE.meta.vertices.toLocaleString();
document.getElementById('backendName').textContent = SCENE.meta.backend || 'unknown';
document.getElementById('qualityName').textContent = SCENE.meta.quality || 'unknown';
document.getElementById('worldName').textContent = SCENE.meta.world;
document.getElementById('warnings').innerHTML = SCENE.meta.warnings.length ? SCENE.meta.warnings.map(x => `<div>${{x}}</div>`).join('') : 'No geometry warnings.';

const positions = new Float32Array(SCENE.positions);
const colorsRaw = new Uint8Array(SCENE.colors);
const ids = new Uint32Array(SCENE.ids);
const colors = new Float32Array(colorsRaw.length);
for (let i=0; i<colorsRaw.length; i++) colors[i] = colorsRaw[i] / 255;
const visible = new Float32Array(ids.length);
visible.fill(1);
let objectVisible = SCENE.objects.map(() => 1);
const materialNames = Object.keys(SCENE.meta.materials).sort();
const materialEnabled = Object.fromEntries(materialNames.map(name => [name, true]));

function fillSelect(id, title, values) {{
  const el = document.getElementById(id);
  el.innerHTML = `<option value="">All ${{title}}</option>` + values.map(v => `<option value="${{v}}">${{v}}</option>`).join('');
}}
fillSelect('solidType', 'solid types', Object.keys(SCENE.meta.solidTypes).sort());

function materialColor(name) {{
  const c = (SCENE.meta.materialColors && SCENE.meta.materialColors[name]) || [128,128,128,200];
  return `rgba(${{c[0]}}, ${{c[1]}}, ${{c[2]}}, ${{Math.max(c[3] / 255, 0.35)}})`;
}}
function renderMaterialToggles() {{
  const box = document.getElementById('materials');
  box.className = 'material-list';
  box.innerHTML = materialNames.map(name => `
    <label class="material-item">
      <input type="checkbox" data-material="${{name}}" ${{materialEnabled[name] ? 'checked' : ''}}>
      <span class="swatch" style="background:${{materialColor(name)}}"></span>
      <span title="${{name}}">${{name}}</span>
      <span>${{SCENE.meta.materials[name].toLocaleString()}}</span>
    </label>
  `).join('');
  box.querySelectorAll('input[data-material]').forEach(input => {{
    input.addEventListener('change', () => {{
      materialEnabled[input.dataset.material] = input.checked;
      applyFilters();
    }});
  }});
}}
renderMaterialToggles();

function shader(type, src) {{
  const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
  return s;
}}
const vs = shader(gl.VERTEX_SHADER, `
attribute vec3 a_pos; attribute vec4 a_col; attribute float a_vis;
uniform mat4 u_mvp; uniform float u_opacity;
varying vec4 v_col; varying vec3 v_pos;
void main() {{
  gl_Position = u_mvp * vec4(a_pos, 1.0);
  v_col = vec4(a_col.rgb, a_col.a * u_opacity * a_vis);
  v_pos = a_pos;
}}`);
const fs = shader(gl.FRAGMENT_SHADER, `
precision mediump float; varying vec4 v_col; varying vec3 v_pos;
uniform float u_phi_cut_enabled; uniform float u_phi_cut_start; uniform float u_phi_cut_width;
void main() {{
  if (v_col.a < 0.01) discard;
  if (u_phi_cut_enabled > 0.5 && u_phi_cut_width > 0.0) {{
    const float two_pi = 6.28318530718;
    float relative_phi = mod(atan(v_pos.y, v_pos.x) - u_phi_cut_start + two_pi, two_pi);
    if (relative_phi <= u_phi_cut_width) discard;
  }}
  gl_FragColor = v_col;
}}`);
const prog = gl.createProgram(); gl.attachShader(prog, vs); gl.attachShader(prog, fs); gl.linkProgram(prog); gl.useProgram(prog);
const locPos = gl.getAttribLocation(prog, 'a_pos');
const locCol = gl.getAttribLocation(prog, 'a_col');
const locVis = gl.getAttribLocation(prog, 'a_vis');
function makeBuffer(data, usage=gl.DYNAMIC_DRAW) {{
  const buf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  gl.bufferData(gl.ARRAY_BUFFER, data, usage);
  return buf;
}}
function setAttrib(loc, buf, size) {{
  gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  gl.enableVertexAttribArray(loc);
  gl.vertexAttribPointer(loc, size, gl.FLOAT, false, 0, 0);
}}
function bindArray(name, data, size) {{
  const loc = gl.getAttribLocation(prog, name);
  const buf = makeBuffer(data);
  setAttrib(loc, buf, size);
  return buf;
}}
const surfacePosBuf = bindArray('a_pos', positions, 3);
const surfaceColBuf = bindArray('a_col', colors, 4);
const visibleBuf = bindArray('a_vis', visible, 1);
const uMvp = gl.getUniformLocation(prog, 'u_mvp');
const uOpacity = gl.getUniformLocation(prog, 'u_opacity');
const uPhiCutEnabled = gl.getUniformLocation(prog, 'u_phi_cut_enabled');
const uPhiCutStart = gl.getUniformLocation(prog, 'u_phi_cut_start');
const uPhiCutWidth = gl.getUniformLocation(prog, 'u_phi_cut_width');
let wireBuffers = null;

const bbox = SCENE.meta.bbox;
const center = [(bbox[0]+bbox[3])/2, (bbox[1]+bbox[4])/2, (bbox[2]+bbox[5])/2];
const radius = Math.hypot(bbox[3]-bbox[0], bbox[4]-bbox[1], bbox[5]-bbox[2]) || 1;
const axisLength = radius * 0.18;
const axisPositions = new Float32Array([
  0,0,0, axisLength,0,0,
  0,0,0, 0,axisLength,0,
  0,0,0, 0,0,axisLength,
]);
const axisColors = new Float32Array([
  0.94,0.16,0.16,1, 0.94,0.16,0.16,1,
  0.10,0.78,0.30,1, 0.10,0.78,0.30,1,
  0.12,0.42,0.96,1, 0.12,0.42,0.96,1,
]);
const axisVisible = new Float32Array(6); axisVisible.fill(1);
const axisPosBuf = makeBuffer(axisPositions, gl.STATIC_DRAW);
const axisColBuf = makeBuffer(axisColors, gl.STATIC_DRAW);
const axisVisBuf = makeBuffer(axisVisible, gl.STATIC_DRAW);
let yaw = -0.65, pitch = 0.45, dist = radius * 1.35, panX = 0, panY = 0;
let cameraUp = [0,1,0];

function mat4mul(a,b) {{
  const o = new Float32Array(16);
  for (let r=0;r<4;r++) for (let c=0;c<4;c++) for (let k=0;k<4;k++) o[c*4+r] += a[k*4+r]*b[c*4+k];
  return o;
}}
function perspective(fovy, aspect, near, far) {{
  const f = 1/Math.tan(fovy/2), nf = 1/(near-far);
  return new Float32Array([f/aspect,0,0,0, 0,f,0,0, 0,0,(far+near)*nf,-1, 0,0,2*far*near*nf,0]);
}}
function lookAt(eye, target, up) {{
  let zx=eye[0]-target[0], zy=eye[1]-target[1], zz=eye[2]-target[2];
  let zl=Math.hypot(zx,zy,zz); zx/=zl; zy/=zl; zz/=zl;
  let xx=up[1]*zz-up[2]*zy, xy=up[2]*zx-up[0]*zz, xz=up[0]*zy-up[1]*zx;
  let xl=Math.hypot(xx,xy,xz); xx/=xl; xy/=xl; xz/=xl;
  let yx=zy*xz-zz*xy, yy=zz*xx-zx*xz, yz=zx*xy-zy*xx;
  return new Float32Array([xx,yx,zx,0, xy,yy,zy,0, xz,yz,zz,0, -(xx*eye[0]+xy*eye[1]+xz*eye[2]), -(yx*eye[0]+yy*eye[1]+yz*eye[2]), -(zx*eye[0]+zy*eye[1]+zz*eye[2]), 1]);
}}
function resize() {{
  const dpr = window.devicePixelRatio || 1;
  const w = Math.floor(canvas.clientWidth*dpr), h = Math.floor(canvas.clientHeight*dpr);
  if (canvas.width !== w || canvas.height !== h) {{ canvas.width = w; canvas.height = h; gl.viewport(0,0,w,h); }}
}}
function themeIsDark() {{ return document.body.classList.contains('dark'); }}
function clearCanvas() {{
  if (themeIsDark()) gl.clearColor(0.031,0.067,0.122,1);
  else gl.clearColor(0.985,0.99,1,1);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
}}
function phiCutRadians() {{
  return {{
    enabled: document.getElementById('phiCutEnabled').checked,
    start: +document.getElementById('phiCutStart').value * Math.PI / 180,
    width: +document.getElementById('phiCutWidth').value * Math.PI / 180,
  }};
}}
function pointIsPhiCut(x, y) {{
  const cut = phiCutRadians();
  if (!cut.enabled || cut.width <= 0) return false;
  const twoPi = 2*Math.PI;
  const relative = ((Math.atan2(y, x) - cut.start) % twoPi + twoPi) % twoPi;
  return relative <= cut.width;
}}
function updateAxisLabels(mvp, enabled) {{
  const endpoints = [
    ['axisLabelX', axisLength, 0, 0],
    ['axisLabelY', 0, axisLength, 0],
    ['axisLabelZ', 0, 0, axisLength],
  ];
  for (const [id, x, y, z] of endpoints) {{
    const label = document.getElementById(id);
    const tx=mvp[0]*x+mvp[4]*y+mvp[8]*z+mvp[12];
    const ty=mvp[1]*x+mvp[5]*y+mvp[9]*z+mvp[13];
    const tw=mvp[3]*x+mvp[7]*y+mvp[11]*z+mvp[15];
    if (!enabled || tw <= 0) {{ label.style.display='none'; continue; }}
    label.style.display='block';
    label.style.left=((tx/tw*0.5+0.5)*canvas.clientWidth)+'px';
    label.style.top=((-ty/tw*0.5+0.5)*canvas.clientHeight)+'px';
  }}
}}
function buildWireBuffers() {{
  if (wireBuffers) return wireBuffers;
  const triCount = positions.length / 9;
  const linePositions = new Float32Array(triCount * 18);
  const lineColors = new Float32Array(triCount * 24);
  const lineIds = new Uint32Array(triCount * 6);
  const lineVisible = new Float32Array(triCount * 6);
  let p = 0, c = 0, n = 0;
  for (let i = 0; i < positions.length; i += 9) {{
    const objId = ids[i / 3];
    const edges = [0,1, 1,2, 2,0];
    for (let e = 0; e < edges.length; e++) {{
      const v = i + edges[e] * 3;
      linePositions[p++] = positions[v];
      linePositions[p++] = positions[v + 1];
      linePositions[p++] = positions[v + 2];
      if (themeIsDark()) {{ lineColors[c++] = 1; lineColors[c++] = 1; lineColors[c++] = 1; lineColors[c++] = 0.42; }}
      else {{ lineColors[c++] = 0.02; lineColors[c++] = 0.04; lineColors[c++] = 0.08; lineColors[c++] = 0.35; }}
      lineIds[n] = objId;
      lineVisible[n] = objectVisible[objId];
      n++;
    }}
  }}
  wireBuffers = {{
    positions: linePositions,
    colors: lineColors,
    ids: lineIds,
    visible: lineVisible,
    posBuf: makeBuffer(linePositions, gl.STATIC_DRAW),
    colBuf: makeBuffer(lineColors),
    visBuf: makeBuffer(lineVisible),
  }};
  return wireBuffers;
}}
function updateWireColors() {{
  if (!wireBuffers) return;
  for (let i = 0; i < wireBuffers.colors.length; i += 4) {{
    if (themeIsDark()) {{
      wireBuffers.colors[i] = 1; wireBuffers.colors[i + 1] = 1; wireBuffers.colors[i + 2] = 1; wireBuffers.colors[i + 3] = 0.42;
    }} else {{
      wireBuffers.colors[i] = 0.02; wireBuffers.colors[i + 1] = 0.04; wireBuffers.colors[i + 2] = 0.08; wireBuffers.colors[i + 3] = 0.35;
    }}
  }}
  gl.bindBuffer(gl.ARRAY_BUFFER, wireBuffers.colBuf);
  gl.bufferData(gl.ARRAY_BUFFER, wireBuffers.colors, gl.DYNAMIC_DRAW);
}}
function updateWireVisibility() {{
  if (!wireBuffers) return;
  for (let i = 0; i < wireBuffers.ids.length; i++) wireBuffers.visible[i] = objectVisible[wireBuffers.ids[i]];
  gl.bindBuffer(gl.ARRAY_BUFFER, wireBuffers.visBuf);
  gl.bufferData(gl.ARRAY_BUFFER, wireBuffers.visible, gl.DYNAMIC_DRAW);
}}
function render() {{
  resize();
  clearCanvas();
  gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); gl.disable(gl.DEPTH_TEST);
  const cx = center[0] + panX, cy = center[1] + panY, cz = center[2];
  const eye = [cx + dist*Math.cos(pitch)*Math.sin(yaw), cy + dist*Math.sin(pitch), cz + dist*Math.cos(pitch)*Math.cos(yaw)];
  const proj = perspective(45*Math.PI/180, canvas.width/canvas.height, radius/500, radius*20);
  const view = lookAt(eye, [cx,cy,cz], cameraUp);
  const mvp = mat4mul(proj, view);
  gl.uniformMatrix4fv(uMvp, false, mvp);
  setAttrib(locPos, surfacePosBuf, 3);
  setAttrib(locCol, surfaceColBuf, 4);
  setAttrib(locVis, visibleBuf, 1);
  gl.uniform1f(uOpacity, +document.getElementById('opacity').value/100);
  const phiCut = phiCutRadians();
  gl.uniform1f(uPhiCutEnabled, phiCut.enabled ? 1 : 0);
  gl.uniform1f(uPhiCutStart, phiCut.start);
  gl.uniform1f(uPhiCutWidth, phiCut.width);
  gl.drawArrays(gl.TRIANGLES, 0, positions.length/3);
  if (document.getElementById('wireframe').checked) {{
    const wire = buildWireBuffers();
    setAttrib(locPos, wire.posBuf, 3);
    setAttrib(locCol, wire.colBuf, 4);
    setAttrib(locVis, wire.visBuf, 1);
    gl.uniform1f(uOpacity, 1.0);
    gl.drawArrays(gl.LINES, 0, wire.positions.length/3);
  }}
  const axesEnabled = document.getElementById('coordinateAxes').checked;
  if (axesEnabled) {{
    setAttrib(locPos, axisPosBuf, 3);
    setAttrib(locCol, axisColBuf, 4);
    setAttrib(locVis, axisVisBuf, 1);
    gl.uniform1f(uOpacity, 1.0);
    gl.uniform1f(uPhiCutEnabled, 0.0);
    gl.drawArrays(gl.LINES, 0, axisPositions.length/3);
  }}
  updateAxisLabels(mvp, axesEnabled);
  requestAnimationFrame(render);
}}
render();

let drag = null;
canvas.addEventListener('mousedown', e => drag = {{x:e.clientX, y:e.clientY, shift:e.shiftKey}});
window.addEventListener('mouseup', () => drag = null);
window.addEventListener('mousemove', e => {{
  if (!drag) return;
  const dx=e.clientX-drag.x, dy=e.clientY-drag.y; drag.x=e.clientX; drag.y=e.clientY;
  if (drag.shift) {{ panX -= dx*dist/900; panY += dy*dist/900; }}
  else {{ yaw += dx*0.006; pitch = Math.max(-1.55, Math.min(1.55, pitch + dy*0.006)); cameraUp = [0,1,0]; }}
}});
canvas.addEventListener('wheel', e => {{ e.preventDefault(); dist = Math.max(radius*0.01, Math.min(radius*50, dist*Math.exp(e.deltaY*0.001))); }}, {{passive:false}});
function setCamera(newYaw, newPitch, up) {{
  yaw=newYaw; pitch=newPitch; cameraUp=up; dist=radius*1.35; panX=0; panY=0;
}}
document.getElementById('reset').onclick = () => setCamera(-0.65, 0.45, [0,1,0]);
document.getElementById('viewXYPos').onclick = () => setCamera(0, 0, [0,1,0]);
document.getElementById('viewXYNeg').onclick = () => setCamera(Math.PI, 0, [0,1,0]);
document.getElementById('viewXZPos').onclick = () => setCamera(0, Math.PI/2, [0,0,1]);
document.getElementById('viewXZNeg').onclick = () => setCamera(0, -Math.PI/2, [0,0,1]);
document.getElementById('viewYZPos').onclick = () => setCamera(Math.PI/2, 0, [0,0,1]);
document.getElementById('viewYZNeg').onclick = () => setCamera(-Math.PI/2, 0, [0,0,1]);
['phiCutStart','phiCutWidth'].forEach(id => document.getElementById(id).addEventListener('input', () => {{
  document.getElementById(id + 'Value').textContent = document.getElementById(id).value + '°';
}}));
document.getElementById('wireframe').addEventListener('change', () => {{ if (document.getElementById('wireframe').checked) buildWireBuffers(); }});
document.getElementById('darkMode').addEventListener('change', e => {{
  document.body.classList.toggle('dark', e.target.checked);
  updateWireColors();
}});
document.getElementById('allMaterials').onclick = () => {{
  materialNames.forEach(name => materialEnabled[name] = true);
  renderMaterialToggles();
  applyFilters();
}};
document.getElementById('noMaterials').onclick = () => {{
  materialNames.forEach(name => materialEnabled[name] = false);
  renderMaterialToggles();
  applyFilters();
}};

function applyFilters() {{
  const q = document.getElementById('search').value.toLowerCase();
  const solid = document.getElementById('solidType').value;
  objectVisible = SCENE.objects.map(o => (materialEnabled[o.material] !== false) && (!solid || o.solidType === solid) && (!q || JSON.stringify(o).toLowerCase().includes(q)) ? 1 : 0);
  for (let i=0; i<ids.length; i++) visible[i] = objectVisible[ids[i]];
  gl.bindBuffer(gl.ARRAY_BUFFER, visibleBuf); gl.bufferData(gl.ARRAY_BUFFER, visible, gl.DYNAMIC_DRAW);
  updateWireVisibility();
}}
['search','solidType'].forEach(id => document.getElementById(id).addEventListener('input', applyFilters));

canvas.addEventListener('click', e => {{
  const rect = canvas.getBoundingClientRect(), x = e.clientX-rect.left, y = e.clientY-rect.top;
  let best=-1, bestD=14*14;
  const aspect = canvas.width/canvas.height;
  const cx = center[0] + panX, cy = center[1] + panY, cz = center[2];
  const eye = [cx + dist*Math.cos(pitch)*Math.sin(yaw), cy + dist*Math.sin(pitch), cz + dist*Math.cos(pitch)*Math.cos(yaw)];
  const mvp = mat4mul(perspective(45*Math.PI/180, aspect, radius/500, radius*20), lookAt(eye, [cx,cy,cz], cameraUp));
  for (let i=0; i<positions.length; i+=30) {{
    const id = ids[i/3]; if (!visible[i/3]) continue;
    const px=positions[i], py=positions[i+1], pz=positions[i+2];
    if (pointIsPhiCut(px, py)) continue;
    const tx=mvp[0]*px+mvp[4]*py+mvp[8]*pz+mvp[12], ty=mvp[1]*px+mvp[5]*py+mvp[9]*pz+mvp[13], tw=mvp[3]*px+mvp[7]*py+mvp[11]*pz+mvp[15];
    if (tw <= 0) continue;
    const sx=(tx/tw*0.5+0.5)*rect.width, sy=(-ty/tw*0.5+0.5)*rect.height;
    const d=(sx-x)*(sx-x)+(sy-y)*(sy-y);
    if (d < bestD) {{ bestD=d; best=id; }}
  }}
  if (best >= 0) document.getElementById('details').textContent = JSON.stringify(SCENE.objects[best], null, 2);
}});
</script>
</body>
</html>
"""


def write_html(scene: Dict[str, object], output: Path) -> None:
    title = html.escape(Path(str(scene["meta"]["source"])).name)  # type: ignore[index]
    output.write_text(
        HTML_TEMPLATE.format(title=title, scene_json=json.dumps(scene, separators=(",", ":"))),
        encoding="utf-8",
    )


def serve_file(path: Path, port: int) -> None:
    os.chdir(path.parent)
    server = ThreadingHTTPServer(("127.0.0.1", port), SimpleHTTPRequestHandler)
    url = f"http://127.0.0.1:{port}/{path.name}"
    print(f"Serving {url}")
    try:
        webbrowser.open(url)
    except Exception:
        pass
    server.serve_forever()


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a standalone HTML/WebGL viewer for GDML geometry.")
    parser.add_argument("gdml", type=Path, help="Input GDML file")
    parser.add_argument("--out", type=Path, default=None, help="Output HTML path")
    parser.add_argument(
        "--backend",
        choices=["auto", "python", "pyg4ometry"],
        default="auto",
        help="Geometry tessellation backend. 'pyg4ometry' gives more accurate boolean solids when installed.",
    )
    parser.add_argument("--quality", choices=sorted(QUALITY_SEGMENTS), default="medium")
    parser.add_argument("--include-containers", action="store_true", help="Draw container/logical volumes too")
    parser.add_argument("--include-world", action="store_true", help="Draw the world volume")
    parser.add_argument("--max-objects", type=int, default=100_000)
    parser.add_argument("--max-vertices", type=int, default=5_000_000)
    parser.add_argument(
        "--debug-pyg4ometry-mesh",
        action="store_true",
        help="Print pyg4ometry mesh internals for one solid and exit.",
    )
    parser.add_argument(
        "--debug-solid",
        default="",
        help="Solid name to inspect with --debug-pyg4ometry-mesh.",
    )
    parser.add_argument("--serve", action="store_true", help="Serve the HTML on localhost after writing it")
    parser.add_argument("--port", type=int, default=8765)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    if not args.gdml.exists():
        print(f"GDML file not found: {args.gdml}", file=sys.stderr)
        return 2
    if args.debug_pyg4ometry_mesh:
        return debug_pyg4ometry_mesh(args.gdml, args.debug_solid or None)

    output = args.out or (Path(tempfile.gettempdir()) / f"{args.gdml.stem}.gdml_viewer.html")
    builder_cls = GDMLSceneBuilder
    if args.backend in {"auto", "pyg4ometry"}:
        try:
            import pyg4ometry  # type: ignore  # noqa: F401
            builder_cls = Pyg4ometrySceneBuilder
            print("Using pyg4ometry backend for solid tessellation.")
        except ModuleNotFoundError:
            if args.backend == "pyg4ometry":
                print(
                    "pyg4ometry is not installed in this Python environment. "
                    "Install it or rerun with --backend python.",
                    file=sys.stderr,
                )
                return 3
            print("pyg4ometry not installed; using pure-Python backend.")

    builder = builder_cls(
        args.gdml,
        quality=args.quality,
        include_containers=args.include_containers,
        include_world=args.include_world,
        max_objects=args.max_objects,
        max_vertices=args.max_vertices,
    )
    scene = builder.build()
    write_html(scene, output)
    meta = scene["meta"]  # type: ignore[index]
    print(f"Wrote {output}")
    print(f"Objects: {meta['objects']:,}  Triangles: {meta['triangles']:,}  Vertices: {meta['vertices']:,}")  # type: ignore[index]
    warnings = meta.get("warnings", [])  # type: ignore[union-attr]
    if warnings:
        print("Warnings:")
        for item in warnings[:8]:
            print(f"  - {item}")
        if len(warnings) > 8:
            print(f"  - ... {len(warnings) - 8} more")
    if args.serve:
        serve_file(output, args.port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
