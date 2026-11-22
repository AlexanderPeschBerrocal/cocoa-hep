#!/usr/bin/env python3
"""
Lazy-loading GDML geometry viewer.

This is a scalable companion to ``gdml_viewer.py``.  It reuses that module's
GDML parsing and solid tessellation, but writes a small HTML application plus
separate binary geometry chunks grouped by detector subsystem and layer.
Chunks are fetched only when their layer is enabled in the browser.

The viewer also implements a shader-based phi cut.  The cut removes (or keeps)
an angular wedge without rebuilding the geometry, which makes tracker and
calorimeter cross-section views inexpensive.

Typical use::

    python gdml_viewer_lazy.py detector.gdml --serve

The generated HTML uses ``fetch()`` for its binary chunks, so open it through
``--serve`` or another HTTP server rather than directly through ``file://``.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import re
import struct
import sys
import tempfile
import webbrowser
from array import array
from collections import Counter, OrderedDict
from dataclasses import dataclass, field
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

try:
    import gdml_viewer as base
except ImportError as exc:  # pragma: no cover - depends on invocation location
    raise SystemExit(
        "gdml_viewer_lazy.py must be placed next to gdml_viewer.py (or that "
        "directory must be on PYTHONPATH)."
    ) from exc


MAGIC = b"GDMLCHN1"
FORMAT_VERSION = 1
HEADER = struct.Struct("<8s6I")
SUBSYSTEM_ORDER = {
    "Tracker": 0,
    "ECAL": 1,
    "Calo support": 2,
    "HCAL": 3,
    "Other": 4,
}


def safe_name(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_.")
    return text[:120] or "chunk"


def canonical_material(name: str) -> str:
    return base.semantic_name(name)


def bbox_union(target_min: List[float], target_max: List[float], lo: Sequence[float], hi: Sequence[float]) -> None:
    for axis in range(3):
        target_min[axis] = min(target_min[axis], float(lo[axis]))
        target_max[axis] = max(target_max[axis], float(hi[axis]))


def classify_chunk(
    volume: base.VolumeDef,
    path: str,
    lo: Sequence[float],
    hi: Sequence[float],
    radial_min: float,
    radial_max: float,
) -> Tuple[str, str]:
    """Return a human-readable subsystem and layer for one placed volume."""

    semantic_volume = base.semantic_name(volume.name)
    semantic_solid = base.semantic_name(volume.solid)
    semantic_path = "/".join(base.semantic_name(part) for part in path.split("/"))
    material = canonical_material(volume.material)
    text = f"{semantic_volume} {semantic_solid} {semantic_path} {material}".lower()

    def calo_layer(prefix: str) -> str:
        match = re.search(rf"{prefix.lower()}(\d+)(?:_(\d+))?", text)
        index = match.group(1) if match else "?"
        if match and match.group(2):
            index += f".{match.group(2)}"
        region = "Endcap" if "endcap" in text else "Barrel"
        side = ""
        if "forward" in text:
            side = " +z"
        elif "back" in text or "backward" in text:
            side = " -z"
        return f"{region}{side} layer {index}"

    if "ecal" in text:
        return "ECAL", calo_layer("ecal")
    if "hcal" in text:
        return "HCAL", calo_layer("hcal")
    if "irongap" in text or "iron_gap" in text:
        region = "Endcap" if "endcap" in text else "Barrel"
        side = ""
        if "forward" in text:
            side = " +z"
        elif "back" in text or "backward" in text:
            side = " -z"
        return "Calo support", f"{region}{side} iron gap"

    tracker_tokens = ("inner", "trk", "pixel", "pix", "strip", "silicon", "g4_si")
    if any(token in text for token in tracker_tokens):
        z_span = float(hi[2]) - float(lo[2])
        radial_span = max(radial_max - radial_min, 1.0e-9)
        z_center = 0.5 * (float(lo[2]) + float(hi[2]))
        radial_center = 0.5 * (radial_min + radial_max)
        is_endcap = "endcap" in text or (abs(z_center) > 1.0 and z_span < 4.0 * radial_span)
        if is_endcap:
            side = "+z" if z_center >= 0.0 else "-z"
            layer = f"Endcap {side} |z|≈{round(abs(z_center))} mm"
        else:
            layer = f"Barrel r≈{round(radial_center)} mm"
        return "Tracker", layer

    label = base.semantic_name(volume.name)
    if material.lower().startswith(("galactic", "g4_galactic")):
        label = "Vacuum/container"
    return "Other", label


@dataclass
class Chunk:
    subsystem: str
    layer: str
    part: int
    positions: array = field(default_factory=lambda: array("f"))
    colors: array = field(default_factory=lambda: array("B"))
    object_count: int = 0
    bbox_min: List[float] = field(default_factory=lambda: [float("inf")] * 3)
    bbox_max: List[float] = field(default_factory=lambda: [float("-inf")] * 3)

    @property
    def vertex_count(self) -> int:
        return len(self.positions) // 3

    @property
    def triangle_count(self) -> int:
        return self.vertex_count // 3

    @property
    def base_key(self) -> str:
        return f"{self.subsystem}\x1f{self.layer}"

    @property
    def key(self) -> str:
        return f"{self.base_key}\x1f{self.part}"


class LazyBuilderMixin:
    """Mixin that changes the original builder from one scene to layer chunks."""

    def __init__(self, *args, chunk_vertices: int, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.chunk_vertices = max(3, int(chunk_vertices))
        self.chunks: "OrderedDict[str, Chunk]" = OrderedDict()
        self.chunk_parts: Dict[str, List[str]] = {}
        self.total_vertices = 0
        self.total_objects = 0
        self.material_counts: Counter[str] = Counter()
        self.global_bbox_min = [float("inf")] * 3
        self.global_bbox_max = [float("-inf")] * 3
        self.stopped_at_limit = False

    def _chunk_for(self, subsystem: str, layer: str, incoming_vertices: int) -> Chunk:
        base_key = f"{subsystem}\x1f{layer}"
        keys = self.chunk_parts.setdefault(base_key, [])
        if keys:
            candidate = self.chunks[keys[-1]]
            if candidate.vertex_count + incoming_vertices <= self.chunk_vertices:
                return candidate
        part = len(keys) + 1
        chunk = Chunk(subsystem=subsystem, layer=layer, part=part)
        self.chunks[chunk.key] = chunk
        keys.append(chunk.key)
        return chunk

    def add_volume(self, volume: base.VolumeDef, matrix: base.Mat4, path: str) -> None:
        tris, _ = self.solid_mesh(volume.solid)
        if not tris:
            return

        positions = array("f")
        lo = [float("inf")] * 3
        hi = [float("-inf")] * 3
        radial_min = float("inf")
        radial_max = 0.0
        for tri in tris:
            for point in tri:
                transformed = base.transform_point(matrix, point)
                positions.extend(transformed)
                radius = math.hypot(transformed[0], transformed[1])
                radial_min = min(radial_min, radius)
                radial_max = max(radial_max, radius)
                for axis in range(3):
                    lo[axis] = min(lo[axis], transformed[axis])
                    hi[axis] = max(hi[axis], transformed[axis])

        vertex_count = len(positions) // 3
        if vertex_count == 0:
            return
        if self.total_objects >= self.max_objects or self.total_vertices + vertex_count > self.max_vertices:
            self.stopped_at_limit = True
            return

        subsystem, layer = classify_chunk(volume, path, lo, hi, radial_min, radial_max)
        chunk = self._chunk_for(subsystem, layer, vertex_count)
        material = canonical_material(volume.material) or "(none)"
        rgba = base.color_for_material(volume.material)

        chunk.positions.extend(positions)
        chunk.colors.extend(array("B", rgba) * vertex_count)
        chunk.object_count += 1
        bbox_union(chunk.bbox_min, chunk.bbox_max, lo, hi)
        bbox_union(self.global_bbox_min, self.global_bbox_max, lo, hi)
        self.material_counts[material] += 1
        self.total_vertices += vertex_count
        self.total_objects += 1

    def traverse(self, volume: base.VolumeDef, matrix: base.Mat4, path: str) -> None:
        if self.total_objects >= self.max_objects or self.total_vertices >= self.max_vertices:
            self.stopped_at_limit = True
            return
        placements = base.children(volume.element, "physvol")
        is_container = bool(placements) or volume.material.startswith("Galactic")
        should_render = (self.include_containers or not is_container) and (
            self.include_world or volume.name != self.world
        )
        if should_render:
            self.add_volume(volume, matrix, path)
        for placement in placements:
            ref_el = base.child(placement, "volumeref")
            if ref_el is None:
                continue
            ref = ref_el.attrib.get("ref", "")
            child_volume = self.volumes.get(ref)
            if child_volume is None:
                self.scene.warn(f"Missing volume reference: {ref}")
                continue
            placement_name = placement.attrib.get("name", ref)
            child_matrix = base.matmul(matrix, base.placement_matrix(placement))
            self.traverse(child_volume, child_matrix, f"{path}/{placement_name}")
            if self.stopped_at_limit:
                return

    def build_chunks(self) -> None:
        self.parse()
        world = self.volumes.get(self.world)
        if world is None:
            raise ValueError(f"World volume {self.world!r} was not found")
        self.traverse(world, base.identity(), self.world)
        if self.global_bbox_min[0] == float("inf"):
            self.global_bbox_min = [-1.0, -1.0, -1.0]
            self.global_bbox_max = [1.0, 1.0, 1.0]
        if self.stopped_at_limit:
            self.scene.warn("Stopped at configured object/vertex limit")


class LazyGDMLSceneBuilder(LazyBuilderMixin, base.GDMLSceneBuilder):
    pass


class LazyPyg4ometrySceneBuilder(LazyBuilderMixin, base.Pyg4ometrySceneBuilder):
    pass


def little_endian_bytes(values: array) -> bytes:
    if sys.byteorder == "little":
        return values.tobytes()
    copied = array(values.typecode, values)
    copied.byteswap()
    return copied.tobytes()


def write_chunk(chunk: Chunk, path: Path) -> int:
    positions = little_endian_bytes(chunk.positions)
    colors = little_endian_bytes(chunk.colors)
    position_offset = HEADER.size
    color_offset = position_offset + len(positions)
    id_offset = color_offset + len(colors)
    total_size = id_offset
    header = HEADER.pack(
        MAGIC,
        FORMAT_VERSION,
        chunk.vertex_count,
        position_offset,
        color_offset,
        id_offset,
        total_size,
    )
    path.write_bytes(header + positions + colors)
    return total_size


def chunk_sort_key(chunk: Chunk) -> Tuple[int, str, int]:
    return SUBSYSTEM_ORDER.get(chunk.subsystem, 99), chunk.layer, chunk.part


def write_assets(
    builder: LazyBuilderMixin,
    html_path: Path,
    default_load: Sequence[str],
) -> dict:
    asset_dir = html_path.parent / f"{html_path.stem}.assets"
    asset_dir.mkdir(parents=True, exist_ok=True)
    chunk_entries = []
    used_names: Counter[str] = Counter()
    for chunk in sorted(builder.chunks.values(), key=chunk_sort_key):
        stem = safe_name(f"{chunk.subsystem}_{chunk.layer}")
        used_names[stem] += 1
        suffix = used_names[stem]
        file_stem = f"{stem}.part{suffix:03d}"
        binary_name = f"{file_stem}.bin"
        byte_size = write_chunk(chunk, asset_dir / binary_name)
        chunk_entries.append(
            {
                "key": chunk.key,
                "baseKey": chunk.base_key,
                "subsystem": chunk.subsystem,
                "layer": chunk.layer,
                "part": chunk.part,
                "url": f"{asset_dir.name}/{binary_name}",
                "objects": chunk.object_count,
                "vertices": chunk.vertex_count,
                "triangles": chunk.triangle_count,
                "bytes": byte_size,
                "bbox": [*chunk.bbox_min, *chunk.bbox_max],
            }
        )

    material_colors = {
        name: list(base.color_for_material(name)) for name in builder.material_counts
    }
    return {
        "format": "gdml-lazy-viewer-v1",
        "source": str(builder.gdml),
        "sourceFile": builder.gdml.name,
        "backend": builder.backend_name,
        "quality": builder.quality,
        "world": builder.world,
        "objects": builder.total_objects,
        "vertices": builder.total_vertices,
        "triangles": builder.total_vertices // 3,
        "bbox": [*builder.global_bbox_min, *builder.global_bbox_max],
        "chunks": chunk_entries,
        "materials": dict(sorted(builder.material_counts.items())),
        "materialColors": material_colors,
        "warnings": [
            f"{message} (x{count})" if count > 1 else message
            for message, count in sorted(builder.scene.warnings.items())
        ],
        "defaultLoad": list(default_load),
    }


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Lazy GDML Viewer — __TITLE__</title>
<style>
:root { color-scheme: light dark; font-family: Inter, system-ui, sans-serif; }
* { box-sizing: border-box; }
html, body { margin: 0; width: 100%; height: 100%; overflow: hidden; }
body { display: grid; grid-template-columns: 360px 1fr; background: #07111f; color: #e6edf7; }
aside { overflow: auto; padding: 14px; background: #0e1b2b; border-right: 1px solid #2b3b50; }
main { position: relative; min-width: 0; }
canvas { display: block; width: 100%; height: 100%; cursor: grab; }
canvas:active { cursor: grabbing; }
h1 { font-size: 17px; margin: 0 0 4px; }
h2 { font-size: 13px; margin: 16px 0 7px; color: #9fc5ff; text-transform: uppercase; letter-spacing: .05em; }
.muted { color: #9cafc6; font-size: 12px; line-height: 1.4; }
.stats { display: grid; grid-template-columns: 1fr 1fr; gap: 6px; margin: 10px 0; }
.stat { background: #13243a; border-radius: 7px; padding: 7px 9px; }
.stat span { display: block; color: #8fa6c0; font-size: 10px; }
.stat b { font-size: 13px; }
.row { display: flex; gap: 7px; align-items: center; margin: 6px 0; }
.row label { flex: 1; font-size: 12px; }
.camera-actions { display: grid; grid-template-columns: repeat(3, 1fr); gap: 5px; }
.camera-actions button { font-size: 11px; padding: 5px 2px; }
.axis-key { display: flex; gap: 10px; margin: 5px 0; font-size: 12px; font-weight: 700; }
.axis-x { color: #ef4444; }
.axis-y { color: #22c55e; }
.axis-z { color: #3b82f6; }
.axis-label { position: absolute; display: none; pointer-events: none; font-size: 13px; font-weight: 800; transform: translate(-50%, -50%); text-shadow: 0 0 3px #07111f, 0 0 3px #07111f; }
button, select, input { accent-color: #65a5ff; }
button { border: 1px solid #3b526e; color: #e9f2ff; background: #172c46; border-radius: 6px; padding: 6px 9px; cursor: pointer; }
button:hover { background: #213c5e; }
input[type=range] { width: 100%; }
.subsystem { margin: 7px 0; border: 1px solid #2c4059; border-radius: 7px; overflow: hidden; }
.subsystem-head { display: flex; align-items: center; gap: 7px; padding: 7px 8px; background: #152942; font-weight: 650; font-size: 13px; }
.layers { padding: 4px 8px 7px 24px; }
.layer { display: grid; grid-template-columns: auto 1fr auto; align-items: center; gap: 7px; padding: 3px 0; font-size: 11px; }
.layer small { color: #8fa6c0; }
.loading { color: #f4c76b; }
.loaded { color: #74d99f; }
.error { color: #ff8585; white-space: pre-wrap; }
.material { display: grid; grid-template-columns: 12px 1fr auto; gap: 7px; align-items: center; font-size: 11px; margin: 3px 0; }
.swatch { width: 11px; height: 11px; border-radius: 2px; }
#overlay { position: absolute; right: 12px; top: 12px; background: #091626dd; border: 1px solid #314761; border-radius: 8px; padding: 8px 10px; font-size: 11px; pointer-events: none; }
#warnings { max-height: 100px; overflow: auto; }
</style>
</head>
<body>
<aside>
  <h1>Lazy GDML Viewer</h1>
  <div class="muted" id="source"></div>
  <div class="stats">
    <div class="stat"><span>Total objects</span><b id="totalObjects"></b></div>
    <div class="stat"><span>Total triangles</span><b id="totalTriangles"></b></div>
    <div class="stat"><span>Loaded chunks</span><b id="loadedChunks">0</b></div>
    <div class="stat"><span>Loaded binary</span><b id="loadedBytes">0 B</b></div>
  </div>
  <div class="row">
    <button id="loadAll">Load all</button>
    <button id="hideAll">Hide all</button>
    <button id="resetView">Reset view</button>
  </div>
  <div class="row">
    <label><input id="releaseHidden" type="checkbox" checked> Unload hidden layers</label>
  </div>

  <h2>View</h2>
  <div class="camera-actions">
    <button id="viewXYPos" title="Look toward the origin along -z">XY +Z</button>
    <button id="viewXYNeg" title="Look toward the origin along +z">XY -Z</button>
    <button id="viewXZPos" title="Look toward the origin along -y">XZ +Y</button>
    <button id="viewXZNeg" title="Look toward the origin along +y">XZ -Y</button>
    <button id="viewYZPos" title="Look toward the origin along -x">YZ +X</button>
    <button id="viewYZNeg" title="Look toward the origin along +x">YZ -X</button>
  </div>
  <div class="row">
    <label><input id="coordinateAxes" type="checkbox" checked> Coordinate axes</label>
  </div>
  <div class="axis-key"><span class="axis-x">X</span><span class="axis-y">Y</span><span class="axis-z">Z</span></div>

  <h2>Subsystems and layers</h2>
  <div id="geometryTree"></div>

  <h2>Phi cutout</h2>
  <div class="row"><label><input id="cutEnabled" type="checkbox"> Enable phi cut</label></div>
  <div class="row"><label>Mode</label><select id="cutMode"><option value="remove">Remove wedge</option><option value="keep">Keep wedge only</option></select></div>
  <label class="muted">Wedge start: <b id="cutStartValue">0°</b></label>
  <input id="cutStart" type="range" min="-180" max="180" value="0" step="1">
  <label class="muted">Wedge width: <b id="cutWidthValue">90°</b></label>
  <input id="cutWidth" type="range" min="1" max="359" value="90" step="1">
  <div class="row"><label>Opacity</label><input id="opacity" type="range" min="10" max="100" value="100"></div>

  <h2>Materials</h2>
  <div id="materials"></div>
  <h2>Warnings</h2>
  <div id="warnings" class="muted"></div>
  <h2>Status</h2>
  <div id="status" class="muted">Select a layer to load it.</div>
</aside>
<main>
  <canvas id="canvas"></canvas>
  <span id="axisLabelX" class="axis-label axis-x">X</span>
  <span id="axisLabelY" class="axis-label axis-y">Y</span>
  <span id="axisLabelZ" class="axis-label axis-z">Z</span>
  <div id="overlay">Drag: rotate · Shift+drag: pan · Wheel: zoom</div>
</main>
<script>
const MANIFEST = __MANIFEST__;
const canvas = document.getElementById('canvas');
const gl = canvas.getContext('webgl', { antialias: true, alpha: false });
if (!gl) throw new Error('WebGL is not available in this browser.');

document.getElementById('source').textContent = `${MANIFEST.sourceFile} · ${MANIFEST.backend}/${MANIFEST.quality}`;
document.getElementById('totalObjects').textContent = MANIFEST.objects.toLocaleString();
document.getElementById('totalTriangles').textContent = MANIFEST.triangles.toLocaleString();
document.getElementById('warnings').innerHTML = MANIFEST.warnings.length ? MANIFEST.warnings.map(x => `<div>${escapeHtml(x)}</div>`).join('') : 'No geometry warnings.';

function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
}
function formatBytes(n) {
  const units=['B','KiB','MiB','GiB']; let i=0, v=n;
  while (v>=1024 && i<units.length-1) { v/=1024; i++; }
  return `${v.toFixed(i ? 1 : 0)} ${units[i]}`;
}

const states = Object.fromEntries(MANIFEST.chunks.map(c => [c.key, { loaded:false, loading:false, visible:false, error:'', buffers:null, bytes:0 }]));
const selectedLayers = new Set();
const layerChecks = new Map();
const grouped = {};
for (const chunk of MANIFEST.chunks) {
  grouped[chunk.subsystem] ??= {};
  grouped[chunk.subsystem][chunk.layer] ??= [];
  grouped[chunk.subsystem][chunk.layer].push(chunk);
}

function buildTree() {
  const root = document.getElementById('geometryTree');
  root.innerHTML = '';
  for (const [subsystem, layers] of Object.entries(grouped)) {
    const box=document.createElement('div'); box.className='subsystem';
    const head=document.createElement('label'); head.className='subsystem-head';
    const subCheck=document.createElement('input'); subCheck.type='checkbox'; subCheck.dataset.subsystem=subsystem;
    head.append(subCheck, document.createTextNode(subsystem)); box.append(head);
    const layerBox=document.createElement('div'); layerBox.className='layers';
    for (const [layer, chunks] of Object.entries(layers)) {
      const key=`${subsystem}\x1f${layer}`;
      const label=document.createElement('label'); label.className='layer';
      const check=document.createElement('input'); check.type='checkbox'; check.dataset.layerKey=key;
      layerChecks.set(key,check);
      const text=document.createElement('span'); text.textContent=layer;
      const count=document.createElement('small');
      count.textContent=`${chunks.reduce((s,c)=>s+c.triangles,0).toLocaleString()} △`;
      label.append(check,text,count); layerBox.append(label);
      check.addEventListener('change', () => setLayer(key, check.checked));
    }
    let changeSequence=0;
    subCheck.addEventListener('change', async () => {
      const checks=[...layerBox.querySelectorAll('input[data-layer-key]')];
      const enabled=subCheck.checked;
      const sequence=++changeSequence;
      for (const check of checks) check.checked=enabled;
      for (const check of checks) {
        if (sequence!==changeSequence) return;
        await setLayer(check.dataset.layerKey,enabled);
      }
      if (sequence===changeSequence) syncSubsystemChecks();
    });
    box.append(layerBox); root.append(box);
  }
}

async function setLayer(baseKey, enabled) {
  if (enabled) selectedLayers.add(baseKey); else selectedLayers.delete(baseKey);
  const chunks=MANIFEST.chunks.filter(c => c.baseKey===baseKey);
  if (enabled) {
    await Promise.all(chunks.map(loadChunk));
    const stillEnabled=selectedLayers.has(baseKey);
    for (const chunk of chunks) if (states[chunk.key].loaded) states[chunk.key].visible=stillEnabled;
  } else {
    for (const chunk of chunks) {
      states[chunk.key].visible=false;
      if (document.getElementById('releaseHidden').checked) unloadChunk(chunk);
    }
  }
  syncSubsystemChecks(); updateStats(); requestRender();
}

async function loadChunk(chunk) {
  const state=states[chunk.key];
  if (state.loaded || state.loading) return;
  state.loading=true; state.error=''; updateStatus(`Loading ${chunk.subsystem} / ${chunk.layer}…`, 'loading');
  try {
    const response=await fetch(chunk.url);
    if (!response.ok) throw new Error(`${response.status} ${response.statusText}`);
    const data=await response.arrayBuffer();
    const view=new DataView(data);
    const magic=String.fromCharCode(...new Uint8Array(data,0,8));
    if (magic!=='GDMLCHN1') throw new Error(`Bad chunk magic in ${chunk.url}`);
    const version=view.getUint32(8,true), vertices=view.getUint32(12,true);
    if (version!==1) throw new Error(`Unsupported chunk version ${version}`);
    const posOff=view.getUint32(16,true), colorOff=view.getUint32(20,true), idOff=view.getUint32(24,true), endOff=view.getUint32(28,true);
    if (endOff!==data.byteLength) throw new Error(`Truncated chunk ${chunk.url}`);
    const positions=new Float32Array(data,posOff,vertices*3);
    const colors=new Uint8Array(data,colorOff,vertices*4);
    state.buffers={
      pos:makeBuffer(positions),
      color:makeBuffer(colors),
      count:vertices,
    };
    state.loaded=true; state.visible=selectedLayers.has(chunk.baseKey); state.bytes=data.byteLength;
    updateStatus(`Loaded ${chunk.subsystem} / ${chunk.layer}.`, 'loaded');
  } catch (error) {
    state.error=String(error); updateStatus(`Could not load ${chunk.url}.\n${error}\nUse --serve or another HTTP server; file:// cannot fetch lazy chunks.`, 'error');
  } finally {
    state.loading=false; updateStats(); requestRender();
  }
}

function unloadChunk(chunk) {
  const state=states[chunk.key];
  if (!state.loaded || !state.buffers) return;
  gl.deleteBuffer(state.buffers.pos); gl.deleteBuffer(state.buffers.color);
  state.loaded=false; state.buffers=null; state.bytes=0;
}
function updateStatus(text, cls='') { const el=document.getElementById('status'); el.textContent=text; el.className=`muted ${cls}`; }
function updateStats() {
  const loaded=Object.values(states).filter(s=>s.loaded);
  document.getElementById('loadedChunks').textContent=`${loaded.length}/${MANIFEST.chunks.length}`;
  document.getElementById('loadedBytes').textContent=formatBytes(loaded.reduce((s,x)=>s+x.bytes,0));
}
function syncSubsystemChecks() {
  document.querySelectorAll('input[data-subsystem]').forEach(sub => {
    const keys=Object.keys(grouped[sub.dataset.subsystem]).map(layer=>`${sub.dataset.subsystem}\x1f${layer}`);
    const n=keys.filter(key=>selectedLayers.has(key)).length;
    sub.checked=n===keys.length; sub.indeterminate=n>0 && n<keys.length;
  });
}

function shader(type, source) {
  const value=gl.createShader(type); gl.shaderSource(value,source); gl.compileShader(value);
  if (!gl.getShaderParameter(value,gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(value));
  return value;
}
const vertexShader=shader(gl.VERTEX_SHADER, `
attribute vec3 a_pos; attribute vec4 a_color;
uniform mat4 u_mvp; varying vec4 v_color; varying vec3 v_world;
void main() { v_world=a_pos; v_color=a_color; gl_Position=u_mvp*vec4(a_pos,1.0); }
`);
const fragmentShader=shader(gl.FRAGMENT_SHADER, `
precision mediump float; varying vec4 v_color; varying vec3 v_world;
uniform float u_opacity, u_cut_enabled, u_cut_start, u_cut_width, u_cut_mode;
const float TWO_PI=6.28318530718;
void main() {
  if (u_cut_enabled>0.5) {
    float phi=atan(v_world.y,v_world.x); if (phi<0.0) phi+=TWO_PI;
    float relative=mod(phi-u_cut_start+TWO_PI,TWO_PI);
    bool inside=relative<=u_cut_width;
    if ((u_cut_mode<0.5 && inside) || (u_cut_mode>0.5 && !inside)) discard;
  }
  gl_FragColor=vec4(v_color.rgb,v_color.a*u_opacity);
}
`);
const program=gl.createProgram(); gl.attachShader(program,vertexShader); gl.attachShader(program,fragmentShader); gl.linkProgram(program);
if (!gl.getProgramParameter(program,gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(program));
const locations={
  pos:gl.getAttribLocation(program,'a_pos'), color:gl.getAttribLocation(program,'a_color'),
  mvp:gl.getUniformLocation(program,'u_mvp'), opacity:gl.getUniformLocation(program,'u_opacity'),
  cutEnabled:gl.getUniformLocation(program,'u_cut_enabled'), cutStart:gl.getUniformLocation(program,'u_cut_start'),
  cutWidth:gl.getUniformLocation(program,'u_cut_width'), cutMode:gl.getUniformLocation(program,'u_cut_mode'),
};
function makeBuffer(values) { const b=gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER,b); gl.bufferData(gl.ARRAY_BUFFER,values,gl.STATIC_DRAW); return b; }

const bbox=MANIFEST.bbox;
const center=[(bbox[0]+bbox[3])/2,(bbox[1]+bbox[4])/2,(bbox[2]+bbox[5])/2];
const radius=Math.hypot(bbox[3]-bbox[0],bbox[4]-bbox[1],bbox[5]-bbox[2])||1;
const axisLength=radius*.18;
const axisPositions=new Float32Array([
  0,0,0, axisLength,0,0,
  0,0,0, 0,axisLength,0,
  0,0,0, 0,0,axisLength,
]);
const axisColors=new Uint8Array([
  239,68,68,255, 239,68,68,255,
  34,197,94,255, 34,197,94,255,
  59,130,246,255, 59,130,246,255,
]);
const axisPosBuffer=makeBuffer(axisPositions);
const axisColorBuffer=makeBuffer(axisColors);
let yaw=-0.65,pitch=0.45,dist=radius*1.35,panX=0,panY=0,renderPending=false;
let cameraUp=[0,1,0];
function mat4mul(a,b) { const o=new Float32Array(16); for(let r=0;r<4;r++)for(let c=0;c<4;c++)for(let k=0;k<4;k++)o[c*4+r]+=a[k*4+r]*b[c*4+k]; return o; }
function perspective(fovy,aspect,near,far) { const f=1/Math.tan(fovy/2),nf=1/(near-far); return new Float32Array([f/aspect,0,0,0,0,f,0,0,0,0,(far+near)*nf,-1,0,0,2*far*near*nf,0]); }
function lookAt(eye,target,up) {
  let zx=eye[0]-target[0],zy=eye[1]-target[1],zz=eye[2]-target[2],zl=Math.hypot(zx,zy,zz); zx/=zl;zy/=zl;zz/=zl;
  let xx=up[1]*zz-up[2]*zy,xy=up[2]*zx-up[0]*zz,xz=up[0]*zy-up[1]*zx,xl=Math.hypot(xx,xy,xz); xx/=xl;xy/=xl;xz/=xl;
  const yx=zy*xz-zz*xy,yy=zz*xx-zx*xz,yz=zx*xy-zy*xx;
  return new Float32Array([xx,yx,zx,0,xy,yy,zy,0,xz,yz,zz,0,-(xx*eye[0]+xy*eye[1]+xz*eye[2]),-(yx*eye[0]+yy*eye[1]+yz*eye[2]),-(zx*eye[0]+zy*eye[1]+zz*eye[2]),1]);
}
function resize() { const dpr=window.devicePixelRatio||1,w=Math.floor(canvas.clientWidth*dpr),h=Math.floor(canvas.clientHeight*dpr); if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;gl.viewport(0,0,w,h);} }
function requestRender() { if(renderPending)return; renderPending=true; requestAnimationFrame(render); }
function updateAxisLabels(mvp,enabled) {
  const endpoints=[
    ['axisLabelX',axisLength,0,0],
    ['axisLabelY',0,axisLength,0],
    ['axisLabelZ',0,0,axisLength],
  ];
  for (const [id,x,y,z] of endpoints) {
    const label=document.getElementById(id);
    const tx=mvp[0]*x+mvp[4]*y+mvp[8]*z+mvp[12];
    const ty=mvp[1]*x+mvp[5]*y+mvp[9]*z+mvp[13];
    const tw=mvp[3]*x+mvp[7]*y+mvp[11]*z+mvp[15];
    if (!enabled||tw<=0) { label.style.display='none'; continue; }
    label.style.display='block';
    label.style.left=((tx/tw*.5+.5)*canvas.clientWidth)+'px';
    label.style.top=((-ty/tw*.5+.5)*canvas.clientHeight)+'px';
  }
}
function render() {
  renderPending=false; resize(); gl.clearColor(0.027,0.067,0.122,1); gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);
  gl.enable(gl.DEPTH_TEST); gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA);
  const opacity=+document.getElementById('opacity').value/100; gl.depthMask(opacity>=0.999);
  const cx=center[0]+panX,cy=center[1]+panY,cz=center[2];
  const eye=[cx+dist*Math.cos(pitch)*Math.sin(yaw),cy+dist*Math.sin(pitch),cz+dist*Math.cos(pitch)*Math.cos(yaw)];
  const mvp=mat4mul(perspective(45*Math.PI/180,canvas.width/canvas.height,radius/1000,radius*30),lookAt(eye,[cx,cy,cz],cameraUp));
  gl.useProgram(program); gl.uniformMatrix4fv(locations.mvp,false,mvp); gl.uniform1f(locations.opacity,opacity);
  const enabled=document.getElementById('cutEnabled').checked;
  const start=(+document.getElementById('cutStart').value*Math.PI/180+Math.PI*2)%(Math.PI*2);
  const width=+document.getElementById('cutWidth').value*Math.PI/180;
  gl.uniform1f(locations.cutEnabled,enabled?1:0); gl.uniform1f(locations.cutStart,start); gl.uniform1f(locations.cutWidth,width);
  gl.uniform1f(locations.cutMode,document.getElementById('cutMode').value==='keep'?1:0);
  for (const chunk of MANIFEST.chunks) {
    const state=states[chunk.key]; if(!state.loaded||!state.visible||!state.buffers)continue;
    gl.bindBuffer(gl.ARRAY_BUFFER,state.buffers.pos); gl.enableVertexAttribArray(locations.pos); gl.vertexAttribPointer(locations.pos,3,gl.FLOAT,false,0,0);
    gl.bindBuffer(gl.ARRAY_BUFFER,state.buffers.color); gl.enableVertexAttribArray(locations.color); gl.vertexAttribPointer(locations.color,4,gl.UNSIGNED_BYTE,true,0,0);
    gl.drawArrays(gl.TRIANGLES,0,state.buffers.count);
  }
  const axesEnabled=document.getElementById('coordinateAxes').checked;
  if (axesEnabled) {
    gl.disable(gl.DEPTH_TEST);
    gl.uniform1f(locations.opacity,1);
    gl.uniform1f(locations.cutEnabled,0);
    gl.bindBuffer(gl.ARRAY_BUFFER,axisPosBuffer); gl.vertexAttribPointer(locations.pos,3,gl.FLOAT,false,0,0);
    gl.bindBuffer(gl.ARRAY_BUFFER,axisColorBuffer); gl.vertexAttribPointer(locations.color,4,gl.UNSIGNED_BYTE,true,0,0);
    gl.drawArrays(gl.LINES,0,6);
  }
  updateAxisLabels(mvp,axesEnabled);
}

let drag=null;
canvas.addEventListener('mousedown',e=>drag={x:e.clientX,y:e.clientY,shift:e.shiftKey});
window.addEventListener('mouseup',()=>drag=null);
window.addEventListener('mousemove',e=>{ if(!drag)return; const dx=e.clientX-drag.x,dy=e.clientY-drag.y;drag.x=e.clientX;drag.y=e.clientY; if(drag.shift){panX-=dx*dist/900;panY+=dy*dist/900;}else{yaw+=dx*.006;pitch=Math.max(-1.45,Math.min(1.45,pitch+dy*.006));cameraUp=[0,1,0];} requestRender(); });
canvas.addEventListener('wheel',e=>{e.preventDefault();dist*=Math.exp(e.deltaY*.001);requestRender();},{passive:false});
window.addEventListener('resize',requestRender);
function setCamera(newYaw,newPitch,up) {
  yaw=newYaw; pitch=newPitch; cameraUp=up; dist=radius*1.35; panX=panY=0; requestRender();
}
document.getElementById('resetView').onclick=()=>setCamera(-.65,.45,[0,1,0]);
document.getElementById('viewXYPos').onclick=()=>setCamera(0,0,[0,1,0]);
document.getElementById('viewXYNeg').onclick=()=>setCamera(Math.PI,0,[0,1,0]);
document.getElementById('viewXZPos').onclick=()=>setCamera(0,Math.PI/2,[0,0,1]);
document.getElementById('viewXZNeg').onclick=()=>setCamera(0,-Math.PI/2,[0,0,1]);
document.getElementById('viewYZPos').onclick=()=>setCamera(Math.PI/2,0,[0,0,1]);
document.getElementById('viewYZNeg').onclick=()=>setCamera(-Math.PI/2,0,[0,0,1]);
document.getElementById('coordinateAxes').addEventListener('change',requestRender);
for (const id of ['cutEnabled','cutMode','cutStart','cutWidth','opacity']) document.getElementById(id).addEventListener('input',()=>{
  document.getElementById('cutStartValue').textContent=`${document.getElementById('cutStart').value}°`;
  document.getElementById('cutWidthValue').textContent=`${document.getElementById('cutWidth').value}°`; requestRender();
});
document.getElementById('loadAll').onclick=async()=>{ for(const check of document.querySelectorAll('input[data-layer-key]'))check.checked=true; for(const key of Object.keys(grouped).flatMap(s=>Object.keys(grouped[s]).map(l=>`${s}\x1f${l}`)))await setLayer(key,true); };
document.getElementById('hideAll').onclick=async()=>{ for(const check of document.querySelectorAll('input[data-layer-key]'))check.checked=false; for(const key of [...selectedLayers])await setLayer(key,false); };

function renderMaterials() {
  const root=document.getElementById('materials');
  root.innerHTML=Object.entries(MANIFEST.materials).map(([name,count])=>{const c=MANIFEST.materialColors[name]||[128,128,128,255];return `<div class="material"><span class="swatch" style="background:rgba(${c[0]},${c[1]},${c[2]},${Math.max(c[3]/255,.35)})"></span><span title="${escapeHtml(name)}">${escapeHtml(name)}</span><span>${count.toLocaleString()}</span></div>`;}).join('');
}

buildTree(); renderMaterials(); updateStats(); requestRender();
for (const requested of MANIFEST.defaultLoad) {
  for (const [subsystem,layers] of Object.entries(grouped)) {
    if (requested.toLowerCase()!==subsystem.toLowerCase() && requested.toLowerCase()!=='all') continue;
    for (const layer of Object.keys(layers)) { const key=`${subsystem}\x1f${layer}`; const check=layerChecks.get(key); if(check)check.checked=true; setLayer(key,true); }
  }
}
</script>
</body>
</html>
"""


def write_html(manifest: dict, output: Path) -> None:
    title = html.escape(str(manifest["sourceFile"]))
    manifest_json = json.dumps(manifest, separators=(",", ":")).replace("</", "<\\/")
    output.write_text(
        HTML_TEMPLATE.replace("__TITLE__", title).replace("__MANIFEST__", manifest_json),
        encoding="utf-8",
    )


def read_existing_manifest(output: Path) -> dict:
    """Read the embedded manifest from an HTML file generated by this viewer."""

    document = output.read_text(encoding="utf-8")
    prefix = "const MANIFEST = "
    suffix = ";\nconst canvas = "
    start = document.find(prefix)
    if start < 0:
        raise ValueError("embedded lazy-viewer manifest was not found")
    start += len(prefix)
    end = document.find(suffix, start)
    if end < 0:
        raise ValueError("end of embedded lazy-viewer manifest was not found")
    manifest = json.loads(document[start:end])
    if not isinstance(manifest, dict) or manifest.get("format") != "gdml-lazy-viewer-v1":
        raise ValueError("unsupported lazy-viewer manifest format")
    return manifest


def serve_file(path: Path, port: int) -> None:
    handler = partial(SimpleHTTPRequestHandler, directory=str(path.parent))
    server = ThreadingHTTPServer(("127.0.0.1", port), handler)
    url = f"http://127.0.0.1:{port}/{path.name}"
    print(f"Serving {url}")
    try:
        webbrowser.open(url)
    except Exception:
        pass
    server.serve_forever()


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a lazy-loading, phi-cut-enabled HTML viewer for GDML geometry.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("gdml", type=Path, help="Input GDML file")
    parser.add_argument("--out", type=Path, default=None, help="Output HTML path")
    parser.add_argument("--backend", choices=("auto", "python", "pyg4ometry"), default="auto")
    parser.add_argument("--quality", choices=sorted(base.QUALITY_SEGMENTS), default="fast")
    parser.add_argument("--include-containers", action="store_true")
    parser.add_argument("--include-world", action="store_true")
    parser.add_argument("--max-objects", type=int, default=500_000)
    parser.add_argument("--max-vertices", type=int, default=500_000_000)
    parser.add_argument(
        "--chunk-vertices",
        type=int,
        default=500_000,
        help="Maximum vertices per binary chunk; large layers are split into parts",
    )
    parser.add_argument(
        "--default-load",
        action="append",
        default=[],
        metavar="SUBSYSTEM",
        help="Subsystem to load when the page opens; repeat or use 'all'",
    )
    parser.add_argument(
        "--reuse-assets",
        action="store_true",
        help="Rewrite and optionally serve an existing output HTML without rebuilding its binary chunks",
    )
    parser.add_argument("--serve", action="store_true", help="Serve the viewer on localhost")
    parser.add_argument("--port", type=int, default=8765)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    if not args.gdml.exists():
        print(f"GDML file not found: {args.gdml}", file=sys.stderr)
        return 2
    output = args.out or (Path(tempfile.gettempdir()) / f"{args.gdml.stem}.lazy.html")
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    if args.reuse_assets:
        if not output.exists():
            print(f"Existing lazy-viewer HTML not found: {output}", file=sys.stderr)
            return 4
        try:
            manifest = read_existing_manifest(output)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            print(f"Could not reuse {output}: {exc}", file=sys.stderr)
            return 4
        if manifest.get("sourceFile") != args.gdml.name:
            print(
                f"Existing viewer source is {manifest.get('sourceFile')!r}, "
                f"not {args.gdml.name!r}.",
                file=sys.stderr,
            )
            return 4
        missing_assets = [
            str(chunk.get("url", ""))
            for chunk in manifest.get("chunks", [])
            if not (output.parent / str(chunk.get("url", ""))).is_file()
        ]
        if missing_assets:
            print(
                f"Could not reuse {output}: {len(missing_assets)} binary chunk(s) are missing; "
                f"first missing asset: {missing_assets[0]}",
                file=sys.stderr,
            )
            return 4
        if args.default_load:
            manifest["defaultLoad"] = list(args.default_load)
        write_html(manifest, output)
        print(f"Rewrote {output} using {len(manifest['chunks']):,} existing binary chunks.")
        if args.serve:
            serve_file(output, args.port)
        return 0

    builder_cls = LazyGDMLSceneBuilder
    if args.backend in {"auto", "pyg4ometry"}:
        try:
            import pyg4ometry  # type: ignore  # noqa: F401

            builder_cls = LazyPyg4ometrySceneBuilder
            print("Using pyg4ometry backend for solid tessellation.")
        except ModuleNotFoundError:
            if args.backend == "pyg4ometry":
                print("pyg4ometry is not installed in this environment.", file=sys.stderr)
                return 3
            print("pyg4ometry not installed; using pure-Python backend.")

    builder = builder_cls(
        args.gdml,
        quality=args.quality,
        include_containers=args.include_containers,
        include_world=args.include_world,
        max_objects=args.max_objects,
        max_vertices=args.max_vertices,
        chunk_vertices=args.chunk_vertices,
    )
    builder.build_chunks()
    manifest = write_assets(builder, output, args.default_load)
    write_html(manifest, output)

    asset_bytes = sum(int(chunk["bytes"]) for chunk in manifest["chunks"])
    print(f"Wrote {output}")
    print(f"Assets: {output.parent / (output.stem + '.assets')}")
    print(
        f"Subsystem chunks: {len(manifest['chunks']):,}  "
        f"Objects: {manifest['objects']:,}  Triangles: {manifest['triangles']:,}  "
        f"Binary: {asset_bytes / (1024 * 1024):.1f} MiB"
    )
    for subsystem in sorted({chunk["subsystem"] for chunk in manifest["chunks"]}, key=lambda x: SUBSYSTEM_ORDER.get(x, 99)):
        entries = [chunk for chunk in manifest["chunks"] if chunk["subsystem"] == subsystem]
        print(
            f"  {subsystem}: {len({entry['layer'] for entry in entries})} layers, "
            f"{sum(entry['triangles'] for entry in entries):,} triangles, {len(entries)} chunks"
        )
    if manifest["warnings"]:
        print("Warnings:")
        for warning in manifest["warnings"][:8]:
            print(f"  - {warning}")
        if len(manifest["warnings"]) > 8:
            print(f"  - ... {len(manifest['warnings']) - 8} more")
    print("Open through HTTP (use --serve); direct file:// loading cannot fetch chunks.")
    if args.serve:
        serve_file(output, args.port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
