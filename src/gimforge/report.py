"""A portable, dependency-free HTML report for inspecting real GIM matrices."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Mapping, Sequence

from .io import as_float


def _json_data(rows: Sequence[Mapping[str, object]]) -> str:
    return json.dumps([dict(row) for row in rows], ensure_ascii=False).replace("</", "<\\/")


def write_report(
    output_path: str | Path,
    *,
    matrix_out: Sequence[Mapping[str, object]],
    members: Sequence[Mapping[str, object]],
    gim_summary: Sequence[Mapping[str, object]],
    edges: Sequence[Mapping[str, object]] = (),
    regions: Sequence[Mapping[str, object]] = (),
    title: str = "GIMForge report",
    conditional_p: float = 1.24741348813236e-8,
    metadata: Mapping[str, object] | None = None,
) -> Path:
    """Write a self-contained interactive report from the conditional matrix.

    The browser never re-fits, clusters, or smooths the result. Cell colour is
    the conditional beta, and the retained-edge border is determined only by
    ``conditional_p``.
    """

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    member_sets: dict[str, dict[str, set[str]]] = {}
    for member in members:
        group = member_sets.setdefault(str(member["gim_id"]), {"SNP": set(), "metabolite": set()})
        group.setdefault(str(member["node_type"]), set()).add(str(member["node_id"]))
    region_gims = {str(item["gim_id"]): str(item["region_id"]) for item in gim_summary}
    heatmap_pairs = {
        (region_gims[gim_id], snp, metabolite)
        for gim_id, nodes in member_sets.items()
        for snp in nodes.get("SNP", set())
        for metabolite in nodes.get("metabolite", set())
        if gim_id in region_gims
    }
    # Keep report.html bounded. The full V_R × M_R matrix remains available in
    # matrix_out.tsv.gz; the report only needs each GIM's Cartesian submatrix.
    edge_lookup = {
        (
            str(row.get("region_id")),
            str(row.get("snp_id")),
            str(row.get("metabolite")),
        ): row
        for row in edges
    }
    report_matrix: list[dict[str, object]] = []
    restored_retained_p = 0
    for source_row in matrix_out:
        row = dict(source_row)
        key = (
            str(row.get("region_id")),
            str(row.get("snp_id")),
            str(row.get("metabolite")),
        )
        retained = edge_lookup.get(key)
        if as_float(row.get("p")) is None and retained is not None:
            if as_float(retained.get("p")) is not None:
                restored_retained_p += 1
            for field in ("beta", "se", "p", "n", "conditioned_on_n", "testable"):
                if row.get(field) in (None, "") and retained.get(field) not in (None, ""):
                    row[field] = retained[field]
        report_matrix.append(row)
    visible_matrix = [
        row
        for row in report_matrix
        if (
            str(row.get("region_id")),
            str(row.get("snp_id")),
            str(row.get("metabolite")),
        )
        in heatmap_pairs
    ]
    payload = {
        "matrix": _json_data(visible_matrix),
        "members": _json_data(members),
        "summary": _json_data(gim_summary),
        "regions": _json_data(regions),
    }
    metadata = dict(metadata or {})
    if restored_retained_p:
        metadata["Legacy matrix compatibility"] = (
            f"Restored {restored_retained_p:,} retained P values from edges.tsv.gz"
        )
    metadata_html = "".join(
        f"<tr><th>{html.escape(str(key))}</th><td>{html.escape(str(value))}</td></tr>"
        for key, value in metadata.items()
    )

    template = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>__TITLE__</title>
  <style>
    :root{--ink:#162e38;--muted:#62747b;--line:#d8e2e5;--paper:#f7fafb;--teal:#086a77;--teal-dark:#07515d;--blue:#2b74b8;--red:#cd4b3f;--neutral:#e9eff0;--cell:42px;--label:166px;--header:172px;color-scheme:light;font-family:Inter,ui-sans-serif,system-ui,-apple-system,sans-serif}
    *{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink)}button,select{font:inherit}
    .masthead{min-height:60px;padding:0 22px;display:flex;align-items:center;gap:22px;border-bottom:1px solid var(--line);background:#fff}
    .brand{display:flex;align-items:center;gap:9px;font-size:18px;font-weight:780;white-space:nowrap}.brand em{color:var(--teal);font-style:normal}
    .brand-mark{width:22px;height:22px;display:flex;align-items:end;gap:3px}.brand-mark i{width:5px;border-radius:4px;background:var(--teal)}.brand-mark i:nth-child(1){height:10px}.brand-mark i:nth-child(2){height:19px;background:var(--red)}.brand-mark i:nth-child(3){height:14px}
    .caption{margin:0;color:var(--muted);font-size:13px;flex:1}.dataset-metrics{display:flex;gap:17px;color:var(--muted);font-size:12px;white-space:nowrap}.dataset-metrics b{color:var(--ink)}
    main{padding:14px 18px 22px}.workbench{display:grid;grid-template-columns:minmax(0,1fr) 340px;gap:12px;max-width:1800px;margin:auto}
    .matrix-panel,.detail-panel,.provenance{border:1px solid var(--line);border-radius:7px;background:#fff;box-shadow:0 2px 7px rgba(21,55,63,.04)}
    .matrix-header{min-height:110px;padding:17px 19px;display:flex;justify-content:space-between;gap:20px;border-bottom:1px solid var(--line)}
    .eyebrow{margin:0 0 3px;color:var(--teal);font-size:10px;font-weight:800;letter-spacing:.12em}.gim-heading h1{margin:0;font-size:25px;letter-spacing:-.03em}.locus{margin:6px 0 0;color:var(--muted);font-size:12px}
    .header-right{display:flex;align-items:center;justify-content:flex-end;gap:22px;flex-wrap:wrap}.module-summary{display:flex;gap:17px}.module-summary span{display:grid;color:var(--muted);font-size:10px;text-align:center}.module-summary b{color:var(--ink);font-size:19px}
    .gim-choice{display:grid;gap:3px;color:var(--muted);font-size:9px;font-weight:750}.gim-choice select{width:min(360px,46vw);padding:7px 9px;border:1px solid #bfcfd3;border-radius:4px;background:#fff;color:var(--ink);font-size:11px}
    .toolbar{min-height:43px;padding:7px 14px;display:flex;align-items:center;gap:14px;border-bottom:1px solid var(--line);background:#fbfdfd}
    .tabs{display:flex;border:1px solid #c8d7da;border-radius:4px;padding:2px}.tabs button{border:0;border-radius:3px;background:transparent;color:#587077;padding:5px 8px;font-size:10px;font-weight:750}.tabs button.active{background:var(--teal-dark);color:#fff}
    .zoom{display:flex;align-items:center;gap:5px;color:var(--muted);font-size:10px}.zoom button{width:25px;height:25px;border:1px solid #c3d3d6;background:#fff;border-radius:3px;color:var(--teal-dark);font-weight:800}.zoom output{width:38px;text-align:center}
    .legend{margin-left:auto;display:flex;align-items:center;gap:4px;color:var(--muted);font-size:9px;white-space:nowrap}.legend i{width:23px;height:9px;border-radius:2px}.legend .negative{background:var(--blue)}.legend .neutral{background:var(--neutral)}.legend .positive{background:var(--red)}.legend b{margin-left:7px;color:#425e65}
    .heatmap-scroll{overflow:auto;padding:17px;min-height:420px;max-height:calc(100vh - 285px)}.heatmap-grid{display:grid;grid-template-columns:var(--label) repeat(var(--cols),var(--cell));grid-template-rows:var(--header) repeat(var(--rows),var(--cell));width:max-content;min-width:100%;align-items:stretch}
    .corner{position:sticky;left:0;z-index:4;display:flex;flex-direction:column;justify-content:end;padding:0 10px 8px 0;background:#fff;border-right:1px solid var(--line);border-bottom:1px solid var(--line);color:var(--muted);font-size:9px}.corner b{color:var(--ink);font-size:10px}
    .column-label{position:relative;border-bottom:1px solid var(--line);background:#fff}.column-label span{position:absolute;left:50%;bottom:8px;width:calc(var(--header) - 20px);transform:rotate(-55deg);transform-origin:left center;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;color:#3c5a62;font-size:11px;font-weight:690}
    .row-label{position:sticky;left:0;z-index:3;display:flex;align-items:center;justify-content:flex-end;padding-right:9px;border-right:1px solid var(--line);border-bottom:1px solid #edf1f2;background:#fff;color:#34515a;font-size:11px;font-weight:720;white-space:nowrap}
    .heat-cell{position:relative;border:1px solid rgba(255,255,255,.72);padding:0;cursor:pointer;color:#17313a;font-size:10px;font-weight:760}.heat-cell.retained{box-shadow:inset 0 0 0 2px #132f38}.heat-cell.selected{z-index:2;outline:3px solid #e2a13f;outline-offset:-2px}.heat-cell.hidden{background:#f7f9f9!important;box-shadow:none;cursor:default}.heat-cell span{filter:drop-shadow(0 1px 0 rgba(255,255,255,.45))}.heat-cell:focus-visible{z-index:3;outline:3px solid #e2a13f}
    .matrix-note{margin:0;padding:9px 14px;border-top:1px solid var(--line);color:var(--muted);font-size:10px}.matrix-note b{color:#405c64}
    .detail-panel{min-height:590px;overflow:hidden}.detail-head{display:flex;justify-content:space-between;align-items:start;padding:17px;border-bottom:1px solid var(--line)}.detail-head h2{margin:0;font-size:18px}.status{padding:4px 7px;border-radius:99px;background:#edf2f3;color:var(--muted);font-size:9px;font-weight:780}.status.retained{background:#dceff0;color:var(--teal-dark)}
    .empty-detail{display:grid;place-items:center;min-height:470px;padding:28px;text-align:center;color:var(--muted);font-size:12px}.detail-body{padding:14px}.selected-pair{display:grid;grid-template-columns:minmax(0,1fr) auto minmax(0,1fr);align-items:center;gap:8px;padding:12px;border:1px solid #cbdadd;border-radius:5px;background:#f7fbfb;color:#31545d;font-size:11px;font-weight:760;word-break:break-word}.selected-pair b{color:var(--red)}
    .effect-grid{display:grid;grid-template-columns:1fr 1fr;margin-top:13px;border:1px solid var(--line);border-radius:5px;overflow:hidden}.effect-grid div{padding:12px;border-right:1px solid var(--line);border-bottom:1px solid var(--line)}.effect-grid div:nth-child(2n){border-right:0}.effect-grid div:nth-last-child(-n+2){border-bottom:0}.effect-grid span{display:block;color:var(--muted);font-size:8px;font-weight:760}.effect-grid b{display:block;margin-top:3px;font-size:14px}.positive{color:var(--red)}.negative{color:var(--blue)}
    .annotation{margin:15px 0 0}.annotation h3{margin:0 0 8px;font-size:11px}.annotation dl{margin:0;border-top:1px solid var(--line)}.annotation dl div{display:grid;grid-template-columns:120px minmax(0,1fr);gap:8px;padding:7px 0;border-bottom:1px solid #edf1f2;font-size:10px}.annotation dt{color:var(--muted)}.annotation dd{margin:0;text-align:right;font-weight:690;word-break:break-word}
    .provenance{max-width:1800px;margin:12px auto 0;overflow:hidden}.provenance summary{cursor:pointer;padding:11px 14px;color:#46636a;font-size:11px;font-weight:760}.provenance table{width:100%;border-collapse:collapse;font-size:10px}.provenance th,.provenance td{padding:7px 14px;border-top:1px solid var(--line);text-align:left}.provenance th{width:230px;color:var(--muted);font-weight:650}
    .empty{color:var(--muted);padding:28px}
    @media(max-width:1050px){.workbench{grid-template-columns:1fr}.detail-panel{min-height:0}.empty-detail{min-height:180px}.heatmap-scroll{max-height:none}.caption,.legend b{display:none}}
    @media(max-width:680px){.masthead{padding:0 12px}.dataset-metrics{display:none}main{padding:0}.matrix-panel,.detail-panel,.provenance{border-radius:0;border-left:0;border-right:0}.matrix-header{flex-direction:column}.header-right{justify-content:start}.gim-choice,.gim-choice select{width:100%}.module-summary{display:none}.toolbar{flex-wrap:wrap}.legend{display:none}.heatmap-scroll{padding:10px}.provenance{margin-top:0}}
  </style>
</head>
<body>
  <header class="masthead">
    <div class="brand"><span class="brand-mark"><i></i><i></i><i></i></span><span>GIM<em>Forge</em></span></div>
    <p class="caption">Conditional genetic influence modules in metabolite genetics</p>
    <div class="dataset-metrics"><span><b id="n-gims">0</b> GIMs</span><span><b id="n-snps">0</b> independent SNPs</span><span><b id="n-metabolites">0</b> metabolites</span></div>
  </header>
  <main>
    <section class="workbench">
      <section class="matrix-panel">
        <header class="matrix-header">
          <div class="gim-heading"><p class="eyebrow">CONDITIONALLY INDEPENDENT GIM</p><h1 id="gim-title">No GIM</h1><p class="locus" id="locus">No significant component</p></div>
          <div class="header-right">
            <div class="module-summary"><span><b id="module-snps">0</b>independent SNPs</span><span><b id="module-metabolites">0</b>metabolites</span><span><b id="module-edges">0</b>retained associations</span></div>
            <label class="gim-choice">CHOOSE GIM<select id="gim-select" aria-label="Choose GIM"></select></label>
          </div>
        </header>
        <div class="toolbar">
          <div class="tabs" aria-label="Association view"><button id="all-tests" class="active">All conditional tests</button><button id="retained-only">Retained GIM associations only</button></div>
          <div class="zoom"><span>Scale</span><button id="zoom-out" aria-label="Decrease heatmap scale">−</button><output id="zoom-label">100%</output><button id="zoom-in" aria-label="Increase heatmap scale">+</button></div>
          <div class="legend"><span>β−</span><i class="negative"></i><i class="neutral"></i><i class="positive"></i><span>β+</span><b>dark border = retained association</b></div>
        </div>
        <div class="heatmap-scroll"><div id="heatmap"></div></div>
        <p class="matrix-note"><b>Rows:</b> conditionally independent SNPs. <b>Columns:</b> metabolites belonging to this GIM. Colour encodes conditional β; dark borders identify associations retained at P ≤ __CONDITIONAL_P_TEXT__. Cell labels are −log<sub>10</sub>(P).</p>
      </section>
      <aside class="detail-panel" id="detail"><div class="empty-detail">Select a heatmap cell to inspect its conditional association.</div></aside>
    </section>
    <details class="provenance"><summary>Run provenance and parameter settings</summary><table><tbody>__METADATA__</tbody></table></details>
  </main>
<script>
const matrix=__MATRIX__;
const members=__MEMBERS__;
const summaryRaw=__SUMMARY__;
const regions=__REGIONS__;
const conditionalP=__CONDITIONAL_P__;
let selectedId=null;
let associationView='all';
let zoom=1;
let selectedCell=null;
const regionById=new Map(regions.map(x=>[String(x.region_id),x]));
const memberByGim=new Map();
for(const row of members){const id=String(row.gim_id);if(!memberByGim.has(id))memberByGim.set(id,{SNP:[],metabolite:[],orders:{}});const group=memberByGim.get(id);const type=String(row.node_type);if(!group[type])group[type]=[];if(!group[type].includes(String(row.node_id)))group[type].push(String(row.node_id));if(type==='SNP')group.orders[String(row.node_id)]=numeric(row.marker_order);}
const rowsByRegion=new Map();
for(const row of matrix){const id=String(row.region_id);if(!rowsByRegion.has(id))rowsByRegion.set(id,[]);rowsByRegion.get(id).push(row);}
function numeric(value){if(value===null||value===undefined||String(value).trim()==='')return null;const n=Number(value);return Number.isFinite(n)?n:null}
function isRetained(row){const p=numeric(row?.p);return p!==null&&p>=0&&p<=conditionalP}
function esc(value){const node=document.createElement('div');node.textContent=String(value??'');return node.innerHTML}
function truthy(value){return value===true||['true','1','yes'].includes(String(value??'').toLowerCase())}
function pFormat(value,testable=null){const p=numeric(value);if(p===0)return'< machine precision';if(p!==null&&p>0)return p.toExponential(2);return truthy(testable)?'result unavailable':'not tested'}
function effectFormat(value){const n=numeric(value);return n===null?'—':`${n>=0?'+':''}${n.toPrecision(4)}`}
function score(value){const p=numeric(value);if(p===null||p<0)return null;return p===0?323.31:-Math.log10(p)}
function regionLabel(id){const r=regionById.get(String(id));if(!r)return String(id);const chr=r.chromosome?`chr${r.chromosome}`:'';const start=numeric(r.start),end=numeric(r.end);return start!==null&&end!==null?`${chr}:${start.toLocaleString()}–${end.toLocaleString()}`:String(id)}
function rowsForGim(item){const group=memberByGim.get(String(item.gim_id))||{SNP:[],metabolite:[]};const snps=new Set(group.SNP||[]),mets=new Set(group.metabolite||[]);return (rowsByRegion.get(String(item.region_id))||[]).filter(x=>snps.has(String(x.snp_id))&&mets.has(String(x.metabolite)))}
const summary=summaryRaw.map(item=>{const rows=rowsForGim(item);const pValues=rows.map(x=>numeric(x.p)).filter(x=>x!==null&&x>=0);return{...item,n_snps:numeric(item.n_snps)??(memberByGim.get(String(item.gim_id))?.SNP.length||0),n_metabolites:numeric(item.n_metabolites)??(memberByGim.get(String(item.gim_id))?.metabolite.length||0),minP:pValues.length?Math.min(...pValues):null,nEdges:rows.filter(isRetained).length}}).sort((a,b)=>b.n_metabolites-a.n_metabolites||b.n_snps-a.n_snps||(a.minP??1)-(b.minP??1)||String(a.gim_id).localeCompare(String(b.gim_id)));
document.querySelector('#n-gims').textContent=summary.length.toLocaleString();
document.querySelector('#n-snps').textContent=new Set(members.filter(x=>x.node_type==='SNP').map(x=>x.node_id)).size.toLocaleString();
document.querySelector('#n-metabolites').textContent=new Set(members.filter(x=>x.node_type==='metabolite').map(x=>x.node_id)).size.toLocaleString();
const select=document.querySelector('#gim-select');
for(const item of summary){const option=document.createElement('option');option.value=String(item.gim_id);option.textContent=`${item.gim_id} · ${item.n_snps} SNP × ${item.n_metabolites} metabolite${item.n_metabolites===1?'':'s'}`;select.appendChild(option)}
function blend(endpoint,t){const base=[233,239,240];const rgb=base.map((v,i)=>Math.round(v+(endpoint[i]-v)*Math.pow(Math.max(0,Math.min(1,t)),.72)));return`rgb(${rgb.join(',')})`}
function heatColor(beta,maxAbs){const b=numeric(beta);if(b===null)return'#f4f6f6';const t=Math.abs(b)/Math.max(maxAbs,.000001);return blend(b>=0?[205,75,63]:[43,116,184],t)}
function openDetail(row,item){
  selectedCell=row;
  const retained=isRetained(row),beta=numeric(row.beta),pScore=score(row.p);
  document.querySelector('#detail').innerHTML=`<header class="detail-head"><div><p class="eyebrow">GIM ASSOCIATION</p><h2>${esc(item.gim_id)}</h2></div><span class="status ${retained?'retained':''}">${retained?'Retained association':'Conditional test'}</span></header><div class="detail-body"><div class="selected-pair"><span>${esc(row.snp_id)}</span><b>×</b><span>${esc(row.metabolite)}</span></div><div class="effect-grid"><div><span>CONDITIONAL β</span><b class="${beta!==null&&beta>=0?'positive':'negative'}">${effectFormat(row.beta)}</b></div><div><span>P VALUE</span><b>${pFormat(row.p,row.testable)}</b></div><div><span>−LOG10 P</span><b>${pScore===null?'—':pScore.toFixed(2)}</b></div><div><span>CONDITIONED ON</span><b>${esc(row.conditioned_on_n??'—')}</b></div></div><section class="annotation"><h3>Conditional test details</h3><dl><div><dt>Region</dt><dd>${esc(regionLabel(item.region_id))}</dd></div><div><dt>Marker order</dt><dd>${esc(row.marker_order??'—')}</dd></div><div><dt>Standard error</dt><dd>${effectFormat(row.se)}</dd></div><div><dt>Sample size</dt><dd>${esc(row.n??'—')}</dd></div><div><dt>Testable</dt><dd>${esc(row.testable??'—')}</dd></div><div><dt>GIM edge threshold</dt><dd>P ≤ ${conditionalP.toExponential(3)}</dd></div></dl></section></div>`;
  renderHeatmap(item,false);
}
function renderHeatmap(item,resetDetail=true){
  const group=memberByGim.get(String(item.gim_id))||{SNP:[],metabolite:[],orders:{}};
  const snps=[...(group.SNP||[])].sort((a,b)=>(group.orders[a]??1e9)-(group.orders[b]??1e9)||a.localeCompare(b));
  const metabolites=[...(group.metabolite||[])];
  const rows=rowsForGim(item),lookup=new Map(rows.map(x=>[`${x.snp_id}\u0000${x.metabolite}`,x]));
  const betas=rows.map(x=>numeric(x.beta)).filter(x=>x!==null),maxAbs=Math.max(...betas.map(Math.abs),.01);
  document.querySelector('#gim-title').textContent=item.gim_id;
  document.querySelector('#locus').innerHTML=`${esc(regionLabel(item.region_id))} · P<sub>min</sub> ${pFormat(item.minP)}`;
  document.querySelector('#module-snps').textContent=snps.length;
  document.querySelector('#module-metabolites').textContent=metabolites.length;
  document.querySelector('#module-edges').textContent=item.nEdges;
  const cellSize=Math.max(22,Math.round(42*zoom)),labelWidth=Math.round(166*Math.max(zoom,.88)),longest=Math.max(12,...metabolites.map(x=>String(x).length)),headerHeight=Math.max(150,Math.min(330,Math.ceil(longest*7.3+26)));
  let out=`<div class="heatmap-grid" style="--cols:${metabolites.length};--rows:${snps.length};--cell:${cellSize}px;--label:${labelWidth}px;--header:${headerHeight}px"><div class="corner"><span>Metabolites →</span><b>Independent SNPs ↓</b></div>`;
  out+=metabolites.map(m=>`<div class="column-label" title="${esc(m)}"><span>${esc(m)}</span></div>`).join('');
  for(const snp of snps){out+=`<div class="row-label">${esc(snp)}${group.orders[snp]!==null&&group.orders[snp]!==undefined?` · ${group.orders[snp]}`:''}</div>`;for(const metabolite of metabolites){const row=lookup.get(`${snp}\u0000${metabolite}`),retained=isRetained(row),hidden=associationView==='retained'&&!retained,active=selectedCell&&String(selectedCell.snp_id)===snp&&String(selectedCell.metabolite)===metabolite;const pScore=row?score(row.p):null;const text=cellSize>=34&&pScore!==null?pScore.toFixed(pScore>=100?0:1):'';const title=row?`${snp} × ${metabolite}\nβ ${effectFormat(row.beta)} · P ${pFormat(row.p,row.testable)}${retained?' · retained association':''}`:'Not tested';out+=`<button class="heat-cell ${retained?'retained':''} ${hidden?'hidden':''} ${active?'selected':''}" ${!row||hidden?'disabled':''} data-snp="${esc(snp)}" data-metabolite="${esc(metabolite)}" style="background:${hidden?'#f7f9f9':heatColor(row?.beta,maxAbs)}" title="${esc(title)}"><span>${hidden?'':text}</span></button>`}}
  out+='</div>';
  document.querySelector('#heatmap').innerHTML=out;
  for(const button of document.querySelectorAll('.heat-cell:not(:disabled)')){button.addEventListener('click',()=>{const row=lookup.get(`${button.dataset.snp}\u0000${button.dataset.metabolite}`);if(row)openDetail(row,item)})}
  if(resetDetail){selectedCell=null;document.querySelector('#detail').innerHTML='<div class="empty-detail">Select a heatmap cell to inspect its conditional association.</div>'}
}
function renderSelected(resetDetail=true){const item=summary.find(x=>String(x.gim_id)===selectedId);if(item)renderHeatmap(item,resetDetail)}
select.addEventListener('change',()=>{selectedId=select.value;renderSelected()});
document.querySelector('#all-tests').addEventListener('click',()=>{associationView='all';document.querySelector('#all-tests').classList.add('active');document.querySelector('#retained-only').classList.remove('active');renderSelected(false)});
document.querySelector('#retained-only').addEventListener('click',()=>{associationView='retained';document.querySelector('#retained-only').classList.add('active');document.querySelector('#all-tests').classList.remove('active');renderSelected(false)});
function setZoom(value){zoom=Math.max(.5,Math.min(1.4,value));document.querySelector('#zoom-label').textContent=`${Math.round(zoom*100)}%`;renderSelected(false)}
document.querySelector('#zoom-out').addEventListener('click',()=>setZoom(zoom-.1));document.querySelector('#zoom-in').addEventListener('click',()=>setZoom(zoom+.1));
if(summary.length){selectedId=String(summary[0].gim_id);select.value=selectedId;renderSelected()}else{document.querySelector('#heatmap').innerHTML='<p class="empty">No significant GIM component was found.</p>'}
</script>
</body>
</html>"""

    replacements = {
        "__TITLE__": html.escape(title),
        "__CONDITIONAL_P__": repr(float(conditional_p)),
        "__CONDITIONAL_P_TEXT__": html.escape(f"{conditional_p:.4g}"),
        "__METADATA__": metadata_html or "<tr><td>No run metadata supplied.</td></tr>",
        "__MATRIX__": payload["matrix"],
        "__MEMBERS__": payload["members"],
        "__SUMMARY__": payload["summary"],
        "__REGIONS__": payload["regions"],
    }
    report = template
    for marker, value in replacements.items():
        report = report.replace(marker, value)
    output_path.write_text(report, encoding="utf-8")
    return output_path
