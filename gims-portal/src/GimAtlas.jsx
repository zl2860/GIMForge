import { useEffect, useMemo, useRef, useState } from 'react'
import './gim-atlas.css'

const number = new Intl.NumberFormat('en-US')
const phewasCache = new Map()

function formatP(value) {
  if (!Number.isFinite(value)) return '—'
  return value < 0.001 ? value.toExponential(2).replace('e-', ' × 10⁻') : value.toFixed(4)
}

function finiteNumber(value) {
  if (value === null || value === undefined || value === '') return null
  const parsed = Number(value)
  return Number.isFinite(parsed) ? parsed : null
}

function scoreLabel(p) {
  if (!Number.isFinite(p) || p <= 0) return '—'
  return String(Math.max(0, Math.ceil(-Math.log10(p))))
}

function formatGimId(id) {
  const match = id?.match(/(?:^|_)region_(\d+)_GIM_(\d+)$/)
  return match ? `GIM ${Number(match[1])}.${Number(match[2])}` : id || 'GIM'
}

function locusNumber(regionId) {
  const match = regionId?.match(/(?:^|_)region_(\d+)$/)
  return match ? String(Number(match[1])) : regionId || '—'
}

function gimNumber(gimId) {
  return gimId?.match(/GIM_(\d+)$/)?.[1]?.replace(/^0+/, '') || '1'
}

function chromosomeOrder(left, right) {
  const leftNumber = Number(left)
  const rightNumber = Number(right)
  if (Number.isFinite(leftNumber) && Number.isFinite(rightNumber)) return leftNumber - rightNumber
  if (Number.isFinite(leftNumber)) return -1
  if (Number.isFinite(rightNumber)) return 1
  return String(left).localeCompare(String(right))
}

function coordinateLabel(region) {
  if (region?.chromosomeGrch38 && region?.startGrch38) {
    const end = region.endGrch38 || region.startGrch38
    return `chr${region.chromosomeGrch38}:${number.format(region.startGrch38)}–${number.format(end)} (GRCh38)`
  }
  if (region?.chromosome && region?.startGrch37) {
    const end = region.endGrch37 || region.startGrch37
    return `chr${region.chromosome}:${number.format(region.startGrch37)}–${number.format(end)} (GRCh37)`
  }
  return 'Coordinate unresolved'
}

function variantCoordinate(snp) {
  if (snp?.grch38?.position) return `chr${snp.grch38.chromosome}:${number.format(snp.grch38.position)} (GRCh38)`
  if (snp?.grch37?.grch37Position) return `chr${snp.grch37.chromosome}:${number.format(snp.grch37.grch37Position)} (GRCh37)`
  return 'Coordinate unresolved'
}

function candidateGenes(entity, snpById) {
  return [...new Set(entity.snps.flatMap((snpId) => snpById.get(snpId)?.geneSymbols || []))]
}

function metaboliteClass(value) {
  if (/^(TAG|DAG|FFA|CE|Cer|GM3|LBPA|LPI|LPC|LPS|Lyso|PA|PC|PE|PI|PS|SM|SPH|Sph)/i.test(value || '')) return { code: 'L', label: 'Lipid', color: '#f6f1ac' }
  if (/(amino|alanine|arginine|aspart|glutam|glycine|histidine|leucine|lysine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|valine)/i.test(value || '')) return { code: 'A', label: 'Amino acid', color: '#e8bfd0' }
  return { code: 'M', label: 'Metabolite', color: '#c9e5c3' }
}

function heatColor(cell, signLimit) {
  if (!cell || !Number.isFinite(cell.beta)) return '#fff'
  const fraction = Math.min(Math.abs(cell.beta) / Math.max(signLimit, 0.01), 1)
  const strength = 0.16 + Math.sqrt(fraction) * 0.84
  const neutral = [242, 242, 245]
  const target = cell.beta < 0 ? [104, 83, 218] : [239, 55, 49]
  const rgb = target.map((channel, index) => Math.round(neutral[index] + (channel - neutral[index]) * strength))
  return `rgb(${rgb.join(', ')})`
}

function representativeCell(entity, metabolite) {
  const direct = entity.heatmap
    .filter((cell) => cell.direct && cell.metabolite === metabolite)
    .sort((left, right) => (left.p ?? 1) - (right.p ?? 1))
  if (direct.length) return direct[0]
  return entity.heatmap
    .filter((cell) => cell.testable && cell.metabolite === metabolite)
    .sort((left, right) => (left.p ?? 1) - (right.p ?? 1))[0]
}

function buildSelection(entity, metabolite, region, snpById) {
  const cell = representativeCell(entity, metabolite)
  return {
    gimId: entity.gimId,
    metabolite,
    entity,
    region,
    cell,
    genes: candidateGenes(entity, snpById),
  }
}

function GlobalGimMatrix({ entities, regions, snpById, selection, onSelect }) {
  const [query, setQuery] = useState('')
  const [metaboliteQuery, setMetaboliteQuery] = useState('')
  const [chromosome, setChromosome] = useState('all')
  const [direction, setDirection] = useState('all')
  const [zoom, setZoom] = useState(1)
  const [viewport, setViewport] = useState({ x: 0, y: 0, width: 20, height: 20 })
  const scrollRef = useRef(null)
  const regionById = useMemo(() => new Map(regions.map((region) => [region.regionId, region])), [regions])

  const allColumns = useMemo(() => entities.map((entity) => {
    const region = regionById.get(entity.regionId) || {}
    const genes = candidateGenes(entity, snpById)
    return {
      entity,
      region,
      genes,
      geneLabel: genes.join('/') || 'Unassigned',
      chromosome: String(region.chromosomeGrch38 || region.chromosome || 'Unresolved'),
      cytoband: region.cytoband || 'Unresolved',
      locus: locusNumber(entity.regionId),
      gim: gimNumber(entity.gimId),
    }
  }).sort((left, right) => Number(left.locus) - Number(right.locus) || Number(left.gim) - Number(right.gim)), [entities, regionById, snpById])

  const chromosomes = useMemo(() => [...new Set(allColumns.map((column) => column.chromosome))].sort(chromosomeOrder), [allColumns])
  const visibleColumns = useMemo(() => {
    const needle = query.trim().toLowerCase()
    const exactGim = needle.match(/^gim\s+(\d+)\.(\d+)$/)
    const exactRegion = needle.match(/^(?:region|locus)\s+(\d+)$/)
    return allColumns.filter((column) => {
      if (chromosome !== 'all' && column.chromosome !== chromosome) return false
      if (!needle) return true
      if (exactGim) return Number(column.locus) === Number(exactGim[1]) && Number(column.gim) === Number(exactGim[2])
      if (exactRegion) return Number(column.locus) === Number(exactRegion[1])
      return [column.geneLabel, column.cytoband, column.entity.gimId, column.entity.regionId, formatGimId(column.entity.gimId), `region ${column.locus}`, `locus ${column.locus}`].join(' ').toLowerCase().includes(needle)
    })
  }, [allColumns, chromosome, query])

  const allMetabolites = useMemo(() => [...new Set(entities.flatMap((entity) => entity.metabolites))].sort((left, right) => {
    const leftClass = metaboliteClass(left).code
    const rightClass = metaboliteClass(right).code
    return leftClass.localeCompare(rightClass) || left.localeCompare(right)
  }), [entities])
  const visibleMetabolites = useMemo(() => {
    const needle = metaboliteQuery.trim().toLowerCase()
    return needle ? allMetabolites.filter((metabolite) => metabolite.toLowerCase().includes(needle)) : allMetabolites
  }, [allMetabolites, metaboliteQuery])

  const globalCells = useMemo(() => {
    const rows = []
    visibleColumns.forEach((column, columnIndex) => {
      column.entity.metabolites.forEach((metabolite) => {
        const rowIndex = visibleMetabolites.indexOf(metabolite)
        if (rowIndex < 0) return
        const cell = representativeCell(column.entity, metabolite)
        if (!cell || (direction === 'positive' && cell.beta < 0) || (direction === 'negative' && cell.beta >= 0)) return
        rows.push({ ...cell, column, columnIndex, rowIndex })
      })
    })
    return rows
  }, [direction, visibleColumns, visibleMetabolites])

  const positiveLimit = Math.max(...globalCells.filter((cell) => cell.beta >= 0).map((cell) => cell.beta), 0.01)
  const negativeLimit = Math.max(...globalCells.filter((cell) => cell.beta < 0).map((cell) => Math.abs(cell.beta)), 0.01)
  const cellSize = Math.max(13, Math.round(21 * zoom))
  const rowLabelWidth = 190
  const headerHeight = 190
  const matrixWidth = Math.max(1, visibleColumns.length * cellSize)
  const matrixHeight = Math.max(1, visibleMetabolites.length * cellSize)

  const regionRuns = useMemo(() => {
    const runs = []
    visibleColumns.forEach((column, index) => {
      const key = column.entity.regionId
      const previous = runs[runs.length - 1]
      if (previous?.key === key) previous.count += 1
      else runs.push({ key, start: index, count: 1, locus: column.locus })
    })
    return runs
  }, [visibleColumns])

  const updateViewport = () => {
    const scroller = scrollRef.current
    if (!scroller) return
    setViewport({
      x: Math.max(scroller.scrollLeft - rowLabelWidth, 0) / cellSize,
      y: Math.max(scroller.scrollTop - headerHeight, 0) / cellSize,
      width: scroller.clientWidth / cellSize,
      height: scroller.clientHeight / cellSize,
    })
  }

  useEffect(updateViewport, [cellSize, headerHeight, matrixHeight, matrixWidth, rowLabelWidth])

  const resetFilters = () => {
    setQuery('')
    setMetaboliteQuery('')
    setChromosome('all')
    setDirection('all')
  }

  return <section className="atlas-overview" aria-label="Global GIM locus by metabolite heatmap">
    <header className="atlas-titlebar">
      <div><p className="atlas-eyebrow">GLOBAL GIM–METABOLITE ATLAS</p><h1>All regions and associated metabolic traits</h1><p>Rows are metabolites; columns are every region/GIM in numerical order, labelled by candidate gene.</p></div>
      <div className="atlas-summary"><span><b>{number.format(visibleColumns.length)}</b> GIM columns</span><span><b>{number.format(visibleMetabolites.length)}</b> metabolites</span><span><b>{number.format(globalCells.length)}</b> displayed associations</span></div>
    </header>
    <div className="atlas-controls">
      <label>Gene or region/GIM<input value={query} onChange={(event) => setQuery(event.target.value)} placeholder="e.g. EGLN1, region 27, GIM 27.2" /></label>
      <label>Metabolite<input value={metaboliteQuery} onChange={(event) => setMetaboliteQuery(event.target.value)} placeholder="e.g. TAG54:5" /></label>
      <label>Chromosome<select value={chromosome} onChange={(event) => setChromosome(event.target.value)}><option value="all">All chromosomes</option>{chromosomes.map((item) => <option key={item} value={item}>{item === 'Unresolved' ? item : `chr${item}`}</option>)}</select></label>
      <label>Effect<select value={direction} onChange={(event) => setDirection(event.target.value)}><option value="all">Both directions</option><option value="positive">β positive</option><option value="negative">β negative</option></select></label>
      <div className="atlas-zoom" aria-label="Global heatmap scale"><span>Scale</span><button onClick={() => setZoom((value) => Math.max(0.65, Number((value - 0.1).toFixed(2))))}>−</button><output>{Math.round(zoom * 100)}%</output><button onClick={() => setZoom((value) => Math.min(1.25, Number((value + 0.1).toFixed(2))))}>+</button></div>
      <button className="atlas-reset" onClick={resetFilters}>Reset</button>
      <div className="atlas-color-key"><span>β−</span><i className="negative"></i><i className="neutral"></i><i className="positive"></i><span>β+</span></div>
    </div>
    <div className="atlas-matrix-area">
      <div className="atlas-matrix-scroller" ref={scrollRef} onScroll={updateViewport}>
        <div className="atlas-matrix-layout" style={{ gridTemplateColumns: `${rowLabelWidth}px ${matrixWidth}px`, gridTemplateRows: `${headerHeight}px ${matrixHeight}px` }}>
          <div className="atlas-corner"><b>Metabolites ↓</b><span>Region/GIM order · candidate gene →</span></div>
          <div className="atlas-x-axis" style={{ width: `${matrixWidth}px` }}>
            <div className="atlas-band-strip">{regionRuns.map((run, index) => <span key={`${run.key}-${run.start}`} className={index % 2 ? 'alternate' : ''} style={{ left: `${run.start * cellSize}px`, width: `${run.count * cellSize}px` }} title={`Region ${run.locus} · ${run.count} GIMs`}><b>R{run.locus}</b></span>)}</div>
            {visibleColumns.map((column, index) => <div className="atlas-column-label" key={column.entity.gimId} style={{ left: `${index * cellSize}px`, width: `${cellSize}px`, height: `${headerHeight - 22}px` }} title={`Region ${column.locus} · ${column.geneLabel} · ${formatGimId(column.entity.gimId)} · ${column.cytoband}`}><span><b>{column.geneLabel}</b><small>R{column.locus}.{column.gim}</small></span></div>)}
          </div>
          <div className="atlas-y-axis" style={{ height: `${matrixHeight}px` }}>{visibleMetabolites.map((metabolite, index) => { const category = metaboliteClass(metabolite); return <div className="atlas-row-label" key={metabolite} style={{ top: `${index * cellSize}px`, height: `${cellSize}px` }} title={`${category.label} · ${metabolite}`}><span>{metabolite}</span><i style={{ background: category.color }}>{category.code}</i></div> })}</div>
          <div className="atlas-cells" style={{ width: `${matrixWidth}px`, height: `${matrixHeight}px`, '--atlas-cell': `${cellSize}px` }}>
            {globalCells.map((cell) => {
              const isSelected = selection?.gimId === cell.column.entity.gimId && selection?.metabolite === cell.metabolite
              const signLimit = cell.beta < 0 ? negativeLimit : positiveLimit
              return <button type="button" className={`atlas-cell ${isSelected ? 'selected' : ''}`} key={`${cell.column.entity.gimId}-${cell.metabolite}`} style={{ left: `${cell.columnIndex * cellSize}px`, top: `${cell.rowIndex * cellSize}px`, width: `${cellSize}px`, height: `${cellSize}px`, backgroundColor: heatColor(cell, signLimit) }} onClick={() => onSelect(buildSelection(cell.column.entity, cell.metabolite, cell.column.region, snpById))} title={`${cell.column.cytoband} · ${cell.column.geneLabel} · ${formatGimId(cell.column.entity.gimId)}\n${cell.metabolite}\nβ ${cell.beta >= 0 ? '+' : ''}${cell.beta.toFixed(3)} · P ${formatP(cell.p)}`}>{cellSize >= 18 && <span>{scoreLabel(cell.p)}</span>}</button>
            })}
          </div>
        </div>
      </div>
      <aside className="atlas-minimap" aria-label="Global heatmap minimap">
        <header><b>Overview</b><span>{visibleColumns.length} × {visibleMetabolites.length}</span></header>
        <svg viewBox={`0 0 ${Math.max(visibleColumns.length, 1)} ${Math.max(visibleMetabolites.length, 1)}`} preserveAspectRatio="none">
          <rect className="minimap-background" x="0" y="0" width={Math.max(visibleColumns.length, 1)} height={Math.max(visibleMetabolites.length, 1)}></rect>
          {globalCells.map((cell) => <rect className={cell.beta < 0 ? 'negative' : 'positive'} key={`${cell.column.entity.gimId}-${cell.metabolite}`} x={cell.columnIndex} y={cell.rowIndex} width="1" height="1"></rect>)}
          <rect className="minimap-viewport" x={viewport.x} y={viewport.y} width={Math.min(viewport.width, visibleColumns.length)} height={Math.min(viewport.height, visibleMetabolites.length)}></rect>
        </svg>
        <p>Drag the main matrix scrollbars to navigate the complete sparse atlas.</p>
      </aside>
    </div>
    <p className="atlas-instruction"><b>Click a coloured square</b> to open its SNP×metabolite GIM and locus annotations. Cell text is rounded −log<sub>10</sub>(P); colour is the conditional β direction and relative magnitude.</p>
  </section>
}

function SmallGimHeatmap({ entity, activeCell, snpById, onSelect }) {
  const directKeys = new Set(entity.heatmap.filter((cell) => cell.direct).map((cell) => `${cell.snpId}|${cell.metabolite}`))
  const cellMap = new Map(entity.heatmap.map((cell) => [`${cell.snpId}|${cell.metabolite}`, cell]))
  const positiveLimit = Math.max(...entity.heatmap.filter((cell) => cell.direct && cell.beta >= 0).map((cell) => cell.beta), 0.01)
  const negativeLimit = Math.max(...entity.heatmap.filter((cell) => cell.direct && cell.beta < 0).map((cell) => Math.abs(cell.beta)), 0.01)
  const cellSize = 21
  return <div className="locus-small-matrix-scroll"><div className="locus-small-matrix" style={{ gridTemplateColumns: `250px repeat(${entity.metabolites.length}, ${cellSize}px)`, gridTemplateRows: `150px repeat(${entity.snps.length}, ${cellSize}px)` }}>
    <div className="locus-small-corner"></div>
    {entity.metabolites.map((metabolite, index) => {
      const category = metaboliteClass(metabolite)
      return <div
        className="locus-small-column"
        key={metabolite}
        style={{ zIndex: entity.metabolites.length - index + 3 }}
        title={`${category.label} · ${metabolite}`}
      ><span><i style={{ background: category.color }}>{category.code}</i>{metabolite}</span></div>
    })}
    {entity.snps.flatMap((snpId) => {
      const snp = snpById.get(snpId)
      return [<div className="locus-small-row" key={`${snpId}-label`}><span>{snpId} <small>{snp?.geneSymbols?.[0] || ''}</small></span><i>{gimNumber(entity.gimId)}</i></div>, ...entity.metabolites.map((metabolite) => {
        const cell = cellMap.get(`${snpId}|${metabolite}`)
        const visible = directKeys.has(`${snpId}|${metabolite}`)
        const selected = activeCell?.snpId === snpId && activeCell?.metabolite === metabolite
        const signLimit = cell?.beta < 0 ? negativeLimit : positiveLimit
        return <button type="button" key={`${snpId}-${metabolite}`} disabled={!visible} className={`locus-small-cell ${visible ? '' : 'blank'} ${selected ? 'selected' : ''}`} style={{ backgroundColor: visible ? heatColor(cell, signLimit) : undefined }} onClick={() => visible && onSelect(cell)} title={visible ? `${snpId} × ${metabolite}\nβ ${cell.beta >= 0 ? '+' : ''}${cell.beta.toFixed(3)} · P ${formatP(cell.p)}` : 'No retained GIM association'}>{visible && <span>{scoreLabel(cell.p)}</span>}</button>
      })]
    })}
  </div></div>
}

function SummaryTab({ entity, activeCell, snpById, onSelect }) {
  const directCells = entity.heatmap.filter((cell) => cell.direct).sort((left, right) => (left.p ?? 1) - (right.p ?? 1))
  const snp = snpById.get(activeCell?.snpId)
  return <div className="locus-tab-body">
    {activeCell && <><div className="locus-selected-pair"><span>{activeCell.snpId}</span><b>×</b><span>{activeCell.metabolite}</span></div><div className="locus-effect-grid"><div><span>CONDITIONAL β</span><b className={activeCell.beta >= 0 ? 'positive' : 'negative'}>{activeCell.beta >= 0 ? '+' : ''}{activeCell.beta.toFixed(3)}</b></div><div><span>P VALUE</span><b>{formatP(activeCell.p)}</b></div><div><span>−LOG10(P)</span><b>{(-Math.log10(activeCell.p)).toFixed(2)}</b></div><div><span>CONDITIONED ON</span><b>{activeCell.conditionedOn ?? '—'} variants</b></div><div><span>GENE</span><b>{snp?.geneSymbols?.join(', ') || 'Unassigned'}</b></div><div><span>COORDINATE</span><b>{variantCoordinate(snp)}</b></div></div></>}
    <div className="locus-table-wrap"><table className="locus-table"><thead><tr><th>GIM</th><th>SNP</th><th>Candidate gene</th><th>MV −log<sub>10</sub>(P)</th><th>MV effect</th><th>Conditioned on</th><th>Metabolite</th></tr></thead><tbody>{directCells.map((cell) => { const itemSnp = snpById.get(cell.snpId); const selected = activeCell?.snpId === cell.snpId && activeCell?.metabolite === cell.metabolite; return <tr className={selected ? 'active' : ''} key={`${cell.snpId}-${cell.metabolite}`}><td>{gimNumber(entity.gimId)}</td><td><button onClick={() => onSelect(cell)}>{cell.snpId}</button></td><td>{itemSnp?.geneSymbols?.join(', ') || '—'}</td><td>{Number.isFinite(cell.p) ? (-Math.log10(cell.p)).toFixed(2) : '—'}</td><td className={cell.beta >= 0 ? 'positive' : 'negative'}>{cell.beta >= 0 ? '+' : ''}{cell.beta.toFixed(3)}</td><td>{cell.conditionedOn ?? '—'}</td><td><button onClick={() => onSelect(cell)}>{cell.metabolite}</button></td></tr> })}</tbody></table></div>
  </div>
}

function GenesTab({ entity, snpById, region }) {
  return <div className="locus-tab-body locus-gene-list">{entity.snps.map((snpId) => {
    const snp = snpById.get(snpId)
    return <article key={snpId}><header><div><b>{snpId}</b><span>{snp?.geneSymbols?.join(', ') || 'No HGNC symbol assigned'}</span></div><small>{variantCoordinate(snp)}</small></header><dl><div><dt>Cytoband</dt><dd>{region?.cytoband || 'Unresolved'}</dd></div><div><dt>Consequence</dt><dd>{snp?.mostSevereConsequence?.replaceAll('_', ' ') || 'Not resolved'}</dd></div><div><dt>Ensembl gene</dt><dd>{snp?.geneIds?.join(', ') || '—'}</dd></div><div><dt>Transcript</dt><dd>{snp?.transcriptIds?.join(', ') || '—'}</dd></div><div><dt>Transcript effect</dt><dd>{snp?.transcriptConsequences?.join(', ').replaceAll('_', ' ') || '—'}</dd></div><div><dt>Alleles / population MAF</dt><dd>{snp?.grch38?.reference || '—'} → {snp?.grch38?.alternates?.join(', ') || '—'}{Number.isFinite(snp?.maf) ? ` · MAF ${snp.maf.toFixed(4)}` : ''}</dd></div><div><dt>Study MAF</dt><dd>{Number.isFinite(snp?.studyMaf) ? `${snp.studyMaf.toFixed(4)} · ${snp.studyMafClass || 'unclassified'}` : 'Not recorded'}</dd></div></dl><footer>{snp?.ucscUrl && <a href={snp.ucscUrl} target="_blank" rel="noreferrer">UCSC locus ↗</a>}{snp?.ensemblUrl && <a href={snp.ensemblUrl} target="_blank" rel="noreferrer">Ensembl variant ↗</a>}</footer></article>
  })}</div>
}

function canonicalGrch38Variant(snp) {
  const chromosome = snp?.grch38?.chromosome
  const position = snp?.grch38?.position
  const reference = snp?.grch38?.reference
  const alternate = snp?.grch38?.alternates?.[0]
  return chromosome && position && reference && alternate ? `${chromosome}:${position}-${reference}-${alternate}` : null
}

function PhewasEvidence({ snp }) {
  const canonical = canonicalGrch38Variant(snp)
  const [state, setState] = useState({ loading: Boolean(canonical), ukb: null, finngen: null })
  useEffect(() => {
    if (!canonical) {
      setState({ loading: false, ukb: null, finngen: null })
      return undefined
    }
    const cached = phewasCache.get(canonical)
    if (cached) {
      setState(cached)
      return undefined
    }
    const controller = new AbortController()
    const query = async (url, key) => {
      const response = await fetch(url, { signal: controller.signal })
      if (!response.ok) throw new Error(String(response.status))
      const payload = await response.json()
      const rows = Array.isArray(payload[key]) ? payload[key] : []
      return rows.filter((row) => { const p = finiteNumber(row.pval); return p !== null && p < 0.05 }).sort((left, right) => finiteNumber(left.pval) - finiteNumber(right.pval))
    }
    setState({ loading: true, ukb: null, finngen: null })
    Promise.allSettled([
      query(`https://pheweb.org/UKB-TOPMed/api/variant/${canonical}`, 'phenos'),
      query(`https://r12.finngen.fi/api/variant/${canonical}`, 'results'),
    ]).then(([ukb, finngen]) => {
      if (controller.signal.aborted) return
      const result = { loading: false, ukb: ukb.status === 'fulfilled' ? ukb.value : { error: true }, finngen: finngen.status === 'fulfilled' ? finngen.value : { error: true } }
      phewasCache.set(canonical, result)
      setState(result)
    })
    return () => controller.abort()
  }, [canonical])

  if (!canonical) return <p className="locus-empty">A canonical GRCh38 REF/ALT representation is required for public PheWAS lookup.</p>
  if (state.loading) return <p className="locus-empty">Retrieving UKB-TOPMed and FinnGen R12 associations for {canonical}…</p>
  const renderSource = (label, rows, href) => <section className="locus-gwas-source"><header><h3>{label}</h3><a href={href} target="_blank" rel="noreferrer">Open cohort page ↗</a></header>{rows?.error ? <p>Source could not be reached from this browser.</p> : rows?.length ? <div className="locus-table-wrap"><table className="locus-table"><thead><tr><th>Phenotype</th><th>Code</th><th>β</th><th>P</th><th>N</th></tr></thead><tbody>{rows.slice(0, 20).map((row, index) => { const beta = finiteNumber(row.beta); const sampleSize = finiteNumber(row.num_samples ?? row.n_sample); return <tr key={`${row.phenocode || row.phenostring}-${index}`}><td>{row.phenostring || row.phenocode || 'Phenotype'}</td><td>{row.phenocode || '—'}</td><td className={beta === null ? '' : beta >= 0 ? 'positive' : 'negative'}>{beta === null ? '—' : `${beta >= 0 ? '+' : ''}${beta.toPrecision(3)}`}</td><td>{formatP(finiteNumber(row.pval))}</td><td>{sampleSize === null ? '—' : number.format(sampleSize)}</td></tr> })}</tbody></table></div> : <p>No nominal association at P &lt; 0.05 was returned.</p>}</section>
  return <div className="locus-tab-body"><p className="locus-caveat">Cohort-specific associations are shown for prioritisation, not as evidence of causality. The top 20 nominal associations per source are displayed without additional re-correction.</p>{renderSource('UKB-TOPMed', state.ukb, `https://pheweb.org/UKB-TOPMed/variant/${canonical}`)}{renderSource('FinnGen R12', state.finngen, `https://r12.finngen.fi/variant/${canonical}`)}</div>
}

function AlphaGenomeTab({ entity, activeCell, snpById, alphaIndex, alphaCoverage, onOpenPrediction }) {
  const alphaBySnp = new Map(alphaIndex.map((row) => [row.snpId, row]))
  return <div className="locus-tab-body"><p className="locus-caveat">AlphaGenome predictions compare ALT with REF. Signed and unsigned scorer definitions remain separate; values are model predictions, not association P values.</p>{alphaCoverage && <div className="locus-alpha-coverage"><b>{number.format(alphaCoverage.scoredSnps || 0)}</b> of {number.format(alphaCoverage.totalInputSnps || 0)} GIM SNPs scored in the local batch.</div>}<div className="locus-alpha-list">{entity.snps.map((snpId) => {
    const summary = alphaBySnp.get(snpId)
    const snp = snpById.get(snpId)
    const contextCell = activeCell?.snpId === snpId ? activeCell : entity.heatmap.find((cell) => cell.direct && cell.snpId === snpId) || { snpId, metabolite: entity.metabolites[0] }
    return <article key={snpId}><header><div><b>{snpId}</b><span>{snp?.geneSymbols?.join(', ') || 'Unassigned gene'}</span></div><small>{variantCoordinate(snp)}</small></header><div><span><b>{number.format(summary?.nTracks || 0)}</b> target tracks</span><span><b>{number.format(summary?.nGastricTissueTracks || 0)}</b> stomach</span><span><b>{number.format(summary?.nGastricCancerTracks || 0)}</b> gastric cancer</span><span><b>{number.format(summary?.nImmuneTracks || 0)}</b> immune</span></div><button disabled={!summary?.nTracks} onClick={() => onOpenPrediction(contextCell)}>Open multimodal prediction <span>{summary?.nTracks ? '↗' : 'pending'}</span></button></article>
  })}</div></div>
}

function MetabolomicsTab({ entity, cancerAnnotations }) {
  return <div className="locus-tab-body"><p className="locus-caveat">Matched LC–MS case–control results are independent observational evidence. Individual-level genotype boxplots are not reconstructed because the portable GIMForge bundle intentionally excludes participant-level data.</p><div className="locus-metabolomics-list">{entity.metabolites.map((metabolite) => {
    const rows = cancerAnnotations[metabolite] || []
    return <article key={metabolite}><header><b>{metabolite}</b><span>{rows.length ? `${rows.length} comparisons` : 'No matched comparison'}</span></header>{rows.length ? <div>{rows.map((row) => <p key={row.comparison} className={row.fdrSignificant ? 'significant' : ''}><span>{row.comparison?.replaceAll('_', ' ')}</span><b className={row.betaCaseMinusOtherSd >= 0 ? 'positive' : 'negative'}>{row.betaCaseMinusOtherSd >= 0 ? '+' : ''}{row.betaCaseMinusOtherSd?.toFixed(2)}</b><small>FDR {formatP(row.fdr)}</small></p>)}</div> : <p>No gastric-cancer differential-metabolomics record was found.</p>}</article>
  })}</div></div>
}

function GimLocusPanel({ selection, snpById, alphaIndex, alphaCoverage, cancerAnnotations, onSelectCell, onOpenPrediction, onClose }) {
  const [tab, setTab] = useState('summary')
  const [activeCell, setActiveCell] = useState(selection.cell)
  useEffect(() => {
    setActiveCell(selection.cell)
    setTab('summary')
  }, [selection])
  const entity = selection.entity
  const region = selection.region
  const activeSnp = snpById.get(activeCell?.snpId || entity.snps[0])
  const chooseCell = (cell) => {
    setActiveCell(cell)
    setTab('summary')
    onSelectCell?.(cell)
  }
  const tabs = [
    ['summary', 'Summary stats'],
    ['genes', 'Genes'],
    ['gwas', 'Other GWAS'],
    ['alpha', 'AlphaGenome'],
    ['metabolomics', 'Metabolomics'],
  ]
  return <aside className="locus-panel" aria-label="Selected GIM and locus annotations">
    <header className="locus-panel-header"><div><p className="atlas-eyebrow">SELECTED LOCUS / GIM</p><h2>{selection.genes.join('/') || 'Unassigned'} · {formatGimId(entity.gimId)}</h2><p>{region?.cytoband || 'Cytoband unresolved'} · {coordinateLabel(region)} · {entity.nSnps} SNPs · {entity.nMetabolites} metabolites</p></div><button onClick={onClose} aria-label="Close locus details">×</button></header>
    <div className="locus-panel-scroll">
      <SmallGimHeatmap entity={entity} activeCell={activeCell} snpById={snpById} onSelect={chooseCell} />
      <p className="locus-selection-note">Selected global association: <b>{selection.metabolite}</b>. Click a coloured SNP×metabolite square for variant-level annotation.</p>
      <nav className="locus-tabs" aria-label="Locus annotation sections">{tabs.map(([value, label]) => <button className={tab === value ? 'active' : ''} key={value} onClick={() => setTab(value)}>{label}</button>)}</nav>
      {tab === 'summary' && <SummaryTab entity={entity} activeCell={activeCell} snpById={snpById} onSelect={chooseCell} />}
      {tab === 'genes' && <GenesTab entity={entity} snpById={snpById} region={region} />}
      {tab === 'gwas' && <PhewasEvidence snp={activeSnp} />}
      {tab === 'alpha' && <AlphaGenomeTab entity={entity} activeCell={activeCell} snpById={snpById} alphaIndex={alphaIndex} alphaCoverage={alphaCoverage} onOpenPrediction={onOpenPrediction} />}
      {tab === 'metabolomics' && <MetabolomicsTab entity={entity} cancerAnnotations={cancerAnnotations} />}
    </div>
  </aside>
}

export default function GimAtlas({ data, snpById, focusGim, onFocusGim, onOpenPrediction }) {
  const [selection, setSelection] = useState(null)
  const regionById = useMemo(() => new Map(data.regions.map((region) => [region.regionId, region])), [data.regions])
  useEffect(() => {
    if (!focusGim || selection?.gimId === focusGim) return
    const entity = data.gimEntities.find((item) => item.gimId === focusGim)
    if (!entity) return
    setSelection(buildSelection(entity, entity.metabolites[0], regionById.get(entity.regionId) || {}, snpById))
  }, [data.gimEntities, focusGim, regionById, selection?.gimId, snpById])
  const selectGlobal = (next) => {
    setSelection(next)
    onFocusGim?.(next.gimId)
  }
  const close = () => {
    setSelection(null)
    onFocusGim?.(null)
  }
  return <section className={`gim-atlas ${selection ? 'has-locus-panel' : ''}`}>
    <GlobalGimMatrix entities={data.gimEntities} regions={data.regions} snpById={snpById} selection={selection} onSelect={selectGlobal} />
    {selection && <GimLocusPanel selection={selection} snpById={snpById} alphaIndex={data.alphaGenomeIndex} alphaCoverage={data.alphaCoverage} cancerAnnotations={data.cancerAnnotations} onOpenPrediction={onOpenPrediction} onClose={close} />}
  </section>
}
