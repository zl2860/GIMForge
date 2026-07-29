import { useEffect, useMemo, useRef, useState } from 'react'
import './App.css'
import GimAtlas from './GimAtlas'

const number = new Intl.NumberFormat('en-US')
const phewasCache = new Map()
const alphaSignalCache = new Map()

function formatP(value) {
  if (!Number.isFinite(value)) return '—'
  return value < 0.001 ? value.toExponential(2).replace('e-', ' × 10⁻') : value.toFixed(4)
}

function finiteNumber(value) {
  if (value === null || value === undefined || value === '') return null
  const parsed = Number(value)
  return Number.isFinite(parsed) ? parsed : null
}

function formatGimId(id) {
  const match = id?.match(/(?:^|_)region_(\d+)_GIM_(\d+)$/)
  return match ? `GIM ${Number(match[1])}.${Number(match[2])}` : id || 'GIM'
}

function chromosomeOrder(left, right) {
  const leftNumber = Number(left)
  const rightNumber = Number(right)
  if (Number.isFinite(leftNumber) && Number.isFinite(rightNumber)) return leftNumber - rightNumber
  if (Number.isFinite(leftNumber)) return -1
  if (Number.isFinite(rightNumber)) return 1
  return left.localeCompare(right)
}

function formatCoordinate(snp) {
  const coordinate = snp?.grch38?.position
  if (coordinate) return `chr${snp.grch38.chromosome}:${number.format(coordinate)} (GRCh38)`
  if (snp?.grch37?.grch37Position) return `chr${snp.grch37.chromosome}:${number.format(snp.grch37.grch37Position)} (GRCh37)`
  return 'Coordinate unresolved'
}

function outputLabel(value) {
  return value?.replaceAll('_', ' ') || 'Prediction'
}

function scoreIsInScope(score, scope) {
  const scopes = Array.isArray(score.scopes) ? score.scopes : []
  return scope === 'target' ? scopes.length > 0 : scopes.includes(scope) || scopes.includes('shared')
}

function predictionScoreValue(score) {
  return Number.isFinite(score.rawScore) ? score.rawScore : score.score
}

function predictionRankingValue(score) {
  return Number.isFinite(score.rankingScore) ? score.rankingScore : null
}

function combinedSplicingScore(scores) {
  const grouped = new Map()
  const components = {
    SPLICE_SITES: 'spliceSites',
    SPLICE_SITE_USAGE: 'spliceSiteUsage',
    SPLICE_JUNCTIONS: 'spliceJunctions',
  }
  scores.forEach((score) => {
    const component = components[score.outputType]
    const value = predictionScoreValue(score)
    if (!component || !Number.isFinite(value)) return
    const key = `${score.alternate || ''}::${score.geneId || score.geneName || 'unassigned'}`
    const item = grouped.get(key) || { alternate: score.alternate, geneId: score.geneId, geneName: score.geneName }
    if (!Number.isFinite(item[component]) || Math.abs(value) > Math.abs(item[component])) item[component] = value
    grouped.set(key, item)
  })
  return [...grouped.values()].map((item) => ({
    ...item,
    score: Math.abs(item.spliceSites || 0) + Math.abs(item.spliceSiteUsage || 0) + Math.abs(item.spliceJunctions || 0) / 5,
  })).sort((left, right) => right.score - left.score)[0] || null
}

function scorerDefinition(score) {
  const scorer = score.variantScorer || ''
  if (scorer.includes('GeneMaskLFCScorer')) return { label: 'Gene-expression log fold-change', formula: 'log(mean ALT + 0.001) − log(mean REF + 0.001)', signed: true }
  if (scorer.includes('CenterMaskScorer') && scorer.includes('DIFF_LOG2_SUM')) return { label: 'Local log₂ fold-change', formula: 'log₂[(ΣALT + 1) / (ΣREF + 1)]', signed: true }
  if (scorer.includes('GeneMaskActiveScorer')) return { label: 'Active-allele expression signal', formula: 'max(mean ALT, mean REF)', signed: false }
  if (scorer.includes('CenterMaskScorer') && scorer.includes('ACTIVE_SUM')) return { label: 'Active-allele local signal', formula: 'max(ΣALT, ΣREF)', signed: false }
  if (scorer.includes('GeneMaskSplicingScorer')) return { label: 'Maximum splice change', formula: 'max(|ALT − REF|) across the gene body', signed: false }
  if (scorer.includes('SpliceJunctionScorer')) return { label: 'Maximum junction log fold-change', formula: 'max(|log ALT − log REF|)', signed: false }
  if (scorer.includes('ContactMapScorer')) return { label: 'Local contact disruption', formula: 'mean(|ALT − REF|) for contacts at the variant bin', signed: false }
  if (scorer.includes('PolyadenylationScorer')) return { label: 'Maximum PAS-usage log fold-change', formula: 'max(|log distal/proximal ALT − log distal/proximal REF|)', signed: false }
  return { label: 'Model score', formula: 'Scorer metadata is being refreshed', signed: false, legacy: true }
}

function downloadAlphaGenomeManifest() {
  const anchor = document.createElement('a')
  anchor.href = '/data/alphagenome_input.tsv'
  anchor.download = 'gims_alphagenome_grch38_input.tsv'
  document.body.appendChild(anchor)
  anchor.click()
  anchor.remove()
}

function downloadGimLocusTable() {
  const anchor = document.createElement('a')
  anchor.href = '/data/gim_locus_metabolites_summary.csv'
  anchor.download = 'gim_locus_metabolites_summary.csv'
  document.body.appendChild(anchor)
  anchor.click()
  anchor.remove()
}

function App() {
  const [data, setData] = useState({ records: [], snps: [], gimEntities: [], regions: [], stats: null, alphaGenomeIndex: [], alphaCoverage: null, cancerAnnotations: {} })
  const [portalView, setPortalView] = useState('gims')
  const [selectedGim, setSelectedGim] = useState(null)
  const [predictionCell, setPredictionCell] = useState(null)
  const [predictionScores, setPredictionScores] = useState([])
  const [predictionStatus, setPredictionStatus] = useState({ loading: false, error: null })
  const predictionRequestRef = useRef(0)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)

  useEffect(() => {
    Promise.all([
      fetch('/data/gim_records.json').then((response) => response.json()),
      fetch('/data/gim_snps.json').then((response) => response.json()),
      fetch('/data/gim_entities.json').then((response) => response.json()),
      fetch('/data/gim_regions.json').then((response) => response.json()),
      fetch('/data/portal_stats.json').then((response) => response.json()),
      fetch('/data/alphagenome_index.json').then((response) => response.ok ? response.json() : []).catch(() => []),
      fetch('/data/alphagenome_coverage.json').then((response) => response.ok ? response.json() : null).catch(() => null),
      fetch('/data/metabolite_cancer_annotations.json').then((response) => response.ok ? response.json() : {}).catch(() => ({})),
    ])
      .then(([records, snps, gimEntities, regions, stats, alphaGenomeIndex, alphaCoverage, cancerAnnotations]) => {
        setData({ records, snps, gimEntities, regions, stats, alphaGenomeIndex, alphaCoverage, cancerAnnotations })
      })
      .catch(() => setError('Portal data is not available. Run npm run portal:build to rebuild the data package.'))
      .finally(() => setLoading(false))
  }, [])

  const snpById = useMemo(() => new Map(data.snps.map((snp) => [snp.snpId, snp])), [data.snps])

  const closePrediction = () => {
    predictionRequestRef.current += 1
    setPredictionCell(null)
    setPredictionScores([])
    setPredictionStatus({ loading: false, error: null })
  }

  const openPrediction = (cell) => {
    const requestId = predictionRequestRef.current + 1
    predictionRequestRef.current = requestId
    setPredictionCell(cell)
    setPredictionScores([])
    setPredictionStatus({ loading: true, error: null })
    fetch(`/data/alphagenome_scores/${encodeURIComponent(cell.snpId)}.json`)
      .then((response) => {
        if (!response.ok) throw new Error(`Prediction file returned ${response.status}`)
        return response.json()
      })
      .then((scores) => {
        if (predictionRequestRef.current !== requestId) return
        setPredictionScores(Array.isArray(scores) ? scores : [])
        setPredictionStatus({ loading: false, error: null })
      })
      .catch(() => {
        if (predictionRequestRef.current !== requestId) return
        setPredictionStatus({ loading: false, error: 'The detailed AlphaGenome track file is not available for this variant.' })
      })
  }

  if (loading) return <main className="loading-page">Loading GIM matrix…</main>
  if (error) return <main className="loading-page error-state">{error}</main>

  return <main className="portal-app">
    <header className="masthead">
      <div className="brand" aria-label="GIMs portal"><span className="brand-mark"><i></i><i></i><i></i></span><span>GIMs <em>Portal</em></span></div>
      <nav className="portal-nav" aria-label="Portal sections"><button className={portalView === 'gims' ? 'active' : ''} onClick={() => setPortalView('gims')}>GIM heatmap</button><button className={portalView === 'variants' ? 'active' : ''} onClick={() => setPortalView('variants')}>Variant explorer</button></nav>
      <p className="masthead-caption">Conditional genetic influence modules in metabolite genetics</p>
      <div className="masthead-metrics" aria-label="Dataset summary"><span><b>{number.format(data.stats?.nGimEntities || data.gimEntities.length)}</b> GIMs</span><span><b>{number.format(data.stats?.nIndependentSnps || data.snps.length)}</b> independent SNPs</span><button onClick={downloadGimLocusTable}>Locus–GIM table ↓</button><button onClick={downloadAlphaGenomeManifest}>AlphaGenome input ↓</button></div>
    </header>

    {portalView === 'gims' && <GimAtlas data={data} snpById={snpById} focusGim={selectedGim} onFocusGim={setSelectedGim} onOpenPrediction={openPrediction} />}
    {portalView === 'variants' && <VariantExplorer snps={data.snps} records={data.records} gimEntities={data.gimEntities} alphaIndex={data.alphaGenomeIndex} alphaCoverage={data.alphaCoverage} onOpenPrediction={(snpId) => { const firstAssociation = data.records.find((record) => record.snpId === snpId); openPrediction({ snpId, metabolite: firstAssociation?.metabolite || 'Variant functional profile' }) }} onOpenGim={(gimId) => { setSelectedGim(gimId); setPortalView('gims') }} />}
    {predictionCell && <VariantFunctionalPredictionDialog cell={predictionCell} snp={snpById.get(predictionCell.snpId)} scores={predictionScores} loading={predictionStatus.loading} error={predictionStatus.error} onClose={closePrediction} />}
  </main>
}

function canonicalGrch38Variant(snp) {
  const chromosome = snp?.grch38?.chromosome
  const position = snp?.grch38?.position
  const reference = snp?.grch38?.reference
  const alternate = snp?.alphaGenome?.input?.alternates?.[0] || snp?.grch38?.studyAlternates?.[0]
  return chromosome && position && reference && alternate ? `${chromosome}:${position}-${reference}-${alternate}` : null
}

function formatEffect(value) {
  if (!Number.isFinite(value)) return '—'
  return `${value >= 0 ? '+' : ''}${value.toPrecision(3)}`
}

function PhewasPlot({ source, result }) {
  const [hovered, setHovered] = useState(null)
  if (result?.error) return <div className="phewas-message">{source} could not be reached for this variant.</div>
  if (!result) return <div className="phewas-message">No association payload returned.</div>
  if (!result.items.length) return <div className="phewas-message">No nominally significant association (P &lt; 0.05) was reported for this variant.</div>
  const effectLimit = Math.max(...result.items.map((row) => Math.abs(finiteNumber(row.beta) ?? 0)), 0.1)
  const significanceLimit = Math.max(...result.items.map((row) => -Math.log10(Math.max(finiteNumber(row.pval) ?? 1, Number.MIN_VALUE))), 2)
  const formattedSignificanceLimit = significanceLimit >= 10 ? significanceLimit.toFixed(0) : significanceLimit.toFixed(1)
  return <>
    <p className="phewas-count">{number.format(result.items.length)} nominally significant associations (P &lt; 0.05; uncorrected) · hover a point for phenotype detail</p>
    <div className="phewas-plot-wrap">
      <div className="phewas-plot" role="group" aria-label={`${source} effect-size and P-value scatter plot`}>
        <span className="phewas-y-label">−log10 P</span><span className="phewas-y-max">{formattedSignificanceLimit}</span><span className="phewas-y-zero">0</span><span className="phewas-x-label">Protective β ← · → Risk β</span><i className="phewas-zero-line"></i><i className="phewas-floor-line"></i>
        {result.items.map((row, index) => {
          const beta = finiteNumber(row.beta)
          const plottedBeta = beta ?? 0
          const p = finiteNumber(row.pval)
          const negLogP = -Math.log10(Math.max(p ?? 1, Number.MIN_VALUE))
          const left = 8 + ((plottedBeta + effectLimit) / (2 * effectLimit)) * 84
          const bottom = 15 + Math.min(negLogP / significanceLimit, 1) * 76
          const sampleSize = row.num_samples ?? row.n_sample
          const caseCount = row.num_cases ?? row.n_case
          const controlCount = row.num_controls ?? row.n_control
          const detail = { phenotype: row.phenostring || row.phenocode || 'Phenotype', phenocode: row.phenocode, beta, p, category: row.category || 'Uncategorised', sampleSize: finiteNumber(sampleSize), caseCount: finiteNumber(caseCount), controlCount: finiteNumber(controlCount), left, bottom }
          const label = `${detail.phenotype}; beta ${formatEffect(beta)}; P ${formatP(detail.p)}`
          return <button key={`${row.phenocode || row.phenostring}-${index}`} className={`phewas-point ${beta === null ? '' : beta >= 0 ? 'positive' : 'negative'}`} style={{ left: `${left}%`, bottom: `${bottom}%` }} onMouseEnter={() => setHovered(detail)} onFocus={() => setHovered(detail)} onClick={() => setHovered(detail)} aria-label={label} aria-describedby={hovered?.phenocode === detail.phenocode ? `${source}-phewas-tooltip` : undefined}></button>
        })}
      </div>
      {hovered && <aside id={`${source}-phewas-tooltip`} className={`phewas-tooltip ${hovered.left > 55 ? 'right' : ''}`} role="status" aria-live="polite" style={{ left: `${hovered.left}%`, top: hovered.bottom > 55 ? '8px' : '46px' }}><b>{hovered.phenotype}</b>{hovered.phenocode && <span>{hovered.phenocode}</span>}<dl><div><dt>β</dt><dd className={hovered.beta === null ? '' : hovered.beta >= 0 ? 'positive' : 'negative'}>{formatEffect(hovered.beta)}</dd></div><div><dt>P</dt><dd>{formatP(hovered.p)}</dd></div><div><dt>Category</dt><dd>{hovered.category}</dd></div>{hovered.sampleSize !== null && <div><dt>N</dt><dd>{number.format(hovered.sampleSize)}</dd></div>}{hovered.caseCount !== null && <div><dt>Cases</dt><dd>{number.format(hovered.caseCount)}{hovered.controlCount !== null ? ` / ${number.format(hovered.controlCount)} controls` : ''}</dd></div>}</dl></aside>}
    </div>
  </>
}

function VariantExplorer({ snps, records, gimEntities, alphaIndex, alphaCoverage, onOpenPrediction, onOpenGim }) {
  const [query, setQuery] = useState('')
  const [chromosome, setChromosome] = useState('all')
  const [priorityFilter, setPriorityFilter] = useState('all')
  const [alphaScope, setAlphaScope] = useState('target')
  const [alphaModality, setAlphaModality] = useState('all')
  const [minimumModalities, setMinimumModalities] = useState(0)
  const [sortMode, setSortMode] = useState('multimodal')
  const [selectedSnpId, setSelectedSnpId] = useState('')
  const [phewas, setPhewas] = useState({ loading: false, ukb: null, finngen: null })

  const variants = useMemo(() => {
    const alphaBySnp = new Map(alphaIndex.map((item) => [item.snpId, item]))
    const bySnp = new Map()
    records.forEach((record) => {
      const item = bySnp.get(record.snpId) || { metabolites: new Set(), regions: new Set(), minP: Infinity, associations: [] }
      item.metabolites.add(record.metabolite)
      item.regions.add(record.regionId)
      item.minP = Math.min(item.minP, record.p)
      item.associations.push(record)
      bySnp.set(record.snpId, item)
    })
    const gimsBySnp = new Map()
    gimEntities.forEach((gim) => gim.snps.forEach((snpId) => {
      const ids = gimsBySnp.get(snpId) || []
      ids.push(gim.gimId)
      gimsBySnp.set(snpId, ids)
    }))
    return snps.map((snp) => {
      const associations = bySnp.get(snp.snpId) || { metabolites: new Set(), regions: new Set(), minP: Infinity, associations: [] }
      return { snp, ...associations, alpha: alphaBySnp.get(snp.snpId), gimIds: [...new Set(gimsBySnp.get(snp.snpId) || [])].sort() }
    }).sort((left, right) => left.snp.snpId.localeCompare(right.snp.snpId, undefined, { numeric: true }))
  }, [alphaIndex, gimEntities, records, snps])

  const chromosomes = useMemo(() => [...new Set(variants.map((item) => String(item.snp.grch38?.chromosome || item.snp.grch37?.chromosome || 'Unresolved')))].sort(chromosomeOrder), [variants])
  const alphaModalities = useMemo(() => [...new Set(alphaIndex.flatMap((item) => item.modalities || []))].sort(), [alphaIndex])
  const filteredVariants = useMemo(() => {
    const needle = query.trim().toLowerCase()
    const rows = variants.filter((item) => {
      const itemChromosome = String(item.snp.grch38?.chromosome || item.snp.grch37?.chromosome || 'Unresolved')
      if (chromosome !== 'all' && itemChromosome !== chromosome) return false
      const consequence = item.snp.mostSevereConsequence || ''
      if (priorityFilter === 'coding' && !/(missense|synonymous|stop_|start_|splice|frameshift|protein_altering)/i.test(consequence)) return false
      if (priorityFilter === 'multiGim' && item.gimIds.length < 2) return false
      if (priorityFilter === 'multiMetabolite' && item.metabolites.size < 2) return false
      const scope = item.alpha?.scopes?.[alphaScope] || {
        nTracks: alphaScope === 'target' ? item.alpha?.nTracks || 0 : 0,
        nModalities: alphaScope === 'target' ? item.alpha?.modalities?.length || 0 : 0,
        modalities: alphaScope === 'target' ? item.alpha?.modalities || [] : [],
        maxAbsRankingScore: alphaScope === 'target' ? item.alpha?.maxAbsRankingScore : null,
      }
      if (alphaScope !== 'target' && !scope.nTracks) return false
      if (alphaModality !== 'all' && !scope.modalities?.includes(alphaModality)) return false
      if ((scope.nModalities || 0) < minimumModalities) return false
      if (!needle) return true
      return [item.snp.snpId, ...item.snp.geneSymbols, ...item.snp.geneIds, ...item.snp.transcriptIds, itemChromosome, item.snp.mostSevereConsequence]
        .filter(Boolean).join(' ').toLowerCase().includes(needle)
    })
    const scopeFor = (item) => item.alpha?.scopes?.[alphaScope] || item.alpha || {}
    return rows.sort((left, right) => {
      const leftScope = scopeFor(left)
      const rightScope = scopeFor(right)
      if (sortMode === 'evidence') return (rightScope.maxAbsRankingScore ?? -1) - (leftScope.maxAbsRankingScore ?? -1) || (rightScope.nModalities || 0) - (leftScope.nModalities || 0)
      if (sortMode === 'splicing') return (rightScope.splicingCombined?.score || 0) - (leftScope.splicingCombined?.score || 0) || (rightScope.maxAbsRankingScore ?? -1) - (leftScope.maxAbsRankingScore ?? -1)
      if (sortMode === 'association') return left.minP - right.minP
      if (sortMode === 'metabolites') return right.metabolites.size - left.metabolites.size || left.minP - right.minP
      return (rightScope.nModalities || 0) - (leftScope.nModalities || 0) || (rightScope.maxAbsRankingScore ?? -1) - (leftScope.maxAbsRankingScore ?? -1) || left.minP - right.minP
    })
  }, [alphaModality, alphaScope, chromosome, minimumModalities, priorityFilter, query, sortMode, variants])

  useEffect(() => {
    if (!filteredVariants.length) return
    if (!filteredVariants.some((item) => item.snp.snpId === selectedSnpId)) setSelectedSnpId(filteredVariants[0].snp.snpId)
  }, [filteredVariants, selectedSnpId])

  const selected = filteredVariants.find((item) => item.snp.snpId === selectedSnpId) || filteredVariants[0] || variants[0]
  const alphaSummary = selected.alpha
  const nPredictedTracks = alphaSummary?.nTracks || 0
  const gastricTissueScoreCount = alphaSummary?.nGastricTissueTracks || 0
  const gastricCancerScoreCount = alphaSummary?.nGastricCancerTracks || 0
  const immuneScoreCount = alphaSummary?.nImmuneTracks || 0
  const canonical = canonicalGrch38Variant(selected?.snp)

  useEffect(() => {
    if (!canonical) {
      setPhewas({ loading: false, ukb: null, finngen: null })
      return undefined
    }
    const cached = phewasCache.get(canonical)
    if (cached) {
      setPhewas(cached)
      return undefined
    }
    const controller = new AbortController()
    const queryPhewas = async (url, payloadKey) => {
      const response = await fetch(url, { signal: controller.signal })
      if (!response.ok) throw new Error(`PheWAS query returned ${response.status}`)
      const payload = await response.json()
      const items = Array.isArray(payload[payloadKey]) ? payload[payloadKey] : []
      return { total: items.length, items: items.filter((item) => { const p = finiteNumber(item.pval); return p !== null && p < 0.05 }).sort((left, right) => finiteNumber(left.pval) - finiteNumber(right.pval)) }
    }
    setPhewas({ loading: true, ukb: null, finngen: null })
    Promise.allSettled([
      queryPhewas(`https://pheweb.org/UKB-TOPMed/api/variant/${canonical}`, 'phenos'),
      queryPhewas(`https://r12.finngen.fi/api/variant/${canonical}`, 'results'),
    ]).then(([ukb, finngen]) => {
      if (controller.signal.aborted) return
      const result = {
        loading: false,
        ukb: ukb.status === 'fulfilled' ? ukb.value : { error: true },
        finngen: finngen.status === 'fulfilled' ? finngen.value : { error: true },
      }
      phewasCache.set(canonical, result)
      setPhewas(result)
    })
    return () => controller.abort()
  }, [canonical])

  if (!selected) return null
  const snp = selected.snp
  const genes = snp.geneSymbols.length ? snp.geneSymbols.join(', ') : snp.geneIds.length ? snp.geneIds.join(', ') : 'No overlapping gene assigned'
  const alt = snp.alphaGenome?.input?.alternates?.join(', ') || snp.grch38?.studyAlternates?.join(', ') || '—'

  return <section className="variant-workbench" aria-label="Variant explorer">
    <aside className="variant-index">
      <header><p className="eyebrow">GIM VARIANT INDEX</p><h1>Variant explorer</h1><p>Search an rsID, gene, transcript, consequence, or chromosome.</p></header>
      <div className="variant-filters"><input value={query} onChange={(event) => setQuery(event.target.value)} placeholder="e.g. rs10001483 or ENSG…" aria-label="Search GIM variants" /><select value={chromosome} onChange={(event) => setChromosome(event.target.value)} aria-label="Filter variant chromosome"><option value="all">All chromosomes</option>{chromosomes.map((item) => <option key={item} value={item}>{item === 'Unresolved' ? item : `chr${item}`}</option>)}</select></div>
      <section className="alpha-priority-panel">
        <header><div><p className="eyebrow">ALPHAGENOME MULTIMODAL PRIORITISATION</p><h2>Filter and rank variants</h2></div><span>Explainable rules</span></header>
        <div className="alpha-priority-grid">
          <label>Biological scope<select value={alphaScope} onChange={(event) => setAlphaScope(event.target.value)}><option value="target">All retained target scopes</option><option value="gastric_tissue">Stomach tissue</option><option value="gastric_cancer">Gastric cancer</option><option value="immune">Immune cells</option></select></label>
          <label>Required modality<select value={alphaModality} onChange={(event) => setAlphaModality(event.target.value)}><option value="all">Any modality</option>{alphaModalities.map((item) => <option key={item} value={item}>{outputLabel(item)}</option>)}</select></label>
          <label>Minimum modalities<select value={minimumModalities} onChange={(event) => setMinimumModalities(Number(event.target.value))}>{[0, 2, 3, 4, 5, 6, 7, 8].map((value) => <option key={value} value={value}>{value ? `${value}+ modalities` : 'No minimum'}</option>)}</select></label>
          <label>Sort by<select value={sortMode} onChange={(event) => setSortMode(event.target.value)}><option value="multimodal">Multimodal breadth</option><option value="evidence">Absolute model quantile</option><option value="splicing">Official merged splicing</option><option value="association">Strongest GIM P value</option><option value="metabolites">Metabolite breadth</option></select></label>
        </div>
        <p>“Multimodal breadth” ranks by distinct output types, then the absolute scorer/track quantile, then the GIM association P value. Strong increases and decreases are treated symmetrically. The splicing option uses AlphaGenome’s recommended merged formula.</p>
      </section>
      <div className="variant-priority-filters" aria-label="Functional variant filters">{[
        ['all', 'All genomic contexts'], ['coding', 'Coding / splice'], ['multiGim', 'Multi-GIM'], ['multiMetabolite', 'Multi-metabolite'],
      ].map(([value, label]) => <button className={priorityFilter === value ? 'active' : ''} onClick={() => setPriorityFilter(value)} key={value}>{label}</button>)}</div>
      <p className="variant-result-count">{number.format(filteredVariants.length)} of {number.format(variants.length)} GIM SNPs</p>
      <div className="variant-list">{filteredVariants.map((item) => {
        const summary = item.alpha?.scopes?.[alphaScope] || item.alpha || {}
        return <button className={item.snp.snpId === selectedSnpId ? 'active' : ''} key={item.snp.snpId} onClick={() => setSelectedSnpId(item.snp.snpId)}>
          <span className="variant-list-title"><b>{item.snp.snpId}</b><i>{summary.nModalities || 0} modalities{Number.isFinite(summary.maxAbsRankingScore) ? ` · |Q| ${summary.maxAbsRankingScore.toFixed(4)}` : ''}</i></span><span>{formatCoordinate(item.snp)}</span><small>{item.snp.geneSymbols.join(', ') || item.snp.geneIds.join(', ') || item.snp.mostSevereConsequence?.replaceAll('_', ' ') || 'unannotated'}</small>
        </button>
      })}</div>
    </aside>
    <article className="variant-detail">
      <header className="variant-hero"><div><p className="eyebrow">GIM-ASSOCIATED VARIANT</p><h1>{snp.snpId}</h1><p>{formatCoordinate(snp)} · {snp.grch38?.reference || '—'} → {alt}</p></div><div className="variant-hero-actions">{snp.ucscUrl && <a href={snp.ucscUrl} target="_blank" rel="noreferrer">UCSC locus ↗</a>}{snp.ensemblUrl && <a href={snp.ensemblUrl} target="_blank" rel="noreferrer">Ensembl variant ↗</a>}</div></header>
      <div className="variant-stat-strip"><span><b>{selected.gimIds.length}</b> GIMs</span><span><b>{selected.metabolites.size}</b> metabolites</span><span><b>{selected.associations.length}</b> conditional tests</span><span><b>{nPredictedTracks}</b> target tracks</span></div>
      <div className="variant-detail-grid">
        <section className="variant-card"><h2>Genomic context</h2><dl className="annotation-list"><div><dt>Gene / locus</dt><dd>{genes}</dd></div><div><dt>Class / alleles</dt><dd>{snp.variantClass || 'Variant'} · {snp.grch38?.reference || '—'} → {alt}</dd></div><div><dt>Consequence</dt><dd>{snp.mostSevereConsequence?.replaceAll('_', ' ') || 'Not resolved'}</dd></div><div><dt>Transcripts</dt><dd>{snp.transcriptIds.join(', ') || 'No transcript overlap'}</dd></div><div><dt>Transcript effects</dt><dd>{snp.transcriptConsequences?.join(', ').replaceAll('_', ' ') || 'No transcript consequence assigned'}</dd></div><div><dt>Population MAF</dt><dd>{Number.isFinite(snp.maf) ? snp.maf.toFixed(4) : 'Not available'}</dd></div><div><dt>Study MAF</dt><dd>{Number.isFinite(snp.studyMaf) ? `${snp.studyMaf.toFixed(4)} · ${snp.studyMafClass || 'unclassified'}` : 'Not recorded'}</dd></div><div><dt>Regulatory</dt><dd>{snp.regulatoryConsequences?.join(', ').replaceAll('_', ' ') || 'No local regulatory consequence assigned'}</dd></div></dl></section>
        <section className="variant-card"><h2>GIM relationships</h2><p className="variant-card-note">Conditional genetic influence modules containing this independent SNP. Select a GIM to return to the complete heatmap and open its locus panel.</p><div className="chip-list">{selected.gimIds.length ? selected.gimIds.map((id) => <button className="gim-link-chip" key={id} onClick={() => onOpenGim(id)}>{formatGimId(id)}</button>) : <em>Not assigned</em>}</div><h3>Associated metabolites</h3><div className="metabolite-list">{[...selected.metabolites].sort().map((metabolite) => <span key={metabolite} title={metabolite}>{metabolite}</span>)}</div></section>
        <section className="variant-card"><h2>Functional prioritisation</h2><p className="variant-card-note">{nPredictedTracks ? <><b>{alphaSummary?.nModalities || alphaSummary?.modalities?.length || 0}</b> modalities across <b>{number.format(nPredictedTracks)}</b> retained AlphaGenome predictions: {gastricTissueScoreCount} stomach-tissue, {gastricCancerScoreCount} gastric-cancer, and {immuneScoreCount} immune-cell tracks.</> : alphaCoverage?.completed ? 'No target-scope AlphaGenome prediction record is available for this variant.' : 'The required 1 MB AlphaGenome scoring batch is pending or still running locally.'}</p>{alphaSummary?.scopes?.target?.splicingCombined?.score ? <p className="variant-card-note"><b>Official merged splicing:</b> {alphaSummary.scopes.target.splicingCombined.score.toPrecision(4)}{alphaSummary.scopes.target.splicingCombined.geneName ? ` · ${alphaSummary.scopes.target.splicingCombined.geneName}` : ''}</p> : null}{alphaSummary?.scopes?.target?.modalityStats?.length ? <div className="modality-evidence-list">{alphaSummary.scopes.target.modalityStats.map((item) => <span key={item.outputType} title={item.topTrack || item.topBiosample || item.outputType}><b>{outputLabel(item.outputType)}</b><i>{item.nTracks} tracks</i><em>{Number.isFinite(item.maxAbsRankingScore) ? `|Q| ${item.maxAbsRankingScore.toFixed(4)}${item.rankingScoreAtMaximum < 0 ? ' ↓' : ' ↑'}` : '|Q| —'}</em></span>)}</div> : null}<button className="functional-prediction-button" disabled={!nPredictedTracks} onClick={() => onOpenPrediction(snp.snpId)}>Inspect multimodal prediction <span>{nPredictedTracks ? '↗' : 'pending'}</span></button><h3>TF binding / motif</h3><p className="variant-card-note">{snp.tfMotif?.status === 'not_loaded' ? 'Inspect ENCODE TF ChIP-seq and JASPAR motif tracks at the UCSC locus; no unverified TF or motif call is shown as a local annotation.' : 'No TF or motif summary is loaded.'}</p></section>
      </div>
      <section className="phewas-section"><header><div><p className="eyebrow">HUMAN GENETIC ASSOCIATIONS</p><h2>PheWAS evidence</h2><p>On-demand public cohort associations for {canonical || snp.snpId}; each point encodes the reported β and P value.</p></div><div className="phewas-links"><a href={`https://pheweb.org/UKB-TOPMed/variant/${canonical || ''}`} target="_blank" rel="noreferrer">UKB-TOPMed ↗</a><a href={`https://r12.finngen.fi/variant/${canonical || ''}`} target="_blank" rel="noreferrer">FinnGen R12 ↗</a></div></header>{phewas.loading ? <p className="phewas-loading">Retrieving public PheWAS summaries…</p> : <div className="phewas-grid"><section><h3>UKB-TOPMed</h3><PhewasPlot source="UKB-TOPMed" result={phewas.ukb} /></section><section><h3>FinnGen R12</h3><PhewasPlot source="FinnGen R12" result={phewas.finngen} /></section></div>}<p className="phewas-caveat">Association results are cohort-specific and do not establish causality. β directions, phenotype definitions, ancestry composition, and sample sizes must be interpreted from the linked source records; displayed P values are not re-corrected for multiple testing.</p></section>
    </article>
  </section>
}

function formatSignal(value) {
  if (!Number.isFinite(value)) return '—'
  const absolute = Math.abs(value)
  if (absolute >= 100) return value.toFixed(1)
  if (absolute >= 10) return value.toFixed(2)
  if (absolute >= 1) return value.toFixed(3)
  return value.toPrecision(3)
}

function signalPath(values, width, height, maxValue, left = 0, top = 0) {
  if (!values?.length) return ''
  return values.map((value, index) => `${index ? 'L' : 'M'}${left + (index / Math.max(values.length - 1, 1)) * width},${top + height - (Number(value) / maxValue) * height}`).join(' ')
}

function PredictedSignalTrack({ track, interval, scopeLabel }) {
  const [hoveredBin, setHoveredBin] = useState(null)
  const maximum = Math.max(...track.ref, ...track.alt, Number.MIN_VALUE)
  const svgWidth = 620
  const svgHeight = 86
  const plotLeft = 45
  const plotRight = 8
  const plotTop = 7
  const plotHeight = 58
  const plotWidth = svgWidth - plotLeft - plotRight
  const variantX = plotLeft + (interval.variantOffset / interval.width) * plotWidth
  const selectedBin = hoveredBin === null ? null : Math.min(Math.max(hoveredBin, 0), track.ref.length - 1)
  const relativeKb = selectedBin === null ? null : ((selectedBin * interval.binSize + interval.binSize / 2 - interval.variantOffset) / 1000)
  const hoverX = selectedBin === null ? null : plotLeft + (selectedBin / Math.max(track.ref.length - 1, 1)) * plotWidth
  const readPointer = (event) => {
    const rect = event.currentTarget.getBoundingClientRect()
    const svgX = ((event.clientX - rect.left) / rect.width) * svgWidth
    const fraction = Math.min(Math.max((svgX - plotLeft) / plotWidth, 0), 1)
    setHoveredBin(Math.round(fraction * Math.max(track.ref.length - 1, 1)))
  }
  return <article className="signal-track"><div><b>{outputLabel(track.modality)}</b><span title={track.name}>{track.name}</span></div><svg viewBox={`0 0 ${svgWidth} ${svgHeight}`} preserveAspectRatio="none" role="img" aria-label={`${track.name}, predicted reference and alternate signal`} onMouseMove={readPointer} onMouseLeave={() => setHoveredBin(null)}><text className="signal-y-title" x="9" y={svgHeight / 2} transform={`rotate(-90 9 ${svgHeight / 2})`}>signal</text><text className="signal-y-tick" x={plotLeft - 5} y={plotTop + 4} textAnchor="end">{formatSignal(maximum)}</text><text className="signal-y-tick" x={plotLeft - 5} y={plotTop + plotHeight} textAnchor="end">0</text><line className="signal-y-axis" x1={plotLeft} x2={plotLeft} y1={plotTop} y2={plotTop + plotHeight}></line><line className="signal-baseline" x1={plotLeft} x2={plotLeft + plotWidth} y1={plotTop + plotHeight} y2={plotTop + plotHeight}></line><line className="signal-variant-line" x1={variantX} x2={variantX} y1={plotTop} y2={plotTop + plotHeight}></line>{hoverX !== null && <line className="signal-hover-line" x1={hoverX} x2={hoverX} y1={plotTop} y2={plotTop + plotHeight}></line>}<path className="signal-ref-path" d={signalPath(track.ref, plotWidth, plotHeight, maximum, plotLeft, plotTop)}></path><path className="signal-alt-path" d={signalPath(track.alt, plotWidth, plotHeight, maximum, plotLeft, plotTop)}></path><text className="signal-x-label" x={plotLeft} y={svgHeight - 4}>−8 kb</text><text className="signal-x-label" x={plotLeft + plotWidth} y={svgHeight - 4} textAnchor="end">+8 kb</text></svg><small>{track.biosample || scopeLabel}{track.biosampleType ? ` · ${track.biosampleType.replaceAll('_', ' ')}` : ''}{selectedBin !== null ? ` · ${relativeKb >= 0 ? '+' : ''}${relativeKb.toFixed(2)} kb · REF ${formatSignal(track.ref[selectedBin])} · ALT ${formatSignal(track.alt[selectedBin])}` : ''}</small></article>
}

function AlphaGenomeSignalTracks({ snpId, scope }) {
  const [state, setState] = useState({ loading: true, profile: null, error: null })
  useEffect(() => {
    if (scope === 'target' || scope === 'gastric_cancer') {
      setState({ loading: false, profile: null, error: null })
      return undefined
    }
    const cached = alphaSignalCache.get(snpId)
    if (cached) {
      setState({ loading: false, profile: cached, error: null })
      return undefined
    }
    let cancelled = false
    setState({ loading: true, profile: null, error: null })
    fetch(`/data/alphagenome_signals/${encodeURIComponent(snpId)}.json`)
      .then((response) => {
        if (!response.ok) throw new Error('Signal profile pending')
        return response.json()
      })
      .then((profile) => {
        if (cancelled) return
        alphaSignalCache.set(snpId, profile)
        setState({ loading: false, profile, error: null })
      })
      .catch(() => !cancelled && setState({ loading: false, profile: null, error: 'Signal profile is being generated for this SNP.' }))
    return () => { cancelled = true }
  }, [snpId, scope])

  if (scope === 'target') return <p className="signal-track-note">Choose Stomach tissue or Immune cells to inspect representative central REF/ALT signal profiles. Scalar prioritisation above uses the full 1 MB context.</p>
  if (scope === 'gastric_cancer') return <p className="signal-track-note">Gastric-cancer evidence is retained in the 1 MB scalar scores when AlphaGenome provides a matching cell-line track. No unrelated cell line is substituted for the local signal plot.</p>
  if (state.loading) return <p className="signal-track-note">Loading actual REF/ALT local signal profiles…</p>
  const scopeRows = scope === 'immune'
    ? [state.profile?.scopes?.tcell, state.profile?.scopes?.bcell].filter(Boolean)
    : [state.profile?.scopes?.gastric_tissue].filter(Boolean)
  const tracks = scopeRows.flatMap((item) => item.tracks || [])
  if (state.error || !tracks.length) return <p className="signal-track-note">{state.error || 'No profile was returned for this biological entity.'}</p>
  const { interval, scopes } = state.profile
  const label = scope === 'immune' ? 'Representative primary T- and B-cell profiles' : scopes.gastric_tissue?.label || 'Stomach tissue'
  return <section className="signal-track-section"><header><div><p className="eyebrow">CENTRAL PREDICTED SIGNAL</p><h3>{label}</h3></div><p>16 kb display window · 128 bp mean bins · <i className="signal-ref"></i>REF <i className="signal-alt"></i>ALT · centre line = variant</p></header><div className="signal-track-list">{tracks.map((track, index) => <PredictedSignalTrack key={`${track.modality}-${track.name}-${index}`} track={track} interval={interval} scopeLabel={label} />)}</div></section>
}

function VariantFunctionalPredictionDialog({ cell, snp, scores, loading, error, onClose }) {
  const [scope, setScope] = useState('gastric_tissue')
  const gastricTissueScores = useMemo(() => scores.filter((score) => scoreIsInScope(score, 'gastric_tissue')), [scores])
  const gastricCancerScores = useMemo(() => scores.filter((score) => scoreIsInScope(score, 'gastric_cancer')), [scores])
  const immuneScores = useMemo(() => scores.filter((score) => scoreIsInScope(score, 'immune')), [scores])
  const displayedScores = scope === 'gastric_tissue' ? gastricTissueScores : scope === 'gastric_cancer' ? gastricCancerScores : scope === 'immune' ? immuneScores : scores
  const modalityGroups = useMemo(() => {
    const groups = new Map()
    displayedScores.forEach((score) => {
      const groupKey = `${score.outputType}::${score.variantScorer || 'legacy'}`
      const rows = groups.get(groupKey) || []
      rows.push(score)
      groups.set(groupKey, rows)
    })
    return [...groups.entries()]
      .map(([groupKey, rows]) => [groupKey, rows.sort((left, right) => Math.abs(predictionRankingValue(right) ?? predictionScoreValue(right)) - Math.abs(predictionRankingValue(left) ?? predictionScoreValue(left)))])
      .sort((left, right) => Math.abs(predictionRankingValue(right[1][0]) ?? predictionScoreValue(right[1][0]) ?? 0) - Math.abs(predictionRankingValue(left[1][0]) ?? predictionScoreValue(left[1][0]) ?? 0))
  }, [displayedScores])
  const modalityCount = new Set(displayedScores.map((score) => score.outputType)).size
  const signedScores = displayedScores.filter((score) => scorerDefinition(score).signed)
  const unsignedScores = displayedScores.filter((score) => !scorerDefinition(score).signed)
  const mergedSplicing = useMemo(() => combinedSplicingScore(displayedScores), [displayedScores])
  const studyAlt = snp?.alphaGenome?.input?.alternates?.join(', ') || snp?.grch38?.studyAlternates?.join(', ') || 'ALT'

  return <div className="prediction-backdrop" role="presentation" onMouseDown={onClose}>
    <section className="prediction-dialog" role="dialog" aria-modal="true" aria-label="Variant functional prediction" onMouseDown={(event) => event.stopPropagation()}>
      <header className="prediction-header">
        <div>
          <p className="eyebrow">ALPHAGENOME VARIANT EFFECT</p>
          <h2>Variant functional prediction</h2>
          <p>{cell.snpId} · {formatCoordinate(snp)} · {snp?.grch38?.reference || 'REF'} → {studyAlt}</p>
          <div className="prediction-scope" aria-label="Prediction entity scope">
            <button className={scope === 'gastric_tissue' ? 'active' : ''} onClick={() => setScope('gastric_tissue')}>Stomach tissue <b>{gastricTissueScores.length}</b></button>
            <button className={scope === 'gastric_cancer' ? 'active' : ''} onClick={() => setScope('gastric_cancer')}>Gastric cancer <b>{gastricCancerScores.length}</b></button>
            <button className={scope === 'immune' ? 'active' : ''} onClick={() => setScope('immune')}>Immune cells <b>{immuneScores.length}</b></button>
            <button className={scope === 'target' ? 'active' : ''} onClick={() => setScope('target')}>All target scopes <b>{scores.length}</b></button>
          </div>
        </div>
        <button className="close-button" onClick={onClose} aria-label="Close variant functional prediction">×</button>
      </header>
      <div className="prediction-summary">
        <div><span>FUNCTIONAL MODALITIES</span><b>{modalityCount}</b></div>
        <div><span>PREDICTED TRACKS</span><b>{displayedScores.length}</b></div>
        <div><span>SCORER DEFINITIONS</span><b>{modalityGroups.length}</b></div>
        <div><span>MERGED SPLICING</span><b>{mergedSplicing?.score ? mergedSplicing.score.toPrecision(3) : '—'}</b></div>
        <p>{signedScores.length ? <><i className="prediction-positive"></i> signed: ALT higher <i className="prediction-negative"></i> ALT lower</> : null}{signedScores.length && unsignedScores.length ? <br /> : null}{unsignedScores.length ? <><i className="prediction-magnitude"></i> unsigned: disruption or activity magnitude</> : null}</p>
      </div>
      <div className="prediction-body">
        <aside className="prediction-explainer">
          <b>What this score means</b>
          <p>Scalar effects use AlphaGenome’s recommended 1 MB context and compare the study ALT with GRCh38 REF. Raw values retain scorer-specific units; |Q| ranks unusual effects within the same scorer and track, while the sign is kept separately. The merged splicing value is max |splice sites| + max |splice-site usage| + max |splice junctions| / 5. None of these values is an association P value.</p>
        </aside>
        {!loading && !error && <AlphaGenomeSignalTracks snpId={cell.snpId} scope={scope} />}
        {loading ? <section className="prediction-empty"><b>Loading the detailed AlphaGenome tracks for {cell.snpId}…</b><p>The GIM heatmap and Variant Explorer use a compact prediction index; only this selected SNP’s complete track set is downloaded.</p></section> : error ? <section className="prediction-empty"><b>{error}</b><p>Try again after the local AlphaGenome batch has published this SNP’s detailed score file.</p></section> : displayedScores.length ? modalityGroups.map(([groupKey, rows]) => {
          const exemplar = rows[0]
          const metric = scorerDefinition(exemplar)
          const groupStrongestScore = Math.max(...rows.map((score) => Math.abs(predictionScoreValue(score))), 0.001)
          return <section className="prediction-modality" key={groupKey}>
          <header><div><h3>{outputLabel(exemplar.outputType)}</h3><p>{metric.label}</p></div><span>{rows.length} tracks</span></header>
          <p className="prediction-formula">{metric.formula}{metric.legacy ? ' — this card will be replaced as the complete rescoring batch publishes.' : ''}</p>
          {rows.map((score, index) => {
            const value = predictionScoreValue(score)
            const percentage = Math.min(Math.abs(value) / groupStrongestScore, 1) * (metric.signed ? 50 : 100)
            const left = metric.signed ? (value >= 0 ? 50 : 50 - percentage) : 0
            const ranking = predictionRankingValue(score)
            const trackLabel = [score.track, score.biosample, score.geneName || score.geneId, score.transcriptionFactor, score.histoneMark].filter(Boolean).join(' · ') || 'Predicted track'
            return <div className="prediction-track" key={`${score.rankInOutput || index}-${trackLabel}`}>
              <div><b title={trackLabel}>{trackLabel}</b><span className={metric.signed ? (value >= 0 ? 'positive' : 'negative') : 'magnitude'}>{metric.signed && value >= 0 ? '+' : ''}{value.toFixed(3)}{Number.isFinite(ranking) ? ` · Q ${ranking >= 0 ? '+' : ''}${ranking.toFixed(4)}` : ''}</span></div>
              <div className={`prediction-bar ${metric.signed ? '' : 'unsigned'}`} aria-label={`${outputLabel(exemplar.outputType)} raw score ${value.toFixed(3)}`}><i className={metric.signed ? (value >= 0 ? 'up' : 'down') : 'magnitude'} style={{ left: `${left}%`, width: `${percentage}%` }}></i></div>
            </div>
          })}
        </section>}) : <section className="prediction-empty"><b>No AlphaGenome track is available for this target scope.</b><p>The portal retains only stomach tissue, gastric-cancer and immune-cell evidence, plus the tissue-independent splice-site component required for merged splicing. It does not substitute unrelated tissues.</p></section>}
      </div>
    </section>
  </div>
}

export default App
