#!/usr/bin/env python3
"""Normalize a portable GIMForge bundle into browser-ready portal data.

The bundle is deliberately self-contained: a server can export it without
internet access, an internet-connected workstation can add Ensembl and
AlphaGenome annotations, and either machine can then build the static portal.
No individual-level genotype, phenotype, covariate, or API credential is stored
in the bundle.
"""

from __future__ import annotations

import csv
import gzip
import json
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


ENSEMBL_URL = "https://rest.ensembl.org"
UCSC_CYTOBAND_URL = (
    "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=cytoBand"
)


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def number(value: str | None, default: float | None = None):
    if value in (None, "", "NA", "NaN"):
        return default
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return int(parsed) if parsed.is_integer() else parsed


def value(row: dict, *names: str, default=None):
    """Return the first populated alias from a GIMForge or legacy result row."""

    for name in names:
        candidate = row.get(name)
        if candidate not in (None, "", "NA", "NaN"):
            return candidate
    return default


def truthy(candidate) -> bool:
    return candidate is True or str(candidate).strip().lower() in {"true", "1", "yes"}


def split_values(value: str | None) -> list[str]:
    if not value:
        return []
    return [part for part in value.split(";") if part]


def optional_table(path: Path) -> list[dict[str, str]]:
    return read_tsv(path) if path.is_file() else []


def table_path(source_dir: Path, stem: str, *, required: bool = False) -> Path | None:
    for suffix in (".tsv.gz", ".tsv"):
        candidate = source_dir / f"{stem}{suffix}"
        if candidate.is_file():
            return candidate
    if required:
        raise FileNotFoundError(
            f"Portable bundle is missing source/{stem}.tsv(.gz). "
            "Run the prepare stage again from a completed GIMForge result directory."
        )
    return None


def read_bundle_manifest(bundle_dir: Path) -> dict:
    path = bundle_dir / "bundle_manifest.json"
    if not path.is_file():
        raise FileNotFoundError(f"Not a GIMs Portal bundle: {path} is missing.")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schemaVersion") != 1:
        raise ValueError(f"Unsupported portal bundle schema: {payload.get('schemaVersion')!r}")
    return payload


def rest_post(path: str, payload: dict, retries: int = 5) -> dict:
    body = json.dumps(payload).encode("utf-8")
    request = Request(
        f"{ENSEMBL_URL}{path}",
        data=body,
        headers={"Accept": "application/json", "Content-Type": "application/json", "User-Agent": "GIMs-portal/1.0"},
        method="POST",
    )
    for attempt in range(retries):
        try:
            with urlopen(request, timeout=90) as response:
                return json.loads(response.read().decode("utf-8"))
        except HTTPError as error:
            retry_after = error.headers.get("Retry-After")
            if error.code not in {429, 500, 502, 503, 504} or attempt == retries - 1:
                raise RuntimeError(f"Ensembl {path}: HTTP {error.code}") from error
            time.sleep(float(retry_after) if retry_after else min(2**attempt, 20))
        except (URLError, TimeoutError, json.JSONDecodeError) as error:
            if attempt == retries - 1:
                raise RuntimeError(f"Ensembl {path}: {error}") from error
            time.sleep(min(2**attempt, 20))
    raise AssertionError("unreachable")


def rest_get(url: str, retries: int = 5) -> dict:
    request = Request(
        url,
        headers={"Accept": "application/json", "User-Agent": "GIMs-portal/1.0"},
    )
    for attempt in range(retries):
        try:
            with urlopen(request, timeout=90) as response:
                return json.loads(response.read().decode("utf-8"))
        except HTTPError as error:
            retry_after = error.headers.get("Retry-After")
            if error.code not in {429, 500, 502, 503, 504} or attempt == retries - 1:
                raise RuntimeError(f"Annotation request {url}: HTTP {error.code}") from error
            time.sleep(float(retry_after) if retry_after else min(2**attempt, 20))
        except (URLError, TimeoutError, json.JSONDecodeError) as error:
            if attempt == retries - 1:
                raise RuntimeError(f"Annotation request {url}: {error}") from error
            time.sleep(min(2**attempt, 20))
    raise AssertionError("unreachable")


def batches(values: list[str], size: int = 50):
    for index in range(0, len(values), size):
        yield values[index : index + size]


def fetch_ensembl_annotations(snp_ids: list[str], cache_path: Path) -> dict[str, dict]:
    cached: dict[str, dict] = {}
    if cache_path.exists():
        cached = json.loads(cache_path.read_text(encoding="utf-8"))
        if set(snp_ids).issubset(cached):
            print(f"Using cached Ensembl annotations for {len(snp_ids):,} variants.")
            return cached

    missing_ids = [snp_id for snp_id in snp_ids if snp_id not in cached]

    def fetch_batch(ids: list[str]) -> dict[str, dict]:
        variation = rest_post("/variation/homo_sapiens", {"ids": ids})
        vep_rows = rest_post(
            "/vep/homo_sapiens/id",
            {"ids": ids, "pick": 1, "canonical": 1, "mane": 1, "regulatory": 1},
        )
        vep: dict[str, dict] = {}
        for row in vep_rows:
            if row.get("id"):
                vep[row["id"]] = row
        return {snp_id: {"variation": variation.get(snp_id, {}), "vep": vep.get(snp_id, {})} for snp_id in ids}

    batch_list = list(batches(missing_ids))
    print(f"Resolving {len(missing_ids):,} uncached variants in {len(batch_list)} Ensembl batches.", flush=True)
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(fetch_batch, ids): index for index, ids in enumerate(batch_list, start=1)}
        for future in as_completed(futures):
            batch_number = futures[future]
            cached.update(future.result())
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(json.dumps(cached, ensure_ascii=False), encoding="utf-8")
            print(f"Completed Ensembl batch {batch_number}/{len(batch_list)}; cache now contains {len(cached):,} variants.", flush=True)

    return cached


def best_hg38_mapping(variation: dict) -> dict:
    mappings = variation.get("mappings") or []
    return next((mapping for mapping in mappings if mapping.get("assembly_name") == "GRCh38"), {})


def parse_alleles(allele_string: str | None) -> tuple[str | None, list[str]]:
    values = (allele_string or "").split("/")
    return (values[0] if values else None, values[1:] if len(values) > 1 else [])


def _complement_allele(allele: str) -> str | None:
    normalized = str(allele or "").upper()
    if not normalized or any(base not in "ACGT" for base in normalized):
        return None
    return normalized.translate(str.maketrans("ACGT", "TGCA"))


def select_study_alternates(
    reference: str | None,
    ensembl_alternates: list[str],
    observed_alleles: list[str],
) -> tuple[list[str], str | None]:
    """Match the two study alleles to the GRCh38 REF orientation.

    PLINK BIM alleles are not guaranteed to use the same strand as Ensembl.
    Only alleles actually observed in the study are returned. This prevents a
    multiallelic rsID from being scored for unrelated Ensembl ALT alleles.
    """

    normalized_reference = str(reference or "").upper()
    normalized_alternates = [
        str(allele).upper() for allele in ensembl_alternates if allele
    ]
    observed = {str(allele).upper() for allele in observed_alleles if allele}
    if not normalized_reference or not normalized_alternates or not observed:
        return [], None

    direct = [
        allele for allele in normalized_alternates if allele in observed
    ]
    if normalized_reference in observed and direct:
        return direct, "forward"

    complemented = {
        value
        for allele in observed
        if (value := _complement_allele(allele)) is not None
    }
    reverse = [
        allele for allele in normalized_alternates if allele in complemented
    ]
    if normalized_reference in complemented and reverse:
        return reverse, "reverse_complement"

    # A rare biallelic study may contain two non-reference alleles at a
    # multiallelic rsID. Both are study-observed and remain valid REF-vs-ALT
    # model comparisons, even though the reference allele was not genotyped.
    if (
        normalized_reference not in observed
        and normalized_reference not in complemented
        and direct
    ):
        return direct, "forward_nonreference_pair"
    if (
        normalized_reference not in observed
        and normalized_reference not in complemented
        and reverse
    ):
        return reverse, "reverse_complement_nonreference_pair"
    return [], None


def clean_gene_consequences(vep: dict) -> tuple[list[str], list[str], list[str], list[str], list[str]]:
    transcript_rows = vep.get("transcript_consequences") or []
    gene_ids = sorted({row["gene_id"] for row in transcript_rows if row.get("gene_id")})
    gene_symbols = sorted({row["gene_symbol"] for row in transcript_rows if row.get("gene_symbol")})
    transcript_ids = sorted({row["transcript_id"] for row in transcript_rows if row.get("transcript_id")})
    consequences = sorted({term for row in transcript_rows for term in row.get("consequence_terms", [])})
    regulatory = sorted(
        {term for row in (vep.get("regulatory_feature_consequences") or []) for term in row.get("consequence_terms", [])}
    )
    return gene_ids, gene_symbols, transcript_ids, consequences, regulatory


def genomic_link(chromosome: str | None, position: int | None) -> str | None:
    if not chromosome or not position:
        return None
    return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr{chromosome}%3A{position - 250}%2D{position + 250}"


def bundle_snp_ids(bundle_dir: Path) -> list[str]:
    source_dir = bundle_dir / "source"
    members_path = table_path(source_dir, "members", required=True)
    members = read_tsv(members_path)
    snp_ids = {
        str(value(row, "node_id"))
        for row in members
        if str(value(row, "node_type", default="")).upper() == "SNP"
        and value(row, "node_id")
    }
    signals_path = table_path(source_dir, "independent_signals")
    if signals_path:
        snp_ids.update(
            str(snp_id)
            for row in read_tsv(signals_path)
            if (snp_id := value(row, "snp_id", "ID", "markername"))
        )
    return sorted(snp_ids)


def annotate_ensembl(bundle_dir: Path) -> dict[str, dict]:
    """Resolve all bundle SNPs and checkpoint the portable annotation cache."""

    read_bundle_manifest(bundle_dir)
    annotation_dir = bundle_dir / "annotations"
    annotation_dir.mkdir(parents=True, exist_ok=True)
    return fetch_ensembl_annotations(
        bundle_snp_ids(bundle_dir),
        annotation_dir / "ensembl_annotation_cache.json",
    )


def annotate_cytobands(bundle_dir: Path) -> dict:
    """Cache canonical GRCh38 cytobands for reproducible offline portal builds."""

    read_bundle_manifest(bundle_dir)
    output = bundle_dir / "annotations" / "cytobands_grch38.json"
    if output.is_file():
        payload = json.loads(output.read_text(encoding="utf-8"))
        if payload.get("assembly") == "GRCh38" and payload.get("chromosomes"):
            print("Using cached UCSC GRCh38 cytobands.")
            return payload

    raw = rest_get(UCSC_CYTOBAND_URL)
    raw_chromosomes = raw.get("cytoBand") or {}
    chromosomes: dict[str, list[dict]] = {}
    for chromosome in [str(value) for value in range(1, 23)] + ["X", "Y"]:
        rows = raw_chromosomes.get(f"chr{chromosome}") or []
        chromosomes[chromosome] = [
            {
                "start": number(row.get("chromStart")),
                "end": number(row.get("chromEnd")),
                "band": str(row.get("name") or ""),
                "stain": str(row.get("gieStain") or ""),
            }
            for row in rows
            if row.get("name")
        ]
    payload = {
        "schemaVersion": 1,
        "assembly": "GRCh38",
        "source": UCSC_CYTOBAND_URL,
        "chromosomes": chromosomes,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    print(
        "Cached "
        f"{sum(len(rows) for rows in chromosomes.values()):,} GRCh38 cytobands."
    )
    return payload


def cytoband_at(
    cytobands: dict, chromosome: str | None, position: int | float | None
) -> str | None:
    if not chromosome or position is None:
        return None
    canonical = str(chromosome).removeprefix("chr")
    for row in (cytobands.get("chromosomes") or {}).get(canonical, []):
        start = number(row.get("start"))
        end = number(row.get("end"))
        if start is not None and end is not None and start <= position < end:
            return f"{canonical}{row.get('band')}"
    return None


def write_alphagenome_input(bundle_dir: Path) -> Path:
    """Create the GRCh38 REF/ALT manifest consumed by AlphaGenome."""

    manifest = read_bundle_manifest(bundle_dir)
    annotation_dir = bundle_dir / "annotations"
    source_variants = _source_variants(
        bundle_dir / "source",
        str(manifest.get("sourceGenomeBuild") or "GRCh37"),
    )
    cache_path = annotation_dir / "ensembl_annotation_cache.json"
    annotations = (
        json.loads(cache_path.read_text(encoding="utf-8"))
        if cache_path.is_file()
        else {}
    )
    rows = []
    for snp_id in bundle_snp_ids(bundle_dir):
        variation = (annotations.get(snp_id) or {}).get("variation") or {}
        mapping = best_hg38_mapping(variation)
        reference, ensembl_alternates = parse_alleles(mapping.get("allele_string"))
        observed_alleles = source_variants.get(snp_id, {}).get(
            "observedAlleles", []
        )
        alternates, orientation = select_study_alternates(
            reference, ensembl_alternates, observed_alleles
        )
        if not mapping or not reference or not alternates:
            continue
        chromosome = mapping.get("seq_region_name")
        position = number(mapping.get("start"))
        if not chromosome or not position:
            continue
        rows.append(
            {
                "snp_id": snp_id,
                "chromosome": f"chr{chromosome}",
                "position": position,
                "reference_bases": reference,
                "alternate_bases": ",".join(alternates),
                "study_observed_alleles": ",".join(observed_alleles),
                "allele_match_orientation": orientation,
                "assembly": "GRCh38",
            }
        )
    output = annotation_dir / "alphagenome_input.tsv"
    annotation_dir.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "snp_id",
                "chromosome",
                "position",
                "reference_bases",
                "alternate_bases",
                "study_observed_alleles",
                "allele_match_orientation",
                "assembly",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    return output


def _source_variants(source_dir: Path, source_build: str) -> dict[str, dict]:
    path = source_dir / "variants.tsv"
    if not path.is_file():
        return {}
    result = {}
    for row in read_tsv(path):
        snp_id = value(row, "snp_id")
        if not snp_id:
            continue
        result[str(snp_id)] = {
            "chromosome": value(row, "chromosome"),
            "grch37Position": number(value(row, "position")),
            "observedAlleles": [
                allele
                for allele in (value(row, "allele1"), value(row, "allele2"))
                if allele
            ],
            "assembly": source_build,
        }
    return result


def _json_write(path: Path, payload, *, indent: int | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=indent),
        encoding="utf-8",
    )


def build_data(bundle_dir: Path, out_dir: Path) -> dict:
    """Build the portal's core JSON/CSV files using only portable bundle data."""

    manifest = read_bundle_manifest(bundle_dir)
    source_dir = bundle_dir / "source"
    annotation_dir = bundle_dir / "annotations"
    matrix_path = table_path(source_dir, "matrix_out", required=True)
    edges_path = table_path(source_dir, "edges", required=True)
    members_path = table_path(source_dir, "members", required=True)
    matrix_rows = read_tsv(matrix_path)
    edges = read_tsv(edges_path)
    members = read_tsv(members_path)
    signals_path = table_path(source_dir, "independent_signals")
    signals = read_tsv(signals_path) if signals_path else edges
    regions_path = table_path(source_dir, "regions")
    region_rows = read_tsv(regions_path) if regions_path else []
    region_summary_path = table_path(source_dir, "region_summary")
    region_summary_rows = (
        read_tsv(region_summary_path) if region_summary_path else []
    )
    region_metabolites_path = table_path(source_dir, "region_metabolites")
    region_metabolite_rows = (
        read_tsv(region_metabolites_path) if region_metabolites_path else []
    )
    differential_path = source_dir / "differential_metabolomics.csv"
    differential_rows = read_csv(differential_path) if differential_path.is_file() else []

    if matrix_rows and not any(
        number(value(row, "p", "P", "p_value")) is not None
        for row in matrix_rows
    ):
        raise ValueError(
            "matrix_out contains no usable P values. This is usually an old "
            "components-only preview; re-run `gimforge components` with the "
            "latest GIMForge or point prepare at the original matrix_out."
        )

    source_build = str(manifest.get("sourceGenomeBuild") or "GRCh37")
    source_variants = _source_variants(source_dir, source_build)
    snp_ids = bundle_snp_ids(bundle_dir)
    cache_path = annotation_dir / "ensembl_annotation_cache.json"
    annotations = (
        json.loads(cache_path.read_text(encoding="utf-8"))
        if cache_path.is_file()
        else {}
    )
    cytoband_path = annotation_dir / "cytobands_grch38.json"
    cytobands = (
        json.loads(cytoband_path.read_text(encoding="utf-8"))
        if cytoband_path.is_file()
        else {}
    )
    cached_ids = set(annotations)
    annotation_status = (
        "resolved"
        if set(snp_ids).issubset(cached_ids)
        else "partial"
        if cached_ids
        else "not_fetched"
    )

    member_regions: dict[str, set[str]] = defaultdict(set)
    study_frequency: dict[str, tuple[float | None, str | None]] = {}
    for row in members:
        if str(value(row, "node_type", default="")).upper() != "SNP":
            continue
        snp_id = str(value(row, "node_id"))
        region_id = value(row, "region_id")
        if region_id:
            member_regions[snp_id].add(str(region_id))
        maf = number(value(row, "maf"))
        maf_class = value(row, "maf_class")
        previous = study_frequency.get(snp_id)
        if previous is None or previous[0] is None:
            study_frequency[snp_id] = (maf, str(maf_class) if maf_class else None)

    signal_by_snp: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in signals:
        snp_id = value(row, "snp_id", "ID", "markername")
        if snp_id:
            signal_by_snp[str(snp_id)].append(row)

    snp_records: list[dict] = []
    for snp_id in snp_ids:
        source = source_variants.get(snp_id, {})
        annotation = annotations.get(snp_id, {})
        variation = annotation.get("variation") or {}
        vep = annotation.get("vep") or {}
        mapping = best_hg38_mapping(variation)
        reference, alternates = parse_alleles(mapping.get("allele_string"))
        study_alternates, allele_orientation = select_study_alternates(
            reference,
            alternates,
            source.get("observedAlleles", []),
        )
        (
            gene_ids,
            gene_symbols,
            transcript_ids,
            transcript_terms,
            regulatory_terms,
        ) = clean_gene_consequences(vep)
        linked_rows = signal_by_snp[snp_id]
        p_values = [
            parsed
            for row in linked_rows
            if (parsed := number(value(row, "p", "P", "p_full", "p_value"))) is not None
        ]
        betas = [
            parsed
            for row in linked_rows
            if (parsed := number(value(row, "beta", "BETA", "beta_full", "effect"))) is not None
        ]
        chromosome = mapping.get("seq_region_name") or source.get("chromosome")
        hg38_position = number(mapping.get("start"))
        study_maf, study_maf_class = study_frequency.get(snp_id, (None, None))
        region_ids = sorted(
            member_regions.get(snp_id, set())
            | {
                str(region)
                for row in linked_rows
                if (region := value(row, "region_id", "region"))
            }
        )
        metabolites = sorted(
            {
                str(metabolite)
                for row in linked_rows
                if (metabolite := value(row, "metabolite", "phenotype", "trait"))
            }
        )
        snp_records.append(
            {
                "snpId": snp_id,
                "grch37": source if source_build.upper().startswith("GRCH37") else None,
                "sourceCoordinate": source or None,
                "grch38": {
                    "chromosome": chromosome,
                    "position": hg38_position,
                    "reference": reference,
                    "alternates": alternates,
                    "studyAlternates": study_alternates,
                    "alleleMatchOrientation": allele_orientation,
                    "assembly": mapping.get("assembly_name"),
                }
                if mapping
                else None,
                "variantClass": variation.get("var_class"),
                "mostSevereConsequence": variation.get("most_severe_consequence")
                or vep.get("most_severe_consequence"),
                "maf": number(variation.get("MAF")),
                "minorAllele": variation.get("minor_allele"),
                "studyMaf": study_maf,
                "studyMafClass": study_maf_class,
                "geneIds": gene_ids,
                "geneSymbols": gene_symbols,
                "transcriptIds": transcript_ids,
                "transcriptConsequences": transcript_terms,
                "regulatoryConsequences": regulatory_terms,
                "regionIds": region_ids,
                "metabolites": metabolites,
                "nAssociations": len(linked_rows),
                "minP": min(p_values) if p_values else None,
                "maxAbsBeta": max((abs(beta) for beta in betas), default=None),
                "ucscUrl": genomic_link(chromosome, hg38_position),
                "ensemblUrl": (
                    f"https://www.ensembl.org/Homo_sapiens/Variation/Explore?v={snp_id}"
                ),
                "alphaGenome": {
                    "status": "awaiting_api_score",
                    "input": {
                        "chromosome": f"chr{chromosome}" if chromosome else None,
                        "position": hg38_position,
                        "reference": reference,
                        "alternates": study_alternates,
                        "observedAlleles": source.get("observedAlleles", []),
                        "alleleMatchOrientation": allele_orientation,
                        "assembly": "GRCh38",
                    },
                },
                "tfMotif": {
                    "status": "not_loaded",
                    "reason": (
                        "Use the linked UCSC locus to inspect ENCODE TF ChIP-seq "
                        "and JASPAR motif tracks."
                    ),
                },
            }
        )

    annotation_by_snp = {row["snpId"]: row for row in snp_records}
    records: list[dict] = []
    for row in signals:
        snp_id = value(row, "snp_id", "ID", "markername")
        metabolite = value(row, "metabolite", "phenotype", "trait")
        region_id = value(row, "region_id", "region")
        if not snp_id or not metabolite or not region_id:
            continue
        snp = annotation_by_snp.get(str(snp_id), {})
        records.append(
            {
                "regionId": str(region_id),
                "snpId": str(snp_id),
                "metabolite": str(metabolite),
                "phenotypeId": value(row, "phenotype_id"),
                "assayGroup": value(row, "assay_group"),
                "beta": number(value(row, "beta", "BETA", "beta_full", "effect")),
                "se": number(value(row, "se", "SE", "se_full", "standard_error")),
                "p": number(value(row, "p", "P", "p_full", "p_value")),
                "n": number(value(row, "n", "N", "n_full", "OBS_CT")),
                "consequence": snp.get("mostSevereConsequence"),
                "geneIds": snp.get("geneIds", []),
                "geneSymbols": snp.get("geneSymbols", []),
                "transcriptIds": snp.get("transcriptIds", []),
            }
        )

    members_by_gim: dict[str, dict[str, list[dict[str, str]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for member in members:
        gim_id = value(member, "gim_id")
        node_type = str(value(member, "node_type", default=""))
        if gim_id and node_type:
            canonical_type = "SNP" if node_type.upper() == "SNP" else "metabolite"
            members_by_gim[str(gim_id)][canonical_type].append(member)
    matrix_by_key = {
        (
            str(value(row, "region_id", "region")),
            str(value(row, "snp_id", "ID", "markername")),
            str(value(row, "metabolite", "phenotype", "trait")),
        ): row
        for row in matrix_rows
    }
    edges_by_gim: dict[str, list[dict[str, str]]] = defaultdict(list)
    edges_by_region: dict[str, list[dict[str, str]]] = defaultdict(list)
    for edge in edges:
        gim_id = value(edge, "gim_id")
        region_id = value(edge, "region_id", "region")
        if gim_id:
            edges_by_gim[str(gim_id)].append(edge)
        if region_id:
            edges_by_region[str(region_id)].append(edge)

    gim_entities: list[dict] = []
    for gim_id, member_nodes in members_by_gim.items():
        all_nodes = member_nodes["SNP"] + member_nodes["metabolite"]
        if not all_nodes:
            continue
        region_id = str(value(all_nodes[0], "region_id", "region"))
        snp_nodes = sorted(
            member_nodes["SNP"],
            key=lambda row: (
                number(value(row, "marker_order", "order"), 99999),
                str(value(row, "node_id")),
            ),
        )
        snp_list = [str(value(row, "node_id")) for row in snp_nodes]
        metabolite_list = [
            str(value(row, "node_id")) for row in member_nodes["metabolite"]
        ]
        direct_pairs = {
            (
                str(value(edge, "snp_id", "ID", "markername")),
                str(value(edge, "metabolite", "phenotype", "trait")),
            )
            for edge in edges_by_gim[gim_id]
        }
        heatmap = []
        for snp_id in snp_list:
            for metabolite in metabolite_list:
                matrix_row = matrix_by_key.get((region_id, snp_id, metabolite), {})
                parsed_p = number(value(matrix_row, "p", "P", "p_value"))
                heatmap.append(
                    {
                        "snpId": snp_id,
                        "metabolite": metabolite,
                        "beta": number(value(matrix_row, "beta", "BETA", "effect")),
                        "p": parsed_p,
                        "testable": truthy(value(matrix_row, "testable"))
                        if "testable" in matrix_row
                        else parsed_p is not None,
                        "conditionedOn": number(
                            value(matrix_row, "conditioned_on_n"), 0
                        ),
                        "maf": number(value(matrix_row, "maf", "MAF")),
                        "mafClass": value(matrix_row, "maf_class"),
                        "direct": (snp_id, metabolite) in direct_pairs,
                    }
                )
        edge_pvalues = [
            parsed
            for edge in edges_by_gim[gim_id]
            if (parsed := number(value(edge, "p", "P", "p_value"))) is not None
        ]
        gim_entities.append(
            {
                "gimId": gim_id,
                "regionId": region_id,
                "snps": snp_list,
                "metabolites": metabolite_list,
                "nSnps": len(snp_list),
                "nMetabolites": len(metabolite_list),
                "nDirectEdges": len(direct_pairs),
                "minP": min(edge_pvalues) if edge_pvalues else None,
                "heatmap": heatmap,
            }
        )
    gim_entities.sort(
        key=lambda entity: (
            -entity["nMetabolites"],
            -entity["nSnps"],
            entity["minP"] if entity["minP"] is not None else 1,
            entity["gimId"],
        )
    )

    gims_by_region: dict[str, list[dict]] = defaultdict(list)
    snps_by_region: dict[str, set[str]] = defaultdict(set)
    metabolites_by_region: dict[str, set[str]] = defaultdict(set)
    for entity in gim_entities:
        gims_by_region[entity["regionId"]].append(entity)
        snps_by_region[entity["regionId"]].update(entity["snps"])
        metabolites_by_region[entity["regionId"]].update(entity["metabolites"])
    sentinels_by_region: dict[str, set[str]] = defaultdict(set)
    for row in region_metabolite_rows:
        region_id = value(row, "region_id")
        if not region_id:
            continue
        if sentinel := value(row, "sentinel_id", "snp_id"):
            sentinels_by_region[str(region_id)].add(str(sentinel))
        if metabolite := value(row, "metabolite", "trait"):
            metabolites_by_region[str(region_id)].add(str(metabolite))
    region_row_by_id = {
        str(value(row, "region_id")): row
        for row in region_rows
        if value(row, "region_id")
    }
    region_summary = {
        str(value(row, "region_id")): row
        for row in region_summary_rows
        if value(row, "region_id")
    }
    region_ids = sorted(
        set(region_row_by_id)
        | set(gims_by_region)
        | {
            str(value(row, "region_id"))
            for row in matrix_rows
            if value(row, "region_id")
        }
    )
    regions_data: list[dict] = []
    for region_id in region_ids:
        row = region_row_by_id.get(region_id, {})
        summary = region_summary.get(region_id, {})
        region_snps = snps_by_region.get(region_id, set())
        source_rows = [
            source_variants[snp_id]
            for snp_id in region_snps
            if snp_id in source_variants
        ]
        positions = [
            item["grch37Position"]
            for item in source_rows
            if item.get("grch37Position") is not None
        ]
        chromosome = value(row, "chromosome")
        if not chromosome and source_rows:
            chromosome = source_rows[0].get("chromosome")
        hg38_rows = [
            annotation_by_snp[snp_id].get("grch38") or {}
            for snp_id in region_snps
            if snp_id in annotation_by_snp
        ]
        hg38_positions = [
            item.get("position")
            for item in hg38_rows
            if item.get("position") is not None
        ]
        hg38_chromosome = next(
            (
                str(item.get("chromosome"))
                for item in hg38_rows
                if item.get("chromosome")
            ),
            str(chromosome) if chromosome else None,
        )
        start_hg38 = min(hg38_positions) if hg38_positions else None
        end_hg38 = max(hg38_positions) if hg38_positions else None
        midpoint_hg38 = (
            int((start_hg38 + end_hg38) / 2)
            if start_hg38 is not None and end_hg38 is not None
            else None
        )
        region_cytobands = sorted(
            {
                band
                for position in hg38_positions
                if (
                    band := cytoband_at(
                        cytobands, hg38_chromosome, position
                    )
                )
            }
        )
        regions_data.append(
            {
                "regionId": region_id,
                "chromosome": chromosome,
                "startGrch37": number(value(row, "start"))
                if value(row, "start")
                else (min(positions) if positions else None),
                "endGrch37": number(value(row, "end"))
                if value(row, "end")
                else (max(positions) if positions else None),
                "chromosomeGrch38": hg38_chromosome,
                "startGrch38": start_hg38,
                "endGrch38": end_hg38,
                "cytoband": cytoband_at(
                    cytobands, hg38_chromosome, midpoint_hg38
                ),
                "cytobands": region_cytobands,
                "cytobandAssembly": "GRCh38" if cytobands else None,
                "sentinels": sorted(
                    set(split_values(value(row, "sentinels")))
                    | sentinels_by_region.get(region_id, set())
                ),
                "metabolites": sorted(
                    set(split_values(value(row, "metabolites")))
                    | metabolites_by_region.get(region_id, set())
                ),
                "nGims": len(gims_by_region.get(region_id, [])),
                "nIndependentSignals": number(
                    value(summary, "n_independent_signals"), 0
                ),
                "nEdges": len(edges_by_region.get(region_id, [])),
                "snpIds": sorted(region_snps),
            }
        )
    region_by_id = {region["regionId"]: region for region in regions_data}

    def annotation_values(snp_id: str, field: str) -> list[str]:
        return annotation_by_snp.get(snp_id, {}).get(field, []) or []

    locus_variant_metabolite_rows: list[dict] = []
    for entity in gim_entities:
        region = region_by_id.get(entity["regionId"], {})
        for snp_id in entity["snps"]:
            locus_variant_metabolite_rows.append(
                {
                    "locus_id": entity["regionId"],
                    "chromosome": region.get("chromosome"),
                    "start_grch37": region.get("startGrch37"),
                    "end_grch37": region.get("endGrch37"),
                    "GIMsID": entity["gimId"],
                    "variantID": snp_id,
                    "gene_symbols": "; ".join(
                        annotation_values(snp_id, "geneSymbols")
                    ),
                    "ensembl_gene_ids": "; ".join(
                        annotation_values(snp_id, "geneIds")
                    ),
                    "transcript_ids": "; ".join(
                        annotation_values(snp_id, "transcriptIds")
                    ),
                    "metabolites": "; ".join(entity["metabolites"]),
                    "n_metabolites": entity["nMetabolites"],
                    "n_variants_in_GIM": entity["nSnps"],
                }
            )
    locus_summary_rows: list[dict] = []
    for region in regions_data:
        region_gims = gims_by_region.get(region["regionId"], [])
        if not region_gims:
            continue
        region_snp_ids = sorted(
            {snp_id for entity in region_gims for snp_id in entity["snps"]}
        )
        locus_summary_rows.append(
            {
                "locus_id": region["regionId"],
                "chromosome": region.get("chromosome"),
                "start_grch37": region.get("startGrch37"),
                "end_grch37": region.get("endGrch37"),
                "GIMsID": "; ".join(entity["gimId"] for entity in region_gims),
                "variantID": " || ".join(
                    f"{entity['gimId']}: {'; '.join(entity['snps'])}"
                    for entity in region_gims
                ),
                "gene_symbols": "; ".join(
                    sorted(
                        {
                            gene
                            for snp_id in region_snp_ids
                            for gene in annotation_values(snp_id, "geneSymbols")
                        }
                    )
                ),
                "ensembl_gene_ids": "; ".join(
                    sorted(
                        {
                            gene
                            for snp_id in region_snp_ids
                            for gene in annotation_values(snp_id, "geneIds")
                        }
                    )
                ),
                "transcript_ids": "; ".join(
                    sorted(
                        {
                            transcript
                            for snp_id in region_snp_ids
                            for transcript in annotation_values(
                                snp_id, "transcriptIds"
                            )
                        }
                    )
                ),
                "metabolites": " || ".join(
                    f"{entity['gimId']}: {'; '.join(entity['metabolites'])}"
                    for entity in region_gims
                ),
                "n_GIMs": len(region_gims),
            }
        )

    gim_metabolites = {
        metabolite for entity in gim_entities for metabolite in entity["metabolites"]
    }
    cancer_annotations: dict[str, list[dict]] = defaultdict(list)
    for row in differential_rows:
        metabolite = row.get("metabolite")
        if metabolite not in gim_metabolites:
            continue
        beta = number(row.get("beta_case_minus_other_sd"))
        fdr = number(row.get("fdr_bh_within_comparison"))
        cancer_annotations[metabolite].append(
            {
                "comparison": row.get("comparison"),
                "caseDefinition": row.get("case_definition"),
                "controlDefinition": row.get("control_definition"),
                "betaCaseMinusOtherSd": beta,
                "direction": (
                    "higher_in_case"
                    if beta is not None and beta > 0
                    else "lower_in_case"
                    if beta is not None and beta < 0
                    else "no_direction"
                ),
                "pValue": number(row.get("p_value")),
                "fdr": fdr,
                "fdrSignificant": bool(fdr is not None and fdr <= 0.05),
                "nCase": number(row.get("n_case")),
                "nOther": number(row.get("n_other")),
                "assayCoverage": row.get("assay_coverage"),
            }
        )
    for rows in cancer_annotations.values():
        rows.sort(
            key=lambda item: (
                not item["fdrSignificant"],
                item["fdr"] if item["fdr"] is not None else 1.0,
            )
        )

    alpha_manifest = [
        {
            "snp_id": snp["snpId"],
            "chromosome": snp["alphaGenome"]["input"]["chromosome"],
            "position": snp["alphaGenome"]["input"]["position"],
            "reference_bases": snp["alphaGenome"]["input"]["reference"],
            "alternate_bases": ",".join(
                snp["alphaGenome"]["input"]["alternates"]
            ),
            "assembly": "GRCh38",
        }
        for snp in snp_records
        if snp["alphaGenome"]["input"]["position"]
        and snp["alphaGenome"]["input"]["reference"]
        and snp["alphaGenome"]["input"]["alternates"]
    ]

    out_dir.mkdir(parents=True, exist_ok=True)
    _json_write(out_dir / "gim_records.json", records)
    _json_write(out_dir / "gim_snps.json", snp_records)
    _json_write(out_dir / "gim_regions.json", regions_data)
    _json_write(out_dir / "gim_entities.json", gim_entities)
    _json_write(
        out_dir / "metabolite_cancer_annotations.json", cancer_annotations
    )
    with (out_dir / "gim_locus_variant_metabolites.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "locus_id",
                "chromosome",
                "start_grch37",
                "end_grch37",
                "GIMsID",
                "variantID",
                "gene_symbols",
                "ensembl_gene_ids",
                "transcript_ids",
                "metabolites",
                "n_metabolites",
                "n_variants_in_GIM",
            ],
        )
        writer.writeheader()
        writer.writerows(locus_variant_metabolite_rows)
    with (out_dir / "gim_locus_metabolites_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "locus_id",
                "chromosome",
                "start_grch37",
                "end_grch37",
                "GIMsID",
                "variantID",
                "gene_symbols",
                "ensembl_gene_ids",
                "transcript_ids",
                "metabolites",
                "n_GIMs",
            ],
        )
        writer.writeheader()
        writer.writerows(locus_summary_rows)
    with (out_dir / "alphagenome_input.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "snp_id",
                "chromosome",
                "position",
                "reference_bases",
                "alternate_bases",
                "assembly",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(alpha_manifest)

    stats = {
        "schemaVersion": 2,
        "generatedAt": datetime.now(timezone.utc).isoformat(),
        "bundleId": manifest.get("bundleId"),
        "sourceBuild": source_build,
        "annotationBuild": "GRCh38 (Ensembl variation mapping)",
        "annotationStatus": annotation_status,
        "nRegions": len(regions_data),
        "nGimAssociations": len(records),
        "nGimEntities": len(gim_entities),
        "nIndependentSnps": len(snp_records),
        "nAlphaGenomeReady": len(alpha_manifest),
        "nAnnotatedSnps": sum(bool(row["grch38"]) for row in snp_records),
        "nConsequences": sum(
            bool(row["mostSevereConsequence"]) for row in snp_records
        ),
        "sources": [
            "GIMForge result tables",
            f"Variant coordinates ({source_build})",
            "Ensembl REST variation and VEP annotations",
            "AlphaGenome GRCh38 input manifest",
        ],
    }
    _json_write(out_dir / "portal_stats.json", stats, indent=2)
    print(
        f"Wrote {len(records):,} GIM associations, {len(snp_records):,} SNPs, "
        f"{len(gim_entities):,} GIMs, and {len(alpha_manifest):,} "
        f"AlphaGenome-ready inputs to {out_dir}"
    )
    return stats
