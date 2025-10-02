#!/usr/bin/env python3
"""
SKIPME: Strategic Knowledge-based In-silico Prediction of Modifiable Exons
"""
import os
from flask import Flask, request, jsonify, render_template
import time
import re
import requests
from typing import Dict, Any, Optional, Tuple, List
from Bio.Seq import Seq

# --- Template Setup ---
# This section will automatically create the necessary HTML files in a 'templates' folder.

def setup_templates():
    """Creates the templates directory and HTML files if they don't exist."""
    if not os.path.exists('templates'):
        os.makedirs('templates')

    # Base template with navigation and footer
    base_html = """
<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{{ title }} - SKIPME</title>
    <style>
        body{font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif;margin:0;background-color:#f8f9fa;color:#333}
        .container{max-width:800px;margin:auto;background:white;padding:2em;border-radius:8px;box-shadow:0 4px 8px rgba(0,0,0,.1); margin-top: 2em; margin-bottom: 2em;}
        h3{color:#005a9c} a{color:#005a9c}
        nav {background-color: #333; padding: 1em; text-align: center;}
        nav a {color: white; margin: 0 15px; text-decoration: none; font-weight: bold;}
        nav a:hover {text-decoration: underline;}
        .disclaimer{font-size:.8em;color:#6c757d;margin-top:2em;text-align:center; padding-top: 1em; border-top: 1px solid #eee;}
        /* Tool-specific styles below */
        form{display:grid;grid-template-columns:1fr;gap:1em;align-items:center;margin-bottom:2em}
        input,button{padding:.75em;border-radius:4px;border:1px solid #ccc;font-size:1em;width:100%;box-sizing:border-box}
        button{background-color:#007bff;color:white;font-weight:bold;cursor:pointer;border:none;width:auto;justify-self:end;padding:.75em 1.5em}
        button:hover{background-color:#0056b3}
        #loader{display:none;text-align:center;padding:1em;font-size:1.2em}
        #results{display:none}.result-header{padding:1em;color:white;border-radius:4px;margin-bottom:1em}.result-header h4{margin:0;font-size:1.5em;text-align:center}.result-block{margin-bottom:1.5em}.result-block h5{margin-bottom:.5em;color:#495057;border-bottom:1px solid #eee;padding-bottom:.3em}.checklist{list-style-type:none;padding:0}.checklist li{margin-bottom:.5em}.pass::before{content:'✔';color:#28a745;margin-right:10px;font-weight:bold}.fail::before{content:'✖';color:#dc3545;margin-right:10px;font-weight:bold}.likely-eligible{background-color:#17a2b8}.unlikely-eligible{background-color:#ffc107;color:#333}.not-eligible{background-color:#dc3545}.unable-to-assess{background-color:#6c757d}.error{background-color:#dc3545}.note{font-style:italic;color:#555;font-size:.9em}
    </style>
</head>
<body>
    <nav>
        <a href="/">Tool</a>
        <a href="/about">About/Methods</a>
        <a href="/cite">How to Cite</a>
    </nav>
    <main class="container">
        {% block content %}{% endblock %}
        <p class="disclaimer">This tool is for informational purposes only and is not for clinical use. Results require manual verification.</p>
    </main>
</body>
</html>
    """
    
    # Index page with the tool
    index_html = """
{% extends "base.html" %}
{% block content %}
<h3>SKIPME: Strategic Knowledge-based In-silico Prediction of Modifiable Exons</h3>
<p>Enter a variant to assess its associated exon for ASO-mediated skipping.</p>
<p style="font-size: 0.9em; color: #333;">This tool specifically assesses eligibility for <strong>exon skipping</strong> to restore a reading frame. For other ASO strategies (e.g., allele-specific knockdown,WT-upregulation) or for variants with different functional effects, please consult the <a href="https://shorturl.at/YqphL" target="_blank" rel="noopener noreferrer"><b>guidelines</b></a>.</p>
<p><b>Example Formats:</b></p>
<ul>
    <li><code>BRAF c.1799T>A</code></li>
    <li><code>NM_015247.4:c.1054G>A</code></li>
</ul>
<form id="assessment-form">
    <label for="query">Variant:</label>
    <input id="query" required placeholder="e.g., NM_015247.4:c.1054G>A">
    <button type="submit">Assess</button>
</form>
<div id="loader">Assessing...</div>
<div id="results"></div>
<script>
document.getElementById('assessment-form').addEventListener('submit',async function(e){e.preventDefault();const t=document.getElementById('results'),n=document.getElementById('loader');t.style.display='none',n.style.display='block';const s={query:document.getElementById('query').value};try{const e=await fetch('/assess',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(s)}),a=await e.json();displayResults(a)}catch(e){displayResults({classification:"Error",reason:"Could not connect to the server."})}finally{n.style.display='none'}});function displayResults(e){const t=document.getElementById('results');t.style.display='block';const n=e.classification.toLowerCase().replace(/ /g,'-');if(e.classification.includes("Error")||e.classification.includes("Unable")){return void(t.innerHTML=`<div class="result-header ${n}"><h4>${e.classification}</h4></div><p>${e.reason}</p>`);}let s=`<div class="result-header ${n}"><h4>${e.classification}</h4></div>`;s+=`<p><strong>Reason:</strong> ${e.reason}</p>`,s+=`<div class="result-block"><h5>Assessment Summary</h5><p class="note">Analysis performed on transcript <strong>${e.transcript_id}</strong>, as determined by the variant's effect.</p><ul><li><strong>Gene:</strong> ${e.gene}</li><li><strong>Query:</strong> ${e.variant}</li><li><strong>Target Exon:</strong> Total #${e.total_exon_number} (Coding #${e.coding_exon_number}) | ID: ${e.exon_id}</li><li><strong>Location:</strong> ${e.location}</li></ul></div>`,s+='<div class="result-block"><h5>Guideline Checks</h5><ul class="checklist">';for(const[t,n]of Object.entries(e.checks))s+=`<li class="${n?'pass':'fail'}">${t}</li>`;s+='</ul></div>',s+=`<div class="result-block"><h5>Evidence from Databases</h5><ul><li><strong>Fraction of Protein:</strong> ${e.frac_cds}</li><li><strong>Overlapping Protein Domains:</strong> ${e.domain_count}</li><li><strong>Pathogenic Variants in Exon (on this transcript):</strong><ul><li>Missense: ${e.pathogenic_variant_counts.missense}</li><li>Nonsense: ${e.pathogenic_variant_counts.nonsense}</li><li>Frameshift: ${e.pathogenic_variant_counts.frameshift}</li><li>In-frame Deletions: ${e.pathogenic_variant_counts.inframe_del}</li><li>Splice Site: ${e.pathogenic_variant_counts.splice}</li></ul></li></ul></div>`,t.innerHTML=s}
</script>
{% endblock %}
    """

    # About page
    about_html = """
{% extends "base.html" %}
{% block content %}
<h3>About & Methods</h3>
<h4>Purpose</h4>
<p>SKIPME (Strategic Knowledge-based In-silico Prediction of Modifiable Exons) is a research tool designed to assess the eligibility of a specific non dominant-inherited exonic variant for therapeutic exon skipping using antisense oligonucleotides (ASOs). It automates the analysis of key criteria to predict whether skipping an exon is likely to restore a functional protein product, and assumes variants are exonic nonsense, missense or frameshift variants with no splice altering effects. For other variants please refer to the <a href="https://shorturl.at/YqphL" target="_blank" rel="noopener noreferrer"><b>guidelines</b></a></p>

<h4>Methodology</h4>
<p>The tool follows a systematic, rules-based workflow:</p>
<ol>
    <li><b>Variant Interpretation:</b> An input exonic variant in HGVS format (e.g., <code>GENE c.123A>G</code>) is submitted to the Ensembl Variant Effect Predictor (VEP) to identify the affected transcript and its genomic coordinates. The tool prioritizes MANE select or canonical transcripts.</li>
    <li><b>Exon Mapping:</b> The variant's genomic position is mapped to a specific exon within the chosen transcript. This is done by coordinate overlap.</li>
    <li><b>Eligibility Assessment:</b> The target exon is evaluated against a series of bioinformatic criteria derived from established guidelines for ASO-mediated exon skipping:</li>
    <ul>
        <li><b>Reading Frame Maintenance:</b> The coding sequence (CDS) length of the exon must be a multiple of three to ensure the downstream reading frame is preserved.</li>
        <li><b>Premature Stop Codons:</b> An in-silico translation is performed on the hypothetical transcript lacking the target exon. If a new premature stop codon is generated, the exon is deemed ineligible.</li>
        <li><b>Functional Importance:</b> The tool checks for overlaps with known protein domains (via Ensembl/InterPro) and flags exons that are mutational hotspots (containing multiple pathogenic missense variants from ClinVar).</li>
        <li><b>Genomic Context:</b> The first and last coding exons are considered ineligible as they are essential for translation initiation and termination. The tool also checks for existing pathogenic splice or in-frame deletion variants within the exon, as their presence suggests that loss of the exon is itself disease-causing.</li>
    </ul>
</ol>

<h4>Data Sources</h4>
<p>SKIPME relies entirely on publicly available data accessed via the <b>Ensembl REST API</b> (GRCh38). This includes gene/transcript structures, sequence data, variant annotations, clinical significance from ClinVar, and protein domain information from databases like CDD and InterPro.</p>

<h4>Limitations</h4>
<p>This is a computational prediction tool intended for <b>research purposes only</b> and is <b>not a substitute for clinical advice</b>. Its predictions are dependent on the accuracy and completeness of the underlying public databases. Experimental validation is required to confirm any in-silico findings.</p>
{% endblock %}
    """

    # Cite page
    cite_html = """
{% extends "base.html" %}
{% block content %}
<h3>How to Cite</h3>
<p>SKIPME is currently being prepared for publication. If you use this tool in your research, please cite this website until a formal publication is available.</p>
<p><strong>(Schober, 2025). SKIPME: Strategic Knowledge-based In-silico Prediction of Modifiable Exons. Retrieved from https://skipme.onrender.com/.</strong></p>
<hr>
<p>Once the software is formally published, we recommend citing both the paper and the specific version of the software used. The code for this tool is open-source and has been archived on Zenodo.</p>
<p><strong>Software DOI:</strong></p>
{% endblock %}
    """

    # Write files with explicit UTF-8 encoding
    with open('templates/base.html', 'w', encoding='utf-8') as f: f.write(base_html)
    with open('templates/index.html', 'w', encoding='utf-8') as f: f.write(index_html)
    with open('templates/about.html', 'w', encoding='utf-8') as f: f.write(about_html)
    with open('templates/cite.html', 'w', encoding='utf-8') as f: f.write(cite_html)


# Run the setup function to ensure templates exist
setup_templates()

# ===== CONFIGURATION =====
ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
# ===== end configuration =====

app = Flask(__name__)

# ----------------- Ensembl Client -----------------
class EnsemblClient:
    def __init__(self, base_url=ENSEMBL_REST, headers=HEADERS, delay=0.1):
        self.base_url = base_url.rstrip('/')
        self.session = requests.Session()
        self.session.headers.update(headers)
        self.delay = delay

    def _get(self, path, params=None, max_retries=5):
        url = f"{self.base_url}{path}"
        backoff = 1.0
        for attempt in range(max_retries):
            time.sleep(self.delay)
            try:
                resp = self.session.get(url, params=params, timeout=30)
                if resp.status_code == 200:
                    try: return resp.json()
                    except ValueError: return resp.text
                elif resp.status_code in (429, 503):
                    wait = float(resp.headers.get('Retry-After', backoff))
                    time.sleep(wait); backoff *= 2
                elif 500 <= resp.status_code < 600:
                    time.sleep(backoff); backoff *= 2
                else:
                    if resp.status_code >= 400 and resp.status_code < 500:
                        return None
            except requests.RequestException:
                if attempt + 1 == max_retries: raise
                time.sleep(backoff); backoff *= 2
        return None

    def lookup_id_expand(self, identifier):
        return self._get(f"/lookup/id/{identifier}", params={'expand': '1'})

    def vep_hgvs(self, hgvs_string):
        return self._get(f"/vep/human/hgvs/{hgvs_string.strip()}", params={'variant_class': 1})

    def get_cds_sequence(self, transcript_id):
        data = self._get(f"/sequence/id/{transcript_id}", params={"type": "cds"})
        return data.get("seq") if isinstance(data, dict) else None
        
    def get_domains(self, protein_id):
        path = f"/overlap/translation/{protein_id}"
        params = {"feature": "protein_feature"}
        all_features = self._get(path, params=params)
        
        if not all_features or not isinstance(all_features, list):
            return []
            
        domain_sources = {'CDD'}
        preliminary_domains = [f for f in all_features if f.get('type') in domain_sources]
        
        unique_interpro_domains = {}
        for feature in preliminary_domains:
            interpro_id = feature.get('interpro')
            if interpro_id and interpro_id not in unique_interpro_domains:
                unique_interpro_domains[interpro_id] = feature
        
        return list(unique_interpro_domains.values())

    def overlap_region_variation(self, chrom, start, end):
        data = self._get(f"/overlap/region/human/{chrom}:{start}-{end}", params={'feature': 'variation'})
        return data if isinstance(data, list) else []

# ----------------- Transcript & Exon Logic -----------------

def choose_best_consequence(consequences: List[Dict[str, Any]], canonical_id: Optional[str] = None) -> Optional[Dict[str, Any]]:
    if not consequences:
        return None

    valid_consequences = consequences
    
    mane_select = [c for c in valid_consequences if c.get('mane_select')]
    if mane_select: return mane_select[0]

    if canonical_id:
        canonical_id_base = canonical_id.split('.')[0]
        for c in valid_consequences:
            if c.get('transcript_id', '').startswith(canonical_id_base):
                return c
    
    coding_consequences = sorted(
        [c for c in valid_consequences if c.get('biotype') == 'protein_coding' and c.get('cds_end')],
        key=lambda c: c['cds_end'] - c['cds_start'] if c.get('cds_start') else -1,
        reverse=True
    )
    if coding_consequences:
        return coding_consequences[0]

    return valid_consequences[0]

def extract_exons_from_transcript(transcript: Dict[str, Any]):
    exons_raw = sorted(transcript.get('Exon', []), key=lambda e: e['start'])
    strand = transcript.get('strand')

    if strand == -1:
        exons_raw.reverse()

    translation_data = transcript.get('Translation', {})
    cds_start, cds_end = translation_data.get('start'), translation_data.get('end')
    seq_region = transcript.get('seq_region_name')
    
    normalized, coding_exon_count = [], 0
    for i, e in enumerate(exons_raw, 1):
        start, end = e['start'], e['end']
        cds_len_of_exon, is_coding = 0, False
        if cds_start and cds_end:
            overlap_start, overlap_end = max(start, cds_start), min(end, cds_end)
            if overlap_end >= overlap_start:
                cds_len_of_exon = overlap_end - overlap_start + 1
                is_coding = True
                coding_exon_count += 1
        
        normalized.append({
            'total_exon_number': i, 
            'coding_exon_number': coding_exon_count if is_coding else None,
            'exon_id': e.get('id'), 
            'start': start, 
            'end': end, 
            'seq_region_name': seq_region,
            'cds_length': cds_len_of_exon
        })
    return normalized

# ----------------- Assessment Logic -----------------

def assess_single_exon(client, original_query, transcript, all_exons, target_exon):
    gene_symbol = transcript.get('Parent')
    transcript_id = transcript.get('id')
    protein_id = transcript.get("Translation", {}).get("id")
    cds_seq = client.get_cds_sequence(transcript_id)
    
    coding_exons = [e for e in all_exons if e['cds_length'] > 0]
    total_coding_exons = len(coding_exons)
    total_cds_len = sum(e['cds_length'] for e in coding_exons)

    if not target_exon.get('coding_exon_number'):
        return {"classification": "Unable to Assess", "reason": f"The variant maps to exon {target_exon['total_exon_number']}, which is non-coding (in the UTR)."}
    
    chrom, start, end = target_exon['seq_region_name'], target_exon['start'], target_exon['end']
    variants_in_region = client.overlap_region_variation(chrom, start, end)
    
    counts = {'missense': 0, 'inframe_del': 0, 'splice': 0, 'nonsense': 0, 'frameshift': 0}
    for v in variants_in_region:
        clclass = classify_variant_clinsig(v.get('clinical_significance'))
        if clclass == "pathogenic":
            conseq = (v.get("consequence_type") or "").lower()
            if "missense" in conseq: counts['missense'] += 1
            elif "inframe_deletion" in conseq: counts['inframe_del'] += 1
            elif "splice_donor" in conseq or "splice_acceptor" in conseq: counts['splice'] += 1
            elif "stop_gained" in conseq: counts['nonsense'] += 1
            elif "frameshift" in conseq: counts['frameshift'] += 1

    exon_cds_len = target_exon['cds_length']
    coding_exon_number = target_exon['coding_exon_number']
    
    cond1_inframe = (exon_cds_len % 3 == 0)
    
    cond2_no_stop = False
    if cds_seq:
        try:
            cds_map, current_pos = {}, 0
            sorted_coding_exons = sorted(coding_exons, key=lambda x: x['coding_exon_number'])
            for exon in sorted_coding_exons:
                length, exon_num = exon['cds_length'], exon['coding_exon_number']
                cds_map[exon_num] = cds_seq[current_pos : current_pos + length]
                current_pos += length
            
            skipped_cds = "".join(cds_map[i] for i in sorted(cds_map.keys()) if i != coding_exon_number)
            if skipped_cds:
                prot = str(Seq(skipped_cds).translate(to_stop=False))
                cond2_no_stop = "*" not in prot[:-1]
        except Exception:
            cond2_no_stop = False

    cond3_not_terminal = (coding_exon_number is not None and coding_exon_number not in (1, total_coding_exons))
    cond4_small = (exon_cds_len / total_cds_len) < 0.10 if total_cds_len > 0 else False
    
    domain_count = 0
    if protein_id:
        domains = client.get_domains(protein_id)
        if domains and isinstance(domains, list):
            cds_pos_start = sum(e['cds_length'] for e in sorted(coding_exons, key=lambda x: x['coding_exon_number']) if e['coding_exon_number'] < coding_exon_number)
            exon_aa_start = (cds_pos_start // 3) + 1
            exon_aa_end = ((cds_pos_start + exon_cds_len -1) // 3) + 1
            for d in domains:
                if d.get('start', 0) <= exon_aa_end and d.get('end', 0) >= exon_aa_start: domain_count += 1

    cond5_no_domain = domain_count == 0
    cond6_missense = counts['missense'] < 3
    cond7_splice = counts['splice'] == 0
    cond8_no_inframe_del = counts['inframe_del'] == 0
    
    classification, reason = "Undetermined", ""
    if not cond3_not_terminal: classification, reason = "Not Eligible", "Exon is the first or last coding exon."
    elif not cond1_inframe: classification, reason = "Not Eligible", "Exon is out-of-frame, which would disrupt the reading frame."
    elif not cond2_no_stop: classification, reason = "Not Eligible", "Skipping this exon is predicted to create a premature stop codon."
    elif not cond7_splice: classification, reason = "Not Eligible", f"Exon contains {counts['splice']} pathogenic splice variant(s), indicating exon loss is pathogenic."
    elif not cond8_no_inframe_del: classification, reason = "Not Eligible", f"Exon contains {counts['inframe_del']} pathogenic in-frame deletion(s)."
    elif not cond6_missense: classification, reason = "Unlikely Eligible", f"Exon is a mutational hotspot with {counts['missense']} pathogenic missense variants."
    elif not cond4_small: classification, reason = "Unlikely Eligible", "Exon constitutes >10% of the protein, risking major functional loss."
    elif not cond5_no_domain: classification, reason = "Unlikely Eligible", f"Exon overlaps with {domain_count} protein domain(s)."
    else: classification, reason = "Likely Eligible", "Exon meets the primary criteria for a skippable exon. To see if there is already an ASO existing look at the <a href=https://generegistry.n1collaborative.org/index.html>N1C-registry</a>."
    
    return {
        "gene": gene_symbol, "variant": original_query, "transcript_id": transcript_id, "exon_id": target_exon['exon_id'],
        "coding_exon_number": coding_exon_number, "total_exon_number": target_exon['total_exon_number'],
        "location": f"chr{target_exon['seq_region_name']}:{target_exon['start']}-{target_exon['end']}",
        "frac_cds": f"{(exon_cds_len / total_cds_len * 100):.2f}%" if total_cds_len > 0 else "N/A",
        "pathogenic_variant_counts": counts, "domain_count": domain_count,
        "checks": {"Is In-Frame": cond1_inframe, "No New Stop Codon": cond2_no_stop, "Not First/Last Exon": cond3_not_terminal,
                   "Is <10% of Protein": cond4_small, "No Domain Overlap": cond5_no_domain, "Low Missense Count (<3)": cond6_missense,
                   "No Pathogenic Splice Variants": cond7_splice, "No Pathogenic In-Frame Deletions": cond8_no_inframe_del},
        "classification": classification, "reason": reason
    }

def classify_variant_clinsig(clinsig_field):
    if clinsig_field is None: return 'other'
    vals = [v.lower() for v in (clinsig_field if isinstance(clinsig_field, list) else [clinsig_field]) if isinstance(v, str)]
    if any('pathogenic' in v for v in vals) and not any('likely' in v for v in vals): return 'pathogenic'
    if any('likely pathogenic' in v for v in vals): return 'pathogenic'
    if any('uncertain' in v for v in vals): return 'VUS'
    return 'other'

def parse_hgvs_query(query: str) -> Tuple[Optional[str], Optional[str]]:
    query = query.strip()
    match = re.search(r'([A-Z0-9\-_.:()]+)\s*:\s*([cgnmp]\..*)', query, re.IGNORECASE)
    if match: return f"{match.group(1)}:{match.group(2)}", match.group(1).split('(')[0]
    match = re.search(r'([A-Z0-9\-_]+)\s+([cgnmp]\..*)', query, re.IGNORECASE)
    if match: return f"{match.group(1)}:{match.group(2)}", None
    return None, None

# ----------------- Main Flask Routes -----------------
@app.route('/')
def index():
    return render_template('index.html', title="Tool")

@app.route('/about')
def about():
    return render_template('about.html', title="About/Methods")

@app.route('/cite')
def cite():
    return render_template('cite.html', title="How to Cite")

@app.route('/assess', methods=['POST'])
def assess_variant():
    data = request.get_json()
    query = data.get('query')
    if not query:
        return jsonify({"classification": "Error", "reason": "Query string is required."}), 400

    client = EnsemblClient()
    try:
        hgvs_query, transcript_id_from_query = parse_hgvs_query(query)
        if not hgvs_query:
            return jsonify({"classification": "Error", "reason": "Invalid input format. Use 'GENE c.123A>G' or 'NM_12345.6:c.123A>G'."}), 400

        vep_data = client.vep_hgvs(hgvs_query)
        if not vep_data or not isinstance(vep_data, list):
            return jsonify({"classification": "Error", "reason": f"VEP analysis failed for '{hgvs_query}'. Variant may be invalid."}), 404

        transcript_consequences = vep_data[0].get('transcript_consequences', [])
        target_consequence = choose_best_consequence(transcript_consequences, transcript_id_from_query)
        if not target_consequence:
            return jsonify({"classification": "Unable to Assess", "reason": f"Variant '{hgvs_query}' does not affect a known transcript."}), 400
        
        definitive_transcript_id = target_consequence['transcript_id']
        transcript_data = client.lookup_id_expand(definitive_transcript_id)
        if not transcript_data:
            return jsonify({"classification": "Error", "reason": f"Could not fetch data for transcript {definitive_transcript_id}."}), 500

        all_exons = extract_exons_from_transcript(transcript_data)
        target_exon = None
        variant_chrom, variant_start, variant_end = vep_data[0]['seq_region_name'], vep_data[0]['start'], vep_data[0]['end']

        for exon in all_exons:
            if exon['seq_region_name'] == variant_chrom and max(variant_start, exon['start']) <= min(variant_end, exon['end']):
                target_exon = exon
                break

        if not target_exon:
            consequence_terms = set(target_consequence.get('consequence_terms', []))
            if 'splice_donor_variant' in consequence_terms:
                adjacent_exons = [ex for ex in all_exons if ex['end'] < variant_start and ex['seq_region_name'] == variant_chrom]
                if adjacent_exons: target_exon = max(adjacent_exons, key=lambda ex: ex['end'])
            elif 'splice_acceptor_variant' in consequence_terms:
                adjacent_exons = [ex for ex in all_exons if ex['start'] > variant_end and ex['seq_region_name'] == variant_chrom]
                if adjacent_exons: target_exon = min(adjacent_exons, key=lambda ex: ex['start'])

        if not target_exon:
            target_exon_id = target_consequence.get('exon_id')
            if target_exon_id:
                target_exon = next((ex for ex in all_exons if ex['exon_id'] == target_exon_id), None)

        if not target_exon and 'exon' in target_consequence:
            try:
                target_exon_rank = int(target_consequence['exon'].split('/')[0])
                target_exon = next((ex for ex in all_exons if ex['total_exon_number'] == target_exon_rank), None)
            except (ValueError, IndexError): pass

        if not target_exon:
            return jsonify({"classification": "Error", "reason": f"Could not map variant to an exon on transcript {definitive_transcript_id}."}), 500

        result = assess_single_exon(client, query, transcript_data, all_exons, target_exon)
        return jsonify(result)

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"classification": "Error", "reason": f"An unexpected server error occurred: {str(e)}"}), 500

if __name__ == '__main__':
    app.run(debug=True)