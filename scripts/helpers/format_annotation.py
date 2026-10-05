import sys
import collections

sqanti_file = sys.argv[1]
orfanage_file = sys.argv[2]
cpat_file = sys.argv[3]
gencode_file = sys.argv[4]
gffcompare_file = sys.argv[5]
out_file = sys.argv[6]

DIVERGENT_MAX_DIST = 1000 # distance to TSS to annotate DT lncRNA

MIN_CDS_NT = 300       # ORFs <100 codons (300 nt) -> not protein-coding
MIN_TRANSCRIPT_NT = 200  # spliced length threshold for lncRNA
MIN_CDS_FRACTION = 0.05   # CDS + its introns must cover >=5% of transcript genomic span
NMD_LAST_JUNCTION_NT = 50 # PTC >50 nt upstream of the last exon-exon junction -> NMD

ISOFORM_CATS = {'full-splice_match', 'incomplete-splice_match',
                'novel_in_catalog', 'novel_not_in_catalog'}
LOCUS_PRIORITY = ('antisense', 'genic', 'genic_intron', 'intergenic')

def parse_attrs(s):
    out = {}
    for tok in s.strip().rstrip(';').split(';'):
        tok = tok.strip()
        if not tok: continue
        if '"' in tok:
            k, _, v = tok.partition(' ')
            out[k.strip()] = v.strip().strip('"')
        elif '=' in tok:
            k, _, v = tok.partition('=')
            out[k.strip()] = v.strip().strip('"')
    return out

gencode = {}
gencode_by_name = {}
all_genes_by_chrom = collections.defaultdict(list)
pc_genes_by_chrom = collections.defaultdict(list)

with open(gencode_file) as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 9 or p[2] != 'gene': continue
        a = parse_attrs(p[8])
        gid = a.get('gene_id', '')
        if not gid: continue
        info = {'id': gid, 'name': a.get('gene_name', gid),
                'type': a.get('gene_type', ''),
                'chrom': p[0], 'start': int(p[3]), 'end': int(p[4]), 'strand': p[6]}
        gencode[gid] = info
        gencode[gid.split('.')[0]] = info
        if info['name']: gencode_by_name[info['name']] = info
        all_genes_by_chrom[p[0]].append(info)
        if info['type'] == 'protein_coding':
            pc_genes_by_chrom[p[0]].append(info)

sqanti = {}
with open(sqanti_file) as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for line in f:
        row = dict(zip(hdr, line.rstrip('\n').split('\t')))
        if 'isoform' in row:
            sqanti[row['isoform']] = row

gffcompare = {}
with open(gffcompare_file) as f:
    for line in f:
        if line.startswith('#'): continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 9 or p[2] != 'transcript': continue
        a = parse_attrs(p[8])
        tid = a.get('transcript_id', '')
        if not tid: continue
        gffcompare[tid] = {'ref_gene_id': a.get('ref_gene_id', ''),
                          'gene_name': a.get('gene_name', ''),
                          'cmp_ref': a.get('cmp_ref', ''),
                          'class_code': a.get('class_code', '')}

transcripts = collections.defaultdict(
    lambda: {'chrom': '', 'strand': '', 'gene_id': '', 'exons': [], 'cds': []})

def load_gtf(path, ignore_tids=None):
    ignore_tids = ignore_tids or set()
    
    with open(path) as f:
        for line in f:
            if line.startswith('#'): continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9: continue
            try: start, end = int(p[3]), int(p[4])
            except ValueError: continue
            a = parse_attrs(p[8])
            tid = a.get('transcript_id', '')
            
            if not tid or tid in ignore_tids: 
                continue
                
            t = transcripts[tid]
            t['chrom'], t['strand'] = p[0], p[6]
            if a.get('gene_id'): t['gene_id'] = a['gene_id']
            if p[2] == 'exon': t['exons'].append((start, end))
            elif p[2] == 'CDS': t['cds'].append((start, end, p[7]))


load_gtf(orfanage_file)
orfanage_tids = set(transcripts.keys())
load_gtf(cpat_file, ignore_tids=orfanage_tids)

# length and ORF filters

def cds_length(t):
    return sum(e - s + 1 for s, e, _ in t['cds'])

def transcript_length(t):
    return sum(e - s + 1 for s, e in t['exons'])

def cds_genomic_span(t):
    if not t['cds']: return 0
    return max(c[1] for c in t['cds']) - min(c[0] for c in t['cds']) + 1

def transcript_genomic_span(t):
    if not t['exons']: return 0
    return max(e[1] for e in t['exons']) - min(e[0] for e in t['exons']) + 1

def is_nmd(t):
    if not t['cds'] or len(t['exons']) < 2: return False
    exons = sorted(t['exons'])
    cmin = min(c[0] for c in t['cds'])
    cmax = max(c[1] for c in t['cds'])
    if t['strand'] == '-':
        last_exon = exons[0]
        utr3 = sum(min(e, cmin - 1) - s + 1 for s, e in exons if s < cmin)
    else:
        last_exon = exons[-1]
        utr3 = sum(e - max(s, cmax + 1) + 1 for s, e in exons if e > cmax)
    last_len = last_exon[1] - last_exon[0] + 1
    return (utr3 - last_len) > NMD_LAST_JUNCTION_NT

def classify_transcript(t):
    if cds_length(t) < MIN_CDS_NT: return 'lncRNA'
    if cds_genomic_span(t) < MIN_CDS_FRACTION * transcript_genomic_span(t): return 'lncRNA'
    if is_nmd(t): return 'nonsense_mediated_decay'
    return 'protein_coding'

def is_protein_coding(t):
    return classify_transcript(t) == 'protein_coding'

for tid in list(transcripts.keys()):
    t = transcripts[tid]
    if not t['exons'] or transcript_length(t) < MIN_TRANSCRIPT_NT:
        del transcripts[tid]

def to_gencode_id(name):
    if not name or name == '.' or name.startswith('novelGene'): return None
    if name in gencode: return gencode[name]['id']
    base = name.split('.')[0]
    if base in gencode: return gencode[base]['id']
    if name in gencode_by_name: return gencode_by_name[name]['id']
    return None

def resolve_fusion_from_gffcompare(tid):
    gc = gffcompare.get(tid, {})
    ref_gid = gc.get('ref_gene_id', '')
    if ref_gid and ref_gid != '-':
        gid = to_gencode_id(ref_gid)
        if gid: return gid
    gname = gc.get('gene_name', '')
    if gname and gname != '-':
        gid = to_gencode_id(gname)
        if gid: return gid
    return None

for tid, t in transcripts.items():
    sq = sqanti.get(tid, {})
    cat = sq.get('structural_category', '').strip().lower()
    assoc = sq.get('associated_gene', '').strip().split(',')[0]

    if cat in ISOFORM_CATS:
        t['final_gene_id'] = to_gencode_id(assoc) or t['gene_id']
    elif cat == 'fusion':
        resolved = resolve_fusion_from_gffcompare(tid)
        if resolved:
            t['final_gene_id'] = resolved
        else:
            first_assoc_gene = assoc.split('_')[0]
            t['final_gene_id'] = to_gencode_id(first_assoc_gene) or t['gene_id']
    elif cat == 'genic':
        host_id = to_gencode_id(assoc)
        if is_protein_coding(t):
            t['final_gene_id'] = host_id or t['gene_id']
        elif host_id:
            t['final_gene_id'] = f"{t['gene_id']}_OT_{host_id.split('.')[0]}"
        else:
            t['final_gene_id'] = t['gene_id']
    else:
        t['final_gene_id'] = t['gene_id']

genes = collections.defaultdict(list)
for tid, t in transcripts.items():
    genes[t['final_gene_id']].append(tid)

gene_coords = {}
for gid, tids in genes.items():
    lo, hi = float('inf'), 0
    chrom = strand = ''
    for tid in tids:
        t = transcripts[tid]
        if not t['exons']: continue
        chrom, strand = t['chrom'], t['strand']
        es = sorted(t['exons'])
        lo = min(lo, es[0][0]); hi = max(hi, es[-1][1])
    if lo != float('inf'):
        gene_coords[gid] = {'chrom': chrom, 'strand': strand, 'start': lo, 'end': hi}

def locus_category(tids):
    cats = {sqanti.get(tid, {}).get('structural_category', '').strip().lower()
            for tid in tids}
    for pri in LOCUS_PRIORITY:
        if pri in cats: return pri
    return 'intergenic'

def sqanti_partner(tids):
    for tid in tids:
        assoc = sqanti.get(tid, {}).get('associated_gene', '').strip()
        if not assoc or assoc == '.' or assoc.startswith('novelGene'): continue
        first = assoc.split(',')[0]
        if first.startswith('ENSG'):
            info = gencode.get(first) or gencode.get(first.split('.')[0])
        else:
            info = gencode_by_name.get(first)
        if info and info['type'] == 'protein_coding' and info['name']:
            return info['name']
    return None

def find_overlap_partner(chrom, s, e, want_strand):
    best_ov, best = 0, None
    for g in pc_genes_by_chrom.get(chrom, []):
        if g['strand'] != want_strand: continue
        if g['end'] < s or g['start'] > e: continue
        ov = min(e, g['end']) - max(s, g['start']) + 1
        if ov > best_ov: best_ov, best = ov, g['name']
    return best

def find_divergent_partner(chrom, s, e, strand):
    tss = s if strand == '+' else e
    best_d, best = DIVERGENT_MAX_DIST + 1, None
    for g in pc_genes_by_chrom.get(chrom, []):
        if g['strand'] == strand: continue
        if strand == '+':
            if g['end'] >= s: continue
            d = tss - g['end']
        else:
            if g['start'] <= e: continue
            d = g['start'] - tss
        if 0 < d <= DIVERGENT_MAX_DIST and d < best_d:
            best_d, best = d, g['name']
    return best

linc_n = 1
as_n = collections.defaultdict(int)
it_n = collections.defaultdict(int)
ot_n = collections.defaultdict(int)
dt_n = collections.defaultdict(int)
final_info = {}

for gid in sorted(genes.keys()):
    tids = genes[gid]
    is_coding = any(is_protein_coding(transcripts[tid]) for tid in tids)

    if gid.startswith('ENSG'):
        info = gencode.get(gid) or gencode.get(gid.split('.')[0])
        final_info[gid] = ({'id': info['id'], 'name': info['name'], 'type': info['type']}
                           if info else
                           {'id': gid, 'name': gid,
                            'type': 'protein_coding' if is_coding else 'lncRNA'})
        continue

    if is_coding:
        final_info[gid] = {'id': gid, 'name': gid, 'type': 'protein_coding'}
        continue

    if gid not in gene_coords:
        final_info[gid] = {'id': gid, 'name': f'CRCN_{gid}', 'type': 'lncRNA'}
        continue

    gc = gene_coords[gid]
    cat = locus_category(tids)
    name = None

    if cat == 'antisense':
        opp = '-' if gc['strand'] == '+' else '+'
        p = sqanti_partner(tids) or find_overlap_partner(gc['chrom'], gc['start'], gc['end'], opp)
        if p:
            as_n[p] += 1
            name = f'CRCN_{p}-AS{as_n[p]}'
    elif cat == 'genic':
        p = sqanti_partner(tids) or find_overlap_partner(gc['chrom'], gc['start'], gc['end'], gc['strand'])
        if p:
            ot_n[p] += 1
            name = f'CRCN_{p}-OT{ot_n[p]}'
    elif cat == 'genic_intron':
        p = sqanti_partner(tids) or find_overlap_partner(gc['chrom'], gc['start'], gc['end'], gc['strand'])
        if p:
            it_n[p] += 1
            name = f'CRCN_{p}-IT{it_n[p]}'
    elif cat == 'intergenic':
        p = find_divergent_partner(gc['chrom'], gc['start'], gc['end'], gc['strand'])
        if p:
            dt_n[p] += 1
            suf = '' if dt_n[p] == 1 else str(dt_n[p])
            name = f'CRCN_{p}-DT{suf}'

    if name is None:
        name = f'CRCN_LINC{str(linc_n).zfill(5)}'
        linc_n += 1

    final_info[gid] = {'id': gid, 'name': name, 'type': 'lncRNA'}

def find_exon_number(exon_map, seg_start, seg_end):
    for es, ee, num in exon_map:
        if es <= seg_start and seg_end <= ee:
            return num
    for es, ee, num in exon_map:
        if seg_start <= ee and seg_end >= es:
            return num
    return exon_map[0][2] if exon_map else 1

out_lines = []
def get_sort_key(gid):
    gc = gene_coords.get(gid)
    return (gc['chrom'], gc['start']) if gc else ('', 0)

out_lines = []
for gid in sorted(genes.keys(), key=get_sort_key):
    tids = genes[gid]
    valid = sorted(tids, key=lambda tid: min(e[0] for e in transcripts[tid]['exons']) if transcripts[tid]['exons'] else 0)
    
    if not valid: continue
    info = final_info.get(gid)
    if not info: continue

    if not gid.startswith('ENSG') and gid in gene_coords:
        gc = gene_coords[gid]
        attrs = f"ID={info['id']};gene_id={info['id']};gene_type={info['type']};gene_name={info['name']}"
        out_lines.append(f"{gc['chrom']}\tCUSTOM\tgene\t{gc['start']}\t{gc['end']}\t.\t{gc['strand']}\t.\t{attrs}")

    for idx, tid in enumerate(valid, 1):
        t = transcripts[tid]
        t_type = classify_transcript(t)
        exons = sorted(set(t['exons']))
        cds = sorted(set(t['cds']))
        strand, chrom = t['strand'], t['chrom']
        t_start, t_end = exons[0][0], exons[-1][1]
        t_name = f"{info['name']}-{idx}"
        base = f"gene_id={info['id']};transcript_id={tid};gene_type={info['type']};gene_name={info['name']};transcript_type={t_type};transcript_name={t_name}"
        out_lines.append(f"{chrom}\tCUSTOM\ttranscript\t{t_start}\t{t_end}\t.\t{strand}\t.\tID={tid};Parent={info['id']};{base}")

        exon_map = []
        for i, ex in enumerate(exons):
            num = i + 1 if strand == '+' else len(exons) - i
            exon_map.append((ex[0], ex[1], num))
            out_lines.append(f"{chrom}\tCUSTOM\texon\t{ex[0]}\t{ex[1]}\t.\t{strand}\t.\tID=exon:{tid}:{num};Parent={tid};{base};exon_number={num}")

        if cds and t_type in ('protein_coding', 'nonsense_mediated_decay'):
            cmin, cmax = cds[0][0], cds[-1][1]
            ordered = cds if strand == '+' else list(reversed(cds))
            
            for c in ordered:
                cs, ce, orig_phase = c[0], c[1], c[2]
                enum = find_exon_number(exon_map, cs, ce)
                out_lines.append(f"{chrom}\tCUSTOM\tCDS\t{cs}\t{ce}\t.\t{strand}\t{orig_phase}\tID=CDS:{tid};Parent={tid};{base};exon_number={enum}")

            utr5, utr3 = [], []
            for ex in exons:
                if strand == '+':
                    if ex[1] < cmin: utr5.append(ex)
                    elif ex[0] < cmin: utr5.append((ex[0], cmin - 1))
                    if ex[0] > cmax: utr3.append(ex)
                    elif ex[1] > cmax: utr3.append((cmax + 1, ex[1]))
                else:
                    if ex[0] > cmax: utr5.append(ex)
                    elif ex[1] > cmax: utr5.append((cmax + 1, ex[1]))
                    if ex[1] < cmin: utr3.append(ex)
                    elif ex[0] < cmin: utr3.append((ex[0], cmin - 1))
            
            for u in utr5:
                enum = find_exon_number(exon_map, u[0], u[1])
                out_lines.append(f"{chrom}\tCUSTOM\tfive_prime_UTR\t{u[0]}\t{u[1]}\t.\t{strand}\t.\tID=UTR5:{tid};Parent={tid};{base};exon_number={enum}")
            for u in utr3:
                enum = find_exon_number(exon_map, u[0], u[1])
                out_lines.append(f"{chrom}\tCUSTOM\tthree_prime_UTR\t{u[0]}\t{u[1]}\t.\t{strand}\t.\tID=UTR3:{tid};Parent={tid};{base};exon_number={enum}")

with open(out_file, 'w') as f:
    f.write('\n'.join(out_lines) + '\n')
