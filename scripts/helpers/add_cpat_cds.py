import sys

cpat_file = sys.argv[1]
gtf_file = sys.argv[2]
out_file = sys.argv[3]
threshold = 0.364

cpat_orfs = {}
with open(cpat_file, 'r') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 11:
            continue
        tid = parts[0]
        orf_strand = parts[3]
        if orf_strand != '+':
            continue
        start = int(parts[5])
        end = int(parts[6])
        prob = float(parts[10])
        if prob >= threshold:
            cpat_orfs[tid] = (start, end)
            
gtf_data = {}
new_gtf_lines = []
with open(gtf_file, 'r') as f:
    for line in f:
        new_gtf_lines.append(line.strip())
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        if parts[2] == 'exon':
            tid = ""
            attrs = parts[8].split(';')
            for attr in attrs:
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    tid = attr.split('"')[1]
                    break
            if tid in cpat_orfs:
                if tid not in gtf_data:
                    gtf_data[tid] = {'strand': parts[6], 'exons': [], 'chr': parts[0], 'attrs': parts[8]}
                gtf_data[tid]['exons'].append((int(parts[3]), int(parts[4])))

for tid, data in gtf_data.items():
    start_t, end_t = cpat_orfs[tid]
    exons = sorted(data['exons'], key=lambda x: x[0])
    strand = data['strand']
    
    if strand == '-':
        exons.reverse()

    current_t_len = 0
    cds_genomic = []
    for ex_start, ex_end in exons:
        ex_len = ex_end - ex_start + 1
        if current_t_len + ex_len >= start_t and current_t_len < end_t:
            overlap_start_t = max(start_t, current_t_len + 1)
            overlap_end_t = min(end_t, current_t_len + ex_len)
            
            offset_start = overlap_start_t - current_t_len - 1
            offset_end = overlap_end_t - current_t_len - 1
            
            if strand == '+':
                g_start = ex_start + offset_start
                g_end = ex_start + offset_end
            else:
                g_start = ex_end - offset_end
                g_end = ex_end - offset_start
            
            cds_genomic.append((min(g_start, g_end), max(g_start, g_end)))
        current_t_len += ex_len

    if strand == '+':
        cds_genomic.sort(key=lambda x: x[0])
    else:
        cds_genomic.sort(key=lambda x: x[0], reverse=True)

    phase = 0
    for g_start, g_end in cds_genomic:
        new_gtf_lines.append(f"{data['chr']}\tCPAT\tCDS\t{g_start}\t{g_end}\t.\t{strand}\t{phase}\t{data['attrs']}")
        length = g_end - g_start + 1
        phase = (3 - ((length - phase) % 3)) % 3

with open(out_file, 'w') as f:
    for line in new_gtf_lines:
        f.write(line + '\n')
