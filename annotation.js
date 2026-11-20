// Annotation reading, shared by the main page and the reads view. Plain script: defines globals.
//   GenBank: reference -> FASTA + features (main page); GenBank feature files (reads view)
//   BED / BED12, GFF3, GTF: annotation tracks in the reads view (parseFeatures)
// Feature records: { start, end (0-based half-open), strand, type, name, info, blocks?, thick?, color? }

// GFF types describing the whole sequence, not a feature: left out of the tracks
const FEAT_SKIP_TYPES = new Set(['region', 'chromosome', 'contig', 'supercontig', 'scaffold', 'source', 'databank_entry']);

// FASTA sequence lines of 80 characters
const wrap80 = (seq) => seq.replace(/(.{80})/g, '$1\n').replace(/\n$/, '');

// GenBank location -> strand + 0-based half-open segments. Handles complement(), join()/order(),
// partial ends (<, >) and single bases; segments on other records (ACC.1:10..20) and
// between-base sites (10^11) are skipped.
function parseGenbankLocation(loc) {
    const strand = /complement\(/.test(loc) ? '-' : '+';
    const clean = loc.replace(/[A-Za-z][\w.]*:[<>]?\d+(?:\.\.[<>]?\d+)?/g, '').replace(/\d+\^\d+/g, '');
    const segs = [...clean.matchAll(/<?(\d+)(?:\.\.>?(\d+))?/g)].map(m => [parseInt(m[1], 10) - 1, parseInt(m[2] || m[1], 10)]);
    return { strand, segs };
}

// One GenBank feature -> feature records for the reads view's annotation track:
// { start, end, strand, type, name, info, blocks? } (blocks = joined segments, as exons)
function genbankFeature(type, loc, q, L) {
    const { strand, segs } = parseGenbankLocation(loc);
    const parts = segs.map(([a, b]) => [Math.max(0, a), Math.min(L, b)]).filter(([a, b]) => b > a).sort((x, y) => x[0] - y[0]);
    if (!parts.length) return [];
    const base = {
        type, strand,
        name: q.gene || q.locus_tag || q.label || q.standard_name || q.product || type,
        info: q.product || q.note || ''
    };
    const covered = parts.reduce((t, [a, b]) => t + b - a, 0);
    const end = Math.max(...parts.map(p => p[1]));
    // a join across the origin of a circular sequence (e.g. join(9000..10000,1..200)): draw the pieces
    if (parts.length > 1 && end - parts[0][0] - covered > L / 2) return parts.map(([a, b]) => ({ ...base, start: a, end: b }));
    const f = { ...base, start: parts[0][0], end };
    if (parts.length > 1) {
        const blocks = [];
        for (const [a, b] of parts) {
            if (blocks.length && a <= blocks[blocks.length - 1]) blocks[blocks.length - 1] = Math.max(blocks[blocks.length - 1], b);
            else blocks.push(a, b);
        }
        if (blocks.length > 2) f.blocks = blocks;
    }
    return [f];
}

// GenBank (one or more records): name from LOCUS, sequence from ORIGIN, features from the
// FEATURES table (except 'source'). Returns { fasta ('' if no record has a sequence),
// annotation: Map record name -> features }.
function parseGenbank(text) {
    const out = [], annotation = new Map();
    let name = null, section = '', seq = [], raw = [], cur = null, qual = null, locusLen = 0;
    const endFeature = () => { if (cur) raw.push(cur); cur = null; qual = null; };
    const flush = () => {
        endFeature();
        if (name) {
            const s = seq.join('').toUpperCase();
            if (s) out.push(`>${name}\n${wrap80(s)}`);
            // features need the record length: the sequence, or the LOCUS line for annotation-only files
            const L = s.length || locusLen;
            const feats = !L ? [] : raw.filter(f => f.type !== 'source').flatMap(f => {
                const q = {};
                for (const [k, v] of f.q) if (!(k in q)) q[k] = v.replace(/^"|"$/g, '').replace(/""/g, '"');
                return genbankFeature(f.type, f.loc, q, L);
            });
            if (feats.length) annotation.set(name, feats);
        }
        name = null; section = ''; seq = []; raw = []; locusLen = 0;
    };
    for (const line of text.split(/\r?\n/)) {
        if (line.startsWith('LOCUS')) {
            flush();
            const t = line.split(/\s+/);
            name = t[1] || `record_${out.length + 1}`;
            locusLen = parseInt(t[2], 10) || 0; // "LOCUS  name  12345 bp ..."
            continue;
        }
        if (line.startsWith('//')) { flush(); continue; }
        if (/^[A-Z]/.test(line)) { // a new top-level keyword: FEATURES, ORIGIN, CONTIG, BASE COUNT, ...
            endFeature();
            section = line.startsWith('FEATURES') ? 'features' : line.startsWith('ORIGIN') ? 'origin' : '';
            continue;
        }
        if (section === 'origin') { seq.push(line.replace(/[^A-Za-z]/g, '')); continue; }
        if (section !== 'features') continue;
        const key = line.match(/^ {5}(\S+)\s+(.*)$/);
        if (key) { endFeature(); cur = { type: key[1], loc: key[2].trim(), q: [] }; continue; }
        if (!cur) continue;
        const body = line.trim();
        const qm = body.match(/^\/(\w+)(?:=(.*))?$/);
        if (qm) { qual = [qm[1], qm[2] ?? '']; if (qm[1] !== 'translation') cur.q.push(qual); continue; }
        if (qual) { if (qual[0] !== 'translation') qual[1] += ' ' + body; }  // continued qualifier value
        else cur.loc += body;                                                   // continued location
    }
    flush(); // last record without a closing //
    return { fasta: out.length ? out.join('\n') + '\n' : '', annotation };
}

// ---- BED / GFF3 / GTF ----
function gffAttrs(s) {
    const a = {};
    for (const part of s.split(';')) {
        const p = part.trim();
        if (!p) continue;
        const m = p.match(/^([^=\s]+)\s*=\s*(.*)$/) || p.match(/^(\S+)\s+"?([^"]*)"?$/); // GFF3 key=value, GTF key "value"
        if (!m) continue;
        try { a[m[1]] = decodeURIComponent(m[2]); } catch { a[m[1]] = m[2]; }
    }
    return a;
}

// GFF colour attribute as IGV reads it: color=#rrggbb, #rgb or r,g,b
function gffColor(v) {
    if (!v) return undefined;
    v = v.trim();
    if (/^#([0-9a-f]{3}|[0-9a-f]{6})$/i.test(v)) return v;
    if (/^\d{1,3},\d{1,3},\d{1,3}$/.test(v)) return `rgb(${v})`;
    return undefined;
}

// Features of one contig from a BED, GFF3, GTF or GenBank file (by content, else file name).
// Returns { features (this contig only, 0-based half-open, clipped to contigLength), others: other seq names }
function parseFeatures(text, filename, contig, contigLength) {
    if (text.trimStart().startsWith('LOCUS')) { // GenBank: its FEATURES tables
        const gb = parseGenbank(text);
        return { features: gb.annotation.get(contig) || [], others: [...gb.annotation.keys()].filter(k => k !== contig) };
    }
    const gff = /\.(gff3?|gtf)(\.gz)?$/i.test(filename) || /^##gff-version/m.test(text.slice(0, 2000));
    const L = contigLength;
    const out = [], seqs = new Set();
    for (const line of text.split(/\r?\n/)) {
        if (!line.trim()) continue;
        if (line.startsWith('##FASTA') || line[0] === '>') break; // GFF3 with embedded sequences
        if (line[0] === '#' || /^(track|browser)\b/.test(line)) continue;
        let f = line.split('\t');
        if (f.length < 3) f = line.trim().split(/\s+/);
        if (f.length < 3) continue;
        seqs.add(f[0]);
        if (f[0] !== contig) continue;
        let feat;
        if (gff) {
            if (f.length < 9 || FEAT_SKIP_TYPES.has(f[2])) continue;
            const a = gffAttrs(f[8]);
            const gtf = !('ID' in a || 'Parent' in a) && 'transcript_id' in a;
            feat = {
                start: parseInt(f[3], 10) - 1, end: parseInt(f[4], 10), strand: f[6], type: f[2],
                name: a.Name || a.gene || a.gene_name || a.locus_tag || a.ID || a.gene_id || a.transcript_id || f[2],
                info: a.product || a.Note || a.description || '',
                color: gffColor(a.color || a.colour),
                // for exon grouping: GFF3 ID/Parent, GTF transcript_id
                id: gtf ? (f[2] === 'transcript' ? a.transcript_id : null) : a.ID,
                parents: gtf ? (f[2] === 'transcript' ? [] : [a.transcript_id]) : (a.Parent ? a.Parent.split(',') : []),
                gtf
            };
        } else {
            feat = { start: parseInt(f[1], 10), end: parseInt(f[2], 10), strand: f[5] || '.', type: 'BED', name: f[3] || '', info: '' };
            if (f[8] && /^\d+,\d+,\d+$/.test(f[8]) && f[8] !== '0,0,0') feat.color = `rgb(${f[8]})`;
            const n = parseInt(f[9], 10);
            if (n > 1 && f[10] && f[11]) { // BED12 blocks (exons)
                const sizes = f[10].split(',').map(Number), starts = f[11].split(',').map(Number);
                feat.blocks = [];
                for (let i = 0; i < n; i++) feat.blocks.push(feat.start + starts[i], feat.start + starts[i] + sizes[i]);
            }
            // BED8+: thickStart/thickEnd = coding part (thickStart == thickEnd: non-coding, all thin)
            const ts = parseInt(f[6], 10), te = parseInt(f[7], 10);
            if (f.length >= 8 && te >= ts) {
                feat.thick = [ts, te];
                if (!feat.blocks) feat.blocks = [feat.start, feat.end]; // single exon: drawn by the block path
            }
        }
        if (!(feat.end > feat.start)) continue;
        feat.start = Math.max(0, feat.start); feat.end = Math.min(L, feat.end);
        if (feat.end > feat.start) out.push(feat);
    }
    seqs.delete(contig);
    return { features: gff ? groupExons(out) : out, others: [...seqs] };
}

// Join transcript parts into their parent (GFF3 Parent=, GTF transcript_id): exons become the
// parent's blocks (or several CDS pieces when there are no exons), the CDS extent becomes its
// thick (coding) part, and the parts themselves are dropped from the track.
// A gene with a single CDS (typical bacterial annotation) is left as it is.
const PART_TYPES = new Set(['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'UTR', '5UTR', '3UTR',
    'start_codon', 'stop_codon', 'intron']);
function groupExons(feats) {
    const byId = new Map();
    for (const f of feats) if (f.id && !PART_TYPES.has(f.type) && !byId.has(f.id)) byId.set(f.id, f);
    const kids = new Map(); // parent id -> { exon: [], CDS: [], all: [] }
    for (const f of feats) {
        if (!PART_TYPES.has(f.type)) continue;
        for (const p of f.parents) {
            if (!kids.has(p)) kids.set(p, { exon: [], CDS: [], all: [] });
            const k = kids.get(p);
            if (f.type === 'exon') k.exon.push(f); else if (f.type === 'CDS') k.CDS.push(f);
            k.all.push(f);
        }
    }
    const drop = new Set(), added = [];
    for (const [pid, k] of kids) {
        const parts = k.exon.length ? k.exon : k.CDS.length > 1 ? k.CDS : null;
        if (!parts) continue;
        let parent = byId.get(pid);
        if (!parent) {
            if (!parts[0].gtf) continue; // GFF3 part whose parent is missing: keep it as is
            // GTF without transcript lines: make the transcript from its parts
            parent = { start: Infinity, end: 0, strand: parts[0].strand, type: 'transcript', name: pid, info: '', color: parts[0].color };
            byId.set(pid, parent); added.push(parent);
        }
        // merge the parts into sorted, non-overlapping blocks
        const iv = parts.map(p => [p.start, p.end]).sort((a, b) => a[0] - b[0]);
        const blocks = [];
        for (const [s, e] of iv) {
            if (blocks.length && s <= blocks[blocks.length - 1]) blocks[blocks.length - 1] = Math.max(blocks[blocks.length - 1], e);
            else blocks.push(s, e);
        }
        parent.blocks = blocks;
        parent.start = Math.min(parent.start, blocks[0]);
        parent.end = Math.max(parent.end, blocks[blocks.length - 1]);
        if (k.CDS.length) {
            parent.thick = [Math.min(...k.CDS.map(c => c.start)), Math.max(...k.CDS.map(c => c.end))];
            if (!parent.info) parent.info = k.CDS.find(c => c.info)?.info || '';
        }
        for (const p of k.all) drop.add(p);
    }
    // features without a colour take their parent's (e.g. an mRNA under a coloured gene)
    for (const f of feats) {
        const p = !f.color && f.parents?.length && byId.get(f.parents[0]);
        if (p?.color) f.color = p.color;
    }
    return feats.filter(f => !drop.has(f)).concat(added);
}
