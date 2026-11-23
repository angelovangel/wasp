// Per-base counts and variant calls, shared by the main page (script.js: the contig-wide variant
// scan) and the reads view (reads.html: coverage colours, fallback scan). Plain script: defines
// globals. Needs MIN_GAP_BP from common.js.

const BASES = ['A', 'C', 'G', 'T'];
const ALLELE_FRAC = 0.2;   // lowest allele frequency the variant scan keeps (the cutoff menu filters these)
const MIN_CLIP_BP = 20;    // soft clips at least this long count as "clipped" for the pile-up marker
const MIN_VAR_DEPTH = 5;   // variant calls: minimum depth
const MIN_ALT_READS = 3;   // variant calls / clip marker: minimum supporting reads

// Per-base counts over [s, e) from all loaded reads (including any hidden from the pile):
//   depth (reads spanning, excl. deletions), del, A/C/G/T, consensus base per position.
// Reference bases come from refSeq (the contig sequence); reads only carry their mismatches.
function countBases(reads, refSeq, s, e) {
    s = Math.max(0, Math.floor(s)); e = Math.min(refSeq.length, Math.ceil(e));
    const n = e - s;
    // Difference arrays (+1 at start, -1 at end) so each read costs O(its features), not O(length)
    const span = new Int32Array(n + 1);     // reads whose aligned blocks cover the base
    const smallDel = new Int32Array(n + 1); // short deletions (inside a block)
    const longDel = new Int32Array(n + 1);  // long deletions (gap between blocks)
    const mm = { A: new Int32Array(n), C: new Int32Array(n), G: new Int32Array(n), T: new Int32Array(n) };
    const addRange = (arr, a, b) => {
        a = Math.max(a, s); b = Math.min(b, e);
        if (a < b) { arr[a - s]++; arr[b - s]--; }
    };
    const clip = new Int32Array(n);   // reads with a long soft clip starting/ending here
    const ins = new Int32Array(n);    // reads with an insertion right after this base
    const insLens = new Map();        // index -> insertion lengths there
    for (const r of reads) {
        if (r.end <= s || r.start >= e) continue;
        if (r.clipL >= MIN_CLIP_BP && r.start >= s && r.start < e) clip[r.start - s]++;
        if (r.clipR >= MIN_CLIP_BP && r.end - 1 >= s && r.end - 1 < e) clip[r.end - 1 - s]++;
        for (let b = 0; b < r.blocks.length; b += 2) addRange(span, r.blocks[b], r.blocks[b + 1]);
        for (let d = 0; d < r.dels.length; d += 2) {
            const p = r.dels[d], len = r.dels[d + 1];
            addRange(len < MIN_GAP_BP ? smallDel : longDel, p, p + len);
        }
        for (let m = 0; m < r.mmPos.length; m++) {
            const p = r.mmPos[m], b = r.mmBase[m];
            if (p >= s && p < e && mm[b]) mm[b][p - s]++;
        }
        for (let k = 0; k < r.ins.length; k += 2) {
            const a = r.ins[k] - 1; // r.ins holds the position the insertion precedes
            if (a < s || a >= e) continue;
            ins[a - s]++;
            const l = insLens.get(a - s);
            if (l) l.push(r.ins[k + 1]); else insLens.set(a - s, [r.ins[k + 1]]);
        }
    }
    const depth = new Int32Array(n), del = new Int32Array(n);
    const cnt = { A: new Int32Array(n), C: new Int32Array(n), G: new Int32Array(n), T: new Int32Array(n) };
    const ref = new Array(n), cons = new Array(n), alt = new Float32Array(n);
    let sp = 0, sd = 0, ld = 0;
    for (let i = 0; i < n; i++) {
        sp += span[i]; sd += smallDel[i]; ld += longDel[i];
        const aligned = Math.max(0, sp - sd);   // reads with a base here
        const dels = sd + ld;                   // reads with a deletion here
        const rb = refSeq[s + i];
        let mmSum = 0;
        for (const b of BASES) { cnt[b][i] = mm[b][i]; mmSum += mm[b][i]; }
        if (cnt[rb]) cnt[rb][i] += Math.max(0, aligned - mmSum); // the rest match the reference
        depth[i] = aligned; del[i] = dels; ref[i] = rb;
        // consensus = most common of A/C/G/T/deletion
        let best = '', bestN = 0;
        for (const b of BASES) if (cnt[b][i] > bestN) { best = b; bestN = cnt[b][i]; }
        if (dels > bestN) best = '-';
        cons[i] = best;
        const tot = aligned + dels;
        alt[i] = tot > 0 && cnt[rb] ? 1 - cnt[rb][i] / tot : 0;
    }
    return { s, e, depth, del, cnt, ref, cons, alt, clip, ins, insLens };
}

// The one rule for calling a variant at index i of per-base counts `c` (from computeCounts),
// used by the table, the coloured coverage bars and the underlined reference bases.
// Returns { alt, n, af } for the most common non-reference allele ('-' = deletion) and
// { insN, insAf } for an insertion after the base; either is null if it doesn't pass
// depth >= MIN_VAR_DEPTH and >= MIN_ALT_READS supporting reads. `minFrac` is the AF cutoff.
function callAt(c, i, minFrac) {
    const tot = c.depth[i] + c.del[i];
    if (tot < MIN_VAR_DEPTH) return { snv: null, ins: null, tot };
    const rb = c.ref[i];
    let alt = null, altN = 0;
    for (const b of BASES) if (b !== rb && c.cnt[b][i] > altN) { alt = b; altN = c.cnt[b][i]; }
    if (c.del[i] > altN) { alt = '-'; altN = c.del[i]; }
    const insN = c.ins[i];
    return {
        tot,
        snv: alt && altN >= MIN_ALT_READS && altN / tot >= minFrac ? { alt, n: altN, af: altN / tot } : null,
        ins: insN >= MIN_ALT_READS && insN / tot >= minFrac ? { n: insN, af: insN / tot } : null,
    };
}

// Candidate variants (AF >= ALLELE_FRAC; the cutoff menu filters further) appended to `out`.
function collectVariants(c, out) {
    for (let i = 0; i < c.depth.length; i++) {
        const { snv, ins, tot } = callAt(c, i, ALLELE_FRAC);
        if (!snv && !ins) continue;
        const row = { pos: c.s + i, ref: c.ref[i], depth: tot,
            A: c.cnt.A[i], C: c.cnt.C[i], G: c.cnt.G[i], T: c.cnt.T[i], del: c.del[i], ins: c.ins[i] };
        if (snv) out.push({ ...row, ...snv });
        // insertion after this base: its own row; length = median of the reads' insertion lengths
        if (ins) {
            const lens = c.insLens.get(i).sort((a, b) => a - b);
            out.push({ ...row, alt: 'ins', insLen: lens[lens.length >> 1], ...ins });
        }
    }
}
