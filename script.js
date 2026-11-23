const CLI = await new Aioli([
    "samtools/1.21", 
    "minimap2/2.22"
], {debug: true});

// --- PARALLEL ALIGNMENT: extra Aioli instances (each its own Web Worker, so its own core) align
// read batches alongside CLI. Each holds a copy of the minimap2 index and one batch (~0.5-1 GB at
// most), so their number is capped by cores, and by device memory when the browser reports it.
const MAX_ALIGNERS = 8;          // workers aligning at once, CLI included
const GB_PER_ALIGNER = 1;        // device memory set aside per worker
const PARALLEL_BATCH_BASES = 50e6; // largest batch when several workers align at once…
const MIN_BATCH_BASES = 5e6;       // …and the smallest (each batch re-loads the index)
const extraAligners = [];        // created on first use and kept for later runs
function alignerCount() {
    const cores = navigator.hardwareConcurrency || 2;
    // Chrome only, rounded and capped at 8 (absent on plain http): a missing value doesn't limit
    const memGB = navigator.deviceMemory || 8;
    return Math.max(1, Math.min(MAX_ALIGNERS, cores - 1, Math.floor(memGB / GB_PER_ALIGNER)));
}
// CLI plus enough extra workers to make `n`; the minimap2 index `mmi` (a path in CLI's file
// system) is copied into each extra one.
async function getAligners(n, mmi) {
    while (extraAligners.length < n - 1) extraAligners.push(await new Aioli(["samtools/1.21", "minimap2/2.22"]));
    const pool = [CLI, ...extraAligners.slice(0, n - 1)];
    if (pool.length > 1) {
        const bytes = await CLI.fs.readFile(mmi);
        for (const a of pool.slice(1)) await a.fs.writeFile(mmi, bytes);
    }
    return pool;
}
async function releaseAligners(pool, mmi) {
    for (const a of pool.slice(1)) { try { await a.fs.unlink(mmi); } catch {} }
}

const NBINS = 300; // number of windows in each contig's depth profile
const DEPTH_CHUNK_BP = 2e6; // samtools depth covers at most this many bases per call, so its output stays far below browser string limits
const DEPTH_LONG_BP = 1e6;  // contigs at least this long are read region by region
const MIN_BIN_BP = 50;      // short contigs get fewer, wider bins (not NBINS bins of a few bp)
const MMI_MAX_REF_BYTES = 50e6; // up to this reference size, a minimap2 index (.mmi) is built once and reused for every batch
const SPARKLINE_COLOR = '#4c9aff';
const MIN_LINE_COLOR = '#c2255c'; // lowest depth per bin

// **NEW GLOBAL VARIABLE** to store max length for proportional scaling
let maxContigLength = 1; // Initialize to 1 to avoid division by zero

// One entry per sample (a FASTQ file or R1/R2 pair, or a BAM):
// { name, bam (path in the in-browser file system), contigs, summary, sampling?, bamFile?, paired, fromBam }
let samples = [];
let referenceFile = null; // the reference fasta File, read again by the reads view
let referenceAnnotation = null; // GenBank reference: { name, byContig: Map contig -> features } for the reads view
const contigSeqCache = new Map();

// Function to reload the page
function resetPage() {
    window.location.reload();
}

// Aggressive cleaning function to ensure keys match
function cleanContigName(name) {
    // Trim whitespace
    return name.trim();
}

// formatBp, niceStep, escapeHtml, downloadBlob, MIN_GAP_BP/MIN_INS_BP, attachQuickTips: common.js

// Thin contig line with tick marks at nice positions
function buildAxisHtml(length, widthPct) {
    const step = Math.max(1, niceStep(length, Math.max(2, Math.round(widthPct / 20))));
    let ticks = '';
    // nice ticks, skipping any too close to the end so they don't collide with the length label
    for (let pos = 0; pos <= length - 0.75 * step; pos += step) {
        const pct = (pos / length) * 100;
        // keep labels inside the line: left-align at start, centre elsewhere
        const shift = pct < 5 ? '0' : '-50%';
        ticks += `<span class="axis-tick" style="left:${pct}%"></span>` +
                 `<span class="axis-label" style="left:${pct}%;transform:translateX(${shift})">${formatBp(pos)}</span>`;
    }
    // final tick: exact contig length, right-aligned to the end of the line
    ticks += `<span class="axis-tick" style="left:100%"></span>` +
             `<span class="axis-label axis-end" style="left:100%;transform:translateX(-100%)">${formatBp(length)}</span>`;
    return `<div class="contig-axis">${ticks}</div>`;
}

// Sparkline + contig line, both scaled to contig length relative to the longest contig
function renderSparklines(root) {
    $(root).find('.profile-track').each(function() {
        const $track = $(this);
        const $spark = $track.find('.sparkline-container');
        const depths = $spark.data('depths');
        const mins = $spark.data('mins');
        const contigLength = $track.data('length');

        const proportionalWidth = Math.max(20, (contigLength / maxContigLength) * 100);
        $track.css('width', `${proportionalWidth}%`);
        $track.append(buildAxisHtml(contigLength, proportionalWidth));

        if (depths && depths.length > 0) {
            const binSize = $track.data('bin') || Math.max(1, Math.ceil(contigLength / NBINS));
            const yMax = Math.max(1, ...depths);
            $spark.sparkline(depths, {
                type: 'line',
                width: '100%',
                height: '50px',
                lineColor: SPARKLINE_COLOR,
                fillColor: SPARKLINE_COLOR + '33',
                spotColor: false,
                minSpotColor: false,
                maxSpotColor: false,
                chartRangeMin: 0,
                chartRangeMax: yMax,
                highlightLineColor: '#6c757d',
                highlightSpotColor: '#6c757d',
                // tooltip shows the window's position range and its mean depth
                tooltipFormatter: (sp, opts, fields) => {
                    const st = fields.x * binSize;
                    const en = Math.min(st + binSize, contigLength);
                    return `${formatBp(st + 1)}–${formatBp(en)}<br>Mean depth: ${fields.y}`;
                }
            });
            // Lowest depth within each bin, drawn as a thin line over the mean
            // so narrow dips (deletions, gaps) aren't averaged away
            if (mins && mins.length === depths.length) {
                $spark.sparkline(mins, {
                    composite: true,
                    type: 'line',
                    lineColor: MIN_LINE_COLOR,
                    lineWidth: 1,
                    fillColor: false,
                    spotColor: false,
                    minSpotColor: false,
                    maxSpotColor: false,
                    highlightSpotColor: MIN_LINE_COLOR,
                    chartRangeMin: 0,
                    chartRangeMax: yMax,
                    tooltipFormatter: (sp, opts, fields) => `<br>Lowest: ${fields.y}`
                });
            }
        } else {
            $spark.text('No depth data');
        }
    });
}

// --- DEPTH PROFILES: mean + lowest depth per bin ---
// `samtools depth` prints one line per covered base, so its output is read in pieces that each
// cover at most DEPTH_CHUNK_BP bases:
//   - long contigs: one call per DEPTH_CHUNK_BP region (-r name:start-end, uses the BAM index)
//   - short contigs, if they add up to <= DEPTH_CHUNK_BP: one call for all of them
//   - short contigs adding up to more (e.g. draft assemblies with many contigs): samtools bedcov
//     over each contig's bins, which gives the mean only (no per-bin minimum)
// Without -a only covered positions are printed, so bins with missing positions have min 0.
// Deletions are not counted as coverage (samtools depth default), so they show as dips.

// Run `samtools <cmd> -o FILE ...` and return FILE as text. Much faster than taking its stdout,
// which reaches JS line by line: use it for commands with big output. If samtools wrote no file
// (an error), returns what exec printed instead, as exec alone would have.
let outFileN = 0;
async function samtoolsText(args) {
    const out = `out${++outFileN}.txt`;
    const log = await CLI.exec("samtools", [args[0], "-o", out, ...args.slice(1)]);
    try { return new TextDecoder().decode(await CLI.fs.readFile(out)); }
    catch { return log; }
    finally { try { await CLI.fs.unlink(out); } catch {} }
}

// Add one `samtools depth` output into the per-contig bin sums (walks it without splitting it)
function addDepthOutput(out, acc) {
    let cur = null, curName = null;
    for (let i = 0; i < out.length;) {
        let nl = out.indexOf('\n', i);
        if (nl === -1) nl = out.length;
        const t1 = out.indexOf('\t', i), t2 = t1 === -1 ? -1 : out.indexOf('\t', t1 + 1);
        if (t1 !== -1 && t2 !== -1 && t2 < nl) {
            const name = out.slice(i, t1);
            if (name !== curName) { curName = name; cur = acc.get(cleanContigName(name)) || null; }
            if (cur) {
                const pos = parseInt(out.slice(t1 + 1, t2), 10) - 1;
                const d = parseInt(out.slice(t2 + 1, nl), 10);
                const b = Math.floor(pos / cur.c.binSize);
                if (b >= 0 && b < cur.sum.length && !isNaN(d)) {
                    cur.sum[b] += d; cur.n[b]++;
                    if (d < cur.min[b]) cur.min[b] = d;
                }
            }
        }
        i = nl + 1;
    }
}

async function depthProfiles(contigs, bam, progress = () => {}) {
    const round1 = x => Math.round(x * 10) / 10;
    const profiles = {}, minProfiles = {};
    const long = contigs.filter(c => c.length >= DEPTH_LONG_BP);
    const short = contigs.filter(c => !long.includes(c));
    const shortBp = short.reduce((t, c) => t + c.length, 0);
    const viaBedcov = shortBp > DEPTH_CHUNK_BP ? short : [];
    const viaDepth = viaBedcov.length ? long : contigs;

    const acc = new Map(viaDepth.map(c => {
        const nb = Math.ceil(c.length / c.binSize);
        return [c.name, { c, sum: new Float64Array(nb), n: new Int32Array(nb), min: new Float64Array(nb).fill(Infinity) }];
    }));
    const totalBp = contigs.reduce((t, c) => t + c.length, 0) || 1;
    let doneBp = 0;

    // short contigs together in one call (all of them, or only these via a BED file)
    if (!viaBedcov.length && short.length) {
        if (long.length) {
            await CLI.fs.writeFile("short.bed", short.map(c => `${c.name}\t0\t${c.length}`).join('\n') + '\n');
            addDepthOutput(await samtoolsText(["depth", "-d", "0", "-b", "short.bed", bam]), acc);
            await CLI.fs.unlink("short.bed");
        } else {
            addDepthOutput(await samtoolsText(["depth", "-d", "0", bam]), acc);
        }
        doneBp += shortBp;
        progress(doneBp / totalBp);
    }
    // long contigs region by region
    for (const c of long) {
        for (let st = 0; st < c.length; st += DEPTH_CHUNK_BP) {
            const en = Math.min(c.length, st + DEPTH_CHUNK_BP);
            addDepthOutput(await samtoolsText(["depth", "-d", "0", "-r", regionOf(c.name, st, en), bam]), acc);
            doneBp += en - st;
            progress(doneBp / totalBp);
        }
    }
    for (const { c, sum, n, min } of acc.values()) {
        const mean = [], low = [];
        for (let b = 0; b < sum.length; b++) {
            const binLen = Math.min(c.binSize, c.length - b * c.binSize);
            mean.push(round1(sum[b] / binLen));
            low.push(n[b] < binLen ? 0 : round1(min[b])); // unreported positions have depth 0
        }
        profiles[c.name] = mean;
        minProfiles[c.name] = low;
    }

    // many short contigs: bedcov over their bins (mean only). -j: don't count deletions as coverage.
    if (viaBedcov.length) {
        const bed = [];
        for (const c of viaBedcov) {
            for (let st = 0; st < c.length; st += c.binSize) bed.push(`${c.name}\t${st}\t${Math.min(st + c.binSize, c.length)}`);
        }
        await CLI.fs.writeFile("windows.bed", bed.join('\n') + '\n');
        let out = await CLI.exec("samtools", ["bedcov", "-j", "windows.bed", bam]);
        if (!/^\S+\t\d+\t\d+\t\d+/.test(out)) {
            // older samtools without -j: fall back to counting deletions as covered
            console.warn("samtools bedcov -j failed, retrying without -j:", out.slice(0, 200));
            out = await CLI.exec("samtools", ["bedcov", "windows.bed", bam]);
        }
        await CLI.fs.unlink("windows.bed");
        for (const line of out.trim().split('\n')) {
            const f = line.split('\t');
            if (f.length < 4) continue;
            const len = parseInt(f[2], 10) - parseInt(f[1], 10);
            const mean = len > 0 ? parseFloat(f[3]) / len : 0;
            (profiles[cleanContigName(f[0])] ??= []).push(round1(mean));
        }
        progress(1);
    }
    return { profiles, minProfiles };
}

// Progress bar: the wasm tools don't report internal progress, so the bar
// advances per pipeline step and an elapsed timer shows it's still working.
let progressTimer = null;

function showProgress(outputDiv) {
    outputDiv.innerHTML = `
        <div class="run-progress">
            <div class="run-progress-text">
                <span id="progress-label" style="white-space: pre-line">Starting…</span>
                <span id="progress-elapsed">0s</span>
            </div>
            <div class="run-progress-track" role="progressbar" aria-label="Analysis progress" aria-valuemin="0" aria-valuemax="100">
                <div id="progress-bar" class="run-progress-fill" style="width: 0%"></div>
            </div>
        </div>`;
    const t0 = Date.now();
    clearInterval(progressTimer);
    progressTimer = setInterval(() => {
        const el = document.getElementById("progress-elapsed");
        if (el) el.textContent = formatDuration(Date.now() - t0);
    }, 1000);
}

function setProgress(pct, label) {
    const bar = document.getElementById("progress-bar");
    if (!bar) return;
    bar.style.width = `${pct}%`;
    bar.parentElement.setAttribute("aria-valuenow", pct);
    document.getElementById("progress-label").textContent = label;
}

function stopProgress() {
    clearInterval(progressTimer);
    progressTimer = null;
}

// Sample name from a fastq file name: strip .fastq/.fq(.gz)
// --- REFERENCE: FASTA is used as is; GenBank, SnapGene .dna, FASTQ and gzipped files are converted to FASTA ---
// wrap80, parseGenbank: annotation.js (shared with the reads view)

// FASTQ (4 lines per record): keep the header and sequence, drop the qualities
function fastqToFasta(text) {
    const lines = text.split(/\r?\n/), out = [];
    for (let i = 0; i + 1 < lines.length; i += 4) {
        if (!lines[i]) { i -= 3; continue; } // skip blank lines
        if (lines[i][0] !== '@') throw new Error(`unexpected FASTQ line ${i + 1}: records must be 4 lines (@name, sequence, +, qualities)`);
        out.push(`>${lines[i].slice(1)}\n${wrap80(lines[i + 1].trim().toUpperCase())}`);
    }
    if (!out.length) throw new Error('no records found in the FASTQ file');
    return out.join('\n') + '\n';
}

async function referenceAsFasta(file) {
    const gz = /\.gz$/i.test(file.name);
    // .gz: unpack into a Blob (the browser may keep it on disk), never into one big string,
    // so large gzipped FASTA files work; only GenBank/FASTQ, which are converted, are read as text
    const src = gz ? await new Response(file.stream().pipeThrough(gunzipStream())).blob() : file;
    // detect the format from the first non-blank line
    const head = new Uint8Array(await src.slice(0, 65536).arrayBuffer());
    const first = new TextDecoder().decode(head).trimStart();
    const base = file.name.replace(/\.gz$/i, '').replace(/\.[^.]+$/, '') || 'reference';
    let fasta, annotation = null;
    if (isSnapGene(head)) { // binary SnapGene .dna: one record, named after the file
        const sg = parseSnapGene(new Uint8Array(await src.arrayBuffer()), base.replace(/\s+/g, '_'));
        if (!sg.fasta) throw new Error('no DNA sequence found in the SnapGene file');
        fasta = sg.fasta;
        annotation = sg.annotation.size ? { name: file.name, byContig: sg.annotation } : null;
    } else if (first[0] === '>') {
        if (!gz) return { file, annotation: null }; // plain FASTA: use the file directly
        return { file: new File([src], `${base}.fasta`, { type: 'text/plain' }), annotation: null };
    } else if (first.startsWith('LOCUS')) {
        const gb = parseGenbank(await src.text());
        if (!gb.fasta) throw new Error('no sequence found in the GenBank file (records need an ORIGIN section)');
        fasta = gb.fasta;
        annotation = gb.annotation.size ? { name: file.name, byContig: gb.annotation } : null;
    } else if (first[0] === '@') {
        fasta = fastqToFasta(await src.text());
    } else {
        throw new Error(`could not recognise the reference format of ${file.name} (expected FASTA, GenBank, SnapGene or FASTQ)`);
    }
    return { file: new File([fasta], `${base}.fasta`, { type: 'text/plain' }), annotation };
}

function sampleNameFromFile(fileName) {
    return fileName.replace(/\.(fastq|fq)(\.gz)?$/i, '').replace(/\.gz$/i, '')
        .replace(/(\.sorted)?\.bam$/i, '') || fileName;
}

// Paired-end files: <name>_R1 / <name>_R2 (also _1/_2, .1/.2, -1/-2, optional _001) with the
// same rest of the name. Returns [{ name, r1, r2 }] for pairs and [{ name, r1 }] for the others,
// in the order the files were given.
function pairFastqs(files) {
    const re = /^(.*?)([._-])(R?)([12])((?:_001)?\.(?:fastq|fq)(?:\.gz)?)$/i;
    const groups = new Map(), out = [];
    for (const f of files) {
        const m = f.name.match(re);
        const key = m ? `${m[1]}|${m[2]}${m[3].toUpperCase()}|${m[5].toLowerCase()}` : `single|${f.name}`;
        if (!groups.has(key)) { groups.set(key, { prefix: m?.[1], files: {}, all: [] }); out.push(groups.get(key)); }
        const g = groups.get(key);
        if (m && !g.files[m[4]]) g.files[m[4]] = f;
        g.all.push(f);
    }
    return out.flatMap(g => g.files['1'] && g.files['2'] && g.all.length === 2
        ? [{ name: g.prefix, r1: g.files['1'], r2: g.files['2'] }]
        : g.all.map(f => ({ name: sampleNameFromFile(f.name), r1: f })));
}

// --- STREAMED SAMPLING: FASTQ files are read in the browser as a stream (gzip unpacked on the
// fly), each read is kept or dropped as it passes (seeded, so runs are repeatable), and the kept
// reads are handed over in batches of about `batchBases` bases. The whole file is never in memory.
const SAMPLE_SEED = 11;
const BATCH_BASES = 100e6;         // bases per alignment batch (FASTQ + SAM of a batch: ~400 MB)
const BATCH_BASES_NO_MMI = 400e6;  // bigger batches when the reference is re-indexed for each one
const ESTIMATE_BYTES = 32 << 20;   // target depth: total bases estimated from the first 32 MB of a file

// Small seeded random number generator (mulberry32): () -> [0, 1)
function seededRandom(seed) {
    return () => {
        seed = seed + 0x6D2B79F5 | 0;
        let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
        t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
        return ((t ^ t >>> 14) >>> 0) / 4294967296;
    };
}

// A FASTQ Blob as a text stream (.gz unpacked); onBytes(n) gets the compressed bytes as they are read
function fastqTextStream(blob, name, onBytes) {
    let st = blob.stream();
    if (onBytes) st = st.pipeThrough(new TransformStream({ transform(chunk, c) { onBytes(chunk.byteLength); c.enqueue(chunk); } }));
    if (/\.gz$/i.test(name)) st = st.pipeThrough(gunzipStream()); // common.js: multi-member gzip
    return st.pipeThrough(new TextDecoderStream());
}

// FASTQ records [header, sequence, '+' line, qualities] from a text stream (4-line records)
async function* fastqRecords(text) {
    const reader = text.getReader();
    let buf = '', rec = [], lineNo = 0;
    for (;;) {
        const { value, done } = await reader.read();
        if (done) break;
        buf += value;
        let start = 0, nl;
        while ((nl = buf.indexOf('\n', start)) !== -1) {
            let line = buf.slice(start, nl);
            start = nl + 1; lineNo++;
            if (line.endsWith('\r')) line = line.slice(0, -1);
            if (!rec.length && !line) continue; // blank line between records
            rec.push(line);
            if (rec.length === 4) {
                if (rec[0][0] !== '@' || rec[2][0] !== '+') throw new Error(`not a 4-line FASTQ record near line ${lineNo - 3}`);
                yield rec;
                rec = [];
            }
        }
        buf = buf.slice(start);
    }
    if (buf.replace(/\r$/, '')) rec.push(buf.replace(/\r$/, ''));
    if (rec.length === 4) yield rec;
    else if (rec.length) throw new Error('the FASTQ file ends in the middle of a record');
}

// Total bases in a FASTQ file: exact for files up to ESTIMATE_BYTES, else estimated from the
// first ESTIMATE_BYTES (bases read there x file size / ESTIMATE_BYTES)
async function estimateFastqBases(file) {
    const exact = file.size <= ESTIMATE_BYTES;
    let bases = 0;
    try {
        for await (const r of fastqRecords(fastqTextStream(exact ? file : file.slice(0, ESTIMATE_BYTES), file.name))) bases += r[1].length;
    } catch (err) {
        if (exact) throw err; // a slice ends mid-record (and mid-gzip): use what was read
    }
    return exact ? bases : bases * file.size / ESTIMATE_BYTES;
}

// Total length of a FASTA file's sequences (bytes outside '>' lines, minus line breaks/spaces)
async function fastaTotalLength(file) {
    const reader = file.stream().getReader();
    let total = 0, inHeader = false, atLineStart = true;
    for (;;) {
        const { value: buf, done } = await reader.read();
        if (done) break;
        for (let i = 0; i < buf.length; i++) {
            const b = buf[i];
            if (atLineStart && b === 62) inHeader = true;   // '>'
            if (b === 10) { inHeader = false; atLineStart = true; continue; }
            atLineStart = false;
            if (!inHeader && b > 32) total++;
        }
    }
    return total;
}

// Stream FASTQ file(s), two = read pairs read in step, keep each read (pair) with probability p,
// and call `await onBatch(data, bases)` for every ~batchBases kept bases (data: one Uint8Array
// per file). onProgress(fractionOfFilesRead, stats). Returns { reads, kept, keptBases }.
async function sampleFastqBatches(files, p, batchBases, onBatch, onProgress = () => {}) {
    const totalBytes = files.reduce((t, f) => t + f.size, 0) || 1;
    let readBytes = 0;
    const iters = files.map(f => fastqRecords(fastqTextStream(f, f.name, n => { readBytes += n; })));
    const rand = seededRandom(SAMPLE_SEED), enc = new TextEncoder();
    let parts = files.map(() => []), pieces = files.map(() => []), partChars = 0, inBatch = 0;
    const stats = { reads: 0, kept: 0, keptBases: 0 };
    const flushParts = () => { // strings -> bytes every ~16 MB, so no string gets large
        parts.forEach((pa, k) => { pieces[k].push(enc.encode(pa.join(''))); parts[k] = []; });
        partChars = 0;
    };
    const flushBatch = async () => {
        flushParts();
        const data = pieces.map(ps => {
            const out = new Uint8Array(ps.reduce((t, x) => t + x.length, 0));
            let o = 0;
            for (const x of ps) { out.set(x, o); o += x.length; }
            return out;
        });
        pieces = files.map(() => []);
        const b = inBatch;
        inBatch = 0;
        await onBatch(data, b);
    };
    for (;;) {
        const recs = await Promise.all(iters.map(it => it.next()));
        if (recs.every(r => r.done)) break;
        if (recs.some(r => r.done)) throw new Error('the R1 and R2 files have different numbers of reads');
        stats.reads++;
        if (p >= 1 || rand() < p) {
            recs.forEach((r, k) => {
                const [h, seq, , qual] = r.value;
                parts[k].push(h, '\n', seq, '\n+\n', qual, '\n');
                partChars += h.length + 2 * seq.length + 5;
                inBatch += seq.length; stats.keptBases += seq.length;
            });
            stats.kept++;
            if (partChars > 16e6) flushParts();
            if (inBatch >= batchBases) await flushBatch();
        }
        if (stats.reads % 5000 === 0) onProgress(readBytes / totalBytes, stats);
    }
    if (inBatch) await flushBatch();
    onProgress(1, stats);
    return stats;
}

// Sample, align and analyse one FASTQ sample (`files`: one File, or R1 + R2 of read pairs).
// sampling: { mode: 'depth', value: target depth } or { mode: 'pct', value: % of reads };
// refBp: total reference length (for the target depth). `progress(frac, label)`: 0..1.
// The reads stream through sampleFastqBatches; each batch is aligned (minimap2) and sorted into
// its own small BAM, and the batch BAMs are merged at the end, so no full-size SAM ever exists.
async function runSample({ name, files, target, preset, sampling, refBp, batchBases, index, aligners = [CLI] }, progress) {
    const tag = `s${index}`;
    let bam = `${tag}.sorted.bam`;
    const paired = files.length === 2;

    // share of reads to keep
    let p = 1, estBases = null;
    if (sampling.mode === 'pct') {
        p = Math.min(1, sampling.value / 100);
    } else {
        progress(0, 'Estimating how much sequence the reads contain…');
        estBases = 0;
        for (const f of files) estBases += await estimateFastqBases(f);
        p = estBases > 0 ? Math.min(1, sampling.value * refBp / estBases) : 1;
    }
    const keepNote = p >= 1 ? 'all reads' : `${(p * 100).toFixed(p < 0.1 ? 2 : 1)}% of reads`;

    // Several workers: batches small enough that each worker gets about two of them
    if (aligners.length > 1) {
        if (estBases === null) {
            progress(0, 'Estimating how much sequence the reads contain…');
            estBases = 0;
            for (const f of files) estBases += await estimateFastqBases(f);
        }
        const kept = p * estBases;
        if (kept > 0) batchBases = Math.max(MIN_BATCH_BASES, Math.min(batchBases, Math.ceil(kept / (2 * aligners.length))));
    }

    // Align + sort one batch on worker `a`; the batch BAM ends up in CLI's file system.
    const alignBatch = async (a, k, data) => {
        const fq = data.map((_, j) => `${tag}.b${k}_${j + 1}.fastq`);
        for (const [j, d] of data.entries()) await a.fs.writeFile(fq[j], d);
        const sam = `${tag}.b${k}.sam`, out = `${tag}.b${k}.bam`;
        await a.exec("minimap2", [
            "-a",
            "-K", "50M", // smaller query batches -> lower peak memory (default 500M)
            "-x", preset,
            ...(paired ? ["--frag=yes"] : []), // two files: mates of read pairs, aligned as fragments
            "-o", sam, target, ...fq
        ]);
        for (const f of fq) await a.fs.unlink(f);
        // -m caps the sort buffer so samtools spills to temp files; -l 1: fast compression (the
        // batch BAMs are only read back by the merge)
        await a.exec("samtools", ["sort", "-m", "200M", "-l", "1", "-o", out, sam]);
        await a.fs.unlink(sam);
        if (a !== CLI) { // move the BAM over to CLI, where it is merged and analysed
            await CLI.fs.writeFile(out, await a.fs.readFile(out));
            await a.fs.unlink(out);
        }
        return out;
    };

    // Batches go to whichever worker is free; reading the FASTQ waits only when all are busy.
    const free = [...aligners], waiting = [], running = new Set();
    const batchBams = [];
    let failure = null, nBatches = 0, nDone = 0;
    const takeWorker = () => free.length ? Promise.resolve(free.pop()) : new Promise(res => waiting.push(res));
    const giveBack = a => { const w = waiting.shift(); if (w) w(a); else free.push(a); };
    const busyNote = () => aligners.length > 1 ? ` · aligning on ${running.size}/${aligners.length} workers` : '';
    const batchNote = () => `\n${nDone}/${nBatches} batches aligned${busyNote()}`; // second line of the label

    let readPct = 0, stats;
    try {
        stats = await sampleFastqBatches(files, p, batchBases, async (data, bases) => {
            if (failure) throw failure;
            const a = await takeWorker();
            if (failure) throw failure;
            const k = nBatches++;
            progress(0.02 + 0.73 * readPct, `Aligning batch ${k + 1} (${formatBp(bases)}) with minimap2…${batchNote()}`);
            const job = alignBatch(a, k, data)
                .then(out => { batchBams[k] = out; nDone++; })
                .catch(err => { failure ||= err; })
                .finally(() => { running.delete(job); giveBack(a); });
            running.add(job);
        }, (frac, st) => {
            readPct = frac;
            progress(0.02 + 0.73 * frac, `Reading and sampling ${paired ? 'read pairs' : 'reads'} (${keepNote}): ` +
                `${Math.round(frac * 100)}% of the file · ${st.kept.toLocaleString()} kept${batchNote()}`);
        });
    } catch (err) {
        await Promise.allSettled(running); // let the workers finish before the error is reported
        throw err;
    }
    if (running.size) progress(0.75, `Finishing alignment of the last ${running.size} batch${running.size > 1 ? 'es' : ''}…`);
    await Promise.all(running);
    if (failure) throw failure;
    if (!batchBams.length) throw new Error(`${name}: no reads ${stats.reads ? 'were sampled (try a higher target)' : 'found in the FASTQ file'}`);

    if (batchBams.length === 1) {
        bam = batchBams[0]; // one batch: it is the BAM
    } else {
        progress(0.76, `Merging ${batchBams.length} batches…`);
        // -l 1: fast compression. The merge is single-threaded and mostly compression; the BAM
        // comes out maybe 20-30% larger than at the default level.
        await CLI.exec("samtools", ["merge", "-f", "-l", "1", "-o", bam, ...batchBams]);
        for (const f of batchBams) await CLI.fs.unlink(f);
    }
    progress(0.78, "Indexing BAM…");
    await CLI.exec("samtools", ["index", bam]);
    const sample = await analyseBam({ name, bam }, (f, label) => progress(0.8 + 0.2 * f, label));
    sample.sampling = {
        mode: sampling.mode, value: sampling.value, p, reads: stats.reads, kept: stats.kept,
        keptBases: stats.keptBases, depth: refBp ? stats.keptBases / refBp : null, estBases
    };
    return sample;
}

// An existing BAM as a sample: index it if no index came with it, then analyse it.
// `bam` is the mounted BAM; `hasIndex` whether its .bai/.csi was mounted next to it.
async function bamSample({ name, bam, file, hasIndex, index }, progress) {
    let path = bam;
    if (!hasIndex) {
        // the mounted files are read-only, so the BAM is copied into memory to be indexed next to it
        progress(0, `Indexing ${file.name} (no .bai given; copying it into memory)…`);
        path = `b${index}.bam`;
        await CLI.fs.writeFile(path, new Uint8Array(await file.arrayBuffer()));
        await CLI.exec("samtools", ["index", path]);
    }
    // The viewer queries regions, which needs a working index: try one on the first contig.
    // A BAM that isn't coordinate-sorted (so can't be indexed), or isn't a BAM, fails here.
    const header = await CLI.exec("samtools", ["view", "-H", path]);
    const first = header.match(/^@SQ\t(?:.*\t)?SN:(\S+)/m)?.[1];
    if (!first) throw new Error(`${file.name}: no reference sequences (@SQ lines) in the BAM header. Is it a BAM file?`);
    const probe = await CLI.exec("samtools", ["view", "-c", path, regionOf(first, 0, 1)]);
    if (!/^\d+\s*$/.test(probe.trim().split('\n').pop() || '')) {
        throw new Error(`${file.name}: could not use its index. Is the BAM coordinate-sorted (samtools sort) and the .bai/.csi made from this file? ${probe.trim().split('\n')[0] || ''}`);
    }
    const sample = await analyseBam({ name, bam: path }, (f, label) => progress(0.1 + 0.9 * f, label));
    sample.bamFile = file; // downloaded as is
    return sample;
}

// Stats, coverage and depth profiles for a sorted, indexed BAM. `progress(frac, label)` 0..1.
async function analyseBam({ name, bam }, progress) {
    // --- read counts (flagstat) ---
    progress(0, "Counting mapped reads…");
    const flagstat = await CLI.exec("samtools", ["flagstat", bam]);
    const grab = (re) => { const m = flagstat.match(re); return m ? parseInt(m[1], 10) : 0; };
    const total = grab(/^(\d+)\s*\+\s*\d+\s*in total/m);
    const secondary = grab(/^(\d+)\s*\+\s*\d+\s*secondary/m);
    const supplementary = grab(/^(\d+)\s*\+\s*\d+\s*supplementary/m);
    const mappedAll = grab(/^(\d+)\s*\+\s*\d+\s*mapped\s+\(/m);
    // flagstat counts secondary/supplementary records as mapped; unmapped records are always primary
    const reads_ = Math.max(0, total - secondary - supplementary);
    const mapped = Math.max(0, mappedAll - secondary - supplementary);
    const unmapped = Math.max(0, total - mappedAll);

    // --- per-contig stats (samtools coverage) ---
    // Columns: 0:#rname 1:startpos 2:endpos 3:numreads 4:covbases 5:coverage 6:meandepth 7:meanbaseq 8:meanmapq
    progress(0.25, "Computing per-contig coverage…");
    const covOutput = await CLI.exec("samtools", ["coverage", "-H", bam]);
    const contigs = [];
    const num = (x) => { const v = parseFloat(x); return isNaN(v) ? 0 : v; };
    for (const line of covOutput.trim().split('\n')) {
        const f = line.split('\t');
        if (f.length < 9) continue;
        contigs.push({
            name: cleanContigName(f[0]),
            length: num(f[2]),
            numreads: num(f[3]),
            coverage: num(f[5]).toFixed(0),
            meandepth: num(f[6]).toFixed(0),
            meanmapq: num(f[8]).toFixed(0),
        });
    }

    // --- binned depth profiles: mean + lowest depth per bin, so narrow dips still show ---
    progress(0.6, "Building depth profiles…");
    for (const c of contigs) {
        c.binSize = Math.max(1, Math.min(c.length, Math.max(MIN_BIN_BP, Math.ceil(c.length / NBINS))));
        maxContigLength = Math.max(maxContigLength, c.length);
    }
    const { profiles, minProfiles } = await depthProfiles(contigs, bam,
        f => progress(0.6 + 0.4 * f, `Building depth profiles… ${Math.round(f * 100)}%`));
    for (const c of contigs) {
        c.profile = profiles[c.name] || [];
        c.minProfile = minProfiles[c.name] || [];
        c.warnings = contigWarnings(c, minProfiles[c.name]);
    }

    return {
        name, bam, contigs,
        summary: {
            reads: reads_, mapped, unmapped,
            mappedPct: reads_ > 0 ? (mapped / reads_ * 100) : 0,
            nWarn: contigs.filter(c => c.warnings.length).length
        }
    };
}

async function run() {
    const outputDiv = document.getElementById("output");
    const resultsDiv = document.getElementById("results");

    showProgress(outputDiv);
    document.getElementById("btn").disabled = true;
    resultsDiv.innerHTML = '';
    samples = [];
    maxContigLength = 1;

    try {
        const q = [...document.getElementById("query").files];
        const r = document.getElementById("ref").files;
        // reads to align (FASTQ), and existing alignments (BAM, with an optional .bai/.csi each)
        const isIndex = f => /\.(bai|csi)$/i.test(f.name);
        const isBam = f => /\.bam$/i.test(f.name);
        const fastqs = q.filter(f => !isIndex(f) && !isBam(f));
        const bams = q.filter(isBam);
        const indexes = q.filter(isIndex);
        const stem = n => n.replace(/\.bam$/i, '');
        const indexFor = b => indexes.find(i => [`${b.name}.bai`, `${b.name}.csi`, `${stem(b.name)}.bai`, `${stem(b.name)}.csi`].includes(i.name));
        const orphanIdx = indexes.filter(i => !bams.some(b => indexFor(b) === i));

        if (fastqs.length === 0 && bams.length === 0) {
            outputDiv.innerHTML = '<div class="alert alert-danger" role="alert">Please select <b>FASTQ</b> read files to align, or <b>BAM</b> files (with their .bai index) to view.</div>';
            return;
        }
        if (fastqs.length && r.length === 0) {
            outputDiv.innerHTML = '<div class="alert alert-danger" role="alert">Aligning FASTQ reads needs a <b>Reference</b>. (For BAM files it is optional: it adds mismatches and variants.)</div>';
            return;
        }
        if (orphanIdx.length) {
            outputDiv.innerHTML = `<div class="alert alert-danger" role="alert">No BAM selected for index file(s) <b>${orphanIdx.map(i => escapeHtml(i.name)).join(', ')}</b>. Name them like the BAM (<i>x.bam.bai</i> or <i>x.bai</i>).</div>`;
            return;
        }

        const sampling = { mode: document.getElementById("sample-mode").value, value: parseFloat(document.getElementById("sample-value").value) };
        if (sampling.mode === 'pct' ? !(sampling.value > 0 && sampling.value <= 100) : !(sampling.value > 0)) {
            outputDiv.innerHTML = `<div class="alert alert-danger" role="alert">${sampling.mode === 'pct'
                ? 'The share of reads must be between 0 and 100%.' : 'The target depth must be above 0.'}</div>`;
            return;
        }
        const preset = document.getElementById("preset").value;

        setProgress(2, "Loading files…");
        // BAMs are mounted (FASTQs are streamed from the File objects instead);
        // each index is mounted under the name samtools looks for: <bam name>.bai / .csi
        const mountList = [];
        for (const b of bams) {
            mountList.push(b);
            const idx = indexFor(b);
            if (idx) mountList.push(new File([idx], `${b.name}${idx.name.toLowerCase().endsWith('.csi') ? '.csi' : '.bai'}`));
        }
        const mounted = mountList.length ? await CLI.mount(mountList) : [];
        const pathOf = f => mounted[mountList.indexOf(f)];

        referenceFile = null; referenceAnnotation = null;
        contigSeqCache.clear();
        let refPath = null, refFile = null;
        if (r.length) {
            setProgress(3, "Reading reference…");
            const conv = await referenceAsFasta(r[0]);
            refFile = conv.file;
            referenceAnnotation = conv.annotation;
            refPath = (await CLI.mount([refFile]))[0];
            referenceFile = refFile;
        }

        // One sample per FASTQ or BAM file; make names unique
        const seen = new Map();
        const uniqueName = (file) => {
            let name = sampleNameFromFile(file.name);
            const k = (seen.get(name) || 0) + 1;
            seen.set(name, k);
            return k > 1 ? `${name}_${k}` : name;
        };
        const jobs = [
            ...pairFastqs(fastqs).map(p => ({
                kind: 'fastq', name: uniqueName({ name: p.name }), files: p.r2 ? [p.r1, p.r2] : [p.r1], paired: !!p.r2
            })),
            ...bams.map(file => ({ kind: 'bam', name: uniqueName(file), file, bam: pathOf(file), hasIndex: !!indexFor(file) })),
        ].map((j, i) => ({ ...j, index: i }));

        // Reads are aligned in batches: build the minimap2 index once and reuse it for all of them.
        // Big references skip the .mmi copy (memory) and are re-indexed per batch, so their batches are bigger.
        // With an index, batches are aligned on several workers at once (see getAligners).
        let target = refPath, batchBases = BATCH_BASES, refBp = 0, aligners = [CLI];
        if (fastqs.length) {
            if (refFile.size <= MMI_MAX_REF_BYTES) {
                setProgress(4, "Indexing reference…");
                await CLI.exec("minimap2", ["-x", preset, "-d", "ref.mmi", refPath]);
                target = "ref.mmi";
                const n = alignerCount();
                if (n > 1) {
                    setProgress(4, `Starting ${n} alignment workers…`);
                    aligners = await getAligners(n, target);
                    batchBases = PARALLEL_BATCH_BASES;
                }
            } else {
                batchBases = BATCH_BASES_NO_MMI;
            }
            if (sampling.mode === 'depth') refBp = await fastaTotalLength(refFile);
        }

        const n = jobs.length;
        for (const job of jobs) {
            const report = (frac, label) => {
                const pct = 5 + 95 * (job.index + frac) / n;
                setProgress(pct, n > 1 ? `Sample ${job.index + 1}/${n} (${job.name}): ${label}` : label);
            };
            const sample = job.kind === 'bam'
                ? await bamSample(job, report)
                : await runSample({ ...job, target, preset, sampling, refBp, batchBases, aligners }, report);
            sample.fromBam = job.kind === 'bam';
            sample.paired = !!job.paired;
            samples.push(sample);
        }
        if (target === "ref.mmi") {
            await releaseAligners(aligners, target);
            await CLI.fs.unlink("ref.mmi");
        }

        if (samples.every(sm => sm.contigs.length === 0)) {
            outputDiv.innerHTML = '<div class="alert alert-warning" role="alert">Analysis completed, but no coverage data was generated. Check your input files.</div>';
            return;
        }

        renderResults();
        outputDiv.innerHTML = ''; // done: remove progress bar
        // Results are shown: hide the inputs and Run button (RESET starts over)
        document.getElementById("inputs").style.display = 'none';
        document.getElementById("btn").style.display = 'none';

    } catch (error) {
        console.error(error);
        outputDiv.innerHTML = `<div class="alert alert-danger" role="alert">An error occurred during execution: <b>${escapeHtml(error.message || String(error) || 'Unknown Aioli/CLI Error')}</b></div>`;
    } finally {
        stopProgress();
        document.getElementById("btn").disabled = false;
    }
}

// --- RESULTS: combined summary + one tab per sample ---
// How a sample's reads were sampled, e.g. "22.4% of 332,046 reads, ~100x (target 100x)"
function samplingNote(sm) {
    const sp = sm.sampling;
    if (!sp) return '';
    const depth = sp.depth ? `~${sp.depth >= 10 ? Math.round(sp.depth) : sp.depth.toFixed(1)}x` : '';
    const target = sp.mode === 'depth' ? ` (target ${sp.value}x${sp.p >= 1 ? ': all reads used, the file has less' : ''})` : '';
    const share = sp.p >= 1 ? `all ${sp.reads.toLocaleString()} reads` : `${(sp.kept / sp.reads * 100).toFixed(1)}% of ${sp.reads.toLocaleString()} reads`;
    return `${share}${depth ? `, ${depth}` : ''}${target}`;
}

function renderResults() {
    const resultsDiv = document.getElementById("results");
    const multi = samples.length > 1;
    const readsLabel = samples.some(sm => sm.sampling?.p < 1) ? 'Reads used' : 'Reads';

    let html = '';
    if (multi) {
        const rows = samples.map((sm, si) => {
            const s = sm.summary;
            return `<tr data-sample="${si}">
                <td class="text-start">${escapeHtml(sm.name)}</td>
                <td>${s.reads.toLocaleString()}</td>
                <td>${s.mapped.toLocaleString()}</td>
                <td>${s.mappedPct.toFixed(1)}%</td>
                <td>${s.unmapped.toLocaleString()}</td>
                <td>${s.nWarn ? `<span class="warn-badge">!</span>${s.nWarn}` : '–'}</td>
            </tr>`;
        }).join('');
        html += `
            <div class="combined-summary">
                <div class="combined-head">
                    <span class="fw-bold">${samples.length} samples</span>
                    <button class="btn btn-outline-secondary btn-sm" id="download-all-csv">All stats CSV</button>
                </div>
                <table class="combined-table">
                    <thead><tr><th class="text-start">Sample</th><th>${readsLabel}</th><th>Mapped</th><th>Mapped %</th><th>Unmapped</th><th>Contigs with warnings</th></tr></thead>
                    <tbody>${rows}</tbody>
                </table>
            </div>
            <ul class="sample-tabs" role="tablist">
                ${samples.map((sm, si) => `<li><button role="tab" data-sample="${si}" class="${si === 0 ? 'active' : ''}">${sm.summary.nWarn ? '<span class="warn-dot"></span>' : ''}${escapeHtml(sm.name)}</button></li>`).join('')}
            </ul>`;
    }
    html += samples.map((sm, si) => samplePanelHtml(sm, si, readsLabel, si === 0)).join('');
    resultsDiv.innerHTML = html;

    showSample(0);
    if (multi) document.getElementById("download-all-csv").addEventListener("click", () => downloadStatsCsv(null));
}

// The contig list draws at most MAX_LIST_CONTIGS rows (each has a sparkline); with more
// contigs, a name filter finds the others. Stats CSV and consensus always cover all contigs.
const MAX_LIST_CONTIGS = 100;

// List rows for the contig indexes `idx` (indexes into sm.contigs)
function contigItemsHtml(sm, si, idx) {
    return idx.map(i => {
        const c = sm.contigs[i];
        const sparklineHtml = `<div class="profile-track" data-length='${c.length}' data-bin='${c.binSize}'><span class="sparkline-container" data-depths='${JSON.stringify(c.profile)}' data-mins='${JSON.stringify(c.minProfile)}'></span></div>`;
        return `
            <li class="depth-item clickable" data-sample="${si}" data-index="${i}">
                <div class="contig-name" data-tip="${escapeHtml(`${c.name}\nclick to view reads`)}">${c.warnings.length ? `<span class="warn-badge" data-tip="${escapeHtml(c.warnings.join('\n'))}">!</span>` : ''}${escapeHtml(c.name)}</div>
                <div class="mean-depth">${c.meandepth}</div>
                <div class="coverage">${c.coverage}%</div>
                <div class="mean-mapq">${c.meanmapq}</div>
                <div class="depth-profile">${sparklineHtml}</div>
            </li>`;
    }).join('');
}

// Indexes of the contigs to list: the first MAX_LIST_CONTIGS, or those matching `query`
function listedContigs(sm, query = '') {
    const q = query.trim().toLowerCase();
    const idx = [];
    for (let i = 0; i < sm.contigs.length && idx.length < MAX_LIST_CONTIGS; i++) {
        if (!q || sm.contigs[i].name.toLowerCase().includes(q)) idx.push(i);
    }
    return idx;
}

function listNote(sm, shown, query) {
    const n = sm.contigs.length;
    if (n <= MAX_LIST_CONTIGS && !query.trim()) return '';
    if (query.trim()) {
        const all = sm.contigs.filter(c => c.name.toLowerCase().includes(query.trim().toLowerCase())).length;
        return all ? `${all.toLocaleString()} contig${all > 1 ? 's' : ''} match${all > MAX_LIST_CONTIGS ? `, showing the first ${MAX_LIST_CONTIGS}` : ''}` : 'no contig matches';
    }
    return `Showing the first ${shown} of ${n.toLocaleString()} contigs. Type a name above to find the others (Stats CSV includes all).`;
}

function samplePanelHtml(sm, si, readsLabel, visible) {
    const s = sm.summary;
    const shownIdx = listedContigs(sm);
    const items = contigItemsHtml(sm, si, shownIdx);
    const many = sm.contigs.length > MAX_LIST_CONTIGS;
    const safe = sm.name.replace(/[^\w.-]+/g, '_');
    return `
        <div class="sample-panel" data-sample="${si}" ${visible ? '' : 'hidden'}>
            <div class="alert alert-info summary-box" role="alert">
                <span class="fw-bold">Sample:</span> ${escapeHtml(sm.name)} &nbsp;|&nbsp;
                <span class="fw-bold">${readsLabel}:</span> ${s.reads.toLocaleString()}${sm.sampling ? ` <span class="text-muted">(${escapeHtml(samplingNote(sm))})</span>` : ''} &nbsp;|&nbsp;
                <span class="fw-bold">Mapped reads:</span> ${s.mapped.toLocaleString()} (${s.mappedPct.toFixed(2)}%) &nbsp;|&nbsp;
                <span class="fw-bold">Unmapped:</span> ${s.unmapped.toLocaleString()}
                ${s.nWarn ? `&nbsp;|&nbsp; <span class="warn-badge">!</span>${s.nWarn} of ${sm.contigs.length} contigs with warnings (hover the ! for details)` : ''}
            </div>
            <div class="sample-downloads">
                <button class="btn btn-success me-2" data-action="bam" data-sample="${si}">Download BAM</button>
                <button class="btn btn-outline-secondary me-2" data-action="csv" data-sample="${si}">Stats CSV</button>
                <button class="btn btn-outline-secondary me-2" data-action="consensus" data-sample="${si}" data-tip="${escapeHtml(`Consensus of each contig from all mapped reads (samtools consensus, simple mode; any MAPQ, secondary/QC-fail/duplicate reads excluded).\nA base is called when ≥ ${CONS_CALL_FRACT * 100}% of the reads carry it, else N; depth < ${MIN_CONS_DEPTH} is N. Deletions are removed, insertions are included.`)}">Consensus FASTA</button>
                ${s.unmapped > 0 ? `<button class="btn btn-outline-secondary" data-action="unmapped" data-sample="${si}">Unmapped FASTQ</button>` : ''}
            </div>
            ${many ? `<input type="search" class="form-control form-control-sm contig-filter" data-sample="${si}"
                placeholder="Find a contig by name (${sm.contigs.length.toLocaleString()} contigs)" aria-label="Find a contig by name">` : ''}
            <div class="depth-list-headers depth-item">
                <div class="contig-name">Contig</div>
                <div class="mean-depth">Depth</div>
                <div class="coverage">Cov (%)</div>
                <div class="mean-mapq">MAPQ</div>
                <div class="depth-profile">Depth Profile</div>
            </div>
            <ul class="depth-list" data-sample="${si}">${items}</ul>
            <p class="list-note text-muted" data-sample="${si}">${listNote(sm, shownIdx.length, '')}</p>
        </div>`;
}

// Contig name filter: redraw the list (at most MAX_LIST_CONTIGS rows) for the matching contigs
let filterTimer = null;
document.getElementById("results").addEventListener("input", (e) => {
    const input = e.target.closest('.contig-filter');
    if (!input) return;
    clearTimeout(filterTimer);
    filterTimer = setTimeout(() => {
        const si = parseInt(input.dataset.sample, 10), sm = samples[si];
        const idx = listedContigs(sm, input.value);
        const panel = input.closest('.sample-panel');
        const ul = panel.querySelector('.depth-list');
        ul.innerHTML = contigItemsHtml(sm, si, idx);
        renderSparklines(ul);
        panel.querySelector('.list-note').textContent = listNote(sm, idx.length, input.value);
    }, 200);
});

// Show one sample's panel; sparklines are drawn the first time a panel is visible
// (jquery.sparkline can't measure widths inside hidden elements)
function showSample(si) {
    document.querySelectorAll('#results .sample-panel').forEach(p => {
        const on = parseInt(p.dataset.sample, 10) === si;
        p.hidden = !on;
        if (on && !p.dataset.rendered) { renderSparklines(p); p.dataset.rendered = '1'; }
    });
    document.querySelectorAll('#results .sample-tabs button').forEach(b => {
        b.classList.toggle('active', parseInt(b.dataset.sample, 10) === si);
    });
    document.querySelectorAll('#results .combined-table tbody tr').forEach(tr => {
        tr.classList.toggle('active', parseInt(tr.dataset.sample, 10) === si);
    });
}

// --- INSTANT TOOLTIPS for the results table (native title tooltips appear only after ~1 s) ---
const tipEl = Object.assign(document.createElement('div'), { className: 'quick-tip' });
document.body.appendChild(tipEl);
attachQuickTips([document.getElementById("results"), document.getElementById("inputs")], tipEl);

// --- CONTIG WARNINGS ---
const WARN_MIN_COVERAGE = 90; // % of bases covered
const WARN_MIN_DEPTH = 5;     // mean depth

// Returns a list of human-readable warnings for one contig
function contigWarnings(c, mins) {
    const w = [];
    const cov = parseFloat(c.coverage), depth = parseFloat(c.meandepth);
    if (c.numreads === 0) return ['no reads mapped'];
    if (cov < WARN_MIN_COVERAGE) w.push(`only ${cov}% of bases covered`);
    if (depth < WARN_MIN_DEPTH) w.push(`low mean depth (${depth}x)`);
    // zero-depth stretches away from the contig ends (first/last bins naturally thin out)
    if (mins && mins.length > 2) {
        let stretches = 0;
        for (let i = 1; i < mins.length - 1; i++) {
            if (mins[i] === 0 && mins[i - 1] !== 0) stretches++;
        }
        if (stretches) w.push(`${stretches} zero-depth stretch${stretches > 1 ? 'es' : ''} inside the contig`);
    }
    return w;
}

// --- DOWNLOADS --- (downloadBlob: common.js)

const csvCell = (v) => /[",\n]/.test(String(v)) ? `"${String(v).replace(/"/g, '""')}"` : String(v);

// Stats CSV for one sample (index) or all samples (null), with a sample column
function downloadStatsCsv(si) {
    const list = si === null ? samples : [samples[si]];
    const header = ['sample', 'contig', 'length', 'reads', 'covered_pct', 'mean_depth', 'mean_mapq', 'warnings'];
    const lines = [header.join(',')];
    for (const sm of list) {
        for (const c of sm.contigs) {
            lines.push([sm.name, c.name, c.length, c.numreads, c.coverage, c.meandepth, c.meanmapq, (c.warnings || []).join('; ')]
                .map(csvCell).join(','));
        }
    }
    const fname = si === null ? 'contig_stats_all_samples.csv' : `${list[0].name.replace(/[^\w.-]+/g, '_')}_contig_stats.csv`;
    downloadBlob(new Blob([lines.join('\n') + '\n'], { type: 'text/csv' }), fname);
}

// The sample's BAM: a BAM given as input is saved as is; an aligned one is read out of the
// in-browser file system only now (keeping a second copy ready would double its memory)
async function downloadBam(si, btn) {
    const sm = samples[si];
    const name = `${sm.name.replace(/[^\w.-]+/g, '_')}.sorted.bam`;
    if (sm.bamFile) return downloadBlob(sm.bamFile, sm.bamFile.name);
    const label = btn.textContent;
    btn.disabled = true;
    btn.textContent = 'Preparing…';
    try {
        downloadBlob(new Blob([await CLI.fs.readFile(sm.bam)], { type: 'application/octet-stream' }), name);
    } catch (err) {
        console.error(err);
        alert(`Could not read the BAM: ${err.message || err}`);
    } finally {
        btn.disabled = false;
        btn.textContent = label;
    }
}

// Unmapped reads as FASTQ. samtools fastq writes reads to a file (-0) so its
// stderr messages can't end up mixed into the FASTQ.
async function downloadUnmapped(si, btn) {
    const sm = samples[si];
    const label = btn.textContent;
    btn.disabled = true;
    btn.textContent = 'Extracting…';
    try {
        // Paired reads (READ1/READ2 flags) go to -1/-2, a mate whose partner is mapped to -s,
        // reads without pair flags to -0. Each non-empty file is downloaded.
        const outs = [["u0.fastq", ""], ["u1.fastq", "_R1"], ["u2.fastq", "_R2"], ["us.fastq", "_singletons"]];
        await CLI.exec("samtools", ["fastq", "-f", "4", "-0", "u0.fastq", "-1", "u1.fastq", "-2", "u2.fastq", "-s", "us.fastq", sm.bam]);
        const base = sm.name.replace(/[^\w.-]+/g, '_');
        let n = 0;
        for (const [f, suffix] of outs) {
            let data = null;
            try { data = await CLI.fs.readFile(f); await CLI.fs.unlink(f); } catch { continue; }
            if (!data?.length) continue;
            downloadBlob(new Blob([data], { type: 'text/plain' }), `${base}_unmapped${suffix}.fastq`);
            n++;
        }
        if (!n) alert('No unmapped reads found.');
    } catch (err) {
        console.error(err);
        alert(`Could not extract unmapped reads: ${err.message || err}`);
    } finally {
        btn.disabled = false;
        btn.textContent = label;
    }
}

// --- CONSENSUS FASTA (samtools consensus, one contig at a time for progress) ---
// Simple mode: per position the base carried by >= CONS_CALL_FRACT of the reads (else N),
// N below MIN_CONS_DEPTH reads; deleted bases are dropped and insertions added.
// Any MAPQ; unmapped, secondary, QC-fail and duplicate reads are excluded (samtools default).
// Works from the BAM alone, so no reference is needed.
const MIN_CONS_DEPTH = 3;
const CONS_CALL_FRACT = 0.5;
const CONS_ARGS = ["consensus", "-a", "-m", "simple", "-d", String(MIN_CONS_DEPTH), "-c", String(CONS_CALL_FRACT),
    "--show-del", "no", "--show-ins", "yes", "-f", "fasta"];

// samtools region for a contig, or for [start, end) of it (0-based). Names containing ':'
// are wrapped in {} as htslib allows, so they can't be misread as name:start-end.
function regionOf(name, start, end) {
    const n = name.includes(':') ? `{${name}}` : name;
    return start === undefined ? n : `${n}:${start + 1}-${end}`;
}

// A progress bar shown under a sample's download buttons while a long task runs:
// returns { set(frac, label), cancelled(), close() }
function taskProgress(btn, title) {
    const box = document.createElement('div');
    box.className = 'run-progress task-progress';
    box.innerHTML = `
        <div class="run-progress-text">
            <span><b>${escapeHtml(title)}</b> <span class="tp-label">starting…</span></span>
            <span><span class="tp-pct">0%</span> · <span class="tp-elapsed">0s</span>
                <button type="button" class="btn btn-link btn-sm p-0 ms-2 tp-cancel">Cancel</button></span>
        </div>
        <div class="run-progress-track" role="progressbar" aria-label="${escapeHtml(title)}" aria-valuemin="0" aria-valuemax="100">
            <div class="run-progress-fill" style="width: 0%"></div>
        </div>`;
    btn.closest('.sample-downloads').after(box);
    let stop = false;
    box.querySelector('.tp-cancel').addEventListener('click', e => {
        stop = true;
        e.target.disabled = true;
        box.querySelector('.tp-label').textContent = 'cancelling…';
    });
    const t0 = Date.now();
    const timer = setInterval(() => {
        box.querySelector('.tp-elapsed').textContent = formatDuration(Date.now() - t0);
    }, 1000);
    return {
        set(frac, label) {
            if (stop) return;
            const pct = Math.max(0, Math.min(100, Math.round(frac * 100)));
            box.querySelector('.run-progress-fill').style.width = `${pct}%`;
            box.querySelector('.run-progress-track').setAttribute('aria-valuenow', pct);
            box.querySelector('.tp-pct').textContent = `${pct}%`;
            if (label) box.querySelector('.tp-label').textContent = label;
        },
        cancelled: () => stop,
        close() { clearInterval(timer); box.remove(); }
    };
}

async function downloadConsensus(si, btn) {
    const sm = samples[si];
    btn.disabled = true;
    const prog = taskProgress(btn, `Consensus FASTA · ${sm.name}`);
    try {
        const fasta = [], skipped = [];
        // progress weighted by contig length
        const total = sm.contigs.reduce((t, c) => t + c.length, 0) || 1;
        const n = sm.contigs.length;
        let done = 0;
        for (const [k, c] of sm.contigs.entries()) {
            if (prog.cancelled()) throw new Error('cancelled');
            prog.set(done / total, `${n > 1 ? `contig ${k + 1}/${n}: ` : ''}${c.name}`);
            const header = `>${c.name}_consensus sample=${sm.name} samtools_consensus=simple min_depth=${MIN_CONS_DEPTH} call_fract=${CONS_CALL_FRACT}`;
            if (!c.numreads) { // no reads: all N
                fasta.push(header, wrap80('N'.repeat(c.length)));
            } else {
                const log = await CLI.exec("samtools", [...CONS_ARGS, "-r", regionOf(c.name), "-o", "consensus.fa", sm.bam]);
                let text;
                try { text = new TextDecoder().decode(await CLI.fs.readFile("consensus.fa")); }
                catch { throw new Error(`samtools consensus failed for ${c.name}: ${String(log).trim().split('\n')[0]}`); }
                await CLI.fs.unlink("consensus.fa");
                const seq = text.split('\n').filter(l => l && l[0] !== '>').join('');
                if (seq) fasta.push(header, wrap80(seq)); else skipped.push(`${c.name}: no consensus produced`);
            }
            done += c.length;
        }
        prog.set(1, 'writing the file…');
        if (fasta.length) downloadBlob(new Blob([fasta.join('\n') + '\n'], { type: 'text/plain' }), `${sm.name.replace(/[^\w.-]+/g, '_')}_consensus.fasta`);
        if (skipped.length) alert(`No consensus for ${skipped.length} contig(s):\n${skipped.slice(0, 5).join('\n')}${skipped.length > 5 ? '\n…' : ''}`);
    } catch (err) {
        if (err?.message !== 'cancelled') {
            console.error(err);
            alert(`Could not build the consensus: ${err.message || err}`);
        }
    } finally {
        prog.close();
        btn.disabled = false;
    }
}

// Clicks in the results: tabs, combined-summary rows, downloads, contig rows (reads view)
document.getElementById("results").addEventListener("click", (e) => {
    const tab = e.target.closest('.sample-tabs button, .combined-table tbody tr');
    if (tab) { showSample(parseInt(tab.dataset.sample, 10)); return; }
    const action = e.target.closest('[data-action]');
    if (action) {
        const si = parseInt(action.dataset.sample, 10);
        if (action.dataset.action === 'csv') downloadStatsCsv(si);
        if (action.dataset.action === 'unmapped') downloadUnmapped(si, action);
        if (action.dataset.action === 'consensus') downloadConsensus(si, action);
        if (action.dataset.action === 'bam') downloadBam(si, action);
        return;
    }
    const item = e.target.closest('.depth-item[data-index]');
    if (item) {
        const sm = samples[parseInt(item.dataset.sample, 10)];
        const contig = sm?.contigs[parseInt(item.dataset.index, 10)];
        if (contig) openReadsView(sm, contig);
    }
});

// --- READS VIEW (opens reads.html in a new tab) ---
const MAX_VIEW_READS = 3000; // above this, reads are randomly subsampled for display
const EXCLUDE_FLAGS = "0x704"; // unmapped, secondary, QC-fail, duplicate

const MAX_SEQ_BP = 300e6;        // longer contigs are shown without mismatches/consensus
const SEQ_CACHE_MAX_BP = 400e6;  // contig sequences kept for reuse (least recently used dropped first)

// One contig's sequence from the reference FASTA, uppercase, as { seq, reason }:
// seq is '' when it can't be used, with `reason` saying why (shown in the reads view).
// The file is streamed as bytes straight into a buffer of the contig's length (1 byte
// per base, no per-line strings), so long chromosomes fit; it stops after the contig.
// Header name = text after '>' up to the first whitespace, as minimap2/samtools use it.
async function readContigSeq(file, name, length) {
    if (!file) return { seq: '', reason: 'no reference file' };
    if (length > MAX_SEQ_BP) return { seq: '', reason: `contig longer than ${formatBp(MAX_SEQ_BP)}` };
    if (contigSeqCache.has(name)) { // move to the most recently used end
        const seq = contigSeqCache.get(name);
        contigSeqCache.delete(name); contigSeqCache.set(name, seq);
        return { seq, reason: '' };
    }
    const out = new Uint8Array(length);
    let n = 0, inTarget = false, found = false, done = false;
    let atLineStart = true, header = null; // header: { parts, len } while reading a '>' line
    const reader = file.stream().getReader();
    while (!done) {
        const { value: buf, done: eof } = await reader.read();
        if (eof) break;
        let i = 0;
        while (i < buf.length) {
            if (header) { // collect the header line (may span chunks)
                const nl = buf.indexOf(10, i), end = nl === -1 ? buf.length : nl;
                if (header.len < 4096) { header.parts.push(buf.slice(i, end)); header.len += end - i; }
                if (nl === -1) break;
                const bytes = new Uint8Array(header.len);
                let o = 0;
                for (const p of header.parts) { bytes.set(p.subarray(0, Math.min(p.length, bytes.length - o)), o); o += p.length; }
                inTarget = new TextDecoder('latin1').decode(bytes).split(/\s/)[0] === name;
                if (inTarget) found = true;
                header = null; atLineStart = true; i = nl + 1;
                continue;
            }
            if (atLineStart && buf[i] === 62) { // '>'
                if (found) { done = true; break; } // the next record: our contig is complete
                header = { parts: [], len: 0 }; i++;
                continue;
            }
            const nl = buf.indexOf(10, i), end = nl === -1 ? buf.length : nl;
            if (inTarget) {
                for (let k = i; k < end; k++) {
                    const b = buf[k];
                    if (b <= 32) continue; // \r, spaces
                    if (n < length) out[n] = b >= 97 && b <= 122 ? b - 32 : b; // uppercase
                    n++;
                }
            }
            atLineStart = nl !== -1;
            i = nl === -1 ? buf.length : nl + 1;
        }
    }
    reader.cancel().catch(() => {});
    if (!found) return { seq: '', reason: `“${name}” not found in the reference file` };
    if (n !== length) return { seq: '', reason: `length differs (${n.toLocaleString()} bp in the reference file, ${length.toLocaleString()} bp in the BAM)` };
    const seq = new TextDecoder('latin1').decode(out);
    contigSeqCache.set(name, seq);
    let total = 0;
    for (const v of contigSeqCache.values()) total += v.length;
    for (const k of contigSeqCache.keys()) { // drop least recently used, but keep this one
        if (total <= SEQ_CACHE_MAX_BP || k === name) break;
        total -= contigSeqCache.get(k).length;
        contigSeqCache.delete(k);
    }
    return { seq, reason: '' };
}

// Parse one SAM line into a compact read record. Only the first 10 fields are split,
// so the long QUAL column is never copied. With the contig sequence `ref`, read bases
// are compared to the reference to find mismatches.
//   blocks: aligned segments [s0,e0,s1,e1,...], split at deletions >= MIN_GAP_BP
//   mmPos/mmBase: mismatch positions and read bases
//   dels: every deletion as [pos,len,...]; ins: every insertion as [pos,len,...]
function parseSamRead(line, ref) {
    const f = line.split('\t', 10);
    if (f.length < 6 || f[5] === '*') return null;
    const flag = parseInt(f[1], 10);
    const start = parseInt(f[3], 10) - 1; // SAM is 1-based
    if (isNaN(flag) || isNaN(start)) return null;
    const seq = ref && f[9] && f[9] !== '*' ? f[9] : null;

    const blocks = [], dels = [], ins = [], mmPos = [];
    let mmBase = '';
    let pos = start, qpos = 0, blockStart = start, clipL = 0, clipR = 0;
    for (const [, lenStr, op] of f[5].matchAll(/(\d+)([MIDNSHP=X])/g)) {
        const len = parseInt(lenStr, 10);
        switch (op) {
            case 'M': case '=': case 'X':
                if (seq) {
                    for (let i = 0; i < len; i++) {
                        const rb = ref.charCodeAt(pos + i), qb = seq.charCodeAt(qpos + i) & 0xDF; // uppercase
                        if (qb !== rb && rb !== 78 && qb !== 78) { mmPos.push(pos + i); mmBase += String.fromCharCode(qb); }
                    }
                }
                pos += len; qpos += len;
                break;
            case 'I':
                ins.push(pos, len);
                qpos += len;
                break;
            case 'D': case 'N':
                if (op === 'D') dels.push(pos, len);
                if (len >= MIN_GAP_BP || op === 'N') { blocks.push(blockStart, pos); blockStart = pos + len; }
                pos += len;
                break;
            case 'S':
                // leading clip if nothing has been aligned yet, otherwise trailing
                if (pos === start && qpos === 0) clipL = len; else clipR = len;
                qpos += len;
                break;
        }
    }
    blocks.push(blockStart, pos);

    return {
        name: f[0], start, end: pos, mapq: parseInt(f[4], 10) || 0,
        reverse: (flag & 0x10) !== 0, supplementary: (flag & 0x800) !== 0,
        // read pairs: flag bits (0x1 paired, 0x2 proper, 0x8 mate unmapped, 0x20 mate reverse,
        // 0x40/0x80 first/second of pair); mate position on this contig (-1: elsewhere), or its contig
        flag, mateStart: f[6] === '=' ? parseInt(f[7], 10) - 1 : -1,
        mateRef: f[6] && f[6] !== '=' && f[6] !== '*' ? f[6] : null, tlen: parseInt(f[8], 10) || 0,
        blocks, dels: Int32Array.from(dels), ins: Int32Array.from(ins),
        mmPos: Int32Array.from(mmPos), mmBase, clipL, clipR
    };
}

async function countReads(bam, region) {
    const out = await CLI.exec("samtools", ["view", "-c", "-F", EXCLUDE_FLAGS, bam, region]);
    return parseInt(out.trim().split('\n').pop(), 10) || 0;
}

// Reads of a contig, or of a region of it ("name:start-end"), randomly subsampled to `limit`.
// `total`: the region's read count, if already known.
async function fetchContigReads(contig, bam, region = regionOf(contig.name), limit = MAX_VIEW_READS, total = null) {
    if (total === null) total = await countReads(bam, region);

    const args = ["view", "-F", EXCLUDE_FLAGS];
    let fraction = 1;
    if (total > limit) {
        fraction = limit / total;
        // samtools -s INT.FRAC: integer part is the seed, fractional part the fraction kept
        args.push("-s", `42${fraction.toFixed(6).slice(1)}`);
    }
    args.push(bam, region);
    const sam = await samtoolsText(args);

    // Reference bases for mismatch calling; refNote says why when they can't be used
    const { seq: ref, reason: refNote } = await readContigSeq(referenceFile, contig.name, contig.length);

    const reads = [];
    for (const line of sam.split('\n')) {
        const r = parseSamRead(line, ref || null);
        if (r) reads.push(r);
    }
    return { total, fraction, reads, ref, refNote };
}

// --- READS FOR A WINDOW: the reads view asks for all reads around what it shows when zoomed in ---
const MAX_WINDOW_READS = 20000;
// Per reads tab and request kind, one fetch at a time; the latest request wins.
//   reads-window: reads around the zoomed-in view, subsampled to MAX_WINDOW_READS
//   reads-scan:   one chunk of the variant scan, called here (scanChunk); the tab asks for the next when done
// Variant-scan chunks span about this many aligned bases (span x the contig's mean depth), so a
// chunk's SAM text stays well under browser string limits. Sized from the depth, not counted.
const SCAN_CHUNK_BASES = 12e6;
const SCAN_MIN_BP = 500;

// One chunk of the variant scan: all reads of [start, end) (end shortened to the chunk budget),
// parsed, filtered like the reads view (MAPQ, supplementary), counted and called here (pileup.js),
// so only the variant rows go back to the tab.
async function scanChunk(contig, bam, m) {
    const start = Math.max(0, Math.floor(m.start));
    const span = Math.max(SCAN_MIN_BP, Math.round(SCAN_CHUNK_BASES / Math.max(1, parseFloat(contig.meandepth) || 1)));
    const end = Math.min(contig.length, Math.ceil(m.end), start + span);
    const sam = await samtoolsText(["view", "-F", EXCLUDE_FLAGS, bam, regionOf(contig.name, start, end)]);
    const { seq: ref, reason } = await readContigSeq(referenceFile, contig.name, contig.length);
    if (!ref) throw new Error(`no reference sequence (${reason})`);
    const reads = [];
    for (const line of sam.split('\n')) {
        const r = parseSamRead(line, ref);
        if (r && r.mapq >= m.minMapq && (m.showSupp || !r.supplementary)) reads.push(r);
    }
    const rows = [];
    collectVariants(countBases(reads, ref, start, end), rows);
    return { start, end, rows };
}

const READ_REQUESTS = { 'reads-window': MAX_WINDOW_READS, 'reads-scan': Infinity };
const windowJobs = new Map(); // `${kind}` -> Map(tab -> { busy, next })

window.addEventListener('message', (e) => {
    const kind = e.data?.type;
    if (e.origin !== window.location.origin || !(kind in READ_REQUESTS)) return;
    if (!windowJobs.has(kind)) windowJobs.set(kind, new Map());
    const jobs = windowJobs.get(kind);
    const job = jobs.get(e.source) || { busy: false, next: null };
    jobs.set(e.source, job);
    job.next = e.data;
    if (!job.busy) runWindowJobs(e.source, job, kind);
});

async function runWindowJobs(win, job, kind) {
    job.busy = true;
    const type = `${kind}-data`;
    while (job.next) {
        const m = job.next;
        job.next = null;
        let reply;
        try {
            const sm = samples.find(x => x.name === m.sample) || (samples.length === 1 ? samples[0] : null);
            const contig = sm?.contigs.find(c => c.name === m.contig);
            if (!contig) throw new Error('contig not found');
            if (kind === 'reads-scan') {
                reply = { type, id: m.id, ...await scanChunk(contig, sm.bam, m) };
            } else {
                const start = Math.max(0, Math.floor(m.start)), end = Math.min(contig.length, Math.ceil(m.end));
                const { total, fraction, reads } = await fetchContigReads(contig, sm.bam, regionOf(contig.name, start, end), READ_REQUESTS[kind]);
                reply = { type, id: m.id, start, end, total, fraction, reads };
            }
        } catch (err) {
            reply = { type, id: m.id, error: String(err?.message || err) };
        }
        if (!win.closed) win.postMessage(reply, window.location.origin);
    }
    job.busy = false;
}

// One pending payload per opened tab, delivered when that tab reports it's ready
const readsViews = new Map(); // window -> { ready, payload }

function deliverReads(win) {
    const v = readsViews.get(win);
    if (v && v.ready && v.payload) {
        win.postMessage({ type: 'reads-data', ...v.payload }, window.location.origin);
        readsViews.delete(win);
    }
}

window.addEventListener('message', (e) => {
    if (e.origin !== window.location.origin || e.data?.type !== 'reads-ready') return;
    const v = readsViews.get(e.source);
    if (v) { v.ready = true; deliverReads(e.source); return; }
    // Unknown tab (e.g. the reads tab was reloaded): fetch again for the contig it names
    const sm = samples.find(x => x.name === e.data.sample) || (samples.length === 1 ? samples[0] : null);
    const contig = sm?.contigs.find(c => c.name === e.data.contig);
    if (contig) sendReadsTo(e.source, sm, contig, true);
});

async function openReadsView(sm, contig) {
    // Open the tab synchronously in the click handler so popup blockers allow it
    const hash = new URLSearchParams({ contig: contig.name, sample: sm.name }).toString();
    const win = window.open(`reads.html#${hash}`, '_blank');
    if (!win) {
        alert('Please allow pop-ups for this page to open the reads view.');
        return;
    }
    sendReadsTo(win, sm, contig, false);
}

async function sendReadsTo(win, sm, contig, ready) {
    readsViews.set(win, { ready, payload: null });
    let payload;
    try {
        const { total, fraction, reads, ref, refNote } = await fetchContigReads(contig, sm.bam);
        payload = {
            contig: { name: contig.name, length: contig.length, meandepth: contig.meandepth, coverage: contig.coverage, binSize: contig.binSize },
            sample: sm.name, profile: contig.profile, minProfile: contig.minProfile, total, fraction, reads, ref, refNote, limit: MAX_VIEW_READS,
            annotation: referenceAnnotation?.byContig.has(contig.name)
                ? { name: referenceAnnotation.name, features: referenceAnnotation.byContig.get(contig.name) } : null
        };
    } catch (err) {
        console.error(err);
        payload = { error: String(err?.message || err) };
    }
    const v = readsViews.get(win);
    if (v) { v.payload = payload; deliverReads(win); }
}

// Attach event listeners
document.getElementById("btn").addEventListener("click", run);
document.getElementById("reset-btn").addEventListener("click", resetPage);

// Colour the file inputs that have files chosen (also after a reload that kept the choice)
for (const id of ["ref", "query"]) {
    const el = document.getElementById(id);
    const mark = () => el.classList.toggle("has-files", el.files.length > 0);
    el.addEventListener("change", mark);
    mark();
}

// --- DEMO DATA: fetch the files in demo/ and put them in the inputs, as if chosen by hand ---
const DEMO = { ref: ["plasmids.gbk"], query: ["sample22.fastq.gz", "sample23.fastq.gz"] };
document.getElementById("demo-btn").addEventListener("click", async (e) => {
    const btn = e.currentTarget, label = btn.textContent;
    btn.disabled = true;
    btn.textContent = "Loading…";
    try {
        for (const [id, names] of Object.entries(DEMO)) {
            const dt = new DataTransfer();
            for (const name of names) {
                const res = await fetch(`demo/${name}`);
                if (!res.ok) throw new Error(`demo/${name}: ${res.status} ${res.statusText}`);
                dt.items.add(new File([await res.blob()], name));
            }
            const el = document.getElementById(id);
            el.files = dt.files;
            el.dispatchEvent(new Event("change"));
        }
        document.getElementById("preset").value = "map-ont";
        btn.textContent = "Demo data loaded";
    } catch (err) {
        btn.disabled = false;
        btn.textContent = label;
        alert(`Could not load the demo data (${err.message}). The page must be served over http(s), not opened as a file.`);
    }
});