// Helpers shared by the main page (script.js) and the reads view (reads.html). Plain script: defines globals.

// Read parsing (main page) and drawing (reads view) must agree on these
const MIN_GAP_BP = 20; // deletions/introns at least this long split a read into blocks (drawn as gaps)
const MIN_INS_BP = 20; // insertions at least this long are marked when zoomed out

// Base pairs, compactly: 950 bp, 12.5 kb, 3 Mb
function formatBp(bp) {
    if (bp >= 1e6) return `${+(bp / 1e6).toFixed(2)} Mb`;
    if (bp >= 1e3) return `${+(bp / 1e3).toFixed(1)} kb`;
    return `${Math.round(bp)} bp`;
}

// "Nice" tick step (1, 2 or 5 x 10^k) giving roughly `target` ticks over `span`
function niceStep(span, target) {
    const raw = span / Math.max(1, target);
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const n = raw / mag;
    return (n < 1.5 ? 1 : n < 3.5 ? 2 : n < 7.5 ? 5 : 10) * mag;
}

// Elapsed time, compactly: 45s, 3m 05s, 1h 02m
function formatDuration(ms) {
    const t = Math.max(0, Math.round(ms / 1000));
    const h = Math.floor(t / 3600), m = Math.floor(t / 60) % 60, sec = t % 60;
    const pad = (n) => String(n).padStart(2, '0');
    if (h) return `${h}h ${pad(m)}m`;
    if (m) return `${m}m ${pad(sec)}s`;
    return `${sec}s`;
}

// Text for HTML content and attribute values
const escapeHtml = (v) => String(v).replace(/[&<>"']/g, ch => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[ch]));

// Save a Blob (or text) as a downloaded file
function downloadBlob(blob, filename) {
    const url = URL.createObjectURL(blob);
    const a = Object.assign(document.createElement('a'), { href: url, download: filename });
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 10000);
}
const downloadText = (text, filename, type = 'text/plain') => downloadBlob(new Blob([text], { type }), filename);

// ---- instant tooltips (the browser's own title tooltips appear only after ~1 s) ----
// Place a tooltip element next to the mouse, kept inside the window
function placeTip(tip, e) {
    tip.style.display = 'block';
    tip.style.left = `${Math.min(e.clientX + 12, window.innerWidth - tip.offsetWidth - 8)}px`;
    tip.style.top = `${e.clientY + 14}px`;
}
// Show `tip` with the data-tip text of the element under the mouse, inside `containers`.
// `className` is added while it shows a data-tip. Returns the hide function.
function attachQuickTips(containers, tip, { className = '' } = {}) {
    const hide = () => { tip.style.display = 'none'; if (className) tip.classList.remove(className); };
    const move = (e) => {
        const el = e.target.closest('[data-tip]');
        if (!el) return hide();
        tip.textContent = el.dataset.tip;
        if (className) tip.classList.add(className);
        placeTip(tip, e);
    };
    for (const c of containers) {
        c.addEventListener('mousemove', move);
        c.addEventListener('mouseleave', hide);
    }
    return hide;
}

// ---- gzip ----
// gzip -> bytes, as a TransformStream. Unlike the browser's DecompressionStream, which stops after
// the first gzip member ("Junk found after end of compressed data"), this reads files made of
// several members: bgzip, nanopore basecaller output, `cat a.gz b.gz`.
// The gzip framing (header, 8-byte trailer, zero padding between members) is read here, and only
// the compressed data of each member is inflated, by pako in raw mode (pako is loaded by the
// page). pako's own multi-member handling is not used: it breaks on some real files ("invalid
// distance too far back" right after a member boundary).
function gunzipStream() {
    let ctrl, inf = null, state = 'header', trailerLeft = 0, pending = null, members = 0;
    const fail = (msg) => { throw new Error(`could not unpack the gzip data (${msg})`); };
    // length of the gzip member header at the start of b, -1 if b doesn't hold all of it yet
    const headerLength = (b) => {
        if (b.length < 10) return -1;
        if (b[0] !== 0x1f || b[1] !== 0x8b || b[2] !== 8) fail(members ? 'bad data between gzip members' : 'not a gzip file');
        const flags = b[3];
        let o = 10;
        if (flags & 4) { if (b.length < o + 2) return -1; o += 2 + (b[o] | b[o + 1] << 8); } // FEXTRA (e.g. bgzip)
        for (const bit of [8, 16]) if (flags & bit) { const z = b.indexOf(0, o); if (z < 0) return -1; o = z + 1; } // FNAME, FCOMMENT
        if (flags & 2) o += 2; // FHCRC
        return b.length >= o ? o : -1;
    };
    return new TransformStream({
        start(c) { ctrl = c; },
        transform(chunk) {
            let data = chunk;
            if (pending) { // an incomplete header from the previous chunk
                data = new Uint8Array(pending.length + chunk.length);
                data.set(pending); data.set(chunk, pending.length);
                pending = null;
            }
            while (data.length) {
                if (state === 'trailer') { // CRC32 + size of the member just inflated
                    const n = Math.min(trailerLeft, data.length);
                    trailerLeft -= n; data = data.subarray(n);
                    if (trailerLeft) break;
                    state = 'header';
                } else if (state === 'header') {
                    let i = 0;
                    while (i < data.length && data[i] === 0) i++; // zero padding between members
                    data = data.subarray(i);
                    if (!data.length) break;
                    const h = headerLength(data);
                    if (h < 0) { pending = data.slice(); break; }
                    data = data.subarray(h);
                    inf = new pako.Inflate({ raw: true });
                    inf.onData = (d) => ctrl.enqueue(d);
                    state = 'body'; members++;
                } else {
                    inf.push(data, false);
                    if (inf.err) fail(inf.msg || 'corrupt data');
                    if (!inf.ended) break;
                    data = data.subarray(data.length - inf.strm.avail_in); // input after this member's data
                    state = 'trailer'; trailerLeft = 8;
                }
            }
        },
        flush() {
            if (!members && !pending) fail('empty or not a gzip file');
            if (state !== 'header' || pending) throw new Error('the gzip file ends early (truncated?)');
        }
    });
}
