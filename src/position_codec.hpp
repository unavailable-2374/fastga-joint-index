// Lossless, independently seekable occurrence compression. All arithmetic on
// positions is modulo 2^64, including orientation and signed-delta zigzag; this
// preserves the full unsigned 8-byte position range without signed overflow.
#ifndef SXGC_POSITION_CODEC_HPP
#define SXGC_POSITION_CODEC_HPP
#include <algorithm>
#include <cstdint>
#include <memory>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>
#include <libdeflate.h>

namespace sxgc_positions {
using Bytes = std::vector<unsigned char>;
constexpr unsigned MAX_ROWS = 4096;
enum Codec : unsigned { LEGACY = 0, COLUMNS = 1, TEMPLATES = 2, PACKED = 3 };
[[noreturn]] inline void fail(const char *message) {
    throw std::runtime_error(std::string("position codec: ") + message);
}
inline void varint(Bytes &b, uint64_t v) {
    while (v >= 128) { b.push_back((v & 127) | 128); v >>= 7; }
    b.push_back(v);
}
inline uint64_t take(const Bytes &b, size_t &at) {
    uint64_t v = 0;
    for (unsigned shift = 0; shift < 64; shift += 7) {
        if (at == b.size()) fail("truncated varint");
        unsigned c = b[at++];
        if (shift == 63 && c > 1) fail("varint overflow");
        v |= uint64_t(c & 127) << shift;
        if (!(c & 128)) return v;
    }
    fail("invalid varint");
}
inline void require(bool condition, const char *message) {
    if (!condition) fail(message);
}
inline void dimensions(unsigned n, unsigned post, unsigned cont) {
    require(n && n <= MAX_ROWS && post && post <= 8 && cont && cont <= 8,
            "invalid dimensions");
}
inline uint64_t number(const unsigned char *p, unsigned width) {
    uint64_t v = 0;
    for (unsigned j = 0; j < width; ++j) v |= uint64_t(p[j]) << (8*j);
    return v;
}
inline uint64_t zigzag(uint64_t delta) {
    return (delta << 1) ^ (uint64_t(0) - (delta >> 63));
}
inline uint64_t unzigzag(uint64_t code) {
    return (code >> 1) ^ (uint64_t(0) - (code & 1));
}
// ---- PACKED codec: fixed-width bit packing of global oriented coordinates ----
// One W-bit value per row: strand * D + start[contig rank] + position. The contigs
// of one strand occupy consecutive ranges of length+1 values (a reverse row's
// position is its end, so 0..length are all possible) and D is their total.
// W = ceil(log2(2D)) follows from the contig lengths alone, so an archive only
// records those. Block bytes are the bit stream itself (little-endian, row i at
// bit i*W) with no entropy stage: positions in key order are a random permutation
// and zlib gains nothing on them (measured).
struct Domain {
    static constexpr unsigned MAX_WIDTH = 57;      // one unaligned 64-bit load per row
    unsigned width = 0, shift = 0;
    bool single = false;                           // no bucket holds two contig starts
    uint64_t strand = 0;                           // D: values per strand
    std::vector<uint64_t> start;                   // ncontig+1 offsets inside one strand
    std::vector<uint32_t> bucket;                  // contig rank holding value b<<shift
    // lengths are indexed by contig rank. Returns null when packing is impossible.
    static std::shared_ptr<const Domain> make(const std::vector<uint64_t> &lengths) {
        if (lengths.empty() || lengths.size() >= UINT32_MAX) return nullptr;
        auto d = std::make_shared<Domain>(); d->start.assign(lengths.size()+1, 0);
        for (size_t r = 0; r < lengths.size(); ++r) {
            if (lengths[r] >= (uint64_t(1) << 56) || d->start[r] >= (uint64_t(1) << 56)) return nullptr;
            d->start[r+1] = d->start[r] + lengths[r] + 1;
        }
        d->strand = d->start.back();
        d->width = 64 - __builtin_clzll(2*d->strand - 1);
        if (d->width > MAX_WIDTH) return nullptr;
        unsigned bits = 64 - __builtin_clzll(d->strand);
        d->shift = bits > 16 ? bits - 16 : 0;
        size_t buckets = ((d->strand - 1) >> d->shift) + 1, r = 0;
        d->bucket.resize(buckets+1);
        for (size_t b = 0; b < buckets; ++b) {
            uint64_t u = uint64_t(b) << d->shift;
            while (d->start[r+1] <= u) ++r;
            d->bucket[b] = r;
        }
        d->bucket[buckets] = lengths.size()-1;
        d->single = true;
        for (size_t b = 0; b < buckets; ++b) d->single &= d->bucket[b+1] - d->bucket[b] <= 1;
        return d;
    }
    unsigned ncontig() const { return start.size()-1; }
    // rank with start[rank] <= u < start[rank+1]; requires u < strand.
    unsigned rank_of(uint64_t u) const {
        size_t b = u >> shift; unsigned lo = bucket[b];
        if (single) return lo + (start[lo+1] <= u);   // branch-free: at most one start inside the bucket
        unsigned hi = bucket[b+1];
        if (hi - lo > 8) return std::upper_bound(start.begin()+lo+1, start.begin()+hi+1, u) - start.begin() - 1;
        while (lo < hi && start[lo+1] <= u) ++lo;
        return lo;
    }
};
inline uint64_t load64(const unsigned char *p) { uint64_t v; memcpy(&v, p, 8); return v; }
inline size_t packed_size(unsigned n, unsigned width) { return (size_t(n)*width + 7) / 8; }
// Packs canonical joint rows (zero id varint, mask and LCP). False when a row does
// not fit (nonzero bookkeeping fields or a coordinate outside the domain); the
// caller then uses the compressed codecs.
inline bool pack_coordinates(const Bytes &raw, unsigned n, unsigned post, unsigned cont,
                             const Domain &dom, Bytes &out) {
    dimensions(n, post, cont);
    const unsigned width = post+cont+3;
    if (raw.size() != size_t(n)*width) return false;
    out.assign(packed_size(n, dom.width) + 8, 0);
    unsigned __int128 acc = 0; unsigned nbits = 0; size_t at = 0;
    const uint64_t strand_bit = uint64_t(1) << (8*cont-1);
    const unsigned char *p = raw.data();
    for (unsigned i = 0; i < n; ++i, p += width) {
        if (p[0] | p[1] | p[2]) return false;
        uint64_t pos = number(p+3, post), rank = number(p+3+post, cont);
        bool reverse = rank & strand_bit; rank &= ~strand_bit;
        if (rank >= dom.ncontig() || pos >= dom.start[rank+1]-dom.start[rank]) return false;
        uint64_t v = (reverse ? dom.strand : 0) + dom.start[rank] + pos;
        acc |= (unsigned __int128)v << nbits; nbits += dom.width;
        if (nbits >= 64) { memcpy(out.data()+at, &acc, 8); at += 8; acc >>= 64; nbits -= 64; }
    }
    if (nbits) memcpy(out.data()+at, &acc, 8);
    out.resize(packed_size(n, dom.width));
    return true;
}
template<unsigned LEAD>
inline void decode_coordinates_impl(const unsigned char *in, size_t length, unsigned n, unsigned post,
                                    unsigned cont, const Domain &dom, unsigned char *out) {
    const unsigned stride = post+cont+2+LEAD, tail = 8/stride + 2;
    const unsigned fast = n > tail ? n - tail : 0;
    const uint64_t mask = (uint64_t(1) << dom.width) - 1, D = dom.strand;
    const uint64_t strand_bit = uint64_t(1) << (8*cont-1);
    const uint64_t *start = dom.start.data();
    const unsigned width = dom.width;
    uint64_t invalid = 0;
    unsigned i = 0;
    // Rows whose 8-byte load and 8-byte stores stay inside both buffers.
    const unsigned load_safe = length >= 8 ? unsigned(std::min<uint64_t>(n, ((length-8)*8)/width + 1)) : 0;
    const unsigned limit = std::min(fast, load_safe);
    for (; i < limit; ++i) {
        size_t bit = size_t(i)*width;
        uint64_t x = (load64(in+(bit>>3)) >> (bit & 7)) & mask;
        invalid |= (x >= 2*D);
        uint64_t reverse = x >= D, u = x - (reverse ? D : 0);
        unsigned rank = dom.rank_of(u < D ? u : D-1); // out-of-domain values are flagged and rejected below
        uint64_t pos = u - start[rank], rk = rank | (reverse ? strand_bit : 0);
        unsigned char *p = out + size_t(i)*stride;
        invalid |= (post < 8 && (pos >> (8*post)) != 0);
        memset(p, 0, LEAD+2);
        memcpy(p+LEAD+2, &pos, 8); memcpy(p+LEAD+2+post, &rk, 8);
    }
    for (; i < n; ++i) {
        size_t bit = size_t(i)*width, byte = bit >> 3; uint64_t x = 0;
        if (byte + 8 <= length) x = load64(in+byte); else memcpy(&x, in+byte, length-byte);
        x = (x >> (bit & 7)) & mask;
        invalid |= (x >= 2*D);
        uint64_t reverse = x >= D, u = x - (reverse ? D : 0);
        unsigned rank = dom.rank_of(u < D ? u : D-1); // out-of-domain values are flagged and rejected below
        uint64_t pos = u - start[rank], rk = rank | (reverse ? strand_bit : 0);
        unsigned char *p = out + size_t(i)*stride;
        invalid |= (post < 8 && (pos >> (8*post)) != 0);
        memset(p, 0, LEAD+2); p += LEAD;
        for (unsigned j = 0; j < post; ++j) p[2+j] = pos >> (8*j);
        for (unsigned j = 0; j < cont; ++j) p[2+post+j] = rk >> (8*j);
    }
    require(!invalid, "packed coordinate outside domain");
}
// Writes n rows of [lead zero bytes][mask 0][lcp 0][post][cont]; lead is 1 for the
// canonical form (zero id varint) and 0 for the reader's csize-wide rows. Whole
// 8-byte stores run while a row has 8 bytes of room behind it; the last rows use
// byte stores.
inline void decode_coordinates(const unsigned char *in, size_t length, unsigned n, unsigned post,
                               unsigned cont, const Domain &dom, unsigned char *out, unsigned lead) {
    dimensions(n, post, cont);
    require(length == packed_size(n, dom.width), "packed block length mismatch");
    require(lead <= 1, "invalid packed lead");
    if (lead) decode_coordinates_impl<1>(in, length, n, post, cont, dom, out);
    else decode_coordinates_impl<0>(in, length, n, post, cont, dom, out);
}
struct Rows {
    std::vector<uint64_t> ids;
    Bytes payload;
};
inline Rows parse(const Bytes &raw, unsigned n, unsigned post, unsigned cont) {
    dimensions(n, post, cont);
    unsigned width = post + cont + 2;
    Rows rows; rows.ids.reserve(n); rows.payload.reserve(size_t(n)*width);
    size_t at = 0; uint64_t id = 0;
    for (unsigned i = 0; i < n; ++i) {
        uint64_t delta = take(raw, at);
        require(!i || delta <= UINT64_MAX-id, "key ID overflow");
        id = i ? id+delta : delta; rows.ids.push_back(id);
        require(raw.size()-at >= width, "short payload");
        rows.payload.insert(rows.payload.end(), raw.begin()+at, raw.begin()+at+width);
        at += width;
    }
    require(at == raw.size(), "trailing payload");
    return rows;
}
inline Bytes canonical(const Rows &rows, unsigned width) {
    Bytes raw; raw.reserve(rows.payload.size()+rows.ids.size());
    uint64_t previous = 0;
    for (size_t i = 0; i < rows.ids.size(); ++i) {
        require(!i || rows.ids[i] >= previous, "decreasing key IDs");
        varint(raw, i ? rows.ids[i]-previous : rows.ids[i]); previous = rows.ids[i];
        raw.insert(raw.end(), rows.payload.begin()+i*width, rows.payload.begin()+(i+1)*width);
    }
    return raw;
}
inline Bytes columns(const Rows &rows, unsigned width) {
    Bytes out; uint64_t previous = 0;
    for (size_t i = 0; i < rows.ids.size(); ++i) {
        varint(out, i ? rows.ids[i]-previous : rows.ids[i]); previous = rows.ids[i];
    }
    for (unsigned j = 0; j < width; ++j)
        for (size_t i = 0; i < rows.ids.size(); ++i) out.push_back(rows.payload[i*width+j]);
    return out;
}
inline Bytes templates(const Rows &rows, unsigned post, unsigned cont) {
    const unsigned width = post+cont+2, n = rows.ids.size();
    std::vector<unsigned> starts, lengths, refs;
    std::vector<Bytes> shapes;
    std::unordered_map<std::string, unsigned> known;
    for (unsigned first = 0; first < n;) {
        unsigned end = first+1;
        while (end < n && rows.ids[end] == rows.ids[first]) ++end;
        starts.push_back(first); lengths.push_back(end-first);
        auto p = rows.payload.data()+size_t(first)*width;
        uint64_t anchor = number(p+2, post);
        if (p[width-1] & 128) anchor = uint64_t(0)-anchor;
        Bytes shape; varint(shape, end-first);
        // The first occurrence has relative position zero, but its complete
        // contig/strand identity still belongs to the shared shape.
        for (unsigned i = first; i < end; ++i) {
            auto q = rows.payload.data()+size_t(i)*width;
            shape.insert(shape.end(), q+2+post, q+width);
        }
        for (unsigned i = first+1; i < end; ++i) {
            auto q = rows.payload.data()+size_t(i)*width;
            uint64_t pos = number(q+2, post);
            if (q[width-1] & 128) pos = uint64_t(0)-pos;
            varint(shape, zigzag(pos-anchor));
        }
        std::string token(reinterpret_cast<const char *>(shape.data()), shape.size());
        auto inserted = known.emplace(std::move(token), shapes.size());
        refs.push_back(inserted.first->second);
        if (inserted.second) shapes.push_back(std::move(shape));
        first = end;
    }
    Bytes out; varint(out, starts.size()); uint64_t previous = 0;
    for (size_t g = 0; g < starts.size(); ++g) {
        uint64_t id = rows.ids[starts[g]];
        varint(out, g ? id-previous : id); varint(out, lengths[g]); previous = id;
    }
    // Masks and LCPs remain lossless, independent of coordinate templates.
    for (unsigned j = 0; j < 2; ++j)
        for (unsigned i = 0; i < n; ++i) out.push_back(rows.payload[size_t(i)*width+j]);
    // Byte planes expose constant high position bytes to the compressor.
    for (unsigned j = 0; j < post; ++j)
        for (unsigned first : starts) out.push_back(rows.payload[size_t(first)*width+2+j]);
    varint(out, shapes.size());
    for (const auto &shape : shapes) {
        varint(out, shape.size()); out.insert(out.end(), shape.begin(), shape.end());
    }
    for (unsigned ref : refs) varint(out, ref);
    return out;
}
inline Bytes unpack(const Bytes &raw, unsigned codec, unsigned n, unsigned post, unsigned cont,
                    const Domain *domain = nullptr) {
    dimensions(n, post, cont);
    if (codec == LEGACY) return raw;
    if (codec == PACKED) {
        require(domain, "packed block without a coordinate domain");
        Bytes out(size_t(n)*(post+cont+3));
        decode_coordinates(raw.data(), raw.size(), n, post, cont, *domain, out.data(), 1);
        return out;
    }
    unsigned width = post+cont+2;
    Rows rows; rows.ids.reserve(n); rows.payload.resize(size_t(n)*width);
    size_t at = 0;
    if (codec == COLUMNS) {
        uint64_t previous = 0;
        for (unsigned i = 0; i < n; ++i) {
            uint64_t delta = take(raw, at);
            require(!i || delta <= UINT64_MAX-previous, "column key ID overflow");
            previous = i ? previous+delta : delta; rows.ids.push_back(previous);
        }
        require(raw.size()-at == rows.payload.size(), "column payload length mismatch");
        for (unsigned j = 0; j < width; ++j)
            for (unsigned i = 0; i < n; ++i) rows.payload[size_t(i)*width+j] = raw[at++];
        return canonical(rows, width);
    }
    require(codec == TEMPLATES, "unknown codec");
    uint64_t groups = take(raw, at);
    require(groups && groups <= n, "invalid group count");
    std::vector<unsigned> lengths; uint64_t previous = 0;
    for (unsigned g = 0; g < groups; ++g) {
        uint64_t delta = take(raw, at), count = take(raw, at);
        require(!g || (delta && delta <= UINT64_MAX-previous), "invalid group key delta");
        previous = g ? previous+delta : delta;
        require(count && count <= n-rows.ids.size(), "invalid group length");
        lengths.push_back(count); rows.ids.insert(rows.ids.end(), count, previous);
    }
    require(rows.ids.size() == n, "group row count mismatch");
    require(raw.size()-at >= size_t(2)*n+groups*post, "truncated coordinate columns");
    for (unsigned j = 0; j < 2; ++j)
        for (unsigned i = 0; i < n; ++i) rows.payload[size_t(i)*width+j] = raw[at++];
    std::vector<uint64_t> anchors(groups);
    for (unsigned j = 0; j < post; ++j)
        for (unsigned g = 0; g < groups; ++g) anchors[g] |= uint64_t(raw[at++]) << (8*j);
    uint64_t shape_count = take(raw, at);
    require(shape_count && shape_count <= groups, "invalid shape count");
    struct Shape { unsigned count; Bytes flags; std::vector<uint64_t> deltas; };
    std::vector<Shape> shapes; unsigned shape_rows = 0;
    for (unsigned s = 0; s < shape_count; ++s) {
        uint64_t bytes = take(raw, at);
        require(bytes && bytes <= raw.size()-at, "truncated shape");
        Bytes encoded(raw.begin()+at, raw.begin()+at+bytes); at += bytes;
        size_t cursor = 0; uint64_t count = take(encoded, cursor);
        require(count && count <= n-shape_rows, "invalid shape rows"); shape_rows += count;
        require(encoded.size()-cursor >= count*cont, "truncated shape contigs");
        Shape shape{unsigned(count), Bytes(encoded.begin()+cursor, encoded.begin()+cursor+count*cont), {0}};
        cursor += count*cont; shape.deltas.reserve(count);
        for (unsigned i = 1; i < count; ++i) shape.deltas.push_back(unzigzag(take(encoded, cursor)));
        require(cursor == encoded.size(), "trailing shape data"); shapes.push_back(std::move(shape));
    }
    unsigned row = 0;
    for (unsigned g = 0; g < groups; ++g) {
        uint64_t ref = take(raw, at); require(ref < shapes.size(), "shape reference out of range");
        const auto &shape = shapes[ref]; require(shape.count == lengths[g], "shape length mismatch");
        uint64_t anchor = anchors[g];
        if (shape.flags[cont-1] & 128) anchor = uint64_t(0)-anchor;
        for (unsigned i = 0; i < shape.count; ++i, ++row) {
            auto p = rows.payload.data()+size_t(row)*width;
            memcpy(p+2+post, shape.flags.data()+size_t(i)*cont, cont);
            uint64_t pos = anchor+shape.deltas[i];
            if (p[width-1] & 128) pos = uint64_t(0)-pos;
            require(post == 8 || pos < (uint64_t(1) << (8*post)), "reconstructed position overflow");
            for (unsigned j = 0; j < post; ++j) p[2+j] = pos >> (8*j);
        }
    }
    require(at == raw.size(), "trailing template data");
    return canonical(rows, width);
}
// Blocks are zlib-format streams written by libdeflate; any zlib reader decodes
// them. Each thread keeps its own compressor and decompressor. Level 1: on chr20
// levels 2 and 4 shrank the index 4-5% but slowed block encoding 25-39%.
#define SXGC_DEFLATE_LEVEL 1
inline Bytes compress(const Bytes &raw) {
    struct Compressor { libdeflate_compressor *c=libdeflate_alloc_compressor(SXGC_DEFLATE_LEVEL); ~Compressor(){libdeflate_free_compressor(c);} };
    thread_local Compressor z;
    if (!z.c) fail("libdeflate compressor allocation failed");
    Bytes packed(libdeflate_zlib_compress_bound(z.c, raw.size()));
    size_t len = libdeflate_zlib_compress(z.c, raw.data(), raw.size(), packed.data(), packed.size());
    if (!len) fail("libdeflate compression failed");
    packed.resize(len); return packed;
}
// Inflates a zlib stream whose decoded size is known exactly; false on any mismatch.
inline bool inflate_exact(unsigned char *out, size_t size, const unsigned char *in, size_t length) {
    struct Decompressor { libdeflate_decompressor *d=libdeflate_alloc_decompressor(); ~Decompressor(){libdeflate_free_decompressor(d);} };
    thread_local Decompressor z;
    if (!z.d) fail("libdeflate decompressor allocation failed");
    return libdeflate_zlib_decompress(z.d, in, length, out, size, nullptr) == LIBDEFLATE_SUCCESS;
}
struct Packed { Bytes raw, packed; unsigned codec; size_t legacy_bytes; };
// With a domain, joint rows are stored bit-packed: codec PACKED, Packed::raw holds the
// stored bytes and Packed::packed stays empty. Other rows use the compressed codecs.
inline Packed pack(const Bytes &raw, unsigned n, unsigned post, unsigned cont,
                   const Domain *domain = nullptr) {
    if (domain) {
        Bytes bits;
        if (pack_coordinates(raw, n, post, cont, *domain, bits)) return Packed{std::move(bits), {}, PACKED, 0};
    }
    auto rows = parse(raw, n, post, cont);
    Packed best{raw, compress(raw), LEGACY, 0}; best.legacy_bytes = best.packed.size();
    for (unsigned codec : {unsigned(COLUMNS), unsigned(TEMPLATES)}) {
        Bytes candidate = codec == COLUMNS ? columns(rows, post+cont+2) :
                          templates(rows, post, cont);
        Bytes packed = compress(candidate);
        if (packed.size() < best.packed.size()) {
            best.raw = std::move(candidate); best.packed = std::move(packed); best.codec = codec;
        }
    }
    return best;
}
} // namespace sxgc_positions
#endif
