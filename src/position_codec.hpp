// Lossless, independently seekable occurrence compression. All arithmetic on
// positions is modulo 2^64, including orientation and signed-delta zigzag; this
// preserves the full unsigned 8-byte position range without signed overflow.
#ifndef SXGC_POSITION_CODEC_HPP
#define SXGC_POSITION_CODEC_HPP
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

namespace sxgc_positions {
using Bytes = std::vector<unsigned char>;
constexpr unsigned MAX_ROWS = 4096;
enum Codec : unsigned { LEGACY = 0, COLUMNS = 1, TEMPLATES = 2 };
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
inline Bytes unpack(const Bytes &raw, unsigned codec, unsigned n, unsigned post, unsigned cont) {
    dimensions(n, post, cont);
    if (codec == LEGACY) return raw;
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
inline Bytes compress(const Bytes &raw) {
    uLongf len = compressBound(raw.size()); Bytes packed(len);
    if (compress2(packed.data(), &len, raw.data(), raw.size(), Z_BEST_SPEED) != Z_OK)
        fail("zlib compression failed");
    packed.resize(len); return packed;
}
struct Packed { Bytes raw, packed; unsigned codec; size_t legacy_bytes; };
inline Packed pack(const Bytes &raw, unsigned n, unsigned post, unsigned cont) {
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
