#include <bitset>
#include <stdint.h>
#include <immintrin.h>


static inline uint64_t compute_shape_mask(uint32_t const shape) {
    uint64_t x = _pdep_u64(shape, 0x5555555555555555ULL);
    return x | (x << 1);
}

static inline constexpr uint32_t bit_length(uint32_t x) {
    return x == 0 ? 0 : 32 - __builtin_clz(x);
}

static inline constexpr uint32_t reverse32(uint32_t x) {
    x = ((x >> 1)  & 0x55555555) | ((x & 0x55555555) << 1);
    x = ((x >> 2)  & 0x33333333) | ((x & 0x33333333) << 2);
    x = ((x >> 4)  & 0x0F0F0F0F) | ((x & 0x0F0F0F0F) << 4);
    x = ((x >> 8)  & 0x00FF00FF) | ((x & 0x00FF00FF) << 8);
    x = (x >> 16) | (x << 16);
    return x;
}

static inline constexpr uint32_t reverse_shape(uint32_t x) {
    assert(x != 0);

    uint32_t len = 32 - __builtin_clz(x);
    return reverse32(x) >> (32 - len);
}


typedef struct {
    unsigned start;
    unsigned end;
    unsigned len;
} run_t;

static inline constexpr run_t find_long_run(uint32_t shape)
{
    run_t best{32, 32, 0};

    unsigned pos = 0;
    uint32_t x = shape;

    while (x) {
        unsigned zeros = std::countr_zero(x);
        pos += zeros;
        x >>= zeros;

        unsigned len = std::countr_zero(~x);

        if (len > best.len) {
            best.start = pos;
            best.len = len;
            best.end = pos + len;
        }

        pos += len;
        x >>= len;
    }

    return best;
}


static inline constexpr unsigned sum_runs(uint32_t shape, unsigned m) {
    unsigned sum = 0;
    while(shape) {
        unsigned start = __builtin_ctz(shape);
        shape >>= start;
        unsigned len = __builtin_ctz(~shape);
        if (len > m)
            sum += len;

        shape >>= len;
    }
    return sum;
}



typedef struct {
    uint32_t value;
    uint64_t mask;
    uint64_t mask_rev;
    size_t weight;
    size_t length;
    size_t kernel_start;
    size_t kernel_end;
} Shape32;

static inline Shape32 shape32_create(uint32_t value)
{
    Shape32 r;

    r.value = value;

    if(value != std::numeric_limits<uint32_t>::max()) {
        r.mask = compute_shape_mask(value);
        r.mask_rev = compute_shape_mask(reverse_shape(value));
        r.weight = __builtin_popcount(value);
        r.length = bit_length(value);
        run_t run = find_long_run(value);
        r.kernel_start = run.start;
        r.kernel_end = run.end;
    }

    // std::cout << "shape: " << std::bitset<32>(value) << "\n";
    // std::cout << "shape mask: " << std::bitset<64>(r.mask) << "\n";
    // std::cout << "shape mask rev: " << std::bitset<64>(r.mask_rev) << "\n";
    // std::cout << "weight: " << r.weight << "\n";
    // std::cout << "length: " << r.length << "\n";
    // std::cout << "kernel start: " << r.kernel_start << "\n";
    // std::cout << "kernel end: " << r.kernel_end << "\n";

    return r;
}
