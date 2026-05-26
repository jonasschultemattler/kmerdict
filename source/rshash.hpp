#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <sux/bits/SimpleSelect.hpp>
#include <gtl/phmap.hpp>
#include "compact_vector.hpp"
#include "EliasFano.hpp"
#include "minimiser_views.hpp"
#include "util.hpp"

using namespace seqan3::literals;
using namespace seqan3::contrib::sdsl;


const uint64_t seed1 = 1;
const uint64_t seed2 = 0x29'6D'BD'33'32'56'8C'64;
const uint64_t seed3 = 0xE5'9A'38'5F'03'76'C9'F6;


struct SkmerInfo {
    size_t seq_id;
    size_t start;
    size_t end;
};

// struct MinimizerInfo32 {
//     uint64_t minimizer_value;
//     uint32_t position;

//     bool compare(const MinimizerInfo32 &x, const MinimizerInfo32 &y) {
//     	return x.minimizer_value < y.minimizer_value;
//     }
// };

// struct MinimizerInfo64 {
//     uint64_t minimizer_value;
//     uint64_t position;

//     bool compare(const MinimizerInfo64 &x, const MinimizerInfo64 &y) {
//     	return x.minimizer_value < y.minimizer_value;
//     }
// };

typedef std::pair<uint64_t, uint32_t> MinimizerInfo32;
typedef std::pair<uint64_t, uint64_t> MinimizerInfo64;

struct RadixTraitsMinimizer32 {
    static const int nBytes = 12;
    int kth_byte(const MinimizerInfo32 &x, int k) {
    	if (k >= 8) return x.first >> ((k - 8) * 8) & 0xFF;
    	return x.second >> (k * 8) & 0xFF;
    }
    bool compare(const MinimizerInfo32 &x, const MinimizerInfo32 &y) {
    	return x.first < y.first;
    }
};

struct RadixTraitsMinimizer64 {
    static const int nBytes = 16;
    int kth_byte(const MinimizerInfo64 &x, int k) {
    	if (k >= 8) return x.first >> ((k - 8) * 8) & 0xFF;
    	return x.second >> (k * 8) & 0xFF;
    }
    bool compare(const MinimizerInfo64 &x, const MinimizerInfo64 &y) {
    	return x.first < y.first;
    }
};



#ifndef RSHASH_HPP
#define RSHASH_HPP


class RSHash
{
private:
    uint64_t level;
    uint64_t k, m1, m_thres1, m2, m_thres2, m3, m_thres3, threshold;
    uint64_t span1, span2, span3;
    uint64_t kmermask, mmermask1, mmermask2, mmermask3;
    uint32_t shape;
    bool loc;
    mixer_64 m_hasher1, m_hasher2, m_hasher3;
    uint64_t no_text_kmers;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> r1, r2, r3;
    bit_vector s1, s2, s3;
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s1_select, s2_select, s3_select;
    bits::compact_vector offsets1, offsets2, offsets3;
    // std::vector<uint32_t> offsets1, offsets2, offsets3;
    gtl::flat_hash_set<uint64_t> hashset;
    // gtl::flat_hash_map<uint64_t, uint64_t> hashmap;
    // gtl::flat_hash_map<uint64_t, bits::compact_vector> hashmap;
    gtl::flat_hash_map<uint64_t, std::vector<uint32_t>> hashmap;
    sux::bits::EliasFano<sux::util::AllocType::MALLOC> endpoints;
    std::vector<uint64_t> text;
    template<int level, typename MinimizerT>
    inline void filter_freq_minimizers(std::vector<MinimizerT> &minimizers,
    std::vector<uint8_t> &counts, size_t &no_minimizers, size_t &no_skmers);
    template<int level, typename MinimizerT>
    inline uint64_t get_minimizers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &, const std::vector<SkmerInfo> &, std::vector<uint64_t> &, std::vector<uint8_t> &);
    template<int level>
    void mark_minimizer_occurences(const size_t, const std::vector<uint8_t> &);
    template<int level>
    void fill_minimizer_offsets(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &, std::vector<size_t> &, std::vector<uint8_t> &, const size_t, const size_t);
    template<int level>
    size_t get_frequent_skmers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &, const std::vector<SkmerInfo> &, std::vector<SkmerInfo> &);
    void fill_ht(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>&, const std::vector<SkmerInfo> &);
    template<int level>
    inline uint64_t find_minimiser(const uint64_t, const uint64_t, size_t &, size_t &);
    template<int level>
    inline void update_minimiser(const uint64_t, const uint64_t, uint64_t&, size_t &, size_t &);
    template<int level>
    inline bool check(const uint64_t, const uint64_t, uint64_t*, const size_t, const size_t, const size_t, const size_t);
    template<int level>
    inline void fill_buffer(uint64_t *, uint64_t *, size_t, size_t, const uint64_t);
    template<int level>
    inline bool check_overlap(uint64_t, uint64_t, uint64_t &, uint64_t &);
    template<int level>
    inline bool check_minimiser_pos(uint64_t *, const uint64_t, const uint64_t, const uint64_t, const size_t, const size_t, const size_t, bool &, uint64_t &, uint64_t &, uint64_t &);
    template<int level>
    inline bool check_minimiser_pos2(uint64_t *, const uint64_t, const uint64_t, const uint64_t, const size_t, const size_t, const size_t, const size_t, bool &, uint64_t &, uint64_t &, uint64_t &);
    template<int level>
    inline bool lookup_buffer(uint64_t *, uint64_t *, const size_t, const uint64_t,  const uint64_t, uint64_t &, const size_t, const size_t, bool &, uint64_t &, uint64_t &);
    inline bool extend_in_text(uint64_t&, uint64_t, uint64_t, bool, const uint64_t, const uint64_t);
    const inline uint64_t get_word64(uint64_t pos);
    const inline uint64_t get_base(uint64_t pos);
    uint64_t streaming_lookup1(const seqan3::bitpacked_sequence<seqan3::dna4>&, uint64_t&);
    uint64_t streaming_lookup2(const seqan3::bitpacked_sequence<seqan3::dna4>&, uint64_t&);
    uint64_t streaming_lookup3(const seqan3::bitpacked_sequence<seqan3::dna4>&, uint64_t&);
    size_t streaming_locate1(const seqan3::bitpacked_sequence<seqan3::dna4>&, std::vector<std::pair<uint64_t, bool>> &, size_t &);
    size_t streaming_locate2(const seqan3::bitpacked_sequence<seqan3::dna4>&, std::vector<std::pair<uint64_t, bool>> &, size_t &);
    size_t streaming_locate3(const seqan3::bitpacked_sequence<seqan3::dna4>&, std::vector<std::pair<uint64_t, bool>> &, size_t &);
    uint64_t lookup1(const std::vector<uint64_t>&);
    uint64_t lookup2(const std::vector<uint64_t>&);
    uint64_t lookup3(const std::vector<uint64_t>&);
    template<int level>
    inline bool report_minimiser_pos(uint64_t *, const uint64_t, const uint64_t, const uint64_t, size_t &, const size_t, const size_t, std::vector<std::pair<uint64_t, bool>> &);
    template<int level>
    inline bool report_minimiser_pos2(uint64_t *, const uint64_t, const uint64_t, const uint64_t, size_t &, const size_t, const size_t, const size_t, std::vector<std::pair<uint64_t, bool>> &);
    template<int level>
    inline void locate_buffer(uint64_t*, uint64_t*, const size_t, const uint64_t, const uint64_t, const size_t, const size_t, std::vector<std::pair<uint64_t, bool>> &, size_t &, size_t &);



public:
    RSHash() : endpoints(std::vector<uint64_t>{}, 1),
        r1(std::vector<uint64_t>{}, 1),
        r2(std::vector<uint64_t>{}, 1),
        r3(std::vector<uint64_t>{}, 1),
        m_hasher1(seed1), m_hasher2(seed2), m_hasher3(seed3)
    {}
    RSHash(uint8_t const k, uint8_t const level, uint8_t const m1, uint8_t const m2, uint8_t const m3,
            uint8_t const m_thres1, uint8_t const m_thres2, uint8_t const m_thres3, uint16_t const threshold, bool const loc, uint32_t const shape)
        : k(k), level(level), m1(m1), m2(m2), m3(m3),
        m_thres1(m_thres1), m_thres2(m_thres2), m_thres3(m_thres3), threshold(threshold), loc(loc), shape(shape),
        span1(k-m1+1), span2(k-m2+1), span3(k-m3+1),
        endpoints(std::vector<uint64_t>{}, 1),
        r1(std::vector<uint64_t>{}, 1),
        r2(std::vector<uint64_t>{}, 1),
        r3(std::vector<uint64_t>{}, 1),
        m_hasher1(seed1), m_hasher2(seed2), m_hasher3(seed3),
        kmermask(compute_mask(2u * k)),
        mmermask1(compute_mask(2u * m1)),
        mmermask2(compute_mask(2u * m2)),
        mmermask3(compute_mask(2u * m3))
    {}
    uint8_t getk() { return k; }
    bool has_locate() { return loc; }
    uint8_t getmaxmthres() { return std::max({m_thres1, m_thres2, m_thres3}); }
    uint64_t number_unitigs() { return endpoints.rank(endpoints.size()); }
    size_t unitig_size(uint64_t unitig_id) { return endpoints.select(unitig_id+1) - endpoints.select(unitig_id) - k + 1; }
    const inline uint64_t access(const uint64_t, const size_t);
    uint64_t lookup(const std::vector<uint64_t>&);
    void build(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>>&);
    uint64_t streaming_lookup(const seqan3::bitpacked_sequence<seqan3::dna4>&, uint64_t&);
    size_t streaming_locate(const seqan3::bitpacked_sequence<seqan3::dna4>&, std::vector<std::pair<uint64_t, bool>> &, size_t &);
    int save(const std::filesystem::path&);
    int load(const std::filesystem::path&);
    void print_info() {
        const size_t N = text.size()*32;
        const size_t offset_width = std::bit_width(N);
        const uint64_t M1 = 1ULL << (m1+m1);
        const uint64_t M2 = 1ULL << (m2+m2);
        const uint64_t M3 = 1ULL << (m3+m3);
        const uint64_t no_minimizers1 = r1.rank(M1);
        const uint64_t no_skmers1 = s1.size();
        const uint64_t no_minimizers2 = r2.rank(M2);
        const uint64_t no_skmers2 = s2.size();
        const uint64_t no_minimizers3 = r3.rank(M3);
        const uint64_t no_skmers3 = s3.size();
        size_t ht_space;
        if(loc) {
            // ht_space = hashmap.capacity()*(sizeof(uint64_t) + sizeof(uint32_t) + 1)*8;
            size_t total = sizeof(hashmap) +  hashmap.capacity()*sizeof(decltype(hashmap)::value_type);
            for (const auto& [k, vec] : hashmap)
                total += vec.capacity() * sizeof(uint32_t);
            ht_space = total*8;
        }
        else
            ht_space = hashset.capacity()*(sizeof(uint64_t) + 1)*8;

        std::cout << "====== report ======\n";
        std::cout << "text length: " << N << "\n";
        std::cout << "textkmers: " << no_text_kmers <<  '\n';
    
        std::cout << "no minimiser1: " << no_minimizers1 << "\n";
        std::cout << "no distinct minimiser1: " << no_skmers1 << "\n";
        std::cout << "avg superkmers1: " << (double) no_skmers1/no_minimizers1 <<  '\n';
        std::cout << "no minimiser2: " << no_minimizers2 << "\n";
        std::cout << "no distinct minimiser2: " << no_skmers2 << "\n";
        std::cout << "avg superkmers2: " << (double) no_skmers2/no_minimizers2 <<  '\n';
        std::cout << "no minimiser3: " << no_minimizers3 << "\n";
        std::cout << "no distinct minimiser3: " << no_skmers3 << "\n";
        std::cout << "avg superkmers3: " << (double) no_skmers3/no_minimizers3 <<  '\n';
        std::cout << "no kmers HT: " << hashmap.size() << " " << (double) hashmap.size()/no_text_kmers*100 << "%\n";

        std::cout << "density r1: " << (double) no_minimizers1/M1*100 << "%\n";
        std::cout << "density r2: " << (double) no_minimizers2/M2*100 << "%\n";
        std::cout << "density r3: " << (double) no_minimizers3/M3*100 << "%\n";
        std::cout << "density s1: " << (double) no_minimizers1/(no_skmers1+1)*100 <<  "%\n";
        std::cout << "density s2: " << (double) no_minimizers2/(no_skmers2+1)*100 <<  "%\n";
        std::cout << "density s3: " << (double) no_minimizers3/(no_skmers3+1)*100 <<  "%\n";
        std::cout << "\nspace per kmer in bit:\n";
        std::cout << "text: " << (double) 2*N/no_text_kmers << "\n";
        std::cout << "endpoints: " << (double) endpoints.bitCount()/no_text_kmers << "\n";
        std::cout << "offsets1: " << (double) no_skmers1*offset_width/no_text_kmers << "\n";
        std::cout << "offsets2: " << (double) no_skmers2*offset_width/no_text_kmers << "\n";
        std::cout << "offsets3: " << (double) no_skmers3*offset_width/no_text_kmers << "\n";
        std::cout << "Hashtable: " << (double) ht_space/no_text_kmers << "\n";
        std::cout << "R_1: " << (double) r1.bitCount()/no_text_kmers << "\n";
        std::cout << "R_2: " << (double) r2.bitCount()/no_text_kmers << "\n";
        std::cout << "R_3: " << (double) r3.bitCount()/no_text_kmers << "\n";
        std::cout << "S_1: " << (double) (no_skmers1+1)/no_text_kmers << "\n";
        std::cout << "S_2: " << (double) (no_skmers2+1)/no_text_kmers << "\n";
        std::cout << "S_3: " << (double) (no_skmers3+1)/no_text_kmers << "\n";
    
        std::cout << "total: " << (double) (no_skmers1*offset_width+no_skmers2*offset_width+no_skmers3*offset_width+2*N+r1.bitCount()+r2.bitCount()+r3.bitCount()+no_skmers1+1+s1_select.bitCount()+no_skmers2+1+s2_select.bitCount()+no_skmers3+1+s3_select.bitCount()+endpoints.bitCount()+ht_space)/no_text_kmers << "\n";
    }
    std::vector<uint64_t> rand_text_kmers(const uint64_t n)
    {
        std::uniform_int_distribution<uint32_t> distr;
        std::mt19937 m_rand(1);
        std::vector<std::uint64_t> kmers;
        kmers.reserve(n);
        const uint64_t l = (text.size()-1)*32;

        const uint64_t no_unitigs = number_unitigs();
        for (uint64_t i = 0; i < n;) {
            const uint64_t offset = distr(m_rand) % l;

            const uint64_t r = endpoints.rank(offset+1);
            const uint64_t next_endpoint = endpoints.select(r);

            if(offset + 64 >= next_endpoint)
                continue;

            const uint64_t kmer = access(0, offset);

            if ((i & 1) == 0)
                kmers.push_back(crc(kmer, k));
            else
                kmers.push_back(kmer);

            i++;
        }
        return kmers;
    }
};

const inline uint64_t RSHash::get_word64(uint64_t pos) {
    uint64_t block = pos >> 5;
    uint64_t shift = (pos & 31) << 1;
    uint64_t lo = text[block];
    uint64_t hi = text[block + 1];

    uint64_t shift_mask = -(shift != 0);
    return (lo >> shift) | ((hi << (64 - shift)) & shift_mask);
}

const inline uint64_t RSHash::get_base(uint64_t pos) {
    return (text[pos >> 5] >> ((pos & 31) << 1)) & 3ULL;
}

const inline uint64_t RSHash::access(const uint64_t unitig_id, const size_t offset) {
    return get_word64(offset) & kmermask;
}


template<int level>
inline void RSHash::fill_buffer(uint64_t *offsets, uint64_t *buffer, size_t p, size_t N, const uint64_t shift)
{
    uint64_t span;
    if constexpr (level == 1)
        span = span1;
    if constexpr (level == 2)
        span = span2;
    if constexpr (level == 3)
        span = span3;
    
    for(size_t i = 0; i < N; i++) {
        if constexpr (level == 1)
            offsets[i] = offsets1.access(p+i) + 1 - span;
        if constexpr (level == 2)
            offsets[i] = offsets2.access(p+i) + 1 - span;
        if constexpr (level == 3)
            offsets[i] = offsets3.access(p+i) + 1 - span;
    }

    for(uint64_t i = 0; i < N; i++) {
        const uint64_t s = offsets[i];

        uint64_t kmer = get_word64(s) & kmermask;
        uint64_t bits = get_word64(s + k);
        *buffer++ = kmer;
        for(uint64_t j=0; j < span-1; j++) {
            uint64_t const next_base = bits & 3ULL;
            bits >>= 2;
            kmer = (kmer >> 2) | (next_base << shift);
            *buffer++ = kmer;
        }

    }
}

template<int level>
inline bool RSHash::check_overlap(uint64_t skmer_pos, uint64_t text_pos, uint64_t &start_pos, uint64_t &end_pos)
{
    size_t span;
    if constexpr (level == 1)
        span = span1;
    if constexpr (level == 2)
        span = span2;
    if constexpr (level == 3)
        span = span3;
    
    const uint64_t r = endpoints.rank(skmer_pos+span);
    start_pos = endpoints.select(r-1, &end_pos);

    return text_pos >= start_pos && text_pos+k-1 < end_pos;
}


template<int level>
inline uint64_t RSHash::find_minimiser(const uint64_t kmer, const uint64_t kmer_rc, size_t &left_minimiser_position, size_t &right_minimiser_position)
{
    uint64_t m, mmermask;
    mixer_64 m_hasher;
    if constexpr (level == 1) {
        m = m1;
        mmermask = mmermask1;
        m_hasher = m_hasher1;
    }
    if constexpr (level == 2) {
        m = m2;
        mmermask = mmermask2;
        m_hasher = m_hasher2;
    }
    if constexpr (level == 3) {
        m = m3;
        mmermask = mmermask3;
        m_hasher = m_hasher3;
    }

    uint64_t mmer = kmer >> 2*(k - m);
    uint64_t mmer_rc = kmer_rc & mmermask;
    uint64_t minimiser = std::min<uint64_t>(m_hasher.hash(mmer) & mmermask, m_hasher.hash(mmer_rc) & mmermask);
    left_minimiser_position = k-m;
    right_minimiser_position = 0;

    for (size_t i = 1; i <= k-m; ++i) {
        mmer = (kmer >> 2*(k-m-i)) & mmermask;
        mmer_rc = (kmer_rc >> 2*i) & mmermask;
        const uint64_t mmerhash = std::min<uint64_t>(m_hasher.hash(mmer) & mmermask, m_hasher.hash(mmer_rc) & mmermask);
        if(mmerhash < minimiser) {
            minimiser = mmerhash;
            left_minimiser_position = k-m-i;
            right_minimiser_position = i;
        }
        else if(mmerhash == minimiser)
            left_minimiser_position = k-m-i;
    }

    return minimiser;
}

template<int level>
inline void RSHash::update_minimiser(const uint64_t kmer, const uint64_t kmer_rc, uint64_t &minimiser, size_t &left_minimiser_position, size_t &right_minimiser_position)
{
    uint64_t m, mmermask;
    mixer_64 m_hasher;
    if constexpr (level == 1) {
        m = m1;
        mmermask = mmermask1;
        m_hasher = m_hasher1;
    }
    if constexpr (level == 2) {
        m = m2;
        mmermask = mmermask2;
        m_hasher = m_hasher2;
    }
    if constexpr (level == 3) {
        m = m3;
        mmermask = mmermask3;
        m_hasher = m_hasher3;
    }

    if(left_minimiser_position-- == 0) {
        minimiser = find_minimiser<level>(kmer, kmer_rc, left_minimiser_position, right_minimiser_position);
        return;
    }

    const uint64_t mmer = kmer >> 2*(k - m);
    const uint64_t mmer_rc = kmer_rc & mmermask;
    const uint64_t mmerhash = std::min<uint64_t>(m_hasher.hash(mmer) & mmermask, m_hasher.hash(mmer_rc) & mmermask);

    if(mmerhash < minimiser) {
        minimiser = mmerhash;
        left_minimiser_position = k - m;
        right_minimiser_position = 0;
        return;
    }
    if(mmerhash == minimiser) {
        right_minimiser_position = 0;
        return;
    }

    right_minimiser_position++;
}

#endif