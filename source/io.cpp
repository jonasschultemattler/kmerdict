#include <bitset>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include "rshash.hpp"

namespace cereal {

template <class Archive>
void save(Archive& ar, const gtl::flat_hash_set<uint64_t>& set) {
    ar(static_cast<std::size_t>(set.size()));
    for (const auto& elem : set) {
        ar(elem);
    }
}

template <class Archive>
void load(Archive& ar, gtl::flat_hash_set<uint64_t>& set) {
    std::size_t size;
    ar(size);
    set.clear();
    set.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        uint64_t elem;
        ar(elem);
        set.insert(std::move(elem));
    }
}


template <class Archive>
void save(Archive& ar, const gtl::flat_hash_map<uint64_t, uint64_t>& map) {
    ar(static_cast<std::size_t>(map.size()));
    for (const auto& [key, val] : map) {
        ar(key, val);
    }
}

template <class Archive>
void load(Archive& ar, gtl::flat_hash_map<uint64_t, uint64_t>& map) {
    std::size_t size;
    ar(size);
    map.clear();
    map.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        uint64_t key;
        uint64_t val;
        ar(key, val);
        map.emplace(key, val);
    }
}

template <class Archive>
void save(Archive& ar, const gtl::flat_hash_map<uint64_t, std::vector<uint32_t>>& map) {
    ar(static_cast<std::size_t>(map.size()));
    for (const auto& [key, val] : map) {
        ar(key, val);
    }
}

template <class Archive>
void load(Archive& ar, gtl::flat_hash_map<uint64_t, std::vector<uint32_t>>& map) {
    std::size_t size;
    ar(size);
    map.clear();
    map.reserve(size);
    for (std::size_t i = 0; i < size; ++i) {
        uint64_t key;
        std::vector<uint32_t> val;
        ar(key, val);
        map.emplace(key, val);
    }
}


} // namespace cereal


template <class Archive>
void save(Archive& ar, const FlatMap& map) {
    ar(static_cast<std::size_t>(map.index.size()));

    for (const auto& [key, entry] : map.index)
        ar(key, entry.offset, entry.size);

    ar(map.values);
}

template <class Archive>
void load(Archive& ar, FlatMap& map) {
    std::size_t size;
    ar(size);

    map.index.clear();
    map.index.reserve(size);

    for (std::size_t i = 0; i < size; ++i) {
        uint64_t key;
        uint32_t offset;
        uint32_t value_count;

        ar(key, offset, value_count);

        map.index.emplace(key, FlatMap::Entry{.offset = offset, .size = value_count});
    }

    ar(map.values);
}


int RSHash::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    cereal::BinaryOutputArchive archive(out);

    archive(k, shape.value, windowmask, windowshift, window_size, shape_overlap_left, shape_overlap_right, overlap,
        shape_mask, shape_mask_rev, kernel_mask, kernel_mask_rev, mmermask1, mmermask2, mmermask3,
        level, m1, m2, m3, m_thres1, m_thres2, m_thres3, threshold, span1, span2, span3,
        s1, s2, s3, s4, s5, endpoints, r1, r2, r3, r4, r5, offsets1, offsets2, offsets3, offsets4, offsets5, text, loc, use_ht, hashset, hashset_rc, hashmap, hashmap_rc);

    out.close();
    return 0;
}

int RSHash::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    cereal::BinaryInputArchive archive(in);

    uint32_t shape_value;
    archive(k, shape_value, windowmask, windowshift, window_size, shape_overlap_left, shape_overlap_right, overlap,
        shape_mask, shape_mask_rev, kernel_mask, kernel_mask_rev, mmermask1, mmermask2, mmermask3,
        level, m1, m2, m3, m_thres1, m_thres2, m_thres3, threshold, span1, span2, span3,
        s1, s2, s3, s4, s5, endpoints, r1, r2, r3, r4, r5, offsets1, offsets2, offsets3, offsets4, offsets5, text, loc, use_ht, hashset, hashset_rc, hashmap, hashmap_rc);
    in.close();

    this->shape = shape32_create(shape_value);

    // std::cout << shape_value << " " << std::bitset<64>(windowmask) << " " << window_size << " ";
    // std::cout << std::bitset<64>(shape_mask) << " " << std::bitset<64>(shape_mask_rev) << " " << shape_overlap_left << " " << shape_overlap_right << " " << std::bitset<64>(kernel_mask) << " " << std::bitset<64>(kernel_mask_rev) << '\n';

    this->s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), s1.size(), 3);
    this->s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), s2.size(), 3);
    this->s3_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s3.data()), s3.size(), 3);
    this->s4_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s4.data()), s4.size(), 3);
    this->s5_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s5.data()), s5.size(), 3);

    initialise_lookupfn();
    if(loc)
        initialise_locatefn();
    
    std::cout << "loaded index...\n";

    return 0;
}