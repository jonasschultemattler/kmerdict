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

} // namespace cereal



int RSHash::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    cereal::BinaryOutputArchive archive(out);

    archive(k, level, m1, m2, m3, m_thres1, m_thres2, m_thres3, s1, s2, s3, endpoints, r1, r2, r3, offsets1, offsets2, offsets3, text, hashmap);

    out.close();
    return 0;
}

int RSHash::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    cereal::BinaryInputArchive archive(in);

    archive(k, level, m1, m2, m3, m_thres1, m_thres2, m_thres3, s1, s2, s3, endpoints, r1, r2, r3, offsets1, offsets2, offsets3, text, hashmap);

    std::cout << "loaded index...\n";
    in.close();

    this->span1 = k - m1 + 1;
    this->span2 = k - m2 + 1;
    this->span3 = k - m3 + 1;
    this->kmermask = compute_mask(2u * k);
    this->mmermask1 = compute_mask(2u * m1);
    this->mmermask2 = compute_mask(2u * m2);
    this->mmermask3 = compute_mask(2u * m3);

    this->s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), s1.size(), 3);
    this->s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), s2.size(), 3);
    this->s3_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s3.data()), s3.size(), 3);

    std::cout << "built rank and select ds...\n";
    
    return 0;
}