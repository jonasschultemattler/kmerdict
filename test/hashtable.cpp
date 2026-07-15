#include <sharg/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include "../source/minimiser_views.hpp"
#include "../source/shape.hpp"
#include <gtl/phmap.hpp>


struct cmd_arguments {
    std::string cmd{};
    std::filesystem::path i{};
    std::filesystem::path q{};
    std::filesystem::path o{};
    uint8_t k{31};
    uint32_t shape{std::numeric_limits<uint32_t>::max()};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.q, sharg::config{.short_id = 'q', .long_id = "query", .description = "provide query file"});
    parser.add_option(args.k, sharg::config{.short_id = 'k', .long_id = "k-mer", .description = "k-mer length"});
    parser.add_option(args.shape, sharg::config{.long_id = "shape", .description = "shape value"});
}

int check_arguments(sharg::parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('i'))
        throw sharg::user_input_error("provide input file.");
    if(!parser.is_option_set('q'))
        throw sharg::user_input_error("provide query file.");

    return 0;
}

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


void load_file(const std::filesystem::path &filepath,
               std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &output)
{
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        seqan3::bitpacked_sequence<seqan3::dna4> seq;
        seq.assign(record.sequence().begin(), record.sequence().end());
        output.push_back(std::move(seq));
    }
}

int main(int argc, char** argv)
{
    sharg::parser parser{"HT", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
        check_arguments(parser, args);
    }
    catch (sharg::parser_error const &ext) {
        return -1;
    }

    std::cout << "loading text...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> text;
    load_file(args.i, text);

    std::cout << "building hashtable...\n";
    std::unordered_set<uint64_t> ht_fwd;
    std::unordered_set<uint64_t> ht_rc;
    // gtl::flat_hash_map<uint64_t, uint64_t> ht;

    Shape32 shape_obj;
    size_t window_size;
    if(args.shape != std::numeric_limits<uint32_t>::max()) {
        shape_obj = shape32_create(args.shape);
        window_size = shape_obj.length;
        // std::cout << shape_value << " " << std::bitset<64>(windowmask) << " " << window_size << " " << overlap << " "
        //         << std::bitset<64>(shape_mask) << " " << std::bitset<64>(shape_mask_rev) << " " << shift_shape << " " << shift_shape_rev << " "
        //         << " " << std::bitset<64>(kmermask) << '\n';
    }
    else
        window_size = args.k;

    for(auto & sequence : text) {
        for(auto && window : sequence | rshash::views::kmerview({.window_size = window_size})) {
            if(args.shape != std::numeric_limits<uint32_t>::max()) {
                uint64_t kmer = _pext_u64(window.value, shape_obj.mask);
                uint64_t kmer_rc = _pext_u64(window.value, shape_obj.mask_rev);
                // ht[kmer]++;
                // ht[kmer_rc]++;
                ht_fwd.insert(kmer);
                ht_rc.insert(kmer_rc);
            }
            else {
                // ht[std::min<uint64_t>(window.value, window.value_rev)]++;
                ht_fwd.insert(std::min<uint64_t>(window.value, window.value_rev));
            }
        }   
    }
     
    std::cout << "loading queries...\n";
    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> queries;
    load_file(args.q, queries);

    std::cout << "querying...\n";
    uint64_t kmers = 0;
    uint64_t found_kmers = 0;
    uint64_t found_positions = 0;

    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    if (args.shape != std::numeric_limits<uint32_t>::max()) {
        for (auto& query : queries) {
            for (auto&& window : query | rshash::views::kmerview({.window_size = window_size})) {
                uint64_t kmer_fwd = _pext_u64(window.value, shape_obj.mask);
                uint64_t kmer_rev = _pext_u64(window.value_rev, shape_obj.mask_rev);
                found_kmers += ht_fwd.contains(kmer_fwd) || ht_rc.contains(kmer_rev);
                // if (auto it = ht.find(kmer_fwd); it != ht.end()) {
                //     found_kmers++;
                //     // found_positions += it->second;
                // }
                // else if (auto it = ht.find(kmer_rev); it != ht.end()) {
                //     found_kmers++;
                //     // found_positions += it->second;
                // }
                // else if (auto it = ht.find(_pext_u64(window.value, shape_obj.mask_rev)); it != ht.end()) {
                //     found_kmers++;
                //     // found_positions += it->second;
                // }
                // else if (auto it = ht.find(_pext_u64(window.value_rev, shape_obj.mask)); it != ht.end()) {
                //     found_kmers++;
                //     // found_positions += it->second;
                // }
                kmers++;
            }
        }
    }
    else {
        for (auto& query : queries) {
            for (auto&& window : query | rshash::views::kmerview({.window_size = args.k})) {
                uint64_t kmer_value = std::min<uint64_t>(window.value, window.value_rev);
                // if (auto it = ht.find(kmer_value); it != ht.end()) {
                //     found_kmers++;
                //     found_positions += it->second;
                // }
                found_kmers += ht_fwd.contains(kmer_value);
                kmers++;
            }
        }
    }
    
    std::chrono::high_resolution_clock::time_point t_stop = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_stop - t_start);
        
    double ns_per_kmer = (double) elapsed.count() / kmers;
        
    std::cout << "==== query report:\n";
    std::cout << "num_kmers = " << kmers << '\n';
    std::cout << "num_positive_kmers = " << found_kmers << " (" << (double) found_kmers/kmers*100 << "%)\n";
    std::cout << "num_positions = " << found_positions << '\n';
    std::cout << "time_per_kmer = " << ns_per_kmer << '\n';
 
    return 0;
}