// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <deque>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <immintrin.h>


namespace rshash
{

struct shape_minimiser_and_positions_parameters
{
    size_t minimiser_size{};
    size_t window_size{};
    uint64_t seed{};
    uint32_t shape{};
};

struct shape_minimiser_and_positions_result
{
    uint64_t minimiser_value;
    size_t range_position;
};

}

namespace rshash::detail
{

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
class shape_minimiser_and_positions : public std::ranges::view_interface<shape_minimiser_and_positions<range_t>>
{
private:
    range_t range{};
    shape_minimiser_and_positions_parameters params{};

    template <bool range_is_const>
    class basic_iterator;

public:
    shape_minimiser_and_positions()
        requires std::default_initializable<range_t>
    = default;
    shape_minimiser_and_positions(shape_minimiser_and_positions const & rhs) = default;
    shape_minimiser_and_positions(shape_minimiser_and_positions && rhs) = default;
    shape_minimiser_and_positions & operator=(shape_minimiser_and_positions const & rhs) = default;
    shape_minimiser_and_positions & operator=(shape_minimiser_and_positions && rhs) = default;
    ~shape_minimiser_and_positions() = default;

    explicit shape_minimiser_and_positions(range_t range, shape_minimiser_and_positions_parameters params) :
        range{std::move(range)},
        params{std::move(params)}
    {}

    basic_iterator<false> begin()
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    basic_iterator<true> begin() const
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    auto end() noexcept
    {
        return std::default_sentinel;
    }

    auto end() const noexcept
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return std::default_sentinel;
    }
};

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
template <bool range_is_const>
class shape_minimiser_and_positions<range_t>::basic_iterator
{
private:
    template <bool>
    friend class basic_iterator;

    using maybe_const_range_t = std::conditional_t<range_is_const, range_t const, range_t>;
    using range_iterator_t = std::ranges::iterator_t<maybe_const_range_t>;

public:
    using difference_type = std::ranges::range_difference_t<maybe_const_range_t>;
    using value_type = shape_minimiser_and_positions_result;
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::conditional_t<std::ranges::forward_range<maybe_const_range_t>,
                                                 std::forward_iterator_tag,
                                                 std::input_iterator_tag>;
    using iterator_concept = iterator_category;

private:
    range_iterator_t range_it{};

    uint64_t kmer_mask{std::numeric_limits<uint64_t>::max()};
    uint64_t seed{};
    mixer_64 m_hasher;
    size_t minimisers_in_window{};
    uint64_t minimiser_size{};
    uint64_t window_size{};

    uint64_t kmer_value{};
    uint64_t kmer_value_rev{};

    size_t range_size{};
    size_t range_position{};

    value_type current{};
    value_type cached{};

    std::deque<uint64_t> kmers_in_window{};
    std::deque<size_t> minimiser_in_window{};

    static inline constexpr uint64_t compute_mask(uint64_t const size)
    {
        assert(size > 0u);
        assert(size <= 64u);

        if(size == 64u)
            return std::numeric_limits<uint64_t>::max();
        else
            return (uint64_t{1u} << (size)) - 1u;
    }

public:
    basic_iterator() = default;
    basic_iterator(basic_iterator const &) = default;
    basic_iterator(basic_iterator &&) = default;
    basic_iterator & operator=(basic_iterator const &) = default;
    basic_iterator & operator=(basic_iterator &&) = default;
    ~basic_iterator() = default;

    basic_iterator(basic_iterator<!range_is_const> const & it)
        requires range_is_const
        :
        range_it{it.range_it},
        kmer_mask{it.kmer_mask},
        kmer_value{it.kmer_value},
        kmer_value_rev{it.kmer_value_rev},
        range_size{it.range_size},
        range_position{it.range_position},
        current{it.current},
        cached{it.cached},
        kmers_in_window{it.kmers_in_window},
        minimisers_in_window{it.minimisers_in_window}
    {}

    basic_iterator(range_iterator_t range_iterator,
                   size_t const range_size,
                   shape_minimiser_and_positions_parameters const & params) :
        range_it{std::move(range_iterator)},
        kmer_mask{compute_mask(2u * params.minimiser_size)},
        range_size{range_size}
    {
        if (range_size < params.window_size)
            range_position = range_size;
        else
            init(params);
    }

    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return lhs.range_it == rhs.range_it;
    }

    friend bool operator==(basic_iterator const & lhs, std::default_sentinel_t const &)
    {
        return lhs.range_position > lhs.range_size;
    }

    basic_iterator & operator++() noexcept
    {
        while (!next_minimiser())
        {}
        return *this;
    }

    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        while (!next_minimiser())
        {}
        return tmp;
    }

    value_type operator*() const noexcept
    {
        return cached;
    }

private:
    enum class pop_first : bool
    {
        no,
        yes
    };

    void rolling_hash()
    {
        uint64_t const new_rank = seqan3::to_rank(*range_it);
        kmer_value = (kmer_value >> 2) | new_rank << 2*(minimiser_size-1);
        kmer_value_rev = ((kmer_value_rev << 2) | new_rank^(0b11)) & kmer_mask;
    }

    template <pop_first pop>
    void next_window()
    {
        ++range_position;
        ++range_it;

        rolling_hash();

        if constexpr (pop == pop_first::yes) {
            kmers_in_window.pop_front();
        }

        kmers_in_window.push_back(std::min<uint64_t>(m_hasher.hash(kmer_value) & kmer_mask, m_hasher.hash(kmer_value_rev) & kmer_mask));
    }

    void find_minimisers_in_window()
    {
        current.minimiser_value = kmers_in_window[0];
        current.range_position = 0;
        for(size_t i = 1; i <= minimisers_in_window; ++i) {
            if(kmers_in_window[i] < current.minimiser_value) {
                current.minimiser_value = kmers_in_window[i];
                current.range_position = i;
            }
        }
        for(size_t i = 0; i <= minimisers_in_window; ++i) {
            if(kmers_in_window[i] == current.minimiser_value)
                minimiser_in_window.push_back(i);
        }
        minimiser_in_window.pop_front();
    }

    void init(shape_minimiser_and_positions_parameters const & params)
    {
        seed = params.seed;
        m_hasher.seed(seed);
        minimiser_size = params.minimiser_size;
        window_size = params.window_size;
        minimisers_in_window = window_size - minimiser_size;

        uint64_t new_rank = seqan3::to_rank(*range_it);
        kmer_value |= new_rank << (2 * (minimiser_size - 1));
        kmer_value_rev |= new_rank^0b11;
        for (size_t i = 1u; i < params.minimiser_size; ++i) {
            ++range_position;
            ++range_it;
            new_rank = seqan3::to_rank(*range_it);
            kmer_value >>= 2;
            kmer_value |= new_rank << (2 * (minimiser_size - 1));
            kmer_value_rev <<= 2;
            kmer_value_rev |= new_rank^0b11;
        }

        kmers_in_window.push_back(std::min<uint64_t>(m_hasher.hash(kmer_value) & kmer_mask, m_hasher.hash(kmer_value_rev) & kmer_mask));

        for (size_t i = minimiser_size; i < window_size; ++i)
            next_window<pop_first::no>();

        find_minimisers_in_window();

        while (!next_minimiser()) {}
    }

    void update_cache() {
        cached.minimiser_value = current.minimiser_value;
        cached.range_position = range_position - window_size + current.range_position;
    }

    bool next_minimiser()
    {
        if(!minimiser_in_window.empty()) {
            cached.range_position = range_position + 1 - window_size + current.range_position;
            current.range_position = minimiser_in_window.front();
            minimiser_in_window.pop_front();
            // std::cout << cached.range_position << " " << range_position << " " << current.range_position;
            // std::cout << " -> " << current.range_position << "\n";
            return true;
        }
        
        if(range_position + 1 >= range_size) {
            ++range_position;
            update_cache();
            return true;
        }

        next_window<pop_first::yes>();

        if(current.range_position == 0) {
            update_cache();
            find_minimisers_in_window();
            return true;
        }

        if(uint64_t new_kmer_hash = kmers_in_window.back(); new_kmer_hash <= current.minimiser_value) {
            update_cache();
            current.minimiser_value = new_kmer_hash;
            current.range_position = minimisers_in_window;
            return true;
        }

        current.range_position--;
        return false;
    }

};


template <std::ranges::viewable_range rng_t>
shape_minimiser_and_positions(rng_t &&, shape_minimiser_and_positions_parameters &&)
    -> shape_minimiser_and_positions<std::views::all_t<rng_t>>;

struct shape_minimiser_and_positions_fn
{
    constexpr auto operator()(shape_minimiser_and_positions_parameters params) const
    {
        return seqan3::detail::adaptor_from_functor{*this, std::move(params)};
    }

    template <std::ranges::range range_t>
    constexpr auto operator()(range_t && range, shape_minimiser_and_positions_parameters params) const
    {
        static_assert(std::same_as<std::ranges::range_value_t<range_t>, seqan3::dna4>, "Only dna4 supported.");
        static_assert(std::ranges::sized_range<range_t>, "Input range must be a std::ranges::sized_range.");

        if (params.minimiser_size == 0u)
            throw std::invalid_argument{"minimiser_size must be > 0."};
        if (params.minimiser_size > 32u)
            throw std::invalid_argument{"minimiser_size must be <= 32."};
        if (params.window_size == 0u)
            throw std::invalid_argument{"window_size must be > 0."};
        if (params.window_size > 32u)
            throw std::invalid_argument{"window_size must be <= 32."};
        if (params.window_size < params.minimiser_size)
            throw std::invalid_argument{"window_size must be >= minimiser_size."};

        return shape_minimiser_and_positions{range, std::move(params)};
    }
};

}

namespace rshash::views
{

inline constexpr auto shape_minimiser_and_positions = rshash::detail::shape_minimiser_and_positions_fn{};

}


namespace rshash
{

struct shape_minimiser_and_occurencess_parameter
{
    size_t minimiser_size{};
    size_t window_size{};
    uint64_t seed{};
    uint32_t shape{};
};

struct shape_minimiser_and_occurencess_result
{
    uint64_t minimiser_value;
    size_t range_position;
    size_t occurrences;
};

}

namespace rshash::detail
{

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
class shape_minimiser_and_occurencess : public std::ranges::view_interface<shape_minimiser_and_occurencess<range_t>>
{
private:
    range_t range{};
    shape_minimiser_and_occurencess_parameter params{};

    template <bool range_is_const>
    class basic_iterator;

public:
    shape_minimiser_and_occurencess()
        requires std::default_initializable<range_t>
    = default;
    shape_minimiser_and_occurencess(shape_minimiser_and_occurencess const & rhs) = default;
    shape_minimiser_and_occurencess(shape_minimiser_and_occurencess && rhs) = default;
    shape_minimiser_and_occurencess & operator=(shape_minimiser_and_occurencess const & rhs) = default;
    shape_minimiser_and_occurencess & operator=(shape_minimiser_and_occurencess && rhs) = default;
    ~shape_minimiser_and_occurencess() = default;

    explicit shape_minimiser_and_occurencess(range_t range, shape_minimiser_and_occurencess_parameter params) :
        range{std::move(range)},
        params{std::move(params)}
    {}

    basic_iterator<false> begin()
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    basic_iterator<true> begin() const
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    auto end() noexcept
    {
        return std::default_sentinel;
    }

    auto end() const noexcept
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return std::default_sentinel;
    }
};

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
template <bool range_is_const>
class shape_minimiser_and_occurencess<range_t>::basic_iterator
{
private:
    template <bool>
    friend class basic_iterator;

    using maybe_const_range_t = std::conditional_t<range_is_const, range_t const, range_t>;
    using range_iterator_t = std::ranges::iterator_t<maybe_const_range_t>;

public:
    using difference_type = std::ranges::range_difference_t<maybe_const_range_t>;
    using value_type = shape_minimiser_and_occurencess_result;
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::conditional_t<std::ranges::forward_range<maybe_const_range_t>,
                                                 std::forward_iterator_tag,
                                                 std::input_iterator_tag>;
    using iterator_concept = iterator_category;

private:
    range_iterator_t range_it{};

    uint64_t kmer_mask{std::numeric_limits<uint64_t>::max()};
    uint64_t kmer_value{};
    uint64_t kmer_value_rev{};
    uint64_t seed{};
    mixer_64 m_hasher;

    uint64_t minimiser_size{};

    size_t range_size{};
    size_t range_position{};

    shape_minimiser_and_occurencess_parameter params{};
    value_type current{}; // range_position -> position in the window
    value_type cached{};  // range_position -> position in the range

    std::deque<uint64_t> kmer_values_in_window{};

    static inline constexpr uint64_t compute_mask(uint64_t const size) {
        assert(size > 0u);
        assert(size <= 64u);

        if (size == 64u)
            return std::numeric_limits<uint64_t>::max();
        else
            return (uint64_t{1u} << (size)) - 1u;
    }

public:
    basic_iterator() = default;
    basic_iterator(basic_iterator const &) = default;
    basic_iterator(basic_iterator &&) = default;
    basic_iterator & operator=(basic_iterator const &) = default;
    basic_iterator & operator=(basic_iterator &&) = default;
    ~basic_iterator() = default;

    basic_iterator(basic_iterator<!range_is_const> const & it)
        requires range_is_const
        :
        range_it{it.range_it},
        kmer_mask{it.kmer_mask},
        kmer_value{it.kmer_value},
        range_size{it.range_size},
        range_position{it.range_position},
        params{it.params},
        current{it.current},
        cached{it.cached},
        kmer_values_in_window{it.kmer_values_in_window}
    {}

    basic_iterator(range_iterator_t range_iterator,
                   size_t const range_size,
                   shape_minimiser_and_occurencess_parameter params) :
        range_it{std::move(range_iterator)},
        kmer_mask{compute_mask(2u * params.minimiser_size)},
        minimiser_size{params.minimiser_size},
        range_size{range_size},
        params{std::move(params)}
    {
        if (range_size < params.window_size)
            range_position = range_size + 1u;
        else
            init();
    }

    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return lhs.range_it == rhs.range_it;
    }

    friend bool operator==(basic_iterator const & lhs, std::default_sentinel_t const &)
    {
        return lhs.range_position > lhs.range_size;
    }

    basic_iterator & operator++() noexcept
    {
        while (!next_minimiser_is_new())
        {}
        return *this;
    }

    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        while (!next_minimiser_is_new())
        {}
        return tmp;
    }

    value_type operator*() const noexcept
    {
        return cached;
    }

private:
    enum class pop_first : bool
    {
        no,
        yes
    };

    void rolling_hash()
    {
        uint64_t const new_rank = seqan3::to_rank(*range_it);
        kmer_value = (kmer_value >> 2) | (new_rank << 2*(minimiser_size-1));
        kmer_value_rev = ((kmer_value_rev << 2) | new_rank^(0b11)) & kmer_mask;
    }

    template <pop_first pop>
    void next_window()
    {
        ++range_position;
        ++range_it;

        rolling_hash();

        if constexpr (pop == pop_first::yes) {
            kmer_values_in_window.pop_front();
        }

        // const uint64_t kmerhash = std::min<uint64_t>(m_hasher.hash(kmer_value), m_hasher.hash(kmer_value_rev)) & kmer_mask;
        const uint64_t kmerhash = std::min<uint64_t>(m_hasher.hash(kmer_value) & kmer_mask, m_hasher.hash(kmer_value_rev) & kmer_mask);

        kmer_values_in_window.push_back(kmerhash);
    }

    void find_minimiser_in_window()
    {
        auto minimiser_it = std::ranges::min_element(kmer_values_in_window, std::less_equal<uint64_t>{});
        current.range_position = std::distance(std::begin(kmer_values_in_window), minimiser_it);
        current.minimiser_value = *minimiser_it;
    }

    void init()
    {
        seed = params.seed;
        m_hasher.seed(seed);
        minimiser_size = params.minimiser_size;

        uint64_t new_rank = seqan3::to_rank(*range_it);
        kmer_value = new_rank << 2*(minimiser_size-1);
        kmer_value_rev = new_rank^(0b11);
        for (size_t i = 1u; i < minimiser_size; ++i) {
            ++range_position;
            ++range_it;

            new_rank = seqan3::to_rank(*range_it);
            kmer_value = (kmer_value >> 2) | (new_rank << 2*(minimiser_size-1));
            kmer_value_rev = ((kmer_value_rev << 2) | new_rank^(0b11)) & kmer_mask;
        }

        // const uint64_t kmerhash = std::min<uint64_t>(m_hasher.hash(kmer_value), m_hasher.hash(kmer_value_rev)) & kmer_mask;
        const uint64_t kmerhash = std::min<uint64_t>(m_hasher.hash(kmer_value) & kmer_mask, m_hasher.hash(kmer_value_rev) & kmer_mask);
        
        kmer_values_in_window.push_back(kmerhash);

        for (size_t i = minimiser_size; i < params.window_size; ++i)
            next_window<pop_first::no>();

        find_minimiser_in_window();

        while (!next_minimiser_is_new())
        {}
    }

    void update_cache()
    {
        cached.minimiser_value = current.minimiser_value;
        cached.range_position = range_position - params.window_size - current.occurrences; // skmer position
        cached.occurrences = current.occurrences + 1;
    }

    bool next_minimiser_is_new()
    {
        // If we reached the end of the range, we are done.
        if (range_position + 1 >= range_size)
        {
            ++range_position;
            update_cache();
            return true;
        }

        next_window<pop_first::yes>();

        // The minimiser left the window.
        if (current.range_position == 0)
        {
            update_cache();
            find_minimiser_in_window();
            current.occurrences = 0;
            return true;
        }

        if (uint64_t new_kmer_hash = kmer_values_in_window.back(); new_kmer_hash < current.minimiser_value)
        {
            update_cache();
            current.minimiser_value = new_kmer_hash;
            current.range_position = kmer_values_in_window.size() - 1u;
            current.occurrences = 0;
            return true;
        }

        --current.range_position;
        ++current.occurrences;
        return false;
    }
};

template <std::ranges::viewable_range rng_t>
shape_minimiser_and_occurencess(rng_t &&, shape_minimiser_and_occurencess_parameter &&)
    -> shape_minimiser_and_occurencess<std::views::all_t<rng_t>>;

struct shape_minimiser_and_occurencess_fn
{
    constexpr auto operator()(shape_minimiser_and_occurencess_parameter params) const
    {
        return seqan3::detail::adaptor_from_functor{*this, std::move(params)};
    }

    template <std::ranges::range range_t>
    constexpr auto operator()(range_t && range, shape_minimiser_and_occurencess_parameter params) const
    {
        static_assert(std::same_as<std::ranges::range_value_t<range_t>, seqan3::dna4>, "Only dna4 supported.");
        static_assert(std::ranges::sized_range<range_t>, "Input range must be a std::ranges::sized_range.");

        if (params.minimiser_size == 0u)
            throw std::invalid_argument{"minimiser_size must be > 0."};
        if (params.minimiser_size > 32u)
            throw std::invalid_argument{"minimiser_size must be <= 32."};
        if (params.window_size < params.minimiser_size)
            throw std::invalid_argument{"window_size must be >= minimiser_size."};

        return shape_minimiser_and_occurencess{std::forward<range_t>(range), std::move(params)};
    }
};

}

namespace rshash::views
{

inline constexpr auto shape_minimiser_and_occurencess = rshash::detail::shape_minimiser_and_occurencess_fn{};

}



namespace rshash
{

struct shapeview_parameters
{
    uint32_t shape{};
};

struct shapeview_result
{
    uint64_t kmer_value;
    uint64_t kmer_value_rev;
};

}

namespace rshash::detail
{

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
class shapeview : public std::ranges::view_interface<shapeview<range_t>>
{
private:
    range_t range{};
    shapeview_parameters params{};

    template <bool range_is_const>
    class basic_iterator;

public:
    shapeview()
        requires std::default_initializable<range_t>
    = default;
    shapeview(shapeview const & rhs) = default;
    shapeview(shapeview && rhs) = default;
    shapeview & operator=(shapeview const & rhs) = default;
    shapeview & operator=(shapeview && rhs) = default;
    ~shapeview() = default;

    explicit shapeview(range_t range, shapeview_parameters params) :
        range{std::move(range)},
        params{std::move(params)}
    {}

    basic_iterator<false> begin()
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    basic_iterator<true> begin() const
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return {std::ranges::begin(range), std::ranges::size(range), params};
    }

    auto end() noexcept
    {
        return std::default_sentinel;
    }

    auto end() const noexcept
        requires std::ranges::view<range_t const> && std::ranges::input_range<range_t const>
    {
        return std::default_sentinel;
    }
};

template <std::ranges::view range_t>
    requires std::ranges::input_range<range_t> && std::ranges::sized_range<range_t>
template <bool range_is_const>
class shapeview<range_t>::basic_iterator
{
private:
    template <bool>
    friend class basic_iterator;

    using maybe_const_range_t = std::conditional_t<range_is_const, range_t const, range_t>;
    using range_iterator_t = std::ranges::iterator_t<maybe_const_range_t>;

public:
    using difference_type = std::ranges::range_difference_t<maybe_const_range_t>;
    using value_type = shapeview_result;
    using pointer = void;
    using reference = value_type;
    using iterator_category = std::conditional_t<std::ranges::forward_range<maybe_const_range_t>,
                                                 std::forward_iterator_tag,
                                                 std::input_iterator_tag>;
    using iterator_concept = iterator_category;

private:
    range_iterator_t range_it{};

    uint64_t kmer_mask{std::numeric_limits<uint64_t>::max()};
    uint64_t shape_length{};
    uint64_t shape_mask{};
    uint64_t shape_mask_rc{};
    uint64_t kmer_value{};
    uint64_t kmer_value_rev{};

    size_t range_size{};
    size_t range_position{};

    value_type current{};

    static inline constexpr uint64_t compute_mask(uint64_t const size) {
        assert(size > 0u);
        assert(size <= 64u);

        if(size == 64u)
            return std::numeric_limits<uint64_t>::max();
        else
            return (uint64_t{1u} << (size)) - 1u;
    }

    static inline constexpr uint64_t compute_shape_mask(uint32_t const shape) {
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

public:
    basic_iterator() = default;
    basic_iterator(basic_iterator const &) = default;
    basic_iterator(basic_iterator &&) = default;
    basic_iterator & operator=(basic_iterator const &) = default;
    basic_iterator & operator=(basic_iterator &&) = default;
    ~basic_iterator() = default;

    basic_iterator(basic_iterator<!range_is_const> const & it)
        requires range_is_const
        :
        range_it{it.range_it},
        kmer_mask{it.kmer_mask},
        range_size{it.range_size},
        range_position{it.range_position},
        current{it.current}
    {}

    basic_iterator(range_iterator_t range_iterator,
                   size_t const range_size,
                   shapeview_parameters const & params) :
        range_it{std::move(range_iterator)},
        shape_mask{compute_shape_mask(params.shape)},
        shape_mask_rc{reverse_shape(shape_mask)},
        shape_length{bit_length(params.shape)},
        kmer_mask{compute_mask(2u * shape_length)},
        range_size{range_size}
    {
        if (range_size < shape_length)
            range_position = range_size;
        else
            init(params);
    }

    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return lhs.range_it == rhs.range_it;
    }

    friend bool operator==(basic_iterator const & lhs, std::default_sentinel_t const &)
    {
        return lhs.range_position == lhs.range_size;
    }

    basic_iterator & operator++() noexcept
    {
        ++range_position;
        ++range_it;
        rolling_hash();
        return *this;
    }

    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        ++range_position;
        ++range_it;
        rolling_hash();
        return tmp;
    }

    value_type operator*() const noexcept
    {
        return current;
    }

private:
    enum class pop_first : bool
    {
        no,
        yes
    };

    void rolling_hash()
    {
        uint64_t const new_rank = seqan3::to_rank(*range_it);
        kmer_value = (kmer_value >> 2) | (new_rank << 2*(shape_length-1));
        kmer_value_rev = ((kmer_value_rev << 2) | (new_rank^0b11)) & kmer_mask;
        // current.shape_value_fwd = _pext_u64(kmer_value, shape_mask);
        // current.shape_value_rc = _pext_u64(kmer_value_rev, shape_mask_rc);
        current.kmer_value = _pext_u64(kmer_value, shape_mask);
        current.kmer_value_rev = _pext_u64(kmer_value_rev, shape_mask_rc);
    }

    void init(shapeview_parameters const & params)
    {
        shape_mask = compute_shape_mask(params.shape);
        shape_mask_rc = reverse_shape(shape_mask);
        shape_length = bit_length(params.shape);

        uint64_t new_rank = seqan3::to_rank(*range_it);
        kmer_value = new_rank << 2*(shape_length-1);
        kmer_value_rev = new_rank^0b11;
        for (size_t i = 1u; i < shape_length; ++i) {
            ++range_position;
            ++range_it;
            rolling_hash();
        }

        // current.shape_value_fwd = _pext_u64(kmer_value, shape_mask);
        // current.shape_value_rc = _pext_u64(kmer_value_rev, shape_mask_rc);
        current.kmer_value = _pext_u64(kmer_value, shape_mask);
        current.kmer_value_rev = _pext_u64(kmer_value_rev, shape_mask_rc);
    }

};


template <std::ranges::viewable_range rng_t>
shapeview(rng_t &&, shapeview_parameters &&)
    -> shapeview<std::views::all_t<rng_t>>;

struct shapeview_fn
{
    constexpr auto operator()(shapeview_parameters params) const
    {
        return seqan3::detail::adaptor_from_functor{*this, std::move(params)};
    }

    template <std::ranges::range range_t>
    constexpr auto operator()(range_t && range, shapeview_parameters params) const
    {
        static_assert(std::same_as<std::ranges::range_value_t<range_t>, seqan3::dna4>, "Only dna4 supported.");
        static_assert(std::ranges::sized_range<range_t>, "Input range must be a std::ranges::sized_range.");

        if (params.shape == 0u)
            throw std::invalid_argument{"shape must be > 0."};

        return shapeview{range, std::move(params)};
    }
};

}

namespace rshash::views
{

inline constexpr auto shapeview = rshash::detail::shapeview_fn{};

}

