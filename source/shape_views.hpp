// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <deque>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>

#include "shape.hpp"



namespace rshash
{

struct shapeview_parameters
{
    uint32_t shape;
};

struct shapeview_result
{
    uint64_t value;
    uint64_t value_rev;
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
    uint64_t window_size{};
    Shape32 shape;

    size_t range_size{};
    size_t range_position{};

    uint64_t kmer_value;
    uint64_t kmer_value_rc;

    value_type current{};

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
        range_size{it.range_size},
        range_position{it.range_position},
        current{it.current}
    {}

    basic_iterator(range_iterator_t range_iterator,
                   size_t const range_size,
                   shapeview_parameters const & params) :
        range_it{std::move(range_iterator)},
        kmer_mask{compute_mask(2u * window_size)},
        range_size{range_size}
    {
        if (range_size < window_size)
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
        kmer_value = (kmer_value >> 2) | (new_rank << 2*(window_size-1));
        kmer_value_rc = ((kmer_value_rc << 2) | (new_rank^0b11)) & kmer_mask; // no & ?
        current.value = _pext_u64(kmer_value, shape.mask);
        current.value_rev = _pext_u64(kmer_value_rc, shape.mask_rev);
    }

    void init(shapeview_parameters const & params)
    {
        shape = shape32_create(params.shape);
        window_size = shape.length;
        kmer_mask = compute_mask(2u * window_size);

        uint64_t new_rank = seqan3::to_rank(*range_it);
        kmer_value = new_rank << 2*(window_size-1);
        kmer_value_rc = new_rank^0b11;
        for (size_t i = 1u; i < window_size; ++i) {
            ++range_position;
            ++range_it;
            rolling_hash();
        }
        
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

        // if (params.window_size == 0u)
        //     throw std::invalid_argument{"window_size must be > 0."};
        // if (params.window_size > 32u)
        //     throw std::invalid_argument{"window_size must be <= 32."};

        return shapeview{range, std::move(params)};
    }
};

}

namespace rshash::views
{

inline constexpr auto shapeview = rshash::detail::shapeview_fn{};

}

