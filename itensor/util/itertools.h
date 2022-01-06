/*******************************************************************************
 * 
 * Copyright (C) 2015-2018, O. Parcollet
 * Copyright (C) 2019, The Simons Foundation
 *   author : N. Wentzell
 * 
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 ******************************************************************************/

#pragma once
#include <tuple>
#include <vector>
#include <iterator>
#include <exception>

namespace triqs::utility {

  template <class Iter, class Value, class Tag = std::forward_iterator_tag, class Reference = Value &, class Difference = std::ptrdiff_t>
  struct iterator_facade;

  /**
   * A helper for the implementation of forward iterators using CRTP
   *
   * @tparam Iter
   * The Iterator Class to be implemented
   * `Iter` is required to have the following member functions
   * - bool equal(Iter other)
   * - value_type [const] [&] dereference()
   * - void increment()
   */
  template <typename Iter, typename Value, typename Reference, typename Difference>
  struct iterator_facade<Iter, Value, std::forward_iterator_tag, Reference, Difference> {

    private:
    Iter &self() { return static_cast<Iter &>(*this); }
    Iter const &self() const { return static_cast<const Iter &>(*this); }

    public:
    using value_type        = Value;
    using reference         = Reference;
    using pointer           = Value *;
    using difference_type   = Difference;
    using iterator_category = std::forward_iterator_tag;

    Iter &operator++() {
      self().increment();
      return self();
    }

    Iter operator++(int) {
      Iter c = self();
      self().increment();
      return c;
    }

    template <typename U> bool operator==(U const &other) const { return self().equal(other); }
    template <typename U> bool operator!=(U const &other) const { return (!operator==(other)); }

    decltype(auto) operator*() const { return self().dereference(); }
    decltype(auto) operator-> () const { return operator*(); }
  };

  /********************* Enumerate Iterator ********************/

  template <typename Iter> struct enum_iter : iterator_facade<enum_iter<Iter>, std::pair<long, typename std::iterator_traits<Iter>::value_type>> {

    Iter it;
    long i = 0;

    enum_iter(Iter it) : it(std::move(it)) {}

    void increment() {
      ++it;
      ++i;
    }

    bool equal(enum_iter const &other) const { return (it == other.it); }

    decltype(auto) dereference() const { return std::tuple<long, decltype(*it)>{i, *it}; }
  };

  /********************* Transform Iterator ********************/

  template <typename Iter, typename L, typename Value = std::invoke_result_t<L,typename std::iterator_traits<Iter>::value_type>>
  struct transform_iter : iterator_facade<transform_iter<Iter, L>, Value> {

    Iter it;
    mutable L lambda;

    transform_iter(Iter it, L lambda) : it(std::move(it)), lambda(std::move(lambda)) {}

    void increment() { ++it; }

    bool equal(transform_iter const &other) const { return (it == other.it); }

    decltype(auto) dereference() const { return lambda(*it); }
  };

  /********************* Zip Iterator ********************/

  template <typename... It> struct zip_iter : iterator_facade<zip_iter<It...>, std::tuple<typename std::iterator_traits<It>::value_type...>> {

    std::tuple<It...> its;
    zip_iter(std::tuple<It...> its) : its(std::move(its)) {}

    template <typename... U> void nop(U &&... u) {} // do nothing...
    template <size_t... Is> void increment_all(std::index_sequence<Is...>) { nop(++std::get<Is>(its)...); }
    void increment() { increment_all(std::index_sequence_for<It...>{}); }

    bool equal(zip_iter const &other) const { return (its == other.its); }

    template <size_t... Is> auto tuple_map_impl(std::index_sequence<Is...>) const {
      return std::tuple<decltype(*std::get<Is>(its))...>(*std::get<Is>(its)...);
    }

    decltype(auto) dereference() const { return tuple_map_impl(std::index_sequence_for<It...>{}); }
  };

  /********************* Product Iterator ********************/

  namespace details {
    // Sentinel_t, used to denote the end of a product range
    template <typename It> struct sentinel_t { It it; };
    template <typename It> sentinel_t<It> make_sentinel(It &&it) { return {it}; }
  } // namespace details

  template <typename... It> struct prod_iter : iterator_facade<prod_iter<It...>, std::tuple<typename std::iterator_traits<It>::value_type...>> {

    std::tuple<It...> its_begin, its_end;
    std::tuple<It...> its = its_begin;

    prod_iter(std::tuple<It...> its_begin, std::tuple<It...> its_end)
       : its_begin(std::move(its_begin)), its_end(std::move(its_end)){}

    template <int N> void _increment() {
      ++std::get<N>(its);
      if constexpr (N < sizeof...(It) - 1) {
        if (std::get<N>(its) == std::get<N>(its_end)) {
          std::get<N>(its) = std::get<N>(its_begin);
          _increment<N + 1>();
        }
      }
    }
    void increment() { _increment<0>(); }

    bool equal(prod_iter const &other) const { return (its == other.its); }
    template <typename U> bool equal(details::sentinel_t<U> const &s) const { return (s.it == std::get<sizeof...(It) - 1>(its)); }

    template <size_t... Is> auto tuple_map_impl(std::index_sequence<Is...>) const {
      return std::tuple<decltype(*std::get<Is>(its))...>(*std::get<Is>(its)...);
    }
    decltype(auto) dereference() const { return tuple_map_impl(std::index_sequence_for<It...>{}); }
  };

  /********************* The Wrapper Classes representing the adapted ranges ********************/

  namespace details {

    template <typename T, typename L> struct transformed {
      T x;
      L lambda;

      using const_iterator = transform_iter<std::decay_t<decltype(std::cbegin(x))>, std::decay_t<L>>;
      using iterator       = const_iterator;

      const_iterator cbegin() const { return {std::cbegin(x), lambda}; }
      const_iterator begin() const { return cbegin(); }

      const_iterator cend() const { return {std::cend(x), lambda}; }
      const_iterator end() const { return cend(); }
    };

    // ---------------------------------------------

    template <typename T> struct enumerated {
      T x;

      using iterator       = enum_iter<std::decay_t<decltype(std::begin(x))>>;
      using const_iterator = enum_iter<std::decay_t<decltype(std::cbegin(x))>>;

      iterator begin() { return std::begin(x); }
      const_iterator cbegin() const { return std::cbegin(x); }
      const_iterator begin() const { return cbegin(); }

      iterator end() { return std::end(x); }
      const_iterator cend() const { return std::cend(x); }
      const_iterator end() const { return cend(); }
    };

    // ---------------------------------------------

    template <typename... T> struct zipped {
      std::tuple<T...> tu; // T can be a ref.

      using seq_t          = std::index_sequence_for<T...>;
      using iterator       = zip_iter<std::decay_t<decltype(std::begin(std::declval<T &>()))>...>;
      using const_iterator = zip_iter<std::decay_t<decltype(std::cbegin(std::declval<T &>()))>...>;

      template <typename... U> zipped(U &&... ranges) : tu{std::forward<U>(ranges)...} {}

      // Apply function to tuple
      template <typename F, size_t... Is> auto tuple_map(F &&f, std::index_sequence<Is...>) {
        return iterator{std::make_tuple(f(std::get<Is>(tu))...)};
      }
      template <typename F, size_t... Is> auto tuple_map(F &&f, std::index_sequence<Is...>) const {
        return const_iterator{std::make_tuple(f(std::get<Is>(tu))...)};
      }

      iterator begin() {
        return tuple_map([](auto &&x) { return std::begin(x); }, seq_t{});
      }
      const_iterator cbegin() const {
        return tuple_map([](auto &&x) { return std::cbegin(x); }, seq_t{});
      }
      const_iterator begin() const { return cbegin(); }

      iterator end() {
        return tuple_map([](auto &&x) { return std::end(x); }, seq_t{});
      }
      const_iterator cend() const {
        return tuple_map([](auto &&x) { return std::cend(x); }, seq_t{});
      }
      const_iterator end() const { return cend(); }
    };

    // ---------------------------------------------

    template <typename... T> struct multiplied {
      std::tuple<T...> tu; // T can be a ref.

      using iterator       = prod_iter<std::decay_t<decltype(std::begin(std::declval<T &>()))>...>;
      using const_iterator = prod_iter<std::decay_t<decltype(std::cbegin(std::declval<T &>()))>...>;

      template <typename... U> multiplied(U &&... ranges) : tu{std::forward<U>(ranges)...} {}

      template <size_t... Is> auto _begin(std::index_sequence<Is...>) {
        return iterator{std::make_tuple(std::begin(std::get<Is>(tu))...), std::make_tuple(std::end(std::get<Is>(tu))...)};
      }

      template <size_t... Is> auto _cbegin(std::index_sequence<Is...>) const {
        return const_iterator{std::make_tuple(std::cbegin(std::get<Is>(tu))...), std::make_tuple(std::cend(std::get<Is>(tu))...)};
      }

      iterator begin() { return _begin(std::index_sequence_for<T...>{}); }
      const_iterator cbegin() const { return _cbegin(std::index_sequence_for<T...>{}); }
      const_iterator begin() const { return cbegin(); }

      auto end() { return make_sentinel(std::end(std::get<sizeof...(T) - 1>(tu))); }
      auto cend() const { return make_sentinel(std::cend(std::get<sizeof...(T) - 1>(tu))); }
      auto end() const { return cend(); }
    };

    template <typename... T> multiplied(T &&...)->multiplied<std::decay_t<T>...>;

  } // namespace details

  /********************* The range adapting functions ********************/

  /**
   * Transform (lazy)applies a unary lambda function to every
   * element of a range. It returns itself a range.
   *
   * @param range The range that the lambda is applied to
   * @param range The lambda to apply to the range
   */
  template <typename T, typename L> auto transform(T &&range, L &&lambda) {
    return details::transformed<T, L>{std::forward<T>(range), std::forward<L>(lambda)};
  }

  /**
   * Lazy-enumerate a range.
   * This function returns itself a range of tuple<int, T>
   *
   * @param range The range to enumerate
   */
  template <typename T> auto enumerate(T &&range) { return details::enumerated<T>{std::forward<T>(range)}; }

  /**
   * Lazy-zip a range.
   * This function returns itself a range of tuple<T...>
   *
   * @param ranges The ranges to zip. Note: They have to be of equal length!
   */
  template <typename... T> auto zip(T &&... ranges) { return details::zipped<T...>{std::forward<T>(ranges)...}; }
  template <typename... T, typename L> auto zip_with(T &&... ranges, L &&lambda) {
    return transform(zip(std::forward<T>(ranges)...), [lambda](std::tuple<T...> t) { return std::apply(lambda, t); });
  }

  /**
   * Lazy-product of multiple ranges. This function returns itself a range of tuple<T...>.
   * Iterating over it will yield all combinations of the different range values.
   * Note: The ranges are incremented beginning with the leftmost range.
   *
   * @tparam T The types of the different ranges
   * @param ranges The ranges to zip. Note: They have to be of equal length!
   */
  template <typename... T> auto product(T &&... ranges) { return details::multiplied<T...>{std::forward<T>(ranges)...}; }


  /********************* Some factory functions ********************/

  namespace details {
    template <typename T, size_t N, size_t... Is> auto make_product_impl(std::array<T, N> &arr, std::index_sequence<Is...>) {
      static_assert(N == sizeof...(Is));
      return product(arr[Is]...);
    }
  } // namespace details

  template <typename T, size_t N> auto make_product(std::array<T, N> &arr) { return details::make_product_impl(arr, std::make_index_sequence<N>{}); }

  /**
   * Create a std::vector from a range
   */
  template <typename R> auto make_vector_from_range(R &&r) {
    std::vector<std::decay_t<decltype(*(std::begin(r)))>> vec;
    for (auto const &x : r) vec.push_back(x);
    return vec;
  }

} // namespace triqs::utility
