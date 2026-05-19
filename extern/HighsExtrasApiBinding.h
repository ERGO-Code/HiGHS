/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsExtrasApiBinding.h
 * @brief Provides useful classes for handling external dependencies
 */
#ifndef HIGHS_EXTRAS_API_BINDING_H_
#define HIGHS_EXTRAS_API_BINDING_H_

#include <tuple>

// provide metadata info for each feature
struct HighsExtrasFeatureInfo {
  HighsExtrasFeatureInfo(const char* provider_ = nullptr,
                         const char* version_ = nullptr,
                         const char* license_ = nullptr, bool enabled_ = false)
      : provider(provider_),
        version(version_),
        license(license_),
        enabled(enabled_) {}

  const char* provider;
  const char* version;
  const char* license;
  const bool enabled;
};

namespace HighsExtras {

template <class... Features>
struct require {};

template <class Methods>
struct feature_api;

// convenience wrapper to access the HighsExtrasApi storage
template <class Family>
struct wrapper_storage {
  template <class Methods>
  static feature_api<Methods>& getApi();

  static const HighsExtrasFeatureInfo* getInfo() { return nullptr; };
};

template <class FeatureFamily, int FeatureIndex>
struct feature_base {
  using family = FeatureFamily;

  static const HighsExtrasFeatureInfo* getInfo() {
    const HighsExtrasFeatureInfo* ptr = wrapper_storage<family>::getInfo();
    return ptr ? ptr + FeatureIndex : nullptr;
  }
};

// captures the function to bind and its function pointer type
// does not store anything itself
template <typename T>
struct method_desc {
  using fnptr_t = T;
};

template <typename T, T value>
struct bound_method_desc : method_desc<T> {
  static typename method_desc<T>::fnptr_t direct() { return value; }
};

// clang wants to link the function, even though it's not used
// so only provide bound method when building highs_extras
#if defined(HIGHS_EXTRAS_LIBRARY_BUILD)
#define HIGHS_API_DESC(fn) bound_method_desc<decltype(&fn), &fn>
#else
#define HIGHS_API_DESC(fn) method_desc<decltype(&fn)>
#endif

// storage for the function pointer, given a method_desc<...>
template <class Desc>
struct method_storage {
  typename Desc::fnptr_t value;
  method_storage() : value(nullptr) {}
};

// builds a struct of function pointers, given tuple<method_desc<...>, ...>
template <class... Desc>
struct feature_api<std::tuple<Desc...>> : method_storage<Desc>... {
  using methods_type = std::tuple<Desc...>;

  // access method by index at compile-time, e.g., api->method<0>(...)
  template <std::size_t Index>
  typename std::tuple_element<Index, methods_type>::type::fnptr_t& method() {
    using desc_type = typename std::tuple_element<Index, methods_type>::type;
    return static_cast<method_storage<desc_type>&>(*this).value;
  }
};

// access function pointer by index, e.g., api::fn<0>()(...)
template <class Family, class Methods>
struct feature_wrapper {
  template <std::size_t Index>
  static typename std::tuple_element<Index, Methods>::type::fnptr_t& fn() {
    return wrapper_storage<Family>::template getApi<Methods>()
        .template method<Index>();
  }
};

// template recursion to set function pointer to internal method
template <class Methods, std::size_t Index, std::size_t Count>
struct bind_methods {
  static void apply(feature_api<Methods>& api) {
    using desc_type = typename std::tuple_element<Index, Methods>::type;
    api.template method<Index>() = desc_type::direct();
    bind_methods<Methods, Index + 1, Count>::apply(api);
  }
};

// specialization to terminate recursion
template <class Methods, std::size_t Count>
struct bind_methods<Methods, Count, Count> {
  static void apply(feature_api<Methods>&) {}
};

// recursively set function pointers to the direct methods
template <class Methods>
void bind_api(feature_api<Methods>& api) {
  bind_methods<Methods, 0, std::tuple_size<Methods>::value>::apply(api);
}

}  // namespace HighsExtras

#endif  // HIGHS_EXTRAS_API_BINDING_H_