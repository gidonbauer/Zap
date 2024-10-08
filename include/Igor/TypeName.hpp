// Copyright 2024 Gidon Bauer <gidon.bauer@rwth-aachen.de>

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:

// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef IGOR_TYPE_NAME_HPP_
#define IGOR_TYPE_NAME_HPP_

#include <cstdlib>
#include <memory>
#include <string>
#include <type_traits>

#ifndef IGOR_NO_CXX_ABI
#include <cxxabi.h>
#endif  // IGOR_NO_CXX_ABI

namespace Igor {

#ifndef IGOR_NO_CXX_ABI

template <typename T>
[[nodiscard]] constexpr auto type_name() -> std::string {
  using namespace std::string_literals;

  int status;
  constexpr auto free_deleter = [](void* p) constexpr noexcept {
    std::free(p);  // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc)
  };
  std::unique_ptr<char, decltype(free_deleter)> name_cstr{
      abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status), free_deleter};

  if (status != 0 || name_cstr == nullptr) {
    switch (status) {
      case -1:
        throw std::runtime_error(
            "Demagleing failed with status -1: A memory allocation failure occurred.");
      case -2:
        throw std::runtime_error("Demagleing failed with status -2: mangled_name is not a valid "
                                 "name under the C++ ABI mangling rules.");
      case -3:
        throw std::runtime_error(
            "Demagleing failed with status -3: One of the arguments is invalid.");
      default:
        throw std::runtime_error("Demagleing failed with unknown status "s +
                                 std::to_string(status) + "."s);
    }
  }

  std::string name{name_cstr.get()};
  if (std::is_volatile_v<std::remove_reference_t<T>>) {
    name += " volatile"s;
  }
  if (std::is_const_v<std::remove_reference_t<T>>) {
    name += " const"s;
  }
  if (std::is_lvalue_reference_v<T>) {
    name += "&"s;
  }
  if (std::is_rvalue_reference_v<T>) {
    name += "&&"s;
  }

  return name;
}

template <typename T>
[[nodiscard]] constexpr auto type_name(T /*ignored*/) -> std::string {
  return type_name<T>();
}

#endif  // IGOR_NO_CXX_ABI
}  // namespace Igor

#endif  // IGOR_TYPE_NAME_HPP_
