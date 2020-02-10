#ifndef COMMON_HPP
#define COMMON_HPP

#include <cstdio>
#include <cstdarg>
#include <exception>

#define noinline_fun __attribute__((noinline))
#define forceinline_fun inline __attribute__((always_inline))
#define flatten_fun __attribute__((flatten))
#define pure_fun __attribute__((pure))
#define hot_fun __attribute__((hot))
#define cold_fun __attribute__((cold))
#define restrict __restrict__
#define likely(expr) __builtin_expect(!!(expr), 1)
#define unlikely(expr) __builtin_expect(!!(expr), 0)
#define prefetchr(addr) __builtin_prefetch((addr), 0)
#define prefetchw(addr) __builtin_prefetch((addr), 1)

#ifndef __has_feature
#    define __has_feature(x) (0)
#endif

#ifdef __has_cpp_attribute
#    if __has_cpp_attribute(noreturn)
#        define noreturn_attr [[noreturn]]
#    endif
#    if __has_cpp_attribute(nodiscard)
#        define nodiscard_attr [[nodiscard]]
#    elif __has_cpp_attribute(gnu::warn_unused_result)
#        define nodiscard_attr [[gnu::warn_unused_result]]
#    endif
#endif

#ifndef noreturn_attr
#    define noreturn_attr __attribute__((noreturn))
#endif

#ifndef nodiscard_attr
#    define nodiscard_attr __attribute__((warn_unused_result))
#endif

#ifdef assert
#    undef assert
#endif

#ifdef NDEBUG
#    define DEBUG 0
#    define assert(...) (static_cast<void>(0))
#    define assume(expr, ...) (likely((expr)) ? static_cast<void>(0) : __builtin_unreachable())
#else
#    define DEBUG 1


/** Handler for assertion/assumption fails
 * Do not throw exception on purpose (directly terminate)
 */
noreturn_attr noinline_fun inline void
abort_message(const char msg[]...)
{
    va_list args;
    va_start(args, msg);
    std::vfprintf(stderr, msg, args);
    va_end(args);
    std::fflush(stderr);
    std::abort();
}


#    define sourceloc_fail(what, msg, ...) abort_message("%s:%u:%s\n\t" what " failed: " msg "\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, ##__VA_ARGS__))
//#    define sourceloc_fail1(what) abort_message("%s:%s:%u: " #what " failed.\n", __FILE__, __PRETTY_FUNCTION__, __LINE__))
//#define __sourceloc_fail_nargs(

#    define assert(expr, ...) (likely((expr)) ? static_cast<void>(0) : sourceloc_fail("Assertion '" #expr "'", ##__VA_ARGS__)
#    define assume(expr, ...) (likely((expr)) ? static_cast<void>(0) : sourceloc_fail("Assumption '" #expr "'", ##__VA_ARGS__)

#endif

#endif // COMMON_HPP
