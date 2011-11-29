#ifndef __ISNAN_H__
#define __ISNAN_H__

#include "wcsconfig.h"

typedef unsigned WCSLIB_INT64 Int64;

#if !defined(U64)
#define U64(u) (* (Int64 *) &(u) )
#endif /* U64 */

#if !defined(isnan64)
#if !defined(_MSC_VER)
#define isnan64(u) \
  ( (( U64(u) & 0x7ff0000000000000LL)  == 0x7ff0000000000000LL)  && ((U64(u) &  0x000fffffffffffffLL) != 0)) ? 1:0
#else
#define isnan64(u) \
  ( (( U64(u) & 0x7ff0000000000000i64) == 0x7ff0000000000000i64)  && ((U64(u) & 0x000fffffffffffffi64) != 0)) ? 1:0
#endif
#endif /* isnan64 */

#if !defined(isinf64)
#if !defined(_MSC_VER)
#define isinf64(u) \
  ( (( U64(u) & 0x7ff0000000000000LL)  == 0x7ff0000000000000LL)  && ((U64(u) &  0x000fffffffffffffLL) == 0)) ? 1:0
#else
#define isinf64(u) \
  ( (( U64(u) & 0x7ff0000000000000i64) == 0x7ff0000000000000i64)  && ((U64(u) & 0x000fffffffffffffi64) == 0)) ? 1:0
#endif
#endif /* isinf64 */

#if !defined(isfinite64)
#if !defined(_MSC_VER)
#define isfinite64(u) \
  ( (( U64(u) & 0x7ff0000000000000LL)  != 0x7ff0000000000000LL)) ? 1:0
#else
#define isfinite64(u) \
  ( (( U64(u) & 0x7ff0000000000000i64) != 0x7ff0000000000000i64)) ? 1:0
#endif
#endif /* isfinite64 */

#endif /* __ISNAN_H__ */
