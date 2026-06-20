// In-house Sobol' sequence generator with Burley hash-based Owen scrambling.
//
// References:
//   * Joe & Kuo (2008). "Constructing Sobol' sequences with better
//     two-dimensional projections." SIAM J. Sci. Comput. 30(5):2635-2654.
//   * Burley (2020). "Practical Hash-based Owen Scrambling."
//     Journal of Computer Graphics Techniques 9(4).
//   * Laine & Karras (2011). "Stratified sampling for stochastic
//     transparency." Computer Graphics Forum 30(4):1197-1204.
//
// Direction numbers (file `sobol_direction_numbers.h`) are the
// public-domain "new-joe-kuo-6.21201" table from
// https://web.maths.unsw.edu.au/~fkuo/sobol/, subsetted to the first
// 250 dimensions (sufficient for realistic sensobol use).

#include <Rcpp.h>
#include <cstdint>
#include <vector>
#include "sobol_direction_numbers.h"

namespace {

constexpr unsigned L = 32u;  // bit width of the Sobol' state

inline std::uint32_t reverse_bits(std::uint32_t x) {
  x = ((x & 0xffff0000u) >> 16) | ((x & 0x0000ffffu) << 16);
  x = ((x & 0xff00ff00u) >>  8) | ((x & 0x00ff00ffu) <<  8);
  x = ((x & 0xf0f0f0f0u) >>  4) | ((x & 0x0f0f0f0fu) <<  4);
  x = ((x & 0xccccccccu) >>  2) | ((x & 0x33333333u) <<  2);
  x = ((x & 0xaaaaaaaau) >>  1) | ((x & 0x55555555u) <<  1);
  return x;
}

// Laine-Karras hash (a permutation of 32-bit integers).
inline std::uint32_t laine_karras_hash(std::uint32_t x, std::uint32_t seed) {
  x += seed;
  x ^= x * 0x6c50b47cu;
  x ^= x * 0xb82f1e52u;
  x ^= x * 0xc7afe638u;
  x ^= x * 0x8d22f6e6u;
  return x;
}

// Burley's hash-based Owen scrambling in base 2.
// Equivalent in distribution to the original Owen scramble, with O(1)
// per-point work and no per-point permutation tables.
inline std::uint32_t owen_scramble(std::uint32_t x, std::uint32_t seed) {
  x = reverse_bits(x);
  x = laine_karras_hash(x, seed);
  x = reverse_bits(x);
  return x;
}

// Compute the L direction integers V[1..L] for dimension d (1-indexed).
//   - d == 1: van der Corput / 1-d Sobol' sequence  => V[i] = 1 << (L - i)
//   - d >= 2: Antonov-Saleev recurrence over Joe-Kuo direction numbers
std::vector<std::uint32_t> compute_V(unsigned d) {
  std::vector<std::uint32_t> V(L + 1, 0u);

  if (d == 1u) {
    for (unsigned i = 1; i <= L; ++i) {
      V[i] = 1u << (L - i);
    }
    return V;
  }

  // d >= 2: look up s, a, m from the embedded table.
  const std::size_t idx = static_cast<std::size_t>(d) - 2u;
  const unsigned s = SOBOL_S[idx];
  const unsigned a = SOBOL_A[idx];

  for (unsigned i = 1; i <= s; ++i) {
    const std::uint32_t mi = SOBOL_M[idx * SOBOL_MAX_S + (i - 1)];
    V[i] = mi << (L - i);
  }

  for (unsigned i = s + 1; i <= L; ++i) {
    V[i] = V[i - s] ^ (V[i - s] >> s);
    for (unsigned k = 1; k + 1 <= s; ++k) {
      const unsigned bit = (a >> (s - 1u - k)) & 1u;
      if (bit) V[i] ^= V[i - k];
    }
  }

  return V;
}

// Position of the rightmost zero bit of n, 1-indexed (so n=0 -> 1, n=1 -> 2,
// n=2 -> 1, n=3 -> 3, n=4 -> 1, ...). This is the Gray-code position used
// in the Antonov-Saleev recurrence X_{n+1} = X_n XOR V[c(n)+1].
inline unsigned rightmost_zero_bit_position(std::uint32_t n) {
  unsigned c = 1;
  while (n & 1u) { n >>= 1; ++c; }
  return c;
}

}  // anonymous namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix sobol_owen_cpp(int N, int dim,
                                   Rcpp::Nullable<int> seed_in,
                                   bool scramble) {

  if (N < 1) Rcpp::stop("'N' must be a positive integer.");
  if (dim < 1 || static_cast<std::size_t>(dim) > SOBOL_MAX_DIM) {
    Rcpp::stop("'dim' must be between 1 and %d (the maximum dimension "
               "supported by the in-house Joe-Kuo direction-number table).",
               static_cast<int>(SOBOL_MAX_DIM));
  }

  // Per-dimension Owen scrambling seed. If seed_in is NULL we draw a
  // random master seed from R's RNG so that the user still gets a
  // *random* Owen-scrambled sequence; if seed_in is supplied the run
  // is fully reproducible.
  std::uint32_t master_seed = 0u;
  if (scramble) {
    if (seed_in.isNotNull()) {
      master_seed = static_cast<std::uint32_t>(Rcpp::as<int>(seed_in));
    } else {
      // R's runif returns [0,1); rescale to a 32-bit integer.
      Rcpp::NumericVector u = Rcpp::runif(1);
      master_seed = static_cast<std::uint32_t>(u[0] * 4294967296.0);
    }
  }

  Rcpp::NumericMatrix out(N, dim);
  const double scale = 1.0 / 4294967296.0;  // 2^-32

  for (int d = 1; d <= dim; ++d) {
    const std::vector<std::uint32_t> V = compute_V(static_cast<unsigned>(d));

    // Derive a per-dimension scramble seed from the master seed.
    // Each dimension must use an *independent* permutation, hence the hash.
    const std::uint32_t dim_seed = scramble
      ? laine_karras_hash(static_cast<std::uint32_t>(d), master_seed)
      : 0u;

    std::uint32_t X = 0u;
    for (int n = 0; n < N; ++n) {
      const unsigned c = rightmost_zero_bit_position(static_cast<std::uint32_t>(n));
      X ^= V[c];
      const std::uint32_t Xs = scramble ? owen_scramble(X, dim_seed) : X;
      out(n, d - 1) = static_cast<double>(Xs) * scale;
    }
  }

  return out;
}
