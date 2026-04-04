# MPInt Optimization Roadmap

## Context

`MPInt` is an arbitrary-precision integer implemented as `class val` with `_digits: Array[U16] val`
(little-endian, base 65536) and `_negative: Bool`. It implements `SignedInteger[MPInt, MPInt]`.

**The central problem:** three multiplication algorithms exist (`mul` schoolbook O(n²), `karatsuba_mul`
O(n^1.585), `fast_mul` FFT O(n log n)) but **`mul()` always uses schoolbook** — there is no automatic
dispatch. All callers in `mpfloat.pony` go through `mul()` even for 500-digit numbers.

**Secondary problem:** `string()` and `from_string()` are O(n²) via repeated `_short_div(10)` /
`_short_mul` loops.

---

**User addition**: Always respect rules in `STYLE_GUIDE.md`.

## Optimizations, Ranked by Performance Impact

### 1. Automatic Multiplication Dispatch ★★★★★
**Easy — Maximum ROI**

Replace `mul()` body with a size-dispatched router:
```
n < ~32 digits  → _mul_schoolbook (current mul body)
n < ~512 digits → _karatsuba_mul  (current karatsuba_mul body, renamed)
n ≥ 512 digits  → _ntt_mul or _fft_mul
```
Rename the three existing algorithm bodies to private helpers. The public `mul()` signature
(part of `SignedInteger`) is unchanged.

- **Impact:** 5–15× for 150-digit numbers (Karatsuba vs schoolbook); 100–1000× for 2000+ digits
  (NTT/FFT vs schoolbook). Benefits all callers automatically: `pow`, `pow_mod`, `isqrt`, all
  `MPFloat` arithmetic.
- **Advantage:** Zero API change. Existing tests cover all three algorithms. Correctness risk is nil.
- **Disadvantage:** Thresholds need empirical tuning; a wrong threshold adds overhead for small numbers.
- **Pony:** No capability issues. All algorithms already use `recover val`.

**User additions**:
- Multiplication algorithm methods must remain public so the client can decide to use a specific algorithm.
- Keep consistant naming with prefix `nul`: `mul_schoolbook`, `mul_karatsuba`, `mul_ntt`, `mul_ftt`
- Find literature for common thresholds.

---

### 2. NTT-Based Multiplication (replacing FFT) ★★★★☆
**Medium — Removes precision limit, enables arbitrary-size exact multiplication**

Implement `_ntt_mul()` using the existing `NTT[USize]` primitive as the convolution engine:
1. Zero-pad both `_digits` arrays to next power of 2 (`m`).
2. Convert U16 limbs to `Array[USize]`.
3. `NTT[USize].transform(a, false)` and `transform(b, false)`.
4. Pointwise: `c(i) = Modular[USize].mul_mod(a(i), b(i), NTT._p())`.
5. `NTT[USize].transform(c, true)`.
6. Carry-propagate back to U16 limbs.

**NTT prime overflow caveat:** On LP64 (Linux/macOS 64-bit), NTT uses the 32-bit prime `998244353`.
Maximum convolution coefficient = n × 65535² ≈ 4.3 × 10⁹ × n, which exceeds `p` for n > 1.
**Fix options:**
- Use the 64-bit prime `2⁶⁴ − 2³² + 1` (currently gated on `llp64` in `ntt.pony` — needs updating
  to also trigger on `lp64`).
- Or three-NTT CRT with three 32-bit primes whose product > max coefficient.

- **Impact:** Exact for arbitrary sizes (no F64 precision limit); ~1.5–3× slower than FFT per call
  but correct. Required for numbers beyond ~1M digits.
- **Advantage:** Reuses existing tested `NTT` code (~60–80 lines of glue). Integer-exact, no rounding.
- **Disadvantage:** LP64 prime selection needs a fix. 32-bit prime alone insufficient for U16 limbs.
- **Pony:** Standard `recover val`; no new issues. `NTT._p()` and `NTT._g()` are `fun` on a primitive.

**User addition**
- `USize` being a variable size unsigned depending on platform, use a fixed size integer like `U32` or `U64` instead.

---

### 3. Divide-and-Conquer String Conversion ★★★★☆
**Hard — Fixes O(n²) → O(n log² n) for `string()` and `from_string()`**

Current: `string()` calls `_short_div(10)` in a loop — O(n) per digit, O(n²) total.

Algorithm: precompute `10^1, 10^2, 10^4, ..., 10^(2^k)` (each an `MPInt`). To convert:
- `(q, r) = divrem(10^(n/2))`
- Recursively convert `q` (high half) and `r` (low half), then concatenate.
Each level does one big `divrem`; O(log n) levels at O(n log n) per level = O(n log² n) total.
Requires fast multiplication to be useful (pairs with Rank 1).

- **Impact:** ~6× for n=500 decimal digits; ~180× for n=10000 digits. Directly improves
  `MPFloat.exact_string()` and `MPFloat.from_string()`.
- **Advantage:** Fundamental complexity improvement; most impactful for large-precision MPFloat work.
- **Disadvantage:** Significant implementation complexity. Need size threshold (slower than linear
  scan below ~50 decimal digits). Powers of 10 table must be computed per call.
- **Pony:** Recursive builder inside `recover iso`; `String iso^` assembled via `append`.

---

### 4. U32 Limbs ★★★☆☆
**Medium — Foundational data structure improvement**

Change `_digits: Array[U16] val` → `Array[U32] val`, base `2³²`.

- Two U32 values multiply to U64 — native in Pony (`a.u64() *~ b.u64()`), no 128-bit needed.
- The U16 constraint came from FFT's F64 precision; with NTT (Rank 2), it disappears.
- Every helper (`_add_arrays`, `_sub_arrays`, `_short_mul`, `_short_div`, `_mul_schoolbook`,
  Algorithm D) must be updated to use U64 intermediates and 32-bit carry masks.
- `BitMap` interop (`_bitwise_op`) needs `from_u32_array`/`to_u32_array` analogues.
- `raw_digits()` output (U8 bytes, big-endian) is unaffected in format; extraction loop changes.

- **Impact:** Halves array length → 2–4× on add/sub/shift (loop count); 4× on schoolbook mul
  (n/2 × n/2 = n²/4 iterations); 2× on NTT array size.
- **Advantage:** Better cache utilization; orthogonal to all other optimizations; natural prerequisite
  for U64 (Rank 6).
- **Disadvantage:** ~60–80 sites in `mpint.pony` to update; mechanical but tedious; high regression
  risk without thorough test coverage. The existing 179 passing tests are the safety net.
- **Pony:** No capability issues. `Array[U32] val` behaves identically to `Array[U16] val`.

---

### 5. Remove Redundant Array Copy in `mul()` ★★☆☆☆
**Easy — Small but free win**

In `mul()` and `fast_mul()`, after the `recover val` block, the result array is copied
into `d2` just to get a right-sized array:
```pony
let d2 = recover val
  let res2 = Array[U16](res.size())
  for x in res.values() do res2.push(x) end  // ← UNNECESSARY
  res2
end
```
`_normalize(res)` already shrinks the array in place. Remove `d2`; use `consume res` directly.

- **Impact:** Eliminates 1 array allocation + O(n) copy per multiplication. ~5–15% throughput
  improvement on multiplication chains.
- **Advantage:** One-line change per operation; zero risk.
- **Pony:** `consume res` inside the `recover val` block produces the `val` array directly.

---

### 6. U64 Limbs ★★★☆☆
**Hard — Maximum efficiency on 64-bit hardware**

Change `_digits: Array[U64] val`, base `2⁶⁴`. Requires U128 intermediates in the schoolbook
inner loop (`a.u128() *~ b.u128()`). Algorithm D's `q_hat` estimation step becomes more complex.
Pony emulates U128 with two 64-bit instructions (~2–4× slower than U64), partially canceling the
limb-count reduction.

- **Impact:** 2× add/sub/shift over U32; 2× NTT array size; schoolbook gain marginal over U32.
- **Recommendation:** Do U32 first (Rank 4); benchmark; only proceed to U64 if data supports it.

**User addition**
- 64-bits platforms are being mainstream, and even Pony on 32-bits systems support `F128` types...
- By the way, rename `_negative` to `_sign` to be consistant with `MPFloat`.

---

### 7. Newton-Raphson Division ★★★☆☆
**Hard — O(n²) → O(n log n) for large division**

Replace Knuth Algorithm D with: compute reciprocal `t ≈ 1/v` via Newton iterations
(`t_{k+1} = t_k × (2 − v × t_k)`), starting from 1 limb, doubling precision each step.
Then `q = trunc(u × t)` with one correction step. With NTT multiplication, division becomes O(n log n).

- **Impact:** ~55× for n=500-digit balanced division. Cascades to `divrem`, `div`, `rem`,
  `gcd`, `isqrt`, string conversion.
- **Disadvantage:** Requires Ranks 1+2 first. Off-by-one errors produce wrong quotients that
  are hard to detect without systematic testing.

---

### 8. Toom-Cook 3-Way ★★☆☆☆
**Medium-Hard — Better intermediate between Karatsuba and NTT**

Reduce n-digit mul to 5 sub-multiplications of size n/3 vs Karatsuba's 3 of size n/2.
O(n^1.465) vs O(n^1.585). Typically beats Karatsuba above ~300 limbs; NTT beats Toom-3
above ~1000–2000 limbs.

- **Impact:** 1.5–2× over Karatsuba for 150–300-digit range; 2–4× for 200–500-digit range.
- **Recommendation:** Low priority if NTT (Rank 2) is implemented; NTT supersedes it at large sizes.

**User addition**
- Method must be public and have a consistant naming with other multiplication algorithms.

---

### 9. Montgomery Multiplication for `pow_mod` ★★☆☆☆
**Hard — 2× speedup for large modular exponentiation**

Precompute `R = 2^k`, `m' = -m⁻¹ mod R`. Each Montgomery multiply computes
`(a × b × R⁻¹) mod m` with shifts+adds only — no division. Saves 1 `divrem` per
`pow_mod` iteration.

- **Impact:** ~2× for `pow_mod` with 1000-digit modulus.
- **Disadvantage:** Only benefits `pow_mod`; complex. Code comment warns against cryptographic use.

---

### 10. Barrett Reduction ★★☆☆☆
**Medium — Amortized cost for repeated `mod m` with same modulus**

Precompute `k = ⌊2^(2n) / m⌋` once; each subsequent `x mod m` costs 2 multiplications
instead of 1 `divrem`.

- **Impact:** 2× for `pow_mod`; no benefit for one-shot `mod`.
- **Pony constraint:** `class val MPInt` cannot cache the precomputed `k` per-instance.
  Must expose as a separate `BarrettContext val` object: `fun pow_mod_barrett(exp, m, ctx: BarrettContext val): MPInt`.
- **Advantage:** Simpler than Montgomery.

---

### 11. Index-Based Karatsuba (Eliminate Sub-Array Copies) ★★☆☆☆
**Medium — Eliminates 4 array allocations per Karatsuba level**

Current Karatsuba creates `this_low`, `this_high`, `that_low`, `that_high` as separate MPInt
objects (each a fresh array). Replace with a private recursive helper that takes
`(array, offset, length)` tuples and accesses limbs by index.

- **Impact:** Saves ~8 allocations for a 500-digit Karatsuba (2 levels). ~5–15% Karatsuba throughput.
- **Pony:** A `val` array can be passed to any context; bounds-checked `_digits(offset+k)?` access.

---

### 12. Lehmer / Binary GCD ★★☆☆☆
**Hard — O(n³) → O(n² log n) or better for GCD**

Lehmer's algorithm uses small-precision leading-limb approximations to extract multiple
Euclidean quotients per big-number step, reducing per-step big-division count from O(n) to
O(n / log n). With fast mul/div, total GCD cost becomes O(n log n × M(n)).

- **Impact:** 3–10× for n > 500 digits.
- **Disadvantage:** GCD is not the primary bottleneck in observed MPFloat workflows.

---

### 13. `_MPInt ref` Proxy ✗
**Do not implement — fundamentally incompatible with Pony capabilities**

A `class val` cannot hold a path to mutable (`ref`) state. Any proxy design either:
- (a) Wraps an `iso` immediately consumed to `val` — equivalent to the current `recover val`
  pattern with no benefit, or
- (b) Changes the API contract away from `class val`, breaking `SignedInteger`.

The `recover val` blocks already provide all the local mutability needed. The cost is
allocation, not mutation semantics. The proxy adds indirection without reducing allocation count.

---

## Implementation Priority

| Phase | Tasks | Effort |
|-------|-------|--------|
| **1 — Quick wins** | Rank 5 (remove redundant copy), Rank 1 (dispatch) | 1–2 days |
| **2 — Data foundation** | Rank 4 (U32 limbs) | 3–5 days |
| **3 — Asymptotic win** | Rank 2 (NTT mul) + fix `lp64` prime in `ntt.pony` | 1–2 weeks |
| **4 — String bottleneck** | Rank 3 (D&C string conv) | 1–2 weeks |
| **5 — Optional** | Ranks 8, 7, 9/10, 12, 6 — based on benchmarks | ongoing |

## Critical Files

- `mathx/mpint.pony` — all changes land here
- `mathx/ntt.pony` — `lp64` prime fix for Rank 2 (line 31: `ifdef ilp32 or lp64`)
- `tests/_tests_mpint.pony` — existing tests cover all 3 mul algorithms; add benchmarks
- `mathx/mpfloat.pony` — main beneficiary; no changes needed (uses `MPInt` via interface)
