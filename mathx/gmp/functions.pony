// Base functions for MPFR


use "lib:mpfr"

use "../../assertx"
use "../../formatx"


type RoundingMode is (RoundingNearest | RoundingNegInf | RoundingPosInf |
                      RoundingZero | RoundingAwayZ | RoundingFaithful )
  """
  The rounding mode that must be applied to operations.

  * MPFR_RNDN: round to nearest, with the even rounding rule (roundTiesToEven in IEEE 754).
  * MPFR_RNDD: round toward negative infinity (roundTowardNegative in IEEE 754).
  * MPFR_RNDU: round toward positive infinity (roundTowardPositive in IEEE 754).
  * MPFR_RNDZ: round toward zero (roundTowardZero in IEEE 754).
  * MPFR_RNDA: round away from zero.
  * MPFR_RNDF: faithful rounding. This feature is currently experimental.

  For full definition: https://www.mpfr.org/mpfr-current/mpfr.html#Rounding
  """

primitive RoundingNearest
    """
    Round to the nearest, with the even rounding rule.
    """
  fun apply(): I32 =>
    0

  fun string(): String =>
    "MPFR_RNDN"

primitive RoundingNegInf
    """
    Round toward negative infinity.
    """
  fun apply(): I32 =>
    3

  fun string(): String =>
    "MPFR_RNDD"

primitive RoundingPosInf
    """
    Round toward positive infinity.
    """
  fun apply(): I32 =>
    2

  fun string(): String =>
    "MPFR_RNDU"

primitive RoundingZero
    """
    Round toward zero.
    """
  fun apply(): I32 =>
    1

  fun string(): String =>
    "MPFR_RNDZ"

primitive RoundingAwayZ
    """
    Round away from zero.
    """
  fun apply(): I32 =>
    4

  fun string(): String =>
    "MPFR_RNDA"

primitive RoundingFaithful
    """
    Faithfull rounding.
    """
  fun apply(): I32 =>
    5

  fun string(): String =>
    "MPFR_RNDF"


primitive MPF
  """
  The `MPF` primitive provides all base functions from the [MPFR](https://www.mpfr.org/)
  library. You can get the maximum performance and memory management when using
  these functions, but you must also manage the memory yourself. The `MPFloat`
  class offers transparent mapping to some of these functions while also
  providing ease of use.
  """

  fun _final() =>
    """
    Free MPFR cache when the functions are no more used.
    """
    @mpfr_free_cache()

  fun get_version(): String iso^ =>
    var pcstring: Pointer[U8] =  @mpfr_get_version()
    let p: String iso = String.from_cstring(pcstring).clone()
    consume p

  fun get_patches(): String iso^ =>
    var pcstring: Pointer[U8] =  @mpfr_get_patches()
    let p: String iso = String.from_cstring(pcstring).clone()
    consume p

  fun buildopt_tls_p(): I32 =>
    @mpfr_buildopt_tls_p()

  fun buildopt_float128_p(): I32 =>
    @mpfr_buildopt_float128_p()

  fun buildopt_decimal_p(): I32 =>
    @mpfr_buildopt_decimal_p()

  fun buildopt_gmpinternals_p(): I32 =>
    @mpfr_buildopt_gmpinternals_p()

  fun buildopt_sharedcache_p(): I32 =>
    @mpfr_buildopt_sharedcache_p()

  fun buildopt_tune_case(): String iso^ =>
    var pcstring: Pointer[U8] =  @mpfr_buildopt_tune_case()
    let p: String iso = String.from_cstring(pcstring).clone()
    consume p

  fun get_emin(): I64 =>
    @mpfr_get_emin()

  fun set_emin(parg0: I64): I32 =>
    @mpfr_set_emin(parg0)

  fun get_emin_min(): I64 =>
    @mpfr_get_emin_min()

  fun get_emin_max(): I64 =>
    @mpfr_get_emin_max()

  fun get_emax(): I64 =>
    @mpfr_get_emax()

  fun set_emax(parg0: I64): I32 =>
    @mpfr_set_emax(parg0)

  fun get_emax_min(): I64 =>
    @mpfr_get_emax_min()

  fun get_emax_max(): I64 =>
    @mpfr_get_emax_max()

  fun set_default_rounding_mode(rnd: RoundingMode) =>
    @mpfr_set_default_rounding_mode(rnd())

  fun get_default_rounding_mode(): RoundingMode =>
    match @mpfr_get_default_rounding_mode()
    | 1 => RoundingZero
    | 2 => RoundingPosInf
    | 3 => RoundingNegInf
    | 4 => RoundingAwayZ
    | 5 => RoundingFaithful
    else
      RoundingNearest
    end

  fun print_rnd_mode(rnd: RoundingMode): String iso^ =>
    var pcstring: Pointer[U8] =  @mpfr_print_rnd_mode(rnd())
    let p: String iso = String.from_cstring(pcstring).clone()
    consume p

  fun clear_flags() =>
    @mpfr_clear_flags()

  fun clear_underflow() =>
    @mpfr_clear_underflow()

  fun clear_overflow() =>
    @mpfr_clear_overflow()

  fun clear_divby0() =>
    @mpfr_clear_divby0()

  fun clear_nanflag() =>
    @mpfr_clear_nanflag()

  fun clear_inexflag() =>
    @mpfr_clear_inexflag()

  fun clear_erangeflag() =>
    @mpfr_clear_erangeflag()

  fun set_underflow() =>
    @mpfr_set_underflow()

  fun set_overflow() =>
    @mpfr_set_overflow()

  fun set_divby0() =>
    @mpfr_set_divby0()

  fun set_nanflag() =>
    @mpfr_set_nanflag()

  fun set_inexflag() =>
    @mpfr_set_inexflag()

  fun set_erangeflag() =>
    @mpfr_set_erangeflag()

  fun underflow_p(): Bool =>
    @mpfr_underflow_p() != 0

  fun overflow_p(): Bool =>
    @mpfr_overflow_p() != 0

  fun divby0_p(): Bool =>
    @mpfr_divby0_p() != 0

  fun nanflag_p(): Bool =>
    @mpfr_nanflag_p() != 0

  fun inexflag_p(): Bool =>
    @mpfr_inexflag_p() != 0

  fun erangeflag_p(): Bool =>
    @mpfr_erangeflag_p() != 0

  fun flags_clear(parg0: U32) =>
    @mpfr_flags_clear(parg0)

  fun flags_set(parg0: U32) =>
    @mpfr_flags_set(parg0)

  fun flags_test(parg0: U32): U32 =>
    @mpfr_flags_test(parg0)

  fun flags_save(): U32 =>
    @mpfr_flags_save()

  fun flags_restore(parg0: U32, parg1: U32) =>
    @mpfr_flags_restore(parg0, parg1)

  fun check_range(mpfr: SMPFr, parg1: I32, rnd: RoundingMode): I32 =>
    @mpfr_check_range(mpfr, parg1, rnd())

  fun init2(mpfr: SMPFr, prec: ILong) =>
    @mpfr_init2(mpfr, prec)

  fun init(mpfr: SMPFr) =>
    @mpfr_init(mpfr)

  fun clear(mpfr: SMPFr) =>
    @mpfr_clear(mpfr)

/*
  // Pony does not support variadic functions
  fun inits2(parg0: I64, parg1: SMPFr, ...) =>
    @mpfr_inits2(parg0, parg1, ...)

  fun inits(mpfr: SMPFr, ...) =>
    @mpfr_inits(mpfr, ...)

  fun clears(mpfr: SMPFr, ...) =>
    @mpfr_clears(mpfr, ...)
*/

  fun prec_round(mpfr: SMPFr, prec: ILong, rnd: RoundingMode): I32 =>
    @mpfr_prec_round(mpfr, prec, rnd())

  fun can_round(mpfr: SMPFr, err: I64, rnd1: RoundingMode, rnd2: RoundingMode, prec: ILong): Bool =>
    @mpfr_can_round(mpfr, err, rnd1(), rnd2(), prec) != 0

  fun min_prec(mpfr: SMPFr): ILong =>
    @mpfr_min_prec(mpfr)

  fun get_exp(mpfr: SMPFr): I64 =>
    @mpfr_get_exp(mpfr)

  fun set_exp(mpfr: SMPFr, parg1: I64): I32 =>
    @mpfr_set_exp(mpfr, parg1)

  fun get_prec(mpfr: SMPFr): ILong =>
    @mpfr_get_prec(mpfr)

  fun set_prec(mpfr: SMPFr, prec: ILong) =>
    @mpfr_set_prec(mpfr, prec)

  fun set_prec_raw(mpfr: SMPFr, prec: ILong) =>
    @mpfr_set_prec_raw(mpfr, prec)

  fun set_default_prec(prec: ILong) =>
    @mpfr_set_default_prec(prec)

  fun get_default_prec(): ILong =>
    @mpfr_get_default_prec()

  fun set_d(mpfr: SMPFr, d: F64, rnd: RoundingMode): I32 =>
    @mpfr_set_d(mpfr, d, rnd())

  fun set_flt(mpfr: SMPFr, f: F32, rnd: RoundingMode): I32 =>
    @mpfr_set_flt(mpfr, f, rnd())

/*
  // No support for F128
  fun set_ld(mpfr: SMPFr, f: F128, rnd: RoundingMode): I32 =>
    @mpfr_set_ld(mpfr, f, rnd())
*/

  fun set_z(mpfr: SMPFr, z: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_set_z(mpfr, z, rnd())

  fun set_z_2exp(mpfr: SMPFr, z: SMPFr, e: I64, rnd: RoundingMode): I32 =>
    @mpfr_set_z_2exp(mpfr, z, e, rnd())

  fun set_nan(mpfr: SMPFr) =>
    @mpfr_set_nan(mpfr)

  fun set_inf(mpfr: SMPFr, sign: I32) =>
    @mpfr_set_inf(mpfr, sign)

  fun set_zero(mpfr: SMPFr, sign: I32) =>
    @mpfr_set_zero(mpfr, sign)

  fun set_f(mpfr: SMPFr, f: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_set_f(mpfr, f, rnd())

  fun cmp_f(op1: SMPFr, op2: SMPFr): I32 =>
    @mpfr_cmp_f(op1, op2)

  fun get_f(mpfr: SMPFr, f: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_get_f(mpfr, f, rnd())

  fun set_si(mpfr: SMPFr, si: ILong, rnd: RoundingMode): I32 =>
    @mpfr_set_si(mpfr, si, rnd())

  fun set_ui(mpfr: SMPFr, ui: ULong, rnd: RoundingMode): I32 =>
    @mpfr_set_ui(mpfr, ui, rnd())

  fun set_sj(mpfr: SMPFr, sj: ISize, rnd: RoundingMode): I32 =>
    @__gmpfr_set_sj(mpfr, sj, rnd())

  fun set_uj(mpfr: SMPFr, uj: USize, rnd: RoundingMode): I32 =>
    @__gmpfr_set_uj(mpfr, uj, rnd())

  fun set_si_2exp(mpfr: SMPFr, si: ILong, e: I64, rnd: RoundingMode): I32 =>
    @mpfr_set_si_2exp(mpfr, si, e, rnd())

  fun set_ui_2exp(mpfr: SMPFr, ui: ULong, e: I64, rnd: RoundingMode): I32 =>
    @mpfr_set_ui_2exp(mpfr, ui, e, rnd())

  fun set_q(mpfr: SMPFr, q: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_set_q(mpfr, q, rnd())

  fun mul_q(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_mul_q(mpfr, op1, op2, rnd())

  fun div_q(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_div_q(mpfr, op1, op2, rnd())

  fun add_q(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_add_q(mpfr, op1, op2, rnd())

  fun sub_q(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_sub_q(mpfr, op1, op2, rnd())

  fun cmp_q(op1: SMPFr, op2: SMPFr): I32 =>
    @mpfr_cmp_q(op1, op2)

  fun get_q(q: SMPFr, f: SMPFr) =>
    @mpfr_get_q(q, f)

  fun set_str(mpfr: SMPFr, s: String, base: I32, rnd: RoundingMode): I32 =>
    @mpfr_set_str(mpfr, s.cstring(), base, rnd())

  fun init_set_str(mpfr: SMPFr, s: String, base: I32, rnd: RoundingMode): I32 =>
    @mpfr_init_set_str(mpfr, s.cstring(), base, rnd())

  fun set4(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode, t: I32): I32 =>
    @mpfr_set4(mpfr, op, rnd(), t)

  fun abs(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_abs(mpfr, op, rnd())

  fun set(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_set(mpfr, op, rnd())

  fun neg(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_neg(mpfr, op, rnd())

  fun signbit(mpfr: SMPFr): I32 =>
    @mpfr_signbit(mpfr)

  fun setsign(mpfr: SMPFr, op: SMPFr, s: I32, rnd: RoundingMode): I32 =>
    @mpfr_setsign(mpfr, op, s, rnd())

  fun copysign(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_copysign(mpfr, op1, op2, rnd())

  fun get_z_2exp(mpfr: SMPFr, op: SMPFr): I64 =>
    @mpfr_get_z_2exp(mpfr, op)

  fun get_flt(mpfr: SMPFr, rnd: RoundingMode): F32 =>
    @mpfr_get_flt(mpfr, rnd())

  fun get_d(mpfr: SMPFr, rnd: RoundingMode): F64 =>
    @mpfr_get_d(mpfr, rnd())

/*
  fun get_ld(mpfr: SMPFr, rnd: RoundingMode): F128 =>
    @mpfr_get_ld(mpfr, rnd())
*/

  fun get_d1(mpfr: SMPFr): F64 =>
    @mpfr_get_d1(mpfr)

  fun get_d_2exp(parg0: Pointer[I64] tag, parg1: SMPFr, rnd: RoundingMode): F64 =>
    @mpfr_get_d_2exp(parg0, parg1, rnd())

/*
  fun get_ld_2exp(parg0: Pointer[I64] tag, parg1: SMPFr, rnd: RoundingMode): F128 =>
    @mpfr_get_ld_2exp(parg0, parg1, rnd())
*/

  fun frexp(parg0: Pointer[I64] tag, parg1: SMPFr, parg2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_frexp(parg0, parg1, parg2, rnd())

  fun get_si(mpfr: SMPFr, rnd: RoundingMode): ILong =>
    @mpfr_get_si(mpfr, rnd())

  fun get_ui(mpfr: SMPFr, rnd: RoundingMode): ULong =>
    @mpfr_get_ui(mpfr, rnd())

  fun get_sj(mpfr: SMPFr, rnd: RoundingMode): ISize =>
    @__gmpfr_mpfr_get_sj(mpfr, rnd())

  fun get_uj(mpfr: SMPFr, rnd: RoundingMode): USize =>
    @__gmpfr_mpfr_get_uj(mpfr, rnd())

  fun get_str_ndigits(base: I32, prec: ILong): USize =>
    @mpfr_get_str_ndigits(base, prec)

  fun get_str(base: I32, prec: ILong, mpfr: SMPFr, rnd: RoundingMode): (String iso^, I64) =>
    let buf_size = MPF.get_str_ndigits(base, prec).max(7)
    let buf: Pointer[U8] = @pony_alloc(@pony_ctx(), buf_size)
    var exponent: I64 = 0
    var pcstring: Pointer[U8] =  @mpfr_get_str(buf, addressof exponent, base, buf_size, mpfr, rnd())
    let p: String iso = String.from_cstring(pcstring).clone()
    (consume p, exponent)

  fun get_z(z: SMPFr, f: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_get_z(z, f, rnd())

  fun free_str(parg0: Pointer[U8] tag) =>
    @mpfr_free_str(parg0)

  fun urandom(mpfr: SMPFr, parg1: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_urandom(mpfr, parg1, rnd())

  fun grandom(mpfr: SMPFr, parg1: SMPFr, parg2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_grandom(mpfr, parg1, parg2, rnd())

  fun nrandom(mpfr: SMPFr, parg1: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_nrandom(mpfr, parg1, rnd())

  fun erandom(mpfr: SMPFr, parg1: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_erandom(mpfr, parg1, rnd())

  fun urandomb(mpfr: SMPFr, parg1: SMPFr): I32 =>
    @mpfr_urandomb(mpfr, parg1)

  fun nextabove(mpfr: SMPFr) =>
    @mpfr_nextabove(mpfr)

  fun nextbelow(mpfr: SMPFr) =>
    @mpfr_nextbelow(mpfr)

  fun nexttoward(mpfr: SMPFr, parg1: SMPFr) =>
    @mpfr_nexttoward(mpfr, parg1)

/*
  // Pony does not support variadic functions
  fun printf(parg0: String, ...): I32 =>
    @mpfr_printf(parg0.cstring(), ...)

  fun asprintf(parg0: Array[String], parg1: String, ...): I32 =>
    @mpfr_asprintf(parg0, parg1.cstring(), ...)

  fun sprintf(parg0: String, parg1: String, ...): I32 =>
    @mpfr_sprintf(parg0.cstring(), parg1.cstring(), ...)

  fun snprintf(parg0: String, parg1: U64, parg2: String, ...): I32 =>
    @mpfr_snprintf(parg0.cstring(), parg1, parg2.cstring(), ...)
*/

  fun pow(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_pow(mpfr, op1, op2, rnd())

  fun pow_si(mpfr: SMPFr, op1: SMPFr, op2: ILong, rnd: RoundingMode): I32 =>
    @mpfr_pow_si(mpfr, op1, op2, rnd())

  fun pow_ui(mpfr: SMPFr, op1: SMPFr, op2: ULong, rnd: RoundingMode): I32 =>
    @mpfr_pow_ui(mpfr, op1, op2, rnd())

  fun ui_pow_ui(mpfr: SMPFr, op1: ULong, op2: ULong, rnd: RoundingMode): I32 =>
    @mpfr_ui_pow_ui(mpfr, op1, op2, rnd())

  fun ui_pow(mpfr: SMPFr, op1: ULong, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_ui_pow(mpfr, op1, op2, rnd())

  fun pow_z(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_pow_z(mpfr, op1, op2, rnd())

  fun sqrt(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_sqrt(mpfr, op, rnd())

  fun sqrt_ui(mpfr: SMPFr, op: ULong, rnd: RoundingMode): I32 =>
    @mpfr_sqrt_ui(mpfr, op, rnd())

  fun rec_sqrt(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_rec_sqrt(mpfr, op, rnd())

  fun add(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_add(mpfr, op1, op2, rnd())

  fun sub(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_sub(mpfr, op1, op2, rnd())

  fun mul(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_mul(mpfr, op1, op2, rnd())

  fun div(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_div(mpfr, op1, op2, rnd())

  fun add_ui(mpfr: SMPFr, op1: SMPFr, op2: ULong, rnd: RoundingMode): I32 =>
    @mpfr_add_ui(mpfr, op1, op2, rnd())

  fun sub_ui(mpfr: SMPFr, op1: SMPFr, op2: ULong, rnd: RoundingMode): I32 =>
    @mpfr_sub_ui(mpfr, op1, op2, rnd())

  fun ui_sub(mpfr: SMPFr, op1: ULong, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_ui_sub(mpfr, op1, op2, rnd())

  fun mul_ui(mpfr: SMPFr, op1: SMPFr, op2: ULong, rnd: RoundingMode): I32 =>
    @mpfr_mul_ui(mpfr, op1, op2, rnd())

  fun div_ui(mpfr: SMPFr, op1: SMPFr, op2: ULong, rnd: RoundingMode): I32 =>
    @mpfr_div_ui(mpfr, op1, op2, rnd())

  fun ui_div(mpfr: SMPFr, op1: ULong, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_ui_div(mpfr, op1, op2, rnd())

  fun add_si(mpfr: SMPFr, op1: SMPFr, op2: ILong, rnd: RoundingMode): I32 =>
    @mpfr_add_si(mpfr, op1, op2, rnd())

  fun sub_si(mpfr: SMPFr, op1: SMPFr, op2: ILong, rnd: RoundingMode): I32 =>
    @mpfr_sub_si(mpfr, op1, op2, rnd())

  fun si_sub(mpfr: SMPFr, op1: ILong, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_si_sub(mpfr, op1, op2, rnd())

  fun mul_si(mpfr: SMPFr, op1: SMPFr, op2: ILong, rnd: RoundingMode): I32 =>
    @mpfr_mul_si(mpfr, op1, op2, rnd())

  fun div_si(mpfr: SMPFr, op1: SMPFr, op2: ILong, rnd: RoundingMode): I32 =>
    @mpfr_div_si(mpfr, op1, op2, rnd())

  fun si_div(mpfr: SMPFr, op1: ILong, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_si_div(mpfr, op1, op2, rnd())

  fun add_d(mpfr: SMPFr, op1: SMPFr, op2: F64, rnd: RoundingMode): I32 =>
    @mpfr_add_d(mpfr, op1, op2, rnd())

  fun sub_d(mpfr: SMPFr, op1: SMPFr, op2: F64, rnd: RoundingMode): I32 =>
    @mpfr_sub_d(mpfr, op1, op2, rnd())

  fun d_sub(mpfr: SMPFr, op1: F64, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_d_sub(mpfr, op1, op2, rnd())

  fun mul_d(mpfr: SMPFr, op1: SMPFr, op2: F64, rnd: RoundingMode): I32 =>
    @mpfr_mul_d(mpfr, op1, op2, rnd())

  fun div_d(mpfr: SMPFr, op1: SMPFr, op2: F64, rnd: RoundingMode): I32 =>
    @mpfr_div_d(mpfr, op1, op2, rnd())

  fun d_div(mpfr: SMPFr, op1: F64, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_d_div(mpfr, op1, op2, rnd())

  fun sqr(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_sqr(mpfr, op, rnd())

  fun const_pi(mpfr: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_const_pi(mpfr, rnd())

  fun const_log2(mpfr: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_const_log2(mpfr, rnd())

  fun const_euler(mpfr: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_const_euler(mpfr, rnd())

  fun const_catalan(mpfr: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_const_catalan(mpfr, rnd())

  fun agm(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_agm(mpfr, op1, op2, rnd())

  fun log(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_log(mpfr, op, rnd())

  fun log2(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_log2(mpfr, op, rnd())

  fun log10(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_log10(mpfr, op, rnd())

  fun log1p(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_log1p(mpfr, op, rnd())

  fun log_ui(mpfr: SMPFr, op: ULong, rnd: RoundingMode): I32 =>
    @mpfr_log_ui(mpfr, op, rnd())

  fun exp(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_exp(mpfr, op, rnd())

  fun exp2(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_exp2(mpfr, op, rnd())

  fun exp10(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_exp10(mpfr, op, rnd())

  fun expm1(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_expm1(mpfr, op, rnd())

  fun eint(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_eint(mpfr, op, rnd())

  fun li2(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_li2(mpfr, op, rnd())

  fun cmp(op1: SMPFr, op2: SMPFr): I32 =>
    @mpfr_cmp(op1, op2)

  fun cmp3(op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_cmp3(op1, op2, rnd())

  fun cmp_d(op1: SMPFr, op2: F64): I32 =>
    @mpfr_cmp_d(op1, op2)

/*
  // Not support for F128
  fun cmp_ld(op1: SMPFr, op2: F128): I32 =>
    @mpfr_cmp_ld(op1, op2)
*/

  fun cmp_ui(op1: SMPFr, op2: ULong): I32 =>
    @mpfr_cmp_ui(op1, op2)

  fun cmp_si(op1: SMPFr, op2: ILong): I32 =>
    @mpfr_cmp_si(op1, op2)

  fun cmp_ui_2exp(op1: SMPFr, op2: ULong, e: I64): I32 =>
    @mpfr_cmp_ui_2exp(op1, op2, e)

  fun cmp_si_2exp(op1: SMPFr, op2: ILong, e: I64): I32 =>
    @mpfr_cmp_si_2exp(op1, op2, e)

  fun cmpabs(op1: SMPFr, op2: SMPFr): I32 =>
    @mpfr_cmpabs(op1, op2)

  fun cmpabs_ui(op1: SMPFr, op2: ULong): I32 =>
    @mpfr_cmpabs_ui(op1, op2)

  fun reldiff(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode) =>
    @mpfr_reldiff(mpfr, op1, op2, rnd())

  fun eq(op1: SMPFr, op2: SMPFr, op3: ULong): I32 =>
    @mpfr_eq(op1, op2, op3)

  fun sgn(mpfr: SMPFr): I32 =>
    @mpfr_sgn(mpfr)

  fun mul_2exp(mpfr: SMPFr, op1: SMPFr, op2: U64, rnd: RoundingMode): I32 =>
    @mpfr_mul_2exp(mpfr, op1, op2, rnd())

  fun div_2exp(mpfr: SMPFr, op1: SMPFr, op2: U64, rnd: RoundingMode): I32 =>
    @mpfr_div_2exp(mpfr, op1, op2, rnd())

  fun mul_2ui(mpfr: SMPFr, op1: SMPFr, op2: ULong, rnd: RoundingMode): I32 =>
    @mpfr_mul_2ui(mpfr, op1, op2, rnd())

  fun div_2ui(mpfr: SMPFr, op1: SMPFr, op2: ULong, rnd: RoundingMode): I32 =>
    @mpfr_div_2ui(mpfr, op1, op2, rnd())

  fun mul_2si(mpfr: SMPFr, op1: SMPFr, op2: ILong, rnd: RoundingMode): I32 =>
    @mpfr_mul_2si(mpfr, op1, op2, rnd())

  fun div_2si(mpfr: SMPFr, op1: SMPFr, op2: ILong, rnd: RoundingMode): I32 =>
    @mpfr_div_2si(mpfr, op1, op2, rnd())

  fun rint(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_rint(mpfr, op, rnd())

  fun roundeven(mpfr: SMPFr, op: SMPFr): I32 =>
    @mpfr_roundeven(mpfr, op)

  fun round(mpfr: SMPFr, op: SMPFr): I32 =>
    @mpfr_round(mpfr, op)

  fun trunc(mpfr: SMPFr, op: SMPFr): I32 =>
    @mpfr_trunc(mpfr, op)

  fun ceil(mpfr: SMPFr, op: SMPFr): I32 =>
    @mpfr_ceil(mpfr, op)

  fun floor(mpfr: SMPFr, op: SMPFr): I32 =>
    @mpfr_floor(mpfr, op)

  fun rint_roundeven(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_rint_roundeven(mpfr, op, rnd())

  fun rint_round(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_rint_round(mpfr, op, rnd())

  fun rint_trunc(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_rint_trunc(mpfr, op, rnd())

  fun rint_ceil(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_rint_ceil(mpfr, op, rnd())

  fun rint_floor(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_rint_floor(mpfr, op, rnd())

  fun frac(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_frac(mpfr, op, rnd())

  fun modf(mpfr1: SMPFr, mpfr2: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_modf(mpfr1, mpfr2, op, rnd())

  fun remquo(mpfr: SMPFr, q: Pointer[I64] tag, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_remquo(mpfr, q, op1, op2, rnd())

  fun remainder(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_remainder(mpfr, op1, op2, rnd())

  fun fmod(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_fmod(mpfr, op1, op2, rnd())

  fun fmodquo(mpfr: SMPFr, q: Pointer[I64] tag, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_fmodquo(mpfr, q, op1, op2, rnd())

  fun fits_ulong_p(mpfr: SMPFr, rnd: RoundingMode): Bool =>
    @mpfr_fits_ulong_p(mpfr, rnd()) != 0

  fun fits_slong_p(mpfr: SMPFr, rnd: RoundingMode): Bool =>
    @mpfr_fits_slong_p(mpfr, rnd()) != 0

  fun fits_uint_p(mpfr: SMPFr, rnd: RoundingMode): Bool =>
    @mpfr_fits_uint_p(mpfr, rnd()) != 0

  fun fits_sint_p(mpfr: SMPFr, rnd: RoundingMode): Bool =>
    @mpfr_fits_sint_p(mpfr, rnd()) != 0

  fun fits_ushort_p(mpfr: SMPFr, rnd: RoundingMode): Bool =>
    @mpfr_fits_ushort_p(mpfr, rnd()) != 0

  fun fits_sshort_p(mpfr: SMPFr, rnd: RoundingMode): Bool =>
    @mpfr_fits_sshort_p(mpfr, rnd()) != 0

  fun fits_uintmax_p(mpfr: SMPFr, rnd: RoundingMode): Bool =>
    @mpfr_fits_uintmax_p(mpfr, rnd()) != 0

  fun fits_intmax_p(mpfr: SMPFr, rnd: RoundingMode): Bool =>
    @mpfr_fits_intmax_p(mpfr, rnd()) != 0

  fun extract(mpfr: SMPFr, parg1: SMPFr, parg2: U32) =>
    @mpfr_extract(mpfr, parg1, parg2)

  fun swap(op1: SMPFr, op2: SMPFr) =>
    @mpfr_swap(op1, op2)

  fun dump(mpfr: SMPFr) =>
    @mpfr_dump(mpfr)

  fun nan_p(mpfr: SMPFr): Bool =>
    @mpfr_nan_p(mpfr) != 0

  fun inf_p(mpfr: SMPFr): Bool =>
    @mpfr_inf_p(mpfr) != 0

  fun number_p(mpfr: SMPFr): Bool =>
    @mpfr_number_p(mpfr) != 0

  fun integer_p(mpfr: SMPFr): Bool =>
    @mpfr_integer_p(mpfr) != 0

  fun zero_p(mpfr: SMPFr): Bool =>
    @mpfr_zero_p(mpfr) != 0

  fun regular_p(mpfr: SMPFr): Bool =>
    @mpfr_regular_p(mpfr) != 0

  fun greater_p(mpfr: SMPFr, parg1: SMPFr): Bool =>
    @mpfr_greater_p(mpfr, parg1) != 0

  fun greaterequal_p(mpfr: SMPFr, parg1: SMPFr): Bool =>
    @mpfr_greaterequal_p(mpfr, parg1) != 0

  fun less_p(mpfr: SMPFr, parg1: SMPFr): Bool =>
    @mpfr_less_p(mpfr, parg1) != 0

  fun lessequal_p(mpfr: SMPFr, parg1: SMPFr): Bool =>
    @mpfr_lessequal_p(mpfr, parg1) != 0

  fun lessgreater_p(mpfr: SMPFr, parg1: SMPFr): Bool =>
    @mpfr_lessgreater_p(mpfr, parg1) != 0

  fun equal_p(mpfr: SMPFr, parg1: SMPFr): Bool =>
    @mpfr_equal_p(mpfr, parg1) != 0

  fun unordered_p(mpfr: SMPFr, parg1: SMPFr): Bool =>
    @mpfr_unordered_p(mpfr, parg1) != 0

  fun atanh(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_atanh(mpfr, op, rnd())

  fun acosh(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_acosh(mpfr, op, rnd())

  fun asinh(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_asinh(mpfr, op, rnd())

  fun cosh(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_cosh(mpfr, op, rnd())

  fun sinh(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_sinh(mpfr, op, rnd())

  fun tanh(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_tanh(mpfr, op, rnd())

  fun sinh_cosh(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_sinh_cosh(mpfr, op1, op2, rnd())

  fun sech(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_sech(mpfr, op, rnd())

  fun csch(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_csch(mpfr, op, rnd())

  fun coth(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_coth(mpfr, op, rnd())

  fun acos(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_acos(mpfr, op, rnd())

  fun asin(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_asin(mpfr, op, rnd())

  fun atan(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_atan(mpfr, op, rnd())

  fun sin(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_sin(mpfr, op, rnd())

  fun sin_cos(mpfr1: SMPFr, mpfr2: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_sin_cos(mpfr1, mpfr2, op, rnd())

  fun cos(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_cos(mpfr, op, rnd())

  fun tan(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_tan(mpfr, op, rnd())

  fun atan2(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_atan2(mpfr, op1, op2, rnd())

  fun sec(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_sec(mpfr, op, rnd())

  fun csc(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_csc(mpfr, op, rnd())

  fun cot(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_cot(mpfr, op, rnd())

  fun hypot(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_hypot(mpfr, op1, op2, rnd())

  fun erf(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_erf(mpfr, op, rnd())

  fun erfc(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_erfc(mpfr, op, rnd())

  fun cbrt(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_cbrt(mpfr, op, rnd())

  fun root(mpfr: SMPFr, op1: SMPFr, op2: ULong, rnd: RoundingMode): I32 =>
    @mpfr_root(mpfr, op1, op2, rnd())

  fun rootn_ui(mpfr: SMPFr, op1: SMPFr, op2: ULong, rnd: RoundingMode): I32 =>
    @mpfr_rootn_ui(mpfr, op1, op2, rnd())

  fun gamma(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_gamma(mpfr, op, rnd())

  fun gamma_inc(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_gamma_inc(mpfr, op1, op2, rnd())

  fun beta(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_beta(mpfr, op1, op2, rnd())

  fun lngamma(mpfr: SMPFr, op1: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_lngamma(mpfr, op1, rnd())

  fun lgamma(mpfr: SMPFr, op1: Pointer[I32] tag, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_lgamma(mpfr, op1, op2, rnd())

  fun digamma(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_digamma(mpfr, op, rnd())

  fun zeta(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_zeta(mpfr, op, rnd())

  fun zeta_ui(mpfr: SMPFr, op: ULong, rnd: RoundingMode): I32 =>
    @mpfr_zeta_ui(mpfr, op, rnd())

  fun fac_ui(mpfr: SMPFr, op: ULong, rnd: RoundingMode): I32 =>
    @mpfr_fac_ui(mpfr, op, rnd())

  fun j0(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_j0(mpfr, op, rnd())

  fun j1(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_j1(mpfr, op, rnd())

  fun jn(mpfr: SMPFr, op1: ILong, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_jn(mpfr, op1, op2, rnd())

  fun y0(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_y0(mpfr, op, rnd())

  fun y1(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_y1(mpfr, op, rnd())

  fun yn(mpfr: SMPFr, op1: ILong, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_yn(mpfr, op1, op2, rnd())

  fun ai(mpfr: SMPFr, op: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_ai(mpfr, op, rnd())

  fun min(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_min(mpfr, op1, op2, rnd())

  fun max(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_max(mpfr, op1, op2, rnd())

  fun dim(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_dim(mpfr, op1, op2, rnd())

  fun mul_z(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_mul_z(mpfr, op1, op2, rnd())

  fun div_z(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_div_z(mpfr, op1, op2, rnd())

  fun add_z(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_add_z(mpfr, op1, op2, rnd())

  fun sub_z(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_sub_z(mpfr, op1, op2, rnd())

  fun z_sub(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_z_sub(mpfr, op1, op2, rnd())

  fun cmp_z(op1: SMPFr, op2: SMPFr): I32 =>
    @mpfr_cmp_z(op1, op2)

  fun fma(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, op3: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_fma(mpfr, op1, op2, op3, rnd())

  fun fms(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, op3: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_fms(mpfr, op1, op2, op3, rnd())

  fun fmma(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, op3: SMPFr, op4: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_fmma(mpfr, op1, op2, op3, op4, rnd())

  fun fmms(mpfr: SMPFr, op1: SMPFr, op2: SMPFr, op3: SMPFr, op4: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_fmms(mpfr, op1, op2, op3, op4, rnd())

  // TODO: Interface will change
  //fun sum(mpfr: SMPFr, tab: Pointer[SMPFr] tag, n: USize, rnd: RoundingMode): I32 =>
  fun sum(mpfr: SMPFr, tab: Array[SMPFr], rnd: RoundingMode): I32 =>
    let ptab = tab.cpointer()
    let n = tab.size()
    @mpfr_sum(mpfr, ptab, n, rnd())

  // TODO: Interface will change
  //fun dot(mpfr: SMPFr, tab1: Pointer[SMPFr] tag, tab2: Pointer[SMPFr] tag, n: USize, rnd: RoundingMode): I32 =>
  fun dot(mpfr: SMPFr, tab1: Array[SMPFr], tab2: Array[SMPFr], rnd: RoundingMode): I32 =>
    let ptab1 = tab1.cpointer()
    let ptab2 = tab2.cpointer()
    let n1 = tab1.size()
    let n2 = tab2.size()
    ifdef debug then
      (n1 == n2) or
        Fail(Format("[MPF.dot] The two arrays don't have the same size: " +
              "{} != {}. Can't calculate dot product.", [n1; n2]))
    end
    @mpfr_dot(mpfr, ptab1, ptab2, n1, rnd())

  fun free_cache() =>
    @mpfr_free_cache()

  fun free_cache2(way: I32) =>
    @mpfr_free_cache2(way)

  fun free_pool() =>
    @mpfr_free_pool()

  fun mp_memory_cleanup(): I32 =>
    @mpfr_mp_memory_cleanup()

  fun subnormalize(mpfr: SMPFr, t: I32, rnd: RoundingMode): I32 =>
    @mpfr_subnormalize(mpfr, t, rnd())

/*
  fun strtofr(mpfr: SMPFr, parg1: String, parg2: Array[String], parg3: I32, parg4: RoundingMode): I32 =>
    @mpfr_strtofr(mpfr, parg1.cstring(), parg2, parg3, parg4)
*/

  fun round_nearest_away_begin(mpfr: SMPFr) =>
    @mpfr_round_nearest_away_begin(mpfr)

  fun round_nearest_away_end(mpfr: SMPFr, rnd: RoundingMode): I32 =>
    @mpfr_round_nearest_away_end(mpfr, rnd())

  fun custom_get_size(prec: ILong): USize =>
    @mpfr_custom_get_size(prec)

  fun custom_init(significand: Pointer[None] tag, prec: ILong) =>
    @mpfr_custom_init(significand, prec)

  fun custom_get_significand(mpfr: SMPFr): Pointer[None] =>
    @mpfr_custom_get_significand(mpfr)

  fun custom_get_exp(mpfr: SMPFr): I64 =>
    @mpfr_custom_get_exp(mpfr)

  fun custom_move(mpfr: SMPFr, new_position: Pointer[None] tag) =>
    @mpfr_custom_move(mpfr, new_position)

  fun custom_init_set(mpfr: SMPFr, kind: I32, e: I64, prec: ILong, significand: Pointer[None] tag) =>
    @mpfr_custom_init_set(mpfr, kind, e, prec, significand)

  fun custom_get_kind(mpfr: SMPFr): I32 =>
    @mpfr_custom_get_kind(mpfr)

  fun total_order_p(op1: SMPFr, op2: SMPFr): I32 =>
    @mpfr_total_order_p(op1, op2)
