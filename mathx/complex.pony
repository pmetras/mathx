// Arithmetic of complex numbers

use "debug"
use "collections"

use "../assertx"
use "../pony_testx"
use "../formatx"



class val Complex[F: (Float & FloatingPoint[F]) = F64]
  is (Equatable[Complex[F]] & Stringable & Formattable & Approximated[Complex[F], F])
  """
  Complex numbers done right for numeric calculations. These functions try to
  avoid some overflows, underflows or loss of precision.

  Default precision for `Complex` is `F64`.

  Most of all operations are also available using unsafe arithmetic operators.
  """
  
  let _re: F
    """
    Real part
    """


  let _im: F
    """
    Imaginary part
    """

  
  new val create(re: F = F.from[ISize](0), im: F = F.from[ISize](0)) =>
    """
    Create a new complex number from real `re` and imaginary `im`parts. If the
    real or imaginary part are not specified, they default to 0.
    Create `re + i * im`.
    """
    _re = re
    _im = im


  new val j() =>
    """
    The `i` complex in literature that is known for `i * i == -1`.

    Name *j* was selected instead of *i* for the reasons:

    - `i` is frequently used as loop variable and to prevent name clash
    - Better visual distinction from other characters like 1 in source
    - Follows some physics conventions using j instead of i.
    """
    _re = F.from[ISize](0)
    _im = F.from[ISize](1)


  new val from_polar(r: F, theta: F) =>
    """
    Create a complex number from polar form `r * exp(i * theta)`.
    `re = r * cos(theta)`, `im = r * sin(theta)`.

    The inverse is `z.abs()` for `r` and `z.arg()` for `theta`.
    """
    _re = r * theta.cos()
    _im = r * theta.sin()


  new val from_string(s: String) ? =>
    """
    Create a complex number interpreting the string `s`.

    **BNF Grammar**

    ```BNF
    complex    ::= ws (pair | polar | cartesian) ws
    pair       ::= '(' ws number ws ',' ws number ws ')'
    polar      ::= number ws '∠' ws number
    cartesian  ::= term (ws term)?
                | term          /* pure real or pure imag */
    term       ::= real_term | imag_term
    real_term  ::= sign? number
    imag_term  ::= sign? imag_coeff
    imag_coeff ::= unit_i                 /* i, j        */
                | unit_i number           /* i3.5, j42   */
                | number unit_i           /* 3.5i, 42j   */
    unit_i     ::= 'i' | 'j'
    sign       ::= '+' | '-'
    number     ::= 'inf' | '@Inf@' | 'nan' | 'NaN' | '@NaN@'
                | digits ('.' digits?)? exp?
                | '.' digits exp?
    digits     ::= [0-9]+
    exp        ::= ('e' | 'E') sign? digits
    ws         ::= ' '*
    ```

    The grammar does not validate that `s` is a valid complex as is accepts
    invalid input like "4i+8i" that is not interpreted as "12i".

    Disambiguation rules:

    - `i` prefix vs `inf`: `i` followed by a digit is imag_coeff (i3.5), but
      `i` followed by `n` is the start of `inf`. The rule is: `unit_i` only
      matches a bare `i` or `j` not followed by an alphanumeric character.
    - Sign of second term: in `a-bi`, the `-` belongs to the second term as
      its sign, not to `a`. The parser must treat `+/-` after a complete
      first term as the start of a new term, not as part of number.
    - Two imaginary terms: `3i+5i` — two `imag_term`s. It is an error
      (only one real, one imaginary allowed).
    - Two real terms: `3+5` — two real_terms with no `i` nor `j`. It is an
      error too.
    - Sign before `i` with `inf` lookahead: `-inf` — sign is `-`, then `inf`
      is a number in a real_term. The `i` of `inf` must be consumed as part
      of the word, not as `unit_i`. Must attempt to match number (including
      `inf` or `nan`) before attempting `unit_i`.

    The following examples are accepted:

    ```pony
    "1.23"                          // Real
    "i"                             // Pure imaginary
    "+i"                            // Positive imaginary
    "-i"                            // Negative imaginary
    "-j"                            // Using j instead of i
    "5+6.78i"                       // i after
    "-6+i854"                       // i before
    "-9.9e12-7.32e5i"               // Scientific numbers
    "inf-5.3j"                      // Infinity
    "1+iinf"                        // Infinity imaginary
    "2-infi"                        // Infinity imaginary
    "infj"                          // Idem but with j
    "nan"                           // Not a number
    "jNaN"                          // Not a number
    "4.2i+7.0"                      // Reverse
    "1.00E+00∠7.85E-01"             // Using polar coordinates
    " -1.23e+03  ∠    7.89e-01   "  // With spaces
    "@NaN@∠3.14"                    // Not a number
    "(5.6, 7.8)"                    // As a pair of numbers
    "(-inf, NaN)"                   // Not a number
    ```
    """
    // Skip leading and trailing whitespace; work on a trimmed val view.
    let t: String val = recover s.clone().>strip() end

    // ── (re, im) pair ─────────────────────────────────────────────────────
    if (t.size() >= 2) and (t.at_offset(0)? == '(') then
      let inner = t.substring(1, (t.size() - 1).isize())
      // find the comma
      var comma: ISize = -1
      for ci in Range(0, inner.size()) do
        if inner.at_offset(ci.isize())? == ',' then
          comma = ci.isize()
          break
        end
      end
      if comma < 0 then
        Debug(Format("[Complex.from_string] Cannot find comma when reading (real, imag) pair from \"{}\"", t))
        error
      end

      let re_str: String val = recover inner.substring(0, comma).>strip() end
      let im_str: String val = recover inner.substring(comma + 1).>strip() end
      if (re_str.size() == 0) or (im_str.size() == 0) then
        Debug(Format("[Complex.from_string] Missing value in ({}, {}) pair", [re_str; im_str]))
        error
      end

      (let rv, _) = _parse_f64(re_str, 0)?
      (let iv, _) = _parse_f64(im_str, 0)?
      _re = F.from[F64](rv)
      _im = F.from[F64](iv)
      return
    end

    // ── polar  rho ∠ theta ────────────────────────────────────────────────
    // ∠ = UTF-8 E2 88 A0 (3 bytes)
    let angle_bytes: Array[U8] val = [0xe2; 0x88; 0xa0]
    var angle_pos: ISize = -1
    var bi: USize = 0
    while bi < t.size() do
      if (bi + 2) < t.size() then
        try
          if (t.at_offset(bi.isize())? == 0xe2)
            and (t.at_offset((bi + 1).isize())? == 0x88)
            and (t.at_offset((bi + 2).isize())? == 0xa0)
          then
            angle_pos = bi.isize()
            break
          end
        end
      end
      bi = bi + 1
    end

    if angle_pos >= 0 then
      let rho_str:   String val = recover t.substring(0, angle_pos).>strip() end
      let theta_str: String val = recover t.substring(angle_pos + 3).>strip() end
      (let rv, _) = _parse_f64(rho_str, 0)?
      (let tv, _) = _parse_f64(theta_str, 0)?
      let rho   = F.from[F64](rv)
      let theta = F.from[F64](tv)
      _re = rho * theta.cos()
      _im = rho * theta.sin()
      return
    end

    // ── cartesian: one or two numeric tokens separated by +/- ────────────
    // An imaginary token ends with i or j (possibly preceded by a number).
    // A real token is a bare number.
    // We parse up to two tokens. Each token is:
    //   [sign] [number] [i/j]
    // where [number] may be absent (meaning 1) and [i/j] marks imaginary.
    //
    // Scan for tokens by consuming the string left-to-right.
    var re_val: F64 = 0.0
    var im_val: F64 = 0.0
    var found_re: Bool = false
    var found_im: Bool = false
    var pos: USize = 0

    var tok: USize = 0
    while (tok < 2) and (pos < t.size()) do
      // Skip spaces.
      while (pos < t.size()) and (t.at_offset(pos.isize())? == ' ') do
        pos = pos + 1
      end
      if pos >= t.size() then
        break
      end

      // Collect optional sign.
      let sign_ch: U8 = t.at_offset(pos.isize())?
      let explicit_sign: Bool = (sign_ch == '+') or (sign_ch == '-')
      let sign: F64 = if sign_ch == '-' then -1.0 else 1.0 end
      if explicit_sign then
        pos = pos + 1
      end

      // Skip spaces after sign.
      while (pos < t.size()) and (t.at_offset(pos.isize())? == ' ') do
        pos = pos + 1
      end
      if pos >= t.size() then
        Debug(Format("[Complex.from_string] Bare sign with nothing after in \"{}\"", t))
        error
      end

      // Check for i/j prefix (e.g. "i854", "-i854", "iinf") or unit i/j ("+i", "-j").
      // Disambiguation: unit_i only when NOT followed by alphanumeric or number-start.
      let next_ch: U8 = t.at_offset(pos.isize())?
      if (next_ch == 'i') or (next_ch == 'j') then
        let after_pos: USize = pos + 1
        let ch2: U8 =
          if after_pos < t.size() then
            t.at_offset(after_pos.isize())?
          else
            0
          end
        let ch2_alpha = ((ch2 >= 'a') and (ch2 <= 'z'))
          or ((ch2 >= 'A') and (ch2 <= 'Z'))
        let ch2_digit = (ch2 >= '0') and (ch2 <= '9')
        let ch2_num_start = (ch2 == '.') or (ch2 == '+') or (ch2 == '-')
        if ch2_digit or ch2_num_start then
          // i/j is an imaginary-unit prefix before digits: "i854", "i.5".
          pos = after_pos
          (let v, let new_pos) = _parse_f64(t, pos)?
          pos = new_pos
          if found_im then
            Debug(Format("[Complex.from_string] Complex \"{}\" can't have more than 1 imaginary part", t))
            error
          end
          im_val = sign * v
          found_im = true
        elseif ch2_alpha then
          // Disambiguate: "inf" (i is part of the word, real term) vs
          // "iinf"/"iNaN" (i is imaginary prefix) vs "jNaN" (j is prefix).
          // Rule: attempt number (including inf/nan) before unit_i.
          // Try parse_f64 from pos first; if that succeeds, it's a real term
          // (e.g. "inf" → +inf). If it fails, i/j is a prefix for imaginary.
          var tried_real = false
          var real_v: F64 = 0.0
          var real_pos: USize = 0
          try
            (let rv, let rp) = _parse_f64(t, pos)?
            real_v = rv
            real_pos = rp
            tried_real = true
          end
          if tried_real then
            // Check trailing i/j (handles "infj" → imaginary inf).
            var np = real_pos
            while (np < t.size()) and (t.at_offset(np.isize())? == ' ') do
              np = np + 1
            end
            let trailing_ij: Bool = if np < t.size() then
                let tc2 = t.at_offset(np.isize())?
                if (tc2 == 'i') or (tc2 == 'j') then
                  np = np + 1
                  true
                else
                  false
                end
              else
                false
              end
            pos = np
            if trailing_ij then
              if found_im then
                Debug(Format("[Complex.from_string] Complex \"{}\" can't have more than 1 imaginary part", t))
                error
              end
              im_val = sign * real_v
              found_im = true
            else
              if found_re then
                Debug(Format("[Complex.from_string] Complex \"{}\" can't have more than 1 real part", t))
                error
              end
              re_val = sign * real_v
              found_re = true
            end
          else
            // parse_f64 from pos failed → i/j is imaginary prefix.
            pos = after_pos
            (let v, let new_pos) = _parse_f64(t, pos)?
            pos = new_pos
            if found_im then
              Debug(Format("[Complex.from_string] Complex \"{}\" can't have more than 1 imaginary part", t))
              error
            end
            im_val = sign * v
            found_im = true
          end
        else
          // Unit imaginary: ±i / ±j (next char is space, sign, end, etc.)
          pos = after_pos
          if found_im then
            Debug(Format("[Complex.from_string] Complex \"{}\" can't have more than 1 imaginary part", t))
            error
          end
          im_val = sign
          found_im = true
        end
      else
        // Try to parse a number.
        (let v, let new_pos) = _parse_f64(t, pos)?
        pos = new_pos
        // Skip spaces after number.
        while (pos < t.size()) and (t.at_offset(pos.isize())? == ' ') do
          pos = pos + 1
        end
        // Check for trailing i/j (e.g. "6.78i").
        let has_ij: Bool = if pos < t.size() then
            let ch = t.at_offset(pos.isize())?
            if (ch == 'i') or (ch == 'j') then
              pos = pos + 1
              true
            else
              false
            end
          else
            false
          end
        if has_ij then
          if found_im then
            Debug(Format("[Complex.from_string] Complex \"{}\" can't have more than 1 imaginary part", t))
            error
          end
          im_val = sign * v
          found_im = true
        else
          if found_re then
            Debug(Format("[Complex.from_string] Complex \"{}\" can't have more than 1 real part", t))
            error
          end
          re_val = sign * v
          found_re = true
        end
      end
      tok = tok + 1
    end

    if pos < t.size() then
      // Skip trailing spaces.
      while (pos < t.size()) and (t.at_offset(pos.isize())? == ' ') do
        pos = pos + 1
      end
      if pos < t.size() then
        let c = t.at_offset(pos.isize())?
        Debug(Format("[Complex.from_string] Unrecognized character '{}' while parsing complex \"{}\"", [c; t]))
        error
      end
    end

    if (not found_re) and (not found_im) then error end
    _re = F.from[F64](re_val)
    _im = F.from[F64](im_val)


  fun tag _parse_f64(src: String val, pos: USize): (F64, USize) ? =>
    """
    Scan one numeric token from `src` starting at byte `pos` (leading spaces
    skipped). Returns `(value, new_pos)` or errors if no number is found.
    Handles: regular decimals with optional exponent, `inf`/`nan` (3-letter
    case-insensitive), and the `@Inf@` / `@NaN@` sentinel forms.
    """
    var p: USize = pos
    // Skip spaces.
    while (p < src.size()) and (src(p)? == ' ') do
      p = p + 1
    end
    if p >= src.size() then
      Debug(Format("[Complex._parse_f64] Can't parse number in \"{}\"", src))
      error
    end
    let start: USize = p
    // Optional sign.
    let ch0 = src(p)?
    if (ch0 == '+') or (ch0 == '-') then
      p = p + 1
    end
    if p >= src.size() then
      Debug(Format("[Complex._parse_f64] Can't parse number after sign in \"{}\"", src))
      error
    end
    // Check for special forms: @Inf@ / @NaN@ or inf / nan.
    let ch1_raw = src(p)?
    let ch1 = ch1_raw or 0x20  // lowercase
    if ch1_raw == '@' then
      // @Inf@ or @NaN@: consume up to closing @.
      p = p + 1  // skip opening @
      while (p < src.size()) and (src(p)? != '@') do
        p = p + 1
      end
      if p < src.size() then // skip closing @
        p = p + 1
      end
    elseif (ch1 == 'i') or (ch1 == 'n') then
      // Consume exactly 3 letter chars (covers "inf"/"nan"; stops before trailing i/j suffix).
      var word_count: USize = 0
      while (p < src.size()) and (word_count < 3) and
            (let c = src(p)? or 0x20
              ((c >= 'a') and (c <= 'z'))) do
        p = p + 1
        word_count = word_count + 1
      end
    else
      // Digits and decimal point.
      while (p < src.size()) and
            (let c = src(p)?
              ((c >= '0') and (c <= '9')) or (c == '.')) do
        p = p + 1
      end
      // Optional exponent: e/E followed by optional sign and digits.
      if p < src.size() then
        let ec = src(p)? or 0x20
        if ec == 'e' then
          p = p + 1
          if p < src.size() then
            let sc = src(p)?
            if (sc == '+') or (sc == '-') then
              p = p + 1
            end
          end
          while (p < src.size()) and
                (let c = src(p)?
                  (c >= '0') and (c <= '9')) do
            p = p + 1
          end
        end
      end
    end
    if p == start then
      Debug(Format("[Complex._parse_f64] Can't parse number from string \"{}\"", src))
      error
    end
    let raw: String val = recover src.substring(start.isize(), p.isize()) end
    // Normalize @Inf@ / @NaN@ to forms strtod recognises (strip @ delimiters).
    let token: String val = recover raw.clone().>replace("@", "") end
    (token.f64()?, p)


  fun string(): String iso^ =>
    """
    Convert the complex to a string for display with the cartesian form
    `re + im * i`. Real and imaginary parts are printed only when non-null.

    There's no special representation for a complex number with `inf` real
    or imaginary parts. When `NaN` real or imaginary, the string is `NaN`.
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    let result = match (_re, _im)
    | if (_re == zero) and (_im == zero) => "0"
    | if (_re == zero) and (_im == one) => "i"
    | if (_re == zero) and (_im == -one) => "-i"
    | if (_re == zero) => _im.string() + "i"
    | if (_im == zero) => _re.string()
    | if (_im == one) => _re.string() + "+i"
    | if (_im == -one) => _re.string() + "-i"
    | if (_im < zero) => _re.string() + _im.string() + "i"
    | if (_im > zero) => _re.string() + "+" + _im.string() + "i"
    else
      "nan"
    end
    result.clone()


  fun format(spec: String = ""): String =>
    """
    Format this complex number according to `spec`.

    The type character controls presentation:
    - `r` — polar form `(ρ∠θ)` where ρ = abs() and θ = arg() in radians.
    - `R` — polar form with uppercase exponent in scientific notation
             (e.g. `(1.00E+00∠7.85E-01)`).
    - anything else (including no type char) — cartesian form `(re+imi)`.

    The numeric sub-spec (precision, sign, `e`/`f`/`g`, width-modifiers) is
    forwarded to each component individually.  Width and alignment are applied
    to the assembled complex string as a whole.

    Examples:
    ```
    Complex[F64](1.0, -2.5).format(".3f")   // "(1.000-2.500i)"
    Complex[F64](1.0,  1.0).format(".4r")   // "(1.4142∠0.7854)"
    Complex[F64](0.0,  0.0).format(".2e")   // "(0.00e+00+0.00e+00i)"
    Complex[F64](1.0,  1.0).format("20.3f") // "     (1.000+1.000i)"
    ```
    """
    let fspec = FormatSpec(spec)

    // Build the inner numeric spec (no width — we apply width to the whole).
    // We strip the type char for component formatting: for polar r/R we use
    // 'e'/'E'; for cartesian we forward the type char as-is (e/f/g/default).
    let tc = fspec.type_char
    let polar = (tc == 'r') or (tc == 'R')
    let upper = tc == 'R'

    // Reconstruct a component spec: sign + optional # + optional z +
    // precision + component type char, but width = 0.
    let cspec: String ref = String
    match fspec.sign
    | SignPlus  => cspec.push('+')
    | SignSpace => cspec.push(' ')
    | SignMinus => cspec.push('-')
    | SignDefault => None
    end
    if fspec.z   then cspec.push('z') end
    if fspec.hash then cspec.push('#') end
    match fspec.precision
    | let p: USize => cspec.push('.'); cspec.append(p.string())
    end
    match fspec.grouping
    | let g: U8 => cspec.push(g)
    end
    if polar then
      cspec.push(if upper then 'E' else 'e' end)
    elseif (tc != 0) and not polar then
      cspec.push(tc)
    end
    let cs: String val = cspec.clone()

    let fmt_component = {(v: F): String => Format("{:" + cs + "}", v)}
    let inner: String =
      if polar then
        let rho   = abs()
        let theta = arg()
        "(" + fmt_component(rho) + "∠" + fmt_component(theta) + ")"
      else
        let re_str = fmt_component(_re)
        let im_abs = _im.abs()
        let im_neg = _im < F.from[ISize](0)
        let im_str = fmt_component(im_abs)
        let sign_str = if im_neg then "-" else "+" end
        "(" + re_str + sign_str + im_str + "i)"
      end

    // Apply width / fill / align to the whole string.
    if fspec.width == 0 then
      inner
    else
      let fill_char = fspec.fill
      let w = fspec.width
      let n = inner.codepoints()
      if n >= w then
        inner
      else
        let pad = w - n
        match fspec.align
        | AlignLeft =>
          inner + (String.from_utf32(fill_char) * pad)
        | AlignCenter =>
          let lpad = pad / 2
          let rpad = pad - lpad
          (String.from_utf32(fill_char) * lpad) + inner + (String.from_utf32(fill_char) * rpad)
        else
          // AlignRight / AlignDefault / AlignNumeric all right-align for complex.
          (String.from_utf32(fill_char) * pad) + inner
        end
      end
    end


  fun real(): F =>
    """
    Get the real part of a complex number.
    """
    _re


  fun imag(): F =>
    """
    Get the imaginary part of a complex number.
    """
    _im


  fun is_real(tol: F = F.epsilon()): Bool =>
    """
    `true` when the complex has no imaginary part, with a tolerance `tol`
    (default to `F.epsilon`).
    """
    ifdef debug then
      (tol >= F.from[ISize](0)) or
        Fail(Format("[Complex.is_real] Precision tolerance ({}) must be positive",
                    tol))
    end
    _im.abs() <= tol


  fun is_imag(abs_tol: F = F.epsilon()): Bool =>
    """
    `true` when the complex has no real part, with an absolute tolerance `abs_tol`
    (default to `F.epsilon)`.
    """
    ifdef debug then
      (abs_tol >= F.from[ISize](0)) or
        Fail(Format("[Complex.is_imag] Precision tolerance ({}) must be positive",
                    abs_tol))
    end
    _re.abs() <= abs_tol


  fun is_null(abs_tol: F = F.epsilon()): Bool =>
    """
    `true` when the complex is equal to 0.0 with an absolute tolerance `abs_tol`
    (default to `F.epsilon)`. The tolerance is applied individually on real and
    imaginary parts and not on the modulus.
    """
    ifdef debug then
      (abs_tol >= F.from[ISize](0)) or
        Fail(Format("[Complex.is_null] Precision tolerance ({}) must be positive",
                    abs_tol))
    end
    (_re.abs() <= abs_tol) and (_im.abs() <= abs_tol)


  fun finite(): Bool =>
    """
    Return `true` when real and imaginary parts of complex are finite.
    """
    _re.finite() and _im.finite()


  fun infinite(): Bool =>
    """
    Return `false` when real or imaginary part of complex is infinite.
    """
    _re.infinite() or _im.infinite()


  fun nan(): Bool =>
    """
    Return `true` when real or imaginary part of complex is NaN (Not a Number).
    """
    _re.nan() or _im.nan()


  fun conj(): Complex[F]^ =>
    """
    Create the conjugate of the complex, `re - i * im`.
    """
    Complex[F](_re, -_im)


  fun conj_unsafe(): Complex[F]^ =>
    """
    Create the conjugate of the complex, `re - i * im` using unsafe arithmetic.
    """
    Complex[F](_re, -~_im)


  fun eq(that: Complex[F]): Bool =>
    """
    Return `true` if `this` and `that` represent the same complex number. The
    two complex numbers are equal if
    `this.real() == that.real() and this.imag() == that.imag()`
    
    Because `Complex` can't be defined as a primitive, there can be multiple
    instances of the same complex number.
    """
    (_re == that._re) and (_im == that._im)


  fun eq_unsafe(that: Complex[F]): Bool =>
    """
    Return `true` if `this` and `that` represent the same complex number, using
    unsafe arithmetic. The two complex numbers are equal if
    `this.real() ==~ that.real() and this.imag() ==~ that.imag()`
    
    Because `Complex` can't be defined as a primitive, there can be multiple
    instances of the same complex number.
    """
    (_re ==~ that._re) and (_im ==~ that._im)


  fun almost_eq(that: box->Complex[F],
                rel_tol: F = F.epsilon().sqrt(),
                abs_tol: F = F.epsilon().sqrt())
               : Bool =>
    """
    Return `true` when `this` and `that` are almost equal with a relative
    precision tolerance of `rel_tol` and absolute tolerance `abs_tol`. Both
    default to `sqrt(F.epsilon())` (~1.49e-08 for `F64`, ~3.45e-04 for
    `F32`). All tolerances are measured using the complex modulus `abs()`.

    Two complex numbers `this` and `that` are considered almost equal when
    `(this - that).abs() <= max(rel_tol * max(this.abs(), that.abs()), abs_tol)`.

    When `this` and `that` moduli are large enough, then they are compared
    using the relative tolerance. One can write that `this` is within 5% of
    `that` with `this.almost_eq(that, 0.05)`. But when `this` and `that` are
    very small, using relative comparison makes no more sense, like when one
    value is 0.0. In that case, the absolute tolerance is used to check that
    both numbers are within that range.

    A complex number that is `NaN` is never equal to another number, according
    to IEEE 754 and so is not "almost equal" either. A complex number that is
    `infinite` is only almost equal to another similar infinite number.

    This function is mostly usefull when writing tests and that rounding and
    other types of errors accumulate, and one needs to accept results within
    predetermined tolerances.

    See Knuth TAOCP 4.2.2.A for approximation of real numbers.
    """
    let zero: F = F.from[F64](0.0)
    let one: F = F.from[F64](1.0)
    ifdef debug then
      (rel_tol >= zero) or
        Fail(Format("[Complex.almost_eq] Relative tolerance ({}) must be positive",
                    rel_tol))
      (rel_tol <= one) or
        Fail(Format("[Complex.almost_eq] Relative tolerance ({}) should be lower than 1.0",
                    rel_tol))
      (abs_tol >= zero) or
        Fail(Format("[Complex.almost_eq] Absolute tolerance ({}) must be positive",
                    abs_tol))
    end

    if nan() or that.nan() then
      false
    elseif infinite() then
      if that.infinite() then
        if (one.copysign(_re) == one.copysign(that._re)) and
          (one.copysign(_im) == one.copysign(that._im)) then
          true
        else
          false
        end
      else
        false
      end
    else
      (this - that).abs() <= (rel_tol * this.abs().max(that.abs())).max(abs_tol)
    end


  fun add(that: box->Complex[F]): Complex[F]^ =>
    """
    Addition of two complex numbers. Result is
    `(this.rea1() + that.real()) + (j() * (this.imag() + that.imag()))`.
    """
    Complex[F](_re + that._re, _im + that._im)


  fun add_unsafe(that: box->Complex[F]): Complex[F]^ =>
    """
    Usafe addition of two complex numbers using usafe arithmetic. Result is
    `(this.rea1() +~ that.real()) + (j() * (this.imag() +~ that.imag()))`.
    """
    Complex[F](_re +~ that._re, _im +~ that._im)


  fun sub(that: box->Complex[F]): Complex[F]^ =>
    """
    Substraction of two complex numbers. Result is
    `(this.rea1() - that.real()) + (j() * (this.imag() - that.imag()))`.
    """
    Complex[F](_re - that._re, _im - that._im)


  fun sub_unsafe(that: box->Complex[F]): Complex[F]^ =>
    """
    Unsafe substraction of two complex numbers using usafe arithmetic. Result is
    `(this.rea1() -~ that.real()) + (j() * (this.imag() -~ that.imag()))`.
    """
    Complex[F](_re -~ that._re, _im -~ that._im)


  fun neg(): Complex[F]^ =>
    """
    Negation of complex.
    """
    Complex[F](-_re, -_im)


  fun neg_unsafe(): Complex[F]^ =>
    """
    Negation of complex, using unsafe arithmetic.
    """
    Complex[F](-~_re, -~_im)


  fun mul(that: box->Complex[F]): Complex[F]^ =>
    """
    Multiplication of two complex numbers. Result is
    `((this.real() * that.real()) - (this.imag() * that.imag())) + (j() * ((this.imag() * that.real()) + (this.real() * that.imag()))`
    but calculation is done in less operations.
    
    This does not report an error if the calculation overflows.
    """
    let ac = _re * that._re
    let bd = _im * that._im
    let abcd = (_re + _im) * (that._re + that._im)
    Complex[F](ac - bd, abcd - ac - bd)


  fun mul_unsafe(that: box->Complex[F]): Complex[F]^ =>
    """
    Unsafe multiplication of two complex numbers using arithmetic. Result is
    `((this.real() *~ that.real()) -~ (this.imag() *~ that.imag())) + (j() * ((this.imag() *~ that.real()) +~ (this.real() *~ that.imag()))`
    but actual calculation is done in less operations.

    This does not report an error if the calculation overflows.
    """
    let ac = _re *~ that._re
    let bd = _im *~ that._im
    let abcd = (_re +~ _im) *~ (that._re +~ that._im)
    Complex[F](ac -~ bd, abcd -~ ac -~ bd)


  fun scale(s: F): Complex[F]^ =>
    """
    Multiply by the real scalar `s`. Equivalent to `mul(Complex(s, 0))` but
    avoids the overhead of a full complex multiplication.
    """
    Complex[F](_re * s, _im * s)


  fun scale_unsafe(s: F): Complex[F]^ =>
    """
    Multiply by the real scalar `s` using unsafe arithmetic.
    """
    Complex[F](_re *~ s, _im *~ s)


  fun abs(): F =>
    """
    Calculate the modulus of complex number. That is `sqrt((real() * real()) + (imag() * imag()))`,
    trying not to overflow.
    """
    let a = _re.abs()
    let b = _im.abs() 
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    if a == zero then
      b
    elseif b == zero then
      a
    elseif a >= b then
      let tmp = b / a
      a * (one + (tmp * tmp)).sqrt()
    else
      let tmp = a / b
      b * (one + (tmp * tmp)).sqrt()
    end
    
    
  fun abs_unsafe(): F =>
    """
    Calculate the unsafe modulus of complex number using unsafe arithmentic.
    That is `sqrt((real() *~ real()) +~ (imag() *~ imag()))`,
    trying not to overflow.
    """
    let a = _re.abs()
    let b = _im.abs()
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    if a ==~ zero then
      b
    elseif b ==~ zero then
      a
    elseif a >=~ b then
      let tmp = b /~ a
      a *~ (one +~ (tmp *~ tmp)).sqrt_unsafe()
    else
      let tmp = a /~ b
      b *~ (one +~ (tmp *~ tmp)).sqrt_unsafe()
    end


  fun abs2(): F =>
    """
    Returns the squared modulus `re² + im²`, avoiding the `sqrt` in `abs()`.
    Satisfies `abs2() == abs() * abs()`.

    Useful when comparing magnitudes (avoids sqrt) or in performance-critical
    inner loops where only relative magnitudes matter.
    """
    (_re * _re) + (_im * _im)


  fun abs2_unsafe(): F =>
    """
    Returns the squared modulus using unsafe arithmetic.
    """
    (_re *~ _re) +~ (_im *~ _im)
    
    
  fun invert(): Complex[F]^ =>
    """
    Inverse of complex, `1 / this`.
    
    This does not report an error if the calculation overflows.
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let inf = F.from[F64](1.0 / 0.0)

    if (_re == zero) and (_im == zero) then
      Complex[F](inf, inf)
    elseif _im == zero then
      Complex[F](one / _re, zero)
    elseif _re == zero then
      Complex[F](zero, -one / _im)
    elseif _re.abs() >= _im.abs() then
      let dc = _im / _re
      let denom = _re + (_im * dc)
      Complex[F](one / denom, -dc / denom)
    else
      let cd = _re / _im
      let denom = (_re * cd) + _im
      Complex[F](cd / denom, -_re / denom)
    end
    
  
  fun invert_unsafe(): Complex[F]^ =>
    """
    Inverse of complex, `1 / this`, using unsafe arithmetic.
    
    This does not report an error if the calculation overflows.
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let inf = F.from[F64](1.0 / 0.0)

    if (_re ==~ zero) and (_im ==~ zero) then
      Complex[F](inf, inf)
    elseif _im ==~ zero then
      Complex[F](one /~ _re, zero)
    elseif _re ==~ zero then
      Complex[F](zero, -one /~ _im)
    elseif _re.abs() >=~ _im.abs() then
      let dc = _im /~ _re
      let denom = _re +~ (_im *~ dc)
      Complex[F](one /~ denom, -~dc /~ denom)
    else
      let cd = _re /~ _im
      let denom = (_re *~ cd) +~ _im
      Complex[F](cd /~ denom, -~_re /~ denom)
    end
    
  
  fun div(that: box->Complex[F]): Complex[F]^ =>
    """
    Division of two complex numbers `this / that`.
    
    This does not report an error if the calculation overflows but an infinite
    complex.
    """
    let zero = F.from[ISize](0)
    let inf = F.from[F64](1.0 / 0.0)

    if (that._re == zero) and (that._im == zero) then
      Complex[F](inf.copysign(_re), inf.copysign(_im))
    elseif that._im == zero then
      Complex[F](_re / that._re, _im / that._re)
    elseif that._re == zero then
      Complex[F](_im / that._im, -_re / that._im)
    elseif that._re.abs() >= that._im.abs() then
      let dc = that._im / that._re
      let denom = that._re + (that._im * dc)
      Complex[F]((_re + (_im * dc)) / denom, (_im - (_re * dc)) / denom)
    else
      let cd = that._re / that._im
      let denom = (that._re * cd) + that._im
      Complex[F](((_re * cd) + _im) / denom, ((_im * cd) - _re) / denom)
    end
    
  
  fun div_unsafe(that: box->Complex[F]): Complex[F]^ =>
    """
    Unsafe division of two complex numbers `this /~ that`.

    It uses usafe arithmetic operations.    
    This does not report an error if the calculation overflows but returns an
    infinite complex.
    """
    let zero = F.from[ISize](0)
    let inf = F.from[F64](1.0 / 0.0)

    if (that._re ==~ zero) and (that._im ==~ zero) then
      Complex[F](inf.copysign(_re), inf.copysign(_im))
    elseif that._im ==~ zero then
      Complex[F](_re /~ that._re, _im /~ that._re)
    elseif that._re.abs() >=~ that._im.abs() then
      let dc = that._im /~ that._re
      let denom = that._re +~ (that._im *~ dc)
      Complex[F]((_re +~ (_im *~ dc)) /~ denom, (_im -~ (_re *~ dc)) /~ denom)
    else
      let cd = that._re /~ that._im
      let denom = (that._re *~ cd) +~ that._im
      Complex[F](((_re *~ cd) +~ _im) /~ denom, ((_im *~ cd) -~ _re) /~ denom)
    end
    
  
  fun sqrt(): Complex[F]^ =>
    """
    Square root of the complex.
    """
    let zero = F.from[ISize](0)
    if (_re == zero) and (_im == zero) then
      return Complex[F](zero, zero)
    end

    let c = _re.abs()
    let d = _im.abs()
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    
    let w =
      if c >= d then
        let dc = d / c
        c.sqrt() * ((one + (one + (dc * dc)).sqrt()) / two).sqrt()
      else
        let cd = c / d
        d.sqrt() * ((cd + (one + (cd * cd)).sqrt()) / two).sqrt()
      end
    
    if _re >= zero then
      Complex[F](w, _im / (two * w))
    elseif (_re < zero) and (_im >= zero) then
      Complex[F](d / (two * w), w)
    else
      Complex[F](d / (two * w), -w)
    end
    
  
  fun sqrt_unsafe(): Complex[F]^ =>
    """
    Unsafe square root of the complex, using unsafe arithmetic operations.
    """
    let zero = F.from[ISize](0)
    if (_re ==~ zero) and (_im ==~ zero) then
      return Complex[F](zero, zero)
    end

    let c = _re.abs()
    let d = _im.abs()
    let one = F.from[ISize](1)
    let two = F.from[ISize](2)
    
    let w =
      if c >=~ d then
        let dc = d /~ c
        c.sqrt() *~ ((one +~ (one +~ (dc *~ dc)).sqrt()) /~ two).sqrt()
      else
        let cd = c /~ d
        d.sqrt() *~ ((cd +~ (one +~ (cd *~ cd)).sqrt()) /~ two).sqrt()
      end
    
    if _re >=~ zero then
      Complex[F](w, _im /~ (two *~ w))
    elseif (_re <~ zero) and (_im >=~ zero) then
      Complex[F](d /~ (two *~ w), w)
    else
      Complex[F](d /~ (two *~ w), -~w)
    end
    
  
  fun arg(): F =>
    """
    The argument (angle) of the complex number, in (-π, π].

    `arg(re + i*im) = atan2(im, re)`
    """
    _im.atan2(_re)
    
    
  fun exp(): Complex[F]^ =>
    """
    Exponential of complex.

    `exp(re + i*im) = e^re * (cos(im) + i*sin(im))`
    """
    let e_re = _re.exp()
    Complex[F](e_re * _im.cos(), e_re * _im.sin())
    
    
  fun exp_unsafe(): Complex[F]^ =>
    """
    Exponential of complex using unsafe arithmetic.

    `exp(re + i*im) = e^re *~ (cos(im) + i*sin(im))`
    """
    let e_re = _re.exp()
    Complex[F](e_re *~ _im.cos(), e_re *~ _im.sin())
 

  fun cos(): Complex[F]^ =>
    """
    Cosine of complex.

    `z.cos() = ((i * z).exp() + (-i * z).exp())) / 2`
    """
    let jz = j() * Complex[F](_re, _im)
    let two = Complex[F](F.from[ISize](2), F.from[ISize](0))
    (jz.exp() + ((-jz).exp())) / two


  fun sin(): Complex[F]^ =>
    """
    Sine of complex.

    `z.sin() = ((i * z).exp() - (-i * z).exp()) / (2 * i)`
    """
    let jz = j() * Complex[F](_re, _im)
    let twoi = Complex[F](F.from[ISize](0), F.from[ISize](2))
    (jz.exp() - ((-jz).exp())) / twoi


  fun tan(): Complex[F]^ =>
    """
    Tangent of complex.

    `z.tan() = z.sin() / z.cos()`
    """
    sin() / cos()


  fun asin(): Complex[F]^ =>
    """
    Inverse sine of complex.

    `asin(z) = -i * log(iz + sqrt(1 - z²))`
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let c_one = Complex[F](one, zero)
    let neg_i = Complex[F](zero, -one)
    // Explicit val copy of this to allow arithmetic (fun receiver is box)
    let z = Complex[F](_re, _im)
    // iz = i * (re + i*im) = -im + i*re
    let iz = Complex[F](-_im, _re)
    let z2 = z * z
    let w = (c_one - z2).sqrt()
    neg_i * (iz + w).log()


  fun acos(): Complex[F]^ =>
    """
    Inverse cosine of complex.

    `acos(z) = -i * log(z + i * sqrt(1 - z²))`
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let c_one = Complex[F](one, zero)
    let neg_i = Complex[F](zero, -one)
    let z = Complex[F](_re, _im)
    let z2 = z * z
    let w = (c_one - z2).sqrt()
    // i * w = (-w._im + i*w._re)
    neg_i * (z + Complex[F](-w._im, w._re)).log()


  fun atan(): Complex[F]^ =>
    """
    Inverse tangent of complex.

    `atan(z) = (i/2) * log((i + z) / (i - z))`
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let half_i = Complex[F](zero, one / F.from[ISize](2))
    let i = Complex[F](zero, one)
    let z = Complex[F](_re, _im)
    half_i * ((i + z) / (i - z)).log()


  fun cosh(): Complex[F]^ =>
    """
    Hyperbolic cosine of complex.

    `z.cosh() = (e^z + e^-z) / 2 = (i * z).cos()`
    """
    let two = Complex[F](F.from[ISize](2), F.from[ISize](0))
    (this.exp() + ((-this).exp())) / two


  fun sinh(): Complex[F]^ =>
    """
    Hyperbolic sine of complex.

    `z.sinh() = (e^z - e^-z) / 2 = (-i * (i * z).sin())`
    """
    let two = Complex[F](F.from[ISize](2), F.from[ISize](0))
    (this.exp() - ((-this).exp())) / two


  fun tanh(): Complex[F]^ =>
    """
    Hyperbolic tangent of complex.

    `z.tanh() = z.sinh() / z.cosh()`
    """
    sinh() / cosh()


  fun asinh(): Complex[F]^ =>
    """
    Inverse hyperbolic sine of complex.

    `asinh(z) = log(z + sqrt(z² + 1))`
    """
    let one = F.from[ISize](1)
    let c_one = Complex[F](one, F.from[ISize](0))
    let z = Complex[F](_re, _im)
    let z2 = z * z
    (z + (z2 + c_one).sqrt()).log()


  fun acosh(): Complex[F]^ =>
    """
    Inverse hyperbolic cosine of complex.

    `acosh(z) = log(z + sqrt(z² - 1))`
    """
    let one = F.from[ISize](1)
    let c_one = Complex[F](one, F.from[ISize](0))
    let z = Complex[F](_re, _im)
    let z2 = z * z
    (z + (z2 - c_one).sqrt()).log()


  fun atanh(): Complex[F]^ =>
    """
    Inverse hyperbolic tangent of complex.

    `atanh(z) = (1/2) * log((1 + z) / (1 - z))`
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)
    let c_one = Complex[F](one, zero)
    let half = Complex[F](one / F.from[ISize](2), zero)
    let z = Complex[F](_re, _im)
    half * ((c_one + z) / (c_one - z)).log()


  fun powi(n: I32): Complex[F]^ =>
    """
    Integer power of complex.
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    // Exponent special values
    if n == 0 then
      return Complex[F](one, zero)
    elseif n == 1 then
      return Complex[F](_re, _im)
    elseif n == -1 then
      return this.invert()
    end

    // Treat real and imaginary as special cases
    if _im == zero then
      return Complex[F](_re.powi(n), zero)
    elseif _re == zero then
      match n % 4
      | -3 => return Complex[F](zero, _im.powi(n))
      | -2 => return Complex[F](-_im.powi(n), zero)
      | -1 => return Complex[F](zero, -_im.powi(n))
      |  0 => return Complex[F](_im.powi(n), zero)
      |  1 => return Complex[F](zero, _im.powi(n))
      |  2 => return Complex[F](-_im.powi(n), zero)
      |  3 => return Complex[F](zero, -_im.powi(n))
      else
        // Can't happen
        return Complex[F](zero, zero)
      end
    end

    // General case
    var i = if n > 0 then n else -n end
    var result = Complex[F](one, zero)
    var x: Complex[F] = Complex[F](_re, _im)
    while i > 0 do
      if (i % 2) == 1 then
        result = result * x
      end
      x = x * x
      i = i / 2
    end

    if n > 0 then
      result
    else
      result.invert()
    end
    
    
  fun powi_unsafe(n: I32): Complex[F]^ =>
    """
    Integer power of complex using unsafe arithmetic.
    """
    let zero = F.from[ISize](0)
    let one = F.from[ISize](1)

    // Exponent special values
    if n == 0 then
      return Complex[F](one, zero)
    elseif n == 1 then
      return Complex[F](_re, _im)
    elseif n == -1 then
      return this.invert()
    end

    // Treat real and imaginary as special cases
    if _im ==~ zero then
      return Complex[F](_re.powi(n), zero)
    elseif _re ==~ zero then
      match n % 4
      | -3 => return Complex[F](zero, _im.powi(n))
      | -2 => return Complex[F](-_im.powi(n), zero)
      | -1 => return Complex[F](zero, -_im.powi(n))
      |  0 => return Complex[F](_im.powi(n), zero)
      |  1 => return Complex[F](zero, _im.powi(n))
      |  2 => return Complex[F](-_im.powi(n), zero)
      |  3 => return Complex[F](zero, -_im.powi(n))
      else
        // Can't happen
        return Complex[F](zero, zero)
      end
    end

    // General case
    var i = if n > 0 then n else -n end
    var result: Complex[F] = Complex[F](one, zero)
    var x: Complex[F] = Complex[F](_re, _im)
    while i > 0 do
      if (i %~ 2) == 1 then
        result = result *~ x
      end
      x = x *~ x
      i = i /~ 2
    end

    if n > 0 then
      result
    else
      result.invert_unsafe()
    end


  fun log(): Complex[F]^ =>
    """
    Natural logarithm of complex (principal value).

    `log(z) = ln|z| + i * arg(z)`
    """
    Complex[F](abs().log(), arg())


  fun log2(): Complex[F]^ =>
    """
    Logarithm base 2 of complex: `log(this) / log(2)`.

    `log2(z).real() = log2(|z|)`, `log2(z).imag() = arg(z) / ln(2)`
    """
    let ln2 = F.from[ISize](2).log()
    Complex[F](abs().log() / ln2, arg() / ln2)


  fun log10(): Complex[F]^ =>
    """
    Logarithm base 10 of complex: `log(this) / log(10)`.

    `log10(z).real() = log10(|z|)`, `log10(z).imag() = arg(z) / ln(10)`
    """
    let ln10 = F.from[ISize](10).log()
    Complex[F](abs().log() / ln10, arg() / ln10)


  fun pow(that: box->Complex[F]): Complex[F]^ =>
    """
    Complex power: `this^that = exp(that * log(this))`.
    """
    (that * log()).exp()


  fun powf(r: F): Complex[F]^ =>
    """
    Power with a real exponent: `this^r = |this|^r * exp(i * r * arg(this))`.

    More efficient than `pow(Complex(r, 0))`. Consistent with polar form:
    the result has modulus `|this|^r` and argument `r * arg(this)`.
    """
    let rho = abs().pow(r)
    let theta = arg() * r
    Complex[F](rho * theta.cos(), rho * theta.sin())


  fun dotp(that: box->Complex[F]): F =>
    """
    Scalar product of `this` and `that`,
    `(this.real() * that.real()) + (this.imag() * that.imag())`
    """
    (_re * that._re) + (_im * that._im)
  
  
  fun dotp_unsafe(that: box->Complex[F]): F =>
    """
    Unsafe scalar product of `this` and `that`,
    `this.real() * that.real() + this.imag() * that.imag()`
    """
    (_re *~ that._re) +~ (_im *~ that._im)

