// Tests for arbitrary-precision integers MPInt
// There are 2 implementations:
//   1. Pure Pony implementation
//   2. Binding to GMP library
// These tests can be used with both implementations, when you change the use
// clause.

// When using GMP/MPFR implementation, change use
use "../mathx"
//use "../mathx/gmp"

use "../pony_testx"

use "collections"
use "random"


class iso _TestMPIntCreate is UnitTest
  """
  Tests on MPInt create
  """

  fun name(): String =>
    "MPInt/create"

  fun apply(h: TestHelper) =>
    h.assert_no_error({() ? => MPInt()? }, "Valid number 0")
    h.assert_no_error({() ? => MPInt("0")? }, "Valid number 0")
    h.assert_no_error({() ? => MPInt("+0")? }, "Valid number 0")
    h.assert_no_error({() ? => MPInt("-0")? }, "Valid number 0")
    h.assert_no_error({() ? => MPInt("-123_456_789")? }, "With separators")
    h.assert_no_error({() ? => MPInt("-123__456__789")? }, "With separators")
    h.assert_no_error({() ? => MPInt("-1_2_3__4_5_6__7_8_9")? }, "With separators")
    h.assert_no_error({() ? => MPInt("-123_456_789@0")? }, "Null exponent")
    h.assert_no_error({() ? => MPInt("-123_456_789@+0")? }, "Null exponent")
    h.assert_no_error({() ? => MPInt("+123_456_789@12")? }, "With exponent")
    h.assert_no_error({() ? => MPInt("-123_456_789e12")? }, "With exponent")
    h.assert_no_error({() ? => MPInt("+123_456_789E12")? }, "With exponent")
    h.assert_no_error({() ? => MPInt("+900d_F00d@+12", 16)? }, "Base 16")
    h.assert_no_error({() ? => MPInt("-1001010010010010111000101110100001110111000101010001101110@+50", 2)? }, "Base 2")
    h.assert_no_error({() ? => MPInt("-1001010010010010111000101110100001110111000101010001101110e+50", 2)? }, "Base 2")
    h.assert_no_error({() ? => MPInt("  -123_456_789")? }, "Leading spaces")
    h.assert_no_error({() ? => MPInt("-e50")? }, "Only exponent without mantissa")

    h.assert_error({() ? => MPInt("-")? }, "Only sign")
    h.assert_error({() ? => MPInt("+")? }, "Only sign")
    h.assert_error({() ? => MPInt("     ")? }, "Empty != 0")
    h.assert_error({() ? => MPInt("123_456_789 ")? }, "Trailing spaces")
    h.assert_error({() ? => MPInt("+-123")? }, "Multiple signs")
    h.assert_error({() ? => MPInt("-12.3")? }, "Float")
    h.assert_error({() ? => MPInt("12e")? }, "Empty exponent")
    h.assert_error({() ? => MPInt("12E+__")? }, "Empty exponent")
    h.assert_error({() ? => MPInt("+900d_F00dE+12", 16)? }, "Wrong exponent")
    h.assert_error({() ? => MPInt("+123_456_789E-12")? }, "Negative exponent")
    h.assert_error({() ? => MPInt("-900d_F00d@-12", 16)? }, "Negative exponent")

    try
      h.assert_true(MPInt("  +123_456_789")?.string() == "123456789")
    else
      h.fail("Wrong value")
    end

    try
      h.assert_true(MPInt("-123e6")?.string() == "-123000000")
    else
      h.fail("Wrong value")
    end

    try
      h.assert_true(MPInt()? == MPInt.from_ilong(), "Valid number 0")
      h.assert_true(MPInt("0")? == MPInt.from_ilong(0), "Valid number 0")
      // Pony parser does not support unary +
      h.assert_true(MPInt("+0")? == MPInt.from_ilong(0), "Valid number 0")
      h.assert_true(MPInt("-0")? == MPInt.from_ilong(-0), "Valid number 0")
      h.assert_true(MPInt("-123_456_789")? == MPInt.from_ilong(-123_456_789), "With separators")
      h.assert_true(MPInt("-123__456__789")? == MPInt.from_ilong(-123_456_789), "With separators")
      h.assert_true(MPInt("-1_2_3__4_5_6__7_8_9")? == MPInt.from_ilong(-1_2_3_4_5_6_7_8_9), "With separators")
      h.assert_true(MPInt("-123_456_789@0")? == MPInt.from_ilong(-123_456_789), "Null exponent")
      h.assert_true(MPInt("-123_456_789@+0")? == MPInt.from_ilong(-123_456_789), "Null exponent")
      h.assert_true(MPInt("+123_456_789@5")? == MPInt.from_ilong(123_456_789_00000), "With exponent")
      h.assert_true(MPInt("-123_456_789e3")? == MPInt.from_ilong(-123_456_789_000), "With exponent")
      h.assert_true(MPInt("+123_456_789E4")? == MPInt.from_ilong(123_456_789_0000), "With exponent")
      h.assert_true(MPInt("+900d_F00d@+2", 16)? == MPInt.from_ilong(0x900d_F00d * 16 * 16), "Base 16")
      h.assert_true(MPInt("-100101001001001011100010111010000111011100010101@+12", 2)? == MPInt.from_ilong(-669116797016952832), "Base 2")
      h.assert_true(MPInt("-100101001001001011100010111010000111011100010101e+12", 2)? == MPInt.from_ilong(-669116797016952832), "Base 2")
      h.assert_true(MPInt("  -123_456_789")? == MPInt.from_ilong(-123_456_789), "Leading spaces")
      h.assert_true(MPInt("-1e10")? == MPInt.from_ilong(-10000000000), "Exponent")
      h.assert_true(MPInt("18795432178911876e21")?.string() == "18795432178911876000000000000000000000", "String representation")
      h.assert_true(MPInt("-37892648294")? == MPInt.from_ilong(-37892648294))
    else
      h.fail("Equalities")
    end

    try
      let rand = Rand()
      for i in Range(0, 1000) do
        let n = rand.ilong()
        let s: String = n.string()
        h.assert_true(MPInt(s)? == MPInt.from_ilong(n), "Random values")
      end
    else
      h.fail("Random values")
    end


class iso _TestMPIntFromILong is UnitTest
  """
  Tests on MPInt from_ilong
  """

  fun name(): String =>
    "MPInt/from_ilong"

  fun apply(h: TestHelper) =>
    let minv = MPInt.from_ilong(ILong.min_value())
    let maxv = MPInt.from_ilong(ILong.max_value())

    try
      h.assert_true((minv == MPInt("-8000_0000", 16)?) or (minv == MPInt("-8000_0000_0000_0000", 16)?) , "Min value")
      h.assert_true((maxv == MPInt("7fff_ffff", 16)?) or (maxv == MPInt("7fff_ffff_ffff_ffff", 16)?) , "Max value")
    else
      h.fail("Min and max values")
    end


class iso _TestMPIntComparisons is UnitTest
  """
  Tests on MPInt comparisons
  """

  fun name(): String =>
    "MPInt/comparisons"

  fun apply(h: TestHelper) =>
    try
      let zero = MPInt.from_ilong(0)
      let zerop = MPInt("+0")?
      let zerom = MPInt("-0")?
      h.assert_true(zero.zero(), "Zero")
      h.assert_true(zerop.zero(), "Zero")
      h.assert_true(zerom.zero(), "Zero")

      let rand = Rand()
      for i in Range(0, 1000) do
        let n = rand.ilong()
        let n' = MPInt.from_ilong(n)
        let m = rand.ilong()
        let m' = MPInt.from_ilong(m)

        h.assert_false((n == m) and (n' != m'), "Random comparisons")
        h.assert_false((n < m) and (n' >= m'), "Random comparisons")
        h.assert_true(n' == n', "Random equality")
        h.assert_true(m' == m', "Random equality")
        if n < m then
          h.assert_true(n' < m', "<")
          h.assert_true(m' > n', ">")
          h.assert_true(n' != m', "!=")
        end
        if n <= m then
          h.assert_true(n' <= m', "<=")
          h.assert_true(m' >= n', ">=")
        end
        if n == m then
          h.assert_true(n' == m', "==")
        end

        let s = rand.u8()
        var x: String val = String(s.usize() * 10)
        var y: String val = String(s.usize() * 10)
        if rand.i8() >= 0 then
          x = "+".clone()
        else
          x = "-".clone()
        end
        if rand.i8() >= 0 then
          y = "+".clone()
        else
          y = "-".clone()
        end
        for j in Range[U8](0, s) do
          x = x + rand.u32().string()
          y = y + rand.u32().string()
        end
        if x == "+" then x = "+1" elseif x == "-" then x = "-9" end
        if y == "+" then y = "+8" elseif y == "-" then y = "-5" end
        let x' = MPInt(x)?
        let y' = MPInt(y)?

        if x(0)? == '+' then
          h.assert_true(x' >= zero, "Positive")
          if y(0)? == '+' then
            h.assert_true(y' >= zero, "Positive")
          else
            h.assert_true(x' >= y', "Greater or equal")
          end
        else
          h.assert_true(x' <= zero, "Negative")
          if y(0)? == '+' then
            h.assert_true(x' <= y', "Less or equal")
          else
            h.assert_true(y' <= zero, "Negative")
          end
        end
      end
    else
      h.fail("Random values")
    end


class iso _TestMPIntAbsComparisons is UnitTest
  """
  Tests on MPInt absolute values comparisons
  """

  fun name(): String =>
    "MPInt/abs_comparisons"

  fun apply(h: TestHelper) =>
    let rand = Rand()
    for i in Range(0, 1000) do
      let v = rand.ilong()
      let a = MPInt.from_ilong(v)
      h.assert_true(a.abs_eq(a), "Absolute value equality")
      let b = MPInt.from_ilong(v + 1)
      h.assert_false(a.abs_eq(b), "Absolute value equality")
    end

    for i in Range(0, 100) do
      let a = rand.ilong()
      let b = rand.ilong()
      let a' = MPInt.from_ilong(a)
      let b' = MPInt.from_ilong(b)
      h.assert_true(not (a.abs() < b.abs()) or (a'.abs_lt(b')), "Absolute value less")
      h.assert_true(not (a.abs() <= b.abs()) or (a'.abs_le(b')), "Absolute value less or equal")
      h.assert_true(not (a.abs() == b.abs()) or (a'.abs_eq(b')), "Absolute value equality")
      h.assert_true(not (a.abs() != b.abs()) or (a'.abs_ne(b')), "Aboslute value difference")
      h.assert_true(not (a.abs() >= b.abs()) or (a'.abs_ge(b')), "Absolute value greater or equal")
      h.assert_true(not (a.abs() > b.abs()) or (a'.abs_gt(b')), "Absolute value greater")
    end


class iso _TestMPIntConversionILong is UnitTest
  """
  Tests on MPInt conversion to ILong
  """

  fun name(): String =>
    "MPInt/ilong"

  fun apply(h: TestHelper) =>
    let rand = Rand()
    for i in Range(0, 1000) do
      let a = rand.ilong()
      let a' = MPInt.from_ilong(a)
      try
        h.assert_true(a == a'.ilong()?, "Conversion to ILong")
      else
        h.fail("Conversion from ILong to ILong")
      end
    end

    try
      let b = MPInt("89334718795878195613438179751988712")?
      h.assert_error({() ? => let b' = b.ilong()? }, "Conversion to large positive ILong")

      let c = MPInt("-19815798458789751981278742378912370")?
      h.assert_error({() ? => let c' = c.ilong()? }, "Conversion to large negative ILong")
    else
      h.fail("Conversion to large ILong")
    end


class iso _TestMPIntMiscellaneous is UnitTest
  """
  Tests on MPInt miscellaneous operations
  """

  fun name(): String =>
    "MPInt/miscellaneous"

  fun apply(h: TestHelper) =>
    let one = MPInt.from_ilong(1)
    h.assert_true((one << 1) == MPInt.from_ilong(65536), "Shift left by 1")
    h.assert_true((one << 2) == MPInt.from_ilong(4294967296), "Shift left by 2")
    h.assert_true((one << 3) == MPInt.from_ilong(281474976710656), "Shift left by 3")
    h.assert_true(((one << 10) >> 10) == one, "Shift left then right")

    let one3 = MPInt.from_ilong(281474976710656)
    let one2 = MPInt.from_ilong(4294967296)
    let one1 = MPInt.from_ilong(65536)
    let zero = MPInt.from_ilong(0)
    h.assert_true((one3 >> 3) == one, "Shift right by 3")
    h.assert_true((one2 >> 2) == one, "Shift right by 2")
    h.assert_true((one1 >> 1) == one, "Shift right by 1")
    h.assert_true((one >> 1) == zero, "Shift right to 0")


class iso _TestMPIntArithmetic is UnitTest
  """
  Tests on MPInt arithmetic operations
  """

  fun name(): String =>
    "MPInt/arithmetic"

  fun apply(h: TestHelper) =>
    let rand = Rand()

    // Addition
    for i in Range(0, 1000) do
      let a = rand.ilong()
      let a' = MPInt.from_ilong(a)
      let b = rand.ilong()
      let b' = MPInt.from_ilong(b)
      try
        let c = a +? b
        let c' = MPInt.from_ilong(c)
        h.assert_true(c' == (a' + b'), "Addition")
      else
        h.log("Overflow with " + a.string() + " + " + b.string() + ", trying other values...")
      end
    end

    try
      h.assert_true((MPInt("+123456")? + MPInt("+1234567890123456789")?) == MPInt("+1234567890123580245")?)
      h.assert_true((MPInt("-123456")? + MPInt("-1234567890123456789")?) == MPInt("-1234567890123580245")?)
      h.assert_true((MPInt("+123456")? + MPInt("-1234567890123456789")?) == MPInt("-1234567890123333333")?)
      h.assert_true((MPInt("-123456")? + MPInt("+1234567890123456789")?) == MPInt("+1234567890123333333")?)
    else
      h.fail("Addition of big and small")
    end

    try
      h.assert_true((MPInt("19286568262396999950202928044648128737703949436999793610613130182105922367448423453159208008561154643649004993410245621815263683558649945634758267696270039303164684506235887452943484121027471056665201277603534728479727382153288066446108997401020223592515566263174704019778189859650264257813887322843943022301697122037047349594826189769800767532590981681753960161163017803752421139870163553766644589335272831995013116481436362323269375288953950813601638251442829597903781688297881241170668619904757069020265297394635546155851093309363569785849294570445068455981948468662793917426384682117599755965209933018194164415994225160507289598509847921")? + MPInt("-89397816278271230459454110395996026663027374858959155341995445324393903731711371299662676953068781699910566197767070398021792037309382471935003522844234398194213451613830240178084330862239920913616552830617626796686980142144475041619332654915691024934414728326254713483694918659373308283669427344580771115999103149093958866604351791680623545304481448200580679477784543921552472823363295771446882280964663265344973021812416638734557654215643369337559645949925843413734820414989440796460578313739641380350952454004469826268119514126221846674571499744339378765148723031976486968205020823203122321431396802517326793155991667238272944727233749")?) == MPInt("19197170446118728719743473934252132711040922062140834455271134736781528463716712081859545331608085861949094427212478551417241891521340563162823264173425804904970471054622057212765399790165231135751584724772917101683040402011143591404489664746104532567581151534848449306294494940990890949530217895499362251185698018887953390728221837978120143987286500233553379481685233259830868667046800257995197707054308168729668143459623945684534817634738307444264078605492903754490046867882891800374208041591017427639914344940631076329582973795237347939174723070700729077216799745630817430458179661294396633643778536215676837622838233493269016653782614172")?, "Addition big numbers +-")
      h.assert_true((MPInt("19286568262396999950202928044648128737703949436999793610613130182105922367448423453159208008561154643649004993410245621815263683558649945634758267696270039303164684506235887452943484121027471056665201277603534728479727382153288066446108997401020223592515566263174704019778189859650264257813887322843943022301697122037047349594826189769800767532590981681753960161163017803752421139870163553766644589335272831995013116481436362323269375288953950813601638251442829597903781688297881241170668619904757069020265297394635546155851093309363569785849294570445068455981948468662793917426384682117599755965209933018194164415994225160507289598509847921")? + MPInt("89397816278271230459454110395996026663027374858959155341995445324393903731711371299662676953068781699910566197767070398021792037309382471935003522844234398194213451613830240178084330862239920913616552830617626796686980142144475041619332654915691024934414728326254713483694918659373308283669427344580771115999103149093958866604351791680623545304481448200580679477784543921552472823363295771446882280964663265344973021812416638734557654215643369337559645949925843413734820414989440796460578313739641380350952454004469826268119514126221846674571499744339378765148723031976486968205020823203122321431396802517326793155991667238272944727233749")?) == MPInt("19375966078675271180662382155044124764366976811858752765955125627430316271180134824458870685514223425348915559608012692213285475595959328106693271219114273701358897957849717693121568451889710977578817830434152355276414362295432541487728330055935914617449980991500958733261884778309637566097556750188523793417696225186141308461430541561481391077895463129954540840640802347673973612693526849538091471616237495260358089503248778962003932943169594182939197897392755441317516508712870681967129198218496710400616249848640015982119212823489791632523866070189407834747097191694770404394589702940802878286641329820711491209150216827745562543237081670")?, "Addition big numbers ++")
      h.assert_true((MPInt("-19286568262396999950202928044648128737703949436999793610613130182105922367448423453159208008561154643649004993410245621815263683558649945634758267696270039303164684506235887452943484121027471056665201277603534728479727382153288066446108997401020223592515566263174704019778189859650264257813887322843943022301697122037047349594826189769800767532590981681753960161163017803752421139870163553766644589335272831995013116481436362323269375288953950813601638251442829597903781688297881241170668619904757069020265297394635546155851093309363569785849294570445068455981948468662793917426384682117599755965209933018194164415994225160507289598509847921")? + MPInt("-89397816278271230459454110395996026663027374858959155341995445324393903731711371299662676953068781699910566197767070398021792037309382471935003522844234398194213451613830240178084330862239920913616552830617626796686980142144475041619332654915691024934414728326254713483694918659373308283669427344580771115999103149093958866604351791680623545304481448200580679477784543921552472823363295771446882280964663265344973021812416638734557654215643369337559645949925843413734820414989440796460578313739641380350952454004469826268119514126221846674571499744339378765148723031976486968205020823203122321431396802517326793155991667238272944727233749")?) == MPInt("-19375966078675271180662382155044124764366976811858752765955125627430316271180134824458870685514223425348915559608012692213285475595959328106693271219114273701358897957849717693121568451889710977578817830434152355276414362295432541487728330055935914617449980991500958733261884778309637566097556750188523793417696225186141308461430541561481391077895463129954540840640802347673973612693526849538091471616237495260358089503248778962003932943169594182939197897392755441317516508712870681967129198218496710400616249848640015982119212823489791632523866070189407834747097191694770404394589702940802878286641329820711491209150216827745562543237081670")?, "Addition big numbers --")
      h.assert_true((MPInt("-19286568262396999950202928044648128737703949436999793610613130182105922367448423453159208008561154643649004993410245621815263683558649945634758267696270039303164684506235887452943484121027471056665201277603534728479727382153288066446108997401020223592515566263174704019778189859650264257813887322843943022301697122037047349594826189769800767532590981681753960161163017803752421139870163553766644589335272831995013116481436362323269375288953950813601638251442829597903781688297881241170668619904757069020265297394635546155851093309363569785849294570445068455981948468662793917426384682117599755965209933018194164415994225160507289598509847921")? + MPInt("89397816278271230459454110395996026663027374858959155341995445324393903731711371299662676953068781699910566197767070398021792037309382471935003522844234398194213451613830240178084330862239920913616552830617626796686980142144475041619332654915691024934414728326254713483694918659373308283669427344580771115999103149093958866604351791680623545304481448200580679477784543921552472823363295771446882280964663265344973021812416638734557654215643369337559645949925843413734820414989440796460578313739641380350952454004469826268119514126221846674571499744339378765148723031976486968205020823203122321431396802517326793155991667238272944727233749")?) == MPInt("-19197170446118728719743473934252132711040922062140834455271134736781528463716712081859545331608085861949094427212478551417241891521340563162823264173425804904970471054622057212765399790165231135751584724772917101683040402011143591404489664746104532567581151534848449306294494940990890949530217895499362251185698018887953390728221837978120143987286500233553379481685233259830868667046800257995197707054308168729668143459623945684534817634738307444264078605492903754490046867882891800374208041591017427639914344940631076329582973795237347939174723070700729077216799745630817430458179661294396633643778536215676837622838233493269016653782614172")?, "Addition big numbers -+")
    else
      h.fail("Addition of big numbers")
    end

    // Substraction
    for i in Range(0, 1000) do
      let a = rand.ilong()
      let a' = MPInt.from_ilong(a)
      let b = rand.ilong()
      let b' = MPInt.from_ilong(b)
      try
        let c = a -? b
        let c' = MPInt.from_ilong(c)
        h.assert_true(c' == (a' - b'), "Substraction")
      else
        h.log("Overflow with " + a.string() + " - " + b.string() + ", trying other values...")
      end
    end

    try
      h.assert_true((MPInt("+654321")? - MPInt("+9876543210987654321")?) == MPInt("-9876543210987000000")?)
      h.assert_true((MPInt("-654321")? - MPInt("-9876543210987654321")?) == MPInt("+9876543210987000000")?)
      h.assert_true((MPInt("+654321")? - MPInt("-9876543210987654321")?) == MPInt("+9876543210988308642")?)
      h.assert_true((MPInt("-654321")? - MPInt("+9876543210987654321")?) == MPInt("-9876543210988308642")?)
    else
      h.fail("Substraction of big and small")
    end

    // Test on boundaries
    h.assert_true((MPInt.from_ilong(ILong.max_value()) + MPInt.from_ilong(1)) == (MPInt.from_ilong(ILong.max_value() - 1) + MPInt.from_ilong(2)), "Max_value + 1")
    h.assert_true((MPInt.from_ilong(ILong.min_value()) - MPInt.from_ilong(1)) == (MPInt.from_ilong(ILong.min_value() + 1) - MPInt.from_ilong(2)), "Min_value - 1")

    // Multiplication
    for i in Range(0, 1000) do
      let a = rand.ilong()
      let a' = MPInt.from_ilong(a)
      let b = rand.ilong()
      let b' = MPInt.from_ilong(b)
      try
        let c = a *? b
        let c' = MPInt.from_ilong(c)
        h.assert_true(c' == (a' * b'), "Multiplication")
      else
        h.log("Overflow with " + a.string() + " * " + b.string() + ", trying other values...")
      end
    end


class iso _TestMPIntKaratsuba is UnitTest
  """
  Tests on MPInt Karatsuba multiplication algorithm
  """

  fun name(): String =>
    "MPInt/karatsuba"

  fun apply(h: TestHelper) =>
    let rand = Rand()
    // We need to generate big numbers with more than 128 base-digits.
    let size1: USize = 130
    let size2: USize = 130

    for _ in Range(0, 100) do
      var big1 = MPInt.from_ilong(rand.ilong())
      for i in Range(1, size1) do
        let n = MPInt.from_ilong(rand.ilong()) << i
        big1 = big1 + n
      end

      var big2 = MPInt.from_ilong(rand.ilong())
      for i in Range(1, size2) do
        let n = MPInt.from_ilong(rand.ilong()) << i
        big2 = big2 + n
      end

      let mult1 = big1 * big2
      h.log("mult1=" + mult1.dump())
      let mult2 = big1.karatsuba_mul(big2)
      h.log("mult2=" + mult2.dump())
      h.assert_true(mult1 == mult2, "Karatsuba multiplication")
      h.log("big1 * big2 = " + big1.string() + " * " + big2.string())
      h.log("mult1=" + mult1.string())
      h.log("mult2=" + mult2.string())
    end


class iso _TestMPIntFastMultiplication is UnitTest
  """
  Tests on MPInt fast multiplication using FFT
  """

  fun name(): String =>
    "MPInt/fast_mult"

  fun apply(h: TestHelper) =>
    let rand = Rand()
    for i in Range(0, 1000) do
      let a = rand.ilong()
      let a' = MPInt.from_ilong(a)
      let b = rand.ilong()
      let b' = MPInt.from_ilong(b)

      h.assert_true((a' * b') == a'.fast_mul(b'), "Fast multiplication")
    end


class iso _TestMPIntFastMultiplicationLarge is UnitTest
  """
  Tests on MPInt fast multiplication using FFT with large numbers
  """

  fun name(): String =>
    "MPInt/fast_mul_large"

  fun apply(h: TestHelper) =>
    let rand = Rand()
    // 128 base-digits is 2^11 * 128 = 2^18 bits
    // Actually base is 2^16. 128 digits * 16 bits = 2048 bits.
    // Let's use 200 digits to be sure we are above max_size = 128
    let size: USize = 200

    for _ in Range(0, 100) do
      var big1 = MPInt.from_ilong(rand.ilong())
      for i in Range(1, size) do
        let n = MPInt.from_ilong(rand.ilong()).abs() << i
        big1 = big1 + n
      end

      var big2 = MPInt.from_ilong(rand.ilong())
      for i in Range(1, size) do
        let n = MPInt.from_ilong(rand.ilong()).abs() << i
        big2 = big2 + n
      end

      let mult1 = big1 * big2
      let mult2 = big1.fast_mul(big2)
      h.assert_true(mult1 == mult2, "Fast multiplication of large numbers")
    end
