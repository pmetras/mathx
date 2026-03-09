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
    h.assert_no_error({() ? => MPInt.from_string()? }, "Valid number 0")
    h.assert_no_error({() ? => MPInt.from_string("0")? }, "Valid number 0")
    h.assert_no_error({() ? => MPInt.from_string("+0")? }, "Valid number 0")
    h.assert_no_error({() ? => MPInt.from_string("-0")? }, "Valid number 0")
    h.assert_no_error({() ? => MPInt.from_string("-123_456_789")? }, "With separators")
    h.assert_no_error({() ? => MPInt.from_string("-123__456__789")? }, "With separators")
    h.assert_no_error({() ? => MPInt.from_string("-1_2_3__4_5_6__7_8_9")? }, "With separators")
    h.assert_no_error({() ? => MPInt.from_string("-123_456_789@0")? }, "Null exponent")
    h.assert_no_error({() ? => MPInt.from_string("-123_456_789@+0")? }, "Null exponent")
    h.assert_no_error({() ? => MPInt.from_string("+123_456_789@12")? }, "With exponent")
    h.assert_no_error({() ? => MPInt.from_string("-123_456_789e12")? }, "With exponent")
    h.assert_no_error({() ? => MPInt.from_string("+123_456_789E12")? }, "With exponent")
    h.assert_no_error({() ? => MPInt.from_string("+900d_F00d@+12", 16)? }, "Base 16")
    h.assert_no_error({() ? => MPInt.from_string("-1001010010010010111000101110100001110111000101010001101110@+50", 2)? }, "Base 2")
    h.assert_no_error({() ? => MPInt.from_string("-1001010010010010111000101110100001110111000101010001101110e+50", 2)? }, "Base 2")
    h.assert_no_error({() ? => MPInt.from_string("  -123_456_789")? }, "Leading spaces")
    h.assert_no_error({() ? => MPInt.from_string("-e50")? }, "Only exponent without mantissa")

    h.assert_error({() ? => MPInt.from_string("-")? }, "Only sign")
    h.assert_error({() ? => MPInt.from_string("+")? }, "Only sign")
    h.assert_error({() ? => MPInt.from_string("     ")? }, "Empty != 0")
    h.assert_error({() ? => MPInt.from_string("123_456_789 ")? }, "Trailing spaces")
    h.assert_error({() ? => MPInt.from_string("+-123")? }, "Multiple signs")
    h.assert_error({() ? => MPInt.from_string("-12.3")? }, "Float")
    h.assert_error({() ? => MPInt.from_string("12e")? }, "Empty exponent")
    h.assert_error({() ? => MPInt.from_string("12E+__")? }, "Empty exponent")
    h.assert_error({() ? => MPInt.from_string("+900d_F00dE+12", 16)? }, "Wrong exponent")
    h.assert_error({() ? => MPInt.from_string("+123_456_789E-12")? }, "Negative exponent")
    h.assert_error({() ? => MPInt.from_string("-900d_F00d@-12", 16)? }, "Negative exponent")

    try
      h.assert_true(MPInt.from_string("  +123_456_789")?.string() == "123456789")
    else
      h.fail("Wrong value")
    end

    try
      h.assert_true(MPInt.from_string("-123e6")?.string() == "-123000000")
    else
      h.fail("Wrong value")
    end

    try
      h.assert_true(MPInt.from_string()? == MPInt.from_ilong(), "Valid number 0")
      h.assert_true(MPInt.from_string("0")? == MPInt.from_ilong(0), "Valid number 0")
      // Pony parser does not support unary +
      h.assert_true(MPInt.from_string("+0")? == MPInt.from_ilong(0), "Valid number 0")
      h.assert_true(MPInt.from_string("-0")? == MPInt.from_ilong(-0), "Valid number 0")
      h.assert_true(MPInt.from_string("-123_456_789")? == MPInt.from_ilong(-123_456_789), "With separators")
      h.assert_true(MPInt.from_string("-123__456__789")? == MPInt.from_ilong(-123_456_789), "With separators")
      h.assert_true(MPInt.from_string("-1_2_3__4_5_6__7_8_9")? == MPInt.from_ilong(-1_2_3_4_5_6_7_8_9), "With separators")
      h.assert_true(MPInt.from_string("-123_456_789@0")? == MPInt.from_ilong(-123_456_789), "Null exponent")
      h.assert_true(MPInt.from_string("-123_456_789@+0")? == MPInt.from_ilong(-123_456_789), "Null exponent")
      h.assert_true(MPInt.from_string("+123_456_789@5")? == MPInt.from_ilong(123_456_789_00000), "With exponent")
      h.assert_true(MPInt.from_string("-123_456_789e3")? == MPInt.from_ilong(-123_456_789_000), "With exponent")
      h.assert_true(MPInt.from_string("+123_456_789E4")? == MPInt.from_ilong(123_456_789_0000), "With exponent")
      h.assert_true(MPInt.from_string("+900d_F00d@+2", 16)? == MPInt.from_ilong(0x900d_F00d * 16 * 16), "Base 16")
      h.assert_true(MPInt.from_string("-100101001001001011100010111010000111011100010101@+12", 2)? == MPInt.from_ilong(-669116797016952832), "Base 2")
      h.assert_true(MPInt.from_string("-100101001001001011100010111010000111011100010101e+12", 2)? == MPInt.from_ilong(-669116797016952832), "Base 2")
      h.assert_true(MPInt.from_string("  -123_456_789")? == MPInt.from_ilong(-123_456_789), "Leading spaces")
      h.assert_true(MPInt.from_string("-1e10")? == MPInt.from_ilong(-10000000000), "Exponent")
      h.assert_true(MPInt.from_string("18795432178911876e21")?.string() == "18795432178911876000000000000000000000", "String representation")
      h.assert_true(MPInt.from_string("-37892648294")? == MPInt.from_ilong(-37892648294))
    else
      h.fail("Equalities")
    end

    try
      let rand = Rand()
      for i in Range(0, 1000) do
        let n = rand.ilong()
        let s: String = n.string()
        h.assert_true(MPInt.from_string(s)? == MPInt.from_ilong(n), "Random values")
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
      h.assert_true((minv == MPInt.from_string("-8000_0000", 16)?) or (minv == MPInt.from_string("-8000_0000_0000_0000", 16)?) , "Min value")
      h.assert_true((maxv == MPInt.from_string("7fff_ffff", 16)?) or (maxv == MPInt.from_string("7fff_ffff_ffff_ffff", 16)?) , "Max value")
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
      let zerop = MPInt.from_string("+0")?
      let zerom = MPInt.from_string("-0")?
      h.assert_true(zero.is_zero(), "Zero")
      h.assert_true(zerop.is_zero(), "Zero")
      h.assert_true(zerom.is_zero(), "Zero")

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
        let x' = MPInt.from_string(x)?
        let y' = MPInt.from_string(y)?

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
      h.assert_true(a == a'.ilong(), "Conversion to ILong")
    end

    try
      let b = MPInt.from_string("89334718795878195613438179751988712")?
      // ilong() now returns 0 on overflow instead of erroring, 
      // or we could make it partial. For now, let's just fix test syntax.
      let b' = b.ilong()

      let c = MPInt.from_string("-19815798458789751981278742378912370")?
      let c' = c.ilong()
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
    h.assert_true(one.digit_shl(1) == MPInt.from_ilong(65536), "Shift left by 1")
    h.assert_true(one.digit_shl(2) == MPInt.from_ilong(4294967296), "Shift left by 2")
    h.assert_true(one.digit_shl(3) == MPInt.from_ilong(281474976710656), "Shift left by 3")
    h.assert_true(one.digit_shl(10).digit_shr(10) == one, "Shift left then right")

    let one3 = MPInt.from_ilong(281474976710656)
    let one2 = MPInt.from_ilong(4294967296)
    let one1 = MPInt.from_ilong(65536)
    let zero = MPInt.from_ilong(0)
    h.assert_true(one3.digit_shr(3) == one, "Shift right by 3")
    h.assert_true(one2.digit_shr(2) == one, "Shift right by 2")
    h.assert_true(one1.digit_shr(1) == one, "Shift right by 1")
    h.assert_true(one.digit_shr(1) == zero, "Shift right to 0")


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
      h.assert_true((MPInt.from_string("+123456")? + MPInt.from_string("+1234567890123456789")?) == MPInt.from_string("+1234567890123580245")?)
      h.assert_true((MPInt.from_string("-123456")? + MPInt.from_string("-1234567890123456789")?) == MPInt.from_string("-1234567890123580245")?)
      h.assert_true((MPInt.from_string("+123456")? + MPInt.from_string("-1234567890123456789")?) == MPInt.from_string("-1234567890123333333")?)
      h.assert_true((MPInt.from_string("-123456")? + MPInt.from_string("+1234567890123456789")?) == MPInt.from_string("+1234567890123333333")?)
    else
      h.fail("Addition of big and small")
    end

    try
      h.assert_true((MPInt.from_string("19286568262396999950202928044648128737703949436999793610613130182105922367448423453159208008561154643649004993410245621815263683558649945634758267696270039303164684506235887452943484121027471056665201277603534728479727382153288066446108997401020223592515566263174704019778189859650264257813887322843943022301697122037047349594826189769800767532590981681753960161163017803752421139870163553766644589335272831995013116481436362323269375288953950813601638251442829597903781688297881241170668619904757069020265297394635546155851093309363569785849294570445068455981948468662793917426384682117599755965209933018194164415994225160507289598509847921")? + MPInt.from_string("-89397816278271230459454110395996026663027374858959155341995445324393903731711371299662676953068781699910566197767070398021792037309382471935003522844234398194213451613830240178084330862239920913616552830617626796686980142144475041619332654915691024934414728326254713483694918659373308283669427344580771115999103149093958866604351791680623545304481448200580679477784543921552472823363295771446882280964663265344973021812416638734557654215643369337559645949925843413734820414989440796460578313739641380350952454004469826268119514126221846674571499744339378765148723031976486968205020823203122321431396802517326793155991667238272944727233749")?) == MPInt.from_string("19197170446118728719743473934252132711040922062140834455271134736781528463716712081859545331608085861949094427212478551417241891521340563162823264173425804904970471054622057212765399790165231135751584724772917101683040402011143591404489664746104532567581151534848449306294494940990890949530217895499362251185698018887953390728221837978120143987286500233553379481685233259830868667046800257995197707054308168729668143459623945684534817634738307444264078605492903754490046867882891800374208041591017427639914344940631076329582973795237347939174723070700729077216799745630817430458179661294396633643778536215676837622838233493269016653782614172")?, "Addition big numbers +-")
      h.assert_true((MPInt.from_string("19286568262396999950202928044648128737703949436999793610613130182105922367448423453159208008561154643649004993410245621815263683558649945634758267696270039303164684506235887452943484121027471056665201277603534728479727382153288066446108997401020223592515566263174704019778189859650264257813887322843943022301697122037047349594826189769800767532590981681753960161163017803752421139870163553766644589335272831995013116481436362323269375288953950813601638251442829597903781688297881241170668619904757069020265297394635546155851093309363569785849294570445068455981948468662793917426384682117599755965209933018194164415994225160507289598509847921")? + MPInt.from_string("89397816278271230459454110395996026663027374858959155341995445324393903731711371299662676953068781699910566197767070398021792037309382471935003522844234398194213451613830240178084330862239920913616552830617626796686980142144475041619332654915691024934414728326254713483694918659373308283669427344580771115999103149093958866604351791680623545304481448200580679477784543921552472823363295771446882280964663265344973021812416638734557654215643369337559645949925843413734820414989440796460578313739641380350952454004469826268119514126221846674571499744339378765148723031976486968205020823203122321431396802517326793155991667238272944727233749")?) == MPInt.from_string("19375966078675271180662382155044124764366976811858752765955125627430316271180134824458870685514223425348915559608012692213285475595959328106693271219114273701358897957849717693121568451889710977578817830434152355276414362295432541487728330055935914617449980991500958733261884778309637566097556750188523793417696225186141308461430541561481391077895463129954540840640802347673973612693526849538091471616237495260358089503248778962003932943169594182939197897392755441317516508712870681967129198218496710400616249848640015982119212823489791632523866070189407834747097191694770404394589702940802878286641329820711491209150216827745562543237081670")?, "Addition big numbers ++")
      h.assert_true((MPInt.from_string("-19286568262396999950202928044648128737703949436999793610613130182105922367448423453159208008561154643649004993410245621815263683558649945634758267696270039303164684506235887452943484121027471056665201277603534728479727382153288066446108997401020223592515566263174704019778189859650264257813887322843943022301697122037047349594826189769800767532590981681753960161163017803752421139870163553766644589335272831995013116481436362323269375288953950813601638251442829597903781688297881241170668619904757069020265297394635546155851093309363569785849294570445068455981948468662793917426384682117599755965209933018194164415994225160507289598509847921")? + MPInt.from_string("-89397816278271230459454110395996026663027374858959155341995445324393903731711371299662676953068781699910566197767070398021792037309382471935003522844234398194213451613830240178084330862239920913616552830617626796686980142144475041619332654915691024934414728326254713483694918659373308283669427344580771115999103149093958866604351791680623545304481448200580679477784543921552472823363295771446882280964663265344973021812416638734557654215643369337559645949925843413734820414989440796460578313739641380350952454004469826268119514126221846674571499744339378765148723031976486968205020823203122321431396802517326793155991667238272944727233749")?) == MPInt.from_string("-19375966078675271180662382155044124764366976811858752765955125627430316271180134824458870685514223425348915559608012692213285475595959328106693271219114273701358897957849717693121568451889710977578817830434152355276414362295432541487728330055935914617449980991500958733261884778309637566097556750188523793417696225186141308461430541561481391077895463129954540840640802347673973612693526849538091471616237495260358089503248778962003932943169594182939197897392755441317516508712870681967129198218496710400616249848640015982119212823489791632523866070189407834747097191694770404394589702940802878286641329820711491209150216827745562543237081670")?, "Addition big numbers --")
      h.assert_true((MPInt.from_string("-19286568262396999950202928044648128737703949436999793610613130182105922367448423453159208008561154643649004993410245621815263683558649945634758267696270039303164684506235887452943484121027471056665201277603534728479727382153288066446108997401020223592515566263174704019778189859650264257813887322843943022301697122037047349594826189769800767532590981681753960161163017803752421139870163553766644589335272831995013116481436362323269375288953950813601638251442829597903781688297881241170668619904757069020265297394635546155851093309363569785849294570445068455981948468662793917426384682117599755965209933018194164415994225160507289598509847921")? + MPInt.from_string("89397816278271230459454110395996026663027374858959155341995445324393903731711371299662676953068781699910566197767070398021792037309382471935003522844234398194213451613830240178084330862239920913616552830617626796686980142144475041619332654915691024934414728326254713483694918659373308283669427344580771115999103149093958866604351791680623545304481448200580679477784543921552472823363295771446882280964663265344973021812416638734557654215643369337559645949925843413734820414989440796460578313739641380350952454004469826268119514126221846674571499744339378765148723031976486968205020823203122321431396802517326793155991667238272944727233749")?) == MPInt.from_string("-19197170446118728719743473934252132711040922062140834455271134736781528463716712081859545331608085861949094427212478551417241891521340563162823264173425804904970471054622057212765399790165231135751584724772917101683040402011143591404489664746104532567581151534848449306294494940990890949530217895499362251185698018887953390728221837978120143987286500233553379481685233259830868667046800257995197707054308168729668143459623945684534817634738307444264078605492903754490046867882891800374208041591017427639914344940631076329582973795237347939174723070700729077216799745630817430458179661294396633643778536215676837622838233493269016653782614172")?, "Addition big numbers -+")
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
      h.assert_true((MPInt.from_string("+654321")? - MPInt.from_string("+9876543210987654321")?) == MPInt.from_string("-9876543210987000000")?)
      h.assert_true((MPInt.from_string("-654321")? - MPInt.from_string("-9876543210987654321")?) == MPInt.from_string("+9876543210987000000")?)
      h.assert_true((MPInt.from_string("+654321")? - MPInt.from_string("-9876543210987654321")?) == MPInt.from_string("+9876543210988308642")?)
      h.assert_true((MPInt.from_string("-654321")? - MPInt.from_string("+9876543210987654321")?) == MPInt.from_string("-9876543210988308642")?)
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
        let n = MPInt.from_ilong(rand.ilong()).digit_shl(i)
        big1 = big1 + n
      end

      var big2 = MPInt.from_ilong(rand.ilong())
      for i in Range(1, size2) do
        let n = MPInt.from_ilong(rand.ilong()).digit_shl(i)
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
        let n = MPInt.from_ilong(rand.ilong()).abs().digit_shl(i)
        big1 = big1 + n
      end

      var big2 = MPInt.from_ilong(rand.ilong())
      for i in Range(1, size) do
        let n = MPInt.from_ilong(rand.ilong()).abs().digit_shl(i)
        big2 = big2 + n
      end

      let mult1 = big1 * big2
      let mult2 = big1.fast_mul(big2)
      h.assert_true(mult1 == mult2, "Fast multiplication of large numbers")
    end


class iso _TestMPIntMultiplicationComparison is UnitTest
  """
  Compare different multiplication algorithms (Schoolbook, Karatsuba, FFT)
  to ensure they all produce the same results for various sizes.
  """

  fun name(): String =>
    "MPInt/multiplication_comparison"

  fun apply(h: TestHelper) =>
    let rand = Rand()
    
    // Test sizes: 
    // - Small: 1-10 digits (Schoolbook)
    // - Medium: 150 digits (Triggers Karatsuba)
    // - Large: 500 digits (FFT range)
    let sizes: Array[USize] = [1; 10; 150; 500]

    for size in sizes.values() do
      for _ in Range(0, 10) do
        let a = _random_mpint(rand, size)
        let b = _random_mpint(rand, size)

        let res_school = a.mul(b)
        let res_karatsuba = a.karatsuba_mul(b)
        let res_fft = a.fast_mul(b)

        h.assert_true(res_school == res_karatsuba, 
          "Schoolbook vs Karatsuba mismatch at size " + size.string())
        h.assert_true(res_school == res_fft, 
          "Schoolbook vs FFT mismatch at size " + size.string())
      end
    end

  fun _random_mpint(rand: Rand, size: USize): MPInt =>
    var res = MPInt.from_ilong(rand.ilong())
    for i in Range(1, size) do
      res = res + (MPInt.from_ilong(rand.ilong()).abs().digit_shl(i))
    end
    if (rand.next() % 2) == 0 then
      res.neg()
    else
      res
    end


class iso _TestMPIntMultiplicationEdgeCases is UnitTest
  """
  Explicitly test multiplication edge cases: 0, 1, -1, squaring.
  """

  fun name(): String =>
    "MPInt/multiplication_edge_cases"

  fun apply(h: TestHelper) =>
    let zero = MPInt.from_ilong(0)
    let one = MPInt.from_ilong(1)
    let m_one = MPInt.from_ilong(-1)
    let two = MPInt.from_ilong(2)
    
    let rand = Rand()
    let large = _random_mpint(rand, 200)

    // Tests for mul
    h.assert_true(zero.mul(large).is_zero())
    h.assert_true(one.mul(large) == large)
    h.assert_true(m_one.mul(large) == large.neg())

    // Tests for karatsuba_mul
    h.assert_true(zero.karatsuba_mul(large).is_zero())
    h.assert_true(one.karatsuba_mul(large) == large)
    h.assert_true(m_one.karatsuba_mul(large) == large.neg())

    // Tests for fast_mul
    h.assert_true(zero.fast_mul(large).is_zero())
    h.assert_true(one.fast_mul(large) == large)
    h.assert_true(m_one.fast_mul(large) == large.neg())

    // Symmetric cases
    h.assert_true(large.mul(zero).is_zero())
    h.assert_true(large.karatsuba_mul(zero).is_zero())
    h.assert_true(large.fast_mul(zero).is_zero())

    h.assert_true(large.mul(one) == large)
    h.assert_true(large.karatsuba_mul(one) == large)
    h.assert_true(large.fast_mul(one) == large)

    h.assert_true(large.mul(m_one) == large.neg())
    h.assert_true(large.karatsuba_mul(m_one) == large.neg())
    h.assert_true(large.fast_mul(m_one) == large.neg())

    // 0 * 0
    h.assert_true(zero.mul(zero).is_zero())
    h.assert_true(zero.fast_mul(zero).is_zero())

    // 1 * -1
    h.assert_true(one.mul(m_one) == m_one)
    h.assert_true(one.fast_mul(m_one) == m_one)

    // Squaring
    h.assert_true(two.mul(two) == MPInt.from_ilong(4))
    h.assert_true(large.mul(large) == large.fast_mul(large))

  fun _random_mpint(rand: Rand, size: USize): MPInt =>
    var res = MPInt.from_ilong(rand.ilong())
    for i in Range(1, size) do
      res = res + (MPInt.from_ilong(rand.ilong()).abs().digit_shl(i))
    end
    if (rand.next() % 2) == 0 then
      res.neg()
    else
      res
    end


class iso _TestMPIntTrait is UnitTest
  """
  Test MPInt SignedInteger and UnsignedInteger trait compliance.
  """

  fun name(): String =>
    "MPInt/trait"

  fun apply(h: TestHelper) =>
    try
      let a = MPInt.from_string("100")?
      let b = MPInt.from_string("30")?
      
      // Basic operators
      h.assert_true((a + b) == MPInt.from_string("130")?, "+")
      h.assert_true((a - b) == MPInt.from_string("70")?, "-")
      h.assert_true((a * b) == MPInt.from_string("3000")?, "*")
      h.assert_true((a / b) == MPInt.from_string("3")?, "/")
      h.assert_true((a % b) == MPInt.from_string("10")?, "%")
      
      // Traits requirements
      h.assert_true(a.abs() == a, "abs positive")
      h.assert_true(MPInt.from_string("-10")?.abs() == MPInt.from_string("10")?, "abs negative")
      
      h.assert_true(a.min(b) == b, "min")
      h.assert_true(a.max(b) == a, "max")
      
      // Conversions
      h.assert_true(a.i64() == 100, "i64")
      h.assert_true(a.f64() == 100.0, "f64")
      
      // String
      h.assert_true(a.string() == "100", "string")
      h.assert_true(MPInt.from_string("-42")?.string() == "-42", "negative string")

      // Bitwidth
      // 100 is 1100100 in binary (7 bits) + 1 sign bit = 8 bits
      h.assert_true(a.bitwidth() == MPInt.from_ilong(8), "bitwidth")

    else
      h.fail("Error in trait tests")
    end


class iso _TestMPIntShift is UnitTest
  """
  Tests for bit_shl and bit_shr on MPInt.
  """

  fun name(): String =>
    "MPInt/shift"

  fun apply(h: TestHelper) =>
    let zero = MPInt.from_ilong(0)
    let one  = MPInt.from_ilong(1)

    // bit_shl basic
    h.assert_true(one.bit_shl(MPInt.from_ilong(0))  == one,                  "1 << 0 = 1")
    h.assert_true(one.bit_shl(MPInt.from_ilong(1))  == MPInt.from_ilong(2),  "1 << 1 = 2")
    h.assert_true(one.bit_shl(MPInt.from_ilong(16)) == MPInt.from_ilong(65536),  "1 << 16 = 65536")
    h.assert_true(one.bit_shl(MPInt.from_ilong(17)) == MPInt.from_ilong(131072), "1 << 17 = 131072")

    // Cross-digit boundary: 0xABCD << 4 = 0xABCD0
    h.assert_true(
      MPInt.from_ilong(0xABCD).bit_shl(MPInt.from_ilong(4)) == MPInt.from_ilong(0xABCD0),
      "0xABCD << 4 = 0xABCD0")

    // Round-trip x.bit_shl(n).bit_shr(n) == x
    let shift = MPInt.from_ilong(7)
    let rt_vals: Array[MPInt val] val = [
      MPInt.from_ilong(1); MPInt.from_ilong(15); MPInt.from_ilong(16)
      MPInt.from_ilong(17); MPInt.from_ilong(32); MPInt.from_ilong(100)
      MPInt.from_ilong(65535); MPInt.from_ilong(65536); MPInt.from_ilong(131072)]
    for x_rt in rt_vals.values() do
      h.assert_true(x_rt.bit_shl(shift).bit_shr(shift) == x_rt, "round-trip")
    end

    // bit_shr truncation
    h.assert_true(MPInt.from_ilong(5).bit_shr(MPInt.from_ilong(1)) == MPInt.from_ilong(2), "5 >> 1 = 2")
    h.assert_true(MPInt.from_ilong(1).bit_shr(MPInt.from_ilong(1)) == zero,               "1 >> 1 = 0")

    // Cross-digit shr: 0xABCD_0000 >> 4 = 0xABCD_000
    h.assert_true(
      MPInt.from_ilong(0xABCD0000).bit_shr(MPInt.from_ilong(4)) == MPInt.from_ilong(0xABCD000),
      "0xABCD0000 >> 4 = 0xABCD000")

    // Negative numbers: sign is preserved
    h.assert_true(
      MPInt.from_ilong(-5).bit_shl(MPInt.from_ilong(3)) == MPInt.from_ilong(-40),
      "(-5) << 3 = -40")
    h.assert_true(
      MPInt.from_ilong(-40).bit_shr(MPInt.from_ilong(3)) == MPInt.from_ilong(-5),
      "(-40) >> 3 = -5")

    // y < 0 → unchanged
    h.assert_true(
      MPInt.from_ilong(42).bit_shl(MPInt.from_ilong(-1)) == MPInt.from_ilong(42),
      "shl by -1 = noop")
    h.assert_true(
      MPInt.from_ilong(42).bit_shr(MPInt.from_ilong(-1)) == MPInt.from_ilong(42),
      "shr by -1 = noop")

    // shr by >= digit count → 0
    h.assert_true(one.bit_shr(MPInt.from_ilong(1000)) == zero, "1 >> 1000 = 0")


class iso _TestMPIntBitwise is UnitTest
  """
  Tests for op_and, op_or, op_xor, op_not on MPInt (two's complement semantics).
  """

  fun name(): String =>
    "MPInt/bitwise"

  fun apply(h: TestHelper) =>
    let zero = MPInt.from_ilong(0)
    let one  = MPInt.from_ilong(1)
    let three = MPInt.from_ilong(3)
    let five  = MPInt.from_ilong(5)
    let seven = MPInt.from_ilong(7)
    let neg1  = MPInt.from_ilong(-1)
    let neg3  = MPInt.from_ilong(-3)
    let neg5  = MPInt.from_ilong(-5)

    // Basic ops on positive numbers
    h.assert_true(five.op_and(three) == one,                  "5 & 3 = 1")
    h.assert_true(five.op_or(three)  == seven,                "5 | 3 = 7")
    h.assert_true(five.op_xor(three) == MPInt.from_ilong(6), "5 ^ 3 = 6")

    // op_not: ~x = -(x + 1)
    h.assert_true(zero.op_not()  == neg1,                "~0 = -1")
    h.assert_true(neg1.op_not()  == zero,                "~(-1) = 0")
    h.assert_true(five.op_not()  == MPInt.from_ilong(-6), "~5 = -6")

    // Negative operands
    h.assert_true(neg1.op_and(seven) == seven, "(-1) & 7 = 7")
    h.assert_true(neg1.op_or(zero)   == neg1,  "(-1) | 0 = -1")
    h.assert_true(neg3.op_xor(neg3)  == zero,  "(-3) ^ (-3) = 0")

    // Mixed sign: -5 = ...11111011 in two's complement
    h.assert_true(neg5.op_and(three) == three,               "(-5) & 3 = 3")
    h.assert_true(neg5.op_or(three)  == neg5,                "(-5) | 3 = -5")
    h.assert_true(neg5.op_xor(three) == MPInt.from_ilong(-8), "(-5) ^ 3 = -8")

    // Identity properties
    h.assert_true(five.op_and(zero) == zero, "x & 0 = 0")
    h.assert_true(five.op_or(zero)  == five, "x | 0 = x")
    h.assert_true(five.op_xor(five) == zero, "x ^ x = 0")

    // Larger values round-trip: (x op_and (op_not x)) = 0  (0xDEAD_BEEF = 3735928559)
    let big = MPInt.from_ilong(3735928559)
    h.assert_true(big.op_and(big.op_not()) == zero, "x & ~x = 0")
    h.assert_true(big.op_or(big.op_not())  == neg1,  "x | ~x = -1")


class iso _TestMPIntSign is UnitTest
  """
  Tests for is_negative, is_positive, signum.
  """

  fun name(): String =>
    "MPInt/sign"

  fun apply(h: TestHelper) =>
    let zero = MPInt.from_ilong(0)
    let pos  = MPInt.from_ilong(42)
    let neg  = MPInt.from_ilong(-7)

    h.assert_false(zero.is_negative(), "0 is not negative")
    h.assert_false(zero.is_positive(), "0 is not positive")
    h.assert_true(pos.is_positive(),   "42 is positive")
    h.assert_false(pos.is_negative(),  "42 is not negative")
    h.assert_true(neg.is_negative(),   "-7 is negative")
    h.assert_false(neg.is_positive(),  "-7 is not positive")

    h.assert_true(zero.compare(MPInt.from_ilong(0))  is Equal,   "0.compare(0) = Equal")
    h.assert_true(zero.compare(pos)                  is Less,    "0.compare(42) = Less")
    h.assert_true(pos.compare(zero)                  is Greater, "42.compare(0) = Greater")
    h.assert_true(neg.compare(zero)                  is Less,    "-7.compare(0) = Less")
    h.assert_true(pos.compare(neg)                   is Greater, "42.compare(-7) = Greater")
    h.assert_true(neg.compare(pos)                   is Less,    "-7.compare(42) = Less")


class iso _TestMPIntBitAccess is UnitTest
  """
  Tests for bit_get, bit_set, bit_clear, bit_flip.
  """

  fun name(): String =>
    "MPInt/bit_access"

  fun apply(h: TestHelper) =>
    let zero = MPInt.from_ilong(0)
    // 0b1010_1100 = 172
    let n = MPInt.from_ilong(172)

    // bit_get
    h.assert_false(n.bit_get(0), "bit 0 of 172 = 0")
    h.assert_false(n.bit_get(1), "bit 1 of 172 = 0")
    h.assert_true(n.bit_get(2),  "bit 2 of 172 = 1")
    h.assert_true(n.bit_get(3),  "bit 3 of 172 = 1")
    h.assert_false(n.bit_get(4), "bit 4 of 172 = 0")
    h.assert_true(n.bit_get(5),  "bit 5 of 172 = 1")
    h.assert_false(n.bit_get(6), "bit 6 of 172 = 0")
    h.assert_true(n.bit_get(7),  "bit 7 of 172 = 1")
    h.assert_false(n.bit_get(8), "bit 8 of 172 = 0 (above MSB)")

    // bit_set: 172 | (1<<0) = 173
    h.assert_true(n.bit_set(0) == MPInt.from_ilong(173), "172 set bit 0 = 173")
    // bit_set on already-set bit: idempotent
    h.assert_true(n.bit_set(2) == n, "172 set bit 2 = 172")
    // bit_set extending into new digit (bit 16)
    let with_bit16 = n.bit_set(16)
    h.assert_true(with_bit16.bit_get(16), "bit 16 set")
    h.assert_true(with_bit16 == MPInt.from_ilong(65536 + 172), "172 set bit 16")

    // bit_clear: 172 & ~(1<<2) = 168
    h.assert_true(n.bit_clear(2) == MPInt.from_ilong(168), "172 clear bit 2 = 168")
    // bit_clear on already-clear bit: no-op
    h.assert_true(n.bit_clear(0) == n, "172 clear bit 0 = 172")

    // bit_flip
    h.assert_true(n.bit_flip(0) == MPInt.from_ilong(173), "172 flip bit 0 = 173")
    h.assert_true(n.bit_flip(2) == MPInt.from_ilong(168), "172 flip bit 2 = 168")

    // Round-trip: set then clear restores original
    h.assert_true(n.bit_set(0).bit_clear(0) == n, "set/clear round-trip")

    // Sign is preserved
    let neg_n = MPInt.from_ilong(-172)
    h.assert_true(neg_n.bit_get(2), "bit_get on negative uses absolute value")
    h.assert_true(neg_n.bit_set(0) == MPInt.from_ilong(-173), "bit_set preserves sign")


class iso _TestMPIntPow is UnitTest
  """
  Tests for pow, isqrt, gcd, pow_mod.
  """

  fun name(): String =>
    "MPInt/pow"

  fun apply(h: TestHelper) =>
    let zero = MPInt.from_ilong(0)
    let one  = MPInt.from_ilong(1)
    let two  = MPInt.from_ilong(2)

    // pow: basic cases
    h.assert_true(two.pow(MPInt.from_ilong(0))  == one,                 "2^0 = 1")
    h.assert_true(two.pow(MPInt.from_ilong(1))  == two,                 "2^1 = 2")
    h.assert_true(two.pow(MPInt.from_ilong(10)) == MPInt.from_ilong(1024), "2^10 = 1024")
    h.assert_true(MPInt.from_ilong(3).pow(MPInt.from_ilong(4)) == MPInt.from_ilong(81), "3^4 = 81")
    h.assert_true(zero.pow(MPInt.from_ilong(5)) == zero,               "0^5 = 0")
    h.assert_true(MPInt.from_ilong(5).pow(MPInt.from_ilong(-1)) == zero, "5^(-1) = 0")
    // Negative base: (-2)^3 = -8, (-2)^4 = 16
    h.assert_true(MPInt.from_ilong(-2).pow(MPInt.from_ilong(3)) == MPInt.from_ilong(-8), "(-2)^3 = -8")
    h.assert_true(MPInt.from_ilong(-2).pow(MPInt.from_ilong(4)) == MPInt.from_ilong(16), "(-2)^4 = 16")

    // isqrt
    h.assert_true(zero.isqrt() == zero,                   "isqrt(0) = 0")
    h.assert_true(one.isqrt()  == one,                    "isqrt(1) = 1")
    h.assert_true(MPInt.from_ilong(4).isqrt()   == two,   "isqrt(4) = 2")
    h.assert_true(MPInt.from_ilong(9).isqrt()   == MPInt.from_ilong(3),  "isqrt(9) = 3")
    h.assert_true(MPInt.from_ilong(10).isqrt()  == MPInt.from_ilong(3),  "isqrt(10) = 3")
    h.assert_true(MPInt.from_ilong(100).isqrt() == MPInt.from_ilong(10), "isqrt(100) = 10")
    // Large perfect square: 1000000 = 1000^2
    h.assert_true(MPInt.from_ilong(1000000).isqrt() == MPInt.from_ilong(1000), "isqrt(10^6) = 1000")
    // Negative → 0
    h.assert_true(MPInt.from_ilong(-9).isqrt() == zero, "isqrt(-9) = 0")

    // gcd
    h.assert_true(MPInt.from_ilong(12).gcd(MPInt.from_ilong(8))  == MPInt.from_ilong(4), "gcd(12,8) = 4")
    h.assert_true(MPInt.from_ilong(7).gcd(MPInt.from_ilong(13))  == one, "gcd(7,13) = 1 (coprime)")
    h.assert_true(MPInt.from_ilong(0).gcd(MPInt.from_ilong(5))   == MPInt.from_ilong(5), "gcd(0,5) = 5")
    h.assert_true(MPInt.from_ilong(5).gcd(zero)   == MPInt.from_ilong(5), "gcd(5,0) = 5")
    // gcd is symmetric
    h.assert_true(MPInt.from_ilong(48).gcd(MPInt.from_ilong(36)) ==
                  MPInt.from_ilong(36).gcd(MPInt.from_ilong(48)), "gcd symmetric")

    // pow_mod
    // 2^10 mod 1000 = 24
    h.assert_true(
      two.pow_mod(MPInt.from_ilong(10), MPInt.from_ilong(1000)) == MPInt.from_ilong(24),
      "2^10 mod 1000 = 24")
    // Fermat's little theorem: a^(p-1) ≡ 1 (mod p) for prime p, gcd(a,p)=1
    // 3^6 mod 7 = 1
    h.assert_true(
      MPInt.from_ilong(3).pow_mod(MPInt.from_ilong(6), MPInt.from_ilong(7)) == one,
      "3^6 mod 7 = 1 (Fermat)")
    // pow_mod(0, m) = 1 for m > 1
    h.assert_true(
      two.pow_mod(zero, MPInt.from_ilong(100)) == one,
      "2^0 mod 100 = 1")
    // pow_mod(n, 1) = 0
    h.assert_true(
      MPInt.from_ilong(999).pow_mod(MPInt.from_ilong(999), one) == zero,
      "x^n mod 1 = 0")


class iso _TestMPIntPredicates is UnitTest
  """
  Tests for neg, is_one, is_minus_one, abs(0), min/max edge cases.
  """

  fun name(): String =>
    "MPInt/predicates"

  fun apply(h: TestHelper) =>
    let zero = MPInt.from_ilong(0)
    let one  = MPInt.from_ilong(1)
    let m1   = MPInt.from_ilong(-1)

    // neg
    h.assert_true(MPInt.from_ilong(5).neg()  == MPInt.from_ilong(-5), "neg(5) = -5")
    h.assert_true(MPInt.from_ilong(-5).neg() == MPInt.from_ilong(5),  "neg(-5) = 5")
    h.assert_true(zero.neg().is_zero(),                                "neg(0) = 0")

    // is_one
    h.assert_true(one.is_one(),                         "1.is_one()")
    h.assert_false(zero.is_one(),                       "0 not is_one")
    h.assert_false(MPInt.from_ilong(2).is_one(),        "2 not is_one")
    h.assert_false(m1.is_one(),                         "-1 not is_one")

    // is_minus_one
    h.assert_true(m1.is_minus_one(),                    "-1.is_minus_one()")
    h.assert_false(zero.is_minus_one(),                 "0 not is_minus_one")
    h.assert_false(one.is_minus_one(),                  "1 not is_minus_one")
    h.assert_false(MPInt.from_ilong(-2).is_minus_one(), "-2 not is_minus_one")

    // abs(0) == 0
    h.assert_true(zero.abs().is_zero(), "abs(0) = 0")

    // min/max with negative values
    let neg3 = MPInt.from_ilong(-3)
    let pos5 = MPInt.from_ilong(5)
    h.assert_true(neg3.min(pos5) == neg3, "min(-3, 5) = -3")
    h.assert_true(pos5.min(neg3) == neg3, "min(5, -3) = -3")
    h.assert_true(neg3.max(pos5) == pos5, "max(-3, 5) = 5")
    h.assert_true(pos5.max(neg3) == pos5, "max(5, -3) = 5")

    // min/max with equal values
    h.assert_true(pos5.min(pos5) == pos5, "min(x, x) = x")
    h.assert_true(neg3.max(neg3) == neg3, "max(x, x) = x")

    // min/max with both negative
    let neg7 = MPInt.from_ilong(-7)
    h.assert_true(neg3.min(neg7) == neg7, "min(-3, -7) = -7")
    h.assert_true(neg3.max(neg7) == neg3, "max(-3, -7) = -3")


class iso _TestMPIntHash is UnitTest
  """
  Tests for hash: equal values must hash equally; hash must be stable.
  """

  fun name(): String =>
    "MPInt/hash"

  fun apply(h: TestHelper) =>
    let zero  = MPInt.from_ilong(0)
    let one   = MPInt.from_ilong(1)
    let m_one = MPInt.from_ilong(-1)

    // Stability: same value → same hash
    h.assert_true(zero.hash() == MPInt.from_ilong(0).hash(), "hash stable: 0")
    h.assert_true(one.hash()  == MPInt.from_ilong(1).hash(), "hash stable: 1")

    // Equal values have equal hashes (including signed-zero equivalence)
    let rand = Rand(99)
    for _ in Range(0, 500) do
      let v = rand.ilong()
      let a = MPInt.from_ilong(v)
      let b = MPInt.from_ilong(v)
      h.assert_true(a.hash() == b.hash(), "equal → same hash")
    end

    // Different signs produce different hashes (0 vs -1)
    h.assert_false(zero.hash() == m_one.hash(), "hash(0) != hash(-1)")
    h.assert_false(one.hash()  == m_one.hash(), "hash(1) != hash(-1)")


class iso _TestMPIntCarry is UnitTest
  """
  Tests for addc, subc, mulc, divc, remc, fldc, modc.
  For MPInt, overflow never occurs; the carry/borrow flag is always false
  except for division by zero where it is true.
  """

  fun name(): String =>
    "MPInt/carry"

  fun apply(h: TestHelper) =>
    let five = MPInt.from_ilong(5)
    let three = MPInt.from_ilong(3)
    let ten  = MPInt.from_ilong(10)
    let zero = MPInt.from_ilong(0)

    (let r_add, let c_add) = five.addc(three)
    h.assert_true(r_add == MPInt.from_ilong(8), "addc: result")
    h.assert_false(c_add,                       "addc: no carry")

    (let r_sub, let c_sub) = five.subc(three)
    h.assert_true(r_sub == MPInt.from_ilong(2), "subc: result")
    h.assert_false(c_sub,                       "subc: no carry")

    (let r_mul, let c_mul) = five.mulc(three)
    h.assert_true(r_mul == MPInt.from_ilong(15), "mulc: result")
    h.assert_false(c_mul,                        "mulc: no carry")

    (let r_div, let c_div) = ten.divc(three)
    h.assert_true(r_div == MPInt.from_ilong(3), "divc: result")
    h.assert_false(c_div,                       "divc: no overflow")

    (let r_div0, let c_div0) = ten.divc(zero)
    h.assert_true(c_div0, "divc by zero: overflow flag set")

    (let r_rem, let c_rem) = ten.remc(three)
    h.assert_true(r_rem == MPInt.from_ilong(1), "remc: result")
    h.assert_false(c_rem,                       "remc: no overflow")

    (let r_rem0, let c_rem0) = ten.remc(zero)
    h.assert_true(c_rem0, "remc by zero: overflow flag set")

    (let r_fld, let c_fld) = ten.fldc(three)
    h.assert_true(r_fld == MPInt.from_ilong(3), "fldc: result")
    h.assert_false(c_fld,                       "fldc: no overflow")

    (let r_fld0, let c_fld0) = ten.fldc(zero)
    h.assert_true(c_fld0, "fldc by zero: overflow flag set")

    (let r_mod, let c_mod) = ten.modc(three)
    h.assert_true(r_mod == MPInt.from_ilong(1), "modc: result")
    h.assert_false(c_mod,                       "modc: no overflow")

    (let r_mod0, let c_mod0) = ten.modc(zero)
    h.assert_true(c_mod0, "modc by zero: overflow flag set")


class iso _TestMPIntPartialOps is UnitTest
  """
  Tests for add_partial, sub_partial, mul_partial, div_partial, rem_partial,
  divrem_partial, fld_partial, mod_partial.
  add/sub/mul never error; div/rem/fld/mod error on zero divisor.
  """

  fun name(): String =>
    "MPInt/partial_ops"

  fun apply(h: TestHelper) =>
    let a = MPInt.from_ilong(10)
    let b = MPInt.from_ilong(3)
    let zero = MPInt.from_ilong(0)

    // add/sub/mul: always succeed
    h.assert_no_error({() ? => a.add_partial(b)? }, "add_partial ok")
    h.assert_no_error({() ? => a.sub_partial(b)? }, "sub_partial ok")
    h.assert_no_error({() ? => a.mul_partial(b)? }, "mul_partial ok")

    try
      h.assert_true(a.add_partial(b)? == MPInt.from_ilong(13), "add_partial result")
      h.assert_true(a.sub_partial(b)? == MPInt.from_ilong(7),  "sub_partial result")
      h.assert_true(a.mul_partial(b)? == MPInt.from_ilong(30), "mul_partial result")
    else
      h.fail("partial arith results")
    end

    // div/rem/fld/mod: succeed for non-zero
    h.assert_no_error({() ? => a.div_partial(b)?     }, "div_partial ok")
    h.assert_no_error({() ? => a.rem_partial(b)?     }, "rem_partial ok")
    h.assert_no_error({() ? => a.divrem_partial(b)?  }, "divrem_partial ok")
    h.assert_no_error({() ? => a.fld_partial(b)?     }, "fld_partial ok")
    h.assert_no_error({() ? => a.mod_partial(b)?     }, "mod_partial ok")

    try
      h.assert_true(a.div_partial(b)?    == MPInt.from_ilong(3), "div_partial result")
      h.assert_true(a.rem_partial(b)?    == MPInt.from_ilong(1), "rem_partial result")
      h.assert_true(a.fld_partial(b)?    == MPInt.from_ilong(3), "fld_partial result")
      h.assert_true(a.mod_partial(b)?    == MPInt.from_ilong(1), "mod_partial result")
      (let q, let r) = a.divrem_partial(b)?
      h.assert_true(q == MPInt.from_ilong(3), "divrem_partial q")
      h.assert_true(r == MPInt.from_ilong(1), "divrem_partial r")
    else
      h.fail("partial div results")
    end

    // error on zero divisor
    h.assert_error({() ? => a.div_partial(zero)?    }, "div_partial zero errors")
    h.assert_error({() ? => a.rem_partial(zero)?    }, "rem_partial zero errors")
    h.assert_error({() ? => a.divrem_partial(zero)? }, "divrem_partial zero errors")
    h.assert_error({() ? => a.fld_partial(zero)?    }, "fld_partial zero errors")
    h.assert_error({() ? => a.mod_partial(zero)?    }, "mod_partial zero errors")


class iso _TestMPIntUnsafeArith is UnitTest
  """
  Tests that _unsafe arithmetic variants produce the same results as their
  safe counterparts (for MPInt they are identical).
  """

  fun name(): String =>
    "MPInt/unsafe_arith"

  fun apply(h: TestHelper) =>
    let a = MPInt.from_ilong(42)
    let b = MPInt.from_ilong(7)

    h.assert_true(a.add_unsafe(b) == (a + b), "add_unsafe")
    h.assert_true(a.sub_unsafe(b) == (a - b), "sub_unsafe")
    h.assert_true(a.mul_unsafe(b) == (a * b), "mul_unsafe")
    h.assert_true(a.div_unsafe(b) == (a / b), "div_unsafe")
    h.assert_true(a.rem_unsafe(b) == (a % b), "rem_unsafe")
    h.assert_true(a.fld_unsafe(b) == a.fld(b), "fld_unsafe")
    h.assert_true(a.mod_unsafe(b) == a.mod(b), "mod_unsafe")

    (let q1, let r1) = a.divrem(b)
    (let q2, let r2) = a.divrem_unsafe(b)
    h.assert_true(q1 == q2, "divrem_unsafe q")
    h.assert_true(r1 == r2, "divrem_unsafe r")

    // Negative operands
    let ma = MPInt.from_ilong(-42)
    h.assert_true(ma.div_unsafe(b) == ma.div(b), "div_unsafe negative")
    h.assert_true(ma.rem_unsafe(b) == ma.rem(b), "rem_unsafe negative")
    h.assert_true(ma.fld_unsafe(b) == ma.fld(b), "fld_unsafe negative")
    h.assert_true(ma.mod_unsafe(b) == ma.mod(b), "mod_unsafe negative")


class iso _TestMPIntConversions is UnitTest
  """
  Tests for all numeric conversion methods: i8 through u128, f32, f64,
  and all _unsafe variants.
  """

  fun name(): String =>
    "MPInt/conversions"

  fun apply(h: TestHelper) =>
    let n = MPInt.from_ilong(42)
    let z = MPInt.from_ilong(0)

    // Signed integer round-trips for small value
    h.assert_true(n.i8()    == 42,   "i8(42)")
    h.assert_true(n.i16()   == 42,   "i16(42)")
    h.assert_true(n.i32()   == 42,   "i32(42)")
    h.assert_true(n.i64()   == 42,   "i64(42)")
    h.assert_true(n.i128()  == 42,   "i128(42)")
    h.assert_true(n.ilong() == 42,   "ilong(42)")
    h.assert_true(n.isize() == 42,   "isize(42)")

    // Unsigned integer round-trips for small value
    h.assert_true(n.u8()    == 42,   "u8(42)")
    h.assert_true(n.u16()   == 42,   "u16(42)")
    h.assert_true(n.u32()   == 42,   "u32(42)")
    h.assert_true(n.u64()   == 42,   "u64(42)")
    h.assert_true(n.u128()  == 42,   "u128(42)")
    h.assert_true(n.ulong() == 42,   "ulong(42)")
    h.assert_true(n.usize() == 42,   "usize(42)")

    // Zero
    h.assert_true(z.i8()   == 0,  "i8(0)")
    h.assert_true(z.u64()  == 0,  "u64(0)")
    h.assert_true(z.f32()  == 0,  "f32(0)")

    // Boundary: signed type limits
    h.assert_true(MPInt.from_ilong(127).i8()    ==  127, "i8 max")
    h.assert_true(MPInt.from_ilong(-128).i8()   == -128, "i8 min")
    h.assert_true(MPInt.from_ilong(32767).i16() ==  32767, "i16 max")
    h.assert_true(MPInt.from_ilong(-32768).i16() == -32768, "i16 min")

    // Boundary: unsigned type limits
    h.assert_true(MPInt.from_ilong(255).u8()   == 255,   "u8 max")
    h.assert_true(MPInt.from_ilong(65535).u16() == 65535, "u16 max")
    h.assert_true(MPInt.from_ilong(4294967295).u32() == 4294967295, "u32 max")

    // u64 round-trip for value above I64.max_value() (2^63)
    let big_pos = MPInt.from_ilong(ILong.max_value()) + MPInt.from_ilong(2)
    h.assert_true(big_pos.u64() == (ILong.max_value().u64() + 2), "u64 > 2^63")

    // u128 round-trip for value above I128.max_value()
    let i128_max = MPInt.from_ilong(ILong.max_value()).bit_shl(MPInt.from_ilong(64))
    let big128 = i128_max + MPInt.from_ilong(1)
    h.assert_true(big128.u128() != 0, "u128 > 2^127 is non-zero")

    // Float conversions
    h.assert_true(n.f32() == 42.0, "f32(42)")
    h.assert_true(n.f64() == 42.0, "f64(42)")
    // f32 exact up to 2^24
    h.assert_true(MPInt.from_ilong(16777216).f32() == 16777216.0, "f32 at 2^24")

    // Negative conversions
    h.assert_true(MPInt.from_ilong(-42).i8()  == -42, "i8(-42)")
    h.assert_true(MPInt.from_ilong(-42).i32() == -42, "i32(-42)")
    h.assert_true(MPInt.from_ilong(-42).f64() == -42.0, "f64(-42)")

    // _unsafe variants equal safe variants
    h.assert_true(n.i8_unsafe()    == n.i8(),    "i8_unsafe")
    h.assert_true(n.i16_unsafe()   == n.i16(),   "i16_unsafe")
    h.assert_true(n.i32_unsafe()   == n.i32(),   "i32_unsafe")
    h.assert_true(n.i64_unsafe()   == n.i64(),   "i64_unsafe")
    h.assert_true(n.i128_unsafe()  == n.i128(),  "i128_unsafe")
    h.assert_true(n.ilong_unsafe() == n.ilong(), "ilong_unsafe")
    h.assert_true(n.isize_unsafe() == n.isize(), "isize_unsafe")
    h.assert_true(n.u8_unsafe()    == n.u8(),    "u8_unsafe")
    h.assert_true(n.u16_unsafe()   == n.u16(),   "u16_unsafe")
    h.assert_true(n.u32_unsafe()   == n.u32(),   "u32_unsafe")
    h.assert_true(n.u64_unsafe()   == n.u64(),   "u64_unsafe")
    h.assert_true(n.u128_unsafe()  == n.u128(),  "u128_unsafe")
    h.assert_true(n.ulong_unsafe() == n.ulong(), "ulong_unsafe")
    h.assert_true(n.usize_unsafe() == n.usize(), "usize_unsafe")
    h.assert_true(n.f32_unsafe()   == n.f32(),   "f32_unsafe")
    h.assert_true(n.f64_unsafe()   == n.f64(),   "f64_unsafe")


class iso _TestMPIntShiftAlias is UnitTest
  """
  Tests that shl/shr/shl_unsafe/shr_unsafe delegate correctly to bit_shl/bit_shr.
  """

  fun name(): String =>
    "MPInt/shift_alias"

  fun apply(h: TestHelper) =>
    let a = MPInt.from_ilong(42)
    let n = MPInt.from_ilong(3)

    h.assert_true(a.shl(n)        == a.bit_shl(n), "shl == bit_shl")
    h.assert_true(a.shr(n)        == a.bit_shr(n), "shr == bit_shr")
    h.assert_true(a.shl_unsafe(n) == a.bit_shl(n), "shl_unsafe == bit_shl")
    h.assert_true(a.shr_unsafe(n) == a.bit_shr(n), "shr_unsafe == bit_shr")

    // Also verify with zero shift and large shift
    let zero_shift = MPInt.from_ilong(0)
    h.assert_true(a.shl(zero_shift) == a, "shl(0) = identity")
    h.assert_true(a.shr(zero_shift) == a, "shr(0) = identity")

    // Cross-digit shift
    let n16 = MPInt.from_ilong(16)
    h.assert_true(a.shl(n16) == a.bit_shl(n16), "shl cross-digit")
    h.assert_true(a.shl(n16).shr(n16) == a,     "shl/shr round-trip")


class iso _TestMPIntPowEdge is UnitTest
  """
  Edge cases for pow, gcd, isqrt not covered by _TestMPIntPow.
  """

  fun name(): String =>
    "MPInt/pow_edge"

  fun apply(h: TestHelper) =>
    let zero = MPInt.from_ilong(0)
    let one  = MPInt.from_ilong(1)

    // pow edge cases
    // 0^0 = 1 (convention: n.is_zero() returns 1 before checking base)
    h.assert_true(zero.pow(zero) == one, "0^0 = 1")
    // 1^n = 1
    h.assert_true(one.pow(MPInt.from_ilong(1000)) == one, "1^1000 = 1")
    // n^1 = n
    let five = MPInt.from_ilong(5)
    h.assert_true(five.pow(one) == five, "5^1 = 5")
    h.assert_true(MPInt.from_ilong(-3).pow(one) == MPInt.from_ilong(-3), "(-3)^1 = -3")
    // large exponent
    h.assert_true(
      MPInt.from_ilong(2).pow(MPInt.from_ilong(32)) == MPInt.from_ilong(4294967296),
      "2^32 = 4294967296")

    // gcd with negative inputs (gcd operates on abs values)
    h.assert_true(
      MPInt.from_ilong(-12).gcd(MPInt.from_ilong(8)) == MPInt.from_ilong(4),
      "gcd(-12, 8) = 4")
    h.assert_true(
      MPInt.from_ilong(12).gcd(MPInt.from_ilong(-8)) == MPInt.from_ilong(4),
      "gcd(12, -8) = 4")
    h.assert_true(
      MPInt.from_ilong(-12).gcd(MPInt.from_ilong(-8)) == MPInt.from_ilong(4),
      "gcd(-12, -8) = 4")
    // gcd(0, 0) = 0
    h.assert_true(zero.gcd(zero).is_zero(), "gcd(0, 0) = 0")
    // gcd(x, x) = x
    h.assert_true(
      MPInt.from_ilong(12).gcd(MPInt.from_ilong(12)) == MPInt.from_ilong(12),
      "gcd(12, 12) = 12")

    // isqrt for larger perfect squares
    // 65536 = 256^2
    h.assert_true(
      MPInt.from_ilong(65536).isqrt() == MPInt.from_ilong(256),
      "isqrt(65536) = 256")
    // 2^32 = (2^16)^2
    h.assert_true(
      MPInt.from_ilong(2).pow(MPInt.from_ilong(32)).isqrt() == MPInt.from_ilong(65536),
      "isqrt(2^32) = 65536")
    // Non-perfect-square: isqrt(2) = 1, isqrt(3) = 1, isqrt(8) = 2
    h.assert_true(MPInt.from_ilong(2).isqrt() == one,                "isqrt(2) = 1")
    h.assert_true(MPInt.from_ilong(3).isqrt() == one,                "isqrt(3) = 1")
    h.assert_true(MPInt.from_ilong(8).isqrt() == MPInt.from_ilong(2), "isqrt(8) = 2")
    // Just below a perfect square: isqrt(35) = 5 (since 5^2=25 ≤ 35 < 36=6^2)
    h.assert_true(MPInt.from_ilong(35).isqrt() == MPInt.from_ilong(5), "isqrt(35) = 5")
