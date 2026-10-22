#import "../template.typ": *

= Floating point arithmetic <cha:rounding_errors>
#minitoc()

#section-star[Introduction] <sec:introduction>

When we study numerical algorithms in the next chapters,
we assume implicitly that the operations involved are performed exactly.
On a computer, however, only a subset of the real numbers can be stored and,
consequently, many arithmetic operations are performed only approximately.
This is the source of the so-called _round-off errors_.
The rest of this chapter is organized as follows.
- In @sec:binary_representation_of_real_numbers,
  we discuss the binary representation of real numbers.

- In @sec:set_of_values,
  we describe the set of floating point numbers that can be represented in the usual floating point formats;

- In @sec:arithmetic_operations_between_floating_point_formats
  we explain how arithmetic operations between floating point numbers behave.
  We insist in particular on the fact that,
  in a calculation involving several successive arithmetic operations,
  the result of each intermediate operation is stored as a floating point number,
  with a possible error.

- In @sec:encoding_of_floating_point_numbers,
  we briefly present how floating point numbers are encoded
  according to the IEEE 754 standard, which is widely adopted today.
  We discuss also the encoding of special values such as #tt[Inf], #tt[-Inf] and #tt[NaN].

- Finally, in @sec:what_about_integer_formats,
  we present the standard integer formats and their encoding.
In order to completely describe computer arithmetic,
one would in principle need to also discuss the conversion mechanisms between different number formats,
as well as a number of edge cases.
// For instance, the reference IEEE 754 standard specifies that NaN == NaN must evaluate to false,
// which is somewhat arbitrary but may be useful in some situations.
A comprehensive discussion of the subject is beyond the scope of this course;
our aim in this chapter is only to introduce the key concepts.

== Binary representation of real numbers <sec:binary_representation_of_real_numbers>

Given any integer number $beta >= 2$, called the _base_,
a real number $x$ can always be expressed as a finite or infinite series of the form
#numbered($ x = plus.minus sum_(k=-n)^infinity a_k beta^(-k), wide a_k in {0, dots, beta - 1}. $) <eq:base_representation>
The number $x$ may then be denoted as $plus.minus (a_(-n) a_(-n+1) dots a_(-1) a_0 . a_1 a_2 dots)_beta$,
where the subscript $beta$ indicates the base.
This numeral system is called the _positional notation_ and is universally used today,
both by humans (usually with $beta=10$) and machines (usually with $beta=2$).
If the base~$beta$ is omitted,
it is always assumed in this course that $beta = 10$ unless otherwise specified
-- this is the _decimal_ representation.
The _digits_ $a_(-n), a_(-n+1), dots$ are also called _bits_ if $beta = 2$.
In computer science, several bases other than 10 are regularly employed,
for example the following:
- Base 2 (binary) is the usual choice for storing numbers on a machine.
  The binary format is convenient because the digits have only two possible values, 0 or 1,
  and so they can be stored using simple electrical circuits with two states.
  We employ the binary notation extensively in the rest of this chapter.
  Notice that, just like multiplying and dividing by~10 is easy in base 10,
  multiplying and dividing by 2 is very simple in base 2:
  these operations amount to shifting all the bits by one position to the left or right,
  respectively.

- Base 16 (hexadecimal) is sometimes convenient to represent numbers in a compact manner.
  In order to represent the values 0-15 with a single digit,
  16 different symbols are required, which are conventionally denoted by ${0,1,2,3,4,5,6,7,8,9,A,B,C,D,E,F}$.
  With this notation, we have $(F F)_16 = (255)_10$, for example.

  The hexadecimal notation is often used in programming languages for describing colors specified by a triplet $(r,g,b)$ of values between 0 and 255,
  corresponding to the primary colors _red_, _green_ and _blue_.
  The number of possible values for each component is~$256 = 16^2$,
  and so only 2 digits are required to represent these in the hexadecimal notation.
  Hexadecimal numbers are also employed in IPv6 addresses,
  which are used to identify computers connected to a network.

=== Conversion between binary and decimal formats
Obtaining the decimal representation of a binary number can be achieved from~@eq:base_representation,
using the decimal representations of the powers of 2.
Since all the positive and negative powers of 2 have a finite decimal representation,
any real number with a finite representation in base 2 has a finite representation also in base 10.
For example, $(0.01)_2 = (0.25)_10$ and $(0.111)_2 = (0.875)_10$.
#example(title: [Converting a binary number to decimal notation])[
  Let us calculate the decimal representation of $x = (0. overline(10))_2$,
  where the horizontal line indicates repetition: $x = (0.101010 dots)_2$.
  By definition, it holds that
  $ x = sum_(k=1)^infinity a_k 2^(-k), $
  where $a_k = 1$ if $k$ is odd and 0 otherwise.
  Thus, the series may be rewritten as
  $ x = sum_(k=0)^infinity 2^(-(2k+1)) = 1/2 sum_(k=0)^infinity (2^(-2))^k. $
  We recognize on the right-hand side a geometric series with common ratio $r = 2^(-2) = 1/4$,
  and so we obtain
  $ x = 1/2 lr((1/(1-r))) = 2/3 = (0. overline(6))_10. $
] <example:converting_binary_to_decimal>

#indent Obtaining the binary representation of a decimal number is more difficult,
because negative powers of 10 have _infinite_ binary representations,
as @exercise:binary_zero_point_one demonstrates.
There is, however, a simple procedure to perform the conversion,
which we present for the specific case of a real number $x$ with decimal representation of the form $x = (0.a_1 dots a_n)_10$.
In this setting,
the bits~$(b_1, b_2, dots)$ in the binary representation of $x = (0.b_1 b_2 b_3 dots)_2$ may be obtained as follows:
#algorithm(
  caption: [Conversion of a number to binary format],
  (0, [$i <- 1$]),
  (0, [*while* $x != 0$ *do*]),
  (1, [$x <- 2x$], "line:after_while"),
  (1, [*if* $x >= 1$ *then*]),
  (2, [$b_i <- 1$]),
  (1, [*else*]),
  (2, [$b_i <- 0$]),
  (1, [*end if*]),
  (1, [$x <- x - b_i$]),
  (1, [$i <- i+1$]),
  (0, [*end while*]),
) <algo:conversion_to_binary>

#example(title: [Converting a decimal number to binary notation])[
  Let us calculate the binary representation of $x = 1/3 = (0. overline(3))_10$.
  We apply~@algo:conversion_to_binary and collate the values of $i$ and $x$ obtained at the beginning of each iteration,
  i.e.~just before~#lineref(<line:after_while>), in the table below.
  #align(center, tabular(
    columns: 3,
    [$i$], [$x$], [Result],
    [$1$], [$1/3$], [0.#strong[0]000#sym.dots],
    [$2$], [$2/3$], [0.0#strong[1]00#sym.dots],
    [$3$], [$1/3$], [0.00#strong[0]0#sym.dots],
  ))
  Since $x$ in the last row is again $1/3$,
  successive bits alternate between 0 and 1,
  and the binary representation of $x$ is given by $(0. overline(01))_2$.
  This is not surprising since $2x = (0. overline(6))_10 = (0. overline(10))_2$,
  as we saw in~@example:converting_binary_to_decimal.
]

== Set of values representable in floating point formats <sec:set_of_values>
We mentioned in the introduction that,
because of memory limitations,
only a subset of the real numbers can be stored exactly in a computer.
Nowadays, the vast majority of programming languages comply with the IEEE 754 standard,
which requires that the set of representable numbers be of the form
#eq-last($
  floating(p, E_min, E_max)
  = Big(\{) & (-1)^s 2^E (b_0 . b_1 b_2 dots b_(p-1))_2 colon \
           & wide s in {0, 1}, b_i in {0, 1} thin "and" thin E_min <= E <= E_max Big(\}).
$) <eq:number_format>
In addition to these, floating number formats provide the special entities #tt[Inf], #tt[-Inf] and #tt[NaN],
the latter being an abbreviation for _Not a Number_.
Three parameters appear in the set definition~@eq:number_format.
The parameter $p in nat_(>0)$ is the number of significant bits (also called the precision),
and $(E_min, E_max) in integer^2$ are respectively the minimum and maximum exponents.
From the precision, the _machine epsilon_ is defined as $epsilon_M = 2^(-(p-1))$;
its significance is discussed in~@sub:relative_error_and_epsilon_machine.

For a number $x in floating(p, E_min, E_max)$,
$s$ is called the _sign_, $E$ is the _exponent_ and
$b_0 . b_1 b_2 dots b_(p-1)$ is the _significand_.
The latter can be divided into a _leading bit_ $b_0$ and the _fraction_ $b_1 b_2 dots b_(p-1)$,
to the right of the binary point.
The most widely used floating point formats are the _single_ and _double precision_ formats,
which are called respectively `Float32` and `Float64` in Julia.
Their parameters,
together with those of the lesser-known half-precision format,
are summarized in~@table:floating_point_formats.
In the rest of this section we use the shorthand notation $floating_16$, $floating_32$ and $floating_64$.
Note that $floating_16 subset floating_32 subset floating_64$.
#figure(
  placement: none,
  kind: table,
  tabular(
    columns: 4,
    [], [Half precision], [Single precision], [Double precision],
    [$p$], [11], [24], [53],
    [$E_min$], [-14], [-126], [-1022],
    [$E_max$], [15], [127], [1023],
  ),
  caption: [
    Floating point formats.
    The first column corresponds to the _half-precision_ format.
    This format,
    which is available through the `Float16` type in Julia,
    is more recent than the single and double precision formats.
    It was introduced in the 2008 revision to the IEEE 754 standard of 1985,
    a revision known as IEEE 754-2008.
  ],
) <table:floating_point_formats>

#remark[
  Some definitions,
  notably that in~@MR2265914[Section 2.5.2],
  include a general base~$beta$ instead of the base 2
  as an additional parameter in the definition of the number format~@eq:number_format.
  Since the binary format ($beta = 2$) is always employed in practice,
  we focus on this case for simplicity in most of this chapter.
]

#remark[
  Given a real number $x in floating(p, E_min, E_max)$,
  the exponent $E$ and significand are generally not uniquely defined.
  For example, the number $2.0 in floating_64$ may be expressed as $(-1)^0 2^1 (1.00 dots 00)_2$ or, equivalently, as $(-1)^0 2^2 (0.100 dots 00)_2$.
] <remark:roundoff_non_uniqueness>

#indent In #tt[Julia], non-integer number literals are interpreted as `Float64` by default,
which can be verified by using the `typeof` function.
For example, the instruction "`a = 0.1`" is equivalent to~"`a = Float64(0.1)`".
In order to define a number of type `Float32`,
the suffix `f0` must be appended to the decimal expansion.
For instance, the instruction "`a = 4.0f0`" defines a floating point number~`a` of type `Float32`;
it is equivalent to writing "`a = Float32(4.0)`".

=== Denormalized floating point numbers <sub:defining_a_unique_representation>
We can decompose the set $floating(p, E_min, E_max)$ in two disjoint parts:
$
  floating(p, E_min, E_max)
  &= Big(\{) (-1)^s 2^E (#text(fill: red, $1$) . b_1 b_2 dots b_(p-1))_2 colon \
  & wide s in {0, 1}, b_i in {0, 1} thin "and" thin E_min <= E <= E_max Big(\}) \
  & wide union Big(\{) (-1)^s 2^#text(fill: red, $E_min$) (#text(fill: red, $0$) . b_1 b_2 dots b_(p-1))_2 colon
   s in {0, 1}, b_i in {0, 1} Big(\}).
$
The numbers in the second set are called _subnormal_ or _denormalized_.

=== Relative error and machine epsilon <sub:relative_error_and_epsilon_machine>
Let $x$ be a nonzero real number and $hat(x)$ be an approximation.
We define the absolute and relative errors of the approximation as follows.
#definition(title: [Absolute and relative error])[
  The absolute error is given by $|x - hat(x)|$,
  whereas the relative error is
  $ (|x - hat(x)|)/(|x|) $
]

#indent The following result establishes a link between the machine $epsilon_M$ and
the relative error between a real number and the closest member of a floating point format.
#proposition[
  Let $x_min$ and $x_max$ denote the smallest and largest _non-denormalized_ positive numbers in a format $F = floating(p, E_min, E_max)$.
  If $x in [-x_max, -x_min] union [x_min, x_max]$,
  then
  #numbered($
    min_(hat(x) in F) (|x - hat(x)|)/(|x|) <= 1/2 2^(-(p-1))
    = 1/2 epsilon_M.
  $) <eq:bound_epsilon_machine>
] <prop:relative_error_and_machine_epsilon>
#proof[
  For simplicity, we assume that $x > 0$.
  Let $n = floor.l log_2(x) floor.r$ and $y := 2^(-n) x$.
  Since $y in [1, 2)$,
  it admits a binary representation of the form
  $(1.b_1 b_2 dots)_2 != (1. overline(1))_2$.
  Thus
  $x = 2^n (1.b_1 b_2 dots)_2,$
  and from the assumption that $x_min <= x <= x_max$
  we deduce that $E_min <= n <= E_max$.
  We now define the number $x_- in F$ by truncating the binary expansion of $x$ as follows:
  $ x_- = 2^n (1.b_1 dots b_(p-1))_2. $
  The distance between $x_-$ and its successor in $F$,
  which we denote by $x_+$, is given by $2^(n-p+1)$.
  Consequently, it holds that
  $ (x_+ - x) + (x-x_-) = x_+ - x_- = 2^(n-p+1). $
  Since both summands on the left-hand side are positive,
  this implies that either $x_+ - x$ or $x - x_-$ is bounded from above by $1/2 2^(n-p+1) <= 1/2 2^(-p+1)x$,
  which concludes the proof.
]

#indent The machine epsilon,
which was defined as $epsilon_M = 2^(-(p-1))$,
coincides with the maximum relative spacing between a non-denormalized floating point number $x$ and its successor in the floating point format,
defined as the smallest number in the format that is strictly larger than~$x$.

@fig:float64_density depicts the density of double-precision floating point numbers,
i.e.~the number of~$floating_64$ members per unit on the real line.
The figure shows that the density decreases as the absolute value of $x$ increases.
We also notice that the density is piecewise constant with discontinuities at powers of 2.
@fig:float64_spacing illustrates the relative spacing between successive floating point numbers.
Although the absolute spacing increases with the absolute value of $x$,
the relative spacing oscillates between~$1/2 epsilon_M$ and~$epsilon_M$.

#figure(
  placement: top,
  image("/figures/float64_density.pdf", width: 80%),
  caption: [
    Absolute spacing between double-precision floating point numbers,
    for $x in floating_64$.
    In this figure,
    $Delta(x)$ denotes the distance between $x$ and its successor in $floating_64$.
  ],
) <fig:float64_density>

#figure(
  placement: top,
  image("/figures/float64_spacing.pdf", width: 65%),
  caption: [
    Relative spacing between successive double-precision floating point numbers in the "normal range".
    The relative spacing oscillates between $1/2 epsilon_M$ and $epsilon_M$.
  ],
) <fig:float64_spacing>

The picture of the relative spacing between successive floating point numbers looks quite different for denormalized numbers.
This is illustrated in @fig:float64_spacing_denormalized,
which shows that the relative spacing increases beyond the machine epsilon in the denormalized range.
Fortunately, in the usual $floating_32$ and $floating_64$ formats,
the transition between denormalized and non-denormalized numbers occurs at such a small value that
it rarely needs worrying about.

#figure(
  placement: top,
  image("/figures/float64_spacing_denormalized.pdf", width: 65%),
  caption: [
    Relative spacing between successive double-precision floating point numbers,
    over a range which includes denormalized numbers.
    The vertical red line indicates the transition from denormalized to non-denormalized numbers.
  ],
) <fig:float64_spacing_denormalized>

#example[
  In Julia,
  the machine epsilon can be obtained using the `eps` function.
  For example, the instruction `eps(Float16)` returns $epsilon_M$ for the half-precision format.
]


== Arithmetic operations between floating point formats <sec:arithmetic_operations_between_floating_point_formats>

Now that we have presented the set of values representable on a computer,
we attempt in this section to understand precisely how arithmetic operations between floating point formats are performed.
The key mechanism governing arithmetic operations on a computer is that of _rounding_,
the action of approximating a real number regarded as infinitely precise by a number in a floating point format $floating(p, E_min, E_max)$.
The IEEE 754 standard stipulates that the default mechanism for rounding a real number $x$,
called _round to nearest_,
should behave as follows:
- *Standard case*:
  The number $x$ is rounded to the _nearest representable number_,
  if this number is unique.
- *Edge case*:
  When there are two equally near representable numbers in the floating point format,
  the one with the least significant bit equal to zero is delivered.
- *Infinities*:
  If the real number $x$ is larger than the largest representable number in the format,
  that is larger than $x_max = 2^(E_max) (2 - epsilon_M)$,
  then there are two cases,
  - If $x < 2^(E_max) (2 - 2^(-p))$, then $x_max$ is delivered;
  - Otherwise, the special value~#tt[Inf] is delivered.
  In other words, $x_max$ is delivered if it would be delivered by following the rules of the first two bullet points
  in a different floating point format with the same precision but a larger exponent $E_max$.
  A similar rule applies for large negative numbers.

When a binary arithmetic operation ($+$, $-$, $times$, $\/$) is performed on floating point numbers in format $floating$,
the result delivered by the computer is obtained by rounding the exact result of the operation according to the rules given above.
In other words,
the arithmetic operation is performed as if the computer first calculated an intermediate exact result,
and then rounded this intermediate result in order to provide a final result in $floating$.

Mathematically,
arithmetic operations between floating point numbers in a given format $floating$ may be formalized by introducing
the rounding operator $fl: real -> floating$ and
by defining, for any binary operation $compose in {+, -, times, \/}$,
the corresponding machine operation
$ hat(compose) colon floating times floating -> floating; (x, y) |-> fl(x compose y). $
We defined this operator for arguments in the same floating point format~$floating$.
If the arguments of a binary arithmetic operation are of different types,
the format of the end result,
known as the _destination format_,
depends on that of the arguments:
as a rule of thumb,
it is given by the most precise among the formats of the arguments.
In addition, recall that a floating point literal whose format is not explicitly specified is rounded to double-precision format
and so, for example, the addition `0.1 + 0.1` produces the result $fl_64 big((fl_64(0.1) + fl_64(0.1)))$,
where $fl_64$ is the rounding operator to the double-precision format.

#example[
  Using the `typeof` function,
  we check that the floating point literal `1.0` is indeed interpreted as a double-precision number:
  ```julia
  julia> a = 1.0; typeof(a)
  Float64
  ```
  When two numbers in different floating point formats are passed to a binary operation,
  the result is in the more precise format.
  ```julia
  julia> typeof(Float16(1) + Float32(1))
  Float32

  julia> typeof(Float32(1) + Float64(1))
  Float64
  ```
]

#indent If a mathematical expression contains several binary arithmetic operations to be performed in succession,
the result of each intermediate calculation is stored in a floating point format dictated by the formats of its argument,
and this floating point number is employed in the next binary operation.
A consequence of this mechanism is that the machine operations $madd$ and $mtimes$ are generally _not associative_.
For example, in general
$ (x madd y) madd z != x madd (y madd z) $
#example[
  Let $x = 1$ and $y = 3 times 2^(-13)$.
  Both of these numbers belong to $floating_16$ and,
  denoting by $hat(+)$ machine addition in~$floating_16$,
  we have
  #numbered($ (x madd y) madd y = 1 $) <eq:non_associativity_1>
  but
  #numbered($ x madd (y madd y) = 1 + 2^(-10). $) <eq:non_associativity_2>
  To explain this somewhat surprising result,
  we begin by writing the normalized representations of~$x$ and $y$ in the $floating_16$ format:
  $
    x &= (-1)^0 times 2^0 times (#tt[1.0000000000])_2 \
    y &= (-1)^0 times 2^(-12) times (#tt[1.1000000000])_2.
  $
  The exact result of the addition $x + y$ is given by $r = 1 + 3 times 2^(-13)$,
  which in binary notation is
  $ r = (1.underbrace(00000000000, "11 zeros")11)_2. $
  Since the length of the significand in the half-precision ($floating_16$) format is only $p = 11$,
  this number is not part of $floating_16$.
  The result of the machine addition $madd$ is therefore obtained by rounding $r$ to the nearest member of~$floating_16$,
  which is 1.
  This reasoning can then be repeated in order to conclude that, indeed,
  $ (x madd y) madd y = x madd y = 1. $

  In order to explain the result of~@eq:non_associativity_2,
  note that the exact result of the addition $y + y$ is $r = 3 times 2^(-12)$,
  which belongs to the floating point format,
  so it also holds that $y madd y = 3 times 2^(-12)$.
  Therefore,
  $ x madd (y madd y) = 1 madd 3 times 2^(-12) = fl_16(1 + 3 times 2^(-12)). $
  The argument of the $floating_16$ rounding operator does not belong to $floating_16$,
  since its binary representation is given by
  $ (1.underbrace(0000000000, "10 zeros")11)_2. $
  This time the nearest member of $floating_16$ is given by $1 + 2^(-10)$.
]

#indent When a numerical computation unexpectedly returns `Inf` or `-Inf`,
we say that an _overflow error_ occurred.
Similarly, _underflow_ is said to occur when a number
is smaller in absolute value than the smallest positive representable number in a floating point format.


== Encoding of floating point numbers~#moreinfo <sec:encoding_of_floating_point_numbers>

Once a number format is specified through parameters $(p, E_min, E_max)$,
the choice of encoding,
i.e.~the machine representation of numbers in this format,
has no bearing on the magnitude and propagation of round-off errors.
Studying encoding is, therefore, not essential for our purposes in this course,
but we opted to cover the topic anyway in the hope that it will help the students build intuition on floating point numbers.
We focus mainly on the single precision format,
but the following discussion applies _mutatis mutandis_ to the double and half-precision formats.
The material in this section is for information purposes only.

We already mentioned in~@remark:roundoff_non_uniqueness that
a number in a floating point format may have several representations.
On a computer, however, a floating point number is always stored in the same manner
(except for the number 0, see @remak:non_uniqueness_0).
The values of the exponent and significand which are selected by the computer,
in the case where there are several possible choices,
are determined from the following rules:
- Either $E > E_min$ and $b_0 = 1$;
- Or $E = E_min$, in which case the leading bit may be 0.
The following result proves that these rules enable to define the exponent and significand of a number in a set of floating point numbers uniquely.
#proposition[
  Assume that
  #numbered($
    (-1)^s 2^E (b_0 . b_1 dots b_(p-1))_2 = (-1)^(tilde(s)) 2^(tilde(E)) (tilde(b)_0 . tilde(b)_1 dots tilde(b)_(p-1))_2,
  $) <eq:equality_floating_point>
  where the parameter sets $(s, E, b_0, dots, b_(p-1))$ and $(tilde(s), tilde(E), tilde(b)_0, dots, tilde(b)_(p-1))$ both satisfy the above rule.
  Then $E = tilde(E)$ and $b_i = tilde(b)_i$ for $i in {0, dots, p-1}$.
] <prop:uniqueness_standard_representation>
#proof[
  We show that $E = tilde(E)$,
  after which the equality of significands follows trivially.
  Let us assume for contradiction that $E > tilde(E)$
  and denote the left and right-hand sides of~@eq:equality_floating_point by~$x$ and~$tilde(x)$, respectively.
  Then~$E > E_min$, implying that $b_0 = 1$ and so $2^E <= |x| < 2^(E+1)$.
  On the other hand, it holds that $|tilde(x)| < 2^(tilde(E)+1)$ regardless of whether $tilde(E) = E_min$ or not.
  Since $E >= tilde(E) + 1$ by assumption,
  we deduce that $|tilde(x)| < 2^E <= |x|$,
  which contradicts the equality $x = tilde(x)$.
]

#indent Now that we have explained how
a unique set of parameters (sign, exponent, significand) can be assigned to any floating point number,
we describe how these parameters are stored on the computer in practice.
As their names suggest, the `Float16`, `Float32` and `Float64` formats use 16, 32 and 64 bits of memory, respectively.
A naive approach for encoding these number formats would be to store the full binary representations of the sign, exponent and significand.

For the `Float32` format,
this approach would require 1 bit for the sign,
8 bits to cover the 254 possible values of the exponent, and 24 bits for the significand,
i.e.~for storing $b_0, dots, b_(p-1)$.
This leads to a total number of 33 bits,
which is one more than is available,
and this is without the special values `NaN`, `Inf` and `-Inf`.
So how are numbers in the $floating_32$ format actually stored?
To answer this question,
we begin with two observations:
- If $E > E_min$,
  then necessarily $b_0 = 1$ in the unique representation of the significand.
  Consequently, the leading bit need not be explicitly specified in this case;
  it is said to be _implicit_.
  We will see,
  as a consequence,
  that $p-1$ instead of $p$ bits are in fact sufficient for the significand.
- In the $floating_32$ format,
  8 bits at minimum need to be reserved for the exponent,
  which enables the representation of $2^8 = 256$ different values,
  but there are only 254 possible values for the exponent.
  This suggests that $256 - 254 = 2$ combinations of the 8 bits can be exploited in order to
  represent the special values #tt[Inf], #tt[-Inf] and #tt[NaN].

Simplifying a little bit,
we may view single precision floating point number as an array of 32 bits as illustrated below:
#align(center, {
  let field(fill, body) = box(fill: fill, stroke: 0.4pt + black, inset: 3pt, body)
  grid(
    columns: 3, column-gutter: 0.1cm, row-gutter: 3pt, align: left,
    field(rgb(84%, 100%, 100%))[Sign],
    field(rgb(64%, 100%, 71%))[#h(1em)Encoded exponent#h(1em)],
    field(lightred)[#h(2em)Encoded significand#h(2em)],
    text(size: 10pt)[1 bit], text(size: 10pt)[8 bits], text(size: 10pt)[23 bits],
  )
})
According to the IEEE 754 standard,
the first bit is the sign $s$,
the next 8 bits $e_0 e_1 dots e_6 e_7$ encode the exponent,
and the last 23 bits $b_1 b_2 dots b_(p-2) b_(p-1)$ encode the significand.
Let us emphasize that when we say "encode the exponent" here,
we just mean that the bits contain information from which the exponent can be uniquely determined,
but we have not yet described how this is achieved.
Let us introduce the integer number $e = (e_0 e_1 dots e_6 e_7)_2$;
that is to say, $0 <= e <= 2^8 -1$ is the integer number whose binary representation
is given by $e_0 e_1 dots e_6 e_7$.
One may determine the exponent and significand of a floating point number from the following rules.
- *Denormalized numbers*:
  If $e = 0$, then the implicit leading bit $b_0$ is zero,
  the fraction is $b_1 b_2 dots b_(p-2) b_(p-1)$, and the exponent is $E = E_min$.
  In other words, using the notation of~@sec:set_of_values,
  we have $x = (-1)^s 2^(E_min) (0.b_1 b_2 dots b_(p-2) b_(p-1))_2$.
  In particular, if $b_1 b_2 dots b_(p-2) b_(p-1) = 00 dots 00$,
  then it holds that $x = 0$.
- *Non-denormalized numbers*:
  If $0 < e < 255$,
  then the implicit leading bit $b_0$ of the significand is $1$
  and the fraction is given by $b_1 b_2 dots b_(p-2) b_(p-1)$.
  The exponent is given by
  $ E = e - "bias" = E_min + e - 1. $
  where the exponent bias for the single and double precision formats are given in~@table:floating_point_formats_encoding.
  In this case $x = (-1)^s 2^(e - "bias") (1.b_1 b_2 dots b_(p-2) b_(p-1))_2$.
  Notice that $E = E_min$ if $e = 1$,
  as in the case of denormalized numbers.
- *Infinities*:
  If $e = 255$ and $b_1 b_2 dots b_(p-2) b_(p-1) = 00 dots 00$,
  then $x = #tt[Inf]$ if $s = 0$ and #tt[-Inf] otherwise.
- *Not a Number*:
  If $e = 255$ and $b_1 b_2 dots b_(p-2) b_(p-1) != 00 dots 00$,
  then $x = #tt[NaN]$.
  Notice that the special value #tt[NaN] can be encoded in many different manners.
  These extra degrees of freedom were reserved for passing information on the reason for the occurrence of $#tt[NaN]$,
  which is usually an indication that something has gone wrong in the calculation.

#figure(
  placement: none,
  kind: table,
  tabular(
    columns: 4,
    [], [Half precision], [Single precision], [Double precision],
    [Exponent bias ($-E_min + 1$)], [$15$], [$127$], [$1023$],
    [Exponent encoding (bits)], [5], [8], [11],
    [Significand encoding (bits)], [10], [23], [52],
  ),
  caption: [Encoding parameters for floating point formats],
) <table:floating_point_formats_encoding>

#remark(title: [Encoding efficiency])[
  With 32 bits,
  at most $2^32$ different numbers could in principle be represented.
  In practice, as we saw in @exercise:cardinality_floating,
  the `Float32` format enables to represent
  $
    (E_max - E_min) 2^p + 2^(p+1) - 1
    = 253 times 2^24 + 2^25 - 1
    = 2^32 - 2^24 - 1
    approx 99.6 % times 2^32,
  $
  different real numbers,
  indicating a very good encoding efficiency.
]

#remark(title: [Nonuniqueness of the floating point representation of 0.0])[
  The sign $s$ is clearly unique for any number in a floating point format,
  except for $0.0$,
  which could in principle be represented as
  $
    (-1)^0 2^(E_min) (0.00 dots 00)_2
    wide "or" wide
    (-1)^1 2^(E_min) (0.00 dots 00)_2.
  $
  In practice, both representations of 0.0 are available on most machines,
  and these behave slightly differently.
  For example $#tt[1/(0.0) = Inf]$ but $#tt[1/(-0.0) = -Inf]$.
] <remak:non_uniqueness_0>

== Integer formats~#moreinfo <sec:what_about_integer_formats>

The machine representation of integer formats is much simpler than that of floating point numbers.
In this short section,
we give a few orders of magnitude for common integer formats and briefly discuss overflow issues.
Programming languages typically provide integer formats based on 16, 32 and 64 bits.
In Julia, these correspond to the types `Int16`, `Int32` and `Int64`,
the latter being the default for integer literals.

The most common encoding for integer numbers,
which is used in Julia, is known as _two's complement_:
a number encoded with $p$ bits given by $b_(p-1) b_(p-2) dots b_0$ corresponds to
$ x = - b_(p-1) 2^(p-1) + sum_(i=0)^(p-2) b_i 2^i. $
This encoding enables to represent uniquely all the integers from $N_min = -2^(p-1)$ to $N_max = 2^(p-1)-1$.
In contrast with floating point formats,
integer formats do not provide special values like `Inf` and `NaN`.
The number delivered by the machine when a calculation exceeds the maximum representable value in the format,
called the _overflow behavior_, generally depends on the programming language.

Since the overflow behavior of integer formats is not universal across programming languages,
a detailed discussion is of little interest.
We only mention that Julia uses a _wraparound_ behavior,
where $N_max + 1$ silently returns $N_min$ and,
similarly, $N_min - 1$ gives $N_max$;
the numbers loop back.
This can lead to unexpected results,
such as `2^64` evaluating to 0.

== Exercises
#exercise[
  Show that if a number $x in real$ admits a finite representation~@eq:base_representation in base~$beta$,
  then it also admits an infinite representation in the same base.
  *Hint:* You may have learned before that $(0. overline(9))_10 = 1$.
] <exercise:non_unique_representation>

#exercise[
  How many digits does it take to represent all the integers from 0 to $10^10 - 1$ in decimal and binary formats?
  What about the hexadecimal format?
]

#exercise[
  Find the decimal representation of $(0.000 overline(1100))_2$.
]

#exercise[
  Find the binary representation of $(0.1)_10$.
] <exercise:binary_zero_point_one>

#compexercise[
  Implement~@algo:conversion_to_binary on a computer and verify that it works.
  Your function should take two arguments: an array of integers #verb("[a_1, ..., a_n]") containing the digits after the decimal point
  and the number~$m$ of bits to return.
  The bits should be _returned_ as an array~#verb("[b_1, ..., b_m]").
]

#compexercise[
  As mentioned above, @algo:conversion_to_binary works only for decimal numbers of the specific form~$x = (0.a_1 dots a_n)_10$.
  Find and implement a similar algorithm for integer numbers.
  More precisely, write a function that takes an integer $n$ as argument and
  returns an array containing the bits of the binary expansion $(b_m dots b_0)_2$ of $n$,
  from the least significant $b_0$ to the most significant $b_m$.
  That is to say, your code should return `[b_0, b_1, ...]`.
  ```julia
  function to_binary(n)
      # Your code comes here
  end

  # Check that it works
  number = 123456789
  bits = to_binary(number)
  pows2 = 2 .^ range(0, length(bits) - 1)
  @assert sum(bits'pows2) == number
  ```
]

#exercise[
  Show that the successor of $1$ in $floating_64$ is $1 + epsilon_64$,
  where $epsilon_64$ is the machine epsilon for the double-precision format.
]

#exercise[
  Write down the values of the smallest and largest, in absolute value,
  positive real numbers representable in the $floating_32$ and $floating_64$ formats.
]

#exercise(title: [Relative error and machine epsilon])[
  Prove that the inequality~@eq:bound_epsilon_machine is sharp.
  More precisely, find $x in real$ such that the inequality is an equality.
] <exercise:machine_epsilon>

#exercise(title: [Cardinality of the set of floating point numbers])[
  Show that, if $E_max >= E_min$, then $floating(p, E_min, E_max)$ contains exactly
  $ (E_max - E_min) 2^p + 2^(p+1) - 1 $
  _distinct real numbers_.
  (In particular, the special values #tt[Inf], #tt[-Inf] and #tt[NaN] are not counted.)
  *Hint:* Count first the numbers with $E > E_min$ and then those with $E = E_min$.
] <exercise:cardinality_floating>

#exercise[
  Calculate the machine epsilon $epsilon_16$ for the $floating_16$ format.
  Write the results of the arithmetic operations $1 madd epsilon_16$ and $1 msub epsilon_16$ in the form
  $ 2^E (1.b_1 dots b_(p-1))_2. $
]

#exercise[
  Let $epsilon_16$ be the machine epsilon for the $floating_16$ format,
  and define $y = 4/3 epsilon_16$.
  What is the relative error between $Delta = (1 + y) - 1$,
  and the machine approximation $hat(Delta) = (1 madd y) msub 1$?
]

#compexercise(title: [Numerical differentiation])[
  Let $f(x) = exp(x)$.
  By definition, the derivative of~$f$ at~$0$ is given by
  $ f'(0) = lim_(delta -> 0) lr(((f(delta) - f(0))/delta)). $
  The expression within brackets on the right-hand side may be used
  with a small but nonzero~$delta$ as an approximation for $f'(0)$.
  Implement this approach using double-precision numbers and the same values for $delta$ as in the table below.
  Explain the results you obtain.
  #align(center, {
    // \arraystretch = 1.5
    set text(style: "normal")
    tabular(
      columns: 4, stretch: 1.5,
      [$delta$], [$epsilon_64/4$], [$epsilon_64/2$], [$epsilon_64$],
      [Approximation of $f'(0)$], [0], [2], [1],
    )
  })
]

#compexercise(title: [Avoiding overflow])[
  Write a code to calculate the weighted average
  $
    S := (sum_(j=0)^J w_j j) / (sum_(j=0)^J w_j),
    wide w_j = exp(j),
    wide J = 1000.
  $
  You may need to first rewrite~$S$ differently.
]

#compexercise[
  Plot the function $x |-> log lr((ee^(ee^x) - 1))$ over the interval $[0, 10]$.
]

#compexercise(title: [Calculating the sample variance])[
  Assume that $(x_n)_(1 <= n <= N)$, with $N = 10^6$, are independent random variables distributed according to
  the uniform distribution $cal(U)(L, L+1)$.
  That is, each~$x_n$ takes a random value uniformly distributed between $L$ and $L+1$ where $L = 10^9$.
  In Julia, these samples can be generated with the following lines of code:
  ```julia
  N, L = 10^6, 10^9
  x = L .+ rand(N)
  ```
  It is well known that the variance of $x_n tilde.op cal(U)(L, L+1)$ is given by~$sigma^2 = 1/12$.
  Numerically, the variance can be estimated from the _sample variance_:
  #numbered($
    s^2 = 1/(N-1) lr((lr((sum_(n=1)^N x_n^2)) - N macron(x)^2)),
    wide macron(x) = 1/N sum_(n=1)^N x_n.
  $) <eq:sample_variance_2>
  Write a computer code to calculate $s^2$ with the best possible accuracy.
  Can you find a formula that enables better accuracy than~@eq:sample_variance_2?
]

#remark[
  In order to estimate the true value of $s^2$ for your samples,
  you can use the `BigFloat` format,
  to which the array $x$ can be converted by using the instruction
  ```julia
  x = BigFloat.(x)
  ```
]

#compexercise[
  Euler proved that
  $ pi^2/6 = lim_(N -> infinity) sum_(n=1)^N 1/n^2. $
  Using the default `Float64` format,
  estimate the error obtained when the series on the right-hand side is truncated after $10^10$ terms.
  Can you rearrange the sum for best accuracy?
]

#exercise[
  Let $x$ and $y$ be positive real numbers in the interval $[2^(-10), 2^10]$
  (so that we do not need to worry about denormalized numbers,
  assuming we are working in single or double precision),
  and let us define the machine addition operator $madd$ for arguments in real numbers as
  $ madd: real times real -> real; (x, y) |-> fl big((fl(x) + fl(y))). $
  Prove the following bound on the relative error between the sum $x+y$ and its machine approximation $x madd y$:
  $
    (big(|(x + y) - (x madd y)|))/(|x + y|)
    <= epsilon_M/2 lr((2 + epsilon_M/2)).
  $
  *Hint:* decompose the numerator as
  $
    (x + y) - (x madd y) = big((x - fl(x))) + big((y - fl(y))) + big((fl(x) + fl(y) - (x madd y))),
  $
  and then use @prop:relative_error_and_machine_epsilon.
]

#exercise[
  Is `Float32(0.1) * Float32(10) == 1` equal to `true` or `false`
  given the default rounding rule defined by the IEEE standard? Explain.
] <exercise:floating_point_illustration>

#solution[
  By default, real numbers are rounded to the nearest floating point number.
  This can be checked in Julia with the command `rounding(Float32)`,
  which prints the default rounding mode.
  The exact binary representation of the real number $x = 0.1$ is
  $
    x
    &= (0.000 overline(1100))_2 \
    &= 2^(-4) times (underbrace(1.10011001100110011001100, "24 bits") overline(1100))_2
  $
  The first task is to determine the member of $floating_32$ that is nearest $x$.
  We have
  $
    x^- := max big(\{) y: y in floating_32 "and" y <= x big(\}) &= 2^(-4) times (1.10011001100110011001100)_2 \
    x^+ := min big(\{) y: y in floating_32 "and" y >= x big(\}) &= 2^(-4) times (1.10011001100110011001101)_2.
  $
  Since the number $(0. overline(1100))_2$ is closer to 1 than to 0,
  the number $x$ is closer to $x^+$ than to $x^-$.
  Therefore, the number obtained when writing `Float32(0.1)` is $x^+$.
  To conclude the exercise,
  we need to calculate~$fl big((10 times x^+))$,
  and to this end we first write the exact binary representation of the real number $10 times x^+ = (1010)_2 times x^+$.
  We have
  $
    (1010)_2 times x^+
    = (1000)_2 times x^+ + (10)_2 times x^+
    &= 2^(-4) times space (1100.11001100110011001101)_2 \
    & + 2^(-4) times #h(.45cm) space (11.0011001100110011001101)_2 \
    &= 2^(-4) times (underbrace(10000.0000000000000000000, "24 bits") #h(-0.17em) 001)_2.
  $
  This can be checked in Julia by writing `bitstring(Float32(0.1) * Float64(10.0))`.
  Clearly, when rounding to the nearest $floating_32$ number,
  the number $2^(-4) (10000)_2 = 1$ is obtained;
  the boolean expression in the question is thus `true`.
]

#remark[
  It should not be inferred from @exercise:floating_point_illustration that `Float32(1/i) * i` is always exact in floating point arithmetic.
  For example `Float32(1/41) * 41` does not evaluate to `1`,
  and neither do `Float16(1/11) * 11` and `Float64(1/49) * 49`.
]

#exercise[
  Given that the default rounding rule specified by the IEEE 754 standard is "round to nearest, tie to even",
  does `Float16(0.1) + Float16(0.2) == Float16(0.3)` evaluate to `true` or `false`?
  Explain.
]

#exercise[
  Explain why `Float32(sqrt(2))^2 - 2` is not zero in Julia.
]
#solution[
  The exact binary representation of $x := sqrt(2)$, found using the `bitstring` in the `Float64` format,
  is given by
  $
    x
    &= (underbrace(1.01101010000010011110011, "24 bits") 001100 dots)_2.
  $
  The first task is to determine the member of $floating_32$ that is nearest $x$.
  We have
  $
    x^- &:= max big(\{) x: x in floating_32 "and" x <= sqrt(2) big(\}) = (underbrace(1.01101010000010011110011, "24 bits"))_2 \
    x^+ &:= min big(\{) x: x in floating_32 "and" x >= sqrt(2) big(\}) = (underbrace(1.01101010000010011110100, "24 bits"))_2,
  $
  and we calculate
  $
    x - x^- &= 2^(-24) (0.01100 dots)_2, \
    x^+ - x &= 2^(-21) big((1 - (0.11001100 dots)_2))
    >= 2^(-21) big((1 - (0.11001101)_2)) =
    2^(-21) lr((0.00110011))_2.
  $
  We deduce that $x - x^- < x^+ - x$,
  and so $fl(x) = x^-$.
  To conclude the exercise,
  we need to show that $fl big(((x^-)^2))$ is not equal to 2.
  The exact binary expansion of $(x^-)^2$ is
  $ (x^-)^2 = (underbrace(1.11111111111111111111111, "24 bits")01101100111111010101001)_2. $
  The member of $floating_32$ nearest this number is
  $ (1.11111111111111111111111)_2 = 2 - 2^(-23), $
  which is precisely the result returned by Julia.
]

#compexercise(title: [Numerical differentiation])[
  Let $f(x) = exp(x)$ and let $d(delta)$ be the approximation of $f'(x)$ obtained from the following piece of code:
  ```julia
  f, x = exp, 1
  d(δ) = (f(x+δ) - f(x))/δ
  ```
  Plot for fixed $x = 1$ the error `abs(d(δ) - exp(x))` as a function of~`δ` in logarithmic scale,
  and explain the result.
] <exercise_differentiation>
#solution[
  We assume that $delta in floating_64$ for simplicity, and we let $epsilon = epsilon_M \/ 2$.
  This is not a restrictive assumption
  as~$delta$ can only take floating point values in computer programs.
  The proof is rather technical, and so it is given for information purposes only.
  We rewrite $d(delta)$ mathematically as
  $ d(delta) = Big((hat(f)(x madd delta) msub hat(f)(x))) mdiv delta, $
  where $hat(f)(x) = fl big((f(x)))$.
  We wish to bound $abs(f'(x) - d(delta))$.
  By the triangle inequality,
  #numbered($
    abs(f'(x) - d(delta))
    <= lr(|f'(x) - (f(x + delta) - f(x))/delta|)
    + lr(|(f(x + delta) - f(x))/delta - Big((hat(f)(x madd delta) msub hat(f)(x))) mdiv delta|).
  $) <eq:triangle_inequality_derivative>
  #paragraph[Bounding the first term in~@eq:triangle_inequality_derivative.][
    By Taylor's theorem, there exists $xi in [x, x+delta]$ such that
    $ f(x + delta) = f(x) + delta f'(x) + delta^2/2 f''(xi). $
    Therefore, the first term in~@eq:triangle_inequality_derivative satisfies
    $
      lr(|f'(x) - (f(x + delta) - f(x))/delta|)
      = delta/2 |f''(xi)|
      = delta/2 |f''(x)| + bigo(delta^2).
    $
  ]

  #paragraph[Bounding the second term in~@eq:triangle_inequality_derivative -- the roundoff error.][
    This is more tedious but not difficult;
    the main ingredient is the triangle inequality.
    Specifically, we will use the bound
    #eq-last($
      &lr(|(f(x + delta) - f(x))/delta - Big((hat(f)(x madd delta) msub hat(f)(x))) mdiv delta|) \
      &wide <= delta^(-1) lr(|f(x+ delta) - f(x madd delta)|)
       + delta^(-1) lr(|f(x madd delta) - hat(f)(x madd delta)|)
       + delta^(-1) lr(|hat(f)(x) - f(x)|) \
      &wide quad + delta^(-1) lr(|Big((hat(f)(x madd delta) - hat(f)(x))) - Big((hat(f)(x madd delta) msub hat(f)(x)))|) \
      &wide quad + lr(|Big((hat(f)(x madd delta) msub hat(f)(x))) \/ delta - Big((hat(f)(x madd delta) msub hat(f)(x))) mdiv delta|).
    $) <eq:numerical_differentiation_big_triangle>
    Note that without absolute values,
    both sides are indeed equal.
    The first three terms on the right-hand side are together a bound from above for
    #numbered($
      lr(|(f(x+delta) - f(x))/delta - (hat(f)(x madd delta) - hat(f)(x))/delta|),
    $) <eq:bound_above_first_three_terms>
    while the two other terms account for the roundoff errors associated with the machine subtraction and division operators,
    respectively.
    We will show that the dominant terms in~@eq:numerical_differentiation_big_triangle are the first three;
    the latter two are negligible in comparison.
    Employing~@prop:relative_error_and_machine_epsilon with $x = f(a)$,
    with $x = a plus.minus b$ and with $x = a\/b$,
    we deduce that the following inequalities are satisfied:
    #numbered($
      &lr(|hat(f)(a) - f(a)|) <= epsilon |f(a)|,
      wide |(a plus.minus b) - (a thin hat(plus.minus) thin b)| <= epsilon abs(a plus.minus b),
      wide abs(a\/b - a mdiv b) <= epsilon |a\/b|.
    $) <eq:inequalities_numerical_differentiation>
    The first two inequalities are valid for all~$(a, b) in floating_64 times floating_64$,
    while the third inequality is valid for all~$(a, b) in floating_64 times floating_64 backslash {0}$.
    The first inequality in~@eq:inequalities_numerical_differentiation
    can be employed in order to bound the second and third term on the right-hand side of~@eq:numerical_differentiation_big_triangle:
    #eq-last($
      delta^(-1) lr(|f(x madd delta) - hat(f)(x madd delta)|)
      + delta^(-1) lr(|hat(f)(x) - f(x)|)
      &<= delta^(-1) epsilon Big((big(|f(x + delta)|) + big(|f(x)|))) \
      &= 2 delta^(-1) epsilon abs(f(x)) + bigo(epsilon).
    $) <eq:num_diff_1>
    Using Taylor's theorem and then the second inequality in~@eq:inequalities_numerical_differentiation (with $+$),
    we then bound the first term on right-hand side of~@eq:numerical_differentiation_big_triangle as
    #eq-last($
      delta^(-1) lr(|f(x + delta) - f(x madd delta)|)
      &= delta^(-1) lr(|f'(xi) big(((x + delta) - (x madd delta)))|) \
      &<= delta^(-1) epsilon |f'(xi)| |x + delta|
      = delta^(-1) epsilon |f'(x) x| + bigo(epsilon).
    $) <eq:num_diff_2>
    Combining~@eq:num_diff_1 and~@eq:num_diff_2,
    we obtain that the expression in~@eq:bound_above_first_three_terms scales as $delta^(-1) epsilon$:
    #numbered($
      lr(|(f(x+delta) - f(x))/delta - (hat(f)(x madd delta) - hat(f)(x))/delta|)
       = bigo(delta^(-1) epsilon).
    $) <eq:numerical_differentiation_term1>
    Next, using the second inequality in~@eq:inequalities_numerical_differentiation (with $-$) and
    then~@eq:numerical_differentiation_term1,
    we bound the fourth term in~@eq:numerical_differentiation_big_triangle as follows:
    #eq-last($
      &delta^(-1) lr(|big((hat(f)(x madd delta) - hat(f)(x))) - big((hat(f)(x madd delta) msub hat(f)(x)))|) \
      &wide <= delta^(-1) epsilon |hat(f)(x madd delta) - hat(f)(x)|
      = delta^(-1) epsilon underbrace(|f(x + delta) - f(x)|, = bigo(delta)) + bigo(delta^(-1) epsilon^2)
      = bigo lr((epsilon + delta^(-1) epsilon^2)).
    $) <eq:numerical_differentiation_term2>
    Notice that this term is smaller than the first three in~@eq:numerical_differentiation_big_triangle when $delta << 1$.
    We deduce by a triangle inequality from~@eq:numerical_differentiation_term1 and~@eq:numerical_differentiation_term2 that
    $
      delta^(-1) lr(|big((f(x + delta) - f(x)))
      - big((hat(f)(x madd delta) msub hat(f)(x)))|)
      &= bigo(delta^(-1) epsilon).
    $
    Finally, using the third inequality in~@eq:inequalities_numerical_differentiation together with this equation,
    we bound the fifth term on the right-hand side of~@eq:numerical_differentiation_big_triangle:
    $
      &lr(|big((hat(f)(x madd delta) msub hat(f)(x))) \/ delta
      - big((hat(f)(x madd delta) msub hat(f)(x))) mdiv delta|) \
      &wide <= epsilon big(|hat(f)(x madd delta) msub hat(f)(x)|) \/ delta
      = epsilon/delta underbrace(big(|f(x + delta) - f(x)|), bigo(delta))
      + thin bigo(delta^(-1) epsilon^2)
      = bigo(epsilon + delta^(-1) epsilon^2).
    $
    This term is also negligible in front of the other dominant contributions to the roundoff error given in~@eq:num_diff_1 and~@eq:num_diff_2.
    Going back to~@eq:numerical_differentiation_big_triangle,
    we conclude that
    $
      lr(|(f(x + delta) - f(x))/delta - Big((hat(f)(x madd delta) msub hat(f)(x))) mdiv delta|)
      <= delta^(-1) epsilon Big((2 |f(x)| + |f'(x)x|)) + bigo(epsilon + delta^(-1) epsilon^2).
    $
  ]

  #paragraph[Concluding the proof.][
    Going back to~@eq:triangle_inequality_derivative,
    we conclude that
    $ |f'(x) - d(delta)| <= delta/2 |f''(x)| + epsilon/delta Big((2 |f(x)| + |f'(x) x|)), $
    up to higher order terms.
    For fixed~$x$,
    the right-hand side is minimized when
    $ delta = sqrt(epsilon) sqrt((2 |f(x)| + |f'(x) x|)/(abs(f''(x)))), $
    which is a well-known formula for the optimal step size in numerical differentiation.
    The error as a function of~$delta$ for $x = 1$ is depicted in~@fig:differentiation.
  ]
]
#figure(
  placement: none,
  image("/figures/numerical_differentiation.pdf", width: 80%),
  caption: [Solution to @exercise_differentiation.],
) <fig:differentiation>

#exercise[
  Explain why `exp(log(Float16(7))) == 7` is `false`.
]
#solution[
  We begin by finding the binary representation of $log(7)$.
  In the $floating_64$ format, the sign and exponent are encoded over 1 and 11 bits respectively,
  and so the fraction is given by the last 52 bits, which can be obtained from the command
  ```julia
  julia> bitstring(log(7))[13:end]
  "1111001000100111001010101110001100100101101001010111"
  ```
  Therefore, since `exponent(log(7))` returns 0, we have
  #numbered($
    log(7) = (underbrace(1.1111001000, "11 bits")10011100101010111000110010010110100101011...)_2
  $) <eq:log7_in_binary>
  The number returned by the command `log(Float16(7))`
  is given by $fl_16 big((log(7)))$,
  where $fl_16$ denotes the half-precision rounding operator.
  Rounding the right-hand side of~@eq:log7_in_binary to 11 bits, we obtain
  $ fl_16 big((log(7))) = (1.1111001001)_2 = 1.9462890625. $
  The number returned by the code `exp(log(Float16(7)))` is
  $
    fl_16 bigg((exp Big((fl_16 big((log(7))))))))
    = fl_16 big((exp(1.9462890625))).
  $
  The rounding operator appears twice on the left-hand side,
  because the computer rounds _after every operation_.
  To explain the result of the rounding operation,
  we begin by calculating the binary expansion of $exp(1.9462890625)$.
  ```julia
  julia> bitstring(exp(1.9462890625))[13:end]
  "1100000000101011011101110000111000100001100010001110"

  julia> exponent(exp(1.9462890625))
  2
  ```
  Therefore,
  $
    exp(1.9462890625)
    = 2^2 (underbrace(1.1100000000, "11 bits")1010110111011100001110001000011000100011...)_2,
  $
  and, rounding to 11 bits,
  we finally obtain
  $
    fl_16 big((exp(1.9462890625)))
    &= 2^2 (1.110000000#text(fill: red, $1$))_2 = 4 lr((1 + 1/2 + 1/4 + 1/2^10)) \
    &= 7 + 1/2^8 = 7.00390625,
  $
  which is different from 7.
]

#exercise[
  Determine the encoding of the following `Float32` numbers:
  - $x_1 = 2.0^(E_min)$
  - $x_2 = - 2.0^(E_min - p + 1) = - 2.0^(-149)$
  - $x_3 = 2.0^(E_max) (2-2^(-p+1))$
  Check your results using the Julia function `bitstring`.
]

#exercise(title: [Summary])[
  True or false?
  + Let $(placeholder)_2$ denote binary representation.
    It holds that
    $(0.1111)_2 + (0.0001)_2 = 1.$
  + It holds that
    $(1000)_2 times (0.001)_2 = 1.$
  + It holds that
    $ (0. overline(1))_3 = 1/2. $
  + In base 16, all the natural numbers from 1 to 200 can be represented using 2 digits.
  + In Julia, `Float64(.1) == Float32(.1)` evaluates to `true`.
  + The spacing (in absolute value) between successive double-precision (`Float64`) floating point numbers is constant.
  + It holds that $(0. overline(10101))_2 = (1.2345)_10$.
  + Machine addition~$madd$ is an associative operation.
    More precisely, given any three double-precision floating point numbers $x$, $y$ and $z$,
    the following equality holds:
    $ (x madd y) madd z = x madd (y madd z). $
  + The machine epsilon is the smallest strictly positive number that can be represented in a floating point format.
  + Let $epsilon$ denote the machine epsilon for the double-precision format.
    Let also~$madd$ and~$mdiv$ denote respectively the machine addition and the machine division operators for the double-precision format.
    It holds that $1 madd (epsilon mdiv 64) = 1$ and that $epsilon mdiv 64 != 0$.
  + Assume that $x in real$ belongs to the double-precision floating point format (that is,
    assume that $x in floating_64$).
    Then $-x in floating_64$.
]

== Discussion and bibliography <sec:discussion_and_bibliograhpy>
This chapter is mostly based on the original 1985 IEEE 754 standard~@ieee754 and the reference book~@MR2265914.
A~significant revision to the 1985 IEEE standard was published in 2008~@ieee2008,
adding for example specifications for the half precision and quad precision formats,
and a minor revision was published in 2019~@ieee2019.
The original IEEE standard and its revisions constitute the authoritative guide on floating point formats.
It was intended to be widely disseminated and is written very clearly and concisely,
but is not available for free online.
Another excellent source for learning about floating point numbers and round-off errors is D.~Goldberg's~paper "_What every computer scientist should know about floating-point arithmetic_"~@goldberg1991every,
freely available online.
