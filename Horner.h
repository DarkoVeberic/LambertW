#ifndef _Horner_h_
#define _Horner_h_


template<typename Float, class Tag, unsigned int order>
struct Polynomial {
  static inline Float Coeff();
  static inline Float Coeff(const Float x);
};


template<typename Float, class Tag, unsigned int order>
struct Horner {
  static Float Recurse(const Float term, const Float x)
  { return Horner<Float, Tag, order-1>::Recurse(term*x + Polynomial<Float, Tag, order>::Coeff(), x); }
  static Float RecurseAlt(const Float term, const Float x)
  { return Horner<Float, Tag, order-1>::RecurseAlt(Polynomial<Float, Tag, order>::Coeff() + x*term, x); }
  static Float Eval(const Float x)
  { return Horner<Float, Tag, order-1>::Recurse(Polynomial<Float, Tag, order>::Coeff(), x); }
  static Float EvalAlt(const Float x)
  { return Horner<Float, Tag, order-1>::RecurseAlt(Polynomial<Float, Tag, order>::Coeff(), x); }
  //
  static Float Recurse(const Float term, const Float x, const Float y)
  { return Horner<Float, Tag, order-1>::Recurse(term*x + Polynomial<Float, Tag, order>::Coeff(y), x, y); }
  static Float Eval(const Float x, const Float y)
  { return Horner<Float, Tag, order-1>::Recurse(Polynomial<Float, Tag, order>::Coeff(y), x, y); }
};


template<typename Float, class Tag>
struct Horner<Float, Tag, 0> {
  static Float Recurse(const Float term, const Float x)
  { return term*x + Polynomial<Float, Tag, 0>::Coeff(); }
  static Float Eval(const Float x)
  { return Polynomial<Float, Tag, 0>::Coeff(); }
  //
  static Float Recurse(const Float term, const Float x, const Float y)
  { return term*x + Polynomial<Float, Tag, 0>::Coeff(y); }
  static Float Eval(const Float x, const Float y)
  { return Polynomial<Float, Tag, 0>::Coeff(y); }
};


#define HORNER_COEFF(_Tag_, _i_, _c_)  \
  template<typename Float> struct Polynomial<Float, _Tag_, _i_> { static Float Coeff() { return Float(_c_); } }
#define HORNER_COEFF2(_Tag_, _i_, _c_y_)  \
  template<typename Float> struct Polynomial<Float, _Tag_, _i_> { static Float Coeff(const Float y) { return Float(_c_y_); } }


#define HORNER0(F, x, c0)                                                (F)(c0)
#define HORNER1(F, x, c1, c0)                                 HORNER0(F, x, (F)(c1)*(x) + (F)(c0)                                )
#define HORNER2(F, x, c2, c1, c0)                             HORNER1(F, x, (F)(c2)*(x) + (F)(c1),                             c0)
#define HORNER3(F, x, c3, c2, c1, c0)                         HORNER2(F, x, (F)(c3)*(x) + (F)(c2),                         c1, c0)
#define HORNER4(F, x, c4, c3, c2, c1, c0)                     HORNER3(F, x, (F)(c4)*(x) + (F)(c3),                     c2, c1, c0)
#define HORNER5(F, x, c5, c4, c3, c2, c1, c0)                 HORNER4(F, x, (F)(c5)*(x) + (F)(c4),                 c3, c2, c1, c0)
#define HORNER6(F, x, c6, c5, c4, c3, c2, c1, c0)             HORNER5(F, x, (F)(c6)*(x) + (F)(c5),             c4, c3, c2, c1, c0)
#define HORNER7(F, x, c7, c6, c5, c4, c3, c2, c1, c0)         HORNER6(F, x, (F)(c7)*(x) + (F)(c6),         c5, c4, c3, c2, c1, c0)
#define HORNER8(F, x, c8, c7, c6, c5, c4, c3, c2, c1, c0)     HORNER7(F, x, (F)(c8)*(x) + (F)(c7),     c6, c5, c4, c3, c2, c1, c0)
#define HORNER9(F, x, c9, c8, c7, c6, c5, c4, c3, c2, c1, c0) HORNER8(F, x, (F)(c9)*(x) + (F)(c8), c7, c6, c5, c4, c3, c2, c1, c0)

#endif
