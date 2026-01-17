// ==========================================================================
// Copyright (C) 2011 Lior Kogan (koganlior1@gmail.com)
// ==========================================================================
// classes defined here:
// CircVal            - circular-value
// CircValTester      - tester for CircVal class
// ==========================================================================

// DRNadler 17-Jan-2026: Replace CircValTypeDef macro with CircValType template.

// DRNadler  2-Jan-2026: Deprecate float equality operators, and suppress GCC warning->error.

// DRNadler 31-Jan-2025: For better portability, replace M_PI with C++ 2020 std::numbers::pi

// DRNadler 31-Jan-2025: To avoid multiple defines when this header is included from multiple C++ source files:
// - static functions must be inline (compiler can decide whether or not to actually inline)
// - CircValTypeDef macro must use static constexpr to force compile-time evaluation and no storage allocation

#pragma once

#include <random>
#include <numbers>         // std::numbers::pi
#include <assert.h>

#include "FPCompare.h"
#include "CircHelper.h"   // Mod

// ==========================================================================
// use this template to define a circular-value type
template <double L_, double H_, double Z_>
struct CircValType
{
    // Prevent instances (no default ctor, no ctors at all in practice)
    CircValType() = delete;

    static constexpr double L = L_;            // range: [L,H)
    static constexpr double H = H_;
    static constexpr double Z = Z_;            // zero-value
    static constexpr double R = H_ - L_;       // range
    static constexpr double R_2 = (H_ - L_) / 2.0; // half range

    static_assert(H_ >  L_, "CircValType: invalid range (require H > L)");
    static_assert(Z_ >= L_, "CircValType: invalid zero (require Z >= L)");
    static_assert(Z_ <  H_, "CircValType: invalid zero (require Z < H)");
};

// define basic circular-value types
using   SignedDegRange = CircValType<           -180.0,                180.0, 0.0>;
using UnsignedDegRange = CircValType<              0.0,                360.0, 0.0>;
using   SignedRadRange = CircValType<-std::numbers::pi,     std::numbers::pi, 0.0>;
using UnsignedRadRange = CircValType<              0.0, 2 * std::numbers::pi, 0.0>;

// for testing only, define some additional circular-value types
using TestRange0 = CircValType<  3.0, 10.0,  5.3>;
using TestRange1 = CircValType< -3.0, 10.0, -3.0>;
using TestRange2 = CircValType< -3.0, 10.0,  9.9>;
using TestRange3 = CircValType<-13.0, -3.0, -5.3>;

// ==========================================================================
// circular value
// Type should be defined using the CircValTypeDef macro
template <typename Type>
class CircVal
{
    double val; // actual value [Type::L, Type::H)

    // ---------------------------------------------
public:
    inline static double GetL() { return Type::L; }
    inline static double GetH() { return Type::H; }
    inline static double GetZ() { return Type::Z; }
    inline static double GetR() { return Type::R; }

    // ---------------------------------------------
    inline static bool IsInRange(double r)
    {
        return (r>=Type::L && r<Type::H);
    }

    // 'wraps' circular-value to [Type::L,Type::H)
    inline static double Wrap(double r)
    {
        // the next lines are for optimization and improved accuracy only
        if (r>=Type::L)
        {
                 if (r< Type::H        ) return r        ;
            else if (r< Type::H+Type::R) return r-Type::R;
        }
        else
                 if (r>=Type::L-Type::R) return r+Type::R;

        // general case
        return Mod(r - Type::L, Type::R) + Type::L;
    }

    // ---------------------------------------------
    // the length of shortest directed walk from c1 to c2
    // return value is in [-Type::R/2, Type::R/2)
    inline static double Sdist(const CircVal& c1, const CircVal& c2)
    {
        double d= c2.val-c1.val;
        if (d <  -Type::R_2) { return d + Type::R; };
        if (d >=  Type::R_2) { return d - Type::R; };
                             { return d          ; };
    }

    // the length of the shortest increasing walk from c1 to c2
    // return value is in [0, Type::R)
    inline static double Pdist(const CircVal& c1, const CircVal& c2)
    {
        return c2.val>=c1.val ? c2.val-c1.val : Type::R-c1.val+c2.val;
    }

    // ---------------------------------------------
    CircVal() : val(Type::Z)
    {
    }

    // construction based on a floating-point value
    // floating-point is wrapped into the range
    // to translate a floating-point such that 0 is mapped to Type::Z, call ToC()
    CircVal(double r) : val(Wrap(r))
    {
    }

    // construction based on a circular value of the same type
    CircVal(const CircVal& c) : val(c.val)
    {
    }

    // construction based on a circular value of another type
    // sample use: CircVal<SignedRadRange> c= c2;   -or-   CircVal<SignedRadRange> c(c2);
    template<typename Type2>
    CircVal(const CircVal<Type2>& c) : val(Wrap(c.Pdist(c.GetZ(), c) * Type::R/Type2::R + Type::Z))
    {
    }

    // ---------------------------------------------
    operator double() const
    {
        return val;
    }

    // ---------------------------------------------
    // assignment from a floating-point value
    // floating-point is wrapped into the range
    // to translate a floating-point such that 0 is mapped to Type::Z, call ToC()
    CircVal& operator= (double r)
    {
        val= Wrap(r);
        return *this;
    }

    // assignment from another type of circular value
    template<typename Type2>
    CircVal& operator= (const CircVal<Type2>& c)
    {
        val= Wrap(c.Pdist(c.GetZ(), c) * Type::R/Type2::R + Type::Z);
        return *this;
    }

    // ---------------------------------------------
    // convert circular-value c to real-value [L-Z,H-Z). Z is converted to 0
    friend double ToR(const CircVal& c) { return c.val - Type::Z; }

    // ---------------------------------------------
    const CircVal  operator+ (                ) const { return val;                                            }
    const CircVal  operator- (                ) const { return Wrap(Type::Z-Sdist(Type::Z,val));               } // return negative circular value
    const CircVal  operator~ (                ) const { return Wrap(val+Type::R_2             );               } // return opposite circular-value

    const CircVal  operator+ (const CircVal& c) const { return Wrap(val+c.val        - Type::Z);               }
    const CircVal  operator- (const CircVal& c) const { return Wrap(val-c.val        + Type::Z);               }
    const CircVal  operator* (const double&  r) const { return Wrap((val-Type::Z)*r  + Type::Z);               }
    const CircVal  operator/ (const double&  r) const { return Wrap((val-Type::Z)/r  + Type::Z);               }

          CircVal& operator+=(const CircVal& c)       { val= Wrap(val+c.val          - Type::Z); return *this; }
          CircVal& operator-=(const CircVal& c)       { val= Wrap(val-c.val          + Type::Z); return *this; }
          CircVal& operator*=(const double&  r)       { val= Wrap((val-Type::Z)*r    + Type::Z); return *this; }
          CircVal& operator/=(const double&  r)       { val= Wrap((val-Type::Z)/r    + Type::Z); return *this; }

          CircVal& operator =(const CircVal& c)       { val= c.val                             ; return *this; }
    #ifdef __GNUC__
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wfloat-equal"
    #endif
    [[deprecated("Equality comparison of float types is risky; only use if checking exact 0 or identically calcd values.")]]
    bool           operator==(const CircVal& c) const { return val == c.val;                                   }
    [[deprecated("Equality comparison of float types is risky; only use if checking exact 0 or identically calcd values.")]]
    bool           operator!=(const CircVal& c) const { return val != c.val;                                   }
    #ifdef __GNUC__
        #pragma GCC diagnostic pop
    #endif

    // note that two circular values can be compared in several different ways.
    // check carefully if this is really what you need!
    bool           operator> (const CircVal& c) const { return val >  c.val;                                   }
    bool           operator>=(const CircVal& c) const { return val >= c.val;                                   }
    bool           operator< (const CircVal& c) const { return val <  c.val;                                   }
    bool           operator<=(const CircVal& c) const { return val <= c.val;                                   }
};

// ==========================================================================
template <typename Type> static double        sin  (const CircVal<Type>& c) { return std::sin(ToR(CircVal<SignedRadRange>(c)));  }
template <typename Type> static double        cos  (const CircVal<Type>& c) { return std::cos(ToR(CircVal<SignedRadRange>(c)));  }
template <typename Type> static double        tan  (const CircVal<Type>& c) { return std::tan(ToR(CircVal<SignedRadRange>(c)));  }
template <typename Type> static CircVal<Type> asin (double r              ) { return CircVal<SignedRadRange>(std::asin (r    )); } // calls copy ctor CircVal(CircVal<SignedRadRange>)
template <typename Type> static CircVal<Type> acos (double r              ) { return CircVal<SignedRadRange>(std::acos (r    )); } // calls copy ctor CircVal(CircVal<SignedRadRange>)
template <typename Type> static CircVal<Type> atan (double r              ) { return CircVal<SignedRadRange>(std::atan (r    )); } // calls copy ctor CircVal(CircVal<SignedRadRange>)
template <typename Type> static CircVal<Type> atan2(double r1, double r2  ) { return CircVal<SignedRadRange>(std::atan2(r1,r2)); } // calls copy ctor CircVal(CircVal<SignedRadRange>)
template <typename Type> static CircVal<Type> ToC  (double r              ) { return CircVal<Type>::Wrap(r + Type::Z);           } // convert real-value r to circular-value in the range. 0 is converted to Type.Z

// ==========================================================================
// tester for CircVal class
template <typename Type>
class CircValTester
{
    // check if 2 circular-values are almost equal
    inline static bool IsCircAlmostEq(const CircVal<Type>& c1, const CircVal<Type>& c2)
    {
        double r1= c1;
        double r2= c2;

        if (::IsAlmostEq(r1, r2))
            return true;

        if (r1 < r2)
            return IsAlmostEq(r1, r2 - Type::R);
        else
            return IsAlmostEq(r1, r2 + Type::R);
    }

    // assert that 2 circular-values are almost equal
    inline static void AssertCircAlmostEq(const CircVal<Type>& c1, const CircVal<Type>& c2)
    {
        assert(IsCircAlmostEq(c1, c2));
    }

    inline static void Test()
    {
        CircVal<Type> ZeroVal= Type::Z;

        // --------------------------------------------------------
        AssertCircAlmostEq(ZeroVal       , -ZeroVal);

        AssertAlmostEq    (sin(ZeroVal)  , 0.      );
        AssertAlmostEq    (cos(ZeroVal)  , 1.      );
        AssertAlmostEq    (tan(ZeroVal)  , 0.      );

        AssertCircAlmostEq(asin<Type>(0.), ZeroVal );
        AssertCircAlmostEq(acos<Type>(1.), ZeroVal );
        AssertCircAlmostEq(atan<Type>(0.), ZeroVal );

        AssertCircAlmostEq(ToC<Type>(0)  , ZeroVal );
        AssertAlmostEq    (ToR(ZeroVal)  , 0.      );

        // --------------------------------------------------------
        std::default_random_engine             rand_engine                 ;
        std::uniform_real_distribution<double> c_uni_dist(Type::L, Type::H);
        std::uniform_real_distribution<double> r_uni_dist(0.     , 1000.  ); // for multiplication,division by real-value
        std::uniform_real_distribution<double> t_uni_dist(-1.    , 1.     ); // for inverse-trigonometric functions

        std::random_device rnd_device;
        rand_engine.seed(rnd_device()); // reseed engine

        for (unsigned i= 10000; i--;)
        {
            CircVal<Type> c1(c_uni_dist(rand_engine)); // random circular value
            CircVal<Type> c2(c_uni_dist(rand_engine)); // random circular value
            CircVal<Type> c3(c_uni_dist(rand_engine)); // random circular value
            double        r (r_uni_dist(rand_engine)); // random real     value [    0, 1000) - for testing *,/ operators
            double        a1(t_uni_dist(rand_engine)); // random real     value [   -1,    1) - for testing asin,acos
            double        a2(t_uni_dist(rand_engine)); // random real     value [-1000, 1000) - for testing atan

            assert            (c1                                 == CircVal<Type>((double)c1)         );

            AssertCircAlmostEq(+c1                                  , c1                               ); // +c         = c
            AssertCircAlmostEq(-(-c1)                               , c1                               ); // -(-c)      = c
            AssertCircAlmostEq(c1 + c2                              , c2 + c1                          ); // c1+c2      = c2+c1
            AssertCircAlmostEq(c1 + (c2 +c3)                        , (c1 + c2) + c3                   ); // c1+(c2+c3) = (c1+c2)+c3
            AssertCircAlmostEq(c1 + -c1                             , ZeroVal                          ); // c+(-c)     = z
            AssertCircAlmostEq(c1 + ZeroVal                         , c1                               ); // c+z        = c

            AssertCircAlmostEq(c1      - c1                         , ZeroVal                          ); // c-c        = z
            AssertCircAlmostEq(c1      - ZeroVal                    , c1                               ); // c-z        = c
            AssertCircAlmostEq(ZeroVal - c1                         , -c1                              ); // z-c        = -c
            AssertCircAlmostEq(c1      - c2                         , -(c2 - c1)                       ); // c1-c2      = -(c2-c1)

            AssertCircAlmostEq(c1 * 0.                              , ZeroVal                          ); // c*0        = 0
            AssertCircAlmostEq(c1 * 1.                              , c1                               ); // c*1        = c
            AssertCircAlmostEq(c1 / 1.                              , c1                               ); // c/1        = c

            AssertCircAlmostEq((c1 * (1./(r+1.))) / (1./(r+1.))     , c1                               ); // (c*r)/r    = c, 0<r<=1
            AssertCircAlmostEq((c1 / (    r+1.) ) * (    r+1. )     , c1                               ); // (c/r)*r    = c,   r>=1

            // --------------------------------------------------------
            AssertCircAlmostEq(~(~c1)                               , c1                               ); // opposite(opposite(c) = c
            AssertCircAlmostEq(c1 - (~c1)                           , ToC<Type>(Type::R/2.)            ); // c - ~c               = r/2+z

            // --------------------------------------------------------
            AssertAlmostEq    (sin(ToR(CircVal<SignedRadRange>(c1))),  sin(c1)                         ); // member func sin
            AssertAlmostEq    (cos(ToR(CircVal<SignedRadRange>(c1))),  cos(c1)                         ); // member func cos
            AssertAlmostEq    (tan(ToR(CircVal<SignedRadRange>(c1))),  tan(c1)                         ); // member func tan

            AssertAlmostEq    (sin(-c1)                             , -sin(c1)                         ); // sin(-c)    = -sin(c)
            AssertAlmostEq    (cos(-c1)                             ,  cos(c1)                         ); // cos(-c)    =  cos(c)
            AssertAlmostEq    (tan(-c1)                             , -tan(c1)                         ); // tan(-c1)   = -tan(c) the error may be large

            AssertAlmostEq    (sin(c1+ToC<Type>(Type::R/4.))        ,  cos(c1)                         ); // sin(c+r/4) =  cos(c)
            AssertAlmostEq    (cos(c1+ToC<Type>(Type::R/4.))        , -sin(c1)                         ); // cos(c+r/4) = -sin(c)
            AssertAlmostEq    (sin(c1+ToC<Type>(Type::R/2.))        , -sin(c1)                         ); // sin(c+r/2) = -sin(c)
            AssertAlmostEq    (cos(c1+ToC<Type>(Type::R/2.))        , -cos(c1)                         ); // cos(c+r/2) = -cos(c)

            AssertAlmostEq    (Sqr(sin(c1))+Sqr(cos(c1))            , 1.                               ); // sin(x)^2+cos(x)^2 = 1

            AssertAlmostEq    (sin(c1)/cos(c1)                      , tan(c1)                          ); // sin(x)/cos(x) = tan(x)

            // --------------------------------------------------------
            AssertCircAlmostEq(asin<Type>(a1)                       , CircVal<SignedRadRange>(asin(a1))); // member func asin
            AssertCircAlmostEq(acos<Type>(a1)                       , CircVal<SignedRadRange>(acos(a1))); // member func acos
            AssertCircAlmostEq(atan<Type>(a2)                       , CircVal<SignedRadRange>(atan(a2))); // member func atan

            AssertCircAlmostEq(asin<Type>(a1) + asin<Type>(-a1)     , ZeroVal                          ); // asin(r)+asin(-r) = z
            AssertCircAlmostEq(acos<Type>(a1) + acos<Type>(-a1)     , ToC<Type>(Type::R/2.)            ); // acos(r)+acos(-r) = r/2+z
            AssertCircAlmostEq(asin<Type>(a1) + acos<Type>( a1)     , ToC<Type>(Type::R/4.)            ); // asin(r)+acos( r) = r/4+z
            AssertCircAlmostEq(atan<Type>(a2) + atan<Type>(-a2)     , ZeroVal                          ); // atan(r)+atan(-r) = z

            // --------------------------------------------------------
            assert            (c1 >  c2                           ==    (c2 <  c1)                     ); // c1> c2 <==>   c2< c1
            assert            (c1 >= c2                           ==    (c2 <= c1)                     ); // c1>=c2 <==>   c2<=c1
            assert            (c1 >= c2                           ==  ( (c1 >  c2) ||  (c1 == c2))     ); // c1>=c2 <==>  (c1> c2)|| (c1==c2)
            assert            (c1 <= c2                           ==  ( (c1 <  c2) ||  (c1 == c2))     ); // c1<=c2 <==>  (c1< c2)|| (c1==c2)
            assert            (c1 >  c2                           ==  (!(c1 == c2) && !(c1 <  c2))     ); // c1> c2 <==> !(c1==c2)&&!(c1< c2)
            assert            (c1 == c2                           ==  (!(c1 >  c2) && !(c1 <  c2))     ); // c1= c2 <==> !(c1> c2)&&!(c1< c2)
            assert            (c1 <  c2                           ==  (!(c1 == c2) && !(c1 >  c2))     ); // c1< c2 <==> !(c1==c2)&&!(c1> c2)
            assert            (!(c1>c2) || !(c2>c3) || (c1>c3)                                         ); // (c1>c2)&&(c2>c3) ==> c1>c3

            // --------------------------------------------------------
            AssertCircAlmostEq(c1                                   , ToC<Type>(ToR( c1)       )       ); //  c1        = ToC(ToR( c1)
            AssertCircAlmostEq(-c1                                  , ToC<Type>(ToR(-c1)       )       ); // -c1        = ToC(ToR(-c1)
            AssertCircAlmostEq(c1 + c2                              , ToC<Type>(ToR(c1)+ToR(c2))       ); // c1+c2      = ToC(ToR(c1)+ToR(c2))
            AssertCircAlmostEq(c1 - c2                              , ToC<Type>(ToR(c1)-ToR(c2))       ); // c1-c2      = ToC(ToR(c1)-ToR(c2))
            AssertCircAlmostEq(c1 * r                               , ToC<Type>(ToR(c1)*r      )       ); // c1*r       = ToC(ToR(c1)*r      )
            AssertCircAlmostEq(c1 / r                               , ToC<Type>(ToR(c1)/r      )       ); // c1/r       = ToC(ToR(c1)/r      )

            // --------------------------------------------------------
        }
    }

public:
    CircValTester()
    {
        Test();
    }
};
