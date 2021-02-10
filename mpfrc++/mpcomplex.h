/*
 * Copyright (c) 2010-2012
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */
/*
Complex class declare for the MPFR
*/

#ifndef __MP_COMPLEX_H__
#define __MP_COMPLEX_H__

#include "mpreal.h"
#include <complex>
#include "mpc.h"

#if defined ___MPACK_BUILD_WITH_GMP___
#include "gmpxx.h"
#include "mpc_class.h"
#endif
#if defined ___MPACK_BUILD_WITH_QD___
#include "qd_complex.h"
#endif
#if defined ___MPACK_BUILD_WITH_DD___
#include "dd_complex.h"
#endif

namespace mpfr {
class mpcomplex {
  private:
   mpc_t mpc;

  public:
  static mpc_rnd_t default_rnd;
  static mp_prec_t default_real_prec;
  static mp_prec_t default_imag_prec;
  static int default_base;
  static int double_bits;

//constructor & deconstructor
    mpcomplex();
    mpcomplex(const mpc_t a);
    mpcomplex(const mpfr_t a, const mpfr_t b);
    mpcomplex(const mpf_t a, const mpf_t b);
    mpcomplex(const char *s, mp_prec_t pr = default_real_prec, mp_prec_t pi = default_imag_prec, mpc_rnd_t mode = default_rnd);
    mpcomplex(const mpcomplex& a);
    mpcomplex(const std::complex<double>& a, mp_prec_t pr = default_real_prec, mp_prec_t pi = default_imag_prec, mpc_rnd_t mode = default_rnd);
    mpcomplex(const mpreal& a, const mpreal& b);
    mpcomplex(const double& a, const double& b, mp_prec_t pr = default_real_prec, mp_prec_t pi = default_imag_prec, mpc_rnd_t mode = default_rnd);
    mpcomplex(const char *s, const char *t, mp_prec_t pr = default_real_prec, mp_prec_t pi = default_imag_prec, mpc_rnd_t mode = default_rnd);

    mpcomplex(const mpreal & a);
    mpcomplex(const mpfr_t a);
    mpcomplex(const mpf_t a);
    mpcomplex(const double a, mp_prec_t pr = default_real_prec, mp_prec_t pi = default_imag_prec, mpc_rnd_t mode = default_rnd);
    ~mpcomplex(); 

    mpcomplex& operator=(const mpcomplex& a);
    mpcomplex& operator=(const mpc_t a);
    mpcomplex& operator=(const std::complex<double> a);
    mpcomplex& operator=(const char* s);

    mpcomplex& operator=(const mpreal& a);
    mpcomplex& operator=(const mpfr_t a);
    mpcomplex& operator=(const double a);

    //+
    mpcomplex& operator+=(const mpcomplex& a);
    mpcomplex& operator+=(const mpc_t a);
    mpcomplex& operator+=(const std::complex<double> a);

    mpcomplex& operator+=(const mpreal& a);
    mpcomplex& operator+=(const mpfr_t a);
    mpcomplex& operator+=(const double a);
    const mpcomplex operator+() const;
    mpcomplex& operator++();
    const mpcomplex operator ++ (int);

    //-
    mpcomplex& operator-=(const mpcomplex& a);
    mpcomplex& operator-=(const mpc_t a);
    mpcomplex& operator-=(const std::complex<double> a);

    mpcomplex& operator-=(const mpreal& a);
    mpcomplex& operator-=(const mpfr_t a);
    mpcomplex& operator-=(const double a);
    const mpcomplex operator-() const;
    mpcomplex& operator--();
    const mpcomplex operator -- (int);

    //*
    mpcomplex& operator*=(const mpcomplex& a);
    mpcomplex& operator*=(const mpc_t a);
    mpcomplex& operator*=(const std::complex<double> a);
    mpcomplex& operator*=(const mpreal& a);
    mpcomplex& operator*=(const mpfr_t a);
    mpcomplex& operator*=(const double a);

    // /
    mpcomplex& operator/=(const mpcomplex& a);
    mpcomplex& operator/=(const mpc_t a);
    mpcomplex& operator/=(const std::complex<double> a);
    mpcomplex& operator/=(const mpreal& a);
    mpcomplex& operator/=(const mpfr_t a);
    mpcomplex& operator/=(const double a);

    //comparison
    friend bool operator == (const mpcomplex& a, const mpcomplex& b);
    friend bool operator == (const mpcomplex& a, const mpreal& b);
    friend bool operator == (const mpreal& a, const mpcomplex& b);
    friend bool operator != (const mpcomplex& a, const mpcomplex& b);
    friend bool operator != (const mpreal& a, const mpcomplex& b);
    friend bool operator != (const mpcomplex& a, const mpreal& b);

    //random
    friend const mpcomplex urandom_c (gmp_randstate_t& state);

    inline mpreal real() //this should not be call by reference, as mpreal is a class contains only pointers.
    {
        mpreal tmp;
        tmp = mpc_realref(mpc);
        return tmp;
    }
    inline const mpreal real() const //this should not be call by reference, as mpreal is a class contains only pointers. 
    {
        mpreal tmp;
        tmp = mpc_realref(mpc);
        return tmp;
    }
    inline mpreal imag() //this should not be call by reference, as mpreal is a class contains only pointers. 
    {
        mpreal tmp;
        tmp = mpc_imagref(mpc);
        return tmp;
    }
    inline const mpreal imag() const //this should not be call by reference, as mpreal is a class contains only pointers. 
    {
        mpreal tmp;
        tmp = mpc_imagref(mpc);
        return tmp;
    }
    inline void real(const mpreal r) //constructor
    {
        mpreal tmp(r);  //required as r is const.
        mpc_real(mpfr_ptr(tmp), mpc, mpreal::default_rnd);
    }
    inline void imag(const mpreal r) //constructor
    {
        mpreal tmp(r);
        mpc_imag(mpfr_ptr(tmp), mpc, mpreal::default_rnd);
    }
    //other functions
    friend const mpreal abs(const mpcomplex& a, mpfr_rnd_t mode);
    friend const mpreal norm(const mpcomplex& a, mpfr_rnd_t mode);
    friend const mpcomplex conj(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpreal arg(const mpcomplex& a, mpfr_rnd_t mode);
    friend const mpcomplex proj(const mpcomplex& a, mpc_rnd_t mode);

    //Powerfunction and Logarithm    
    friend const mpcomplex sqrt(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex pow(const mpcomplex& a, const mpcomplex& b,mpc_rnd_t mode);
    friend const mpcomplex exp(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex log(const mpcomplex& a, mpc_rnd_t mode);

    //trigonometric functions
    friend const mpcomplex sin(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex cos(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex tan(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex sinh(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex cosh(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex tanh(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex asin(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex acos(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex atan(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex asinh(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex acosh(const mpcomplex& a, mpc_rnd_t mode);
    friend const mpcomplex atanh(const mpcomplex& a, mpc_rnd_t mode);

// Type Conversion operators
    operator std::complex<double>() const;
    operator std::string() const;
    std::string to_string(size_t n = 0, int b = default_base, mpc_rnd_t mode = default_rnd) const;
// Set/Get instance properties
   static void	set_default_prec(mp_prec_t prec);
   mp_prec_t get_prec() const;
   mp_prec_t get_prec_re() const;
   mp_prec_t get_prec_im() const;
   void set_prec(mp_prec_t prec, mpc_rnd_t rnd_mode);
   void set_prec2(mp_prec_t pr, mp_prec_t pi, mpc_rnd_t rnd_mode);

#if defined ___MPACK_BUILD_WITH_GMP___
mpcomplex(const mpc_class& a);
mpcomplex& operator=(const mpc_class& a);
#endif
#if defined ___MPACK_BUILD_WITH_QD___
mpcomplex(const qd_complex& a);
mpcomplex& operator=(const qd_complex& a);
#endif
#if defined ___MPACK_BUILD_WITH_DD___
mpcomplex(const dd_complex& a);
mpcomplex& operator=(const dd_complex& a);
#endif
#if defined ___MPACK_BUILD_WITH___FLOAT128___
mpcomplex(const std::complex<__float128>& a);
mpcomplex& operator=(const std::complex<__float128>& a);
#endif
#if defined ___MPACK_BUILD_WITH___DOUBLE___
mpcomplex(const std::complex<double>& a);
mpcomplex& operator=(const std::complex<double>& a);
#endif
};

//+ addition
const mpcomplex operator+(const mpcomplex& a, const mpcomplex& b);
const mpcomplex operator+(const mpcomplex& a, const std::complex<double> b);
const mpcomplex operator+(const mpcomplex& a, const char* b);

const mpcomplex operator+(const std::complex<double> a, const mpcomplex& b);
const mpcomplex operator+(const char* a, const mpcomplex& b);

const mpcomplex operator+(const mpcomplex& a, const mpreal b);
const mpcomplex operator+(const mpcomplex& a, const double b);
const mpcomplex operator+(const mpcomplex& a, const int b);

const mpcomplex operator+(const mpreal a, const mpcomplex& b);
const mpcomplex operator+(const double a, const mpcomplex& b);
const mpcomplex operator+(const int a, const mpcomplex& b);
//- subtraction
const mpcomplex operator-(const mpcomplex& a, const mpcomplex& b);
const mpcomplex operator-(const mpcomplex& a, const std::complex<double>& b);
const mpcomplex operator-(const mpcomplex& a, const char* b);

const mpcomplex operator-(const std::complex<double>& a, const mpcomplex& b);
const mpcomplex operator-(const char* a, const mpcomplex& b);

const mpcomplex operator-(const mpcomplex& a, const mpreal b);
const mpcomplex operator-(const mpcomplex& a, const double b);
const mpcomplex operator-(const mpcomplex& a, const int b);

const mpcomplex operator-(const mpreal a, const mpcomplex& b);
const mpcomplex operator-(const double a, const mpcomplex& b);
const mpcomplex operator-(const int a, const mpcomplex& b);

//* multiplication
const mpcomplex operator*(const mpcomplex& a, const mpcomplex& b);
const mpcomplex operator*(const mpcomplex& a, const mpreal& b);
const mpcomplex operator*(const mpreal& a, const mpcomplex& b);

/// division
const mpcomplex operator/(const mpcomplex& a, const mpcomplex& b);
const mpcomplex operator/(const mpcomplex& a, const mpreal& b);
const mpcomplex operator/(const mpreal& a, const mpcomplex& b);

inline void mpcomplex::set_default_prec(mp_prec_t prec)
{ 
	default_real_prec = prec;
	default_imag_prec = prec;
}

inline mp_prec_t mpcomplex::get_prec() const
{
	return mpc_get_prec(mpc);
}

inline mp_prec_t mpcomplex::get_prec_re() const
{
      mp_prec_t pr, pi;
      mpc_get_prec2(&pr, &pi, mpc);
      return pr;
}
inline mp_prec_t mpcomplex::get_prec_im() const
{
      mp_prec_t pr, pi;
      mpc_get_prec2(&pr, &pi, mpc);
      return pi;
}

inline void mpcomplex::set_prec(mp_prec_t prec, mpc_rnd_t rnd_mode)
{
    mpcomplex tmp (*this);
    mpc_init2(mpc, prec);
    mpc_set(mpc, tmp.mpc, default_rnd);
}

inline void mpcomplex::set_prec2(mp_prec_t pr, mp_prec_t pi, mpc_rnd_t rnd_mode)
{
    mpcomplex tmp (*this);
    mpc_init3(mpc, pr, pi);
    mpc_set(mpc, tmp.mpc, default_rnd);
}

//constructor and deconstructor
inline mpcomplex::mpcomplex()
{
   mpc_init3(mpc, default_real_prec, default_imag_prec);
   mpc_set_ui(mpc,0UL,default_rnd);
}

inline mpcomplex::mpcomplex(const mpc_t a) 
{
   mp_prec_t pr, pi;
   mpc_get_prec2(&pr, &pi, a);
   mpc_init3(mpc, pr, pi);
   mpc_set(mpc, a, default_rnd);
}

inline mpcomplex::mpcomplex(const mpfr_t a, const mpfr_t b)
{
   mp_prec_t pr, pi;
   pr = mpfr_get_prec(a);  pi = mpfr_get_prec(b);
   mpc_init3(mpc, pr, pi);
   mpc_set_fr_fr(mpc, a, b, default_rnd);
}

inline mpcomplex::mpcomplex(const mpf_t a, const mpf_t b)
{
   mp_prec_t pr, pi;
   pr = mpf_get_prec(a);  pi = mpf_get_prec(b);
   mpc_init3(mpc, pr, pi);
   mpc_set_f_f(mpc, a, b, default_rnd);
}

inline mpcomplex::mpcomplex(const char *s, mp_prec_t pr, mp_prec_t pi, mpc_rnd_t mode)
{
    mpc_init3(mpc, pr, pi);
    mpc_set_str(mpc, (char *)s, default_base, mode); 
}

inline mpcomplex::mpcomplex(const mpcomplex& a) 
{
   mp_prec_t pr, pi;
   mpc_get_prec2(&pr, &pi, a.mpc);
   mpc_init3(mpc, pr, pi);
   mpc_set(mpc, a.mpc, default_rnd);
}

inline mpcomplex::mpcomplex(const std::complex<double>& a, mp_prec_t pr, mp_prec_t pi, mpc_rnd_t mode)
{
   mpc_init3(mpc, pr, pi);
   mpc_set_d_d(mpc, a.real(), a.imag(), default_rnd);
}

inline mpcomplex::mpcomplex(const mpreal& a, const mpreal& b)
{
   mp_prec_t pr, pi;
   mpreal tmp1(a), tmp2(b);
   pr = a.get_prec();  pi = b.get_prec();
   mpc_init3(mpc, pr, pi);
   mpc_set_fr_fr(mpc, (mpfr_ptr)tmp1, (mpfr_ptr)tmp2, default_rnd);
}

inline mpcomplex::mpcomplex(const double& a, const double& b, mp_prec_t pr, mp_prec_t pi, mpc_rnd_t mode)
{
   mpc_init3(mpc, default_real_prec, default_imag_prec);
   mpc_set_d_d(mpc, a, b,default_rnd);
}

inline mpcomplex::mpcomplex(const char *s, const char *t, mp_prec_t pr, mp_prec_t pi, mpc_rnd_t mode)
{
   mpfr_t a, b;
   mpfr_init2(a, pr); mpfr_init2(b, pi);
   mpfr_set_str(a, s, default_base, MPC_RND_RE(mode));
   mpfr_set_str(b, t, default_base, MPC_RND_IM(mode));
   mpc_init3(mpc, pr, pi);
   mpc_set_fr_fr(mpc, a, b, default_rnd);
   mpfr_clear(a); mpfr_clear(b);
}

inline mpcomplex::mpcomplex(const mpreal& a)
{
   mpreal tmp(a);
   mpc_init3(mpc, default_real_prec, default_imag_prec);
   mpc_set_fr(mpc, (mpfr_ptr)(tmp), default_rnd);
}

inline mpcomplex::mpcomplex(const mpfr_t a)
{
   mp_prec_t pr;
   pr = mpfr_get_prec(a);
   mpc_init2(mpc, pr);
   mpc_set_fr(mpc, a, default_rnd);
}

inline mpcomplex::mpcomplex(const mpf_t a)
{
   mp_prec_t pr;
   pr = mpf_get_prec(a);
   mpc_init2(mpc, pr);
   mpc_set_f(mpc, a, default_rnd);
}

inline mpcomplex::mpcomplex(const double a, mp_prec_t pr, mp_prec_t pi, mpc_rnd_t mode)
{
   mpc_init3(mpc, pr, pi);
   mpc_set_d(mpc, a, default_rnd);
}

inline mpcomplex::~mpcomplex() 
{ 
   mpc_clear(mpc);
}                           

inline mpcomplex& mpcomplex::operator=(const mpcomplex& a)
{
	if (this!= &a) mpc_set(mpc,a.mpc,default_rnd);
	return *this;
}

inline mpcomplex& mpcomplex::operator=(const mpc_t a)
{
	mpc_set(mpc,a,default_rnd);
	return *this;
}

inline mpcomplex& mpcomplex::operator=(const std::complex<double> a)
{
    mpc_set_d_d(mpc, a.real(), a.imag(), default_rnd);
    return *this;
}

inline mpcomplex& mpcomplex::operator=(const char* s)
{
     mpc_init3(mpc, default_real_prec, default_imag_prec);
     mpc_set_str(mpc, (char *)s,default_base,default_rnd);				
     return *this;
}

inline mpcomplex& mpcomplex::operator=(const mpreal& a)
{
   mpc_set_fr(mpc, (mpfr_ptr)(&a), default_rnd);
   return *this;
}

inline mpcomplex& mpcomplex::operator=(const mpfr_t a)
{
   mpc_set_fr(mpc, a, default_rnd);
   return *this;
}

inline mpcomplex& mpcomplex::operator=(const double a)
{
   mpc_set_d(mpc, a, default_rnd);
   return *this;
}
// + Addition
inline mpcomplex& mpcomplex::operator+=(const mpcomplex& a)
{
	mpc_add(mpc,mpc,a.mpc,default_rnd);
	return *this;
}

inline mpcomplex& mpcomplex::operator+=(const mpc_t a)
{
	*this += mpcomplex(a);
	return *this;
}

inline mpcomplex& mpcomplex::operator+=(const std::complex<double> a)
{
	return *this += mpcomplex(a);	
}

inline mpcomplex& mpcomplex::operator+=(const mpreal& a)
{
     mpc_add_fr(mpc, mpc, (mpfr_ptr)(&a), default_rnd);
     return *this;
}

inline mpcomplex& mpcomplex::operator+=(const mpfr_t a)
{
     mpc_add_fr(mpc, mpc, a, default_rnd);
     return *this;
}


inline mpcomplex& mpcomplex::operator+=(const double a)
{
     mpreal p = a;
     mpc_add_fr(mpc, mpc, (mpfr_ptr)(&p), default_rnd);
     return *this;
}

inline const mpcomplex mpcomplex::operator+()const
{
	return mpcomplex(*this);
}

inline mpcomplex& mpcomplex::operator++() 
{
	*this += 1.0;
	return *this;
}

inline const mpcomplex mpcomplex::operator++ (int)
{
	mpcomplex x(*this);
	*this += 1.0;
	return x;
}

// + Addition again
inline const mpcomplex operator+(const mpcomplex& a, const mpcomplex& b)
{
  mpcomplex tmp;
  if (!a.get_prec() == 0 && !b.get_prec() == 0) {
      if(a.get_prec()>b.get_prec()) { tmp=a; tmp += b; return tmp;} 
      else {tmp=b; tmp += a; return tmp;} 
  } else {
    tmp=a;
    tmp.set_prec2(std::max(a.get_prec_re(), b.get_prec_re()), std::max(a.get_prec_im(), b.get_prec_im()), mpcomplex::default_rnd);
    return tmp += b;
  }
}

inline const mpcomplex operator+(const mpcomplex& a, const std::complex<double> b)
{
  return mpcomplex(a) += b;
}

inline const mpcomplex operator+(const mpcomplex& a, const char *b)
{
  return mpcomplex(b) += a;
}

inline const mpcomplex operator+(const std::complex<double> a, const mpcomplex& b)
{
  return mpcomplex(b) += a;
}

inline const mpcomplex operator+(const char *a, const mpcomplex& b)
{
  return mpcomplex(a) += b;
}


inline const mpcomplex operator+(const mpcomplex& a, const mpreal b)
{
  mpcomplex tmp(a);
  tmp.set_prec2(std::max(a.get_prec_re(), b.get_prec()), a.get_prec_im(), mpcomplex::default_rnd);
  return tmp += b;
}

inline const mpcomplex operator+(const mpreal a, const mpcomplex& b)
{
  mpcomplex tmp(b);
  tmp.set_prec2(std::max(b.get_prec_re(), a.get_prec()), b.get_prec_im(), mpcomplex::default_rnd);
  return tmp += a;
}

inline const mpcomplex operator+(const double a, const mpcomplex& b)
{
  mpcomplex tmp(b);
  return tmp += a;
}

inline const mpcomplex operator+(const int a, const mpcomplex& b)
{
  mpcomplex tmp(b);
  return tmp += a;
}

// - Subtraction
inline mpcomplex& mpcomplex::operator-=(const mpcomplex& a)
{
        mpc_sub(mpc,mpc,a.mpc,default_rnd);
        return *this;
}

inline mpcomplex& mpcomplex::operator-=(const mpc_t a)
{
	*this -= mpcomplex(a);
	return *this;
}

inline mpcomplex& mpcomplex::operator-=(const std::complex<double> a)
{
	return *this -= mpcomplex(a);	
}

inline mpcomplex& mpcomplex::operator-=(const mpreal& a)
{
     mpc_sub_fr(mpc, mpc, (mpfr_ptr)(&a), default_rnd);
     return *this;
}

inline mpcomplex& mpcomplex::operator-=(const mpfr_t a)
{
     mpc_sub_fr(mpc, mpc, a, default_rnd);
     return *this;
}


inline mpcomplex& mpcomplex::operator-=(const double a)
{
     mpreal p = a;
     mpc_sub_fr(mpc, mpc, (mpfr_ptr)(&p), default_rnd);
     return *this;
}

inline const mpcomplex mpcomplex::operator-()const
{
     mpcomplex u(*this);
     mpc_neg (u.mpc, u.mpc, default_rnd);
     return u;
}

inline mpcomplex& mpcomplex::operator--() 
{
	*this -= 1.0;
	return *this;
}

inline const mpcomplex mpcomplex::operator-- (int)
{
	mpcomplex x(*this);
	*this -= 1.0;
	return x;
}

inline const mpcomplex operator-(const mpcomplex& a, const mpcomplex& b)
{
  if (!a.get_prec() == 0 && !b.get_prec() == 0) {
    if(a.get_prec()>b.get_prec()) return mpcomplex(a) -= b;
    else return -(mpcomplex(b) -= a);
  } else {
    mpcomplex tmp(a);
    tmp.set_prec2(std::max(a.get_prec_re(), b.get_prec_re()), std::max(a.get_prec_im(), b.get_prec_im()), mpcomplex::default_rnd);
    return tmp -= b;
  }
}

inline const mpcomplex operator-(const mpcomplex& a, const std::complex<double>& b)
{
    return -(mpcomplex(b) -= a);
}

inline const mpcomplex operator-(const mpcomplex& a, const char* b)
{
       return a-mpcomplex(b);
}

inline const mpcomplex operator-(const std::complex<double>& a, const mpcomplex& b)
{
    return mpcomplex(a) -= b;
}

inline const mpcomplex operator-(const char* a, const mpcomplex& b)
{
        return mpcomplex(a)-b;
}

inline const mpcomplex operator-(const mpcomplex& a, const mpreal b)
{
  if (!a.get_prec() == 0 && !b.get_prec() == 0) {
    if(a.get_prec()>b.get_prec()) return mpcomplex(a) -= b;
    else return -(mpcomplex(b) -= a);
  } else {
    mpcomplex tmp(a);
    tmp.set_prec2(std::max(a.get_prec_re(), b.get_prec()), a.get_prec_im(), mpcomplex::default_rnd);
    return tmp -= b;
  }
}

inline const mpcomplex operator-(const mpcomplex& a, const double b)
{
        return mpcomplex(a)-mpreal(b);
}

inline const mpcomplex operator-(const mpcomplex& a, const int b)
{
        return mpcomplex(a)-mpreal(b);
}

inline const mpcomplex operator-(const mpreal a, const mpcomplex& b)
{
        return mpcomplex(a)-mpcomplex(b);
}

inline const mpcomplex operator-(const double a, const mpcomplex& b)
{
        return mpcomplex(a)-mpcomplex(b);
}

inline const mpcomplex operator-(const int a, const mpcomplex& b)
{
        return mpcomplex(a)-mpcomplex(b);
}

// * multiplication
inline mpcomplex& mpcomplex::operator*=(const mpcomplex& a)
{
	mpc_mul(mpc,mpc,a.mpc,default_rnd);
	return *this;
}

inline mpcomplex& mpcomplex::operator*=(const mpc_t a)
{
	*this *= mpcomplex(a);
	return *this;
}

inline mpcomplex& mpcomplex::operator*=(const std::complex<double> a)
{
	return *this *= mpcomplex(a);	
}

inline mpcomplex& mpcomplex::operator*=(const mpreal& a)
{
     mpc_mul_fr(mpc, mpc, (mpfr_ptr)(&a), default_rnd);
     return *this;
}

inline mpcomplex& mpcomplex::operator*=(const mpfr_t a)
{
     mpc_mul_fr(mpc, mpc, a, default_rnd);
     return *this;
}


inline mpcomplex& mpcomplex::operator*=(const double a)
{
     mpreal p(a);
     mpc_add_fr(mpc, mpc, (mpfr_ptr)(&p), default_rnd);
     return *this;
}

inline const mpcomplex operator*(const mpcomplex& a, const mpcomplex& b)
{
  if (!a.get_prec() == 0 && !b.get_prec() == 0) {
    if(a.get_prec()>b.get_prec()) return mpcomplex(a) *= b;
    else return mpcomplex(b) *= a;
  } else {
    mpcomplex tmp(a);
    tmp.set_prec2(std::max(a.get_prec_re(), b.get_prec_re()), std::max(a.get_prec_im(), b.get_prec_im()), mpcomplex::default_rnd);
    return tmp *= b;
  }
}

inline const mpcomplex operator*(const mpcomplex& a, const mpreal& b)
{
   mpcomplex tmp(a);
   tmp.set_prec2(std::max(a.get_prec_re(), b.get_prec()), a.get_prec_im(), mpcomplex::default_rnd);
   return tmp *= mpcomplex(b);
}

inline const mpcomplex operator*(const mpreal& a, const mpcomplex& b)
{
   mpcomplex tmp(a);
   tmp.set_prec2(std::max(b.get_prec_re(), a.get_prec()), b.get_prec_im(), mpcomplex::default_rnd);
   return tmp *= b;
}

// / division
inline mpcomplex& mpcomplex::operator/=(const mpcomplex& a)
{
	mpc_div(mpc,mpc,a.mpc,default_rnd);
	return *this;
}

inline mpcomplex& mpcomplex::operator/=(const mpc_t a)
{
	*this /= mpcomplex(a);
	return *this;
}

inline mpcomplex& mpcomplex::operator/=(const std::complex<double> a)
{
	return *this /= mpcomplex(a);	
}

inline mpcomplex& mpcomplex::operator/=(const mpreal& a)
{
     mpc_div_fr(mpc, mpc, (mpfr_ptr)(&a), default_rnd);
     return *this;
}

inline mpcomplex& mpcomplex::operator/=(const mpfr_t a)
{
     mpc_div_fr(mpc, mpc, a, default_rnd);
     return *this;
}

inline mpcomplex& mpcomplex::operator/=(const double a)
{
     mpreal p(a);
     mpc_div_fr(mpc, mpc, (mpfr_ptr)(&p), default_rnd);
     return *this;
}

inline bool operator == (const mpcomplex& a, const mpcomplex& b)
{
	return (mpc_cmp(a.mpc,b.mpc)==0);
}

inline bool operator == (const mpcomplex& a, const mpreal& b)
{
        mpcomplex c(b);
	return (mpc_cmp(a.mpc, c.mpc)==0);
}

inline bool operator == (const mpreal& a, const mpcomplex& b)
{
        mpcomplex c(a);
	return (mpc_cmp(c.mpc, b.mpc)==0);
}

inline bool operator != (const mpcomplex& a, const mpcomplex& b)
{
	return (!mpc_cmp(a.mpc, b.mpc)==0);
}

inline bool operator != (const mpcomplex& a, const mpreal& b)
{
        mpcomplex c(b);
	return (!mpc_cmp(a.mpc, c.mpc)==0);
}

inline bool operator != (const mpreal& a, const mpcomplex& b)
{
        mpcomplex c(a);
	return (!mpc_cmp(b.mpc, c.mpc)==0);
}

inline mpcomplex::operator std::complex<double>() const
{
  mpreal re, im;
  std::complex<double> tmp;
  re = (*this).real();
  im = (*this).imag();
  tmp.real(mpfr_get_d(re,MPC_RND_RE(default_rnd)));
  tmp.imag(mpfr_get_d(im,MPC_RND_IM(default_rnd)));
  return tmp;
}

inline const mpcomplex operator/(const mpcomplex& a, const mpcomplex& b)
{
  mpcomplex tmp(a);
  tmp.set_prec2(std::max(a.get_prec_re(), b.get_prec_re()), std::max(a.get_prec_im(), b.get_prec_im()), mpcomplex::default_rnd);
  return tmp /= b;
}

inline const mpcomplex operator/(const mpreal& a, const mpcomplex& b)
{
  mpcomplex tmp(a);
  tmp.set_prec2(std::max(a.get_prec(), b.get_prec_re()), std::max(a.get_prec(), b.get_prec_im()), mpcomplex::default_rnd);
  return tmp /= b;
}

inline const mpcomplex operator/(const mpcomplex& a, const mpreal& b)
{
  mpcomplex tmp(a);
  tmp.set_prec2(std::max(a.get_prec_re(), b.get_prec()), std::max(a.get_prec_im(), b.get_prec()), mpcomplex::default_rnd);
  return tmp /= b;
}

inline const mpreal abs(const mpcomplex& a, mpfr_rnd_t rnd_mode = mpreal::default_rnd)
{
        mpreal x;
        mpcomplex y(a);
	mpc_abs((mpfr_ptr)(&x), y.mpc, rnd_mode);
	return x;
}

inline const mpreal norm(const mpcomplex& a, mpfr_rnd_t rnd_mode = mpreal::default_rnd)
{
        mpreal x;
        mpcomplex y(a);
	mpc_norm((mpfr_ptr)(&x), y.mpc, rnd_mode);
	return x;
}

inline const mpcomplex conj(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
        mpcomplex y(a);
	mpc_conj(y.mpc, y.mpc, rnd_mode);
	return y;
}

inline const mpreal arg(const mpcomplex& a, mpfr_rnd_t rnd_mode = mpreal::default_rnd)
{
        mpreal x(a.real());
        mpcomplex y(a);
	mpc_arg((mpfr_ptr)(&x), y.mpc, rnd_mode);
	return x;
}

inline const mpcomplex proj(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
        mpcomplex y(a);
	mpc_proj(y.mpc, y.mpc, rnd_mode);
	return y;
}

inline const mpcomplex sqrt(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_sqrt(x.mpc, x.mpc, rnd_mode);
	return x;
}

inline const mpcomplex pow(const mpcomplex& a, const mpcomplex& b, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_pow(x.mpc, x.mpc, b.mpc, rnd_mode);
	return x;
}

inline const mpcomplex exp(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_exp(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex log(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_log(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex sin(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_sin(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex cos(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_cos(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex tan(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_tan(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex sinh(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_sinh(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex cosh(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_cosh(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex tanh(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_tanh(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex asin(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_asin(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex acos(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_acos(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex atan(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_atan(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex asinh(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_asinh(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex acosh(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_acosh(x.mpc, a.mpc, rnd_mode);
	return x;
}

inline const mpcomplex atanh(const mpcomplex& a, mpc_rnd_t rnd_mode = mpcomplex::default_rnd)
{
	mpcomplex x(a);
	mpc_atanh(x.mpc, a.mpc, rnd_mode);
	return x;
}
inline const mpcomplex urandom_c (gmp_randstate_t& state)
{
	mpcomplex x;
	mpc_urandom(x.mpc,state);
	return x;
}

#if defined ___MPACK_BUILD_WITH_GMP___
inline mpcomplex::mpcomplex(const mpc_class& a)
{
   mp_prec_t pr, pi;
   pr = mpf_get_prec(a.real().get_mpf_t());  pi = mpf_get_prec(a.imag().get_mpf_t());
   mpc_init3(mpc, pr, pi);
   mpc_set_f_f(mpc, a.real().get_mpf_t(), a.imag().get_mpf_t(), default_rnd);
}

inline const mpcomplex operator-(const mpcomplex& a, const mpc_class& b)
{
   return mpcomplex(b) -= a;
}

inline const mpcomplex operator-(const mpc_class& a, const mpcomplex& b)
{
   return -(mpcomplex(a) -= b);
}

inline mpc_class cast2mpc_class(const mpcomplex &b)
{
//mpcomplex -> mpfr, mpfr -> mpf, mpf -> mpc_class
//I have to rewrite soon.
  mpf_t  tmpre, tmpim;
  mpfr_t tmpre2, tmpim2;
  mpcomplex a(b);
  mpreal are, aim;
  mp_prec_t pr, pi;

  pr = a.get_prec_re();
  pi = a.get_prec_im();
  are=a.real();
  aim=a.imag();

  mpfr_init2(tmpre2, pr); mpfr_init2(tmpim2, pi);
  mpfr_set(tmpre2, (mpfr_ptr)are, MPC_RND_RE(mpcomplex::default_rnd));
  mpfr_set(tmpim2, (mpfr_ptr)aim, MPC_RND_IM(mpcomplex::default_rnd));

  mpf_init2 (tmpre, pr); mpf_init2 (tmpim, pr);
  mpfr_get_f(tmpre, tmpre2, MPC_RND_RE(mpcomplex::default_rnd));  mpfr_get_f(tmpim, tmpim2, MPC_RND_IM(mpcomplex::default_rnd));

  mpc_class tmp(tmpre, tmpim);

  mpf_clear(tmpre);
  mpf_clear(tmpim);

  mpfr_clear(tmpre2);
  mpfr_clear(tmpim2);
  return tmp;
}     
inline mpcomplex& mpcomplex::operator=(const mpc_class& a)
{
   mpcomplex tmp(a);
   *this = tmp; 
   return *this;
}

#endif

#if defined ___MPACK_BUILD_WITH_QD___
inline mpcomplex::mpcomplex(const qd_complex& a)
{
   mpfr_t mp_real, mp_imag;
   mpc_init3(mpc, default_real_prec, default_imag_prec);

   mpfr_init2(mp_real, default_real_prec);
   mpfr_set_d (mp_real, a.real().x[0],  MPC_RND_RE(default_rnd));
   mpfr_add_d (mp_real, mp_real, a.real().x[1], MPC_RND_RE(default_rnd));
   mpfr_add_d (mp_real, mp_real, a.real().x[2], MPC_RND_RE(default_rnd));
   mpfr_add_d (mp_real, mp_real, a.real().x[3], MPC_RND_RE(default_rnd));

   mpfr_init2(mp_imag, default_imag_prec);
   mpfr_set_d (mp_imag, a.imag().x[0], MPC_RND_IM(default_rnd));
   mpfr_add_d (mp_imag, mp_imag, a.imag().x[1], MPC_RND_IM(default_rnd));
   mpfr_add_d (mp_imag, mp_imag, a.imag().x[2], MPC_RND_IM(default_rnd));
   mpfr_add_d (mp_imag, mp_imag, a.imag().x[3], MPC_RND_IM(default_rnd));
   
   mpc_set_fr_fr(mpc, mp_real, mp_imag, default_rnd);
   mpfr_clear(mp_imag);
   mpfr_clear(mp_real);
}
inline qd_complex cast2qd_complex(const mpcomplex &b)
{
  std::complex<double> p[4];
  mpcomplex c (b);
  qd_complex q;

  p[0] = c;
  c = c - p[0];
  p[1] = c;
  c = c - p[1];
  p[2] = c;
  c = c - p[2];
  p[3] = c;
//do not optimize
  q.real().x[0] = p[0].real();
  q.imag().x[0] = p[0].imag();
  q.real().x[1] = p[1].real();
  q.imag().x[1] = p[1].imag();
  q.real().x[2] = p[2].real();
  q.imag().x[2] = p[2].imag();
  q.real().x[3] = p[3].real();
  q.imag().x[3] = p[3].imag();
  return q;
}     
inline const mpcomplex operator-(const mpcomplex& a, const qd_complex& b)
{
   return mpcomplex(b) -= a;
}

inline const mpcomplex operator-(const qd_complex& a, const mpcomplex& b)
{
   return -(mpcomplex(a) -= b);
}

inline mpcomplex& mpcomplex::operator=(const qd_complex& a)
{
   mpcomplex tmp(a);
   *this = tmp; 
   return *this;
}

#endif

#if defined ___MPACK_BUILD_WITH_DD___
inline mpcomplex::mpcomplex(const dd_complex& a)
{
   mpfr_t mp_real, mp_imag;
   mpc_init3(mpc, default_real_prec, default_imag_prec);

   mpfr_init2(mp_real, default_real_prec);
   mpfr_set_d (mp_real, a.real().x[0],  MPC_RND_RE(default_rnd));
   mpfr_add_d (mp_real, mp_real, a.real().x[1], MPC_RND_RE(default_rnd));

   mpfr_init2(mp_imag, default_imag_prec);
   mpfr_set_d (mp_imag, a.imag().x[0], MPC_RND_IM(default_rnd));
   mpfr_add_d (mp_imag, mp_imag, a.imag().x[1], MPC_RND_IM(default_rnd));
   
   mpc_set_fr_fr(mpc, mp_real, mp_imag, default_rnd);
   mpfr_clear(mp_imag);
   mpfr_clear(mp_real);
}
inline dd_complex cast2dd_complex(const mpcomplex &b)
{
  std::complex<double> p[2];
  mpcomplex c (b);
  dd_complex q;

  p[0] = c;
  c = c - p[0];
  p[1] = c;
  c = c - p[1];
  q.real().x[0] = p[0].real();
  q.imag().x[0] = p[0].imag();
  q.real().x[1] = p[1].real();
  q.imag().x[1] = p[1].imag();

  return q;
}     
inline const mpcomplex operator-(const mpcomplex& a, const dd_complex& b)
{
   return mpcomplex(b) -= a;
}

inline const mpcomplex operator-(const dd_complex& a, const mpcomplex& b)
{
   return -(mpcomplex(a) -= b);
}

inline mpcomplex& mpcomplex::operator=(const dd_complex& a)
{
   mpcomplex tmp(a);
   *this = tmp; 
   return *this;
}

#endif

#if defined ___MPACK_BUILD_WITH___FLOAT128___
inline mpcomplex::mpcomplex(const std::complex<__float128>& a)
{
   mpfr_t mp_real, mp_imag;
   mpc_init3(mpc, default_real_prec, default_imag_prec);

   mpfr_init2(mp_real, default_real_prec);
   mpfr_set_float128 (&mp_real, a.real());

   mpfr_init2(mp_imag, default_imag_prec);
   mpfr_set_float128 (&mp_imag, a.imag());
   
   mpc_set_fr_fr(mpc, mp_real, mp_imag, default_rnd);
   mpfr_clear(mp_imag);
   mpfr_clear(mp_real);
}
inline std::complex<__float128> cast2complex__float128(const mpcomplex &b)
{
  std::complex<__float128> q;
  mpreal re_tmp, im_tmp; 
  re_tmp = b.real();
  im_tmp = b.imag();
  q.real() = mpfr_get_float128((mpfr_ptr)(re_tmp));
  q.imag() = mpfr_get_float128((mpfr_ptr)(im_tmp));
  return q;
}     
inline const mpcomplex operator-(const mpcomplex& a, const std::complex<__float128>& b)
{
   return mpcomplex(b) -= a;
}

inline const mpcomplex operator-(const std::complex<__float128>& a, const mpcomplex& b)
{
   return -(mpcomplex(a) -= b);
}

inline mpcomplex& mpcomplex::operator=(const std::complex<__float128>& a)
{
   mpcomplex tmp(a);
   *this = tmp; 
   return *this;
}

#endif

}

#endif /* __MP_COMPLEX_H__ */
