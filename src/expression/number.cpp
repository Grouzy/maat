#include "maat/number.hpp"
#include "boost/numeric/conversion/cast.hpp"

namespace maat
{

using namespace maat::serial;

// Interpret src as a signed mpz value and puts it in res
// src must be more than 64 bits, mpz_t is initialized
void mpz_init_force_signed(boost::multiprecision::uint512_t& res, const Number& src)
{
    if (!src.is_mpz())
        throw expression_exception("mpz_force_signed(): shouldn't be called with regular Number!");
    
    int bit = boost::multiprecision::bit_test(src.mpz_, src.size-1);
    if (bit == 0)
        res.assign(src.mpz_); // Unsigned, keep the same value
    else
    {
        boost::multiprecision::uint512_t tmp;
        boost::multiprecision::bit_set(tmp, src.size);
        tmp -= src.mpz_;
        res = tmp;
        res.backend().negate();
    }
}

// TODO(boyan): is setting mpz to zero causing memory allocation here ???
Number::Number(): size(0), cst_(-1), mpz_(0){}

Number::Number(size_t bits): size(bits), cst_(0), mpz_(0){}

Number::Number(size_t bits, const std::string& value, int base): size(bits), cst_(0)
{
    set_mpz(value, base);
}

Number::~Number(){}

uid_t Number::class_uid() const {return ClassId::NUMBER;}

void Number::dump(Serializer& s) const
{
    bool _is_mpz = is_mpz();
    s << bits(size) << bits(_is_mpz);
    if (_is_mpz)
    {
        // If mpz, dump as string
        std::stringstream ss;
        print(ss, true); // Decimal = True
        s << ss.str();
    }
    else
        s << bits(cst_);
}

void Number::load(Deserializer& d)
{
    bool _is_mpz;
    d >> bits(size) >> bits(_is_mpz);
    if (_is_mpz)
    {
        std::string val;
        d >> val;
        set_mpz(val, 10); // Base 10 because dump() prints with decimal = True
    }
    else
    {
        cst_t val;
        d >> bits(val);
        set_cst(val);
    }
}

void Number::adjust_mpz()
{
    boost::multiprecision::uint512_t tmp;
    if (!is_mpz())
        return;

    tmp = mpz_;
    mpz_ = boost::multiprecision::uint512_t(0);

    // Copy bit by bit
    for (unsigned int i = 0; i < size; i++)
    {
        if (boost::multiprecision::bit_test(tmp, i))
        {
            boost::multiprecision::bit_set(mpz_, i);
        }
        else
        {
            boost::multiprecision::bit_unset(mpz_, i);
        }
    }
}

cst_t __number_cst_mask(size_t size)
{
    if( size == sizeof(cst_t)*8 )
        return (cst_t)-1;
    else
        return ((ucst_t)1<<(ucst_t)size)-1; 
}

ucst_t __number_cst_unsign_trunc(size_t size, cst_t c)
{
    if( size >= sizeof(cst_t)*8)
    {
        return c;
    }
    return (ucst_t)__number_cst_mask(size) & (ucst_t)c;
}

cst_t __number_cst_sign_extend(size_t size, cst_t val)
{
    if( size == sizeof(cst_t)*8)
    {
        return val;
    }
    if( ((ucst_t)1<<((ucst_t)size-1)) & (ucst_t)val )
    {
        // Negative, set higher bits to 1
        val = ((ucst_t)0xffffffffffffffff<< size) | val;
    }
    else
    {
        // Positive, set higher bits to 0
        val = ((((ucst_t)1<<size)-1) & val);
    }
    return val;
}

Number::Number(size_t bits, cst_t value): size(bits)
{
    if (bits > 64)
        set_mpz(value);
    else
        cst_ = __number_cst_sign_extend(size, value);
}

/// Set the number to simple value 'val'
void Number::set_cst(cst_t val)
{
    // Truncate/extend value if needed
    cst_ = __number_cst_sign_extend(size, val);
}

void Number::set(cst_t val)
{
    cst_ = val;
    if (is_mpz())
    {
        mpz_ = boost::multiprecision::uint512_t(val);
        adjust_mpz();
    }
}

cst_t Number::get_cst() const
{
    if (!is_mpz())
        return cst_;
    else
    {
        cst_t res = 0;
        for (int i = (sizeof(cst_t)*8) -1; i >= 0; i--)
        {
            res = (res<<1) + boost::multiprecision::bit_test(mpz_, i);
        }
        return res;
    }
}

ucst_t Number::get_ucst() const
{
    if (!is_mpz())
        return __number_cst_unsign_trunc(size, cst_);
    else
    {
        cst_t res = 0;
        for (int i = (sizeof(cst_t)*8) -1; i >= 0; i--)
        {
            res = (res<<1) + boost::multiprecision::bit_test(mpz_, i);
        }
        return __number_cst_unsign_trunc(size, res);
    }
}

/// Set the number to multiprecision value 'val'
void Number::set_mpz(cst_t val)
{
    mpz_ = boost::multiprecision::uint512_t(val);
    adjust_mpz();
}

void Number::set_mpz(std::string_view val, int base)
{
    if (base < 2 || base > 62)
        throw expression_exception("Number::set_mpz() needs a base between 2 and 62");

    std::string boost_parsable;
    switch(base)
    {
        using namespace std::string_literals;
        case 2:
            {
                if(val.starts_with("0b"))
                {
                    boost_parsable = val;
                    break;
                }

                boost_parsable = "0b" + std::string(val);
                break;
            }
        case 16:
            {
                if(val.starts_with("0x"))
                {
                    boost_parsable = val;
                    break;
                }

                boost_parsable = "0x" + std::string(val);
                break;
            }
        case 10:
            {
                boost_parsable = val;
                break;
            }
        default:
            {
                assert(false);
            }
    }

    mpz_ = boost::multiprecision::uint512_t(boost_parsable);
    adjust_mpz();
}

void Number::set_neg(const Number& n)
{
    size = n.size;
    if (n.size <= 64)
        set_cst(-1 * n.cst_);
    else
    {
        mpz_ = n.mpz_;
        mpz_.backend().negate();
        adjust_mpz();
    }
}

void Number::set_nan(const Number& n)
{
    size = n.size;
    if (n.size <= 64)
        set_cst(boost::math::isnan((fcst_t)n.cst_));
    else
    {
        // TODO use 128bit float
        mpz_ = boost::math::isnan(n.mpz_.convert_to<double>());
        adjust_mpz();
    }
}

void Number::set_not(const Number& n)
{
    size = n.size;
    if (n.size <= 64)
        set_cst(~(ucst_t)(n.cst_));
    else
    {
        mpz_ = ~n.mpz_;
        adjust_mpz();
    }
}

void Number::set_add(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
        set_cst(n1.cst_ + n2.cst_);
    else
    {
        mpz_ = n1.mpz_ + n2.mpz_;
        adjust_mpz();
    }
}

void Number::set_sub(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
        set_cst(n1.cst_ - n2.cst_);
    else
    {
        mpz_ = n1.mpz_ - n2.mpz_;
        adjust_mpz();
    }
}

void Number::set_and(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
        set_cst((ucst_t)n1.cst_ & (ucst_t)n2.cst_);
    else
    {
        mpz_ = n1.mpz_ & n2.mpz_;
    }
}

void Number::set_mul(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
    {
        set_cst((ucst_t)n1.cst_ * (ucst_t)n2.cst_);
    }
    else
    {
        mpz_ = n1.mpz_ * n2.mpz_;
        adjust_mpz();
    }
}

void Number::set_fmul(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
    {
        set_cst((fcst_t)n1.cst_ * (fcst_t)n2.cst_);
    }
    else
    {
        mpz_ = n1.mpz_ * n2.mpz_;
        adjust_mpz();
    }
}

void Number::set_xor(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
    {
        set_cst(n1.cst_ ^ n2.cst_);
    }
    else
    {
        mpz_ = n1.mpz_ ^ n2.mpz_;
        adjust_mpz();
    }
}

void Number::set_rem(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
    {
        set_cst(__number_cst_unsign_trunc(n1.size, n1.cst_) % __number_cst_unsign_trunc(n2.size, n2.cst_));
    }
    else
    {
        mpz_ = n1.mpz_ % n2.mpz_;
        adjust_mpz();
    }
}

void Number::set_srem(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
    {
        set_cst(__number_cst_sign_extend(n1.size, n1.cst_) % __number_cst_sign_extend(n2.size, n2.cst_));
    }
    else
    {
        boost::multiprecision::uint512_t tmp1, tmp2;
        mpz_init_force_signed(tmp1, n1);
        mpz_init_force_signed(tmp2, n2);
        mpz_ = tmp1 % tmp2;
        adjust_mpz();
    }
}

// Helper for exponentiation on integers
ucst_t uint_pow(ucst_t base, ucst_t exp)
{
    ucst_t result = 1;
    while (exp > 0)
    {
        if (exp % 2)
           result *= base;
        exp /= 2;
        base *= base;
    }
    return result;
}

void Number::set_exp(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
    {
        set_cst(uint_pow(n1.get_ucst(), n2.get_ucst()));
    }
    else
    {
        mpz_ = 1;
        for (boost::multiprecision::uint512_t i = 0; i < n2.mpz_; i++)
        {
            mpz_ *= n1.mpz_;
        }

        adjust_mpz();
    }
}

void Number::set_shl(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
    {
        cst_t tmp;
        if (n2.cst_ >= n1.size)
            tmp = 0;
        else
            tmp = ((ucst_t)n1.cst_) << ((ucst_t)n2.cst_);
        set_cst(tmp);
    }
    else
    {
        mpz_ = n1.mpz_ << n2.get_cst();
        adjust_mpz();
    }
}

void Number::set_shr(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
    {
        cst_t tmp;
        if (n2.cst_ >= n1.size)
            tmp = 0;
        else
            tmp = n1.get_ucst() >> n2.get_ucst();
        set_cst(tmp);
    }
    else
    {
        mpz_ = n1.mpz_ >> n2.get_cst();
        adjust_mpz();
    }
}

void Number::set_sar(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
    {
        cst_t tmp;
        if (n2.cst_ >= n1.size)
        {
            if( n1.cst_ & (0x1 << (n1.size-1)))
                tmp = 0xffffffffffffffff;
            else
                tmp = 0;
        }
        else
        {
            tmp = n1.get_cst() >> n2.get_ucst();
        }
        set_cst(tmp);
    }
    else
    {
        mpz_ = 0;
        uint64_t shift = n2.mpz_.convert_to<uint64_t>();
        uint64_t i;
        // Copy bits
        for (i = 0; i < size - shift; i++)
        {
            if (boost::multiprecision::bit_test(n1.mpz_, i + shift) == 1)
                boost::multiprecision::bit_set(mpz_, i);
            else
                boost::multiprecision::bit_unset(mpz_, i);
        }
        // Set the shifted mask to 0 or 0xffff....
        if (boost::multiprecision::bit_test(n1.mpz_, n1.size-1) == 1)
            for (i = 0; i < shift; i++)
                boost::multiprecision::bit_set(mpz_, size-1-i);
        else 
            for (i = 0; i < shift; i++)
                boost::multiprecision::bit_unset(mpz_, size-1-i);
        // Adjust
        adjust_mpz();
    }
}

void Number::set_or(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
    {
        set_cst(n1.cst_ | n2.cst_);
    }
    else
    {
        mpz_ = n1.mpz_ | n2.mpz_;
    }
}

void Number::set_sdiv(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
    {
        set_cst(
            __number_cst_sign_extend(n1.size, n1.get_cst()) /
            __number_cst_sign_extend(n2.size, n2.get_cst())
        );
    }
    else
    {
        
        boost::multiprecision::uint512_t tmp1, tmp2;
        mpz_init_force_signed(tmp1, n1);
        mpz_init_force_signed(tmp2, n2);
        mpz_ = tmp1 / tmp2;
        adjust_mpz();
    }
}

void Number::set_div(const Number& n1, const Number& n2)
{
    size = n1.size;
    if (size <= 64)
    {
        ucst_t t1 = n1.get_ucst();
        ucst_t t2 = n2.get_ucst();
        set_cst(t1 / t2);
    }
    else
    {
        // TODO: this is signed division, not unsigned ???
        mpz_ = n1.mpz_ / n2.mpz_;
        adjust_mpz();
    }
}

void Number::set_extract(const Number& n, unsigned int high, unsigned int low)
{
    cst_t tmp;
    size_t tmp_size = high - low + 1;
    if (n.size <= 64)
    {
        ucst_t mask;
        if( high == 63 )
            mask = 0xffffffffffffffff;
        else
            mask = (((cst_t)1 << (high+1))-1);

        tmp =  ((ucst_t)n.cst_ & mask) >> (ucst_t)low;
        size = tmp_size;
        set_cst(tmp);
    }
    else
    {
        boost::multiprecision::uint512_t tmp;
        tmp.assign(0);
        // Copy bit by bit
        for (unsigned int i = 0; i < tmp_size; i++)
        {
            if (boost::multiprecision::bit_test(n.mpz_, i+low) == 1)
                boost::multiprecision::bit_set(tmp, i);
            else
                boost::multiprecision::bit_unset(tmp, i);
        }

        size = tmp_size;
        mpz_ = boost::multiprecision::uint512_t(tmp);

        // adjust_mpz(); no need to adjust, we set bits manually
        // If result size on 64 bits or less, transform into cst, not mpz 
        if (this->size <= 64)
        {
            set_cst(mpz_.convert_to<uint64_t>());
        }
    }
}

void Number::set_concat(const Number& n1, const Number& n2)
{
    size_t tmp_size = n1.size + n2.size; // Use tmp size because *this might be n1 or n2

    if (tmp_size <= 64)
    {
        cst_t tmp = n2.cst_;
        // Mask higher bits before doing OR
        tmp &= (((ucst_t)1<<(ucst_t)n2.size)-1);
        // Do OR to set higher part
        tmp |= (ucst_t)n1.cst_ << (ucst_t)n2.size;
        size = tmp_size;
        set_cst(tmp);
    }
    else
    {
        boost::multiprecision::uint512_t tmp;
        // Set higher (set then shift)
        if (n1.is_mpz())
        {
            tmp = n1.mpz_;
        }
        else
        {
            tmp = n1.get_ucst();
        }

        tmp <<= n2.size;
        // Set lower
        if (n2.is_mpz())
        {
            tmp |= n2.mpz_;
        }
        else
        {
            tmp |= n2.get_ucst();
        }

        mpz_ = tmp;
        size = tmp_size;
        adjust_mpz();
    }
}

void Number::set_popcount(int dest_size, const Number& n)
{
    size = dest_size;
    ucst_t res = 0;
    if (n.size <= 64)
    {
        for (int i = 0; i < n.size; i++)
        {
            res += (n.cst_ >> i) & 1;
        }
    }
    else
    {
        for (int i = 0; i < n.size; i++)
        {
            res += boost::multiprecision::bit_test(n.mpz_, i);
        }
    }

    // Assign res
    if (size <= 64)
        set_cst(res);
    else
        set_mpz(res);
}

void Number::set_zext(int ext_size, const Number& n)
{
    this->size = ext_size;
    if (ext_size <= 64)
    {
        cst_t tmp =  ((ucst_t)__number_cst_unsign_trunc(n.size, n.cst_));
        set_cst(tmp);
    }
    else
    {
        if (n.is_mpz())
            mpz_ = n.mpz_;
        else
            mpz_ = n.get_ucst();
        // Extend higher bits to zero
        for (unsigned int i = n.size; i < ext_size; i++)
        {
            boost::multiprecision::bit_unset(mpz_, i);
        }
    }
}

void Number::set_sext(int ext_size, const Number& n)
{
    this->size = ext_size;
    if (ext_size <= 64)
    {
        cst_t tmp =  ((ucst_t)__number_cst_unsign_trunc(n.size, n.cst_));
        if (tmp & (1U << (n.size-1)))
        {
            // hsb is 1 add mask
            tmp |= (__number_cst_mask(ext_size - n.size) << n.size);
        }
        set_cst(tmp);
    }
    else
    {
        if (n.is_mpz())
            mpz_ = n.mpz_;
        else
            mpz_ = n.get_ucst();
        // Extend higher bits
        bool hsb_set = boost::multiprecision::bit_test(mpz_, n.size-1);
        for (unsigned int i = n.size; i < ext_size; i++)
        {
            if (hsb_set)
                boost::multiprecision::bit_set(mpz_, i);
            else
                boost::multiprecision::bit_unset(mpz_, i);
        }
        adjust_mpz();
    }
}

void Number::set_mask(int mask_size)
{
    if (size <= 64)
    {
        set_cst(__number_cst_mask(mask_size));
    }
    else
    {
        for (unsigned int i = 0; i < mask_size; i++)
        {
            boost::multiprecision::bit_set(mpz_, i);
        }
    }
}

void Number::set_overwrite(const Number& n1, const Number& n2, int lb)
{
    ucst_t mask;
    ucst_t res;

    if (n2.size + lb > n1.size)
        throw expression_exception("Number::set_overwrite(): src number is too big to fit in dest!");
    if (n2.size == n1.size)
    {
        *this = n2;
        return;
    }

    if (n1.size <= 64)
    {
        mask = ~ ((ucst_t)__number_cst_mask(n2.size) << (ucst_t)lb);
        res = (n1.cst_ & mask) | (__number_cst_unsign_trunc(n2.size, n2.cst_) << lb);

        // Set size at the end because n2 might be a reference to 'this'
        this->size = n1.size;
        set_cst(res);
    }
    else
    {
        // Make copies in case n1 or n2 is a reference to 'this'
        boost::multiprecision::uint512_t tmp = n1.mpz_;
        boost::multiprecision::uint512_t tmp2 = n2.is_mpz() ? n2.mpz_ : n2.get_ucst();
        for (int i = 0; i < n2.size; i++)
        {
            if (boost::multiprecision::bit_test(tmp2, i) == 1)
            {
                boost::multiprecision::bit_set(tmp, i + lb);
            }
            else
            {
                boost::multiprecision::bit_unset(tmp, i + lb);
            }
        }
        mpz_ = tmp;
        this->size = n1.size;
    }
}

bool Number::sless_than(const Number& other) const
{
    if (size <= 64)
    {
        return (cst_t)cst_ < (cst_t)other.cst_;
    }
    else
    {
        // mpz_cmp returns a positive value if op1 > op2, 
        // zero if op1 = op2, or a negative value if op1 < op2
        boost::multiprecision::uint512_t s1, s2;
        mpz_init_force_signed(s1, *this);
        mpz_init_force_signed(s2, other);
        return boost::multiprecision::isless(s1, s2);
    }
}

bool Number::slessequal_than(const Number& other) const
{
    if (size <= 64)
    {
        return (cst_t)cst_ <= (cst_t)other.cst_;
    }
    else
    {
        return sless_than(other) || equal_to(other);
    }
}

bool Number::less_than(const Number& other) const
{
    if (size <= 64)
    {
        return (ucst_t)cst_ < (ucst_t)other.cst_;
    }
    else
    {   
        return mpz_ < other.mpz_;
    }
}

bool Number::fless_than(const Number& other) const
{
    if (size <= 64)
    {
        return (fcst_t)cst_ < (fcst_t)other.cst_;
    }
    else
    {   
        return mpz_ < other.mpz_;
    }
}

bool Number::lessequal_than(const Number& other) const
{
    if (size <= 64)
    {
        return (ucst_t)cst_ <= (ucst_t)other.cst_;
    }
    else
    {   
        return less_than(other) || equal_to(other);
    }
}

bool Number::equal_to(const Number& other) const
{
    if (size <= 64)
    {
        return (cst_t)cst_ == (cst_t)other.cst_;
    }
    else
    {
        // mpz_cmp returns a positive value if op1 > op2, 
        // zero if op1 = op2, or a negative value if op1 < op2
        return mpz_.compare(other.mpz_) == 0;
    }
}

bool Number::fequal_to(const Number& other) const
{
    if (size <= 64)
    {
        return (fcst_t)cst_ == (fcst_t)other.cst_;
    }
    else
    {
        // mpz_cmp returns a positive value if op1 > op2, 
        // zero if op1 = op2, or a negative value if op1 < op2
        return mpz_.compare(other.mpz_) == 0;
    }
}

bool Number::is_null() const
{
    if (size <= 64)
        return cst_ == 0;
    else
        return mpz_.is_zero();
}

bool Number::is_mpz() const
{
    return size > 64;
}

std::ostream& operator<<(std::ostream& os, const Number& n)
{
    n.print(os);
    return os;
}

const char* __hex_format = "%Zx";
const char* __dec_format = "%Zd";

void Number::print(std::ostream& os, bool decimal) const
{
    if (is_mpz())
    {
        if (!decimal)
            os << std::hex << "0x";
        os << mpz_;
    }
    else
        os << std::showbase << get_ucst() << std::noshowbase;
}

} // namespace maat
