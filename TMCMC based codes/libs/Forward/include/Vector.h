#ifndef VECTOR_H
#define VECTOR_H


template <typename ValueType1, typename ValueType2>
struct LargerType
{
	static const bool value = (sizeof(ValueType1) >= sizeof(ValueType2));

	typedef typename std::conditional<value,ValueType1,ValueType2>::type type;
};

template <typename ValueType>
struct CVector{

	ValueType x;
	ValueType y;
	ValueType z;

	CVector(const ValueType &xe=ValueType(), const ValueType &ye=ValueType(), const ValueType &ze=ValueType())
		: x(xe), y(ye), z(ze) {}

	template <typename ValueType2>
	CVector(const CVector<ValueType2> &rhs) : x((ValueType)rhs.x), y((ValueType)rhs.y), z((ValueType)rhs.z) {}

	template <typename ValueType2>
	inline CVector& operator=(const CVector<ValueType2> &rhs)
	{
		x = (ValueType)rhs.x;
		y = (ValueType)rhs.y;
		z = (ValueType)rhs.z;
		return (*this);
	}

	template <typename ValueType2>
	inline CVector& operator+=(const CVector<ValueType2> &rhs)
	{
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return (*this);
	}

	template <typename ValueType2>
	inline CVector& operator+=(const ValueType2 &rhs)
	{
		x += rhs;
		y += rhs;
		z += rhs;
		return (*this);
	}

	template <typename ValueType2>
	inline CVector& operator-=(const CVector<ValueType2> &rhs)
	{
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		return (*this);
	}

	template <typename ValueType2>
	inline CVector& operator-=(const ValueType2 &rhs)
	{
		x -= rhs;
		y -= rhs;
		z -= rhs;
		return (*this);
	}

	template <typename ValueType2>
	inline CVector& operator*=(const ValueType2 &rhs)
	{
		x *= rhs;
		y *= rhs;
		z *= rhs;
		return (*this);
	}

	template <typename ValueType2>
	inline CVector& operator/=(const ValueType2 &rhs)
	{
		x /= rhs;
		y /= rhs;
		z /= rhs;
		return (*this);
	}

	inline ValueType length()
	{
		ValueType tmp = std::sqrt( x*x + y*y + z*z );
		return tmp;
	}

};

// define two operator "+"
template <typename ValueType1, typename ValueType2>
inline CVector<typename LargerType<ValueType1,ValueType2>::type> operator+(const CVector<ValueType1> &first, const CVector<ValueType2> &second)
{
	CVector<typename LargerType<ValueType1,ValueType2>::type> tmp = first;
	tmp += second;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CVector<typename LargerType<ValueType1,ValueType2>::type> operator+(const CVector<ValueType1> &first, const ValueType2 &second)
{
	CVector<typename LargerType<ValueType1,ValueType2>::type> tmp = first;
	tmp += second;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CVector<typename LargerType<ValueType1,ValueType2>::type> operator+(const ValueType1 &first, const CVector<ValueType2> &second)
{
	CVector<typename LargerType<ValueType1,ValueType2>::type> tmp = second;
	tmp += first;
	return tmp;
}

// define two operator "-"
template <typename ValueType1, typename ValueType2>
inline CVector<typename LargerType<ValueType1,ValueType2>::type> operator-(const CVector<ValueType1> &first, const CVector<ValueType2> &second)
{
	CVector<typename LargerType<ValueType1,ValueType2>::type> tmp = first;
	tmp -= second;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CVector<typename LargerType<ValueType1,ValueType2>::type> operator-(const CVector<ValueType1> &first, const ValueType2 &second)
{
	CVector<typename LargerType<ValueType1,ValueType2>::type> tmp = first;
	tmp -= second;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CVector<typename LargerType<ValueType1,ValueType2>::type> operator-(const ValueType1 &first, const CVector<ValueType2> &second)
{
	CVector<typename LargerType<ValueType1,ValueType2>::type> tmp = second;
	tmp -= first;
	return tmp;
}

// define two operator "*"
template <typename ValueType1, typename ValueType2>
inline CVector<typename LargerType<ValueType1,ValueType2>::type> operator*(const CVector<ValueType1> &first, const ValueType2 &second)
{
	CVector<typename LargerType<ValueType1,ValueType2>::type> tmp = first;
	tmp *= second;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CVector<typename LargerType<ValueType1,ValueType2>::type> operator*(const ValueType1 &first, const CVector<ValueType2> &second)
{
	CVector<typename LargerType<ValueType1,ValueType2>::type> tmp = second;
	tmp *= first;
	return tmp;
}

// define two operator "/"
template <typename ValueType1, typename ValueType2>
inline CVector<typename LargerType<ValueType1,ValueType2>::type> operator/(const CVector<ValueType1> &first, const ValueType2 &second)
{
	CVector<typename LargerType<ValueType1,ValueType2>::type> tmp = first;
	tmp /= second;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CVector<typename LargerType<ValueType1,ValueType2>::type> operator/(const ValueType1 &first, const CVector<ValueType2> &second)
{
	CVector<typename LargerType<ValueType1,ValueType2>::type> tmp = second;
	tmp /= first;
	return tmp;
}

template <typename ValueType>
class CDyadic
{
public:

	ValueType xx, xy, xz;
	ValueType yx, yy, yz;
	ValueType zx, zy, zz;

	CDyadic(const ValueType &xxe=ValueType(), const ValueType &xye=ValueType(), const ValueType &xze=ValueType(),
			const ValueType &yxe=ValueType(), const ValueType &yye=ValueType(), const ValueType &yze=ValueType(),
			const ValueType &zxe=ValueType(), const ValueType &zye=ValueType(), const ValueType &zze=ValueType())
		: xx(xxe), xy(xye), xz(xze), yx(yxe), yy(yye), yz(yze), zx(zxe), zy(zye), zz(zze)  {}

	template <typename ValueType2>
	CDyadic(const CDyadic<ValueType2> &rhs) : xx((ValueType)rhs.xx), xy((ValueType)rhs.xy), xz((ValueType)rhs.xz), 
											  yx((ValueType)rhs.yx), yy((ValueType)rhs.yy), yz((ValueType)rhs.yz), 
											  zx((ValueType)rhs.zx), zy((ValueType)rhs.zy), zz((ValueType)rhs.zz) {}

	template <typename ValueType2>
	inline CDyadic& operator=(const CDyadic<ValueType2> &rhs)
	{
		xx = (ValueType)rhs.xx;		xy = (ValueType)rhs.xy;		xz = (ValueType)rhs.xz;
		yx = (ValueType)rhs.yx;		yy = (ValueType)rhs.yy;		yz = (ValueType)rhs.yz;
		zx = (ValueType)rhs.zx;		zy = (ValueType)rhs.zy;		zz = (ValueType)rhs.zz;
		return (*this);
	}

	template <typename ValueType2>
	inline CDyadic& operator+=(const CDyadic<ValueType2> &rhs)
	{
		xx += rhs.xx;		xy += rhs.xy;		xz += rhs.xz;
		yx += rhs.yx;		yy += rhs.yy;		yz += rhs.yz;
		zx += rhs.zx;		zy += rhs.zy;		zz += rhs.zz;
		return (*this);
	}

	template <typename ValueType2>
	inline CDyadic& operator+=(const ValueType2 &rhs)
	{
		xx += rhs;		xy += rhs;		xz += rhs;
		yx += rhs;		yy += rhs;		yz += rhs;
		zx += rhs;		zy += rhs;		zz += rhs;
		return (*this);
	}

	template <typename ValueType2>
	inline CDyadic& operator-=(const CDyadic<ValueType2> &rhs)
	{
		xx -= rhs.xx;		xy -= rhs.xy;		xz -= rhs.xz;
		yx -= rhs.yx;		yy -= rhs.yy;		yz -= rhs.yz;
		zx -= rhs.zx;		zy -= rhs.zy;		zz -= rhs.zz;
		return (*this);
	}

	template <typename ValueType2>
	inline CDyadic& operator-=(const ValueType2 &rhs)
	{
		xx -= rhs;		xy -= rhs;		xz -= rhs;
		yx -= rhs;		yy -= rhs;		yz -= rhs;
		zx -= rhs;		zy -= rhs;		zz -= rhs;
		return (*this);
	}

	template <typename ValueType2>
	inline CDyadic& operator*=(const ValueType2 &rhs)
	{
		xx *= rhs;		xy *= rhs;		xz *= rhs;
		yx *= rhs;		yy *= rhs;		yz *= rhs;
		zx *= rhs;		zy *= rhs;		zz *= rhs;
		return (*this);
	}

	template <typename ValueType2>
	inline CDyadic& operator/=(const ValueType2 &rhs)
	{
		xx /= rhs;		xy /= rhs;		xz /= rhs;
		yx /= rhs;		yy /= rhs;		yz /= rhs;
		zx /= rhs;		zy /= rhs;		zz /= rhs;
		return (*this);
	}

};

// define two operator "+"
template <typename ValueType1, typename ValueType2>
inline CDyadic<typename LargerType<ValueType1,ValueType2>::type> operator+(const CDyadic<ValueType1> &first, const CDyadic<ValueType2> &second)
{
	CDyadic<typename LargerType<ValueType1,ValueType2>::type> tmp = first;
	tmp += second;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CDyadic<typename LargerType<ValueType1,ValueType2>::type> operator+(const CDyadic<ValueType1> &first, const ValueType2 &second)
{
	CDyadic<typename LargerType<ValueType1,ValueType2>::type> tmp = first;
	tmp += second;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CDyadic<typename LargerType<ValueType1,ValueType2>::type> operator+(const ValueType1 &first, const CDyadic<ValueType2> &second)
{
	CDyadic<typename LargerType<ValueType1,ValueType2>::type> tmp = second;
	tmp += first;
	return tmp;
}

// define two operator "-"
template <typename ValueType1, typename ValueType2>
inline CDyadic<typename LargerType<ValueType1,ValueType2>::type> operator-(const CDyadic<ValueType1> &first, const CDyadic<ValueType2> &second)
{
	CDyadic<typename LargerType<ValueType1,ValueType2>::type> tmp = first;
	tmp -= second;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CDyadic<typename LargerType<ValueType1,ValueType2>::type> operator-(const CDyadic<ValueType1> &first, const ValueType2 &second)
{
	CDyadic<typename LargerType<ValueType1,ValueType2>::type> tmp = first;
	tmp -= second;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CDyadic<typename LargerType<ValueType1,ValueType2>::type> operator-(const ValueType1 &first, const CDyadic<ValueType2> &second)
{
	CDyadic<typename LargerType<ValueType1,ValueType2>::type> tmp = second;
	tmp -= first;
	return tmp;
}

// define two operator "*"
template <typename ValueType1, typename ValueType2>
inline CDyadic<typename LargerType<ValueType1,ValueType2>::type> operator*(const CDyadic<ValueType1> &first, const ValueType2 &second)
{
	CDyadic<typename LargerType<ValueType1,ValueType2>::type> tmp = first;
	tmp *= second;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CDyadic<typename LargerType<ValueType1,ValueType2>::type> operator*(const ValueType1 &first, const CDyadic<ValueType2> &second)
{
	CDyadic<typename LargerType<ValueType1,ValueType2>::type> tmp = second;
	tmp *= first;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CDyadic<typename LargerType<ValueType1,ValueType2>::type> operator*(const CDyadic<ValueType1> &first, const CDyadic<ValueType2> &second)
{
	CDyadic<typename LargerType<ValueType1,ValueType2>::type> tmp;
	tmp.xx = first.xx * second.xx + first.xy * second.yx + first.xz * second.zx;
	tmp.xy = first.xx * second.xy + first.xy * second.yy + first.xz * second.zy;
	tmp.xz = first.xx * second.xz + first.xy * second.yz + first.xz * second.zz;
	tmp.yx = first.yx * second.xx + first.yy * second.yx + first.yz * second.zx;
	tmp.yy = first.yx * second.xy + first.yy * second.yy + first.yz * second.zy;
	tmp.yz = first.yx * second.xz + first.yy * second.yz + first.yz * second.zz;
	tmp.zx = first.zx * second.xx + first.zy * second.yx + first.zz * second.zx;
	tmp.zy = first.zx * second.xy + first.zy * second.yy + first.zz * second.zy;
	tmp.zz = first.zx * second.xz + first.zy * second.yz + first.zz * second.zz;
	return tmp;
}

// define dyadic-vector multiplication
template <typename ValueType1, typename ValueType2>
inline CVector<typename LargerType<ValueType1,ValueType2>::type> operator*( const CDyadic<ValueType1> &first,  const CVector<ValueType2> &second)
{
	CVector<typename LargerType<ValueType1,ValueType2>::type> tmp;
	tmp.x = first.xx * second.x + first.xy * second.y + first.xz * second.z;
	tmp.y = first.yx * second.x + first.yy * second.y + first.yz * second.z;
	tmp.z = first.zx * second.x + first.zy * second.y + first.zz * second.z;
	return tmp;
}

// define two operator "/"
template <typename ValueType1, typename ValueType2>
inline CDyadic<typename LargerType<ValueType1,ValueType2>::type> operator/(const CDyadic<ValueType1> &first, const ValueType2 &second)
{
	CDyadic<typename LargerType<ValueType1,ValueType2>::type> tmp = first;
	tmp /= second;
	return tmp;
}

template <typename ValueType1, typename ValueType2>
inline CDyadic<typename LargerType<ValueType1,ValueType2>::type> operator/(const ValueType1 &first, const CDyadic<ValueType2> &second)
{
	CDyadic<typename LargerType<ValueType1,ValueType2>::type> tmp = second;
	tmp /= first;
	return tmp;
}

#endif