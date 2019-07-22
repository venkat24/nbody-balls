#pragma once

#include <cmath>
#include <iostream>
#include <utility>

namespace physics {

/**
 * Class for 2D vectors
 */
template <typename T> class Vector {
  public:
	Vector();
	Vector(T x, T y);

	/**
	 * Alias value for a null vector
	 */
	static Vector<T> null;

	/**
	 * Equal to operator for vectors
	 *
	 * @param[in]  rhs   The vector to be compared against
	 *
	 * @return     true if the vectors are equal, else false
	 */
	bool operator==(const Vector<T> &rhs) const;

	/**
	 * Not equal to operator for vectors
	 *
	 * @param[in]  rhs   The vector to be compared against
	 *
	 * @return     true if the vectors are not equal, else false
	 */
	bool operator!=(const Vector<T> &rhs) const;

	/**
	 * Addition operator for vectors
	 *
	 * @param[in]  rhs   The vector to be added with
	 *
	 * @return     Sum of the vectors
	 */
	Vector<T> operator+(const Vector<T> &rhs) const;

	/**
	 * Minus operator for vectors
	 * Subtracts the paramater passed from the vector
	 *
	 * @param[in]  rhs   The vector that subtracts from
	 *
	 * @return     Difference of the vectors
	 */
	Vector<T> operator-(const Vector<T> &rhs) const;

	/**
	 * Scalar Addition operator
	 * vec.x + scalar, vec.y + scalar
	 *
	 * @param[in]  scalar  The scalar to add
	 *
	 * @return     The new vector
	 */
	template <typename T2> Vector<T> operator+(const T2 &scalar) const;

	/**
	 * Scalar Subtraction operator
	 * vec.x - scalar, vec.y - scalar
	 *
	 * @param[in]  scalar  The scalar to subtract
	 *
	 * @return     The new vector
	 */
	template <typename T2> Vector<T> operator-(const T2 &scalar) const;

	/**
	 * Scalar multiplication operator
	 *
	 * @param[in]  scalar  The factor to multiply by
	 *
	 * @return     The scaled vector
	 */
	template <typename T2> Vector<T> operator*(const T2 &scalar) const;

	/**
	 * Scalar division operator
	 *
	 * @param[in]  scalar  The factor to divide by
	 *
	 * @return     The scaled vector
	 */
	template <typename T2> Vector<T> operator/(const T2 &scalar) const;

	template <typename T2>
	friend std::ostream &operator<<(std::ostream &ostream,
	                                const Vector<T2> &vector);

	/**
	 * Dot product of vectors
	 *
	 * @param[in]  rhs   The vector with which the operation is performed
	 *
	 * @return     The result of the dot product
	 */
	double dot(const Vector<T> &rhs) const;

	/**
	 * The magnitude of the vector
	 *
	 * @return     The magnitude of the vector
	 */
	double magnitude() const;

	/**
	 * The polar direction of the vector, in radians
	 *
	 * @return     The direction of the vector
	 */
	double direction() const;

	/**
	 * Return Vector as polar coordinates as a pair of (r, θ)
	 *
	 * @return std::pair<double, double> r, theta
	 */
	std::pair<double, double> as_polar() const;

	/**
	 * Set x and y from a polar pair of (r, θ)
	 * 
	 * @param r radius
	 * @param theta angle θ
	 */
	void set_as_polar(double r, double theta);

	/**
	 * Distance between this and another vector
	 *
	 * @param[in]  other  The other vector
	 *
	 * @return     The distance between the two vectors
	 */
	double distance(const Vector<T> &other) const;

	/**
	 * Calculates floor of the members of the vector
	 *
	 * @return     The result Vector with floored values
	 */
	Vector<T> floor() const;

	/**
	 * Calculates ceiling of the members of the vector
	 *
	 * @return     The result Vector with ceiled values
	 */
	Vector<T> ceil() const;

	/**
	 * Return a new Vector<double_t> from another Vector type
	 *
	 * @param vec A numeric Vector type
	 * @return Vector<double_t>
	 */
	Vector<double_t> to_double() const;

	/**
	 * Return a new Vector<int64_t> from another Vector type. Floor values.
	 *
	 * @param vec A numeric Vector type
	 * @return Vector<int64_t>
	 */
	Vector<int64_t> to_int() const;

	/**
	 * Return true if this vector has any value other than null
	 *
	 * @return true if value is Vector<T>::null
	 * @return false otherwise
	 */
	explicit operator bool() const;

	T x;
	T y;
};

template <typename T> Vector<T>::Vector() : x(), y() {}

template <typename T> Vector<T>::Vector(T x, T y) : x(x), y(y) {}

template <typename T> bool Vector<T>::operator==(const Vector<T> &rhs) const {
	return (x == rhs.x && y == rhs.y);
}

template <typename T> bool Vector<T>::operator!=(const Vector<T> &rhs) const {
	return (x != rhs.x || y != rhs.y);
}

template <typename T>
Vector<T> Vector<T>::operator+(const Vector<T> &rhs) const {
	return Vector(x + rhs.x, y + rhs.y);
}

template <typename T>
Vector<T> Vector<T>::operator-(const Vector<T> &rhs) const {
	return Vector(x - rhs.x, y - rhs.y);
}

template <typename T>
template <typename T2>
Vector<T> Vector<T>::operator+(const T2 &scalar) const {
	return Vector(x + scalar, y + scalar);
}

template <typename T>
template <typename T2>
Vector<T> Vector<T>::operator-(const T2 &scalar) const {
	return Vector(x - scalar, y - scalar);
}

template <typename T>
template <typename T2>
Vector<T> Vector<T>::operator*(const T2 &scalar) const {
	return Vector(x * scalar, y * scalar);
}

template <typename T>
template <typename T2>
Vector<T> Vector<T>::operator/(const T2 &scalar) const {
	return Vector(x / scalar, y / scalar);
}

template <typename T>
std::ostream &operator<<(std::ostream &ostream, const Vector<T> &vector) {
	ostream << "(" << vector.x << ", " << vector.y << ")";
	return ostream;
}

template <typename T> double Vector<T>::dot(const Vector<T> &rhs) const {
	return (x * rhs.x + y * rhs.y);
}

template <typename T> double Vector<T>::magnitude() const {
	return sqrt(x * x + y * y);
}

template <typename T> double Vector<T>::direction() const {
	return atan2(y, x);
}

template <typename T>
void Vector<T>::set_as_polar(double r, double theta) {
	x = r * cos(theta);
	y = r * sin(theta);
}

template <typename T> std::pair<double, double> Vector<T>::as_polar() const {
	double r = magnitude();
	double theta = direction();
	return std::make_pair(r, theta);
}

template <typename T> double Vector<T>::distance(const Vector<T> &other) const {
	return sqrt(pow(x - other.x, 2) + pow(y - other.y, 2));
}

template <typename T> Vector<T> Vector<T>::floor() const {
	return Vector<T>(std::floor(x), std::floor(y));
}

template <typename T> Vector<T> Vector<T>::ceil() const {
	return Vector<T>(std::ceil(x), std::ceil(y));
}

template <typename T> Vector<double_t> Vector<T>::to_double() const {
	return Vector<double_t>(static_cast<double_t>(x), static_cast<double_t>(y));
}

template <typename T> Vector<int64_t> Vector<T>::to_int() const {
	return Vector<int64_t>(static_cast<int64_t>(x), static_cast<int64_t>(y));
}

template <typename T> Vector<T>::operator bool() const {
	if (*this == Vector<T>::null) {
		return false;
	}
	return true;
}

} // namespace physics

using Vector = physics::Vector<double>;
