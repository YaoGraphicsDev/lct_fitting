#pragma once

#include <initializer_list>
#include <array>
#include <iostream>

template<typename T, int N>
struct Vector {
	typedef T EleType;

	Vector() {
		this->data.fill(0.0f);
	}
	Vector(EleType v) {
		this->data.fill(v);
	}
	Vector(std::initializer_list<EleType> init) {
		size_t n = std::min(init.size(), data.size());
		std::copy_n(init.begin(), n, data.begin());
	}
	EleType& operator[](std::size_t index) {
		return data[index];
	}
	const EleType& operator[](std::size_t index) const {
		return data[index];
	}
	Vector operator+(const Vector& other) const {
		Vector result;
		for (int i = 0; i < N; ++i) {
			result.data[i] = this->data[i] + other.data[i];
		}
		return result;
	}
	Vector operator-(const Vector& other) const {
		Vector result;
		for (int i = 0; i < N; ++i) {
			result.data[i] = this->data[i] - other.data[i];
		}
		return result;
	}
	Vector operator*(EleType scalar) const {
		Vector result;
		for (int i = 0; i < N; ++i) {
			result.data[i] = this->data[i] * scalar;
		}
		return result;
	}
	Vector operator/(EleType scalar) const {
		Vector result;
		for (int i = 0; i < N; ++i) {
			result.data[i] = this->data[i] / scalar;
		}
		return result;
	}
	static Vector lerp(const Vector& v1, const Vector& v2, EleType t) {
		return v1 + (v2 - v1) * t;
	}
	static EleType dot(const Vector& v1, const Vector& v2) {
		EleType result = 0.0f;
		for (int i = 0; i < N; ++i) {
			result += v1[i] * v2[i];
		}
		return result;
	}
	EleType length() const {
		return sqrtf(dot(*this, *this));
	}

	friend std::istream& operator>>(std::istream& is, Vector<T, N>& vec) {
		for (size_t i = 0; i < N; ++i) {
			is >> vec.data[i];
		}
		return is;
	}

	friend std::ostream& operator<<(std::ostream& os, const Vector<T, N>& vec) {
		for (int i = 0; i < N; ++i) {
			os << vec.data[i];
			if (i < N - 1) os << ", ";
		}
		return os;
	}
	std::array<EleType, N> data;
};