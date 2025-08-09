#pragma once

#include <initializer_list>
#include <array>
#include <iostream>

template<int N>
struct Vector {
	Vector() {
		this->data.fill(0.0f);
	}
	Vector(float v) {
		this->data.fill(v);
	}
	Vector(std::initializer_list<float> init) {
		size_t n = std::min(init.size(), data.size());
		std::copy_n(init.begin(), n, data.begin());
	}
	float& operator[](std::size_t index) {
		return data[index];
	}
	const float& operator[](std::size_t index) const {
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
	Vector operator*(float scalar) const {
		Vector result;
		for (int i = 0; i < N; ++i) {
			result.data[i] = this->data[i] * scalar;
		}
		return result;
	}
	Vector operator/(float scalar) const {
		Vector result;
		for (int i = 0; i < N; ++i) {
			result.data[i] = this->data[i] / scalar;
		}
		return result;
	}
	static Vector lerp(const Vector& v1, const Vector& v2, float t) {
		return v1 + (v2 - v1) * t;
	}
	static float dot(const Vector& v1, const Vector& v2) {
		float result = 0.0f;
		for (int i = 0; i < N; ++i) {
			result += v1[i] * v2[i];
		}
		return result;
	}
	float length() const {
		return sqrtf(dot(*this, *this));
	}
	
	friend std::ostream& operator<<(std::ostream& os, const Vector<N>& vec) {
		for (int i = 0; i < N; ++i) {
			os << vec.data[i];
			if (i < N - 1) os << ", ";
		}
		return os;
	}
	std::array<float, N> data;
};