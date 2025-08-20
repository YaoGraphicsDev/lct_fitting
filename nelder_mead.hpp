#pragma once

#include "vector_n.hpp"

#include <iostream>
#include <algorithm>
#include <array>
#include <list>
#include <functional>
#include <cassert>

template<int DIM>
struct NelderMead{

	using VecN = Vector<float, DIM>;
	static constexpr int SimplexSize = DIM + 1;
	struct Vertex {
		VecN x;
		double f;
	};
	using Simplex = std::array<Vertex, SimplexSize>;

	static void ascend_sort(typename Simplex::iterator begin, typename Simplex::iterator end) {
		std::sort(begin, end, [&](Vertex& v1, Vertex& v2) {
			return v1.f < v2.f;
		});
	}

	static VecN centroid(typename Simplex::iterator begin, typename Simplex::iterator end) {
		VecN c(0.0f);
		for (Simplex::iterator iter = begin; iter < end; ++iter) {
			c = c + iter->x;
		}
		c = c / float(end - begin);

		return c;
	}

	static void shrink(typename  Simplex::iterator begin, typename  Simplex::iterator end, const VecN& x_o, std::function<double(const VecN&)> obj_func) {
		for (Simplex::iterator iter = begin; iter < end; ++iter) {
			iter->x = VecN::lerp(iter->x, x_o, 0.5f);
			iter->f = obj_func(iter->x);
		}
	}

	// pass objective function as parameter
	struct Opti {
		VecN x;
		float error;
		int n_iter;
	};

	static Opti optimize(
		std::function<double(VecN)> obj_func,
		VecN start_x,
		float delta,
		float tolerance,
		int max_iter) {

		Simplex simplex;

		// initialize starting vertices
		simplex[0].x = start_x;
		for (int i = 1; i < simplex.size(); ++ i) {
			simplex[i].x = start_x;
			simplex[i].x[i - 1] += delta; // spread simplex vertices across dimensions
		}

		// one round of evaluation
		for (Vertex& vert : simplex) {
			vert.f = obj_func(vert.x);
		}

		float error = 0.0f;
		int iter = 0;
		VecN x_opt;
		for (iter = 0; iter < max_iter; ++iter) {
			if (iter == max_iter - 5) {
				int a = 0;
			}

			ascend_sort(simplex.begin(), simplex.end());

			Vertex& hi = simplex.back();
			Vertex& lo = simplex.front();
			assert(simplex.size() >= 3);
			Vertex& nh = *(simplex.end() - 2);

			// TODO: check tolerance, break
			double a = glm::abs(lo.f);
			double b = glm::abs(hi.f);
			error = fabs(a - b);
			// std::cout << "error = " << error << std::endl;
			if (error < tolerance) {
				x_opt = lo.x;
				return { lo.x, error, iter };
			}

			VecN x_o = centroid(simplex.begin(), simplex.end() - 1);
			
			// reflect
			VecN x_r = VecN::lerp(hi.x, x_o, 2.0f);
			double f_r = obj_func(x_r);
			if (f_r >= lo.f && f_r < nh.f) {
				hi.x = x_r;
				hi.f = f_r;
			}
			else if (f_r < lo.f) {
				// expansion
				VecN x_e = VecN::lerp(x_o, x_r, 2.0f);
				double f_e = obj_func(x_e);
				if (f_e < f_r) {
					hi.x = x_e;
					hi.f = f_e;
				}
				else {
					hi.x = x_r;
					hi.f = f_r;
				}
			}
			else {
				// contraction
				assert(f_r >= nh.f);
				if (f_r < hi.f) {
					VecN x_c = VecN::lerp(x_o, x_r, 0.5f);
					double f_c = obj_func(x_c);
					if (f_c < f_r) {
						hi.x = x_c;
						hi.f = f_c;
					}
					else {
						shrink(simplex.begin() + 1, simplex.end(), lo.x, obj_func);
					}
				}
				else {
					VecN x_c = VecN::lerp(x_o, hi.x, 0.5f);
					double f_c = obj_func(x_c);
					if (f_c < hi.f) {
						hi.x = x_c;
						hi.f = f_c;
					}
					else {
						shrink(simplex.begin() + 1, simplex.end(), lo.x, obj_func);
					}
				}
			}
		}
		
		// should not run out of iterations before convergence
		assert(false);
		std::cout << "Out of iterations. Error = " << error << ", tolerance = " << tolerance << std::endl;
		
		ascend_sort(simplex.begin(), simplex.end());
		x_opt = simplex.front().x;
		return {x_opt, error, iter};
	}
};

