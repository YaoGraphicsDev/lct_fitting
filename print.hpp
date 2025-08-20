#include "brdf.h"
#include "glm/glm.hpp"
#include "lct.hpp"
#include "fstream"

void print_lct(const std::string& name, LCT& lct, int n = 128) {
	std::ofstream ofs;
	ofs.open(name, std::fstream::out);
	if (!ofs.is_open()) {
		std::cout << "Cannot open file " << name << std::endl;
		return;
	}

	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= n; ++j) {
			float theta = glm::half_pi<float>() * glm::clamp((float)i / n, 0.0f, 1.0f);
			float phi = glm::two_pi<float>() * glm::clamp((float)j / n, 0.0f, 1.0f);
			glm::vec3 l = UVA::spherical(theta, phi);
			LCT::EvalResult lct_eval = lct.eval(l);

			ofs << glm::clamp((float)i / n, 0.0f, 1.0f) * 0.5f << ", ";
			ofs << glm::clamp((float)j / n, 0.0f, 1.0f) * 2.0f << ", ";
			ofs << lct_eval.brdfxcos_l << std::endl;
		}
	}

	std::cout << "lct data written to " << name << std::endl;
}

void print_ggx(const std::string& name, BRDF& brdf, glm::vec3 v, int n = 256) {
	std::ofstream ofs;
	ofs.open(name, std::fstream::out);
	if (!ofs.is_open()) {
		std::cout << "Cannot open file " << name << std::endl;
		return;
	}

	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= n; ++j) {
			float theta = glm::half_pi<float>() * glm::clamp((float)i / n, 0.0f, 1.0f);
			float phi = glm::two_pi<float>() * glm::clamp((float)j / n, 0.0f, 1.0f);
			glm::vec3 l = UVA::spherical(theta, phi);
			BRDF::EvalResult brdf_eval = brdf.eval(v, l);

			ofs << glm::clamp((float)i / n, 0.0f, 1.0f) * 0.5f << ", ";
			ofs << glm::clamp((float)j / n, 0.0f, 1.0f) * 2.0f << ", ";
			ofs << brdf_eval.brdfxcos_l << std::endl;
		}
	}

	std::cout << "ggx brdf data written to " << name << std::endl;
}

//void print_diff(BRDF& brdf, LCT& lct, glm::vec3 v, int n = 128) {
//	for (int i = 0; i <= n; ++i) {
//		for (int j = 0; j <= n; ++j) {
//			float theta = glm::half_pi<float>() * glm::clamp((float)i / n, 0.0f, 1.0f);
//			float phi = glm::two_pi<float>() * glm::clamp((float)j / n, 0.0f, 1.0f);
//			glm::vec3 l = UVA::spherical(theta, phi);
//			LCT::EvalResult lct_eval = lct.eval(l);
//			BRDF::EvalResult brdf_eval = brdf.eval(v, l);
//			
//			std::cout << glm::clamp((float)i / n, 0.0f, 1.0f) * 0.5f << ", ";
//			std::cout << glm::clamp((float)j / n, 0.0f, 1.0f) * 2.0f << ", ";
//			std::cout << brdf_eval.brdfxcos_l - lct_eval.brdfxcos_l << std::endl;
//		}
//	}
//}


#pragma once
