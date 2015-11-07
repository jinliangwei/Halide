#include "Halide.h"

#include <iostream>

const int kRank = 8;
const int kHeight = 4;
const int kWidth = 4;
const float step_size = 0.1;

int main(int argc, char *argv[]) {

  Halide::Func data_mat;
  Halide::Var row, col, ki, iter;
  //Halide::Expr e = row + col;
  // initialize data matrix
  data_mat(row, col) = 1.0f;

  // initialize model matrix L, R
  Halide::Func L, R;
  L(row, ki) = 1.0f;
  R(ki, col) = 1.0f;

  Halide::Func approx_mat;
  approx_mat(row, col) = 0.0f;

  Halide::RDom k(0, kRank);
  approx_mat(row, col) = approx_mat(row, col) + L(row, k) * R(k, col);

  Halide::Func diff_mat;
  diff_mat(row, col) = approx_mat(row, col) - data_mat(row, col);

  Halide::Func Ln, Rn;
  Ln(row, ki) = 1.0f;
  Rn(ki, col) = 1.0f;
  Halide::RDom c(0, kWidth);
  Ln(row, ki) = L(row, ki) - step_size * diff_mat(row, c) * R(ki, c);
  Halide::RDom r(0, kHeight);
  Rn(ki, col) = R(ki, col) - step_size * diff_mat(r, col) * L(r, ki);

  Halide::Image<float> Ln_output = Ln.realize(kHeight, kRank);
  Halide::Image<float> Rn_output = Rn.realize(kRank, kWidth);

  for (int i = 0; i < Ln_output.height(); i++) {
    std::cout << i << " - ";
    for (int j = 0; j < Ln_output.width(); j++) {
      std::cout << j <<":" << Ln_output(i, j) << " ";
    }
    std::cout << std::endl;
  }

  for (int i = 0; i < Rn_output.height(); i++) {
    std::cout << i << " - ";
    for (int j = 0; j < Rn_output.width(); j++) {
      std::cout << j <<":" << Rn_output(i, j) << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "Success!" << std::endl;
  return 0;
}
