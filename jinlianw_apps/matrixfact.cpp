#include "Halide.h"

#include <iostream>

const int kRank = 100;
const int kHeight = 16;
const int kWidth = 16;
const float step_size = 0.1;

const int kIters = 5;

int main(int argc, char *argv[]) {

  Halide::Func data_mat;
  Halide::Var row, col, ki, iter;
  //Halide::Expr e = row + col;
  // initialize data matrix
  data_mat(row, col) = 0.0f;

  // initialize model matrix L, R
  Halide::Func L, R;
  // Normally, this should be random initialization,
  // but this is not available in Halide.
  L(row, ki, iter) = 1.0f;
  R(ki, col, iter) = 1.0f;

  Halide::Func approx_mat;
  Halide::RDom k(0, kRank);
  approx_mat(row, col, iter) += L(row, k, iter) * R(k, col, iter);

  Halide::Func diff_mat;
  diff_mat(row, col, iter) = approx_mat(row, col, iter) - data_mat(row, col);

  Halide::RDom c(0, kWidth), i(1, kIters);
  L(row, ki, i) = L(row, ki, i - 1) - step_size * diff_mat(row, c, i - 1) * R(ki, c, i - 1);
  Halide::RDom r(0, kHeight);
  R(ki, col, i) = R(ki, col, i - 1) - step_size * diff_mat(r, col, i - 1) * L(r, ki, i - 1);

  L.realize(kHeight, kRank, kIters);
  R.realize(kRank, kWidth, kIters);
  std::cout << "Success!" << std::endl;
  return 0;
}
