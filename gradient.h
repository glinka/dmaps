#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include <Eigen/Dense>

typedef Eigen::VectorXd Vector;

/**
 * \class Gradient
 *
 * Customizable class that computes the gradient of some functional
 */

class Gradient {
 public:
  Gradient() {}
  ~Gradient() {}
  /**
   * Calculates the gradient of some function, in this case \f$f(x,y) = x^2 + y^2\f$
   *
   * \param x vector at which to evaluate the gradient
   * \returns vector containing gradient at x, in this case \f$\nabla f(x) = \(x[0], x[1]\) \f$
   */
  Vector operator()(const Vector& x) const {
    return x;
  }
};

#endif /* _GRADIENT_H_ */
