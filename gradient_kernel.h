#ifndef _GRADIENT_KERNEL_H_
#define _GRADIENT_KERNEL_H_

#include <math.h>
#include <Eigen/Dense>
#include "gradient.h" /// requires user-defined gradient functor

typedef Eigen::VectorXd Vector;

/**
 * \class Kernel_Function
 *
 * Provides a kernel for use in DMAPS described in Lafon's thesis. It uses gradient information to produce embeddings that are constant on level sets of \f$f(x)\f$. The kernel is given by \f$k(x,y)=e^{\frac{-\| x-y \|^2}{\epsilon} - \frac{<\nabla f(x), x-y>^2}{\epsilon^2}}\f$
 */

class Kernel_Function {
 public:
  /// constructor
 Kernel_Function(const double epsilon, const Gradient gradient = Gradient()):_epsilon(epsilon), _gradient(gradient) {}
  /// copy constructor
 Kernel_Function(const Kernel_Function& gk):_epsilon(gk._epsilon)  {}
  /// move constructor
 Kernel_Function(Kernel_Function&& gk): _epsilon(std::move(gk._epsilon)) {}
  /* no assignment operator, only const members */
  ~Kernel_Function() {}
  /**
   * Calculates the anisotropic, gradient kernel between two vectors
   *
   * \param x1 first vector
   * \param x2 second vector
   * \returns \f$k(x,y)=e^{\frac{-\| x-y \|^2}{\epsilon} - \frac{<\nabla f(x), x-y>^2}{\epsilon^2}}\f$
   */
  double operator()(const Vector& x1, const Vector& x2) const {
    return std::exp(-(x1-x2).norm()/_epsilon - pow(_gradient(x1).dot(x1-x2),2)/(_epsilon*_epsilon));
  }
 private:
  const double _epsilon; ///< DMAPS parameter defining a points neighborhood: only those points within approximately distance _epsilon will be considered neighbors
  const Gradient _gradient;
};

#endif /* _GRADIENT_KERNEL_H_ */
