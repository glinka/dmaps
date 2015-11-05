#ifndef _GAUSSIAN_KERNEL_H_
#define _GAUSSIAN_KERNEL_H_

/**
 * \class Kernel_Function
 *
 * Provides a Gaussian kernel for use in DMAPS, given by \f$k(x,y)=e^{\frac{-\| x-y \|^2}{\epsilon^2}}\f$. This is a standard choice when datapoints are themselves vectors in \f$R^n\f$.
 */


class Kernel_Function {
 public:
  /// constructor
 Kernel_Function(const double epsilon):_epsilon(epsilon) {}
  /// copy constructor
 Kernel_Function(const Kernel_Function& gk):_epsilon(gk._epsilon)  {}
  /// move constructor
 Kernel_Function(Kernel_Function&& gk): _epsilon(std::move(gk._epsilon)) {}
  /* no assignment operator, only const members */
  ~Kernel_Function() {}
  /**
   * Calculates the Gaussian kernel between two vectors
   *
   * \param x1 first vector
   * \param x2 second vector
   * \returns \f$k(x,y)=e^{\frac{-\| x-y \|^2}{\epsilon^2}}\f$
   */
  double operator()(const std::vector<double>& x1, const std::vector<double>& x2) const {
    int n = x1.size();
    double norm = 0;
    for(int i = 0; i < n; i++) {
      norm += std::pow(x1[i] - x2[i], 2);
    }
    return std::exp(-norm/(_epsilon*_epsilon));
  }
 private:
  const double _epsilon; ///< DMAPS parameter defining a points neighborhood: only those points within approximately distance _epsilon will be considered neighbors
};


#endif /* _GAUSSIAN_KERNEL_H_ */


