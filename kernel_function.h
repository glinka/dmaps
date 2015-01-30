#ifndef KERNEL_FUNCTION_H
#define KERNEL_FUNCTION_H

template <typename T>
class Kernel_Function {
 public:
  virtual double kernel(const T& x1, const T& x2) const = 0;
};

#endif
