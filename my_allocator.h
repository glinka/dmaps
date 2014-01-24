#ifndef MY_ALLOCATOR_H
#define MY_ALLOCATOR_H

#include <limits>
#include <memory>

namespace extensions {
  template <typename T>
  class my_allocator : public std::allocator<T> {
  public:
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T value_type;
  public:
    my_allocator(const value_type* firstval);
    pointer allocate(size_type n, allocator<U>::const_pointer hint=0);
  private:
    const T* firstval;
    bool IS_FIRST_ALLOC;
  };
}

#endif
