#include "my_allocator.h"

my_allocator::my_allocator(const value_type* firstval): firstval(firstval), IS_FIRST_ALLOC(false) {};
}

pointer my_allocator::allocate(size_type n, allocator<U>::const_pointer hint=0) {
  if(IS_FIRST_ALLOC) {
    return firstval;
    IS_FIRST_ALLOC = false;
  }
  else {
    return std::allocator<T>::allocate(n, hint);
  }
}

void my_allocator::deallocate(pointer p) {
  if(p != firstval) {
    std::allocator<T>::deallocate(pointer);
  }
}
