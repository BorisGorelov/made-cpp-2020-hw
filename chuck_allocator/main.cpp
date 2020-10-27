#include <iostream>
#include <limits>
#include <memory>
#include <vector>

template<class T>
class ChuckAllocator {
 public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = size_t;
    using difference_type = ptrdiff_t;

    template<class U>
    struct rebind {
        typedef ChuckAllocator<U> other;
    };

    ChuckAllocator() = default;

    template<class U>
    ChuckAllocator(const ChuckAllocator<U> &other) {}

    ~ChuckAllocator() = default;

    pointer allocate(size_type numObjects) {
        alloc_counter += numObjects;
        return static_cast<pointer>(operator new(sizeof(T) * numObjects));
    }

    void deallocate(pointer p, size_type numOfObjs) {
        delete(p);
    }

    size_type get_alloc_counter() const {
        return alloc_counter;
    }

    template <class U, class... Args>
    void construct(U* obj, Args&& ...args) {
        new(obj) U(std::forward<Args>(args)...);
    }

    template<class U>
    void destroy(U* obj) {
        obj->~T();
    }

 private:
    static size_type alloc_counter;
};

template<class T>
typename ChuckAllocator<T>::size_type ChuckAllocator<T>::alloc_counter = 0;


int main() {
    std::vector<int, ChuckAllocator<int>> v(10);
    return 0;
}
