#pragma once
#include <utility>
#include <memory>
#include <cassert>

namespace task {

template <class T>
class UniquePtr {
 private:
  T* data_ptr = nullptr;

 public:
  using pointer = T*;
  UniquePtr() = default;
  UniquePtr(const UniquePtr& obj) = delete;
  UniquePtr& operator=(const UniquePtr& obj) = delete;
  explicit UniquePtr(pointer obj_ptr): data_ptr(obj_ptr) {}
  UniquePtr(UniquePtr&& another) {
    data_ptr = another.data_ptr;
    another.data_ptr = nullptr;
  }

  void swap(UniquePtr& other) noexcept {
    std::swap(this->data_ptr, other.data_ptr);
  }

  UniquePtr& operator=(UniquePtr&& obj) {
    data_ptr = obj.data_ptr;
    obj.data_ptr = nullptr;
    return *this;
  }

  pointer get() const noexcept {
    return data_ptr;
  }

  T& operator*() {return *get();}
  pointer operator->() {return get();}

  pointer release() noexcept {
    pointer p = get();
    data_ptr = nullptr;
    return p;
  }

  void reset(T* ptr) {
      delete data_ptr;
      *this = UniquePtr<T>(ptr);
  }

  ~UniquePtr() {
    delete data_ptr;
  }
};

template<class T>
struct Control_block {
  T* data_ptr;
  int32_t shared_counter;
  int32_t weak_counter;

  Control_block(T* ptr):
    data_ptr(ptr),
    shared_counter(1),
    weak_counter(1) {}
};

template<class T>
class WeakPtr;
template <class T>
class SharedPtr {
  friend class WeakPtr<T>;
 private:
  Control_block<T>* block;

 public:
  SharedPtr<T>() {
    block = new Control_block<T>(nullptr);
  }

  SharedPtr<T>(T* ptr) {
    block = new Control_block<T>(ptr);
  }

  SharedPtr<T>(const SharedPtr<T>& another) {
    block = another.block;
    (block->shared_counter)++;
  }

  SharedPtr<T>(const WeakPtr<T>& w_p):SharedPtr<T>() {
    if (!w_p.expired()) {
      block = w_p.block;
      (block->shared_counter)++;
    }
  }

  SharedPtr<T>(SharedPtr<T>&& another) {
    block = std::move(another.block);
    another.block = nullptr;
  }

  SharedPtr<T>& operator=(const SharedPtr<T>& another) {
    if (&another == this)
      return *this;
    (another.block->shared_counter)++;
    (block->shared_counter)--;
    if (use_count() == 0)
      this->~SharedPtr();
    block = another.block;
    return *this;
  }

  bool operator==(const SharedPtr<T>& another) {
    return get() == another.get();
  }

  void swap(SharedPtr<T>& another) {
    std::swap(block, another.block);
  }

  SharedPtr<T>& operator=(SharedPtr<T>&& another) {
    if (&another == this)
      return *this;
    another.swap(*this);
    return *this;
  }

  ~SharedPtr<T>() {
    if (!this->block)
      return;
    (block->shared_counter)--;
    if (use_count() == 0) {
      delete block->data_ptr;
      block->data_ptr = nullptr;
    }
    if (use_count() == 0 && block->weak_counter == 0)
      delete block;
  }

  void reset(T* ptr = nullptr) {
    SharedPtr<T>(ptr).swap(*this);
  }

  T* get() const noexcept {
    return block->data_ptr;
  }

  int32_t use_count() const noexcept {
    return block->shared_counter;
  }

  T& operator*() {return *get();}
  T* operator->() {return get();}
};

template <class T>
class WeakPtr {
  friend class SharedPtr<T>;
 private:
  Control_block<T>* block;

 public:
  WeakPtr<T>(): block(nullptr) {}
  WeakPtr<T>(const SharedPtr<T>& sh_p) {
    block = sh_p.block;
    (block->weak_counter)++;
  }

  WeakPtr<T>(const WeakPtr<T>& another) {
    block = another.block;
    (block->weak_counter)++;
  }

  WeakPtr<T>(WeakPtr<T>&& another): WeakPtr<T>() {
    *this = std::move(another);
  }

  WeakPtr<T>& operator=(const WeakPtr<T>& another) {
    if (this == &another)
      return *this;
    ~WeakPtr<T>();
    block = another.block;
    (block->weak_counter)++;
    return *this;
  }

  WeakPtr<T>& operator=(WeakPtr<T>&& another) {
    if (this == &another)
      return *this;
    block = another.block;
    another.block = nullptr;
    return *this;
  }

  WeakPtr<T>& operator=(const SharedPtr<T>& sh_p) {
    block = sh_p.block;
    (block->weak_counter)++;
    return *this;
  }

  int32_t use_count() const noexcept {
    if (block)
      return block->shared_counter;
    else
      return 0;
  }

  SharedPtr<T> lock() const {
      return (expired() ? SharedPtr<T>() : SharedPtr<T>(*this));
  }

  bool expired() const noexcept {
    return use_count() == 0;
  }

  ~WeakPtr<T>() {
    if (!block)
      return;
    (block->weak_counter)--;
    if (use_count() == 0 && block->weak_counter == 0)
      delete block;
  }
};

// Your template function definitions may go here...

}  // namespace task

#include "smart_pointers.tpp"
