
#ifndef LINALG_FILE_VECTOR_HXX

#include "Allocator.hxx"

namespace linalg
{
    
  //! Class for basic vectors
  template<class T, class Allocator = 
	   typename DefaultAllocator<T>::allocator>
  class Vector
  {
  protected:
    //! Number of elements
    size_t m_;
    //! Pointer to stored elements
    T* data_;
    
  public:
    Vector();
    explicit Vector(size_t n);
    Vector(const Vector<T, Allocator>& A);
    
    ~Vector();
    void Clear();
    
    int GetM() const;
    size_t GetSize() const;
    
    T* GetData() const;
    void* GetDataVoid() const;
    
    void Reallocate(size_t);
    void Resize(size_t);
    
    void SetData(size_t, T*);
    void Nullify();
    
    T& operator()(size_t);
    const T& operator()(size_t) const;
    
    Vector<T, Allocator>& operator=(const Vector<T, Allocator>&);
    Vector<T, Allocator>& operator*=(const T&);

    void PushBack(const T& x);
    void PushBack(const Vector<T, Allocator>& x);
    
    void Zero();
    void Fill(const T& x);
    void FillRand();

    void Write(string FileName, bool with_size = true) const;
    void Write(ostream& FileStream, bool with_size = true) const;    
    void WriteText(const string&) const;
    void WriteText(ostream&) const;
    
    void ReadText(const string&);
    void ReadText(istream&);
    
  };

  template<class T>
  void GetRand(T& x);

  template<class T>
  void GetRand(complex<T>& x);

  template<class T>
  T DotProd(const Vector<T>& x, Vector<T>& y);

  template<class T>
  T Norm2(const Vector<T>& x);

  template<class T>
  T Norm2(const Vector<complex<T> >& x);
  
  template<class T>
  void Add(const T& alpha, const Vector<T>& x, Vector<T>& y);
  
  template<class T, class Allocator>
  ostream& operator<<(ostream& out, const Vector<T, Allocator>&);
  
}



#define LINALG_FILE_VECTOR_HXX
#endif
