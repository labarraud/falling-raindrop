
#ifndef LINALG_FILE_ALLOCATOR_HXX

namespace linalg
{
  using namespace std;
  
  //! Allocator using malloc/free operators
  template<class T>
  class MallocAlloc
  {
  public:
    static T* allocate(size_t num);
    static void deallocate(T* data, size_t num);
    static void* reallocate(T* data, size_t num);
    static void* resize(T* data, size_t old_n, size_t num);
    static void memoryset(T* data, char c, size_t num);
    static void memorycpy(T* datat, T* datas, size_t num);
  };


  //! Allocator using new/delete operators
  template<class T>
  class NewAlloc
  {
  public:
    static T* allocate(size_t num);
    static void deallocate(T* data, size_t num);    
    static void* reallocate(T* data, size_t num);
    static void* resize(T* data, size_t old_n, size_t num);
    static void memoryset(T* data, char c, size_t num);
    static void memorycpy(T* datat, T* datas, size_t num);
  };
  
  
  //! Default allocator is NewAlloc
  template<class T>
  class DefaultAllocator
  {
  public:
    typedef NewAlloc<T> allocator;
  };

  
  // For basic types, MallocAlloc is prefered
  
  template<>
  class DefaultAllocator<bool>
  {
  public:
    typedef MallocAlloc<bool> allocator;
  };


  template<>
  class DefaultAllocator<int>
  {
  public:
    typedef MallocAlloc<int> allocator;
  };


  template<>
  class DefaultAllocator<unsigned>
  {
  public:
    typedef MallocAlloc<unsigned> allocator;
  };


  template<>
  class DefaultAllocator<float>
  {
  public:
    typedef MallocAlloc<float> allocator;
  };


  template<>
  class DefaultAllocator<double>
  {
  public:
    typedef MallocAlloc<double> allocator;
  };


  template<>
  class DefaultAllocator<complex<float> >
  {
  public:
    typedef MallocAlloc<complex<float> > allocator;
  };


  template<>
  class DefaultAllocator<complex<double> >
  {
  public:
    typedef MallocAlloc<complex<double> > allocator;
  };


  class Error
  {
  protected:
    //! Message describing the exception type.
    string description_;
    //! Name of the function in which the error occurred.
    string function_;
    //! A comment about the error.
    string comment_;

  public:
    Error(string function = "", string comment = "");
    Error(string description, string function, string comment);
    virtual ~Error();

    void CoutWhat();
  };


  class NoMemory: public Error
  {
  public:
    NoMemory(string function = "", string comment = "");
  };


  class WrongIndex: public Error
  {
  public:
    WrongIndex(string function = "", string comment = "");
  };

  class IOError: public Error
  {
  public:
    IOError(string function = "", string comment = "");
  };

}

#define LINALG_FILE_ALLOCATOR_HXX
#endif
