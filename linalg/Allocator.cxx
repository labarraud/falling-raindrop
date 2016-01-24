
#ifndef LINALG_FILE_ALLOCATOR_CXX

#include "Allocator.hxx"

namespace linalg
{
  
  /* MallocAlloc */
  
  
  template <class T>
  inline T* MallocAlloc<T>::allocate(size_t num)
  {
    return static_cast<T*>( malloc(num * sizeof(T)) );
  }

  template <class T>
  inline void MallocAlloc<T>::deallocate(T* data, size_t num)
  {
    free(data);
  }

  template <class T>
  inline void* MallocAlloc<T>::reallocate(T* data, size_t num)
  {
    return realloc(reinterpret_cast<void*>(data), num * sizeof(T));
  }

  template <class T>
  inline void* MallocAlloc<T>::resize(T* data, size_t old_n, size_t num)
  {
    return realloc(reinterpret_cast<void*>(data), num * sizeof(T));
  }

  template <class T>
  inline void MallocAlloc<T>::memoryset(T* data, char c, size_t num)
  {
    memset(reinterpret_cast<void*>(data), c, num);
  }

  template <class T>
  inline void
  MallocAlloc<T>::memorycpy(T* datat, T* datas, size_t num)
  {
    memcpy(reinterpret_cast<void*>(datat), reinterpret_cast<void*>(datas),
	   num * sizeof(T));
  }
  
  
  /* NewAlloc */


  template <class T>
  inline T* NewAlloc<T>::allocate(size_t num)
  {
    return static_cast<T*>(new T[num]);
  }

  template <class T>
  inline void NewAlloc<T>::deallocate(T* data, size_t num)
  {
    delete [] data;
  }

  template <class T>
  inline void* NewAlloc<T>::reallocate(T* data, size_t num)
  {
    if (data != NULL)
      delete [] data;

    return (new T[num]);
  }

  template <class T>
  inline void* NewAlloc<T>::resize(T* data, size_t old_n, size_t num)
  {
    T* new_data = new T[num];
    for (size_t i = 0; i < min(old_n, num); i++)
      new_data[i] = data[i];
    
    if (data != NULL)
      delete [] data;

    return new_data;
  }

  template <class T>
  inline void NewAlloc<T>::memoryset(T* data, char c, size_t num)
  {
    memset(reinterpret_cast<void*>(data), c, num);
  }

  template <class T>
  inline void
  NewAlloc<T>::memorycpy(T* datat, T* datas, size_t num)
  {
    for (size_t i = 0; i < num; i++)
      datat[i] = datas[i];
  }


  /* Error */
  
  
  Error::Error(string function, string comment):
    description_("ERROR!\nAn undefined error occurred"),
    function_(function), comment_(comment)
  {
    this->CoutWhat();
    abort();
  }


  Error::Error(string description, string function, string comment):
    description_("ERROR!\n" + description),
    function_(function), comment_(comment)
  {
    this->CoutWhat();
    abort();
  }

  
  Error::~Error()
  {
  }


  void Error::CoutWhat()
  {
    cout << this->description_ << " in " << this->function_
	 << ".\n" << this->comment_ << endl;
  }
  

  NoMemory::NoMemory(string function, string comment):
    Error("Out of memory", function, comment)
  {
  }
  

  WrongIndex::WrongIndex(string function, string comment):
    Error("Index out of range", function, comment)
  {
  }


  IOError::IOError(string function, string comment):
    Error("Error while performing an I/O operation", function, comment)
  {
  }

}

#define LINALG_FILE_ALLOCATOR_CXX
#endif
